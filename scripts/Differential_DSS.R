# ============================================================================ #
# SCRIPT:
#   Differential methylation analysis (DSS) for one chromosome
# DESCRIPTION:
#   Loads BSBolt CGmap files for one chromosome, builds a BSseq object, keeps
#   CpGs covered in all samples and calls DMLs/DMRs between two groups of a
#   metadata variable. diff = mean(test_group) - mean(reference group).
# USAGE (standalone):
#   Rscript Differential_DSS.R \
#     --metadata metadata.tsv --samples_col SampleID \
#     --indir <BSBolt 03_CallMethylation dir> --outdir <out dir> \
#     --chrom chr1 --var KRAS --test_group MUT [--min_cov 5]
# ============================================================================ #

# Redirect output to the log file when running under Snakemake
if (exists("snakemake") && length(snakemake@log) != 0) {
  dir.create(dirname(snakemake@log[[1]]), recursive = TRUE, showWarnings = FALSE)
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log, append = TRUE)
  sink(log, append = TRUE, type = "message")
}

# Load packages (install if missing: CRAN first, Bioconductor as fallback)
load_pkgs <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(p, update = FALSE, ask = FALSE)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}
load_pkgs(c("optparse", "tidyverse", "data.table", "bsseq", "DSS"))

set.seed(42)


# ============================================================================ #
#                         Setup command line arguments                         #
# ============================================================================ #

if (exists("snakemake")) {
  opt <- list(
    metadata    = snakemake@input[["metadata"]],
    samples     = unlist(snakemake@params[["samples"]]),
    samples_col = snakemake@params[["samples_col"]],
    indir       = snakemake@params[["indir"]],
    outdir      = snakemake@params[["outdir"]],
    chrom       = snakemake@wildcards[["chrom"]],
    var         = snakemake@params[["var"]],
    test_group  = snakemake@params[["test_group"]],
    min_cov     = snakemake@params[["min_cov"]]
  )
  out_bsseq <- snakemake@output[["bsseq"]]
  out_dml   <- snakemake@output[["dmls"]]
  out_dmr   <- snakemake@output[["dmrs"]]
} else {
  option_list <- list(
    make_option(c("-m", "--metadata"),    type = "character"),
    make_option(c("-S", "--samples_col"), type = "character", default = "SampleID"),
    make_option(c("-d", "--indir"),       type = "character"),
    make_option(c("-o", "--outdir"),      type = "character"),
    make_option(c("-c", "--chrom"),       type = "character"),
    make_option(c("-v", "--var"),         type = "character"),
    make_option(c("-t", "--test_group"),  type = "character"),
    make_option(c("--min_cov"),           type = "integer", default = 5)
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  # Standalone: use every sample directory found in indir
  opt$samples <- list.dirs(opt$indir, full.names = FALSE, recursive = FALSE)
  out_bsseq <- file.path(opt$outdir, paste0("BSseq_object_", opt$chrom, ".rds"))
  out_dml   <- file.path(opt$outdir, paste0("DML_results_", opt$var, "_", opt$chrom, ".csv"))
  out_dmr   <- file.path(opt$outdir, paste0("DMR_results_", opt$var, "_", opt$chrom, ".csv"))
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)


# ============================================================================ #
#                         Metadata and group definition                        #
# ============================================================================ #

metadata <- read_tsv(opt$metadata, show_col_types = FALSE)
for (col in c(opt$samples_col, opt$var)) {
  if (!col %in% colnames(metadata)) stop(paste0("Column '", col, "' not found in metadata."))
}

# Keep samples processed by the pipeline with a non-missing grouping value
metadata <- metadata %>%
  mutate(across(all_of(c(opt$samples_col, opt$var)), as.character)) %>%
  filter(.data[[opt$samples_col]] %in% opt$samples, !is.na(.data[[opt$var]]))

groups <- unique(metadata[[opt$var]])
if (length(groups) != 2) {
  stop(paste0("Variable '", opt$var, "' must have exactly 2 non-missing levels; found: ",
              paste(groups, collapse = ", ")))
}
if (!opt$test_group %in% groups) {
  stop(paste0("test_group '", opt$test_group, "' is not a level of '", opt$var, "'."))
}
ref_group <- setdiff(groups, opt$test_group)

sample_ids <- metadata[[opt$samples_col]]
group1 <- sample_ids[metadata[[opt$var]] == opt$test_group]  # test
group2 <- sample_ids[metadata[[opt$var]] == ref_group]       # reference
cat("Chromosome:", opt$chrom, "\n")
cat("Test group (", opt$test_group, "):", length(group1), "samples\n")
cat("Reference group (", ref_group, "):", length(group2), "samples\n")


# ============================================================================ #
#                         Load CGmap files (one chromosome)                    #
# ============================================================================ #

# BSBolt CGmap columns: chr, nucleotide, pos, context, subcontext, meth_level, meth_reads, total_reads
load_bsbolt_cgmap <- function(sample_id) {
  file_path <- file.path(opt$indir, sample_id, paste0(sample_id, "_meth.CGmap.gz"))
  if (!file.exists(file_path)) stop(paste("File not found:", file_path))
  cat("Processing sample:", sample_id, "\n")

  # Stream-filter by chromosome and CG context to avoid loading the whole genome
  cmd <- sprintf("zcat '%s' | awk -F'\\t' '$1 == \"%s\" && $4 == \"CG\"'", file_path, opt$chrom)
  tmp <- suppressWarnings(fread(cmd = cmd, header = FALSE, sep = "\t"))

  if (nrow(tmp) == 0) {
    return(data.frame(chr = character(), pos = numeric(), N = numeric(), X = numeric()))
  }
  data.frame(chr = tmp$V1, pos = as.numeric(tmp$V3),
             N = as.numeric(tmp$V8), X = as.numeric(tmp$V7))
}

samples_list <- lapply(sample_ids, load_bsbolt_cgmap)

# Helpers to write empty (header-only) results so the workflow can continue
empty_dml <- data.frame(chr = character(), pos = integer(), mu1 = numeric(), mu2 = numeric(),
                        diff = numeric(), diff.se = numeric(), stat = numeric(),
                        phi1 = numeric(), phi2 = numeric(), pval = numeric(), fdr = numeric(),
                        postprob.overThreshold = numeric())
empty_dmr <- data.frame(chr = character(), start = integer(), end = integer(),
                        length = integer(), nCG = integer(), meanMethy1 = numeric(),
                        meanMethy2 = numeric(), diff.Methy = numeric(), areaStat = numeric())

if (any(vapply(samples_list, nrow, integer(1)) == 0)) {
  cat("At least one sample has no CpGs on", opt$chrom, "- writing empty results.\n")
  saveRDS(NULL, out_bsseq)
  write.csv(empty_dml, out_dml, row.names = FALSE)
  write.csv(empty_dmr, out_dmr, row.names = FALSE)
  quit(save = "no", status = 0)
}


# ============================================================================ #
#                    BSseq object and differential testing                     #
# ============================================================================ #

BSobj <- makeBSseqData(samples_list, sample_ids)

# Keep CpGs with coverage >= min_cov in every sample
BSobj <- BSobj[which(rowSums(getCoverage(BSobj) >= opt$min_cov) == ncol(BSobj)), ]
cat("CpGs covered >=", opt$min_cov, "x in all samples:", nrow(BSobj), "\n")
saveRDS(BSobj, file = out_bsseq)

if (nrow(BSobj) == 0) {
  cat("No CpGs left after coverage filtering - writing empty results.\n")
  write.csv(empty_dml, out_dml, row.names = FALSE)
  write.csv(empty_dmr, out_dmr, row.names = FALSE)
  quit(save = "no", status = 0)
}

# group1 = test, group2 = reference -> diff > 0 means hypermethylated in test
dmlTest <- DMLtest(BSobj, group1 = group1, group2 = group2, smoothing = TRUE)

dmls <- callDML(dmlTest, p.threshold = 0.05, delta = 0.1)
dmrs <- callDMR(dmlTest, p.threshold = 0.05, delta = 0.1)

# callDML/callDMR return NULL when nothing is significant
if (is.null(dmls) || nrow(dmls) == 0) dmls <- empty_dml
if (is.null(dmrs) || nrow(dmrs) == 0) dmrs <- empty_dmr

write.csv(dmls, out_dml, row.names = FALSE)
write.csv(dmrs, out_dmr, row.names = FALSE)

cat("DMLs:", nrow(dmls), "| DMRs:", nrow(dmrs), "\n")
cat("Differential methylation analysis completed successfully.\n")
