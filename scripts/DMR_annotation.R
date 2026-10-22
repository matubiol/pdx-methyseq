# ============================================================================ #
# SCRIPT:
#   DMR filtering, annotation and functional enrichment
# DESCRIPTION:
#   Filters genome-wide DSS DMRs, annotates them with annotatr (CpG context and
#   gene features), plots summary distributions, runs GO BP / KEGG enrichment
#   on genes with promoter DMRs (split by direction) and exports a BED file
#   for HOMER motif analysis.
#   Direction is relative to the test group (DSS diff = test - reference).
# USAGE (standalone):
#   Rscript DMR_annotation.R \
#     --dmrs DMR_results_<var>_all.csv --gen_ver hg38 --test MUT \
#     --outdir <output_directory> [--out_ann <annotated_csv>]
# ============================================================================ #

# Redirect output to the log file when running under Snakemake
if (exists("snakemake") && length(snakemake@log) != 0) {
  dir.create(dirname(snakemake@log[[1]]), recursive = TRUE, showWarnings = FALSE)
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log, append = TRUE)
  sink(log, append = TRUE, type = "message")
}

# Load packages (install if missing via BiocManager, which also handles CRAN)
load_pkgs <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(p, update = FALSE, ask = FALSE)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}


# ============================================================================ #
#                         Setup command line arguments                         #
# ============================================================================ #

if (exists("snakemake")) {
  load_pkgs("optparse")
  opt <- list(
    dmrs    = snakemake@input[["dmrs"]],
    gen_ver = snakemake@params[["gen_ver"]],
    test    = snakemake@params[["test"]],
    outdir  = snakemake@params[["outdir"]],
    out_ann = snakemake@output[["dmrs_ann"]]
  )
} else {
  load_pkgs("optparse")
  option_list <- list(
    make_option(c("-i", "--dmrs"),    type = "character", help = "Path to DMR results CSV file"),
    make_option(c("-g", "--gen_ver"), type = "character", default = "hg38",
                help = "Genome version (hg38 or hg19)"),
    make_option(c("-t", "--test"),    type = "character", help = "Test group name"),
    make_option(c("-o", "--outdir"),  type = "character", help = "Output directory"),
    make_option(c("-a", "--out_ann"), type = "character", default = NULL,
                help = "Annotated DMR table (default: <outdir>/DMR_annotated.csv)")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$out_ann)) opt$out_ann <- file.path(opt$outdir, "DMR_annotated.csv")
}

if (!opt$gen_ver %in% c("hg38", "hg19")) stop("gen_ver must be 'hg38' or 'hg19'.")
txdb_pkg <- paste0("TxDb.Hsapiens.UCSC.", opt$gen_ver, ".knownGene")

load_pkgs(c("tidyverse", "scales", "GenomicRanges", "annotatr", "clusterProfiler",
            "org.Hs.eg.db", txdb_pkg))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

test_group <- opt$test
lab_hyper <- paste0("Hypermethylated_in_", test_group)
lab_hypo  <- paste0("Hypomethylated_in_", test_group)
dir_colors <- setNames(c("#E64B35", "#4DBBD5"), c(lab_hyper, lab_hypo))


# ============================================================================ #
#                         Import and filter DMRs                               #
# ============================================================================ #

dmrs_raw <- read_csv(opt$dmrs, show_col_types = FALSE)
cat("Total DMRs before filtering:", nrow(dmrs_raw), "\n")

dmrs_sig <- dmrs_raw %>%
  filter(
    nCG >= 5,                 # minimum CpGs per DMR
    abs(diff.Methy) >= 0.1,   # >= 10% methylation difference
    abs(areaStat) >= 100      # sufficient effect size
  ) %>%
  mutate(
    direction = if_else(diff.Methy > 0, lab_hyper, lab_hypo),
    start = as.integer(start),
    end = as.integer(end)     # DSS may write coordinates as floats
  )

cat("DMRs after filtering:", nrow(dmrs_sig), "\n")
cat("Hypermethylated in", test_group, ":", sum(dmrs_sig$direction == lab_hyper), "\n")
cat("Hypomethylated in",  test_group, ":", sum(dmrs_sig$direction == lab_hypo), "\n")

write_csv(dmrs_sig, file.path(opt$outdir, "DMR_filtered.csv"))

if (nrow(dmrs_sig) == 0) {
  cat("No DMRs passed the filters - writing an empty annotation table and stopping.\n")
  write_csv(dmrs_sig, opt$out_ann)
  quit(save = "no", status = 0)
}


# ============================================================================ #
#                      Convert DMRs to GRanges and annotate                    #
# ============================================================================ #

dmrs_gr <- GRanges(
  seqnames = dmrs_sig$chr,
  ranges = IRanges(start = dmrs_sig$start, end = dmrs_sig$end),
  strand = "*",
  diff.Methy = dmrs_sig$diff.Methy,
  areaStat = dmrs_sig$areaStat,
  direction = dmrs_sig$direction,
  nCG = dmrs_sig$nCG
)

annot_types <- paste0(opt$gen_ver, c(
  "_cpgs",        # CpG islands, shores, shelves, open sea
  "_basicgenes"   # 1-5kb upstream, promoters, 5'UTR, exons, introns, 3'UTR
))
annotations <- build_annotations(genome = opt$gen_ver, annotations = annot_types)

dmrs_annotated <- annotate_regions(regions = dmrs_gr, annotations = annotations,
                                   ignore.strand = TRUE, quiet = FALSE)
dmrs_ann_df <- as.data.frame(dmrs_annotated)

write_csv(dmrs_ann_df, opt$out_ann)
cat("Annotated DMR table saved:", opt$out_ann, "\n")


# ============================================================================ #
#                                Summary plots                                 #
# ============================================================================ #

# 1. Distribution across genomic features
p_features <- dmrs_ann_df %>%
  filter(!is.na(annot.type)) %>%
  mutate(feature = case_when(
    str_detect(annot.type, "promoter") ~ "Promoter",
    str_detect(annot.type, "1to5kb")   ~ "1-5 kb upstream",
    str_detect(annot.type, "5UTR")     ~ "5' UTR",
    str_detect(annot.type, "exon")     ~ "Exon",
    str_detect(annot.type, "intron")   ~ "Intron",
    str_detect(annot.type, "3UTR")     ~ "3' UTR",
    str_detect(annot.type, "island")   ~ "CpG Island",
    str_detect(annot.type, "shore")    ~ "CpG Shore",
    str_detect(annot.type, "shelf")    ~ "CpG Shelf",
    str_detect(annot.type, "inter")    ~ "Open Sea",
    TRUE ~ "Other"
  )) %>%
  count(direction, feature) %>%
  ggplot(aes(x = reorder(feature, n), y = n, fill = direction)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = dir_colors) +
  theme_minimal(base_size = 13) +
  labs(title = "DMR distribution across genomic features",
       x = NULL, y = "Number of DMRs", fill = "Direction")
ggsave(file.path(opt$outdir, "DMR_genomic_features.svg"), p_features, width = 9, height = 6)

# 2. CpG context distribution
p_cpg <- dmrs_ann_df %>%
  filter(str_detect(annot.type, "cpg")) %>%
  mutate(cpg_context = case_when(
    str_detect(annot.type, "island") ~ "CpG Island",
    str_detect(annot.type, "shore")  ~ "CpG Shore",
    str_detect(annot.type, "shelf")  ~ "CpG Shelf",
    TRUE ~ "Open Sea"
  )) %>%
  count(direction, cpg_context) %>%
  ggplot(aes(x = cpg_context, y = n, fill = direction)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = dir_colors) +
  theme_minimal(base_size = 13) +
  labs(title = "DMR distribution across CpG contexts",
       x = NULL, y = "Proportion", fill = "Direction")
ggsave(file.path(opt$outdir, "DMR_cpg_context.svg"), p_cpg, width = 7, height = 5)

# 3. Methylation difference distribution
p_diff <- ggplot(dmrs_sig, aes(x = diff.Methy, fill = direction)) +
  geom_histogram(bins = 40, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = dir_colors) +
  theme_minimal(base_size = 13) +
  labs(title = paste0("Distribution of methylation differences (", test_group, " - reference)"),
       x = "diff.Methy", y = "Count", fill = "Direction")
ggsave(file.path(opt$outdir, "DMR_diff_distribution.svg"), p_diff, width = 8, height = 5)


# ============================================================================ #
#             Extract genes in promoters for enrichment                        #
# ============================================================================ #

promoter_dmrs <- dmrs_ann_df %>%
  filter(str_detect(annot.type, "promoter"), !is.na(annot.symbol)) %>%
  dplyr::select(seqnames, start, end, diff.Methy, areaStat, direction,
                nCG, annot.symbol, annot.gene_id) %>%
  distinct(annot.symbol, direction, .keep_all = TRUE)

write_csv(promoter_dmrs, file.path(opt$outdir, "DMR_promoter_genes.csv"))
cat("Genes with DMRs in promoters:", nrow(promoter_dmrs), "\n")

get_genes <- function(lab) {
  promoter_dmrs %>% filter(direction == lab) %>%
    pull(annot.gene_id) %>% na.omit() %>% as.character() %>% unique()
}
genes_hyper <- get_genes(lab_hyper)
genes_hypo  <- get_genes(lab_hypo)

# Universe: all genes with an Entrez ID in org.Hs.eg.db
universe <- keys(org.Hs.eg.db, keytype = "ENTREZID")


# ============================================================================ #
#                     GO / KEGG enrichment (clusterProfiler)                   #
# ============================================================================ #

run_enrichment <- function(gene_ids, label, outdir, universe) {
  if (length(gene_ids) < 5) {
    cat("Not enough genes for enrichment in:", label, "\n")
    return(invisible(NULL))
  }

  ego <- enrichGO(gene = gene_ids, universe = universe, OrgDb = org.Hs.eg.db,
                  ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("No significant GO terms for:", label, "\n")
  } else {
    write_csv(as.data.frame(ego), file.path(outdir, paste0("GO_BP_", label, ".csv")))
    p <- dotplot(ego, showCategory = 20, font.size = 10) +
      labs(title = paste0("GO Biological Process - ", label)) +
      theme_minimal(base_size = 11)
    ggsave(file.path(outdir, paste0("GO_BP_dotplot_", label, ".svg")), p, width = 10, height = 8)
  }

  # KEGG needs internet access to the KEGG REST API; failures are not fatal
  ekegg <- tryCatch(
    enrichKEGG(gene = gene_ids, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05),
    error = function(e) { cat("KEGG enrichment failed for", label, ":", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    write_csv(as.data.frame(ekegg), file.path(outdir, paste0("KEGG_", label, ".csv")))
    p_kegg <- dotplot(ekegg, showCategory = 20, font.size = 10) +
      labs(title = paste0("KEGG Pathways - ", label)) +
      theme_minimal(base_size = 11)
    ggsave(file.path(outdir, paste0("KEGG_dotplot_", label, ".svg")), p_kegg, width = 10, height = 8)
  }

  invisible(ego)
}

run_enrichment(genes_hyper, lab_hyper, opt$outdir, universe)
run_enrichment(genes_hypo,  lab_hypo,  opt$outdir, universe)


# ============================================================================ #
#                Export BED file for HOMER motif analysis                      #
# ============================================================================ #

# BED6: chr, start (0-based), end, name, score, strand
dmrs_sig %>%
  mutate(start = start - 1L,
         name = paste0(direction, "_", row_number()),
         score = round(abs(areaStat)),
         strand = "+") %>%
  dplyr::select(chr, start, end, name, score, strand) %>%
  write_tsv(file.path(opt$outdir, "DMR_filtered.bed"), col_names = FALSE)

cat("\nAnnotation analysis completed. Results saved in:", opt$outdir, "\n")
