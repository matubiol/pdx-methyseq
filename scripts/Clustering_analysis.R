# ============================================================================ #
# SCRIPT:
#   Clustering analysis for methylation data
# DESCRIPTION:
#   PCA on highly variable CpGs of the imputed BSBolt matrix, optimal k
#   selection by silhouette (k-means, hierarchical, PAM), PERMANOVA validation
#   and 2D/3D PCA plots.
# USAGE (standalone):
#   Rscript Clustering_analysis.R \
#     --metadata <metadata_file> \
#     --matrix <methyl_matrix_imputed.tsv> \
#     --outdir <output_directory> \
#     --samples_col <sample_ID_column_in_metadata> \
#     [--batch_var <metadata_column_for_ComBat, empty to skip>] \
#     [-k <number_of_clusters, default: best>]
# ============================================================================ #

# Redirect output to the log file when running under Snakemake
if (exists("snakemake") && length(snakemake@log) != 0) {
  dir.create(dirname(snakemake@log[[1]]), recursive = TRUE, showWarnings = FALSE)
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log, append = TRUE)
  sink(log, append = TRUE, type = "message")
}

# Load packages (install if missing: CRAN/Bioconductor via BiocManager, GitHub for pairwiseAdonis)
load_pkgs <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (p == "pairwiseAdonis") {
        if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
        remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
      } else {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(p, update = FALSE, ask = FALSE)
      }
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}
load_pkgs(c("tidyverse", "PCAtools", "cluster", "plotly", "factoextra", "ggsci",
            "vegan", "pairwiseAdonis", "htmlwidgets", "sva", "optparse"))

set.seed(42)


# ============================================================================ #
#                         Setup command line arguments                         #
# ============================================================================ #

if (exists("snakemake")) {
  opt <- list(
    metadata    = snakemake@input[["metadata"]],
    matrix      = snakemake@input[["matrix"]],
    k           = snakemake@params[["k"]],
    batch_var   = snakemake@params[["batch_var"]],
    samples_col = snakemake@params[["samples_col"]],
    outdir      = snakemake@params[["outdir"]]
  )
} else {
  option_list <- list(
    make_option(c("-m", "--metadata"),    type = "character"),
    make_option(c("-i", "--matrix"),      type = "character"),
    make_option(c("-k"),                  type = "character", default = "best"),
    make_option(c("-b", "--batch_var"),   type = "character", default = ""),
    make_option(c("-s", "--samples_col"), type = "character", default = "SampleID"),
    make_option(c("-o", "--outdir"),      type = "character")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
}

opt$k <- as.character(opt$k)
use_best_k <- is.null(opt$k) || opt$k %in% c("", "best")
outdir_k <- file.path(opt$outdir, if (use_best_k) "k_best" else paste0("k_", opt$k))
dir.create(outdir_k, showWarnings = FALSE, recursive = TRUE)


# ============================================================================ #
#                         Import data and filter sites                         #
# ============================================================================ #

metadata <- read_tsv(opt$metadata, show_col_types = FALSE)
matrix <- read_tsv(opt$matrix, show_col_types = FALSE)

if (!opt$samples_col %in% colnames(metadata)) {
  stop(paste0("Column '", opt$samples_col, "' not found in metadata."))
}
metadata[[opt$samples_col]] <- as.character(metadata[[opt$samples_col]])

# First column holds the CpG site IDs; the rest are samples.
# Strip any path / "_meth.CGmap.gz" suffix in case labels came from file names.
site_col <- colnames(matrix)[1]
clean_id <- function(x) sub("_meth(\\.CGmap(\\.gz)?)?$", "", basename(x))
colnames(matrix)[-1] <- clean_id(colnames(matrix)[-1])

missing <- setdiff(colnames(matrix)[-1], metadata[[opt$samples_col]])
if (length(missing) > 0) {
  stop(paste("Samples in the matrix not present in the metadata:",
             paste(missing, collapse = ", ")))
}

# Samples as rows, sites as columns
data_pca <- matrix %>%
  column_to_rownames(site_col) %>%
  as.matrix() %>%
  t()

# Filter out sites with NAs and constant sites
nas_per_site <- colSums(is.na(data_pca))
data_no_nas <- data_pca[, nas_per_site == 0, drop = FALSE]
sd_sites <- apply(data_no_nas, 2, sd)
data_pca_final <- data_no_nas[, sd_sites > 0, drop = FALSE]

cat("Initial sites:", ncol(data_pca), "\n")
cat("Sites with NAs removed:", sum(nas_per_site > 0), "\n")
cat("Constant sites removed:", sum(sd_sites == 0), "\n")


# ============================================================================ #
#                    Batch correction (optional) via ComBat                    #
# ============================================================================ #
# Runs on all retained sites, before variance-based site selection.

if (!is.null(opt$batch_var) && nzchar(opt$batch_var)) {
  batch <- metadata[[opt$batch_var]][match(rownames(data_pca_final), metadata[[opt$samples_col]])]
  if (any(is.na(batch))) {
    stop(paste0("Some samples have no matching '", opt$batch_var, "' in metadata."))
  }
  # ComBat expects sites as rows, samples as columns
  data_combat <- ComBat(dat = t(data_pca_final), batch = batch, mod = NULL,
                        par.prior = TRUE, prior.plots = FALSE)
  data_pca_final <- t(data_combat)
  cat("Batch correction (ComBat) applied: batch =", opt$batch_var, "\n")
} else {
  cat("No batch_var provided - skipping ComBat correction.\n")
}


# ============================================================================ #
#                    Filter sites by excess variance (chi-square)              #
# ============================================================================ #
# Keeps sites whose variance significantly exceeds the median site variance.

site_var <- apply(data_pca_final, 2, var)
var_pool <- median(site_var)
n_samples <- nrow(data_pca_final)
chisq_stat <- (n_samples - 1) * site_var / var_pool
p_vals <- pchisq(chisq_stat, df = n_samples - 1, lower.tail = FALSE)
q_vals <- p.adjust(p_vals, method = "BH")
data_pca_final <- data_pca_final[, q_vals < 0.05, drop = FALSE]

cat("Sites retained after variance filtering (q < 0.05):", ncol(data_pca_final), "\n")
if (ncol(data_pca_final) < 3) stop("Too few variable sites left for PCA.")


# ============================================================================ #
#                              Run PCA                                         #
# ============================================================================ #

pca_res <- prcomp(data_pca_final, scale. = TRUE)

pca_df <- as.data.frame(pca_res$x)
pca_df[[opt$samples_col]] <- rownames(pca_df)
pca_df <- pca_df %>%
  left_join(metadata, by = opt$samples_col) %>%
  mutate(across(where(is.character),
                ~ factor(if_else(is.na(.), "Unknown", .),
                         levels = c(unique(na.omit(.)), "Unknown"))))

# Number of PCs to retain: scree elbow (at least 2)
variance <- (pca_res$sdev)^2 / sum(pca_res$sdev^2) * 100
elbow <- max(2, findElbowPoint(variance))
cat(paste0("Using PC1-PC", elbow, ", explaining ",
           round(sum(summary(pca_res)$importance[2, 1:elbow]) * 100, 2),
           "% of total variance cumulatively\n"))


# ============================================================================ #
#                       Find optimal k via silhouette loop                     #
# ============================================================================ #

pca_scores <- pca_res$x[, 1:elbow, drop = FALSE]
dist_pca <- dist(pca_scores)

k_range <- 2:min(10, nrow(pca_scores) - 1)
sil_kmeans <- sil_hclust <- sil_pam <- numeric(length(k_range))

# Hierarchical tree computed once; only the cut point (k) changes
hc_full <- hclust(dist_pca, method = "ward.D2")

for (i in seq_along(k_range)) {
  k <- k_range[i]
  set.seed(42)
  km <- kmeans(pca_scores, centers = k, nstart = 25)
  sil_kmeans[i] <- summary(silhouette(km$cluster, dist_pca))$avg.width
  sil_hclust[i] <- summary(silhouette(cutree(hc_full, k = k), dist_pca))$avg.width
  sil_pam[i] <- pam(dist_pca, k = k)$silinfo$avg.width
}

sil_df <- tibble(
  k = rep(k_range, 3),
  method = rep(c("kmeans", "hclust", "pam"), each = length(k_range)),
  avg_sil = c(sil_kmeans, sil_hclust, sil_pam)
)

# Best method/k by silhouette; optional manual k override (method chosen at that k)
if (use_best_k) {
  best_row <- sil_df[which.max(sil_df$avg_sil), ]
} else {
  best_k_manual <- as.integer(opt$k)
  if (!best_k_manual %in% k_range) stop(paste("k must be within", min(k_range), "-", max(k_range)))
  sub_df <- sil_df[sil_df$k == best_k_manual, ]
  best_row <- sub_df[which.max(sub_df$avg_sil), ]
}
best_method <- best_row$method
best_k <- best_row$k
cat(paste0("Best method: ", best_method, " | k: ", best_k,
           " | Silhouette: ", round(best_row$avg_sil, 3), "\n"))

svg(file.path(outdir_k, "Scree_plot.svg"), width = 9, height = 5)
print(
  fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 20),
           linecolor = "#2E9FDF", barfill = "steelblue",
           xlab = "Principal component",
           ylab = "Percentage of variance explained") +
    theme_minimal(base_size = 12)
)
dev.off()

svg(file.path(outdir_k, "Silhouette_k_selection.svg"), width = 8, height = 5)
print(
  ggplot(sil_df, aes(x = k, y = avg_sil, color = method)) +
    geom_line() + geom_point(size = 2) +
    scale_x_continuous(breaks = k_range) +
    scale_color_jco() +
    theme_minimal(base_size = 12) +
    labs(title = "Optimal k selection by silhouette",
         x = "Number of clusters (k)", y = "Average silhouette width", color = "Method")
)
dev.off()


# ============================================================================ #
#                    Clustering and silhouette analyses                        #
# ============================================================================ #

make_clusters <- function(pca_scores, dist, method = c("kmeans", "hclust", "pam"),
                          k, hclust_m = "ward.D2", seed = 42) {
  method <- match.arg(method)
  clusters <- switch(method,
    kmeans = { set.seed(seed); kmeans(pca_scores, centers = k, nstart = 25)$cluster },
    hclust = cutree(hclust(dist, method = hclust_m), k = k),
    pam    = pam(dist, k = k)$clustering
  )
  clusters <- as.integer(clusters)

  sil <- silhouette(clusters, dist)
  row.names(sil) <- rownames(as.matrix(dist))
  sil_avg <- round(summary(sil)$avg.width, 3)

  print(
    fviz_silhouette(sil, label = TRUE, print.summary = FALSE) +
      scale_fill_jco() + scale_color_jco() +
      theme_minimal(base_size = 12) +
      labs(title = paste0("Method: ", method, " | k = ", k, " | Avg silhouette: ", sil_avg)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  )
  clusters
}

svg(file.path(outdir_k, "Silhouette_kmeans.svg"), width = 9, height = 5)
pca_df$Cluster_km <- make_clusters(pca_scores, dist_pca, "kmeans", best_k)
dev.off()

svg(file.path(outdir_k, "Silhouette_hclust.svg"), width = 9, height = 5)
pca_df$Cluster_hc <- make_clusters(pca_scores, dist_pca, "hclust", best_k)
dev.off()

svg(file.path(outdir_k, "Silhouette_PAM.svg"), width = 9, height = 5)
pca_df$Cluster_pam <- make_clusters(pca_scores, dist_pca, "pam", best_k)
dev.off()

# Primary cluster assignment (categorical)
pca_df$Cluster_best <- factor(pca_df[[paste0("Cluster_",
  switch(best_method, kmeans = "km", hclust = "hc", pam = "pam"))]])

write_tsv(pca_df %>% dplyr::select(all_of(opt$samples_col), starts_with("Cluster_")),
          file.path(outdir_k, "Cluster_assignments.tsv"))


# ============================================================================ #
#                       Statistical validation (PERMANOVA)                     #
# ============================================================================ #

# Cluster_best must be a factor: numeric labels would be treated as a covariate
permanova_res <- adonis2(data_pca_final ~ Cluster_best, data = pca_df,
                         permutations = 999, method = "euclidean")
writeLines(capture.output(print(permanova_res)), con = file.path(outdir_k, "adonis.txt"))

pairwise_res <- pairwise.adonis(data_pca_final, factors = pca_df$Cluster_best,
                                sim.method = "euclidean", p.adjust.m = "BH")
write.table(pairwise_res, file.path(outdir_k, "pairwiseAdonis.tsv"),
            sep = "\t", row.names = FALSE)


# ============================================================================ #
#                              PCA Visualization                               #
# ============================================================================ #

pc_var <- round(100 * summary(pca_res)$importance[2, 1:3], 1)

svg(file.path(outdir_k, "PCA_2D_plot.svg"),
    width = 6 * (1 + pc_var[1] / 100), height = 5 * (1 + pc_var[2] / 100))
print(
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster_best)) +
    geom_point(size = 3) +
    scale_color_jco() +
    theme_minimal() +
    labs(caption = paste0("*Method: ", best_method, " | k = ", best_k),
         color = "Cluster",
         x = paste0("PC1 (", pc_var[1], "%)"),
         y = paste0("PC2 (", pc_var[2], "%)")) +
    theme(plot.caption = element_text(hjust = 0))
)
dev.off()

p <- plot_ly(pca_df,
             x = ~PC1, y = ~PC2, z = ~PC3,
             color = ~Cluster_best,
             colors = pal_jco()(10)[seq_len(nlevels(pca_df$Cluster_best))],
             text = ~paste("Sample:", .data[[opt$samples_col]]),
             type = "scatter3d", mode = "markers",
             marker = list(size = 5, opacity = 0.8)) %>%
  layout(
    legend = list(title = list(text = "<b>Cluster</b>", font = list(size = 18)),
                  font = list(size = 16)),
    scene = list(
      xaxis = list(title = paste0("PC1 (", pc_var[1], "%)")),
      yaxis = list(title = paste0("PC2 (", pc_var[2], "%)")),
      zaxis = list(title = paste0("PC3 (", pc_var[3], "%)"))
    )
  )
saveWidget(p, file.path(outdir_k, "PCA_3D_plot.html"), selfcontained = TRUE)
unlink(file.path(outdir_k, "PCA_3D_plot_files"), recursive = TRUE, force = TRUE)

cat("Clustering and PCA analyses completed successfully. Results saved in:", outdir_k, "\n")
