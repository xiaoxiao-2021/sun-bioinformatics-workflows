args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Expected one config.yml argument")
cfg <- yaml::read_yaml(normalizePath(args[1], mustWork = TRUE))
project_dir <- normalizePath(cfg$project_dir, mustWork = TRUE)
metadata <- read.delim(file.path(project_dir, cfg$sample_metadata), check.names = FALSE, stringsAsFactors = FALSE)
dat <- readxl::read_excel(file.path(project_dir, cfg$input_file), sheet = 1)

# Build expression matrix in metadata order
sample_cols <- metadata$sample
expr <- as.matrix(dat[, sample_cols, drop = FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- as.character(dat[[cfg$gene_id_col]])
expr_log2 <- log2(expr + 1)
group <- factor(metadata$group, levels = c(cfg$control_group, cfg$case_group))
qc_dir <- file.path(project_dir, "results", cfg$dataset_id, "QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# Compare the configured filter with a two-fold stricter expression threshold
qc_cutoffs <- unique(c(cfg$expression_cutoff, cfg$expression_cutoff * 2))
for (cutoff in qc_cutoffs) {
  keep <- rowSums(expr >= cutoff) >= cfg$min_samples
  filtered <- expr_log2[keep, , drop = FALSE]
  tag <- gsub("\\.", "p", format(cutoff, trim = TRUE, scientific = FALSE))
  if (nrow(filtered) < 2) stop("Too few genes remain for QC at cutoff ", cutoff)
  cat("QC cutoff", cutoff, ":", nrow(filtered), "genes\n")

  png(file.path(qc_dir, paste0("boxplot_filter_", tag, ".png")), width = 1400, height = 900, res = 150)
  boxplot(filtered, las = 2, col = "#8FBBD9", ylab = "log2(expression + 1)", main = paste("Filter >=", cutoff, "in", cfg$min_samples, "samples"))
  dev.off()

  # PCA colored from metadata groups
  pca <- prcomp(t(filtered), scale. = FALSE)
  variance <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  pca_df <- data.frame(sample = sample_cols, group = group, PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  p <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, color = group, label = sample)) +
    ggplot2::geom_point(size = 4) + ggplot2::geom_text(vjust = -0.8, show.legend = FALSE) +
    ggplot2::labs(x = paste0("PC1 (", variance[1], "%)"), y = paste0("PC2 (", variance[2], "%)"), color = "Group") +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(qc_dir, paste0("PCA_filter_", tag, ".png")), p, width = 7, height = 5, dpi = 150)

  # Sample correlation
  cor_mat <- cor(filtered, method = "pearson")
  write.table(cor_mat, file.path(qc_dir, paste0("correlation_filter_", tag, ".tsv")), sep = "\t", quote = FALSE, col.names = NA)
  annotation_col <- data.frame(Group = group, row.names = sample_cols)
  pheatmap::pheatmap(cor_mat, annotation_col = annotation_col, filename = file.path(qc_dir, paste0("correlation_filter_", tag, ".png")), width = 7, height = 6)
}
