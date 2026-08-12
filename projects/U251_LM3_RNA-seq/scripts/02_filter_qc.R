# ============================================================
# 02_filter_qc.R
# Purpose:
#   Compare two low-expression filtering thresholds:
#   1) expression >= 0.5 in at least 3 samples
#   2) expression >= 1   in at least 3 samples
#
# QC:
#   - Boxplot
#   - PCA
# ============================================================


# =========================
# 1. Load packages
# =========================

library(readxl)
library(ggplot2)


# =========================
# 2. Input
# =========================

input_file <- "datasets/raw/基因表达水平检测结果-U251.xlsx"

sample_cols <- c(
  "OE1", "OE2", "OE3",
  "NC1", "NC2", "NC3"
)

group <- factor(
  c("OE", "OE", "OE", "NC", "NC", "NC"),
  levels = c("NC", "OE")
)


# =========================
# 3. Read expression data
# =========================

dat <- read_excel(input_file)

expr <- as.matrix(dat[, sample_cols])

storage.mode(expr) <- "numeric"

rownames(expr) <- dat$gene_id


# =========================
# 4. log2 transformation
# =========================

expr_log2 <- log2(expr + 1)


# =========================
# 5. Filtering
# =========================

# Strategy 1:
# expression >= 0.5 in at least 3 samples

keep_05 <- rowSums(expr >= 0.5) >= 3


# Strategy 2:
# expression >= 1 in at least 3 samples

keep_1 <- rowSums(expr >= 1) >= 3


expr_log2_05 <- expr_log2[
  keep_05,
  ,
  drop = FALSE
]

expr_log2_1 <- expr_log2[
  keep_1,
  ,
  drop = FALSE
]


# Check dimensions

cat("\nOriginal matrix:\n")
print(dim(expr_log2))

cat("\nFilter >= 0.5 in at least 3 samples:\n")
print(dim(expr_log2_05))

cat("\nFilter >= 1 in at least 3 samples:\n")
print(dim(expr_log2_1))


# =========================
# 6. Output directory
# =========================

dir.create(
  "results/QC",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 7. Boxplot
# ============================================================


# -------------------------
# 7.1 Filter >= 0.5
# -------------------------

png(
  "results/QC/boxplot_filter_0.5.png",
  width = 1200,
  height = 900,
  res = 150
)

boxplot(
  expr_log2_05,
  las = 2,
  main = "Expression >= 0.5 in at least 3 samples",
  ylab = "log2(expression + 1)"
)

dev.off()


# -------------------------
# 7.2 Filter >= 1
# -------------------------

png(
  "results/QC/boxplot_filter_1.png",
  width = 1200,
  height = 900,
  res = 150
)

boxplot(
  expr_log2_1,
  las = 2,
  main = "Expression >= 1 in at least 3 samples",
  ylab = "log2(expression + 1)"
)

dev.off()


# ============================================================
# 8. PCA: filter >= 0.5
# ============================================================

pca_05 <- prcomp(
  t(expr_log2_05),
  center = TRUE,
  scale. = FALSE
)


# Percentage of variance explained by each PC

variance_05 <- (
  pca_05$sdev^2 /
  sum(pca_05$sdev^2)
) * 100


# Build PCA dataframe

pca_df_05 <- data.frame(
  sample = colnames(expr_log2_05),
  group = group,
  PC1 = pca_05$x[, 1],
  PC2 = pca_05$x[, 2]
)


# Plot

p_05 <- ggplot(
  pca_df_05,
  aes(
    x = PC1,
    y = PC2,
    color = group,
    fill = group
  )
) +

  # Group ellipse
  stat_ellipse(
    geom = "polygon",
    level = 0.95,
    alpha = 0.15,
    linewidth = 0.8
  ) +

  # Sample points
  geom_point(
    size = 4
  ) +

  # Sample names
  geom_text(
    aes(label = sample),
    vjust = -1,
    size = 4,
    show.legend = FALSE
  ) +

  # Group colors
  scale_color_manual(
    values = c(
      "NC" = "#4C78A8",
      "OE" = "#E45756"
    )
  ) +

  scale_fill_manual(
    values = c(
      "NC" = "#4C78A8",
      "OE" = "#E45756"
    )
  ) +

  labs(
    title = "PCA - filter >= 0.5",
    x = paste0(
      "PC1 (",
      round(variance_05[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_05[2], 1),
      "%)"
    ),
    color = "Group",
    fill = "Group"
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "right"
  )


# Save PCA

ggsave(
  filename = "results/QC/PCA_filter_0.5.png",
  plot = p_05,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 9. PCA: filter >= 1
# ============================================================

pca_1 <- prcomp(
  t(expr_log2_1),
  center = TRUE,
  scale. = FALSE
)


# Percentage of variance explained

variance_1 <- (
  pca_1$sdev^2 /
  sum(pca_1$sdev^2)
) * 100


# Build PCA dataframe

pca_df_1 <- data.frame(
  sample = colnames(expr_log2_1),
  group = group,
  PC1 = pca_1$x[, 1],
  PC2 = pca_1$x[, 2]
)


# Plot

p_1 <- ggplot(
  pca_df_1,
  aes(
    x = PC1,
    y = PC2,
    color = group,
    fill = group
  )
) +

  stat_ellipse(
    geom = "polygon",
    level = 0.95,
    alpha = 0.15,
    linewidth = 0.8
  ) +

  geom_point(
    size = 4
  ) +

  geom_text(
    aes(label = sample),
    vjust = -1,
    size = 4,
    show.legend = FALSE
  ) +

  scale_color_manual(
    values = c(
      "NC" = "#4C78A8",
      "OE" = "#E45756"
    )
  ) +

  scale_fill_manual(
    values = c(
      "NC" = "#4C78A8",
      "OE" = "#E45756"
    )
  ) +

  labs(
    title = "PCA - filter >= 1",
    x = paste0(
      "PC1 (",
      round(variance_1[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_1[2], 1),
      "%)"
    ),
    color = "Group",
    fill = "Group"
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "right"
  )


# Save PCA

ggsave(
  filename = "results/QC/PCA_filter_1.png",
  plot = p_1,
  width = 8,
  height = 6,
  dpi = 300
)


# =========================
# 10. Finished
# =========================

cat("\nQC finished.\n")

cat(
  "\nPCA >= 0.5: PC1 =",
  round(variance_05[1], 1),
  "%, PC2 =",
  round(variance_05[2], 1),
  "%\n"
)

cat(
  "PCA >= 1: PC1 =",
  round(variance_1[1], 1),
  "%, PC2 =",
  round(variance_1[2], 1),
  "%\n"
)

library(pheatmap)

# ============================================================
# 10. Sample correlation heatmap
# ============================================================

annotation_col <- data.frame(
  Group = group
)

rownames(annotation_col) <- sample_cols

cor_05 <- cor(
  expr_log2_05,
  method = "pearson"
)

print(round(cor_05, 3))

pheatmap(
  cor_05,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_col,
  annotation_col = annotation_col,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Sample correlation - filter >= 0.5",
  filename = "results/QC/correlation_filter_0.5.png",
  width = 7,
  height = 6
)

cor_1 <- cor(
  expr_log2_1,
  method = "pearson"
)

print(round(cor_1, 3))

pheatmap(
  cor_1,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_col,
  annotation_col = annotation_col,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Sample correlation - filter >= 1",
  filename = "results/QC/correlation_filter_1.png",
  width = 7,
  height = 6
)

write.table(
  cor_05,
  "results/QC/correlation_filter_0.5.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  cor_1,
  "results/QC/correlation_filter_1.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)