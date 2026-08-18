args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("03_visualization.R requires one config path")
for (package in c("yaml", "ggplot2", "pheatmap")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing R package: ", package)
}

config <- yaml::read_yaml(normalizePath(args[[1]], mustWork = TRUE))
project_dir <- normalizePath(config$project_dir, mustWork = TRUE)
processed_dir <- file.path(project_dir, "datasets", "processed")
figure_dir <- file.path(project_dir, "figures", "DE")
result_dir <- file.path(project_dir, "results", "DE")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
data <- utils::read.delim(
  file.path(processed_dir, paste0(config$dataset_id, "_all_proteins.tsv")),
  check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA")
)
metadata <- NULL
if (!is.null(config$sample_metadata)) {
  metadata <- utils::read.delim(file.path(project_dir, config$sample_metadata), check.names = FALSE, stringsAsFactors = FALSE)
}
flags <- readRDS(file.path(processed_dir, paste0(config$dataset_id, "_QC_flags.rds")))

# Encode every result-changing threshold in filenames so parameter sweeps do not
# overwrite previous plots or their audit tables (for example 0.05 -> 0p05).
cutoff_tag <- function(value) {
  text <- format(value, scientific = FALSE, trim = TRUE, digits = 15)
  if (grepl("\\.", text)) text <- sub("0+$", "", text)
  text <- sub("\\.$", "", text)
  gsub("\\.", "p", text)
}
formal_de_tag <- paste0("q", cutoff_tag(config$deg_qvalue_cutoff), "_logFC", cutoff_tag(config$deg_logFC_cutoff))
exploratory_de_tag <- paste0("p", cutoff_tag(config$deg_pvalue_cutoff), "_logFC", cutoff_tag(config$deg_logFC_cutoff))
volcano_tag <- paste0(
  "q", cutoff_tag(config$volcano_qvalue_cutoff),
  "_logFC", cutoff_tag(config$volcano_logFC_cutoff),
  "_label", as.integer(config$volcano_label_n)
)
heatmap_tag <- paste0(formal_de_tag, "_top", as.integer(config$heatmap_top_n))

data$logFC <- suppressWarnings(as.numeric(data$logFC))
data$P.Value <- suppressWarnings(as.numeric(data$P.Value))
data$Q.Value <- suppressWarnings(as.numeric(data$Q.Value))
formal_de <- data[
  is.finite(data$Q.Value) & data$Q.Value < config$deg_qvalue_cutoff &
    is.finite(data$logFC) & abs(data$logFC) >= config$deg_logFC_cutoff,
  , drop = FALSE
]
exploratory_de <- data[
  is.finite(data$P.Value) & data$P.Value < config$deg_pvalue_cutoff &
    is.finite(data$logFC) & abs(data$logFC) >= config$deg_logFC_cutoff,
  , drop = FALSE
]
formal_de$Direction <- ifelse(formal_de$logFC > 0, "Up", "Down")
exploratory_de$Direction <- ifelse(exploratory_de$logFC > 0, "Up", "Down")
utils::write.table(
  formal_de, file.path(result_dir, paste0(config$dataset_id, "_DE_formal_", formal_de_tag, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
utils::write.table(
  exploratory_de, file.path(result_dir, paste0(config$dataset_id, "_DE_exploratory_", exploratory_de_tag, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
positive_q <- data$Q.Value[is.finite(data$Q.Value) & data$Q.Value > 0]
q_floor <- if (length(positive_q)) min(positive_q) / 10 else .Machine$double.xmin
data$plot_q <- ifelse(is.finite(data$Q.Value), pmax(data$Q.Value, q_floor), NA_real_)
data$significance <- "Not significant"
formal <- is.finite(data$Q.Value) & data$Q.Value < config$volcano_qvalue_cutoff &
  is.finite(data$logFC) & abs(data$logFC) >= config$volcano_logFC_cutoff
data$significance[formal & data$logFC > 0] <- "Up"
data$significance[formal & data$logFC < 0] <- "Down"
data$significance <- factor(data$significance, levels = c("Down", "Not significant", "Up"))
data$plot_label <- ifelse(!is.na(data$SYMBOL) & nzchar(data$SYMBOL), data$SYMBOL, data$PROTEIN_ID)

pick_labels <- function(direction) {
  candidates <- which(data$significance == direction & !is.na(data$plot_label) & nzchar(data$plot_label))
  if (!length(candidates)) return(integer())
  order_values <- order(data$Q.Value[candidates], -abs(data$logFC[candidates]), na.last = TRUE)
  head(candidates[order_values], as.integer(config$volcano_label_n))
}
label_rows <- c(pick_labels("Up"), pick_labels("Down"))
label_table <- data[label_rows, c("PROTEIN_ID", "SYMBOL", "DESCRIPTION", "logFC", "P.Value", "Q.Value", "significance"), drop = FALSE]
utils::write.table(
  label_table,
  file.path(result_dir, paste0(config$dataset_id, "_volcano_labeled_", volcano_tag, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

volcano <- ggplot2::ggplot(data, ggplot2::aes(x = logFC, y = -log10(plot_q), color = significance)) +
  ggplot2::geom_point(alpha = 0.75, size = 1.7, na.rm = TRUE) +
  ggplot2::geom_vline(xintercept = c(-1, 1) * config$volcano_logFC_cutoff, linetype = "dashed", color = "grey45") +
  ggplot2::geom_hline(yintercept = -log10(config$volcano_qvalue_cutoff), linetype = "dashed", color = "grey45") +
  ggplot2::geom_text(
    data = data[label_rows, , drop = FALSE], ggplot2::aes(label = plot_label),
    color = "black", size = 3, check_overlap = TRUE, vjust = -0.6, show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = c(Down = "#377EB8", `Not significant` = "#BDBDBD", Up = "#E64B35")) +
  ggplot2::labs(
    title = paste0(config$dataset_id, ": vendor q-value volcano"),
    subtitle = sprintf("Formal: q < %g and |log2FC| >= %g", config$volcano_qvalue_cutoff, config$volcano_logFC_cutoff),
    x = paste0("Vendor log2 fold change (", config$case_group, " - ", config$control_group, ")"),
    y = expression(-log[10](vendor~q-value)), color = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = "top")
volcano_stem <- file.path(figure_dir, paste0(config$dataset_id, "_volcano_vendor_", volcano_tag))
ggplot2::ggsave(paste0(volcano_stem, ".png"), volcano, width = 8, height = 6, dpi = 300)
ggplot2::ggsave(paste0(volcano_stem, ".pdf"), volcano, width = 8, height = 6)

if (!isTRUE(flags$abundance_scale_ok)) {
  cat("Heatmap skipped: abundance scale did not pass the vendor-logFC agreement check.\n")
  quit(save = "no", status = 0L)
}

heatmap_rows <- which(
  is.finite(data$Q.Value) & data$Q.Value < config$deg_qvalue_cutoff &
    is.finite(data$logFC) & abs(data$logFC) >= config$deg_logFC_cutoff
)
heatmap_rows <- head(heatmap_rows[order(-abs(data$logFC[heatmap_rows]))], as.integer(config$heatmap_top_n))
if (length(heatmap_rows) < 2L) {
  cat("Heatmap skipped: fewer than two proteins pass the formal vendor thresholds.\n")
  quit(save = "no", status = 0L)
}

sample_ids <- metadata$sample_id
matrix_values <- as.matrix(data[heatmap_rows, sample_ids, drop = FALSE])
storage.mode(matrix_values) <- "numeric"
row_sd <- apply(matrix_values, 1L, stats::sd, na.rm = TRUE)
keep <- is.finite(row_sd) & row_sd > 0 & rowSums(is.finite(matrix_values)) >= 2L
matrix_values <- matrix_values[keep, , drop = FALSE]
heatmap_rows <- heatmap_rows[keep]
if (nrow(matrix_values) < 2L) {
  cat("Heatmap skipped: fewer than two selected proteins have variable finite abundance.\n")
  quit(save = "no", status = 0L)
}
z_matrix <- t(scale(t(matrix_values)))
labels <- ifelse(!is.na(data$SYMBOL[heatmap_rows]) & nzchar(data$SYMBOL[heatmap_rows]), data$SYMBOL[heatmap_rows], data$PROTEIN_ID[heatmap_rows])
rownames(z_matrix) <- make.unique(labels)
annotation <- data.frame(Group = metadata$group, row.names = metadata$sample_id, check.names = FALSE)
audit <- data[heatmap_rows, c("PROTEIN_ID", "SYMBOL", "DESCRIPTION", "logFC", "P.Value", "Q.Value"), drop = FALSE]
audit$heatmap_label <- rownames(z_matrix)
utils::write.table(
  audit, file.path(result_dir, paste0(config$dataset_id, "_actual_heatmap_proteins_", heatmap_tag, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

draw_heatmap <- function() {
  pheatmap::pheatmap(
    z_matrix, annotation_col = annotation, cluster_rows = TRUE, cluster_cols = TRUE,
    border_color = NA, fontsize_row = 8, main = paste0(config$dataset_id, " formal DE proteins (row Z-score)"),
    color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
    silent = TRUE
  )
}
heatmap_stem <- file.path(figure_dir, paste0(config$dataset_id, "_heatmap_formal_", heatmap_tag))
grDevices::png(paste0(heatmap_stem, ".png"), width = 8, height = 7, units = "in", res = 300)
grid::grid.newpage(); grid::grid.draw(draw_heatmap()$gtable); grDevices::dev.off()
grDevices::pdf(paste0(heatmap_stem, ".pdf"), width = 8, height = 7)
grid::grid.newpage(); grid::grid.draw(draw_heatmap()$gtable); grDevices::dev.off()
cat("Volcano and heatmap outputs written to:", figure_dir, "\n")
