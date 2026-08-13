args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Expected one config.yml argument")
cfg <- yaml::read_yaml(normalizePath(args[1], mustWork = TRUE))
project_dir <- normalizePath(cfg$project_dir, mustWork = TRUE)
metadata <- read.delim(file.path(project_dir, cfg$sample_metadata), check.names = FALSE, stringsAsFactors = FALSE)
de_dir <- file.path(project_dir, "results", cfg$dataset_id, "DE")
figure_dir <- file.path(project_dir, "figures", cfg$dataset_id, "DE")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- paste(cfg$dataset_id, cfg$case_group, "vs", cfg$control_group, sep = "_")
lfc_tag <- format(cfg$TREAT_lfc_cutoff, trim = TRUE, scientific = FALSE)
fdr_tag <- format(cfg$TREAT_FDR_cutoff, trim = TRUE, scientific = FALSE)
deg_fdr_tag <- format(cfg$deg_FDR_cutoff, trim = TRUE, scientific = FALSE)
deg_lfc_tag <- format(cfg$deg_logFC_cutoff, trim = TRUE, scientific = FALSE)
pvalue_tag <- format(cfg$pvalue_cutoff, trim = TRUE, scientific = FALSE)

# Volcano plots use the same standard limma all-genes result
volcano_all <- read.delim(file.path(de_dir, paste0(prefix, "_limma_all_genes_annotated.tsv")), check.names = FALSE)

# Save separate FDR and nominal P-value volcano plots
save_volcano <- function(p_column, cutoff, statistic_label, file_tag) {
  volcano_data <- volcano_all[
    !is.na(volcano_all$logFC) & !is.na(volcano_all[[p_column]]),
    , drop = FALSE
  ]
  volcano_data$change <- "NS"
  volcano_data$change[
    volcano_data[[p_column]] < cutoff &
      volcano_data$logFC >= cfg$volcano_logFC_cutoff
  ] <- "Up"
  volcano_data$change[
    volcano_data[[p_column]] < cutoff &
      volcano_data$logFC <= -cfg$volcano_logFC_cutoff
  ] <- "Down"
  volcano_plot_data <- volcano_data[volcano_data[[p_column]] < 1, , drop = FALSE]
  volcano_plot_data$minus_log10_significance <- -log10(
    pmax(volcano_plot_data[[p_column]], .Machine$double.xmin)
  )

  label_available <- !is.na(volcano_plot_data$SYMBOL) & volcano_plot_data$SYMBOL != ""
  label_up <- volcano_plot_data[volcano_plot_data$change == "Up" & label_available, , drop = FALSE]
  label_down <- volcano_plot_data[volcano_plot_data$change == "Down" & label_available, , drop = FALSE]
  label_up <- head(label_up[order(label_up[[p_column]]), , drop = FALSE], cfg$volcano_label_n)
  label_down <- head(label_down[order(label_down[[p_column]]), , drop = FALSE], cfg$volcano_label_n)
  label_genes <- rbind(label_up, label_down)

  volcano <- ggplot2::ggplot(
    volcano_plot_data,
    ggplot2::aes(logFC, minus_log10_significance, color = change)
  ) +
    ggplot2::geom_point(alpha = .55, size = 1.2) +
    ggplot2::geom_vline(xintercept = c(-cfg$volcano_logFC_cutoff, cfg$volcano_logFC_cutoff), linetype = 2, color = "grey50") +
    ggplot2::geom_hline(yintercept = -log10(cutoff), linetype = 2, color = "grey50") +
    ggrepel::geom_text_repel(data = label_genes, ggplot2::aes(label = SYMBOL), size = 3, max.overlaps = Inf) +
    ggplot2::scale_color_manual(values = c(Down = "#4DBBD5", NS = "#BDBDBD", Up = "#E64B35")) +
    ggplot2::labs(
      title = paste(cfg$dataset_id, cfg$case_group, "vs", cfg$control_group, "-", statistic_label),
      subtitle = paste0(statistic_label, " < ", cutoff, "; |logFC| >= ", cfg$volcano_logFC_cutoff),
      x = "log2 fold change", y = paste0("-log10(", statistic_label, ")"), color = NULL
    ) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(figure_dir, paste0(prefix, "_volcano_", file_tag, ".png")), volcano, width = 7, height = 6, dpi = 300)
  ggplot2::ggsave(file.path(figure_dir, paste0(prefix, "_volcano_", file_tag, ".pdf")), volcano, width = 7, height = 6)
  cat(statistic_label, "volcano classes:\n"); print(table(volcano_data$change))
}

unlink(file.path(figure_dir, paste0(prefix, "_volcano.", c("png", "pdf"))))
save_volcano("adj.P.Val", cfg$volcano_FDR_cutoff, "FDR", paste0("FDR", format(cfg$volcano_FDR_cutoff, trim = TRUE)))
save_volcano("P.Value", cfg$pvalue_cutoff, "P.Value", paste0("P", pvalue_tag))

# Select downstream DEG source
if (cfg$downstream_deg_method == "limma") {
  deg_file <- file.path(
    de_dir,
    paste0(prefix, "_limma_FDR", deg_fdr_tag, "_logFC", deg_lfc_tag, "_DEG_annotated.tsv")
  )
} else if (cfg$downstream_deg_method == "limma_pvalue") {
  deg_file <- file.path(
    de_dir,
    paste0(prefix, "_limma_P", pvalue_tag, "_logFC", deg_lfc_tag, "_DEG_annotated.tsv")
  )
} else {
  deg_file <- file.path(
    de_dir,
    paste0(prefix, "_TREAT_lfc", lfc_tag, "_FDR", fdr_tag, "_DEG_annotated.tsv")
  )
}
cat("Heatmap DEG source:", cfg$downstream_deg_method, "\n")
deg <- read.delim(deg_file, check.names = FALSE)
existing_heatmaps <- list.files(figure_dir, full.names = TRUE)
existing_heatmaps <- existing_heatmaps[
  startsWith(basename(existing_heatmaps), paste0(prefix, "_top")) &
    grepl("_DEG_heatmap\\.(png|pdf)$", basename(existing_heatmaps))
]
unlink(existing_heatmaps)

# Skip heatmap if too few genes
if (!nrow(deg)) {
  cat("No DEG available for heatmap. Heatmap skipped.\n")
} else if (nrow(deg) == 1) {
  cat("Only one DEG available. Clustered heatmap skipped.\n")
} else {
  expr <- readRDS(file.path(project_dir, "datasets", "processed", paste0(cfg$dataset_id, "_bulk_workflow_expression_log2_filtered.rds")))
  deg <- deg[deg$ENSEMBL %in% rownames(expr), , drop = FALSE]
  deg <- head(deg[order(abs(deg$logFC), decreasing = TRUE), , drop = FALSE], cfg$heatmap_top_n)

  if (nrow(deg) < 2) {
    cat("Fewer than two DEG remain after expression matching. Clustered heatmap skipped.\n")
  } else {
    heat_expr <- expr[deg$ENSEMBL, metadata$sample, drop = FALSE]
    labels <- ifelse(is.na(deg$SYMBOL) | deg$SYMBOL == "", deg$ENSEMBL, deg$SYMBOL)
    rownames(heat_expr) <- make.unique(labels)
    heat_expr <- heat_expr[apply(heat_expr, 1, sd) > 0, , drop = FALSE]

    if (nrow(heat_expr) < 2) {
      cat("Fewer than two DEG remain after SD filtering. Clustered heatmap skipped.\n")
    } else {
      # Row Z-score for heatmap and metadata-based sample annotation
      heat_z <- t(scale(t(heat_expr)))
      annotation_col <- data.frame(Group = factor(metadata$group, levels = c(cfg$control_group, cfg$case_group)), row.names = metadata$sample)
      annotation_colors <- list(Group = setNames(c("#4DBBD5", "#F28E85"), c(cfg$control_group, cfg$case_group)))
      heatmap_base <- file.path(figure_dir, paste0(prefix, "_top", nrow(heat_z), "_DEG_heatmap"))
      pheatmap::pheatmap(heat_z, annotation_col = annotation_col, annotation_colors = annotation_colors, cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = TRUE, border_color = NA, filename = paste0(heatmap_base, ".png"), width = 8, height = max(6, 2 + .18 * nrow(heat_z)))
      pheatmap::pheatmap(heat_z, annotation_col = annotation_col, annotation_colors = annotation_colors, cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = TRUE, border_color = NA, filename = paste0(heatmap_base, ".pdf"), width = 8, height = max(6, 2 + .18 * nrow(heat_z)))
      cat("Heatmap genes:", nrow(heat_z), "\n")
    }
  }
}
