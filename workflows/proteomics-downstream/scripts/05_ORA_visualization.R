args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("05_ORA_visualization.R requires one config path")
for (package in c("yaml", "ggplot2")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing R package: ", package)
}

config <- yaml::read_yaml(normalizePath(args[[1]], mustWork = TRUE))
project_dir <- normalizePath(config$project_dir, mustWork = TRUE)
result_dir <- file.path(project_dir, "results", "ORA")
figure_dir <- file.path(project_dir, "figures", "ORA")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

for (evidence in c("formal", "exploratory")) {
  metric <- if (evidence == "formal") "p.adjust" else "pvalue"
  cutoff <- if (evidence == "formal") config$enrichment_qvalue_cutoff else config$enrichment_pvalue_cutoff
  evidence_title <- if (evidence == "formal") {
    paste0("formal adjusted evidence (BH FDR < ", cutoff, ")")
  } else {
    paste0("exploratory nominal evidence (P < ", cutoff, ")")
  }
  for (database in c("GO_BP", "KEGG")) {
    for (direction in c("UP", "DOWN")) {
      stem <- paste(config$dataset_id, "ORA", evidence, database, direction, sep = "_")
      input_file <- file.path(result_dir, paste0(stem, "_all.tsv"))
      if (!file.exists(input_file)) {
        cat("ORA plot skipped; table absent:", input_file, "\n")
        next
      }
      result <- utils::read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
      if (!nrow(result) || !metric %in% names(result)) {
        cat("ORA plot skipped; no pathways:", stem, "\n")
        next
      }
      result[[metric]] <- suppressWarnings(as.numeric(result[[metric]]))
      result$Count <- suppressWarnings(as.numeric(result$Count))
      result <- result[is.finite(result[[metric]]) & result[[metric]] < cutoff, , drop = FALSE]
      if (!nrow(result)) {
        cat("ORA plot skipped; no", evidence_title, "pathways:", database, direction, "\n")
        next
      }
      result <- head(result[order(result[[metric]], -result$Count), , drop = FALSE], as.integer(config$show_category_n))
      result$plot_description <- factor(make.unique(result$Description), levels = rev(make.unique(result$Description)))
      result$minus_log10 <- -log10(pmax(result[[metric]], .Machine$double.xmin))
      title <- paste(config$dataset_id, database, direction, evidence_title, sep = " | ")

      dot <- ggplot2::ggplot(result, ggplot2::aes(x = Count, y = plot_description, size = Count, color = minus_log10)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_gradient(low = "#4DBBD5", high = "#E64B35") +
        ggplot2::labs(title = title, x = "Foreground Entrez genes", y = NULL, color = paste0("-log10(", metric, ")"), size = "Count") +
        ggplot2::theme_bw(base_size = 10) + ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
      bar <- ggplot2::ggplot(result, ggplot2::aes(x = minus_log10, y = plot_description, fill = minus_log10)) +
        ggplot2::geom_col() + ggplot2::scale_fill_gradient(low = "#4DBBD5", high = "#E64B35") +
        ggplot2::labs(title = title, x = paste0("-log10(", metric, ")"), y = NULL, fill = NULL) +
        ggplot2::theme_bw(base_size = 10) + ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
      for (plot_type in c("dotplot", "barplot")) {
        plot <- if (plot_type == "dotplot") dot else bar
        output_stem <- file.path(figure_dir, paste0(stem, "_", plot_type))
        ggplot2::ggsave(paste0(output_stem, ".png"), plot, width = 9, height = 6, dpi = 300)
        ggplot2::ggsave(paste0(output_stem, ".pdf"), plot, width = 9, height = 6)
      }
    }
  }
}
