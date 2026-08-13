args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Expected one config.yml argument")
cfg <- yaml::read_yaml(normalizePath(args[1], mustWork = TRUE))
project_dir <- normalizePath(cfg$project_dir, mustWork = TRUE)
result_dir <- file.path(project_dir, "results", cfg$dataset_id, "enrichment")
figure_dir <- file.path(project_dir, "figures", cfg$dataset_id, "enrichment")
plot_data_dir <- file.path(result_dir, "plot_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_data_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- paste(cfg$dataset_id, cfg$case_group, "vs", cfg$control_group, sep = "_")

# Draw horizontal dotplots for each significant ORA result
make_dotplot <- function(database, direction) {
  analysis_name <- paste(database, direction, sep = "_")
  output_base <- paste(prefix, analysis_name, sep = "_")
  input_file <- file.path(result_dir, paste0(output_base, "_FDR", cfg$enrichment_FDR_cutoff, ".tsv"))
  png_file <- file.path(figure_dir, paste0(output_base, "_dotplot.png"))
  pdf_file <- file.path(figure_dir, paste0(output_base, "_dotplot.pdf"))
  plot_data_file <- file.path(plot_data_dir, paste0(output_base, "_plot_terms.tsv"))
  if (!file.exists(input_file) || file.info(input_file)$size == 0) {
    unlink(c(png_file, pdf_file, plot_data_file))
    cat("No significant enrichment result for", analysis_name, ". Plot skipped.\n")
    return(invisible(NULL))
  }
  enrichment <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(enrichment)) {
    unlink(c(png_file, pdf_file, plot_data_file))
    cat("No significant enrichment result for", analysis_name, ". Plot skipped.\n")
    return(invisible(NULL))
  }
  required <- c("Description", "GeneRatio", "Count", "p.adjust")
  if (!all(required %in% names(enrichment))) {
    if (!nrow(enrichment)) {
      unlink(c(png_file, pdf_file, plot_data_file))
      cat("No significant enrichment result for", analysis_name, ". Plot skipped.\n")
      return(invisible(NULL))
    }
    stop("Missing enrichment plot column(s): ", paste(setdiff(required, names(enrichment)), collapse = ", "))
  }
  ratio_parts <- strsplit(enrichment$GeneRatio, "/", fixed = TRUE)
  enrichment$GeneRatio_numeric <- vapply(ratio_parts, function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
  enrichment$FDR <- as.numeric(enrichment$p.adjust)
  enrichment$minus_log10_FDR <- -log10(pmax(enrichment$FDR, .Machine$double.xmin))
  enrichment <- enrichment[is.finite(enrichment$GeneRatio_numeric) & is.finite(enrichment$FDR) & enrichment$FDR < cfg$enrichment_FDR_cutoff, ]
  if (!nrow(enrichment)) {
    unlink(c(png_file, pdf_file, plot_data_file))
    cat("No significant enrichment result for", analysis_name, ". Plot skipped.\n")
    return(invisible(NULL))
  }

  # Show at most configured N; when fewer are significant, show all
  plot_data <- head(enrichment[order(enrichment$FDR), ], cfg$show_category_n)
  write.table(plot_data, plot_data_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(GeneRatio_numeric, Description, size = Count, color = minus_log10_FDR)) +
    ggplot2::geom_point(alpha = .85) +
    ggplot2::scale_color_gradient(low = "#4DBBD5", high = "#E64B35") +
    ggplot2::labs(title = paste(gsub("_", " ", analysis_name), "ORA"), x = "GeneRatio", y = NULL, color = "-log10(FDR)") +
    ggplot2::theme_classic(base_size = 12)
  plot_height <- max(3.5, 2.5 + .30 * nrow(plot_data))
  ggplot2::ggsave(png_file, p, width = 9, height = plot_height, dpi = 300)
  ggplot2::ggsave(pdf_file, p, width = 9, height = plot_height)
  cat(analysis_name, ": plotted", nrow(plot_data), "terms\n")
}

make_dotplot("GO_BP", "UP")
make_dotplot("GO_BP", "DOWN")
make_dotplot("KEGG", "UP")
make_dotplot("KEGG", "DOWN")
