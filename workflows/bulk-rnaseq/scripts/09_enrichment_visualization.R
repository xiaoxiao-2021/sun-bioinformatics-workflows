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
enrichment_p_tag <- format(cfg$enrichment_pvalue_cutoff, trim = TRUE, scientific = FALSE)
enrichment_fdr_tag <- format(cfg$enrichment_FDR_cutoff, trim = TRUE, scientific = FALSE)
enrichment_gene_logFC_cutoff <- cfg$enrichment_gene_logFC_cutoff
if (is.null(enrichment_gene_logFC_cutoff)) {
  enrichment_gene_logFC_cutoff <- cfg$deg_logFC_cutoff
}
if (is.null(enrichment_gene_logFC_cutoff)) {
  stop("Missing both enrichment_gene_logFC_cutoff and fallback deg_logFC_cutoff")
}
enrichment_gene_lfc_tag <- format(
  enrichment_gene_logFC_cutoff, trim = TRUE, scientific = FALSE
)
pvalue_tag <- format(cfg$pvalue_cutoff, trim = TRUE, scientific = FALSE)
deg_fdr_tag <- format(cfg$deg_FDR_cutoff, trim = TRUE, scientific = FALSE)
lfc_tag <- format(cfg$TREAT_lfc_cutoff, trim = TRUE, scientific = FALSE)
fdr_tag <- format(cfg$TREAT_FDR_cutoff, trim = TRUE, scientific = FALSE)

if (cfg$downstream_deg_method == "limma") {
  ora_base <- paste0(
    prefix, "_ORA_geneFDR", deg_fdr_tag,
    "_logFC", enrichment_gene_lfc_tag
  )
  foreground_label <- paste0(
    "Gene FDR < ", cfg$deg_FDR_cutoff,
    "; |logFC| >= ", enrichment_gene_logFC_cutoff
  )
} else if (cfg$downstream_deg_method == "limma_pvalue") {
  ora_base <- paste0(
    prefix, "_ORA_geneP", pvalue_tag,
    "_logFC", enrichment_gene_lfc_tag
  )
  foreground_label <- paste0(
    "Gene P < ", cfg$pvalue_cutoff,
    "; |logFC| >= ", enrichment_gene_logFC_cutoff
  )
} else {
  ora_base <- paste0(
    prefix, "_ORA_TREAT_lfc", lfc_tag, "_FDR", fdr_tag
  )
  foreground_label <- paste0(
    "TREAT |logFC| > ", cfg$TREAT_lfc_cutoff,
    "; FDR < ", cfg$TREAT_FDR_cutoff
  )
}

# Draw either formal FDR or exploratory nominal P enrichment
make_dotplot <- function(database, direction, result_type) {
  analysis_name <- paste(database, direction, sep = "_")
  analysis_label <- if (database == "GO_BP") "GO Biological Process" else "KEGG Pathways"
  direction_label <- if (direction == "UP") "Upregulated Genes" else "Downregulated Genes"
  output_base <- paste(ora_base, analysis_name, sep = "_")

  if (result_type == "FDR") {
    result_tag <- paste0("FDR", enrichment_fdr_tag)
    metric_column <- "p.adjust"
    metric_label <- "-log10(FDR)"
    cutoff <- cfg$enrichment_FDR_cutoff
    subtitle <- paste0("FDR < ", cutoff, "\n", foreground_label)
    missing_message <- "No FDR enrichment result for"
  } else {
    result_tag <- paste0("P", enrichment_p_tag)
    metric_column <- "pvalue"
    metric_label <- "-log10(P.Value)"
    cutoff <- cfg$enrichment_pvalue_cutoff
    subtitle <- paste0(
      "Exploratory: nominal P < ", cutoff, "\n", foreground_label
    )
    missing_message <- "No nominal enrichment result for"
  }

  input_file <- file.path(result_dir, paste0(output_base, "_", result_tag, ".tsv"))
  png_file <- file.path(figure_dir, paste0(output_base, "_", result_tag, "_dotplot.png"))
  pdf_file <- file.path(figure_dir, paste0(output_base, "_", result_tag, "_dotplot.pdf"))
  plot_data_file <- file.path(plot_data_dir, paste0(output_base, "_", result_tag, "_plot_terms.tsv"))

  if (!file.exists(input_file) || file.info(input_file)$size == 0) {
    unlink(c(png_file, pdf_file, plot_data_file))
    cat(missing_message, analysis_name, ". Plot skipped.\n")
    return(invisible(NULL))
  }
  enrichment <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(enrichment)) {
    unlink(c(png_file, pdf_file, plot_data_file))
    cat(missing_message, analysis_name, ". Plot skipped.\n")
    return(invisible(NULL))
  }

  required <- c("Description", "GeneRatio", "Count", metric_column)
  if (!all(required %in% names(enrichment))) {
    stop("Missing enrichment plot column(s): ", paste(setdiff(required, names(enrichment)), collapse = ", "))
  }

  ratio_parts <- strsplit(as.character(enrichment$GeneRatio), "/", fixed = TRUE)
  enrichment$GeneRatio_numeric <- vapply(
    ratio_parts,
    function(x) as.numeric(x[1]) / as.numeric(x[2]),
    numeric(1)
  )
  enrichment$plot_metric <- as.numeric(enrichment[[metric_column]])
  enrichment$minus_log10_metric <- -log10(
    pmax(enrichment$plot_metric, .Machine$double.xmin)
  )
  enrichment <- enrichment[
    is.finite(enrichment$GeneRatio_numeric) &
      is.finite(enrichment$plot_metric) & enrichment$plot_metric < cutoff,
    , drop = FALSE
  ]
  if (!nrow(enrichment)) {
    unlink(c(png_file, pdf_file, plot_data_file))
    cat(missing_message, analysis_name, ". Plot skipped.\n")
    return(invisible(NULL))
  }

  # Show at most configured N, ordered by the selected pathway statistic
  plot_data <- head(enrichment[order(enrichment$plot_metric), , drop = FALSE], cfg$show_category_n)
  write.table(plot_data, plot_data_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(GeneRatio_numeric, Description, size = Count, color = minus_log10_metric)
  ) +
    ggplot2::geom_point(alpha = .85) +
    ggplot2::scale_color_gradient(low = "#4DBBD5", high = "#E64B35") +
    ggplot2::labs(
      title = paste(analysis_label, "-", direction_label),
      subtitle = subtitle,
      x = "GeneRatio", y = NULL, color = metric_label
    ) +
    ggplot2::theme_classic(base_size = 12)
  plot_height <- max(3.5, 2.5 + .30 * nrow(plot_data))
  ggplot2::ggsave(png_file, p, width = 9, height = plot_height, dpi = 300)
  ggplot2::ggsave(pdf_file, p, width = 9, height = plot_height)
  cat(analysis_name, result_tag, ": plotted", nrow(plot_data), "terms\n")
}

# Remove legacy untagged dotplots from earlier workflow versions
legacy_dotplots <- list.files(
  figure_dir,
  pattern = paste0("^", prefix, "_(GO_BP|KEGG)_(UP|DOWN)_dotplot\\.(png|pdf)$"),
  full.names = TRUE
)
unlink(legacy_dotplots)
legacy_plot_data <- list.files(
  plot_data_dir,
  pattern = paste0("^", prefix, "_(GO_BP|KEGG)_(UP|DOWN)_plot_terms\\.tsv$"),
  full.names = TRUE
)
unlink(legacy_plot_data)

for (database in c("GO_BP", "KEGG")) {
  for (direction in c("UP", "DOWN")) {
    # Plot FDR enrichment
    make_dotplot(database, direction, "FDR")
    # Plot nominal P enrichment
    make_dotplot(database, direction, "P")
  }
}
