args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript workflows/bulk-rnaseq/run_pipeline.R path/to/config.yml")
}

config_file <- normalizePath(args[1], mustWork = TRUE)
config <- yaml::read_yaml(config_file)
required <- c(
  "dataset_id", "project_dir", "input_file", "sample_metadata", "gene_id_col",
  "case_group", "control_group", "expression_cutoff", "min_samples",
  "limma_FDR_cutoff", "deg_FDR_cutoff", "deg_logFC_cutoff", "pvalue_cutoff",
  "TREAT_lfc_cutoff", "TREAT_FDR_cutoff", "downstream_deg_method",
  "volcano_FDR_cutoff", "volcano_logFC_cutoff", "volcano_label_n",
  "heatmap_top_n", "heatmap_gene_filter", "enrichment_pvalue_cutoff",
  "enrichment_FDR_cutoff",
  "min_GS_size", "max_GS_size",
  "show_category_n", "gsea_pvalue_cutoff", "gsea_FDR_cutoff",
  "gsea_min_GS_size", "gsea_max_GS_size", "gsea_show_category_n",
  "gsea_curve_n"
)
missing_keys <- required[!vapply(required, function(x) !is.null(config[[x]]), logical(1))]
if (length(missing_keys)) stop("Missing config key(s): ", paste(missing_keys, collapse = ", "))
if (!grepl("^[A-Za-z0-9._-]+$", config$dataset_id)) stop("dataset_id contains unsafe characters")
if (config$case_group == config$control_group) stop("case_group and control_group must differ")
if (!config$downstream_deg_method %in% c("limma", "limma_pvalue", "treat")) {
  stop("downstream_deg_method must be 'limma', 'limma_pvalue', or 'treat'")
}
if (!config$heatmap_gene_filter %in% c("all", "annotated", "protein_coding")) {
  stop("heatmap_gene_filter must be 'all', 'annotated', or 'protein_coding'")
}
if (!is.null(config$volcano_label_filter) &&
    !config$volcano_label_filter %in% c("all", "annotated", "protein_coding")) {
  stop("volcano_label_filter must be 'all', 'annotated', or 'protein_coding'")
}
if (config$gsea_pvalue_cutoff <= 0 || config$gsea_pvalue_cutoff > 1 ||
    config$gsea_FDR_cutoff <= 0 || config$gsea_FDR_cutoff > 1) {
  stop("GSEA P-value and FDR cutoffs must be in (0, 1]")
}
if (config$gsea_min_GS_size < 1 ||
    config$gsea_max_GS_size < config$gsea_min_GS_size) {
  stop("GSEA gene-set size limits are invalid")
}
if (config$gsea_show_category_n < 1 || config$gsea_curve_n < 1) {
  stop("GSEA plot counts must be positive")
}

project_dir <- normalizePath(config$project_dir, mustWork = TRUE)
log_dir <- file.path(project_dir, "logs", config$dataset_id)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, "pipeline.log")
writeLines(paste("Pipeline started:", format(Sys.time())), log_file)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
workflow_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE))
scripts <- sprintf("%02d_%s.R", 1:11, c(
  "data_check", "filter_qc", "data_clean", "limma_DE", "gene_annotation",
  "add_annotation", "visualization", "enrichment", "enrichment_visualization",
  "GSEA", "GSEA_visualization"
))
labels <- c(
  "Data check", "Filter QC", "Data cleaning", "limma / TREAT",
  "Gene annotation", "Add annotation", "Visualization", "Enrichment",
  "Enrichment visualization", "GSEA", "GSEA visualization"
)

for (i in seq_along(scripts)) {
  cat(sprintf("[%02d/11] %s\n", i, labels[i]))
  cat(sprintf("\n[%02d/11] %s\n", i, labels[i]), file = log_file, append = TRUE)
  output <- system2(
    "Rscript",
    c(file.path(workflow_dir, "scripts", scripts[i]), config_file),
    stdout = TRUE,
    stderr = TRUE
  )
  cat(paste(output, collapse = "\n"), "\n", file = log_file, append = TRUE)
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  if (!identical(status, 0L)) {
    message("FAILED at step ", scripts[i])
    message("See log: ", log_file)
    quit(save = "no", status = status)
  }
}

cat("Pipeline completed successfully\n")
cat(paste("Pipeline completed:", format(Sys.time()), "\n"), file = log_file, append = TRUE)
cat("Log:", log_file, "\n")
