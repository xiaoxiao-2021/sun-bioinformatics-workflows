args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript workflows/proteomics-downstream/run_pipeline.R path/to/config.yml")
}
if (!requireNamespace("yaml", quietly = TRUE)) stop("Missing R package: yaml")

config_file <- normalizePath(args[[1]], mustWork = TRUE)
config <- yaml::read_yaml(config_file)
required <- c(
  "dataset_id", "project_dir", "input_file", "input_sheet",
  "sample_metadata", "gene_id_col", "case_group", "control_group", "protein_id_col",
  "gene_symbol_col", "description_col", "unique_peptides_col",
  "vendor_logFC_col", "vendor_pvalue_col", "vendor_qvalue_col",
  "expression_already_log2", "deg_pvalue_cutoff", "deg_qvalue_cutoff",
  "deg_logFC_cutoff", "volcano_qvalue_cutoff", "volcano_logFC_cutoff",
  "volcano_label_n", "heatmap_top_n", "enrichment_gene_logFC_cutoff",
  "enrichment_pvalue_cutoff", "enrichment_qvalue_cutoff", "min_GS_size",
  "max_GS_size", "show_category_n", "gsea_pvalue_cutoff",
  "gsea_FDR_cutoff", "gsea_min_GS_size", "gsea_max_GS_size",
  "gsea_show_category_n", "gsea_draw_curves", "gsea_curve_n", "run_gsea"
)
missing_keys <- required[!required %in% names(config)]
if (length(missing_keys)) stop("Missing config key(s): ", paste(missing_keys, collapse = ", "))
if (!grepl("^[A-Za-z0-9._-]+$", config$dataset_id)) stop("dataset_id contains unsafe characters")
if (identical(config$case_group, config$control_group)) stop("case_group and control_group must differ")
if (!is.logical(config$gsea_draw_curves) || length(config$gsea_draw_curves) != 1L || is.na(config$gsea_draw_curves)) {
  stop("gsea_draw_curves must be true or false")
}
if (!is.logical(config$run_gsea) || length(config$run_gsea) != 1L || is.na(config$run_gsea)) {
  stop("run_gsea must be true or false")
}

project_dir <- normalizePath(config$project_dir, mustWork = TRUE)
log_dir <- file.path(project_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, "pipeline.log")
writeLines(paste("Pipeline started:", format(Sys.time())), log_file)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
workflow_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE))
scripts <- c(
  "01_import_and_standardize.R", "02_data_check.R", "03_visualization.R",
  "04_ORA.R", "05_ORA_visualization.R"
)
labels <- c(
  "Import and standardize", "Data checks", "Vendor-statistic visualization",
  "ORA", "ORA visualization"
)
if (isTRUE(config$run_gsea)) {
  scripts <- c(scripts, "06_GSEA.R", "07_GSEA_visualization.R")
  labels <- c(labels, "GSEA", "GSEA visualization")
} else {
  cat("GSEA skipped by config: run_gsea = false\n")
  cat("GSEA skipped by config: run_gsea = false\n", file = log_file, append = TRUE)
}

for (i in seq_along(scripts)) {
  header <- sprintf("[%d/%d] %s", i, length(scripts), labels[[i]])
  cat(header, "\n")
  cat("\n", header, "\n", file = log_file, append = TRUE, sep = "")
  output <- suppressWarnings(system2(
    "Rscript",
    c(file.path(workflow_dir, "scripts", scripts[[i]]), config_file),
    stdout = TRUE,
    stderr = TRUE
  ))
  if (length(output)) {
    cat(paste(output, collapse = "\n"), "\n", file = log_file, append = TRUE)
  }
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  if (!identical(as.integer(status), 0L)) {
    cat("FAILED STEP: ", scripts[[i]], "\n", sep = "", file = stderr())
    if (length(output)) cat(paste(output, collapse = "\n"), "\n", file = stderr())
    cat("Log: ", log_file, "\n", sep = "", file = stderr())
    quit(save = "no", status = as.integer(status))
  }
}

cat("Pipeline completed successfully\n")
cat("Pipeline completed: ", format(Sys.time()), "\n", file = log_file, append = TRUE, sep = "")
cat("Log:", log_file, "\n")
