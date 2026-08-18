args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("02_data_check.R requires one config path")
if (!requireNamespace("yaml", quietly = TRUE)) stop("Missing R package: yaml")

config <- yaml::read_yaml(normalizePath(args[[1]], mustWork = TRUE))
project_dir <- normalizePath(config$project_dir, mustWork = TRUE)
processed_dir <- file.path(project_dir, "datasets", "processed")
data_file <- file.path(processed_dir, paste0(config$dataset_id, "_all_proteins.tsv"))
if (!file.exists(data_file)) stop("Standardized protein table not found: ", data_file)

data <- utils::read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
metadata <- NULL
if (!is.null(config$sample_metadata)) {
  metadata_file <- file.path(project_dir, config$sample_metadata)
  if (!file.exists(metadata_file)) stop("Sample metadata not found: ", metadata_file)
  metadata <- utils::read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(metadata), c("sample_id", "group"))) stop("Metadata columns must be sample_id and group")
  if (anyDuplicated(metadata$sample_id)) stop("Metadata sample_id values are duplicated")
  if (!all(metadata$sample_id %in% names(data))) stop("Not all metadata samples exist in the standardized table")
  if (!all(metadata$group %in% c(config$case_group, config$control_group))) {
    stop("Metadata contains group(s) other than the configured case/control groups")
  }
}

sample_ids <- if (is.null(metadata)) character() else metadata$sample_id
case_samples <- if (is.null(metadata)) character() else metadata$sample_id[metadata$group == config$case_group]
control_samples <- if (is.null(metadata)) character() else metadata$sample_id[metadata$group == config$control_group]
if (!is.null(metadata) && (length(case_samples) < 2L || length(control_samples) < 2L)) {
  stop("At least two biological replicates are required in each configured group")
}

summary_rows <- list()
add_metric <- function(section, metric, value) {
  summary_rows[[length(summary_rows) + 1L]] <<- data.frame(
    section = section, metric = metric, value = as.character(value),
    stringsAsFactors = FALSE
  )
}
add_metric("proteins", "total_rows", nrow(data))
add_metric("proteins", "unique_nonmissing_PROTEIN_ID", length(unique(stats::na.omit(data$PROTEIN_ID))))
add_metric("proteins", "duplicated_PROTEIN_ID_rows", sum(duplicated(data$PROTEIN_ID) & !is.na(data$PROTEIN_ID)))
add_metric("annotation", "missing_SYMBOL", sum(is.na(data$SYMBOL) | !nzchar(data$SYMBOL)))
add_metric("annotation", "missing_ENTREZID", sum(is.na(data$ENTREZID) | !nzchar(data$ENTREZID)))
add_metric("annotation", "duplicated_ENTREZID_rows", sum(duplicated(data$ENTREZID) & !is.na(data$ENTREZID)))
proteins_with_symbol <- sum(!is.na(data$SYMBOL) & nzchar(data$SYMBOL))
proteins_mapped_entrez <- sum(!is.na(data$ENTREZID) & grepl("^[0-9]+$", data$ENTREZID))
mapping_percentage <- if (nrow(data)) 100 * proteins_mapped_entrez / nrow(data) else 0
mapping_percentage_among_symbol <- if (proteins_with_symbol) 100 * proteins_mapped_entrez / proteins_with_symbol else 0
add_metric("mapping", "total_proteins", nrow(data))
add_metric("mapping", "proteins_with_SYMBOL", proteins_with_symbol)
add_metric("mapping", "proteins_mapped_to_ENTREZID", proteins_mapped_entrez)
add_metric("mapping", "mapping_percentage", signif(mapping_percentage, 5))
add_metric("mapping", "mapping_percentage_among_SYMBOL", signif(mapping_percentage_among_symbol, 5))
add_metric("mapping", "unmapped_proteins", nrow(data) - proteins_mapped_entrez)
add_metric("mapping", "duplicated_ENTREZID_count", sum(duplicated(data$ENTREZID) & !is.na(data$ENTREZID)))

valid_entrez <- !is.na(data$ENTREZID) & grepl("^[0-9]+$", data$ENTREZID)
nonmissing_entrez <- !is.na(data$ENTREZID) & nzchar(data$ENTREZID)
entrez_format_fraction <- if (any(nonmissing_entrez)) mean(valid_entrez[nonmissing_entrez]) else 0
unique_valid_entrez <- length(unique(data$ENTREZID[valid_entrez]))
# Enrichment is considered auditable only when most non-empty IDs are Entrez-like
# and the mapped universe is large enough to support configured gene-set sizes.
entrez_valid_for_enrichment <- entrez_format_fraction >= 0.80 &&
  unique_valid_entrez >= max(10L, as.integer(config$min_GS_size))
add_metric("annotation", "entrez_like_fraction_among_nonmissing", signif(entrez_format_fraction, 4))
add_metric("annotation", "unique_valid_ENTREZID", unique_valid_entrez)
add_metric("annotation", "entrez_valid_for_enrichment", entrez_valid_for_enrichment)

# Convert in place so non-syntactic sample IDs such as "ADCY7-1" are retained.
# Rebuilding with as.data.frame(lapply(...)) would silently change '-' to '.'.
abundance <- data[sample_ids]
abundance[] <- lapply(abundance, function(x) suppressWarnings(as.numeric(x)))
for (sample_id in sample_ids) {
  values <- abundance[[sample_id]]
  finite <- is.finite(values)
  add_metric("sample", paste0(sample_id, ":group"), metadata$group[match(sample_id, metadata$sample_id)])
  add_metric("sample", paste0(sample_id, ":valid_n"), sum(finite))
  add_metric("sample", paste0(sample_id, ":missing_fraction"), signif(mean(!finite), 4))
  add_metric("sample", paste0(sample_id, ":range"), if (any(finite)) paste(signif(range(values[finite]), 5), collapse = " to ") else "NA")
}
if (is.null(metadata)) {
  add_metric("metadata", "sample_metadata_available", FALSE)
  add_metric("abundance", "sample_level_abundance_available", FALSE)
  n_comparable <- 0L
  correlation <- NA_real_
  median_abs_difference <- NA_real_
  direction_concordance <- NA_real_
  abundance_scale_ok <- FALSE
} else {
  add_metric("metadata", "sample_metadata_available", TRUE)
  add_metric("metadata", paste0(config$case_group, "_replicates"), length(case_samples))
  add_metric("metadata", paste0(config$control_group, "_replicates"), length(control_samples))
  add_metric("abundance", "total_NA_or_Inf", sum(!is.finite(as.matrix(abundance))))
  case_mean <- rowMeans(abundance[case_samples], na.rm = FALSE)
  control_mean <- rowMeans(abundance[control_samples], na.rm = FALSE)
  sample_delta <- case_mean - control_mean
  vendor_logfc <- suppressWarnings(as.numeric(data$logFC))
  comparable <- is.finite(sample_delta) & is.finite(vendor_logfc)
  n_comparable <- sum(comparable)
  correlation <- if (n_comparable >= 3L) stats::cor(sample_delta[comparable], vendor_logfc[comparable]) else NA_real_
  median_abs_difference <- if (n_comparable) stats::median(abs(sample_delta[comparable] - vendor_logfc[comparable])) else NA_real_
  nonzero <- comparable & sample_delta != 0 & vendor_logfc != 0
  direction_concordance <- if (any(nonzero)) mean(sign(sample_delta[nonzero]) == sign(vendor_logfc[nonzero])) else NA_real_
  # These explicit thresholds catch raw/intensity-scale input without modifying it.
  abundance_scale_ok <- isTRUE(config$expression_already_log2) && n_comparable >= 10L &&
    is.finite(correlation) && correlation >= 0.90 &&
    is.finite(median_abs_difference) && median_abs_difference <= 0.35 &&
    is.finite(direction_concordance) && direction_concordance >= 0.90
}
add_metric("scale_check", "compared_rows", n_comparable)
add_metric("scale_check", "pearson_correlation", signif(correlation, 5))
add_metric("scale_check", "median_absolute_difference", signif(median_abs_difference, 5))
add_metric("scale_check", "direction_concordance", signif(direction_concordance, 5))
add_metric("scale_check", "abundance_scale_ok", abundance_scale_ok)

qc_summary <- do.call(rbind, summary_rows)
summary_file <- file.path(processed_dir, paste0(config$dataset_id, "_QC_summary.tsv"))
utils::write.table(qc_summary, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
flags <- list(
  entrez_valid_for_enrichment = entrez_valid_for_enrichment,
  abundance_scale_ok = abundance_scale_ok,
  case_samples = case_samples,
  control_samples = control_samples,
  scale_correlation = correlation,
  scale_median_absolute_difference = median_abs_difference,
  scale_direction_concordance = direction_concordance
)
saveRDS(flags, file.path(processed_dir, paste0(config$dataset_id, "_QC_flags.rds")))

if (is.null(metadata)) {
  warning("sample_metadata is null; Volcano and ORA can run, but heatmap and sample-level limma/GSEA are unavailable.")
} else if (abundance_scale_ok) {
  cat("Sample abundance is consistent with log2-scale values; no additional log2 transform will be applied.\n")
} else {
  warning("Abundance values do not agree strongly with vendor logFC. Heatmap and limma/GSEA will be skipped or stopped; values were not transformed.")
}
if (!entrez_valid_for_enrichment) {
  warning("GeneID values are not sufficiently Entrez-like or mapped; ORA/GSEA will stop with a clear error. Volcano visualization may continue.")
}
cat("QC summary:", summary_file, "\n")
