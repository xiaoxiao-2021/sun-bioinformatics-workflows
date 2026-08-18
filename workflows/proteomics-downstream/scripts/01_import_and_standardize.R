args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("01_import_and_standardize.R requires one config path")
for (package in c("yaml", "readxl")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing R package: ", package)
}

config <- yaml::read_yaml(normalizePath(args[[1]], mustWork = TRUE))
project_dir <- normalizePath(config$project_dir, mustWork = TRUE)
input_file <- file.path(project_dir, config$input_file)
if (!file.exists(input_file)) stop("Input Excel file not found: ", input_file)
if (!config$input_sheet %in% readxl::excel_sheets(input_file)) {
  stop("Excel sheet not found: ", config$input_sheet)
}

metadata <- NULL
if (!is.null(config$sample_metadata)) {
  metadata_file <- file.path(project_dir, config$sample_metadata)
  if (!file.exists(metadata_file)) stop("Sample metadata not found: ", metadata_file)
  metadata <- utils::read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(metadata), c("sample_id", "group"))) {
    stop("Metadata must contain exactly two columns in this order: sample_id, group")
  }
  if (anyNA(metadata$sample_id) || any(!nzchar(metadata$sample_id)) || anyDuplicated(metadata$sample_id)) {
    stop("Metadata sample_id values must be non-empty and unique")
  }
  if (!all(c(config$case_group, config$control_group) %in% metadata$group)) {
    stop("Metadata does not contain both configured case and control groups")
  }
}

raw <- as.data.frame(
  readxl::read_excel(input_file, sheet = config$input_sheet, .name_repair = "minimal"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!nrow(raw)) stop("The configured Excel sheet contains no protein rows")
if (anyDuplicated(names(raw))) stop("The Excel sheet contains duplicated column names")

source_columns <- c(
  config$protein_id_col, config$gene_symbol_col, config$description_col,
  config$unique_peptides_col, config$vendor_logFC_col,
  config$vendor_pvalue_col, config$vendor_qvalue_col
)
sample_ids <- if (is.null(metadata)) character() else metadata$sample_id
gene_id_available <- !is.null(config$gene_id_col) &&
  length(config$gene_id_col) == 1L && nzchar(config$gene_id_col) &&
  config$gene_id_col %in% names(raw)
if (gene_id_available) source_columns <- c(source_columns, config$gene_id_col)
required_columns <- unique(c(source_columns, sample_ids))
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns)) {
  stop("Required Excel column(s) missing: ", paste(missing_columns, collapse = ", "))
}

as_numeric_checked <- function(x, label) {
  value <- suppressWarnings(as.numeric(x))
  bad <- !is.na(x) & nzchar(trimws(as.character(x))) & is.na(value)
  if (any(bad)) warning(sum(bad), " non-numeric value(s) became NA in ", label)
  value
}
clean_entrez <- function(x) {
  value <- trimws(as.character(x))
  value[value %in% c("", "NA", "N/A", "na")] <- NA_character_
  sub("\\.0+$", "", value)
}

map_symbols_to_entrez <- function(symbols) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("Missing R package required for SYMBOL -> ENTREZID mapping: AnnotationDbi")
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("Missing R package required for SYMBOL -> ENTREZID mapping: org.Hs.eg.db")
  mapped <- rep(NA_character_, length(symbols))
  usable <- !is.na(symbols) & nzchar(symbols)
  if (!any(usable)) return(mapped)
  keys <- unique(symbols[usable])
  symbol_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = keys,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  mapped[usable] <- clean_entrez(unname(symbol_map[symbols[usable]]))
  mapped
}

entrez_values <- if (gene_id_available) {
  clean_entrez(raw[[config$gene_id_col]])
} else {
  map_symbols_to_entrez(trimws(as.character(raw[[config$gene_symbol_col]])))
}

# Keep every vendor column, then append a stable downstream schema.
standardized <- data.frame(
  PROTEIN_ID = trimws(as.character(raw[[config$protein_id_col]])),
  SYMBOL = trimws(as.character(raw[[config$gene_symbol_col]])),
  DESCRIPTION = as.character(raw[[config$description_col]]),
  UNIQUE_PEPTIDES = as_numeric_checked(raw[[config$unique_peptides_col]], config$unique_peptides_col),
  ENTREZID = entrez_values,
  logFC = as_numeric_checked(raw[[config$vendor_logFC_col]], config$vendor_logFC_col),
  P.Value = as_numeric_checked(raw[[config$vendor_pvalue_col]], config$vendor_pvalue_col),
  Q.Value = as_numeric_checked(raw[[config$vendor_qvalue_col]], config$vendor_qvalue_col),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
standardized$SYMBOL[!nzchar(standardized$SYMBOL)] <- NA_character_
standardized$PROTEIN_ID[!nzchar(standardized$PROTEIN_ID)] <- NA_character_
for (sample_id in sample_ids) {
  raw[[sample_id]] <- as_numeric_checked(raw[[sample_id]], sample_id)
}

overlap <- intersect(names(raw), names(standardized))
if (length(overlap)) stop("Input already contains reserved standardized column(s): ", paste(overlap, collapse = ", "))
all_proteins <- cbind(raw, standardized)
output_dir <- file.path(project_dir, "datasets", "processed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, paste0(config$dataset_id, "_all_proteins.tsv"))
utils::write.table(all_proteins, output_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

cat("Imported", nrow(all_proteins), "unscreened protein rows\n")
cat("Preserved", ncol(raw), "original columns and added 8 standardized columns\n")
cat("ENTREZID source:", if (gene_id_available) config$gene_id_col else "SYMBOL mapped with org.Hs.eg.db", "\n")
cat("Standardized table:", output_file, "\n")
