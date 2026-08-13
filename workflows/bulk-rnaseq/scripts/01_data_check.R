args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Expected one config.yml argument")
cfg <- yaml::read_yaml(normalizePath(args[1], mustWork = TRUE))
project_dir <- normalizePath(cfg$project_dir, mustWork = TRUE)
input_file <- file.path(project_dir, cfg$input_file)
metadata_file <- file.path(project_dir, cfg$sample_metadata)

# Read sample metadata
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)
metadata <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("sample", "group") %in% names(metadata))) stop("Metadata must contain sample and group columns")
if (anyNA(metadata$sample) || any(metadata$sample == "") || anyDuplicated(metadata$sample)) {
  stop("Metadata sample values must be non-empty and unique")
}
if (!cfg$case_group %in% metadata$group) stop("case_group not found in metadata: ", cfg$case_group)
if (!cfg$control_group %in% metadata$group) stop("control_group not found in metadata: ", cfg$control_group)
unexpected_groups <- setdiff(unique(metadata$group), c(cfg$case_group, cfg$control_group))
if (length(unexpected_groups)) stop("Metadata contains group(s) outside this two-group comparison: ", paste(unexpected_groups, collapse = ", "))
group_n <- table(metadata$group)
if (any(group_n[c(cfg$case_group, cfg$control_group)] < 2)) stop("case and control groups each require at least 2 samples")
if (cfg$min_samples < 1 || cfg$min_samples > nrow(metadata)) stop("min_samples is outside the sample count")

# Read expression matrix and check sample columns
if (!file.exists(input_file)) stop("Expression file not found: ", input_file)
dat <- readxl::read_excel(input_file, sheet = 1)
if (!cfg$gene_id_col %in% names(dat)) stop("gene_id_col not found: ", cfg$gene_id_col)
expression_samples <- setdiff(names(dat), cfg$gene_id_col)
missing_expression <- setdiff(metadata$sample, expression_samples)
missing_metadata <- setdiff(expression_samples, metadata$sample)
if (length(missing_expression)) stop("Metadata sample(s) absent from expression matrix: ", paste(missing_expression, collapse = ", "))
if (length(missing_metadata)) stop("Expression sample(s) absent from metadata: ", paste(missing_metadata, collapse = ", "))

# Match sample order with metadata
sample_cols <- metadata$sample
expr <- as.matrix(dat[, sample_cols, drop = FALSE])
suppressWarnings(storage.mode(expr) <- "numeric")
gene_id <- as.character(dat[[cfg$gene_id_col]])
if (anyNA(gene_id) || any(gene_id == "")) stop("Gene IDs contain missing or empty values")
if (anyDuplicated(gene_id)) stop("Duplicated gene IDs: ", sum(duplicated(gene_id)))
if (anyNA(expr)) stop("Expression matrix contains NA or non-numeric values")
if (any(!is.finite(expr))) stop("Expression matrix contains non-finite values")
if (any(expr < 0)) stop("Expression matrix contains negative values before log2 transformation")

cat("Data dimensions:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("Sample order:", paste(colnames(expr), collapse = ", "), "\n")
cat("Groups:\n"); print(group_n)
cat("NA count:", sum(is.na(expr)), "\n")
cat("Duplicated gene IDs:", sum(duplicated(gene_id)), "\n")
cat("Zero proportion:", mean(expr == 0), "\n")
cat("Expression range:", paste(range(expr), collapse = " to "), "\n")
cat("Non-zero quantiles:\n"); print(quantile(expr[expr > 0], c(0, .25, .5, .75, .9, .95, .99, 1)))
cat("Genes detected in each number of samples:\n"); print(table(rowSums(expr > 0)))
cat("Genes passing configured filter:", sum(rowSums(expr >= cfg$expression_cutoff) >= cfg$min_samples), "\n")
