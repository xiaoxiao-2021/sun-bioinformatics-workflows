args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Expected one config.yml argument")
cfg <- yaml::read_yaml(normalizePath(args[1], mustWork = TRUE))
project_dir <- normalizePath(cfg$project_dir, mustWork = TRUE)
metadata <- read.delim(file.path(project_dir, cfg$sample_metadata), check.names = FALSE, stringsAsFactors = FALSE)
dat <- readxl::read_excel(file.path(project_dir, cfg$input_file), sheet = 1)

# Build expression matrix in metadata order
sample_cols <- metadata$sample
expr <- as.matrix(dat[, sample_cols, drop = FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- as.character(dat[[cfg$gene_id_col]])

# Low-expression filtering
keep <- rowSums(expr >= cfg$expression_cutoff) >= cfg$min_samples
expr_filtered <- expr[keep, , drop = FALSE]
if (!nrow(expr_filtered)) stop("No genes remain after low-expression filtering")

# log2 transformation
expr_clean <- log2(expr_filtered + 1)
processed_dir <- file.path(project_dir, "datasets", "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
rds_file <- file.path(processed_dir, paste0(cfg$dataset_id, "_bulk_workflow_expression_log2_filtered.rds"))
tsv_file <- file.path(processed_dir, paste0(cfg$dataset_id, "_bulk_workflow_expression_log2_filtered.tsv"))
saveRDS(expr_clean, rds_file)
write.table(data.frame(ENSEMBL = rownames(expr_clean), expr_clean, check.names = FALSE), tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Retained", nrow(expr_clean), "of", nrow(expr), "genes\n")
cat("Saved:", rds_file, "\n")
