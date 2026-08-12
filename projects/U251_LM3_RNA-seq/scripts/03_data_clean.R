# ============================================================
# 03_data_clean.R
#
# Purpose:
#   Generate the final cleaned expression matrix for
#   downstream differential expression analysis.
#
# Filtering rule:
#   Expression >= 0.5 in at least 3 samples
#
# Transformation:
#   log2(expression + 1)
# ============================================================


# =========================
# 1. Load packages
# =========================

library(readxl)


# =========================
# 2. Parameters
# =========================

input_file <- "datasets/raw/基因表达水平检测结果-U251.xlsx"

sample_cols <- c(
  "OE1", "OE2", "OE3",
  "NC1", "NC2", "NC3"
)

expression_cutoff <- 0.5
min_samples <- 3


# =========================
# 3. Read raw data
# =========================

dat <- read_excel(input_file)


# =========================
# 4. Basic checks
# =========================

# Check whether all required sample columns exist
stopifnot(
  all(sample_cols %in% colnames(dat))
)

# Check gene_id column
stopifnot(
  "gene_id" %in% colnames(dat)
)

# Check duplicated gene IDs
stopifnot(
  !anyDuplicated(dat$gene_id)
)


# =========================
# 5. Build expression matrix
# =========================

expr <- as.matrix(
  dat[, sample_cols]
)

# Excel columns were originally read as character,
# so explicitly convert them to numeric
storage.mode(expr) <- "numeric"

rownames(expr) <- dat$gene_id


# Check whether numeric conversion produced NA
stopifnot(
  !anyNA(expr)
)


cat(
  "\nOriginal expression matrix:\n"
)

print(
  dim(expr)
)


# =========================
# 6. Low-expression filtering
# =========================

# Keep genes with expression >= 0.5
# in at least 3 of the 6 samples

keep <- rowSums(
  expr >= expression_cutoff
) >= min_samples


expr_filtered <- expr[
  keep,
  ,
  drop = FALSE
]


cat(
  "\nFiltering rule:\n"
)

cat(
  "Expression >=",
  expression_cutoff,
  "in at least",
  min_samples,
  "samples\n"
)


cat(
  "\nGenes retained:\n"
)

print(
  sum(keep)
)


cat(
  "\nGenes removed:\n"
)

print(
  sum(!keep)
)


cat(
  "\nFiltered matrix:\n"
)

print(
  dim(expr_filtered)
)


# =========================
# 7. log2 transformation
# =========================

expr_clean <- log2(
  expr_filtered + 1
)


cat(
  "\nExpression range after log2 transformation:\n"
)

print(
  range(expr_clean)
)


# =========================
# 8. Create output directory
# =========================

dir.create(
  "datasets/processed",
  recursive = TRUE,
  showWarnings = FALSE
)


# =========================
# 9. Save R object
# =========================

saveRDS(
  expr_clean,
  file = "datasets/processed/U251_expression_log2_filtered.rds"
)


# =========================
# 10. Save readable TSV
# =========================

expr_output <- data.frame(
  gene_id = rownames(expr_clean),
  expr_clean,
  check.names = FALSE
)


write.table(
  expr_output,
  file = "datasets/processed/U251_expression_log2_filtered.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# =========================
# 11. Final checks
# =========================

cat(
  "\nFinal clean matrix:\n"
)

print(
  dim(expr_clean)
)


cat(
  "\nFinal sample names:\n"
)

print(
  colnames(expr_clean)
)


cat(
  "\nData cleaning completed successfully.\n"
)