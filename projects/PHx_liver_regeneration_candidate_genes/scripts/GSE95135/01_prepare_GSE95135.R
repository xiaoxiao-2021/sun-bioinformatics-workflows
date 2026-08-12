#!/usr/bin/env Rscript

# ============================================================
# GSE95135：整理候选基因表达数据
#
# 这个脚本只做数据整理，不进行统计检验：
#   1. 从作者提供的 log2 RPKM 矩阵中提取候选基因；
#   2. 从样本编号识别 Control、Sham、Post-PH 和时间；
#   3. 保存成长表，供后续画图和比较使用。
# ============================================================

# 需要分析的基因。以后增加基因时，需要同时在 gene_reference 中补充信息。
genes_to_analyse <- c("Adcy7", "Piezo1", "Cxcl2", "Ifrd1")

# 当前基因名、矩阵中的基因名及 Ensembl ID。
# GSE95135 使用较早的 Fam38a 表示现在的 Piezo1。
gene_reference <- data.frame(
  current_symbol = c("Adcy7", "Piezo1", "Cxcl2", "Ifrd1"),
  matrix_symbol  = c("Adcy7", "Fam38a", "Cxcl2", "Ifrd1"),
  ensembl_id     = c(
    "ENSMUSG00000031659",
    "ENSMUSG00000014444",
    "ENSMUSG00000058427",
    "ENSMUSG00000001627"
  ),
  stringsAsFactors = FALSE
)

# 自动确定课题根目录：脚本位于 scripts/GSE95135/ 下。
get_project_root <- function() {
  all_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", all_args, value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- normalizePath(
      sub("^--file=", "", file_arg),
      winslash = "/",
      mustWork = TRUE
    )
    return(normalizePath(
      file.path(dirname(script_path), "..", ".."),
      winslash = "/",
      mustWork = TRUE
    ))
  }
  normalizePath(".", winslash = "/", mustWork = TRUE)
}

# 将样本名（例如 X10H3）拆分为组别、时间和重复编号。
parse_sample_codes <- function(sample_codes) {
  matches <- regexec("^([CSX])(\\d+)([HW])(\\d+)$", sample_codes, perl = TRUE)
  parts <- regmatches(sample_codes, matches)

  if (any(lengths(parts) != 5)) {
    bad_codes <- sample_codes[lengths(parts) != 5]
    stop("无法识别这些样本编号：", paste(bad_codes, collapse = ", "))
  }

  prefix <- vapply(parts, `[`, character(1), 2)
  time_value <- as.integer(vapply(parts, `[`, character(1), 3))
  time_unit <- vapply(parts, `[`, character(1), 4)
  replicate <- as.integer(vapply(parts, `[`, character(1), 5))

  data.frame(
    sample_code = sample_codes,
    group = unname(c(C = "Control", S = "Sham", X = "Post-PH")[prefix]),
    time_value = time_value,
    time_unit = ifelse(time_unit == "H", "hour", "week"),
    time_label = ifelse(time_unit == "H", paste(time_value, "h"), paste(time_value, "wk")),
    time_hours = ifelse(time_unit == "H", time_value, time_value * 168L),
    replicate = replicate,
    stringsAsFactors = FALSE
  )
}

# 逐行读取 gzip 矩阵，找到全部目的基因后立即停止。
# 这样即使当前下载文件只包含矩阵前半部分，也不会把它误用于全基因组分析。
extract_requested_rows <- function(matrix_path, requested_reference) {
  connection <- gzfile(matrix_path, open = "rt")
  on.exit(close(connection), add = TRUE)

  header_line <- readLines(connection, n = 1, warn = FALSE)
  if (length(header_line) == 0) stop("表达矩阵为空：", matrix_path)
  header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]

  found <- list()
  repeat {
    line <- readLines(connection, n = 1, warn = FALSE)
    if (length(line) == 0) break
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3) next

    for (i in seq_len(nrow(requested_reference))) {
      current <- requested_reference$current_symbol[i]
      accepted_symbols <- tolower(c(
        requested_reference$current_symbol[i],
        requested_reference$matrix_symbol[i]
      ))
      is_match <- fields[1] == requested_reference$ensembl_id[i] ||
        tolower(fields[2]) %in% accepted_symbols

      if (is_match) {
        if (length(fields) != length(header)) {
          stop(current, " 行的列数为 ", length(fields),
               "，但表头列数为 ", length(header), "。")
        }
        found[[current]] <- fields
      }
    }

    if (all(requested_reference$current_symbol %in% names(found))) break
  }

  missing <- setdiff(requested_reference$current_symbol, names(found))
  if (length(missing) > 0) {
    stop("矩阵中未找到：", paste(missing, collapse = ", "))
  }

  list(header = header, rows = found)
}

project_root <- get_project_root()
raw_dir <- file.path(project_root, "datasets", "GSE95135", "raw")
processed_dir <- file.path(project_root, "datasets", "GSE95135", "processed")
metadata_dir <- file.path(project_root, "datasets", "GSE95135", "metadata")
figure_dir <- file.path(project_root, "results", "GSE95135", "figures")
table_dir <- file.path(project_root, "results", "GSE95135", "tables")
log_dir <- file.path(project_root, "results", "GSE95135", "logs")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

matrix_candidates <- list.files(
  raw_dir,
  pattern = "^GSE95135_Rib_et_al\\.RPKM_log2.*\\.csv\\.gz$",
  full.names = TRUE
)
if (length(matrix_candidates) == 0) stop("没有找到 GSE95135 log2 RPKM 矩阵。")
matrix_path <- matrix_candidates[which.max(file.info(matrix_candidates)$size)]

requested_reference <- gene_reference[
  gene_reference$current_symbol %in% genes_to_analyse,
  ,
  drop = FALSE
]
extracted <- extract_requested_rows(matrix_path, requested_reference)
sample_codes <- extracted$header[-seq_len(3)]
sample_metadata <- parse_sample_codes(sample_codes)

expression_list <- lapply(genes_to_analyse, function(current_symbol) {
  row <- extracted$rows[[current_symbol]]
  values <- as.numeric(row[-seq_len(3)])
  if (length(values) != nrow(sample_metadata)) {
    stop(current_symbol, " 的表达值数量与样本数量不一致。")
  }
  data.frame(
    sample_metadata,
    gene_id = row[1],
    gene_symbol_current = current_symbol,
    gene_symbol_matrix = row[2],
    gene_type = row[3],
    log2_rpkm = values,
    stringsAsFactors = FALSE
  )
})
expression <- do.call(rbind, expression_list)

write.table(
  expression,
  file.path(processed_dir, "candidate_genes_log2RPKM_long.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  unique(sample_metadata),
  file.path(metadata_dir, "matrix_sample_metadata.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

for (gene in genes_to_analyse) {
  gene_data <- expression[expression$gene_symbol_current == gene, , drop = FALSE]
  write.table(
    gene_data,
    file.path(processed_dir, paste0(gene, "_log2RPKM_long.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

message("数据整理完成：", nrow(expression), " 行")
message("输入矩阵：", matrix_path)
message("输出目录：", processed_dir)

