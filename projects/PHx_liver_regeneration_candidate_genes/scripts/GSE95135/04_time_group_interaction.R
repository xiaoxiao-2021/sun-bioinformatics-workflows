#!/usr/bin/env Rscript

# ============================================================
# GSE95135：检验“组别 × 时间”交互作用
#
# 问题：Post-PH 与 Sham 的差异是否会随时间改变？
#
# 仅使用两组都实际测量的时间点，不补齐缺失值。
# 对每个基因比较两个线性模型：
#   简化模型：表达量 ~ 组别 + 时间
#   完整模型：表达量 ~ 组别 * 时间
#
# 交互P值较小表示两组轨迹形状可能不同。
# 由于输入是log2 RPKM且样本量较小，本结果属于探索性证据。
# ============================================================

get_project_root <- function() {
  all_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", all_args, value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
    return(normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE))
  }
  normalizePath(".", winslash = "/", mustWork = TRUE)
}

test_one_gene <- function(gene_data) {
  post_times <- unique(gene_data$time_hours[gene_data$group == "Post-PH"])
  sham_times <- unique(gene_data$time_hours[gene_data$group == "Sham"])
  matched_times <- sort(intersect(post_times, sham_times))

  model_data <- gene_data[
    gene_data$group %in% c("Post-PH", "Sham") &
      gene_data$time_hours %in% matched_times,
    ,
    drop = FALSE
  ]
  model_data$group <- factor(model_data$group, levels = c("Sham", "Post-PH"))
  model_data$time_factor <- factor(model_data$time_hours)

  reduced_model <- stats::lm(log2_rpkm ~ group + time_factor, data = model_data)
  full_model <- stats::lm(log2_rpkm ~ group * time_factor, data = model_data)
  model_comparison <- stats::anova(reduced_model, full_model)

  data.frame(
    gene = gene_data$gene_symbol_current[1],
    matched_times = paste(matched_times, collapse = ","),
    n_samples = nrow(model_data),
    numerator_df = model_comparison$Df[2],
    denominator_df = stats::df.residual(full_model),
    interaction_f = model_comparison$F[2],
    interaction_p = model_comparison$`Pr(>F)`[2],
    stringsAsFactors = FALSE
  )
}

project_root <- get_project_root()
input_path <- file.path(
  project_root, "datasets", "GSE95135", "processed",
  "candidate_genes_log2RPKM_long.tsv"
)
if (!file.exists(input_path)) stop("请先运行 01_prepare_GSE95135.R")

expression <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
gene_list <- split(expression, expression$gene_symbol_current)
interaction_results <- do.call(rbind, lapply(gene_list, test_one_gene))
interaction_results$interaction_fdr <- stats::p.adjust(
  interaction_results$interaction_p,
  method = "BH"
)
interaction_results <- interaction_results[order(interaction_results$interaction_fdr), , drop = FALSE]

table_dir <- file.path(project_root, "results", "GSE95135", "tables", "R")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(table_dir, "time_by_group_interaction_R.tsv")
write.table(
  interaction_results,
  output_path,
  sep = "\t", row.names = FALSE, quote = FALSE
)

message("组别×时间交互检验完成。")
message("结果表：", output_path)

