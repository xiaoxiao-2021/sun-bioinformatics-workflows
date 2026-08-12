#!/usr/bin/env Rscript

# ============================================================
# GSE95135：在相同时间点比较 Post-PH 与 Sham
#
# 输出内容：
#   mean_difference_log2 = Post-PH均值 - Sham均值
#   fold_difference      = 2^(mean_difference_log2)
#   p_value              = Welch t检验（探索性）
#   fdr                  = 对全部基因和时间点进行BH校正
#
# 注意：作者提供的是 log2 RPKM，不是原始read counts；
# 因此这里的P值仅供候选筛选，不能替代DESeq2/edgeR原始计数分析。
# ============================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("缺少 ggplot2。请先运行：install.packages('ggplot2')")
}

get_project_root <- function() {
  all_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", all_args, value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
    return(normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE))
  }
  normalizePath(".", winslash = "/", mustWork = TRUE)
}

compare_one_time_point <- function(part) {
  post <- part$log2_rpkm[part$group == "Post-PH"]
  sham <- part$log2_rpkm[part$group == "Sham"]
  if (length(post) < 2 || length(sham) < 2) return(NULL)

  test_result <- stats::t.test(post, sham, var.equal = FALSE)
  difference <- mean(post) - mean(sham)
  data.frame(
    gene = part$gene_symbol_current[1],
    time_hours = part$time_hours[1],
    n_post_ph = length(post),
    n_sham = length(sham),
    mean_post_ph = mean(post),
    mean_sham = mean(sham),
    mean_difference_log2 = difference,
    fold_difference = 2^difference,
    p_value = test_result$p.value,
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
analysis_data <- expression[
  expression$group %in% c("Post-PH", "Sham") &
    expression$time_hours >= 1 & expression$time_hours <= 72,
  ,
  drop = FALSE
]

pieces <- split(
  analysis_data,
  list(analysis_data$gene_symbol_current, analysis_data$time_hours),
  drop = TRUE
)
comparison_list <- lapply(pieces, compare_one_time_point)
comparison_list <- Filter(Negate(is.null), comparison_list)
comparison <- do.call(rbind, comparison_list)
comparison$fdr <- stats::p.adjust(comparison$p_value, method = "BH")
comparison <- comparison[order(comparison$gene, comparison$time_hours), , drop = FALSE]

table_dir <- file.path(project_root, "results", "GSE95135", "tables")
log_dir <- file.path(project_root, "results", "GSE95135", "logs")
figure_dir <- file.path(project_root, "results", "GSE95135", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
  comparison,
  file.path(table_dir, "matched_time_PHx_vs_sham_R.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

comparison$segment_id <- ave(
  comparison$time_hours,
  comparison$gene,
  FUN = function(x) cumsum(c(TRUE, diff(x) > 12))
)
comparison$gene <- factor(comparison$gene, levels = c("Adcy7", "Piezo1", "Cxcl2", "Ifrd1"))

effect_plot <- ggplot2::ggplot(
  comparison,
  ggplot2::aes(
    x = time_hours,
    y = mean_difference_log2,
    group = interaction(gene, segment_id)
  )
) +
  ggplot2::geom_hline(yintercept = 0, colour = "#6B7280", linetype = 2) +
  ggplot2::geom_line(colour = "#7C3AED", linewidth = 0.9) +
  ggplot2::geom_point(colour = "#7C3AED", size = 2.8) +
  ggplot2::facet_wrap(~gene, scales = "free_y", ncol = 2) +
  ggplot2::scale_x_continuous(breaks = c(1, 4, 10, 20, 48)) +
  ggplot2::labs(
    title = "Matched-time Post-PH versus Sham differences (GSE95135)",
    subtitle = "Positive values indicate higher expression after partial hepatectomy",
    x = "Time after surgery (hours)",
    y = "Post-PH minus Sham (log2 RPKM)",
    caption = "Only time points measured in both groups are included; gaps longer than 12 h are not connected."
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(
  file.path(figure_dir, "matched_time_PHx_vs_sham_effect_R.png"),
  effect_plot, width = 10, height = 7, dpi = 300
)
ggplot2::ggsave(
  file.path(figure_dir, "matched_time_PHx_vs_sham_effect_R.pdf"),
  effect_plot, width = 10, height = 7
)

message("匹配时点比较完成。")
message("结果表：", file.path(table_dir, "matched_time_PHx_vs_sham_R.tsv"))

