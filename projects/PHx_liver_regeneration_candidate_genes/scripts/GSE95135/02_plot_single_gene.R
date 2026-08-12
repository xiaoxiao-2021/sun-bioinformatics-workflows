#!/usr/bin/env Rscript

# ============================================================
# GSE95135：绘制一个目的基因的 PHx / Sham 时序图
#
# 最简单的用法：
#   Rscript scripts/GSE95135/02_plot_single_gene.R Piezo1
#
# 如果不在命令中提供基因名，就使用下面的默认基因。
# ============================================================

target_gene <- "Piezo1"       # 可改成 Adcy7、Cxcl2、Ifrd1 或其他已提取基因
maximum_time_hours <- 72       # 只展示术后 72 h 以内
maximum_line_gap <- 12         # 相邻实测点相差超过 12 h 时不连线

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("缺少 ggplot2。请先运行：install.packages('ggplot2')")
}

command_line <- commandArgs(trailingOnly = TRUE)
if (length(command_line) >= 1 && nzchar(command_line[1])) {
  target_gene <- command_line[1]
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

# 对每个组别、每个时间点计算 n、均值、标准差和标准误。
summarise_time_points <- function(data) {
  pieces <- split(data, list(data$group, data$time_hours), drop = TRUE)
  result <- lapply(pieces, function(part) {
    data.frame(
      group = part$group[1],
      time_hours = part$time_hours[1],
      n = nrow(part),
      mean_log2_rpkm = mean(part$log2_rpkm),
      sd_log2_rpkm = stats::sd(part$log2_rpkm),
      sem_log2_rpkm = stats::sd(part$log2_rpkm) / sqrt(nrow(part)),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, result)
  result[order(result$group, result$time_hours), , drop = FALSE]
}

# 如果两个相邻实测点距离过大，给它们分配不同的线段编号，不进行插值连线。
add_line_segments <- function(data, maximum_gap) {
  pieces <- split(data, data$group)
  result <- lapply(pieces, function(part) {
    part <- part[order(part$time_hours), , drop = FALSE]
    part$segment_id <- cumsum(c(TRUE, diff(part$time_hours) > maximum_gap))
    part
  })
  do.call(rbind, result)
}

project_root <- get_project_root()
input_path <- file.path(
  project_root, "datasets", "GSE95135", "processed",
  "candidate_genes_log2RPKM_long.tsv"
)
if (!file.exists(input_path)) {
  stop("没有找到整理后的表达表。请先运行 01_prepare_GSE95135.R")
}

expression <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
gene_data <- expression[expression$gene_symbol_current == target_gene, , drop = FALSE]
if (nrow(gene_data) == 0) stop("整理后的数据中没有基因：", target_gene)

focus <- gene_data[
  (gene_data$group == "Control" & gene_data$time_hours == 0) |
    (gene_data$group %in% c("Post-PH", "Sham") &
       gene_data$time_hours >= 1 & gene_data$time_hours <= maximum_time_hours),
  ,
  drop = FALSE
]
summary_data <- summarise_time_points(focus)

# 0 h Control 是两条手术后曲线共同的起点。
baseline <- summary_data[
  summary_data$group == "Control" & summary_data$time_hours == 0,
  ,
  drop = FALSE
]
line_data <- do.call(rbind, lapply(c("Post-PH", "Sham"), function(current_group) {
  group_data <- summary_data[summary_data$group == current_group, , drop = FALSE]
  baseline_copy <- baseline
  baseline_copy$group <- current_group
  rbind(baseline_copy, group_data)
}))
line_data <- add_line_segments(line_data, maximum_line_gap)

focus$legend_group <- ifelse(focus$group == "Control", "Control 0 h", focus$group)
line_data$legend_group <- line_data$group
plot_colours <- c(
  "Control 0 h" = "#6B7280",
  "Post-PH" = "#D1495B",
  "Sham" = "#277DA1"
)

plot_object <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = focus,
    ggplot2::aes(x = time_hours, y = log2_rpkm, colour = legend_group),
    position = ggplot2::position_jitter(width = 0.20, height = 0),
    size = 2.5, alpha = 0.90
  ) +
  ggplot2::geom_errorbar(
    data = line_data,
    ggplot2::aes(
      x = time_hours,
      ymin = mean_log2_rpkm - sem_log2_rpkm,
      ymax = mean_log2_rpkm + sem_log2_rpkm,
      colour = legend_group
    ),
    width = 0.55, linewidth = 0.65
  ) +
  ggplot2::geom_line(
    data = line_data,
    ggplot2::aes(
      x = time_hours, y = mean_log2_rpkm,
      colour = legend_group,
      group = interaction(group, segment_id)
    ),
    linewidth = 1.05
  ) +
  ggplot2::geom_point(
    data = line_data,
    ggplot2::aes(x = time_hours, y = mean_log2_rpkm, colour = legend_group),
    size = 3.4
  ) +
  ggplot2::scale_colour_manual(values = plot_colours, name = NULL) +
  ggplot2::scale_x_continuous(
    breaks = c(0, 1, 4, 10, 20, 28, 36, 44, 48, 72),
    limits = c(-0.6, maximum_time_hours + 1)
  ) +
  ggplot2::labs(
    title = paste0(target_gene, " expression during mouse liver regeneration (GSE95135)"),
    subtitle = "Post-PH and Sham are joined only across closely spaced measured time points",
    x = "Time after surgery (hours)",
    y = paste0(target_gene, " expression (log2 RPKM)"),
    caption = paste(
      "Small points: individual samples; large points: means; error bars: SEM.",
      "Intervals longer than 12 h are intentionally not connected."
    )
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "top",
    axis.text.x = ggplot2::element_text(colour = "#374151"),
    axis.text.y = ggplot2::element_text(colour = "#374151")
  )

figure_dir <- file.path(project_root, "results", "GSE95135", "figures")
table_dir <- file.path(project_root, "results", "GSE95135", "tables")
log_dir <- file.path(project_root, "results", "GSE95135", "logs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

output_stem <- paste0(target_gene, "_PHx_vs_sham_0_72h_R")
ggplot2::ggsave(file.path(figure_dir, paste0(output_stem, ".png")), plot_object,
                width = 11, height = 6.7, dpi = 300)
ggplot2::ggsave(file.path(figure_dir, paste0(output_stem, ".pdf")), plot_object,
                width = 11, height = 6.7)
write.table(
  summary_data,
  file.path(table_dir, paste0(target_gene, "_timecourse_summary_R.tsv")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

message("绘图完成：", target_gene)
message("图形目录：", figure_dir)
message("汇总表目录：", table_dir)

