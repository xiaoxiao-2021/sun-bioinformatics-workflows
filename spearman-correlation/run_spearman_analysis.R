# ============================================================
# Spearman correlation workflow：清晰版启动调用脚本
# ============================================================

# ============================================================
# 0. 路径设置
# ============================================================

getwd()

workflow_file <- paste0(
  "spearman-correlation/R/",
  "run_spearman_workflow.R"
)

feature_file <- paste0(
  "spearman-correlation/examples/",
  "feature_df_example.csv"
)

target_file <- paste0(
  "spearman-correlation/examples/",
  "target_df_example.csv"
)

# 可选："csv" 或 "rds"
input_format <- "csv"


# ============================================================
# 1. Spearman 分析参数
# ============================================================

# ------------------------------------------------------------
# 1.1 输入结构
# ------------------------------------------------------------

feature_col <- "feature"
sample_col <- "Sample"


# ------------------------------------------------------------
# 1.2 分析对象
# ------------------------------------------------------------

selected_features <- NULL

target_cols <- c(
  "Score_A"
)


# ------------------------------------------------------------
# 1.3 全样本与组内分析
# ------------------------------------------------------------

run_all_samples <- TRUE
run_group_analysis <- TRUE

group_col <- "Group"
selected_groups <- NULL


# ------------------------------------------------------------
# 1.4 统计阈值
# ------------------------------------------------------------

min_n <- 5
r_cutoff <- 0.5
p_cutoff <- 0.05
padj_cutoff <- 0.05


# ------------------------------------------------------------
# 1.5 样本匹配
# ------------------------------------------------------------

strict_sample_match <- FALSE


# ------------------------------------------------------------
# 1.6 保存设置
# ------------------------------------------------------------

save_results <- FALSE

outdir <- "D:/bioinformatics-results/spearman_analysis"


# ------------------------------------------------------------
# 1.7 运行信息
# ------------------------------------------------------------

verbose <- TRUE


# ============================================================
# 2. 绘图参数
# ============================================================

# ------------------------------------------------------------
# 2.1 总览图
# ------------------------------------------------------------

make_overview_plot <- TRUE
label_top_n_each <- 0


# ------------------------------------------------------------
# 2.2 单 feature-target 图
# ------------------------------------------------------------

run_single_feature_plot <- TRUE

single_feature_name <- "Gene_Pos"
single_target_name <- "Score_A"
single_group_name <- "All"

# "scatter"：连续型 target
# "ordinal"：有序等级 target
single_plot_type <- "scatter"

# 只在 scatter 模式下有意义
single_add_lm <- TRUE


# ============================================================
# 3. 其他可选操作
# ============================================================

run_result_filter_examples <- TRUE
run_manual_validation <- TRUE

validation_feature <- "Gene_Pos"
validation_target <- "Score_A"
validation_group <- "All"


# ============================================================
# 4. 读取输入数据
# ============================================================

if (identical(input_format, "csv")) {

  feature_df <- read.csv(
    feature_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  target_df <- read.csv(
    target_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

} else if (identical(input_format, "rds")) {

  feature_df <- readRDS(feature_file)
  target_df <- readRDS(target_file)

} else {

  stop(
    "input_format 只能设置为 'csv' 或 'rds'。",
    call. = FALSE
  )
}


# ============================================================
# 5. 输入数据检查
# ============================================================

dim(feature_df)
dim(target_df)

head(feature_df)
head(target_df)

str(feature_df)
str(target_df)


# ============================================================
# 6. 加载 workflow 函数
# ============================================================

source(
  workflow_file
)


# ============================================================
# 7. 启动 Spearman workflow
# ============================================================

spearman_result <- run_spearman_workflow(
  feature_df = feature_df,
  target_df = target_df,

  feature_col = feature_col,
  sample_col = sample_col,

  selected_features = selected_features,
  target_cols = target_cols,

  run_all_samples = run_all_samples,
  run_group_analysis = run_group_analysis,

  group_col = group_col,
  selected_groups = selected_groups,

  min_n = min_n,
  r_cutoff = r_cutoff,
  p_cutoff = p_cutoff,
  padj_cutoff = padj_cutoff,

  strict_sample_match = strict_sample_match,

  outdir = outdir,
  save_results = save_results,

  make_overview_plot = make_overview_plot,
  label_top_n_each = label_top_n_each,

  verbose = verbose
)


# ============================================================
# 8. 查看主要结果
# ============================================================

spearman_result
spearman_result$summary

head(
  spearman_result$results
)

spearman_result$significant_p
spearman_result$significant_padj


# ============================================================
# 9. 查看质控信息
# ============================================================

spearman_result$sample_report

spearman_result$features_used
spearman_result$targets_used
spearman_result$groups_used

table(
  spearman_result$results$status,
  useNA = "ifany"
)


# ============================================================
# 10. 查看总览图
# ============================================================

if (isTRUE(make_overview_plot)) {

  names(
    spearman_result$overview_plots
  )

  if (length(spearman_result$overview_plots) > 0L) {
    print(
      spearman_result$overview_plots[[1]]
    )
  }
}


# ============================================================
# 11. 可选：绘制单 feature-target 图
# ============================================================

if (isTRUE(run_single_feature_plot)) {

  p_feature <- plot_spearman_feature(
    workflow_result = spearman_result,

    feature_name = single_feature_name,
    target_name = single_target_name,
    group_name = single_group_name,

    plot_type = single_plot_type,
    add_lm = single_add_lm
  )

  print(
    p_feature
  )
}


# ============================================================
# 12. 可选：筛选结果
# ============================================================

if (isTRUE(run_result_filter_examples)) {

  feature_result <- subset(
    spearman_result$results,
    feature == single_feature_name
  )

  feature_result

  group_result <- subset(
    spearman_result$results,
    group == single_group_name
  )

  head(
    group_result
  )

  significant_result <- subset(
    spearman_result$results,
    status == "OK" &
      significant_by_padj
  )

  head(
    significant_result
  )
}


# ============================================================
# 13. 可选：手工验证
# ============================================================

if (isTRUE(run_manual_validation)) {

  target_df_check <- spearman_result$target_df_aligned

  if (!identical(validation_group, "All")) {

    group_col_check <- spearman_result$parameters$group_col

    target_df_check <- target_df_check[
      !is.na(target_df_check[[group_col_check]]) &
        target_df_check[[group_col_check]] == validation_group,
      ,
      drop = FALSE
    ]
  }

  samples_check <- target_df_check[
    [sample_col]
  ]

  x_check <- as.numeric(
    spearman_result$feature_mat[
      validation_feature,
      samples_check
    ]
  )

  y_check <- target_df_check[
    [validation_target]
  ]

  complete_check <- (
    stats::complete.cases(x_check, y_check) &
      is.finite(x_check) &
      is.finite(y_check)
  )

  manual_result <- stats::cor.test(
    x = x_check[complete_check],
    y = y_check[complete_check],
    method = "spearman",
    exact = FALSE
  )

  manual_result

  workflow_result_check <- subset(
    spearman_result$results,
    feature == validation_feature &
      target == validation_target &
      group == validation_group
  )

  workflow_result_check
}


# ============================================================
# 运行结束
# ============================================================
