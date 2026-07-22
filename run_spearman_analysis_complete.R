# ============================================================
# Spearman correlation workflow：启动调用脚本
# ============================================================
#
# 使用方法：
# 1. 将本文件放在仓库根目录：
#    sun-bioinformatics-workflows/run_spearman_analysis.R
# 2. 在 RStudio / VS Code 中把工作目录设为仓库根目录
# 3. 修改“输入文件路径”和“分析参数”
# 4. 运行本脚本
#
# 函数定义文件：
# spearman-correlation/R/run_spearman_workflow.R
# ============================================================


# ------------------------------------------------------------
# 0. 检查当前工作目录
# ------------------------------------------------------------

getwd()

# 当前目录中应当能够找到：
# spearman-correlation/R/run_spearman_workflow.R


# ------------------------------------------------------------
# 1. 输入文件路径
# ------------------------------------------------------------

# 示例数据路径
feature_file <- paste0(
  "spearman-correlation/examples/",
  "feature_df_example.csv"
)

target_file <- paste0(
  "spearman-correlation/examples/",
  "target_df_example.csv"
)

# 使用真实数据时，替换为实际路径，例如：
#
# feature_file <- "D:/project/data/feature_df.csv"
# target_file  <- "D:/project/data/target_df.csv"


# ------------------------------------------------------------
# 2. 读取输入数据
# ------------------------------------------------------------

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

# 基础检查
dim(feature_df)
dim(target_df)

head(feature_df)
head(target_df)

str(feature_df)
str(target_df)


# ------------------------------------------------------------
# 3. 加载 Spearman workflow 函数
# ------------------------------------------------------------

source(
  "spearman-correlation/R/run_spearman_workflow.R"
)

# source() 只加载函数，不会自动开始分析。
# 加载后应当能够使用：
#
# run_spearman_workflow()
# plot_spearman_overview()
# plot_spearman_feature()


# ------------------------------------------------------------
# 4. 设置本次分析参数
# ------------------------------------------------------------

# feature_df 中保存 feature 名称的列
feature_col <- "feature"

# target_df 中保存样本名称的列
sample_col <- "Sample"


# 4.1 选择 feature
#
# NULL：
# 分析 feature_df 中全部 feature
#
# 字符向量：
# 只分析指定 feature

selected_features <- NULL

# 示例：
# selected_features <- c(
#   "Gene_Pos",
#   "Gene_Neg"
# )


# 4.2 选择 target
#
# NULL：
# 自动分析 target_df 中除 sample/group 外的全部 numeric 列
#
# 字符向量：
# 只分析指定 target

target_cols <- c(
  "Score_A"
)

# 示例：
# target_cols <- c(
#   "Score_A",
#   "Score_B",
#   "ALT"
# )


# 4.3 控制全样本和组内分析

# TRUE：运行全部共同样本分析
run_all_samples <- TRUE

# TRUE：按照 group 分别运行组内分析
run_group_analysis <- TRUE


# 4.4 设置分组参数

# target_df 中保存分组信息的列
group_col <- "Group"

# NULL：
# 自动分析 group_col 中全部非缺失组
#
# 字符向量：
# 只分析指定组

selected_groups <- NULL

# 示例：
# selected_groups <- c(
#   "IA",
#   "ENEG"
# )


# 4.5 统计阈值

# 每次相关性计算所需的最低有效样本数
min_n <- 5

# 相关系数绝对值阈值
r_cutoff <- 0.5

# 原始 p 值阈值
p_cutoff <- 0.05

# BH 校正后 p 值阈值
padj_cutoff <- 0.05


# 4.6 样本匹配策略

# FALSE：
# 报告不匹配样本，但继续使用共同样本分析
#
# TRUE：
# 只要存在不匹配样本就停止

strict_sample_match <- FALSE


# 4.7 输出设置

# FALSE：
# 结果只保存在 R 内存对象中，不写入硬盘
#
# TRUE：
# 自动导出表格、图片、对齐数据和 RDS 文件

save_results <- FALSE

# save_results = TRUE 时才会真正使用 outdir
# 建议正式结果保存在代码仓库外部

outdir <- "D:/bioinformatics-results/spearman_analysis"


# 4.8 绘图设置

# TRUE：
# 在返回结果中生成每个 target × group 的总览图

make_overview_plot <- TRUE

# 每个正、负相关方向标注多少个 feature
# 0 表示不标注

label_top_n_each <- 0


# 4.9 运行信息

# TRUE：
# 在控制台显示样本匹配、分析数量和结果摘要

verbose <- TRUE


# ------------------------------------------------------------
# 5. 启动 Spearman workflow
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# 6. 查看主要结果
# ------------------------------------------------------------

# 直接打印结果对象
spearman_result

# 分析摘要
spearman_result$summary

# 完整相关性结果
head(
  spearman_result$results
)

# 按原始 p 值筛选的结果
spearman_result$significant_p

# 按 BH 校正后 p 值筛选的结果
spearman_result$significant_padj


# ------------------------------------------------------------
# 7. 查看质控信息
# ------------------------------------------------------------

# 样本匹配报告
spearman_result$sample_report

# 本次实际分析的 feature
spearman_result$features_used

# 本次实际分析的 target
spearman_result$targets_used

# 本次实际分析的组
spearman_result$groups_used

# 检查不同状态的结果数量
table(
  spearman_result$results$status,
  useNA = "ifany"
)


# ------------------------------------------------------------
# 8. 查看相关性总览图
# ------------------------------------------------------------

# 查看图的名称
names(
  spearman_result$overview_plots
)

# 显示第一张图
if (length(spearman_result$overview_plots) > 0) {
  print(
    spearman_result$overview_plots[[1]]
  )
}

# 也可以按名称查看，例如：
#
# print(
#   spearman_result$overview_plots[
#     ["All__Score_A"]
#   ]
# )


# ------------------------------------------------------------
# 9. 可选：绘制单个 feature-target 散点图
# ------------------------------------------------------------

# 修改为实际存在的 feature、target 和 group。
# 示例数据中可以使用 Gene_Pos、Score_A、All。

p_feature <- plot_spearman_feature(
  workflow_result = spearman_result,
  feature_name = "Gene_Pos",
  target_name = "Score_A",
  group_name = "All",
  add_lm = TRUE
)

print(p_feature)


# ------------------------------------------------------------
# 10. 可选：筛选指定 feature 或 group
# ------------------------------------------------------------

# 查看某个 feature 在 All、IA、ENEG 中的全部结果
subset(
  spearman_result$results,
  feature == "Gene_Pos"
)

# 只查看 IA 组结果
subset(
  spearman_result$results,
  group == "IA"
)

# 只查看成功计算且校正后显著的结果
subset(
  spearman_result$results,
  status == "OK" &
    significant_by_padj
)


# ------------------------------------------------------------
# 11. 可选：手工验证一个相关性结果
# ------------------------------------------------------------

# 以下示例验证：
# Gene_Pos vs Score_A，All 样本

samples_check <- (
  spearman_result$target_df_aligned[
    [sample_col]
  ]
)

x_check <- as.numeric(
  spearman_result$feature_mat[
    "Gene_Pos",
    samples_check
  ]
)

y_check <- spearman_result$target_df_aligned[
  ["Score_A"]
]

manual_result <- cor.test(
  x = x_check,
  y = y_check,
  method = "spearman",
  exact = FALSE
)

manual_result

subset(
  spearman_result$results,
  feature == "Gene_Pos" &
    target == "Score_A" &
    group == "All"
)


# ============================================================
# 运行结束
# ============================================================
