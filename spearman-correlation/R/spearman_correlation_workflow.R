# ============================================================
# Spearman correlation workflow
# ============================================================
# 功能：
# 1. 接收两个标准 data.frame：
#    feature_df：第一列/指定列为 feature 名称，其余列为样本
#    target_df ：每一行为一个样本，包含 sample、target、group 等列
# 2. 检查并统一样本
# 3. 支持全样本分析和组内分析
# 4. 支持全部或指定 feature 分析
# 5. 支持全部或指定 target 批量分析
# 6. 支持全部或指定 group 分析
# 7. 输出 rho、pvalue、padj、n、direction 等结果
# 8. 可选输出相关性总览图和单 feature 散点图
#
# 重要：
# 本脚本规定输入接口为 data.frame，
# 但内部计算会将 feature 数值部分转换为 numeric matrix，
# 最终进入 cor.test() 的仍然是两个 numeric vector。
# ============================================================


# ============================================================
# 0. 加载 R 包
# ============================================================

# 第一次使用时可取消下面一行注释安装所需 R 包
# install.packages(c("dplyr", "tibble", "purrr", "readr", "ggplot2", "ggrepel"))

library(dplyr)
library(tibble)
library(purrr)
library(readr)
library(ggplot2)
library(ggrepel)


# ============================================================
# 1. 准备输入对象
# ============================================================

# ------------------------------------------------------------
# 【必须准备】feature_df
# ------------------------------------------------------------
# 标准格式：
# 1. 每一行代表一个 feature
# 2. feature 名称保存在 feature_col 指定的列中
# 3. 其余列全部为样本
# 4. 所有样本列必须为 numeric
#
# 示例：
#
# feature   S01   S02   S03
# Gene_A    4.2   5.1   3.8
# Gene_B    2.3   2.8   4.1
#
# 若原始对象是 numeric matrix，且行名是 feature：
#
# feature_df <- tibble::rownames_to_column(
#   as.data.frame(feature_mat, check.names = FALSE),
#   var = "feature"
# )
#
# 若从 CSV 读取：
#
# feature_df <- read.csv(
#   "feature_df.csv",
#   check.names = FALSE,
#   stringsAsFactors = FALSE
# )
#
# check.names = FALSE 很重要：
# 可避免 R 自动修改样本名，例如把 "-" 改成 "."


# ------------------------------------------------------------
# 【必须准备】target_df
# ------------------------------------------------------------
# 标准格式：
# 1. 每一行代表一个样本
# 2. sample_col 指定的列保存样本名
# 3. target_cols 指定需要分析的一个或多个目标变量
# 4. group_col 可选，用于组内分析
#
# 示例：
#
# sample   Score_A   Score_B   group
# S01      0.42      1.15      IA
# S02      0.81      0.94      IA
# S03      1.10      0.75      ENEG
#
# 若从 CSV 读取：
#
# target_df <- read.csv(
#   "target_df.csv",
#   check.names = FALSE,
#   stringsAsFactors = FALSE
# )


# ============================================================
# 2. 参数区
# ============================================================

# ------------------------------------------------------------
# 【必须检查】列名参数
# ------------------------------------------------------------

# feature_df 中保存 feature 名称的列名
feature_col <- "feature"

# target_df 中保存样本名称的列名
sample_col <- "sample"

# ------------------------------------------------------------
# 【按需修改】选择需要分析的 feature
# ------------------------------------------------------------
#
# 模式 1：分析 feature_df 中全部 feature
selected_features <- NULL
#
# 模式 2：只分析指定 feature
# selected_features <- c(
#   "RIPOR2",
#   "TRAF1",
#   "BCL11B"
# )
#
# 说明：
# NULL 表示全量分析全部 feature；
# 字符向量表示仅分析列出的 feature。
# 程序会按照 selected_features 中给出的顺序提取 feature。

# ------------------------------------------------------------
# 【按需修改】选择需要分析的 target
# ------------------------------------------------------------
#
# 模式 1：分析 target_df 中全部可用 numeric target
# target_cols <- NULL
#
# 模式 2：只分析指定 target
target_cols <- c(
  "T_Cell_Exhaustion_Progenitor_like"
  # "Score_B",
  # "ALT"
)
#
# 说明：
# NULL 表示自动选择 target_df 中除 sample_col 和 group_col
# 之外的全部 numeric 列；
# 字符向量表示仅分析列出的 target。

# ------------------------------------------------------------
# 【按需修改】是否进行全样本分析
# ------------------------------------------------------------

# TRUE：运行所有共同样本的相关性分析
# FALSE：不运行全样本分析
run_all_samples <- TRUE

# ------------------------------------------------------------
# 【按需修改】是否进行组内分析
# ------------------------------------------------------------

# TRUE：按照 group_col 进行组内相关性分析
# FALSE：不进行组内分析
run_group_analysis <- TRUE

# target_df 中保存分组信息的列名
# 仅当 run_group_analysis = TRUE 时使用
group_col <- "Group"

# 需要分析的组
#
# 模式 1：分析 group_col 中全部非缺失组
# selected_groups <- NULL
#
# 模式 2：只分析指定组
selected_groups <- c(
  "IA",
  "ENEG"
)
#
# 说明：
# NULL 表示自动分析 group_col 中全部非缺失组；
# 字符向量表示仅分析列出的组。

# ------------------------------------------------------------
# 【按需修改】统计参数
# ------------------------------------------------------------

# 最小有效样本数
# 每一对 feature-target 在去除缺失值后，n 小于此值时不计算
min_n <- 5

# 相关系数筛选阈值
r_cutoff <- 0.5

# 原始 p 值筛选阈值
p_cutoff <- 0.05

# BH 校正后 p 值筛选阈值
padj_cutoff <- 0.05

# ------------------------------------------------------------
# 【按需修改】样本不匹配时的处理方式
# ------------------------------------------------------------

# TRUE：
# 只要 feature_df 和 target_df 存在不匹配样本，就立即停止
#
# FALSE：
# 显示不匹配样本，但继续使用共同样本分析
#
# 探索性分析通常可设为 FALSE
# 正式分析前建议确认不匹配样本的原因
strict_sample_match <- FALSE

# ------------------------------------------------------------
# 【按需修改】输出路径
# ------------------------------------------------------------

outdir <- "results/spearman_correlation"

table_dir <- file.path(outdir, "tables")
plot_dir  <- file.path(outdir, "plots")
input_dir <- file.path(outdir, "aligned_inputs")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 【按需修改】画图参数
# ------------------------------------------------------------

# 是否绘制每个 target × group 的相关性总览图
make_overview_plot <- TRUE

# 每个方向标注多少个基因
# 0 表示不标注
# 基因较多时建议设为 0、5 或 10
label_top_n_each <- 0


# ============================================================
# 3. 输入对象与列名检查
# ============================================================

# ------------------------------------------------------------
# 3.1 检查对象类型
# ------------------------------------------------------------

if (!exists("feature_df")) {
  stop("未找到 feature_df，请先读取或创建 feature_df。")
}

if (!exists("target_df")) {
  stop("未找到 target_df，请先读取或创建 target_df。")
}

if (!is.data.frame(feature_df)) {
  stop("feature_df 必须是 data.frame。")
}

if (!is.data.frame(target_df)) {
  stop("target_df 必须是 data.frame。")
}

# ------------------------------------------------------------
# 3.2 清理列名和关键名称中的首尾空格
# ------------------------------------------------------------

colnames(feature_df) <- trimws(colnames(feature_df))
colnames(target_df)  <- trimws(colnames(target_df))

# ------------------------------------------------------------
# 3.3 检查必须存在的列
# ------------------------------------------------------------

if (!feature_col %in% colnames(feature_df)) {
  stop(
    "feature_df 中不存在 feature_col 指定的列：",
    feature_col
  )
}

if (!sample_col %in% colnames(target_df)) {
  stop(
    "target_df 中不存在 sample_col 指定的列：",
    sample_col
  )
}

if (run_group_analysis) {
  if (!group_col %in% colnames(target_df)) {
    stop(
      "run_group_analysis = TRUE，但 target_df 中不存在 group_col：",
      group_col
    )
  }
}

# ------------------------------------------------------------
# 解析 target_cols
# ------------------------------------------------------------
#
# target_cols = NULL：
# 自动选择除 sample_col 和 group_col 外的全部 numeric 列
#
# target_cols = 字符向量：
# 仅分析指定 target

if (is.null(target_cols)) {

  excluded_target_cols <- sample_col

  if (group_col %in% colnames(target_df)) {
    excluded_target_cols <- c(
      excluded_target_cols,
      group_col
    )
  }

  candidate_target_cols <- setdiff(
    colnames(target_df),
    excluded_target_cols
  )

  if (length(candidate_target_cols) == 0) {
    stop(
      "target_df 中除 sample/group 外没有可作为 target 的列。"
    )
  }

  candidate_numeric_check <- vapply(
    target_df[, candidate_target_cols, drop = FALSE],
    is.numeric,
    logical(1)
  )

  target_cols <- candidate_target_cols[
    candidate_numeric_check
  ]

  if (length(target_cols) == 0) {
    stop(
      "target_cols = NULL，但未检测到可用的 numeric target 列。"
    )
  }

  message(
    "target_cols = NULL，将分析全部 numeric target：",
    paste(target_cols, collapse = ", ")
  )

} else {

  target_cols <- trimws(
    as.character(target_cols)
  )

  if (
    length(target_cols) == 0 ||
    any(is.na(target_cols)) ||
    any(target_cols == "")
  ) {
    stop(
      "target_cols 不能为空；请设置为 NULL 或有效的字符向量。"
    )
  }

  if (anyDuplicated(target_cols) > 0) {
    stop(
      "target_cols 中存在重复 target：",
      paste(
        unique(target_cols[duplicated(target_cols)]),
        collapse = ", "
      )
    )
  }

  missing_target_cols <- setdiff(
    target_cols,
    colnames(target_df)
  )

  if (length(missing_target_cols) > 0) {
    stop(
      "target_df 中缺少以下 target 列：",
      paste(missing_target_cols, collapse = ", ")
    )
  }

  invalid_target_cols <- intersect(
    target_cols,
    c(sample_col, group_col)
  )

  if (length(invalid_target_cols) > 0) {
    stop(
      "以下列不能作为 target：",
      paste(invalid_target_cols, collapse = ", "),
      "。sample/group 列只用于样本匹配或分组。"
    )
  }
}

# ------------------------------------------------------------
# 3.4 清理 feature、sample 和 group 内容中的首尾空格
# ------------------------------------------------------------

feature_df[[feature_col]] <- trimws(
  as.character(feature_df[[feature_col]])
)

target_df[[sample_col]] <- trimws(
  as.character(target_df[[sample_col]])
)

if (run_group_analysis) {
  target_df[[group_col]] <- trimws(
    as.character(target_df[[group_col]])
  )
}

# ------------------------------------------------------------
# 3.5 检查空名称和缺失名称
# ------------------------------------------------------------

if (
  any(is.na(feature_df[[feature_col]])) ||
  any(feature_df[[feature_col]] == "")
) {
  stop("feature_df 的 feature 名称中存在 NA 或空字符串。")
}

if (
  any(is.na(target_df[[sample_col]])) ||
  any(target_df[[sample_col]] == "")
) {
  stop("target_df 的 sample 名称中存在 NA 或空字符串。")
}

# ------------------------------------------------------------
# 3.6 检查重复 feature 和重复 sample
# ------------------------------------------------------------

duplicated_features <- unique(
  feature_df[[feature_col]][duplicated(feature_df[[feature_col]])]
)

if (length(duplicated_features) > 0) {
  stop(
    "feature_df 中存在重复 feature，请先合并或去重：",
    paste(head(duplicated_features, 20), collapse = ", ")
  )
}

duplicated_target_samples <- unique(
  target_df[[sample_col]][duplicated(target_df[[sample_col]])]
)

if (length(duplicated_target_samples) > 0) {
  stop(
    "target_df 中存在重复 sample，每个样本必须只对应一行：",
    paste(head(duplicated_target_samples, 20), collapse = ", ")
  )
}

# feature_df 中除 feature_col 外，其余列均被视为样本列
feature_samples <- setdiff(
  colnames(feature_df),
  feature_col
)

if (length(feature_samples) == 0) {
  stop("feature_df 中没有检测到样本列。")
}

duplicated_feature_samples <- unique(
  feature_samples[duplicated(feature_samples)]
)

if (length(duplicated_feature_samples) > 0) {
  stop(
    "feature_df 中存在重复样本列名：",
    paste(duplicated_feature_samples, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3.7 检查数值列
# ------------------------------------------------------------

# feature_df 中所有样本列都必须为 numeric
feature_numeric_check <- vapply(
  feature_df[, feature_samples, drop = FALSE],
  is.numeric,
  logical(1)
)

if (!all(feature_numeric_check)) {
  non_numeric_feature_cols <- names(
    feature_numeric_check[!feature_numeric_check]
  )

  stop(
    "feature_df 中以下样本列不是 numeric：",
    paste(non_numeric_feature_cols, collapse = ", "),
    "\n请检查是否存在 '-', 'unknown', 空字符串或单位字符。"
  )
}

# target_cols 必须为 numeric
target_numeric_check <- vapply(
  target_df[, target_cols, drop = FALSE],
  is.numeric,
  logical(1)
)

if (!all(target_numeric_check)) {
  non_numeric_target_cols <- names(
    target_numeric_check[!target_numeric_check]
  )

  stop(
    "target_df 中以下 target 列不是 numeric：",
    paste(non_numeric_target_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3.8 根据 selected_features 选择分析范围
# ------------------------------------------------------------
#
# selected_features = NULL：
# 分析 feature_df 中全部 feature
#
# selected_features = 字符向量：
# 仅分析指定 feature

if (is.null(selected_features)) {

  feature_df_selected <- feature_df

  message(
    "selected_features = NULL，将分析全部 feature；数量：",
    nrow(feature_df_selected)
  )

} else {

  selected_features <- trimws(
    as.character(selected_features)
  )

  if (
    length(selected_features) == 0 ||
    any(is.na(selected_features)) ||
    any(selected_features == "")
  ) {
    stop(
      "selected_features 不能为空；请设置为 NULL 或有效字符向量。"
    )
  }

  if (anyDuplicated(selected_features) > 0) {
    stop(
      "selected_features 中存在重复 feature：",
      paste(
        unique(
          selected_features[
            duplicated(selected_features)
          ]
        ),
        collapse = ", "
      )
    )
  }

  missing_selected_features <- setdiff(
    selected_features,
    feature_df[[feature_col]]
  )

  if (length(missing_selected_features) > 0) {
    stop(
      "以下 selected_features 不存在于 feature_df：",
      paste(
        missing_selected_features,
        collapse = ", "
      )
    )
  }

  feature_df_selected <- feature_df[
    match(
      selected_features,
      feature_df[[feature_col]]
    ),
    ,
    drop = FALSE
  ]

  message(
    "仅分析指定 feature；数量：",
    nrow(feature_df_selected)
  )
}

# ------------------------------------------------------------
# 3.9 解析并检查 selected_groups
# ------------------------------------------------------------

if (run_group_analysis) {

  print(
    table(
      target_df[[group_col]],
      useNA = "ifany"
    )
  )

  available_groups <- unique(
    target_df[[group_col]][
      !is.na(target_df[[group_col]]) &
        target_df[[group_col]] != ""
    ]
  )

  if (length(available_groups) == 0) {
    stop(
      "group_col 中没有可用的非缺失分组。"
    )
  }

  if (is.null(selected_groups)) {

    selected_groups <- available_groups

    message(
      "selected_groups = NULL，将分析全部组：",
      paste(selected_groups, collapse = ", ")
    )

  } else {

    selected_groups <- trimws(
      as.character(selected_groups)
    )

    if (
      length(selected_groups) == 0 ||
      any(is.na(selected_groups)) ||
      any(selected_groups == "")
    ) {
      stop(
        "selected_groups 不能为空；请设置为 NULL 或有效字符向量。"
      )
    }

    if (anyDuplicated(selected_groups) > 0) {
      stop(
        "selected_groups 中存在重复组：",
        paste(
          unique(
            selected_groups[
              duplicated(selected_groups)
            ]
          ),
          collapse = ", "
        )
      )
    }

    missing_groups <- setdiff(
      selected_groups,
      available_groups
    )

    if (length(missing_groups) > 0) {
      stop(
        "selected_groups 中以下组不存在于 target_df：",
        paste(missing_groups, collapse = ", ")
      )
    }
  }
}


# ============================================================
# 4. 样本匹配与顺序统一
# ============================================================

# target_df 中的样本名称
target_samples <- target_df[[sample_col]]

# ------------------------------------------------------------
# 4.1 使用 setdiff() 检查不匹配样本
# ------------------------------------------------------------

# target_df 中存在，但 feature_df 中不存在的样本
only_in_target <- setdiff(
  target_samples,
  feature_samples
)

# feature_df 中存在，但 target_df 中不存在的样本
only_in_feature <- setdiff(
  feature_samples,
  target_samples
)

message(
  "target_df 中有、feature_df 中没有的样本数：",
  length(only_in_target)
)

if (length(only_in_target) > 0) {
  print(only_in_target)
}

message(
  "feature_df 中有、target_df 中没有的样本数：",
  length(only_in_feature)
)

if (length(only_in_feature) > 0) {
  print(only_in_feature)
}

# strict_sample_match = TRUE 时，不允许任何样本不匹配
if (
  strict_sample_match &&
  (length(only_in_target) > 0 || length(only_in_feature) > 0)
) {
  stop(
    "检测到不匹配样本。请检查样本名称，",
    "或将 strict_sample_match 设置为 FALSE 后使用共同样本。"
  )
}

# ------------------------------------------------------------
# 4.2 使用 intersect() 获取共同样本
# ------------------------------------------------------------

# intersect(x, y) 返回 x 和 y 共有的元素
# 返回顺序以第一个向量 feature_samples 的顺序为基础
common_samples <- intersect(
  feature_samples,
  target_samples
)

message("共同样本数：", length(common_samples))

if (length(common_samples) < min_n) {
  stop(
    "共同样本数为 ",
    length(common_samples),
    "，小于 min_n = ",
    min_n,
    "，无法继续分析。"
  )
}

# ------------------------------------------------------------
# 4.3 按照 common_samples 同时整理两个数据框
# ------------------------------------------------------------

# feature_df：
# 保留 feature 列，并按照 common_samples 的顺序排列样本列
feature_df_aligned <- feature_df_selected[
  ,
  c(feature_col, common_samples),
  drop = FALSE
]

# target_df：
# match() 返回 common_samples 在 target_df$sample 中的位置
# 从而按照 common_samples 的顺序重新排列 target_df 的行
target_df_aligned <- target_df[
  match(common_samples, target_df[[sample_col]]),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 4.4 最终确认样本顺序完全一致
# ------------------------------------------------------------

sample_order_ok <- identical(
  colnames(feature_df_aligned)[-1],
  target_df_aligned[[sample_col]]
)

print(sample_order_ok)

if (!sample_order_ok) {
  stop("样本顺序统一失败，请检查样本名称。")
}

# 保存整理后的标准输入，便于 Debug 和复现
write_csv(
  feature_df_aligned,
  file.path(input_dir, "feature_df_aligned.csv")
)

write_csv(
  target_df_aligned,
  file.path(input_dir, "target_df_aligned.csv")
)

# 保存本次实际使用的 feature、target 和 group 列表
write_csv(
  tibble(
    feature = feature_df_aligned[[feature_col]]
  ),
  file.path(input_dir, "features_used.csv")
)

write_csv(
  tibble(
    target = target_cols
  ),
  file.path(input_dir, "targets_used.csv")
)

if (run_group_analysis) {
  write_csv(
    tibble(
      group = selected_groups
    ),
    file.path(input_dir, "groups_used.csv")
  )
}


# ============================================================
# 5. 将 feature_df 转为内部计算矩阵
# ============================================================

# 数据框适合输入、检查和保存
# numeric matrix 更适合批量按行计算
feature_mat <- as.matrix(
  feature_df_aligned[
    ,
    common_samples,
    drop = FALSE
  ]
)

# 前面已检查所有样本列为 numeric
# 这里强制设置为 double，确保适合数值计算
storage.mode(feature_mat) <- "double"

# 将 feature 名称设置为矩阵行名
rownames(feature_mat) <- feature_df_aligned[[feature_col]]

# 检查矩阵结构
print(dim(feature_mat))
print(feature_mat[seq_len(min(3, nrow(feature_mat))),
                  seq_len(min(3, ncol(feature_mat))),
                  drop = FALSE])


# ============================================================
# 6. 定义单个 target × group 的相关性分析函数
# ============================================================

run_one_target_group <- function(
    feature_mat,
    target_sub,
    target_col,
    sample_col,
    group_label,
    min_n = 5
) {

  # 当前分析所使用的样本
  samples_use <- target_sub[[sample_col]]

  # 按当前样本顺序提取 feature 矩阵
  feature_sub <- feature_mat[
    ,
    samples_use,
    drop = FALSE
  ]

  # 最终确认矩阵列顺序与 target_sub 样本顺序完全一致
  stopifnot(
    identical(
      colnames(feature_sub),
      samples_use
    )
  )

  # 提取当前 target 数值向量
  y <- target_sub[[target_col]]

  # 对 feature_sub 的每一行分别进行相关性计算
  result_list <- lapply(
    seq_len(nrow(feature_sub)),
    function(i) {

      feature_name <- rownames(feature_sub)[i]

      # 当前 feature 在所有样本中的数值向量
      x <- as.numeric(feature_sub[i, ])

      # 只保留 x 和 y 同时不缺失的样本
      complete_idx <- complete.cases(x, y)

      x_complete <- x[complete_idx]
      y_complete <- y[complete_idx]

      # 实际参与本次计算的样本数
      n_complete <- sum(complete_idx)

      # 样本数不足时不计算
      if (n_complete < min_n) {
        return(
          tibble(
            feature = feature_name,
            target = target_col,
            group = group_label,
            rho = NA_real_,
            pvalue = NA_real_,
            n = n_complete,
            status = "insufficient_n"
          )
        )
      }

      # feature 为常数时无法计算相关性
      if (length(unique(x_complete)) < 2) {
        return(
          tibble(
            feature = feature_name,
            target = target_col,
            group = group_label,
            rho = NA_real_,
            pvalue = NA_real_,
            n = n_complete,
            status = "constant_feature"
          )
        )
      }

      # target 为常数时无法计算相关性
      if (length(unique(y_complete)) < 2) {
        return(
          tibble(
            feature = feature_name,
            target = target_col,
            group = group_label,
            rho = NA_real_,
            pvalue = NA_real_,
            n = n_complete,
            status = "constant_target"
          )
        )
      }

      # 执行 Spearman correlation
      # exact = FALSE 可避免 tied ranks 导致的 exact p-value 警告
      test_result <- tryCatch(
        suppressWarnings(
          cor.test(
            x = x_complete,
            y = y_complete,
            method = "spearman",
            exact = FALSE
          )
        ),
        error = function(e) e
      )

      # 捕获极少数计算错误，避免整个批量流程中断
      if (inherits(test_result, "error")) {
        return(
          tibble(
            feature = feature_name,
            target = target_col,
            group = group_label,
            rho = NA_real_,
            pvalue = NA_real_,
            n = n_complete,
            status = "calculation_error"
          )
        )
      }

      tibble(
        feature = feature_name,
        target = target_col,
        group = group_label,
        rho = unname(test_result$estimate),
        pvalue = test_result$p.value,
        n = n_complete,
        status = "OK"
      )
    }
  )

  result_df <- bind_rows(result_list)

  # ----------------------------------------------------------
  # 在当前 target × group 内进行 BH 多重检验校正
  # ----------------------------------------------------------

  result_df$padj <- NA_real_

  valid_idx <- which(
    result_df$status == "OK" &
      !is.na(result_df$pvalue)
  )

  if (length(valid_idx) > 0) {
    result_df$padj[valid_idx] <- p.adjust(
      result_df$pvalue[valid_idx],
      method = "BH"
    )
  }

  # ----------------------------------------------------------
  # 补充方向和绝对相关系数
  # ----------------------------------------------------------

  result_df <- result_df %>%
    mutate(
      abs_rho = abs(rho),
      direction = case_when(
        rho > 0 ~ "Positive",
        rho < 0 ~ "Negative",
        TRUE ~ NA_character_
      )
    ) %>%
    select(
      feature,
      target,
      group,
      rho,
      abs_rho,
      pvalue,
      padj,
      n,
      direction,
      status
    ) %>%
    arrange(pvalue)

  return(result_df)
}


# ============================================================
# 7. 批量运行相关性分析
# ============================================================

message(
  "本次实际分析 feature 数量：",
  nrow(feature_mat)
)

message(
  "本次实际分析 target 数量：",
  length(target_cols)
)

if (run_group_analysis) {
  message(
    "本次实际分析 group 数量：",
    length(selected_groups)
  )
}

all_result_list <- list()
result_index <- 1

# ------------------------------------------------------------
# 7.1 全样本分析
# ------------------------------------------------------------

if (run_all_samples) {

  message("开始运行全样本相关性分析。")

  for (target_name in target_cols) {

    result_one <- run_one_target_group(
      feature_mat = feature_mat,
      target_sub = target_df_aligned,
      target_col = target_name,
      sample_col = sample_col,
      group_label = "All",
      min_n = min_n
    )

    all_result_list[[result_index]] <- result_one
    result_index <- result_index + 1
  }
}

# ------------------------------------------------------------
# 7.2 组内分析
# ------------------------------------------------------------

if (run_group_analysis) {

  message("开始运行组内相关性分析。")

  for (group_name in selected_groups) {

    # 只保留当前组样本
    target_group <- target_df_aligned[
      target_df_aligned[[group_col]] == group_name,
      ,
      drop = FALSE
    ]

    message(
      "当前组：",
      group_name,
      "；原始样本数：",
      nrow(target_group)
    )

    # 当前组总样本数小于 min_n 时跳过
    if (nrow(target_group) < min_n) {
      warning(
        "组 ",
        group_name,
        " 的样本数小于 min_n，已跳过。"
      )
      next
    }

    for (target_name in target_cols) {

      result_one <- run_one_target_group(
        feature_mat = feature_mat,
        target_sub = target_group,
        target_col = target_name,
        sample_col = sample_col,
        group_label = group_name,
        min_n = min_n
      )

      all_result_list[[result_index]] <- result_one
      result_index <- result_index + 1
    }
  }
}

if (length(all_result_list) == 0) {
  stop(
    "没有产生任何结果。请检查 run_all_samples、",
    "run_group_analysis 和 selected_groups 参数。"
  )
}

cor_results <- bind_rows(all_result_list)


# ============================================================
# 8. 增加显著性筛选标记
# ============================================================

cor_results <- cor_results %>%
  mutate(
    # 探索性筛选：
    # 同时满足 |rho| >= r_cutoff 且 pvalue < p_cutoff
    significant_by_p = (
      status == "OK" &
        !is.na(rho) &
        !is.na(pvalue) &
        abs_rho >= r_cutoff &
        pvalue < p_cutoff
    ),

    # 更严格筛选：
    # 同时满足 |rho| >= r_cutoff 且 padj < padj_cutoff
    significant_by_padj = (
      status == "OK" &
        !is.na(rho) &
        !is.na(padj) &
        abs_rho >= r_cutoff &
        padj < padj_cutoff
    )
  )

# 探索性显著结果
cor_sig_p <- cor_results %>%
  filter(significant_by_p) %>%
  arrange(group, target, pvalue)

# FDR 显著结果
cor_sig_padj <- cor_results %>%
  filter(significant_by_padj) %>%
  arrange(group, target, padj)


# ============================================================
# 9. 导出结果表
# ============================================================

# ------------------------------------------------------------
# 9.1 导出全部结果
# ------------------------------------------------------------

write_csv(
  cor_results,
  file.path(table_dir, "Spearman_all_results.csv")
)

# ------------------------------------------------------------
# 9.2 导出按原始 p 值筛选的结果
# ------------------------------------------------------------

write_csv(
  cor_sig_p,
  file.path(
    table_dir,
    paste0(
      "Spearman_sig_absR",
      r_cutoff,
      "_p",
      p_cutoff,
      ".csv"
    )
  )
)

# ------------------------------------------------------------
# 9.3 导出按 padj 筛选的结果
# ------------------------------------------------------------

write_csv(
  cor_sig_padj,
  file.path(
    table_dir,
    paste0(
      "Spearman_sig_absR",
      r_cutoff,
      "_padj",
      padj_cutoff,
      ".csv"
    )
  )
)

# ------------------------------------------------------------
# 9.4 按 target × group 分开导出
# ------------------------------------------------------------

safe_filename <- function(x) {
  gsub(
    pattern = "[^A-Za-z0-9._-]+",
    replacement = "_",
    x = x
  )
}

result_split <- split(
  cor_results,
  interaction(
    cor_results$group,
    cor_results$target,
    drop = TRUE
  )
)

walk(
  result_split,
  function(df_one) {

    group_name  <- unique(df_one$group)
    target_name <- unique(df_one$target)

    filename <- paste0(
      safe_filename(group_name),
      "__",
      safe_filename(target_name),
      "__all_results.csv"
    )

    write_csv(
      df_one,
      file.path(table_dir, filename)
    )
  }
)

# ------------------------------------------------------------
# 9.5 输出简单统计汇总
# ------------------------------------------------------------

result_summary <- cor_results %>%
  group_by(group, target) %>%
  summarise(
    total_features = n(),
    valid_results = sum(status == "OK"),
    significant_by_p = sum(significant_by_p, na.rm = TRUE),
    significant_by_padj = sum(significant_by_padj, na.rm = TRUE),
    positive_by_p = sum(
      significant_by_p & direction == "Positive",
      na.rm = TRUE
    ),
    negative_by_p = sum(
      significant_by_p & direction == "Negative",
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(result_summary)

write_csv(
  result_summary,
  file.path(table_dir, "Spearman_result_summary.csv")
)


# ============================================================
# 10. 可选：相关性总览图
# ============================================================

plot_correlation_overview <- function(
    result_df,
    target_name,
    group_name,
    r_cutoff = 0.5,
    p_cutoff = 0.05,
    label_top_n_each = 0
) {

  plot_df <- result_df %>%
    filter(
      target == target_name,
      group == group_name,
      status == "OK",
      !is.na(rho),
      !is.na(pvalue)
    ) %>%
    mutate(
      # pvalue 极接近 0 时避免 -log10(0) 得到 Inf
      neg_log10_p = -log10(
        pmax(
          pvalue,
          .Machine$double.xmin
        )
      ),
      correlation_class = case_when(
        rho >= r_cutoff &
          pvalue < p_cutoff ~ "Positive",
        rho <= -r_cutoff &
          pvalue < p_cutoff ~ "Negative",
        TRUE ~ "Not significant"
      )
    )

  p <- ggplot(
    plot_df,
    aes(
      x = rho,
      y = neg_log10_p,
      color = correlation_class
    )
  ) +
    geom_point(
      size = 1.7,
      alpha = 0.8
    ) +
    geom_vline(
      xintercept = c(-r_cutoff, r_cutoff),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed"
    ) +
    scale_color_manual(
      values = c(
        "Negative" = "#3B6FB6",
        "Not significant" = "grey75",
        "Positive" = "#D94B3D"
      )
    ) +
    labs(
      title = paste0(
        group_name,
        " within-group correlation overview"
      ),
      subtitle = paste0(
        "Target: ",
        target_name
      ),
      x = "Spearman correlation coefficient",
      y = "-log10(p value)",
      color = "Correlation"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      ),
      panel.grid.minor = element_blank()
    )

  # 仅在 label_top_n_each > 0 时标注基因
  if (label_top_n_each > 0) {

    label_df <- bind_rows(
      plot_df %>%
        filter(correlation_class == "Positive") %>%
        arrange(pvalue, desc(rho)) %>%
        slice_head(n = label_top_n_each),

      plot_df %>%
        filter(correlation_class == "Negative") %>%
        arrange(pvalue, rho) %>%
        slice_head(n = label_top_n_each)
    )

    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(label = feature),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.3,
        point.padding = 0.2,
        show.legend = FALSE
      )
  }

  return(p)
}

if (make_overview_plot) {

  plot_plan <- cor_results %>%
    distinct(group, target)

  for (i in seq_len(nrow(plot_plan))) {

    group_name  <- plot_plan$group[i]
    target_name <- plot_plan$target[i]

    p_overview <- plot_correlation_overview(
      result_df = cor_results,
      target_name = target_name,
      group_name = group_name,
      r_cutoff = r_cutoff,
      p_cutoff = p_cutoff,
      label_top_n_each = label_top_n_each
    )

    print(p_overview)

    file_prefix <- paste0(
      safe_filename(group_name),
      "__",
      safe_filename(target_name),
      "__correlation_overview"
    )

    ggsave(
      filename = file.path(
        plot_dir,
        paste0(file_prefix, ".png")
      ),
      plot = p_overview,
      width = 7,
      height = 5.5,
      dpi = 300
    )

    ggsave(
      filename = file.path(
        plot_dir,
        paste0(file_prefix, ".pdf")
      ),
      plot = p_overview,
      width = 7,
      height = 5.5
    )
  }
}


# ============================================================
# 11. 可选：单个 feature 的相关性散点图
# ============================================================

plot_feature_scatter <- function(
    feature_name,
    target_name,
    group_name = "All",
    feature_mat,
    target_df_aligned,
    sample_col,
    group_col = NULL
) {

  # 检查 feature 是否存在
  if (!feature_name %in% rownames(feature_mat)) {
    stop(
      "feature_mat 中不存在 feature：",
      feature_name
    )
  }

  # 检查 target 是否存在
  if (!target_name %in% colnames(target_df_aligned)) {
    stop(
      "target_df_aligned 中不存在 target：",
      target_name
    )
  }

  # 根据 group_name 选择样本
  if (group_name == "All") {

    target_sub <- target_df_aligned

  } else {

    if (is.null(group_col)) {
      stop("组内画图时必须提供 group_col。")
    }

    target_sub <- target_df_aligned[
      target_df_aligned[[group_col]] == group_name,
      ,
      drop = FALSE
    ]
  }

  samples_use <- target_sub[[sample_col]]

  scatter_df <- tibble(
    sample = samples_use,
    feature_value = as.numeric(
      feature_mat[
        feature_name,
        samples_use
      ]
    ),
    target_value = target_sub[[target_name]]
  ) %>%
    filter(
      complete.cases(
        feature_value,
        target_value
      )
    )

  if (nrow(scatter_df) < min_n) {
    stop("完整样本数小于 min_n，无法绘图。")
  }

  cor_result <- cor.test(
    x = scatter_df$feature_value,
    y = scatter_df$target_value,
    method = "spearman",
    exact = FALSE
  )

  rho_value <- unname(cor_result$estimate)
  p_value   <- cor_result$p.value

  p <- ggplot(
    scatter_df,
    aes(
      x = target_value,
      y = feature_value
    )
  ) +
    geom_point(
      size = 2.5,
      alpha = 0.85
    ) +
    # 该直线用于展示样本点的总体趋势
    # 图中报告的统计结果仍然是 Spearman rho
    geom_smooth(
      method = "lm",
      se = TRUE
    ) +
    annotate(
      geom = "text",
      x = Inf,
      y = Inf,
      hjust = 1.1,
      vjust = 1.3,
      label = paste0(
        "Spearman rho = ",
        round(rho_value, 3),
        "\np = ",
        signif(p_value, 3),
        "\nn = ",
        nrow(scatter_df)
      )
    ) +
    labs(
      title = paste0(
        feature_name,
        " vs ",
        target_name
      ),
      subtitle = paste0(
        "Group: ",
        group_name
      ),
      x = target_name,
      y = paste0(
        feature_name,
        " value"
      )
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      )
    )

  return(p)
}


# ============================================================
# 12. 单 feature 散点图使用示例
# ============================================================

# 【按需取消注释并替换参数】
#
# p_gene <- plot_feature_scatter(
#   feature_name = "RIPOR2",
#   target_name = "T_Cell_Exhaustion_Progenitor_like",
#   group_name = "IA",
#   feature_mat = feature_mat,
#   target_df_aligned = target_df_aligned,
#   sample_col = sample_col,
#   group_col = group_col
# )
#
# print(p_gene)
#
# ggsave(
#   filename = file.path(
#     plot_dir,
#     "IA__RIPOR2__scatter.png"
#   ),
#   plot = p_gene,
#   width = 5,
#   height = 4,
#   dpi = 300
# )


# ============================================================
# 13. 最终对象说明
# ============================================================

# feature_df_selected：
# 根据 selected_features 选择后的 feature 数据框
#
# feature_df_aligned：
# feature 已选择，且样本已经匹配并统一顺序的数据框
#
# target_df_aligned：
# 样本已经匹配并统一顺序的 target 数据框
#
# feature_mat：
# 用于内部批量计算的 numeric matrix
#
# cor_results：
# 全部相关性结果
#
# cor_sig_p：
# 满足 |rho| >= r_cutoff 且 pvalue < p_cutoff 的结果
#
# cor_sig_padj：
# 满足 |rho| >= r_cutoff 且 padj < padj_cutoff 的结果
#
# result_summary：
# 每个 target × group 的结果数量汇总


# ============================================================
# 14. 保存 R 对象
# ============================================================

save(
  feature_df_selected,
  feature_df_aligned,
  target_df_aligned,
  feature_mat,
  target_cols,
  selected_groups,
  cor_results,
  cor_sig_p,
  cor_sig_padj,
  result_summary,
  file = file.path(
    outdir,
    "Spearman_workflow_results.RData"
  )
)

message("Spearman correlation workflow 运行完成。")
message("结果目录：", normalizePath(outdir))
