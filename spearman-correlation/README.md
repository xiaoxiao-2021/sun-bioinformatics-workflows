# Spearman Correlation Workflow

A reusable R workflow for batch Spearman rank correlation analysis between multiple features and one or more target variables.

该工作流用于批量分析多个 `feature` 与一个或多个 `target` 之间的 Spearman 秩相关关系，支持全样本分析、组内分析、样本匹配、缺失值处理、多重检验校正、结果筛选和基础可视化。


本流程中路径可根据自身实际存储路径进行改动！！ 流程中所显示的路径仅为试验时所用！！ お疲れ様(●ˇ∀ˇ●)
---

## Repository structure

```text
sun-bioinformatics-workflows/
└── spearman-correlation/
    ├── README.md
    ├── run_spearman_analysis.R
    ├── R/
    │   └── run_spearman_workflow.R
    ├── examples/
    │   ├── feature_df_example.csv
    │   ├── target_df_example.csv
    │   └── run_function_example.R
    └── figures/
```

### Function file

```text
spearman-correlation/R/run_spearman_workflow.R
```

该文件只定义可复用函数：

```r
run_spearman_workflow()
plot_spearman_overview()
plot_spearman_feature()
```

执行：

```r
source(
  "spearman-correlation/R/run_spearman_workflow.R"
)
```

只会加载函数，不会自动读取数据或开始分析。

### Analysis launcher

```text
spearman-correlation/run_spearman_analysis.R
```

该文件是实际分析入口，负责：

```text
读取输入数据
→ 加载函数
→ 设置分析参数
→ 调用 workflow
→ 查看与验证结果
```

分析新数据时，主要修改启动脚本，不需要反复修改函数文件。

---

## Required R packages

核心统计分析使用 base R。

绘图需要：

```r
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
```

在总览图中自动标注 feature 时还需要：

```r
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
```

---

## Input data

Workflow 使用两个标准化 `data.frame`：

```r
feature_df
target_df
```

### `feature_df`

要求：

- 每一行代表一个 feature
- 一列保存 feature 名称
- 其余列代表样本
- 样本列必须为 numeric
- feature 名称不能重复
- 样本列名不能重复

示例：

| feature | S01 | S02 | S03 | S04 |
|---|---:|---:|---:|---:|
| Gene_A | 4.2 | 5.1 | 3.8 | 4.6 |
| Gene_B | 2.3 | 2.8 | 4.1 | 3.7 |
| Gene_C | 7.4 | 6.9 | 8.0 | 7.2 |

如果原始数据是 numeric matrix，并且 feature 名称保存在行名中：

```r
feature_df <- data.frame(
  feature = rownames(feature_mat),
  feature_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
```

### `target_df`

要求：

- 每一行代表一个样本
- 一列保存样本名称
- 一个或多个 numeric 列保存 target
- 可选分组列用于组内分析
- 每个样本只能出现一次

示例：

| Sample | Score_A | Group |
|---|---:|---|
| S01 | 0.42 | IA |
| S02 | 0.81 | IA |
| S03 | 1.10 | ENEG |
| S04 | 0.35 | ENEG |

一次明确的评分分析建议直接指定单个 target：

```r
target_cols <- "Score_A"
```

`target_cols <- NULL` 会自动选择除样本列和分组列之外的全部 numeric 列。

---

# Quick start

## 1. 设置工作目录

建议将工作目录设置为整个 GitHub 仓库根目录：

```r
setwd(
  "D:/Git/sun-bioinformatics-workflows"
)
```

---

## 2. 打开启动脚本

打开：

```text
spearman-correlation/run_spearman_analysis.R
```

每次分析主要修改该文件，而不是直接修改：

```text
spearman-correlation/R/run_spearman_workflow.R
```

---

## 3. 修改输入文件路径

模拟数据：

```r
feature_file <- paste0(
  "spearman-correlation/examples/",
  "feature_df_example.csv"
)

target_file <- paste0(
  "spearman-correlation/examples/",
  "target_df_example.csv"
)
```

真实数据示例：

```r
feature_file <- "D:/project_A/data/feature_df.csv"
target_file  <- "D:/project_A/data/target_df.csv"
```

读取数据：

```r
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
```

---

## 4. 设置分析参数

### Input columns

```r
feature_col <- "feature"
sample_col <- "Sample"
```

列名区分大小写。

### Feature selection

分析全部 feature：

```r
selected_features <- NULL
```

只分析指定 feature：

```r
selected_features <- c(
  "RIPOR2",
  "TOX",
  "TCF7"
)
```

### Target selection

分析一个指定评分：

```r
target_cols <- "Score_A"
```

### Whole-cohort and group-wise analysis

```r
run_all_samples <- TRUE
run_group_analysis <- TRUE
```

- `run_all_samples = TRUE`：分析全部共同样本，结果组名为 `All`
- `run_group_analysis = TRUE`：分别进行组内分析

### Group settings

```r
group_col <- "Group"
selected_groups <- NULL
```

当 `selected_groups <- NULL` 时，会自动分析全部非缺失、非空分组。

例如真实分组为 `IA` 和 `ENEG`，且两个分析开关均为 `TRUE`，最终运行：

```text
All
IA
ENEG
```

### Statistical thresholds

```r
min_n <- 5
r_cutoff <- 0.5
p_cutoff <- 0.05
padj_cutoff <- 0.05
```

参数含义：

- `min_n`：每次相关性计算所需的最少完整样本数
- `r_cutoff`：相关系数绝对值阈值
- `p_cutoff`：原始 p 值阈值，用于探索性候选筛选
- `padj_cutoff`：BH 校正后 p 值阈值，用于 FDR 控制后的高置信筛选

### Significance levels and intended use

Workflow 同时保留原始 p 值和 BH 校正后的 p 值，两类结果的证据等级和用途不同。

```text
探索性候选：
|rho| >= r_cutoff 且 pvalue < p_cutoff

FDR 控制后的高置信结果：
|rho| >= r_cutoff 且 padj < padj_cutoff
```

- `pvalue < p_cutoff` 表示达到名义显著性，适合候选 feature 发现、与其他证据取交集、富集分析和后续验证。
- `padj < padj_cutoff` 表示在同一 `target × group` 分析范围内通过 BH 多重检验校正，可作为更严格的高置信结果。
- `r_cutoff` 用于限制效应量，避免仅凭较小 p 值筛选相关强度很弱的结果。
- 当组内样本量较少而 feature 数量较多时，`padj < 0.05` 可能非常严格，甚至没有结果通过筛选。
- 未通过 `padj` 阈值不等同于证明不存在相关性，只表示当前数据不足以支持 FDR 控制后的显著结论。

建议始终同时保留并报告：

```text
rho
pvalue
padj
n
```

### Sample matching

```r
strict_sample_match <- FALSE
```

- `FALSE`：报告不匹配样本，并使用共同样本继续分析
- `TRUE`：发现任何不匹配样本时停止

### Output settings

测试时：

```r
save_results <- FALSE
```

结果只保存在当前 R 会话中的 `spearman_result`，不会写入硬盘。

正式分析时：

```r
save_results <- TRUE
outdir <- "D:/project_A/results/spearman_analysis"
```

建议将正式结果保存到代码仓库外部。

### Plot settings

```r
make_overview_plot <- TRUE
label_top_n_each <- 0
verbose <- TRUE
```

---

## 5. Start the analysis

启动脚本先加载函数：

```r
source(
  "spearman-correlation/R/run_spearman_workflow.R"
)
```

真正启动分析的是：

```r
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
```

---

# Sample matching

函数会：

```text
提取两个输入中的样本名
→ 使用 setdiff() 报告不匹配样本
→ 使用 intersect() 提取共同样本
→ 使用 match() 统一 target_df 顺序
→ 检查最终样本顺序
```

运行后检查：

```r
spearman_result$sample_report
```

其中：

```r
spearman_result$sample_report$sample_order_ok
```

正常应为：

```r
TRUE
```

样本顺序必须一致，因为 `cor.test()` 按向量位置计算，不会自动根据样本名匹配。

---

# Statistical analysis

对每个 `feature × target × group`，函数运行：

```r
cor.test(
  x = x_complete,
  y = y_complete,
  method = "spearman",
  exact = FALSE
)
```

BH 多重检验校正在每个独立的：

```text
target × group
```

内部进行。

---

# Returned result object

主函数返回：

```r
spearman_result
```

它是一个 `spearman_workflow_result` list，包含：

```text
call
parameters
sample_report
features_used
targets_used
groups_used
feature_df_aligned
target_df_aligned
feature_mat
results
significant_p
significant_padj
summary
overview_plots
```

## `spearman_result$summary`

每个 `target × group` 的结果摘要。

主要字段：

| Field | Description |
|---|---|
| `group` | `All` 或具体分组 |
| `target` | target 名称 |
| `total_features` | 分析的 feature 总数 |
| `valid_results` | 正常完成计算的结果数 |
| `significant_by_p` | 同时满足相关系数阈值与原始 p 值阈值的探索性候选数量 |
| `significant_by_padj` | 同时满足相关系数阈值与 BH 校正 p 值阈值的高置信结果数量 |
| `positive_by_p` | 满足效应量与原始 p 值阈值的正相关候选数量 |
| `negative_by_p` | 满足效应量与原始 p 值阈值的负相关候选数量 |

## `spearman_result$results`

完整结果表。

| Field | Description |
|---|---|
| `feature` | Feature name |
| `target` | Target variable |
| `group` | `All` 或具体分组 |
| `rho` | Spearman correlation coefficient |
| `abs_rho` | `rho` 的绝对值 |
| `pvalue` | 原始 p 值 |
| `padj` | BH 校正后的 p 值 |
| `n` | 实际参与计算的完整样本数 |
| `direction` | `Positive`、`Negative` 或 `Zero` |
| `status` | 计算状态 |
| `significant_by_p` | 是否同时满足相关系数阈值与原始 p 值阈值，用于探索性候选筛选 |
| `significant_by_padj` | 是否同时满足相关系数阈值与 BH 校正 p 值阈值，用于 FDR 控制后的高置信筛选 |

可能的 `status`：

| Status | Meaning |
|---|---|
| `OK` | 正常完成计算 |
| `insufficient_n` | 有效样本数小于 `min_n` |
| `constant_feature` | feature 没有变化 |
| `constant_target` | target 没有变化 |
| `calculation_error` | 计算失败 |

## Candidate and FDR-controlled results

### Exploratory candidates based on nominal p value

```r
spearman_result$significant_p
```

筛选条件：

```text
|rho| >= r_cutoff
pvalue < p_cutoff
```

该结果用于探索性候选 feature 筛选。由于相关性检验通常同时覆盖大量 feature，该结果未控制全体检验的假发现率，不能直接表述为“BH 校正后显著相关结果”。

推荐表述：

```text
达到预设相关系数阈值且名义 p 值小于 0.05 的探索性相关候选
```

### FDR-controlled results

```r
spearman_result$significant_padj
```

筛选条件：

```text
|rho| >= r_cutoff
padj < padj_cutoff
```

该结果通过 BH 多重检验校正，表示在当前 `target × group` 分析范围内达到设定 FDR 阈值的高置信相关结果。

推荐表述：

```text
经过 BH 多重检验校正后显著相关的 feature
```

### Recommended use of both result levels

对于小样本组内分析，推荐同时保留两层结果：

```text
significant_p
→ 用于探索性候选发现、与差异表达结果取交集、富集分析和后续验证

significant_padj
→ 用于判断是否存在通过全 feature 多重检验校正的高置信结果
```

与差异表达基因取交集可以增加候选结果的生物学针对性，但不会使未经校正的相关性 p 值自动成为 FDR 显著结果。

---

# Viewing results

查看摘要：

```r
spearman_result$summary
```

查看完整结果：

```r
head(
  spearman_result$results
)
```

查看基于原始 p 值的探索性候选：

```r
spearman_result$significant_p
```

查看 BH 校正后的高置信结果：

```r
spearman_result$significant_padj
```

查看某个 feature：

```r
subset(
  spearman_result$results,
  feature == "RIPOR2"
)
```

查看某个组：

```r
subset(
  spearman_result$results,
  group == "IA"
)
```

---

# Correlation overview plots

查看图名：

```r
names(
  spearman_result$overview_plots
)
```

显示第一张图：

```r
print(
  spearman_result$overview_plots[[1]]
)
```

每个点代表一个 feature：

- x 轴：Spearman `rho`
- y 轴：`-log10(p value)`
- 虚线：rho 和原始 p 值筛选阈值

总览图主要用于观察相关方向、相关强度和名义 p 值分布。图中越过原始 p 值虚线的点仍属于探索性候选，是否达到 FDR 控制后的显著性应查看结果表中的 `padj` 或 `significant_by_padj`。

---

# Single-feature scatter plot

使用：

```r
plot_spearman_feature()
```

连续型 target 的散点图示例：

```r
p_gene <- plot_spearman_feature(
  workflow_result = spearman_result,
  feature_name = "Gene_Pos",
  target_name = "Score_A",
  group_name = "All",
  plot_type = "scatter",
  add_lm = TRUE
)

print(p_gene)
```

有序等级 target 的箱线图加抖动点示例：

```r
p_ordinal <- plot_spearman_feature(
  workflow_result = spearman_result,
  feature_name = "Gene_Pos",
  target_name = "Scheuer_S",
  group_name = "All",
  plot_type = "ordinal",
  add_lm = FALSE
)

print(p_ordinal)
```

每个点代表一个样本。图中显示：

- Spearman `rho`
- 原始 p 值
- 有效样本数

`plot_type = "scatter"` 适合连续型 target；`plot_type = "ordinal"` 适合 0、1、2、3、4 等有序等级 target。

`add_lm = TRUE` 只在线性散点图中添加趋势线用于可视化，不改变 Spearman 统计结果；在 `ordinal` 模式中不会用于统计计算。

---

# Output files

当：

```r
save_results <- TRUE
```

时，结果保存到 `outdir`：

```text
outdir/
├── aligned_inputs/
│   ├── feature_df_aligned.csv
│   ├── target_df_aligned.csv
│   ├── features_used.csv
│   ├── targets_used.csv
│   └── groups_used.csv
├── plots/
│   ├── All__Score_A__correlation_overview.png
│   ├── All__Score_A__correlation_overview.pdf
│   └── ...
├── tables/
│   ├── Spearman_all_results.csv
│   ├── Spearman_sig_absR0.5_p0.05.csv
│   ├── Spearman_sig_absR0.5_padj0.05.csv
│   ├── Spearman_result_summary.csv
│   └── group__target__all_results.csv
└── Spearman_workflow_results.rds
```

其中：

```text
Spearman_sig_absR0.5_p0.05.csv
```

保存基于原始 p 值的探索性候选；

```text
Spearman_sig_absR0.5_padj0.05.csv
```

保存经过 BH 多重检验校正后的高置信结果。

重新读取完整结果对象：

```r
saved_result <- readRDS(
  "D:/project_A/results/spearman_analysis/Spearman_workflow_results.rds"
)
```

---

# Recommended validation

运行前：

```r
str(feature_df)
str(target_df)
```

运行后：

```r
spearman_result$sample_report$sample_order_ok

table(
  spearman_result$results$status,
  useNA = "ifany"
)
```

手工验证一个结果：

```r
samples_check <- spearman_result$target_df_aligned[
  [sample_col]
]

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
```

与结果表对应行比较：

```r
subset(
  spearman_result$results,
  feature == "Gene_Pos" &
    target == "Score_A" &
    group == "All"
)
```

---

# Interpretation notes

结果解释应综合考虑：

- 相关方向
- 相关强度
- 原始 p 值
- BH 校正后的 p 值
- 有效样本数
- 散点或有序等级分布
- 生物学合理性
- 独立数据或实验验证

## Exploratory and confirmatory interpretation

本工作流同时输出名义 p 值结果和 BH 校正后的结果。

- `pvalue < 0.05` 且达到相关系数阈值的结果，应表述为探索性相关候选。
- `padj < 0.05` 且达到相关系数阈值的结果，可表述为经过 BH 多重检验校正后的显著相关结果。
- 当组内样本量较少、feature 数量较多时，BH 校正后的筛选可能非常严格。
- 未通过 `padj` 阈值不等同于证明不存在相关性，而表示当前数据不足以支持 FDR 控制后的显著结论。
- 探索性结果适合用于候选基因筛选、与差异表达结果取交集、功能富集和后续验证，但应同时报告 `rho`、`pvalue`、`padj` 和有效样本数 `n`。
- 与差异表达基因取交集可以增加候选结果的生物学针对性，但不会使未经校正的相关性自动成为 FDR 显著结果。

例如：

```text
上调差异基因
∩
组内正相关探索性候选
```

更准确的解释是：

```text
同时具有组间差异表达证据和组内探索性相关证据的候选 feature
```

如果相关性结果同时满足 `padj < padj_cutoff`，才可以进一步描述为同时具备差异表达证据与 FDR 校正后相关证据。

Correlation does not establish causation.

---

# Common issues

## Unmatched sample names

常见原因：

- 前后空格
- 字母大小写不同
- 分隔符不同
- R 自动修改列名
- 某个输入缺少样本

## Non-numeric columns

常见原因：

- `"-"`
- `"unknown"`
- `"not detected"`
- 空字符串
- 数值中混入单位

## Duplicate samples or features

- `target_df` 中每个样本只能出现一次
- `feature_df` 中每个 feature 名称只能出现一次
- feature 样本列名不能重复

## Insufficient sample size

某个分组样本数小于 `min_n` 时，该组会被跳过。

某个 feature 因缺失值导致有效样本数不足时，状态为：

```text
insufficient_n
```

## Constant variables

对应状态：

```text
constant_feature
constant_target
```

## Tied ranks

Workflow 使用：

```r
exact = FALSE
```

计算近似 p 值。

---

# Recommended project use

建议代码仓库与真实分析项目分开管理：

```text
GitHub repository
├── reusable function
├── analysis launcher template
├── example data
└── documentation

Local analysis project
├── real data
├── copied analysis launcher
└── result files
```

实际分析时，可以将 `run_spearman_analysis.R` 复制到本地项目目录，再修改：

```text
feature_file
target_file
target_cols
group_col
selected_groups
save_results
outdir
```

---

# Detailed documentation

[Spearman correlation workflow](https://xiaoxiao-2021.github.io/sun-lab-life-astro-v2/bioinformatics/spearman-correlation-workflow/)

---

# Development note

This workflow was developed with AI-assisted code drafting and was reviewed, tested, and maintained by the repository author.

The current function version has passed a minimal simulated-data test covering:

- Whole-cohort analysis
- Within-group analysis
- Constant-feature detection
- Result summary generation
- Overview plot generation
- Single-feature scatter plotting
- Manual `cor.test()` verification

---

# License

This project is released under the MIT License.
