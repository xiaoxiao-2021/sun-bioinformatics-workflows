# Spearman Correlation Workflow

这是一个配置驱动、可重复运行、适合命令行和服务器执行的 Spearman 秩相关分析 workflow。它用于批量计算多个 feature 与一个或多个 target 之间的相关性，支持全样本分析、组内分析、候选结果筛选、方向拆分、可选四象限后处理、绘图、手工复核和 RDS 断点复用。

`run_spearman_analysis.R` 是唯一正式命令行入口。workflow 不清空工作区、不依赖 RStudio Global Environment、不调用 `setwd()`，也不要求终端当前目录固定。旧的 `run_spearman_workflow()` 函数接口仅为已有 R 代码保留兼容性；新项目应使用 YAML 配置和正式入口。

## 统计语义

对每个 `feature × target × group` 组合，workflow 在完整且有限的成对观测上运行：

```r
cor.test(x, y, method = "spearman", exact = FALSE)
```

- `rho` 是 Spearman 秩相关系数。
- `pvalue` 使用 `exact = FALSE` 得到近似 P 值，这也避免 tied ranks 导致 exact test 不可用。
- `padj` 在每个独立的 `target × group` 分析范围内，对该范围的有效 feature P 值执行 Benjamini–Hochberg 校正。
- 探索性候选的边界是 `abs(rho) >= r_cutoff` 且 `pvalue < p_cutoff`。
- FDR 筛选的边界是 `abs(rho) >= r_cutoff` 且 `padj < padj_cutoff`。
- `filters.include_positive` 和 `filters.include_negative` 只控制是否额外导出相应方向的筛选表，不改变完整结果、P 值、BH 校正或主筛选结果。

相关不等于因果。解释结果时应同时报告 `rho`、`pvalue`、`padj`、有效样本数 `n`、分组和数据质量。

## 输入数据

两个输入都支持 `csv`、`tsv` 或 `rds`。RDS 文件中必须保存一个 `data.frame`。

### 表达矩阵

宽格式，每行一个 feature：

| feature | S01 | S02 | S03 |
|---|---:|---:|---:|
| Gene_A | 4.2 | 5.1 | 3.8 |
| Gene_B | 2.3 | 2.8 | 4.1 |

要求：

- `columns.feature` 指定 feature 名称列。
- 其余被用作样本的列名是样本 ID，必须唯一。
- feature 名称必须唯一、非缺失且非空。
- 样本数据必须是真正的数值型数据；字符数字、单位、`unknown` 或 `-` 不会被静默当作数值。

### Score/metadata 表

长格式，每行一个样本：

| Sample | Score_A | Group |
|---|---:|---|
| S01 | 0.42 | Group_A |
| S02 | 0.81 | Group_A |
| S03 | 1.10 | Group_B |

要求：

- `columns.sample` 指定样本 ID 列，ID 必须唯一、非缺失且非空。
- `analysis.target_columns` 指定一个或多个数值型 target 列。
- `columns.group` 指定可选分组列；启用组内分析时该列必须存在。
- 表达矩阵和 metadata 的样本按名称匹配并显式重排，不依赖原始列/行顺序。

`analysis.strict_sample_match: false` 时会报告不匹配样本并仅使用交集；设为 `true` 时任何不匹配都会终止运行。

## 环境

不要把依赖安装到 `base` 环境。`environment.yml` 定义的环境名是 `spearman-r`。

首次创建：

```bash
cd "$WORKFLOW_REPO/spearman-correlation"
micromamba create -f environment.yml
micromamba activate spearman-r
```

环境已存在时更新：

```bash
cd "$WORKFLOW_REPO/spearman-correlation"
micromamba update -n spearman-r -f environment.yml
micromamba activate spearman-r
```

主要依赖包括 `r-base`、`r-dplyr`、`r-tibble`、`r-purrr`、`r-readr`、`r-ggplot2`、`r-ggrepel`、`r-yaml` 和 `r-optparse`。

## 快速开始

直接运行仓库示例：

```bash
micromamba activate spearman-r
Rscript "$WORKFLOW_REPO/spearman-correlation/run_spearman_analysis.R" \
  --config "$WORKFLOW_REPO/spearman-correlation/config/config.example.yml"
```

推荐把真实分析配置和数据放在仓库外：

```bash
mkdir -p "$ANALYSIS_DIR"/{data,results}
cp "$WORKFLOW_REPO/spearman-correlation/config/config.example.yml" \
  "$ANALYSIS_DIR/config.yml"
```

编辑 `$ANALYSIS_DIR/config.yml` 后运行：

```bash
micromamba activate spearman-r
Rscript "$WORKFLOW_REPO/spearman-correlation/run_spearman_analysis.R" \
  --config "$ANALYSIS_DIR/config.yml"
```

在 WSL/Linux 中，先把变量指向实际 Linux/WSL 路径，例如：

```bash
export WORKFLOW_REPO=/path/to/sun-bioinformatics-workflows
export ANALYSIS_DIR=/path/to/analysis-project
```

不要把 Windows 盘符路径直接写进面向 Linux 的共享配置；WSL 下应使用 `/mnt/...` 形式或纯 Linux 路径。

## 路径规则与推荐目录

YAML 中所有相对输入路径和 `output_dir` 都以该 YAML 文件所在目录为基准解析，而不是以调用 `Rscript` 时的当前目录为基准。绝对路径保持为绝对路径。因此，同一个配置可以从任意终端目录重复执行。

推荐结构：

```text
$WORKFLOW_REPO/
└── spearman-correlation/
    ├── config/config.example.yml
    ├── examples/
    ├── R/
    ├── environment.yml
    └── run_spearman_analysis.R

$ANALYSIS_DIR/
├── config.yml
├── data/
│   ├── expression.csv
│   └── metadata.csv
└── results/                 # workflow 自动创建
```

若 `$ANALYSIS_DIR/config.yml` 使用上述结构，可写：

```yaml
inputs:
  expression_matrix: data/expression.csv
  score_metadata: data/metadata.csv
output_dir: results/spearman
```

## 完整配置示例

下面展示当前 schema 的全部字段。`null` 可表示不限制 feature/group；`analysis.target_columns: null` 会按 metadata 原顺序选择除 sample/group 外的全部数值列。启用依赖具体项目名称的模块前，应替换其中的 feature、target 和 group。

```yaml
project_id: spearman_example

inputs:
  expression_matrix: ../examples/feature_df_example.csv
  score_metadata: ../examples/target_df_example.csv
  expression_format: csv
  metadata_format: csv

output_dir: ../results/spearman_example

columns:
  feature: feature
  sample: Sample
  group: Group

analysis:
  selected_features: null
  target_columns:
    - Score_A
  run_all_samples: true
  run_group_analysis: true
  selected_groups: null
  correlation_method: spearman
  min_n: 5
  strict_sample_match: false

thresholds:
  r_cutoff: 0.5
  p_cutoff: 0.05
  padj_cutoff: 0.05

filters:
  include_positive: true
  include_negative: true

quadrant:
  enabled: false
  x_target: null
  x_group: null
  y_target: null
  y_group: null
  significance: padj
  require_significant_both: true

outputs:
  aligned_inputs: true
  combined_tables: true
  per_analysis_tables: true
  filtered_tables: true
  result_rds: true

intermediate:
  enabled: true
  reuse_existing: true

plots:
  overview:
    enabled: true
    label_top_n_each: 0
    width: 7
    height: 5.5
    dpi: 300
    formats:
      - png
      - pdf
  single_feature:
    enabled: true
    width: 5
    height: 4
    dpi: 300
    formats:
      - png
      - pdf
    items:
      - feature: Gene_Pos
        target: Score_A
        group: All
        type: scatter
        add_lm: true
  quadrant:
    enabled: false
    width: 7
    height: 5.5
    dpi: 300
    formats:
      - png
      - pdf

validation:
  manual:
    enabled: true
    items:
      - feature: Gene_Pos
        target: Score_A
        group: All

logging:
  verbose: true
```

## 配置字段

配置会递归拒绝未知字段，拼写错误不会被静默忽略。

| 字段 | 类型/允许值 | 作用 |
|---|---|---|
| `project_id` | 非空字符串 | 项目标识，也用于日志文件名；不要写个人信息。 |
| `inputs.expression_matrix` | 路径 | 表达矩阵文件。相对路径以配置目录为基准。 |
| `inputs.score_metadata` | 路径 | score/metadata 文件。相对路径以配置目录为基准。 |
| `inputs.expression_format` | `csv` / `tsv` / `rds` | 表达矩阵格式。 |
| `inputs.metadata_format` | `csv` / `tsv` / `rds` | metadata 格式。 |
| `output_dir` | 路径 | 结果根目录；相对路径以配置目录为基准。 |
| `columns.feature` | 非空字符串 | 表达矩阵中的 feature 名称列。 |
| `columns.sample` | 非空字符串 | metadata 中的样本 ID 列。 |
| `columns.group` | 字符串或 `null` | 分组列；启用组内分析时必须配置并存在。 |
| `analysis.selected_features` | 字符串列表或 `null` | 只分析指定 feature；`null` 表示全部。 |
| `analysis.target_columns` | 字符串列表或 `null` | 要分析的数值型 target 列；`null` 自动选择除 sample/group 外的全部数值列。 |
| `analysis.run_all_samples` | 布尔值 | 是否生成组名为 `All` 的全样本分析。 |
| `analysis.run_group_analysis` | 布尔值 | 是否分别运行组内分析。 |
| `analysis.selected_groups` | 字符串列表或 `null` | 指定组；`null` 表示所有有效非空组。 |
| `analysis.correlation_method` | `spearman` | 相关方法。当前最终版只允许 Spearman，以保留统计定义。 |
| `analysis.min_n` | 整数，至少 3 | 每次相关计算所需的最少完整样本数。 |
| `analysis.strict_sample_match` | 布尔值 | 是否要求两个输入中的样本集合完全一致。 |
| `thresholds.r_cutoff` | 数值，`0` 到 `1` | 相关强度筛选阈值，使用 `>=`。 |
| `thresholds.p_cutoff` | 数值，`0 < p <= 1` | 原始 P 值阈值，使用严格 `<`。 |
| `thresholds.padj_cutoff` | 数值，`0 < p <= 1` | BH 校正后 P 值阈值，使用严格 `<`。 |
| `filters.include_positive` | 布尔值 | 是否导出正相关方向筛选表。 |
| `filters.include_negative` | 布尔值 | 是否导出负相关方向筛选表。 |
| `quadrant.enabled` | 布尔值 | 是否运行独立的双分析四象限后处理；默认 `false`。 |
| `quadrant.x_target` | 字符串 | X 轴所用分析的 target；启用时必填。 |
| `quadrant.x_group` | 字符串 | X 轴所用分析的 group，例如 `All`；启用时必填。 |
| `quadrant.y_target` | 字符串 | Y 轴所用分析的 target；启用时必填。 |
| `quadrant.y_group` | 字符串 | Y 轴所用分析的 group；启用时必填。 |
| `quadrant.significance` | `p` / `padj` / `none` | 四象限 `selected` 标记采用哪种显著性定义。 |
| `quadrant.require_significant_both` | 布尔值 | 是否仅将两轴都满足所选显著性定义的 feature 标为 `selected`。 |
| `outputs.aligned_inputs` | 布尔值 | 是否导出匹配、重排后的输入及使用清单。 |
| `outputs.combined_tables` | 布尔值 | 是否导出完整结果和汇总表。 |
| `outputs.per_analysis_tables` | 布尔值 | 是否按 `group × target` 导出结果。 |
| `outputs.filtered_tables` | 布尔值 | 是否导出 P/P-adjusted 候选表及启用的方向表。 |
| `outputs.result_rds` | 布尔值 | 是否保存完整 workflow 结果对象。 |
| `intermediate.enabled` | 布尔值 | 是否保存阶段性 RDS。 |
| `intermediate.reuse_existing` | 布尔值 | 是否在签名匹配时复用已有阶段结果；仅在 `enabled: true` 时生效。 |
| `plots.overview.enabled` | 布尔值 | 是否为每个 `group × target` 生成相关总览图。 |
| `plots.overview.label_top_n_each` | 非负整数 | 每个方向标注的 top feature 数；`0` 不标注。 |
| `plots.overview.width` | 正数 | 总览图宽度，单位为英寸。 |
| `plots.overview.height` | 正数 | 总览图高度，单位为英寸。 |
| `plots.overview.dpi` | 正整数 | 位图分辨率。 |
| `plots.overview.formats` | `png` / `pdf` 列表 | 总览图格式。 |
| `plots.single_feature.enabled` | 布尔值 | 是否生成配置指定的单 feature 图。 |
| `plots.single_feature.width` | 正数 | 单 feature 图宽度，单位为英寸。 |
| `plots.single_feature.height` | 正数 | 单 feature 图高度，单位为英寸。 |
| `plots.single_feature.dpi` | 正整数 | 位图分辨率。 |
| `plots.single_feature.formats` | `png` / `pdf` 列表 | 单 feature 图格式。 |
| `plots.single_feature.items` | 对象列表 | 要绘制的 feature/target/group 组合；模块启用时至少一项。 |
| `plots.single_feature.items[].feature` | 字符串 | feature 名称。 |
| `plots.single_feature.items[].target` | 字符串 | target 名称。 |
| `plots.single_feature.items[].group` | 字符串 | `All` 或实际分组值。 |
| `plots.single_feature.items[].type` | `scatter` / `ordinal` | 连续 target 散点图或有序等级箱线图加抖动点。 |
| `plots.single_feature.items[].add_lm` | 布尔值 | 散点图是否添加线性趋势线；只影响展示，不参与 Spearman 统计。 |
| `plots.quadrant.enabled` | 布尔值 | 是否为已启用的四象限后处理生成图；四象限表与此开关相互独立。 |
| `plots.quadrant.width` | 正数 | 四象限图宽度，单位为英寸。 |
| `plots.quadrant.height` | 正数 | 四象限图高度，单位为英寸。 |
| `plots.quadrant.dpi` | 正整数 | 四象限位图分辨率。 |
| `plots.quadrant.formats` | `png` / `pdf` 列表 | 四象限图格式。 |
| `validation.manual.enabled` | 布尔值 | 是否对指定结果重新调用独立 `cor.test()` 并严格比对。 |
| `validation.manual.items` | 对象列表 | 手工复核条目；模块启用时至少一项。 |
| `validation.manual.items[].feature` | 字符串 | 要复核的 feature。 |
| `validation.manual.items[].target` | 字符串 | 要复核的 target。 |
| `validation.manual.items[].group` | 字符串 | `All` 或实际分组值。 |
| `logging.verbose` | 布尔值 | 是否同时把 INFO 日志打印到终端；警告和错误始终可见。 |

至少应启用 `analysis.run_all_samples` 或 `analysis.run_group_analysis` 之一。启用依赖条目的模块时，其 `items` 不能为空。

## 四象限后处理

四象限是新增的、独立的结果后处理，默认关闭。它从已完成的主结果中选取两个明确的 `target × group` 分析作为 X/Y 轴，按两个 `rho` 的符号分类：

- Q1：X 正、Y 正。
- Q2：X 负、Y 正。
- Q3：X 负、Y 负。
- Q4：X 正、Y 负。
- 任一轴 `rho == 0` 时标为 `Axis`。

`significance` 和 `require_significant_both` 只决定四象限结果中的 `selected` 标记。该模块不会重新计算或修改主结果、相关系数、P 值、BH 校正、方向筛选或阈值定义。`quadrant.enabled` 控制表格后处理，`plots.quadrant.enabled` 独立控制是否绘图，因此可以只导出四象限表。

## 输出

启用全部模块时，目录大致如下：

```text
output_dir/
├── aligned_inputs/
│   ├── feature_df_aligned.csv
│   ├── target_df_aligned.csv
│   ├── features_used.csv
│   ├── targets_used.csv
│   └── groups_used.csv
├── intermediate/
│   ├── 01_inputs.rds
│   ├── 02_prepared_data.rds
│   └── 03_correlation_results.rds
├── logs/
│   └── <project_id>.log
├── plots/
│   ├── <group>__<target>__correlation_overview.<format>
│   ├── <group>__<feature>__<target>__<type>.<format>
│   └── Spearman_quadrant.<format>
├── tables/
│   ├── Spearman_all_results.csv
│   ├── Spearman_result_summary.csv
│   ├── Spearman_sig_absR<r>_p<p>.csv
│   ├── Spearman_sig_absR<r>_padj<padj>.csv
│   ├── Spearman_positive_sig_absR<r>_p<p>.csv
│   ├── Spearman_positive_sig_absR<r>_padj<padj>.csv
│   ├── Spearman_negative_sig_absR<r>_p<p>.csv
│   ├── Spearman_negative_sig_absR<r>_padj<padj>.csv
│   ├── <group>__<target>__all_results.csv
│   ├── Spearman_quadrant_results.csv
│   └── Spearman_manual_validation.csv
└── Spearman_workflow_results.rds
```

原有的完整结果、摘要、阈值筛选、逐分析结果、对齐输入、总览图和 `Spearman_workflow_results.rds` 命名结构尽量保持兼容。新增的是 `logs/`、`intermediate/`、正负方向表、四象限结果/图、手工验证表和配置指定的单 feature 图；关闭对应模块时不会生成相应文件。

完整结果表保留 `feature`、`target`、`group`、`rho`、`abs_rho`、`pvalue`、`padj`、`n`、`direction`、`status`、`significant_by_p` 和 `significant_by_padj` 等字段。常见 `status` 包括 `OK`、`insufficient_n`、`constant_feature`、`constant_target` 和 `calculation_error`。

## 从中间结果继续

设置：

```yaml
intermediate:
  enabled: true
  reuse_existing: true
```

workflow 会依次保存输入、整理后数据和相关结果三个 RDS。后续绘图、导出或可选后处理失败后，修正配置或环境并使用同一配置再次运行，匹配的前序阶段会被复用。

复用不是仅按文件名判断：

- 原始输入内容、解析后的输入路径或输入格式变化，会使相关 checkpoint 失效并重新读取。
- 列映射、feature/target/group 选择、全样本/组内开关、相关方法、`min_n` 或严格样本匹配设置变化，会使整理和/或相关阶段失效。
- 阈值、方向导出、四象限、绘图、手工验证、普通输出和日志开关不改变相关计算签名，因此可以复用已有相关结果并重新后处理。
- 更换 `output_dir` 会使用新目录中的 checkpoint。
- 每个 checkpoint 还包含该阶段及其上游计算函数的源码 MD5。输入读取函数变化会使三个阶段失效，数据整理函数变化会使整理与相关阶段失效，相关核心变化只需重算相关阶段；后处理、绘图或导出代码变化不会迫使重新读取和整理原始数据。
- 损坏、结构异常或签名不匹配的 RDS 不会被盲目复用；日志会说明原因并重新计算。
- `intermediate.enabled: false` 时既不保存也不复用 checkpoint。

如果需要无条件重算，可将 `reuse_existing` 设为 `false`，或把旧 `intermediate/` 移出结果目录后重新运行。

## 常见报错

### 输入文件不存在或不可读

检查 `inputs.expression_matrix`、`inputs.score_metadata`，并记住相对路径相对于 YAML 所在目录。错误会指出具体配置字段和解析后的文件路径。

### 未知配置字段或类型错误

schema 会拒绝未知键、非法布尔值、非法格式、越界阈值和不完整模块条目。检查缩进和字段拼写。`min_n` 至少为 3，`r_cutoff` 必须在 `[0, 1]`，P 值阈值必须在 `(0, 1]`。

### 必需列不存在

检查 `columns.feature`、`columns.sample`、`columns.group` 和 `analysis.target_columns` 的大小写。错误会指出缺失列及对应输入文件/参数。

### 重复或缺失样本/feature

metadata 样本 ID、表达矩阵样本列名和 feature 名称必须唯一。空字符串、`NA` 和重复名称都会被明确报告，不能依赖 R 自动修复。

### 样本无法匹配

检查空格、大小写、分隔符和导入时列名变化。非严格模式只分析交集并记录差异；严格模式会立即报错。共同样本不足时无法进行相关分析。

### 非数值矩阵或 target

去除单位、`unknown`、`-`、逗号格式数字等字符内容，并在上游明确转换。workflow 不会为通过验证而静默改变数据含义。

### n=1、样本不足、全 NA 或零方差

- 整体或请求分组没有足够样本时，会给出包含分组和 `min_n` 的错误或警告。
- 单个 feature/target 的完整配对数不足时，结果状态为 `insufficient_n`。
- 全 `NA` 会导致完整配对数不足。
- 零方差 feature/target 分别标记为 `constant_feature` / `constant_target`，不会产生有效相关系数。
- metadata 中缺失或空分组不会被当作有效分组；日志会说明跳过或失败原因。

### 绘图失败但相关表已算完

确认 `ggplot2`/`ggrepel` 可加载、条目名称存在且绘图目录可写。保留 checkpoint 后重新运行，不需要重新读取和整理全部原始数据。

### 服务器任务显示失败

正式入口会捕获错误并返回非零退出码。配置成功读取、输出目录已确定后的运行错误会写入 `logs/<project_id>.log`；配置文件不存在、YAML 解析失败或配置校验失败发生在日志创建前，只输出到标准错误。不要只依赖交互式 R 会话状态。

## GitHub 与本地数据隔离

应提交：

- `run_spearman_analysis.R`
- `R/` 中的通用模块
- `config/config.example.yml`
- `environment.yml`
- `README.md`
- `.gitignore`
- 少量、匿名、可公开且受 Git 管理的 `examples/` 示例数据与示例脚本

不应提交：

- 真实项目配置，如 `config.local.yml`
- 真实表达矩阵、metadata 或其他 `data/`
- `results/`、`logs/`、`intermediate/`
- `*.RData`、`*.rds`、`.Rhistory`、`.Rproj.user/`
- 任何真实绝对路径、个人信息、受限数据或数据集专用编号

不要因为 `.gitignore` 新规则删除已经由 Git 正常管理的必要示例文件。真实项目配置可以位于仓库外任意位置，并通过 `--config` 传入。

## 兼容函数接口

已有 R 脚本仍可 source `R/run_spearman_workflow.R` 并调用旧函数 API。该模式不会清空工作区或自动开始分析，但只作为兼容层维护；配置验证、日志、标准退出码、外部路径解析和 checkpoint 的完整运行方式以 `run_spearman_analysis.R --config ...` 为准。

## License

This project is released under the MIT License.
