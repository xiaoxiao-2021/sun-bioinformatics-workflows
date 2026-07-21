# Spearman Correlation Workflow

A reusable R workflow for batch Spearman rank correlation analysis between multiple features and one or more target variables.

该工作流用于批量分析多个 `feature` 与一个或多个 `target` 之间的 Spearman 秩相关关系，支持全样本分析、组内分析、样本匹配、多重检验校正、结果筛选和基础可视化。

---

## Features

- Standardized `data.frame` input
- Automatic sample matching and order alignment
- Whole-cohort correlation analysis
- Within-group correlation analysis
- Multiple target variables
- Missing-value handling
- Constant-variable detection
- Spearman `rho` and p-value calculation
- Benjamini–Hochberg multiple-testing correction
- Significant feature screening
- Correlation overview plots
- Single-feature scatter plots
- CSV, PNG, PDF and RData output

---

## Repository structure

```text
spearman-correlation/
├── README.md
├── R/
│   └── spearman_correlation_workflow.R
├── examples/
└── figures/
```

The core workflow is stored in:

```text
R/spearman_correlation_workflow.R
```

---

## Input data

The workflow uses two standardized `data.frame` objects:

```r
feature_df
target_df
```

Although the external input is standardized as `data.frame`, the workflow converts the feature values to a numeric matrix internally. The actual statistical calculation performed by `cor.test()` uses two matched numeric vectors.

### `feature_df`

Required structure:

- Each row represents one feature
- One column stores feature names
- All other columns represent samples
- Sample columns must be numeric
- Feature names must not be duplicated
- Sample column names must not be duplicated

Example:

| feature | S01 | S02 | S03 | S04 |
|---|---:|---:|---:|---:|
| Gene_A | 4.2 | 5.1 | 3.8 | 4.6 |
| Gene_B | 2.3 | 2.8 | 4.1 | 3.7 |
| Gene_C | 7.4 | 6.9 | 8.0 | 7.2 |

If the original data are stored as a numeric matrix and feature names are stored in row names:

```r
feature_df <- tibble::rownames_to_column(
  as.data.frame(feature_mat, check.names = FALSE),
  var = "feature"
)
```

### `target_df`

Required structure:

- Each row represents one sample
- One column stores sample names
- One or more numeric columns store target variables
- An optional group column can be used for within-group analysis
- Each sample must appear only once

Example:

| sample | Score_A | Score_B | Group |
|---|---:|---:|---|
| S01 | 0.42 | 1.15 | IA |
| S02 | 0.81 | 0.94 | IA |
| S03 | 1.10 | 0.75 | ENEG |
| S04 | 0.35 | 1.34 | ENEG |

Possible target variables include:

- ssGSEA scores
- pathway scores
- immune infiltration scores
- clinical indices
- laboratory measurements
- phenotype scores
- other continuous or ordinal variables

---

## Required R packages

```r
install.packages(
  c(
    "dplyr",
    "tibble",
    "purrr",
    "readr",
    "ggplot2",
    "ggrepel"
  )
)
```

---

## Quick start

### 1. Prepare the input objects

Before running the workflow, create or import:

```r
feature_df
target_df
```

For example:

```r
feature_df <- read.csv(
  "feature_df.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

target_df <- read.csv(
  "target_df.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
```

Using `check.names = FALSE` helps prevent R from automatically modifying sample names.

### 2. Modify the parameter section

Open:

```text
R/spearman_correlation_workflow.R
```

Check and modify the following parameters.

#### Parameters that must be checked for each dataset

```r
feature_col <- "feature"
sample_col <- "sample"

target_cols <- c(
  "T_Cell_Exhaustion_Progenitor_like"
)

group_col <- "Group"

selected_groups <- c(
  "IA",
  "ENEG"
)
```

Column names are case-sensitive.

For example:

```text
sample
Sample
```

are treated as different column names.

#### Analysis switches

```r
run_all_samples <- TRUE
run_group_analysis <- TRUE
```

- `run_all_samples = TRUE`: run correlation analysis using all matched samples
- `run_group_analysis = TRUE`: run correlation analysis within the groups listed in `selected_groups`

#### Statistical parameters

```r
min_n <- 5
r_cutoff <- 0.5
p_cutoff <- 0.05
padj_cutoff <- 0.05
```

#### Sample-matching behavior

```r
strict_sample_match <- FALSE
```

- `FALSE`: report unmatched samples and continue using common samples
- `TRUE`: stop the workflow when any unmatched sample is detected

#### Output directory

```r
outdir <- "results/spearman_correlation"
```

### 3. Run the workflow

From the repository root directory:

```r
source(
  "spearman-correlation/R/spearman_correlation_workflow.R"
)
```

---

## Sample matching

The workflow first extracts sample names from both input objects.

```r
feature_samples <- setdiff(
  colnames(feature_df),
  feature_col
)

target_samples <- target_df[[sample_col]]
```

It then uses `setdiff()` to identify unmatched samples:

```r
only_in_target <- setdiff(
  target_samples,
  feature_samples
)

only_in_feature <- setdiff(
  feature_samples,
  target_samples
)
```

It uses `intersect()` to obtain common samples:

```r
common_samples <- intersect(
  feature_samples,
  target_samples
)
```

Both input objects are then reordered according to the same sample vector.

Before correlation analysis, the workflow confirms:

```r
identical(
  colnames(feature_df_aligned)[-1],
  target_df_aligned[[sample_col]]
)
```

The result must be:

```r
TRUE
```

This step is essential because `cor.test()` compares vector positions and does not automatically match samples by name.

---

## Statistical analysis

For each feature-target pair, the workflow performs:

```r
cor.test(
  x = x_complete,
  y = y_complete,
  method = "spearman",
  exact = FALSE
)
```

Where:

- `x_complete` is the numeric vector of one feature
- `y_complete` is the matched numeric target vector
- `exact = FALSE` allows approximate p-value calculation when tied ranks are present

The workflow reports:

| Field | Description |
|---|---|
| `feature` | Feature name |
| `target` | Target variable |
| `group` | Analysis group |
| `rho` | Spearman correlation coefficient |
| `abs_rho` | Absolute correlation coefficient |
| `pvalue` | Raw p-value |
| `padj` | BH-adjusted p-value |
| `n` | Number of complete samples |
| `direction` | Positive or negative correlation |
| `status` | Calculation status |

Multiple-testing correction is performed within each:

```text
target × group
```

using:

```r
p.adjust(
  pvalue,
  method = "BH"
)
```

---

## Output objects

The main R objects generated by the workflow are:

```text
feature_df_aligned
target_df_aligned
feature_mat
cor_results
cor_sig_p
cor_sig_padj
result_summary
```

### `feature_df_aligned`

Feature data after sample matching and order alignment.

### `target_df_aligned`

Target data after sample matching and order alignment.

### `feature_mat`

Numeric matrix used internally for batch calculation.

### `cor_results`

Complete correlation results.

### `cor_sig_p`

Results satisfying:

```text
|rho| >= r_cutoff
pvalue < p_cutoff
```

### `cor_sig_padj`

Results satisfying:

```text
|rho| >= r_cutoff
padj < padj_cutoff
```

### `result_summary`

Summary of valid and significant results for each `target × group`.

---

## Output files

Results are saved by default to:

```text
results/spearman_correlation/
```

Output structure:

```text
results/spearman_correlation/
├── aligned_inputs/
│   ├── feature_df_aligned.csv
│   └── target_df_aligned.csv
├── plots/
│   ├── *_correlation_overview.png
│   └── *_correlation_overview.pdf
├── tables/
│   ├── Spearman_all_results.csv
│   ├── Spearman_sig_*.csv
│   ├── Spearman_result_summary.csv
│   └── group__target__all_results.csv
└── Spearman_workflow_results.RData
```

---

## Correlation overview plot

The workflow can generate one overview plot for each `target × group`.

In this plot:

- Each point represents one feature
- The x-axis represents Spearman `rho`
- The y-axis represents `-log10(p value)`
- Positive, negative and non-significant results are displayed separately
- Dashed lines represent the selected `rho` and p-value thresholds

---

## Single-feature scatter plot

The workflow also provides:

```r
plot_feature_scatter()
```

Example:

```r
p_gene <- plot_feature_scatter(
  feature_name = "RIPOR2",
  target_name = "T_Cell_Exhaustion_Progenitor_like",
  group_name = "IA",
  feature_mat = feature_mat,
  target_df_aligned = target_df_aligned,
  sample_col = sample_col,
  group_col = group_col
)

print(p_gene)
```

Each point represents one sample.

The figure reports:

- Spearman `rho`
- p-value
- effective sample size

---

## Recommended validation

Before interpreting the results, check:

```r
str(feature_df)
str(target_df)
```

Confirm sample alignment:

```r
identical(
  colnames(feature_df_aligned)[-1],
  target_df_aligned[[sample_col]]
)
```

Inspect unmatched samples:

```r
only_in_target
only_in_feature
```

It is also recommended to manually verify one feature:

```r
x <- as.numeric(
  feature_mat[
    "RIPOR2",
    target_df_aligned[[sample_col]]
  ]
)

y <- target_df_aligned[
  ["T_Cell_Exhaustion_Progenitor_like"]
]

cor.test(
  x,
  y,
  method = "spearman",
  exact = FALSE
)
```

The manually calculated result should agree with the corresponding row in `cor_results`.

---

## Interpretation notes

A positive Spearman coefficient means that the feature and target tend to increase together.

A negative coefficient means that one variable tends to decrease as the other increases.

The interpretation should consider:

- Correlation direction
- Correlation strength
- Raw and adjusted p-values
- Effective sample size
- Scatterplot distribution
- Biological plausibility
- Independent validation

Correlation does not establish causation.

---

## Common issues

### Unmatched sample names

Possible causes:

- Leading or trailing spaces
- Different letter case
- Different separators
- R automatically modified sample names
- Missing samples in one input table

### Non-numeric columns

Possible causes:

- `"-"`
- `"unknown"`
- `"not detected"`
- Empty strings
- Unit characters

### Duplicate samples or features

Each sample in `target_df` must appear only once.

Duplicated features should be merged or otherwise resolved before running the workflow.

### Insufficient sample size

Within-group analysis reduces the number of available samples.

Results based on small groups should be interpreted cautiously.

### Constant variables

A feature or target with no variation cannot be used to calculate correlation.

### Tied ranks

Repeated values may produce tied ranks. The workflow uses:

```r
exact = FALSE
```

to calculate an approximate p-value.

---

## Detailed documentation

A detailed Chinese explanation of the workflow is available here:

[Spearman correlation workflow](https://xiaoxiao-2021.github.io/sun-lab-life-astro-v2/bioinformatics/spearman-correlation-workflow/)

---

## License

This project is released under the MIT License.


## Development note

This workflow was developed with AI-assisted code drafting and was reviewed, tested, and maintained by the repository author.