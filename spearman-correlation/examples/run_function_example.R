# ============================================================
# Example: run the function-based Spearman workflow
# ============================================================

# Run this file from the repository root.

# 1. Read example input data
feature_df <- read.csv(
  "spearman-correlation/examples/feature_df_example.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

target_df <- read.csv(
  "spearman-correlation/examples/target_df_example.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 2. Load function definitions
source(
  "spearman-correlation/R/run_spearman_workflow.R"
)

# 3. Run the workflow
example_result <- run_spearman_workflow(
  feature_df = feature_df,
  target_df = target_df,

  feature_col = "feature",
  sample_col = "Sample",

  # NULL = use all features
  selected_features = NULL,

  # Replace this with the target column in your example data
  target_cols = c("Score_A"),

  run_all_samples = TRUE,
  run_group_analysis = TRUE,

  group_col = "Group",

  # NULL = automatically use all non-missing groups
  selected_groups = NULL,

  min_n = 5,
  r_cutoff = 0.5,
  p_cutoff = 0.05,
  padj_cutoff = 0.05,

  strict_sample_match = FALSE,

  outdir = "spearman-correlation/examples/example_results",
  save_results = FALSE,
  make_overview_plot = TRUE,
  label_top_n_each = 0,
  verbose = TRUEFALSE
)

# 4. Inspect returned objects
example_result
head(example_result$results)
example_result$summary

# 5. Optional: draw one feature-target scatter plot
#
# p_example <- plot_spearman_feature(
#   workflow_result = example_result,
#   feature_name = "Gene_A",
#   target_name = "Score_A",
#   group_name = "All"
# )
#
# print(p_example)
