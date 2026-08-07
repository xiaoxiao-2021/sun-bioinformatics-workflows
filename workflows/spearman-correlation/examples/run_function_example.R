# Example of the retained function API.
#
# Formal command-line analyses should use:
# Rscript run_spearman_analysis.R --config config/config.example.yml

.example_script_path <- function() {
  file_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  candidate <- if (length(file_arg) == 1L) {
    sub("^--file=", "", file_arg)
  } else {
    tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  }

  if (is.null(candidate)) {
    stop(
      "Unable to locate run_function_example.R.",
      call. = FALSE
    )
  }

  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}

example_file <- .example_script_path()
workflow_dir <- dirname(dirname(example_file))

feature_df <- read.csv(
  file.path(
    workflow_dir,
    "examples",
    "feature_df_example.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
target_df <- read.csv(
  file.path(
    workflow_dir,
    "examples",
    "target_df_example.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

source(
  file.path(
    workflow_dir,
    "R",
    "run_spearman_workflow.R"
  ),
  encoding = "UTF-8"
)

example_result <- run_spearman_workflow(
  feature_df = feature_df,
  target_df = target_df,
  feature_col = "feature",
  sample_col = "Sample",
  selected_features = NULL,
  target_cols = "Score_A",
  run_all_samples = TRUE,
  run_group_analysis = TRUE,
  group_col = "Group",
  selected_groups = NULL,
  min_n = 5,
  r_cutoff = 0.5,
  p_cutoff = 0.05,
  padj_cutoff = 0.05,
  strict_sample_match = FALSE,
  save_results = FALSE,
  make_overview_plot = TRUE,
  label_top_n_each = 0,
  verbose = TRUE
)

print(example_result)

p_example <- plot_spearman_feature(
  workflow_result = example_result,
  feature_name = "Gene_Pos",
  target_name = "Score_A",
  group_name = "All",
  plot_type = "scatter",
  add_lm = TRUE
)
print(p_example)
