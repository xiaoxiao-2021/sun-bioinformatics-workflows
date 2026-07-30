# Spearman workflow orchestration and backward-compatible function API.
#
# Sourcing this file loads all workflow modules relative to this file. It does
# not read data, change the working directory, or start an analysis.

.spearman_locate_source_file <- function(expected_name) {
  frames <- sys.frames()
  if (length(frames) == 0L) {
    return(NULL)
  }

  for (i in rev(seq_along(frames))) {
    candidate <- frames[[i]]$ofile
    if (is.null(candidate) ||
        !is.character(candidate) ||
        length(candidate) != 1L ||
        is.na(candidate) ||
        !nzchar(candidate)) {
      next
    }
    candidate <- tryCatch(
      normalizePath(
        candidate,
        winslash = "/",
        mustWork = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(candidate) &&
        identical(basename(candidate), expected_name)) {
      return(candidate)
    }
  }

  NULL
}

.spearman_this_file <- .spearman_locate_source_file(
  "run_spearman_workflow.R"
)

if (is.null(.spearman_this_file)) {
  stop(
    "Unable to locate R/run_spearman_workflow.R while sourcing it. ",
    "Source the file by path rather than pasting it into a session.",
    call. = FALSE
  )
}

.spearman_r_dir <- dirname(.spearman_this_file)
.spearman_modules <- c(
  "runtime_utils.R",
  "config.R",
  "input_output.R",
  "analysis_core.R",
  "postprocessing.R",
  "visualization.R"
)

for (.spearman_module in .spearman_modules) {
  .spearman_module_path <- file.path(
    .spearman_r_dir,
    .spearman_module
  )
  if (!file.exists(.spearman_module_path)) {
    stop(
      "Workflow module is missing: ",
      .spearman_module_path,
      call. = FALSE
    )
  }
  source(
    .spearman_module_path,
    local = environment(),
    encoding = "UTF-8"
  )
}

.spearman_function_code_md5 <- function(
    function_names,
    source_environment
) {
  missing_functions <- function_names[
    !vapply(
      function_names,
      exists,
      logical(1),
      envir = source_environment,
      mode = "function",
      inherits = FALSE
    )
  ]
  if (length(missing_functions) > 0L) {
    stop(
      "Cannot calculate checkpoint code signature; missing function(s): ",
      paste(missing_functions, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  definitions <- lapply(
    function_names,
    function(function_name) {
      function_object <- get(
        function_name,
        envir = source_environment,
        mode = "function",
        inherits = FALSE
      )
      list(
        name = function_name,
        formals = formals(function_object),
        body = body(function_object)
      )
    }
  )

  signature_file <- tempfile("spearman-code-signature-")
  on.exit(unlink(signature_file), add = TRUE)
  writeBin(
    serialize(definitions, NULL, version = 2L),
    signature_file
  )
  unname(as.character(tools::md5sum(signature_file)))
}

.spearman_compute_code_signatures <- function(
    source_environment
) {
  input_functions <- c(
    ".spearman_read_one_input",
    "read_spearman_inputs"
  )
  preparation_functions <- c(
    ".spearman_core_label",
    ".spearman_core_preview",
    ".spearman_core_validate_logger",
    ".spearman_core_log",
    ".spearman_core_warn",
    ".spearman_core_validate_flag",
    ".spearman_core_validate_min_n",
    ".spearman_core_scalar_name",
    ".spearman_core_trim_columns",
    ".spearman_core_character_selection",
    ".spearman_core_warn_degenerate_inputs",
    "prepare_spearman_data"
  )
  correlation_functions <- c(
    ".spearman_core_one_target_group",
    ".spearman_core_validate_prepared",
    "run_spearman_correlations"
  )

  stage_functions <- list(
    inputs = input_functions,
    prepared = c(
      input_functions,
      preparation_functions
    ),
    correlations = c(
      input_functions,
      preparation_functions,
      correlation_functions
    )
  )

  vapply(
    stage_functions,
    .spearman_function_code_md5,
    character(1),
    source_environment = source_environment,
    USE.NAMES = TRUE
  )
}

.spearman_workflow_code_signatures <-
  .spearman_compute_code_signatures(environment())

rm(
  .spearman_locate_source_file,
  .spearman_this_file,
  .spearman_r_dir,
  .spearman_modules,
  .spearman_module,
  .spearman_module_path
)

.spearman_console_logger <- function(verbose) {
  force(verbose)
  function(level = "INFO", ...) {
    level <- toupper(as.character(level)[1L])
    if (isTRUE(verbose) || level %in% c("WARN", "ERROR")) {
      message(
        "[",
        level,
        "] ",
        paste0(..., collapse = "")
      )
    }
    invisible(NULL)
  }
}

.spearman_assemble_result <- function(
    prepared,
    screened,
    config,
    call,
    reused_stages = character(0)
) {
  result <- list(
    call = call,
    project_id = config$project_id,
    config = config,
    parameters = list(
      feature_col = config$columns$feature,
      sample_col = config$columns$sample,
      selected_features = config$analysis$selected_features,
      target_cols = prepared$targets_used,
      run_all_samples = config$analysis$run_all_samples,
      run_group_analysis = config$analysis$run_group_analysis,
      group_col = config$columns$group,
      selected_groups = prepared$groups_used,
      min_n = config$analysis$min_n,
      correlation_method = config$analysis$correlation_method,
      r_cutoff = config$thresholds$r_cutoff,
      p_cutoff = config$thresholds$p_cutoff,
      padj_cutoff = config$thresholds$padj_cutoff,
      strict_sample_match = config$analysis$strict_sample_match
    ),
    sample_report = prepared$sample_report,
    features_used = prepared$features_used,
    targets_used = prepared$targets_used,
    groups_used = prepared$groups_used,
    feature_df_aligned = prepared$feature_df_aligned,
    target_df_aligned = prepared$target_df_aligned,
    feature_mat = prepared$feature_mat,
    results = screened$results,
    significant_p = screened$significant_p,
    significant_padj = screened$significant_padj,
    directional_results = screened$directional,
    summary = build_spearman_summary(screened$results),
    quadrant = NULL,
    manual_validation = NULL,
    overview_plots = list(),
    single_feature_plots = list(),
    quadrant_plot = NULL,
    run_info = list(
      completed_at = format(
        Sys.time(),
        "%Y-%m-%dT%H:%M:%S%z"
      ),
      reused_stages = reused_stages
    )
  )

  class(result) <- c("spearman_workflow_result", "list")
  result
}

#' Run a configuration-driven Spearman analysis project.
#'
#' @param config Normalized configuration returned by
#'   `read_spearman_config()`.
#' @param logger Logger function with signature `logger(level, ...)`.
#'
#' @return A `spearman_workflow_result`.
run_spearman_project <- function(config, logger) {
  if (!is.list(config) || is.null(config$output_dir)) {
    stop(
      "config must be returned by read_spearman_config().",
      call. = FALSE
    )
  }
  if (!is.function(logger)) {
    stop("logger must be a function.", call. = FALSE)
  }

  if (!dir.exists(config$output_dir) &&
      !dir.create(
        config$output_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )) {
    stop(
      "Unable to create output directory: ",
      config$output_dir,
      call. = FALSE
    )
  }

  logger(
    "INFO",
    "Starting project '",
    config$project_id,
    "'."
  )

  input_paths <- c(
    expression_matrix = config$inputs$expression_matrix,
    score_metadata = config$inputs$score_metadata
  )
  missing_inputs <- input_paths[!file.exists(input_paths)]
  if (length(missing_inputs) > 0L) {
    stop(
      "Configured input file(s) do not exist: ",
      paste(
        paste0(
          names(missing_inputs),
          "=",
          unname(missing_inputs)
        ),
        collapse = "; "
      ),
      call. = FALSE
    )
  }

  input_md5 <- as.character(tools::md5sum(input_paths))
  names(input_md5) <- names(input_paths)
  signature <- .spearman_config_signature(config, input_md5)

  intermediate_dir <- file.path(
    config$output_dir,
    "intermediate"
  )
  checkpoint_paths <- c(
    inputs = file.path(intermediate_dir, "01_inputs.rds"),
    prepared = file.path(intermediate_dir, "02_prepared_data.rds"),
    correlations = file.path(
      intermediate_dir,
      "03_correlation_results.rds"
    )
  )
  reuse <- (
    isTRUE(config$intermediate$enabled) &&
      isTRUE(config$intermediate$reuse_existing)
  )
  reused_stages <- character(0)

  inputs <- .spearman_load_checkpoint(
    checkpoint_paths[["inputs"]],
    .spearman_stage_signature(signature, "inputs"),
    "inputs",
    reuse,
    logger
  )
  if (is.null(inputs)) {
    inputs <- read_spearman_inputs(config, logger)
    .spearman_save_checkpoint(
      inputs,
      .spearman_stage_signature(signature, "inputs"),
      checkpoint_paths[["inputs"]],
      "inputs",
      config$intermediate$enabled,
      logger
    )
  } else {
    reused_stages <- c(reused_stages, "inputs")
  }

  prepared <- .spearman_load_checkpoint(
    checkpoint_paths[["prepared"]],
    .spearman_stage_signature(signature, "prepared"),
    "prepared",
    reuse,
    logger
  )
  if (is.null(prepared)) {
    prepared <- prepare_spearman_data(
      feature_df = inputs$feature_df,
      target_df = inputs$target_df,
      feature_path = config$inputs$expression_matrix,
      metadata_path = config$inputs$score_metadata,
      feature_col = config$columns$feature,
      sample_col = config$columns$sample,
      group_col = config$columns$group,
      selected_features = config$analysis$selected_features,
      target_cols = config$analysis$target_columns,
      run_group_analysis = config$analysis$run_group_analysis,
      selected_groups = config$analysis$selected_groups,
      min_n = config$analysis$min_n,
      strict_sample_match = config$analysis$strict_sample_match,
      logger = logger
    )
    .spearman_save_checkpoint(
      prepared,
      .spearman_stage_signature(signature, "prepared"),
      checkpoint_paths[["prepared"]],
      "prepared",
      config$intermediate$enabled,
      logger
    )
  } else {
    reused_stages <- c(reused_stages, "prepared")
  }

  correlations <- .spearman_load_checkpoint(
    checkpoint_paths[["correlations"]],
    .spearman_stage_signature(signature, "correlations"),
    "correlations",
    reuse,
    logger
  )
  if (is.null(correlations)) {
    correlations <- run_spearman_correlations(
      prepared = prepared,
      run_all_samples = config$analysis$run_all_samples,
      run_group_analysis = config$analysis$run_group_analysis,
      min_n = config$analysis$min_n,
      correlation_method = config$analysis$correlation_method,
      logger = logger
    )
    .spearman_save_checkpoint(
      correlations,
      .spearman_stage_signature(signature, "correlations"),
      checkpoint_paths[["correlations"]],
      "correlations",
      config$intermediate$enabled,
      logger
    )
  } else {
    reused_stages <- c(reused_stages, "correlations")
  }

  screened <- apply_spearman_thresholds(
    correlations,
    config$thresholds$r_cutoff,
    config$thresholds$p_cutoff,
    config$thresholds$padj_cutoff,
    config$filters$include_positive,
    config$filters$include_negative
  )

  result <- .spearman_assemble_result(
    prepared,
    screened,
    config,
    match.call(),
    reused_stages
  )

  result$quadrant <- run_spearman_quadrant(
    result$results,
    config$quadrant,
    logger
  )
  result$manual_validation <- run_spearman_manual_validation(
    result,
    config$validation$manual,
    logger
  )

  plot_bundle <- create_spearman_plots(
    result,
    config,
    logger
  )
  result$overview_plots <- plot_bundle$overview
  result$single_feature_plots <- plot_bundle$single_feature
  result$quadrant_plot <- plot_bundle$quadrant

  save_spearman_plots(plot_bundle, config, logger)
  export_spearman_results(result, config, logger)

  logger(
    "INFO",
    "Spearman workflow completed for project '",
    config$project_id,
    "'."
  )
  result
}

#' Run the function API retained for existing scripts.
#'
#' This preserves the previous public arguments and statistical definitions.
#' New command-line analyses should use `run_spearman_analysis.R --config ...`.
run_spearman_workflow <- function(
    feature_df,
    target_df,
    feature_col = "feature",
    sample_col = "sample",
    selected_features = NULL,
    target_cols = NULL,
    run_all_samples = TRUE,
    run_group_analysis = FALSE,
    group_col = NULL,
    selected_groups = NULL,
    min_n = 5L,
    r_cutoff = 0.5,
    p_cutoff = 0.05,
    padj_cutoff = 0.05,
    strict_sample_match = FALSE,
    outdir = "results/spearman_correlation",
    save_results = FALSE,
    make_overview_plot = TRUE,
    label_top_n_each = 0L,
    verbose = TRUE
) {
  logger <- .spearman_console_logger(verbose)
  prepared <- prepare_spearman_data(
    feature_df = feature_df,
    target_df = target_df,
    feature_path = "<feature_df argument>",
    metadata_path = "<target_df argument>",
    feature_col = feature_col,
    sample_col = sample_col,
    group_col = group_col,
    selected_features = selected_features,
    target_cols = target_cols,
    run_group_analysis = run_group_analysis,
    selected_groups = selected_groups,
    min_n = min_n,
    strict_sample_match = strict_sample_match,
    logger = logger
  )
  correlations <- run_spearman_correlations(
    prepared,
    run_all_samples,
    run_group_analysis,
    min_n,
    "spearman",
    logger
  )
  screened <- apply_spearman_thresholds(
    correlations,
    r_cutoff,
    p_cutoff,
    padj_cutoff,
    TRUE,
    TRUE
  )

  config <- list(
    project_id = "function_api",
    output_dir = outdir,
    inputs = list(
      expression_matrix = "<feature_df argument>",
      score_metadata = "<target_df argument>",
      expression_format = "data.frame",
      metadata_format = "data.frame"
    ),
    columns = list(
      feature = feature_col,
      sample = sample_col,
      group = group_col
    ),
    analysis = list(
      selected_features = selected_features,
      target_columns = target_cols,
      run_all_samples = run_all_samples,
      run_group_analysis = run_group_analysis,
      selected_groups = selected_groups,
      correlation_method = "spearman",
      min_n = as.integer(min_n),
      strict_sample_match = strict_sample_match
    ),
    thresholds = list(
      r_cutoff = r_cutoff,
      p_cutoff = p_cutoff,
      padj_cutoff = padj_cutoff
    ),
    filters = list(
      include_positive = TRUE,
      include_negative = TRUE
    ),
    quadrant = list(enabled = FALSE),
    outputs = list(
      aligned_inputs = TRUE,
      combined_tables = TRUE,
      per_analysis_tables = TRUE,
      filtered_tables = TRUE,
      result_rds = TRUE
    ),
    intermediate = list(
      enabled = FALSE,
      reuse_existing = FALSE
    ),
    plots = list(
      overview = list(
        enabled = make_overview_plot,
        label_top_n_each = as.integer(label_top_n_each),
        width = 7,
        height = 5.5,
        dpi = 300L,
        formats = c("png", "pdf")
      ),
      single_feature = list(
        enabled = FALSE,
        width = 5,
        height = 4,
        dpi = 300L,
        formats = c("png", "pdf"),
        items = list()
      ),
      quadrant = list(
        enabled = FALSE,
        width = 7,
        height = 5.5,
        dpi = 300L,
        formats = c("png", "pdf")
      )
    ),
    validation = list(
      manual = list(enabled = FALSE, items = list())
    ),
    logging = list(verbose = verbose)
  )

  result <- .spearman_assemble_result(
    prepared,
    screened,
    config,
    match.call()
  )
  plot_bundle <- create_spearman_plots(result, config, logger)
  result$overview_plots <- plot_bundle$overview

  if (isTRUE(save_results)) {
    save_spearman_plots(plot_bundle, config, logger)
    export_spearman_results(result, config, logger)
  }

  if (isTRUE(verbose)) {
    print(result$summary)
    logger("INFO", "Spearman workflow completed.")
  }
  result
}

print.spearman_workflow_result <- function(x, ...) {
  cat("<spearman_workflow_result>\n", sep = "")
  cat("Project: ", x$project_id, "\n", sep = "")
  cat("Features: ", length(x$features_used), "\n", sep = "")
  cat(
    "Targets: ",
    paste(x$targets_used, collapse = ", "),
    "\n",
    sep = ""
  )
  cat(
    "Analyses: ",
    paste(unique(x$results$group), collapse = ", "),
    "\n",
    sep = ""
  )
  cat(
    "Matched samples: ",
    length(x$sample_report$common_samples),
    "\n\n",
    sep = ""
  )
  print(x$summary)
  invisible(x)
}
