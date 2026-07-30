# Input readers and result exporters.

.spearman_read_one_input <- function(path, format, parameter_name) {
  if (!file.exists(path)) {
    stop(
      "Input file configured by '",
      parameter_name,
      "' does not exist: ",
      path,
      call. = FALSE
    )
  }
  if (file.access(path, mode = 4L) != 0L) {
    stop(
      "Input file configured by '",
      parameter_name,
      "' is not readable: ",
      path,
      call. = FALSE
    )
  }

  value <- tryCatch(
    switch(
      format,
      csv = utils::read.csv(
        path,
        check.names = FALSE,
        stringsAsFactors = FALSE
      ),
      tsv = utils::read.delim(
        path,
        check.names = FALSE,
        stringsAsFactors = FALSE
      ),
      rds = readRDS(path),
      stop(
        "Unsupported input format '",
        format,
        "' for ",
        parameter_name,
        ".",
        call. = FALSE
      )
    ),
    error = function(e) {
      stop(
        "Failed to read ",
        parameter_name,
        " from ",
        path,
        ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (!is.data.frame(value)) {
    stop(
      "Input '",
      parameter_name,
      "' must contain a data.frame; file ",
      path,
      " produced object class: ",
      paste(class(value), collapse = "/"),
      ".",
      call. = FALSE
    )
  }

  as.data.frame(
    value,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

read_spearman_inputs <- function(config, logger) {
  paths <- c(
    expression_matrix = config$inputs$expression_matrix,
    score_metadata = config$inputs$score_metadata
  )

  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    detail <- paste(
      paste0(names(missing), "=", unname(missing)),
      collapse = "; "
    )
    stop("Configured input file(s) do not exist: ", detail, call. = FALSE)
  }

  input_md5 <- as.character(tools::md5sum(paths))
  names(input_md5) <- names(paths)

  logger(
    "INFO",
    "Reading expression matrix: ",
    config$inputs$expression_matrix
  )
  feature_df <- .spearman_read_one_input(
    config$inputs$expression_matrix,
    config$inputs$expression_format,
    "inputs.expression_matrix"
  )

  logger(
    "INFO",
    "Reading score/metadata table: ",
    config$inputs$score_metadata
  )
  target_df <- .spearman_read_one_input(
    config$inputs$score_metadata,
    config$inputs$metadata_format,
    "inputs.score_metadata"
  )

  logger(
    "INFO",
    "Loaded expression matrix ",
    nrow(feature_df),
    " x ",
    ncol(feature_df),
    " and metadata ",
    nrow(target_df),
    " x ",
    ncol(target_df),
    "."
  )

  list(
    feature_df = feature_df,
    target_df = target_df,
    input_md5 = input_md5
  )
}

.spearman_atomic_save_rds <- function(value, path) {
  parent <- dirname(path)
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create output directory: ", parent, call. = FALSE)
  }

  tmp <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = parent
  )
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(value, tmp, version = 3L)

  if (file.exists(path) && !file.remove(path)) {
    stop("Unable to replace output file: ", path, call. = FALSE)
  }
  if (!file.rename(tmp, path)) {
    stop("Unable to finalize output file: ", path, call. = FALSE)
  }
  invisible(path)
}

.spearman_export_aligned_inputs <- function(result, directory) {
  .spearman_write_csv(
    result$feature_df_aligned,
    file.path(directory, "feature_df_aligned.csv")
  )
  .spearman_write_csv(
    result$target_df_aligned,
    file.path(directory, "target_df_aligned.csv")
  )
  .spearman_write_csv(
    data.frame(
      feature = result$features_used,
      stringsAsFactors = FALSE
    ),
    file.path(directory, "features_used.csv")
  )
  .spearman_write_csv(
    data.frame(
      target = result$targets_used,
      stringsAsFactors = FALSE
    ),
    file.path(directory, "targets_used.csv")
  )

  if (!is.null(result$parameters$group_col) &&
      length(result$groups_used) > 0L) {
    .spearman_write_csv(
      data.frame(
        group = result$groups_used,
        stringsAsFactors = FALSE
      ),
      file.path(directory, "groups_used.csv")
    )
  }
}

.spearman_export_combined_tables <- function(result, directory) {
  .spearman_write_csv(
    result$results,
    file.path(directory, "Spearman_all_results.csv")
  )
  .spearman_write_csv(
    result$summary,
    file.path(directory, "Spearman_result_summary.csv")
  )
}

.spearman_export_filtered_tables <- function(result, config, directory) {
  r_text <- .spearman_format_threshold(
    config$thresholds$r_cutoff
  )
  p_text <- .spearman_format_threshold(
    config$thresholds$p_cutoff
  )
  padj_text <- .spearman_format_threshold(
    config$thresholds$padj_cutoff
  )

  .spearman_write_csv(
    result$significant_p,
    file.path(
      directory,
      paste0(
        "Spearman_sig_absR",
        r_text,
        "_p",
        p_text,
        ".csv"
      )
    )
  )
  .spearman_write_csv(
    result$significant_padj,
    file.path(
      directory,
      paste0(
        "Spearman_sig_absR",
        r_text,
        "_padj",
        padj_text,
        ".csv"
      )
    )
  )

  directional_names <- names(result$directional_results)
  for (name in directional_names) {
    suffix <- if (grepl("_padj$", name)) {
      paste0("_padj", padj_text)
    } else {
      paste0("_p", p_text)
    }
    direction <- sub("_(p|padj)$", "", name)

    .spearman_write_csv(
      result$directional_results[[name]],
      file.path(
        directory,
        paste0(
          "Spearman_",
          direction,
          "_sig_absR",
          r_text,
          suffix,
          ".csv"
        )
      )
    )
  }
}

.spearman_export_per_analysis <- function(result, directory) {
  combinations <- unique(
    result$results[, c("group", "target"), drop = FALSE]
  )
  labels <- paste(
    combinations$group,
    combinations$target,
    sep = "__"
  )
  .spearman_assert_unique_filenames(
    labels,
    "per-analysis result tables"
  )

  for (i in seq_len(nrow(combinations))) {
    group_name <- combinations$group[i]
    target_name <- combinations$target[i]
    one <- result$results[
      result$results$group == group_name &
        result$results$target == target_name,
      ,
      drop = FALSE
    ]
    filename <- paste0(
      .spearman_safe_filename(group_name),
      "__",
      .spearman_safe_filename(target_name),
      "__all_results.csv"
    )
    .spearman_write_csv(one, file.path(directory, filename))
  }
}

export_spearman_results <- function(result, config, logger) {
  outdir <- config$output_dir
  table_dir <- file.path(outdir, "tables")
  input_dir <- file.path(outdir, "aligned_inputs")

  if (!dir.exists(outdir) &&
      !dir.create(outdir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create output directory: ", outdir, call. = FALSE)
  }

  if (isTRUE(config$outputs$aligned_inputs)) {
    .spearman_export_aligned_inputs(result, input_dir)
    logger("INFO", "Exported aligned inputs to ", input_dir)
  }

  if (isTRUE(config$outputs$combined_tables)) {
    .spearman_export_combined_tables(result, table_dir)
  }

  if (isTRUE(config$outputs$filtered_tables)) {
    .spearman_export_filtered_tables(result, config, table_dir)
  }

  if (isTRUE(config$outputs$per_analysis_tables)) {
    .spearman_export_per_analysis(result, table_dir)
  }

  if (!is.null(result$quadrant)) {
    .spearman_write_csv(
      result$quadrant,
      file.path(table_dir, "Spearman_quadrant_results.csv")
    )
  }

  if (!is.null(result$manual_validation)) {
    .spearman_write_csv(
      result$manual_validation,
      file.path(table_dir, "Spearman_manual_validation.csv")
    )
  }

  if (isTRUE(config$outputs$result_rds)) {
    path <- file.path(outdir, "Spearman_workflow_results.rds")
    .spearman_atomic_save_rds(result, path)
    logger("INFO", "Saved workflow result object: ", path)
  }

  logger("INFO", "Result export completed: ", outdir)
  invisible(outdir)
}
