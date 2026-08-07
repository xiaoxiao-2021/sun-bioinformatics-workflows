# Core data preparation and correlation functions.
#
# This file intentionally contains function definitions only. It can be sourced
# safely from an interactive session, Rscript, or a server-side runner without
# changing the working directory or attaching packages.

.spearman_core_label <- function(path, fallback) {
  if (is.character(path) &&
      length(path) == 1L &&
      !is.na(path) &&
      nzchar(trimws(path))) {
    paste0(fallback, " [", path, "]")
  } else {
    fallback
  }
}

.spearman_core_preview <- function(x, n = 20L) {
  x <- unique(as.character(x))
  shown <- utils::head(x, n)
  suffix <- if (length(x) > length(shown)) {
    paste0(" ... (+", length(x) - length(shown), " more)")
  } else {
    ""
  }
  paste0(paste(shown, collapse = ", "), suffix)
}

.spearman_core_validate_logger <- function(logger) {
  if (!is.null(logger) && !is.function(logger)) {
    stop("logger must be NULL or a function.", call. = FALSE)
  }
  invisible(TRUE)
}

.spearman_core_log <- function(logger, level = "INFO", ...) {
  if (!is.null(logger)) {
    logger(level, ...)
  }
  invisible(NULL)
}

.spearman_core_warn <- function(logger, ...) {
  text <- paste0(..., collapse = "")
  warning(text, call. = FALSE)
  invisible(NULL)
}

.spearman_core_validate_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be a single TRUE/FALSE value.", call. = FALSE)
  }
  invisible(TRUE)
}

.spearman_core_validate_min_n <- function(min_n) {
  if (!is.numeric(min_n) ||
      length(min_n) != 1L ||
      is.na(min_n) ||
      !is.finite(min_n) ||
      min_n != floor(min_n) ||
      min_n < 3) {
    stop(
      "min_n must be one finite integer greater than or equal to 3; received: ",
      paste(min_n, collapse = ", "),
      call. = FALSE
    )
  }
  as.integer(min_n)
}

.spearman_core_scalar_name <- function(x, name, allow_null = FALSE) {
  if (allow_null && is.null(x)) {
    return(NULL)
  }
  if (length(x) != 1L || is.na(x)) {
    stop(name, " must be one non-empty column name.", call. = FALSE)
  }
  value <- trimws(as.character(x))
  if (!nzchar(value)) {
    stop(name, " must be one non-empty column name.", call. = FALSE)
  }
  value
}

.spearman_core_trim_columns <- function(x, source_label) {
  trimmed_names <- trimws(as.character(colnames(x)))

  if (length(trimmed_names) != ncol(x) ||
      any(is.na(trimmed_names)) ||
      any(!nzchar(trimmed_names))) {
    stop(
      source_label,
      " contains an NA or empty column name after trimming whitespace.",
      call. = FALSE
    )
  }

  duplicated_names <- unique(trimmed_names[duplicated(trimmed_names)])
  if (length(duplicated_names) > 0L) {
    stop(
      source_label,
      " contains duplicate column names after trimming whitespace: ",
      .spearman_core_preview(duplicated_names),
      call. = FALSE
    )
  }

  colnames(x) <- trimmed_names
  x
}

.spearman_core_character_selection <- function(x, name) {
  values <- trimws(as.character(x))
  if (length(values) == 0L ||
      any(is.na(values)) ||
      any(!nzchar(values))) {
    stop(name, " must not contain NA or empty values.", call. = FALSE)
  }
  duplicated_values <- unique(values[duplicated(values)])
  if (length(duplicated_values) > 0L) {
    stop(
      name,
      " contains duplicate values after trimming whitespace: ",
      .spearman_core_preview(duplicated_values),
      call. = FALSE
    )
  }
  values
}

.spearman_core_warn_degenerate_inputs <- function(
    feature_mat,
    target_df,
    target_cols,
    feature_label,
    metadata_label,
    logger
) {
  feature_all_missing <- character(0)
  feature_constant <- character(0)

  for (i in seq_len(nrow(feature_mat))) {
    values <- as.numeric(feature_mat[i, ])
    finite_values <- values[!is.na(values) & is.finite(values)]
    if (length(finite_values) == 0L) {
      feature_all_missing <- c(feature_all_missing, rownames(feature_mat)[i])
    } else if (length(unique(finite_values)) < 2L) {
      feature_constant <- c(feature_constant, rownames(feature_mat)[i])
    }
  }

  if (length(feature_all_missing) > 0L) {
    .spearman_core_warn(
      logger,
      feature_label,
      " has features with no finite matched values: ",
      .spearman_core_preview(feature_all_missing)
    )
  }
  if (length(feature_constant) > 0L) {
    .spearman_core_warn(
      logger,
      feature_label,
      " has constant features across their finite matched values: ",
      .spearman_core_preview(feature_constant)
    )
  }

  target_all_missing <- character(0)
  target_constant <- character(0)
  for (target_name in target_cols) {
    values <- target_df[[target_name]]
    finite_values <- values[!is.na(values) & is.finite(values)]
    if (length(finite_values) == 0L) {
      target_all_missing <- c(target_all_missing, target_name)
    } else if (length(unique(finite_values)) < 2L) {
      target_constant <- c(target_constant, target_name)
    }
  }

  if (length(target_all_missing) > 0L) {
    .spearman_core_warn(
      logger,
      metadata_label,
      " has targets with no finite matched values: ",
      .spearman_core_preview(target_all_missing)
    )
  }
  if (length(target_constant) > 0L) {
    .spearman_core_warn(
      logger,
      metadata_label,
      " has constant targets across their finite matched values: ",
      .spearman_core_preview(target_constant)
    )
  }

  invisible(NULL)
}

#' Validate, select, and align inputs for Spearman analysis.
#'
#' Relative paths are resolved before this function is called. The path
#' arguments are retained solely to make validation messages actionable.
prepare_spearman_data <- function(
    feature_df,
    target_df,
    feature_path,
    metadata_path,
    feature_col,
    sample_col,
    group_col,
    selected_features,
    target_cols,
    run_group_analysis,
    selected_groups,
    min_n,
    strict_sample_match,
    logger
) {
  .spearman_core_validate_logger(logger)
  .spearman_core_validate_flag(run_group_analysis, "run_group_analysis")
  .spearman_core_validate_flag(strict_sample_match, "strict_sample_match")
  min_n <- .spearman_core_validate_min_n(min_n)

  feature_label <- .spearman_core_label(feature_path, "expression matrix")
  metadata_label <- .spearman_core_label(metadata_path, "score/metadata table")
  feature_col <- .spearman_core_scalar_name(feature_col, "feature_col")
  sample_col <- .spearman_core_scalar_name(sample_col, "sample_col")
  group_col <- .spearman_core_scalar_name(
    group_col,
    "group_col",
    allow_null = !run_group_analysis
  )

  if (!is.data.frame(feature_df)) {
    stop(feature_label, " must be read as a data.frame.", call. = FALSE)
  }
  if (!is.data.frame(target_df)) {
    stop(metadata_label, " must be read as a data.frame.", call. = FALSE)
  }

  feature_df <- as.data.frame(
    feature_df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  target_df <- as.data.frame(
    target_df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  feature_df <- .spearman_core_trim_columns(feature_df, feature_label)
  target_df <- .spearman_core_trim_columns(target_df, metadata_label)

  if (!feature_col %in% colnames(feature_df)) {
    stop(
      feature_label,
      " does not contain feature_col '",
      feature_col,
      "'. Available columns: ",
      .spearman_core_preview(colnames(feature_df)),
      call. = FALSE
    )
  }
  if (!sample_col %in% colnames(target_df)) {
    stop(
      metadata_label,
      " does not contain sample_col '",
      sample_col,
      "'. Available columns: ",
      .spearman_core_preview(colnames(target_df)),
      call. = FALSE
    )
  }
  if (run_group_analysis && !group_col %in% colnames(target_df)) {
    stop(
      metadata_label,
      " does not contain group_col '",
      group_col,
      "'. Available columns: ",
      .spearman_core_preview(colnames(target_df)),
      call. = FALSE
    )
  }

  feature_df[[feature_col]] <- trimws(as.character(feature_df[[feature_col]]))
  target_df[[sample_col]] <- trimws(as.character(target_df[[sample_col]]))
  if (!is.null(group_col) && group_col %in% colnames(target_df)) {
    target_df[[group_col]] <- trimws(as.character(target_df[[group_col]]))
  }

  if (nrow(feature_df) == 0L) {
    stop(feature_label, " contains zero feature rows.", call. = FALSE)
  }
  if (any(is.na(feature_df[[feature_col]])) ||
      any(!nzchar(feature_df[[feature_col]]))) {
    stop(
      feature_label,
      " column '",
      feature_col,
      "' contains NA or empty feature identifiers after trimming.",
      call. = FALSE
    )
  }
  duplicated_features <- unique(
    feature_df[[feature_col]][duplicated(feature_df[[feature_col]])]
  )
  if (length(duplicated_features) > 0L) {
    stop(
      feature_label,
      " column '",
      feature_col,
      "' contains duplicate feature identifiers after trimming: ",
      .spearman_core_preview(duplicated_features),
      call. = FALSE
    )
  }

  if (nrow(target_df) == 0L) {
    stop(metadata_label, " contains zero sample rows.", call. = FALSE)
  }
  if (any(is.na(target_df[[sample_col]])) ||
      any(!nzchar(target_df[[sample_col]]))) {
    stop(
      metadata_label,
      " column '",
      sample_col,
      "' contains NA or empty sample identifiers after trimming.",
      call. = FALSE
    )
  }
  duplicated_samples <- unique(
    target_df[[sample_col]][duplicated(target_df[[sample_col]])]
  )
  if (length(duplicated_samples) > 0L) {
    stop(
      metadata_label,
      " column '",
      sample_col,
      "' contains duplicate sample identifiers after trimming: ",
      .spearman_core_preview(duplicated_samples),
      call. = FALSE
    )
  }

  feature_samples <- setdiff(colnames(feature_df), feature_col)
  if (length(feature_samples) == 0L) {
    stop(
      feature_label,
      " has no sample columns besides feature_col '",
      feature_col,
      "'.",
      call. = FALSE
    )
  }

  feature_numeric <- vapply(
    feature_df[, feature_samples, drop = FALSE],
    is.numeric,
    logical(1)
  )
  if (!all(feature_numeric)) {
    stop(
      feature_label,
      " has non-numeric sample columns: ",
      .spearman_core_preview(names(feature_numeric)[!feature_numeric]),
      call. = FALSE
    )
  }

  excluded_targets <- sample_col
  if (!is.null(group_col) && group_col %in% colnames(target_df)) {
    excluded_targets <- c(excluded_targets, group_col)
  }

  if (is.null(target_cols)) {
    candidate_targets <- setdiff(colnames(target_df), excluded_targets)
    candidate_numeric <- vapply(
      target_df[, candidate_targets, drop = FALSE],
      is.numeric,
      logical(1)
    )
    target_cols <- candidate_targets[candidate_numeric]
    if (length(target_cols) == 0L) {
      stop(
        metadata_label,
        " has no eligible numeric target columns after excluding '",
        paste(excluded_targets, collapse = "', '"),
        "'. Set target_cols explicitly or correct the input types.",
        call. = FALSE
      )
    }
    .spearman_core_log(
      logger,
      "INFO",
      "Automatically selected numeric targets in metadata order: ",
      paste(target_cols, collapse = ", ")
    )
  } else {
    target_cols <- .spearman_core_character_selection(target_cols, "target_cols")
    missing_targets <- setdiff(target_cols, colnames(target_df))
    if (length(missing_targets) > 0L) {
      stop(
        metadata_label,
        " is missing configured target columns: ",
        .spearman_core_preview(missing_targets),
        call. = FALSE
      )
    }
    invalid_targets <- intersect(target_cols, excluded_targets)
    if (length(invalid_targets) > 0L) {
      stop(
        "These configured target_cols are reserved identifier/group columns: ",
        .spearman_core_preview(invalid_targets),
        call. = FALSE
      )
    }
  }

  target_numeric <- vapply(
    target_df[, target_cols, drop = FALSE],
    is.numeric,
    logical(1)
  )
  if (!all(target_numeric)) {
    stop(
      metadata_label,
      " has configured target columns that are not numeric: ",
      .spearman_core_preview(names(target_numeric)[!target_numeric]),
      call. = FALSE
    )
  }

  if (is.null(selected_features)) {
    feature_df_selected <- feature_df
  } else {
    selected_features <- .spearman_core_character_selection(
      selected_features,
      "selected_features"
    )
    missing_features <- setdiff(selected_features, feature_df[[feature_col]])
    if (length(missing_features) > 0L) {
      stop(
        feature_label,
        " does not contain configured selected_features: ",
        .spearman_core_preview(missing_features),
        call. = FALSE
      )
    }
    feature_df_selected <- feature_df[
      match(selected_features, feature_df[[feature_col]]),
      ,
      drop = FALSE
    ]
  }
  if (nrow(feature_df_selected) == 0L) {
    stop("Feature selection produced zero features.", call. = FALSE)
  }

  groups_used <- character(0)
  if (run_group_analysis) {
    group_values <- target_df[[group_col]]
    missing_group <- is.na(group_values) | !nzchar(group_values)
    if (any(missing_group)) {
      .spearman_core_warn(
        logger,
        metadata_label,
        " column '",
        group_col,
        "' has ",
        sum(missing_group),
        " missing/empty group value(s); those samples are excluded from group analyses."
      )
    }

    available_groups <- unique(group_values[!missing_group])
    if (length(available_groups) == 0L) {
      stop(
        metadata_label,
        " column '",
        group_col,
        "' has no usable group values.",
        call. = FALSE
      )
    }

    if (is.null(selected_groups)) {
      groups_used <- available_groups
    } else {
      groups_used <- .spearman_core_character_selection(
        selected_groups,
        "selected_groups"
      )
      missing_groups <- setdiff(groups_used, available_groups)
      if (length(missing_groups) > 0L) {
        stop(
          metadata_label,
          " column '",
          group_col,
          "' does not contain configured selected_groups: ",
          .spearman_core_preview(missing_groups),
          call. = FALSE
        )
      }
    }

  }

  target_samples <- target_df[[sample_col]]
  only_in_target <- setdiff(target_samples, feature_samples)
  only_in_feature <- setdiff(feature_samples, target_samples)

  .spearman_core_log(
    logger,
    "INFO",
    "Sample matching: ",
    length(feature_samples),
    " expression sample(s), ",
    length(target_samples),
    " metadata sample(s), ",
    length(only_in_feature),
    " expression-only and ",
    length(only_in_target),
    " metadata-only."
  )

  if (strict_sample_match &&
      (length(only_in_target) > 0L || length(only_in_feature) > 0L)) {
    details <- c(
      if (length(only_in_feature) > 0L) {
        paste0("expression-only: ", .spearman_core_preview(only_in_feature))
      },
      if (length(only_in_target) > 0L) {
        paste0("metadata-only: ", .spearman_core_preview(only_in_target))
      }
    )
    stop(
      "strict_sample_match is TRUE, but samples do not match between ",
      feature_label,
      " and ",
      metadata_label,
      ". ",
      paste(details, collapse = "; "),
      call. = FALSE
    )
  }

  # Preserve the expression-matrix column order, matching base::intersect()
  # in the canonical implementation.
  common_samples <- feature_samples[feature_samples %in% target_samples]
  if (length(common_samples) < min_n) {
    stop(
      "Only ",
      length(common_samples),
      " matched sample(s) were found between ",
      feature_label,
      " and ",
      metadata_label,
      "; min_n is ",
      min_n,
      ".",
      call. = FALSE
    )
  }

  feature_df_aligned <- feature_df_selected[
    ,
    c(feature_col, common_samples),
    drop = FALSE
  ]
  target_row_index <- match(common_samples, target_df[[sample_col]])
  if (anyNA(target_row_index)) {
    stop(
      "Internal sample alignment failed for ",
      metadata_label,
      "; match() returned NA.",
      call. = FALSE
    )
  }
  target_df_aligned <- target_df[target_row_index, , drop = FALSE]
  sample_order_ok <- identical(
    common_samples,
    target_df_aligned[[sample_col]]
  )
  if (!sample_order_ok) {
    stop(
      "Internal sample alignment failed: expression and metadata sample order differ.",
      call. = FALSE
    )
  }

  feature_mat <- as.matrix(
    feature_df_aligned[, common_samples, drop = FALSE]
  )
  storage.mode(feature_mat) <- "double"
  rownames(feature_mat) <- feature_df_aligned[[feature_col]]

  .spearman_core_warn_degenerate_inputs(
    feature_mat = feature_mat,
    target_df = target_df_aligned,
    target_cols = target_cols,
    feature_label = feature_label,
    metadata_label = metadata_label,
    logger = logger
  )

  .spearman_core_log(
    logger,
    "INFO",
    "Prepared ",
    nrow(feature_mat),
    " feature(s), ",
    length(target_cols),
    " target(s), and ",
    length(common_samples),
    " matched sample(s)."
  )

  prepared <- list(
    feature_df_aligned = feature_df_aligned,
    target_df_aligned = target_df_aligned,
    feature_mat = feature_mat,
    features_used = feature_df_aligned[[feature_col]],
    targets_used = target_cols,
    groups_used = groups_used,
    sample_report = list(
      only_in_target = only_in_target,
      only_in_feature = only_in_feature,
      common_samples = common_samples,
      sample_order_ok = sample_order_ok
    ),
    feature_col = feature_col,
    sample_col = sample_col,
    group_col = group_col,
    source_paths = list(
      feature = feature_path,
      metadata = metadata_path
    ),
    parameters = list(
      feature_col = feature_col,
      sample_col = sample_col,
      group_col = group_col,
      selected_features = selected_features,
      target_cols = target_cols,
      selected_groups = groups_used,
      run_group_analysis = run_group_analysis,
      min_n = min_n,
      strict_sample_match = strict_sample_match
    )
  )
  class(prepared) <- c("spearman_prepared_data", "list")
  prepared
}

.spearman_core_one_target_group <- function(
    feature_mat,
    target_sub,
    target_col,
    sample_col,
    group_label,
    min_n,
    correlation_method
) {
  samples_use <- as.character(target_sub[[sample_col]])
  missing_matrix_samples <- setdiff(samples_use, colnames(feature_mat))
  if (length(missing_matrix_samples) > 0L) {
    stop(
      "Prepared feature matrix is missing samples required for group '",
      group_label,
      "', target '",
      target_col,
      "': ",
      .spearman_core_preview(missing_matrix_samples),
      call. = FALSE
    )
  }

  feature_sub <- feature_mat[, samples_use, drop = FALSE]
  if (!identical(colnames(feature_sub), samples_use)) {
    stop(
      "Internal sample order differs for group '",
      group_label,
      "', target '",
      target_col,
      "'.",
      call. = FALSE
    )
  }

  y <- target_sub[[target_col]]
  n_features <- nrow(feature_sub)
  feature_name <- rownames(feature_sub)
  rho <- rep(NA_real_, n_features)
  pvalue <- rep(NA_real_, n_features)
  n_complete <- integer(n_features)
  status <- rep(NA_character_, n_features)

  for (i in seq_len(n_features)) {
    x <- as.numeric(feature_sub[i, ])
    complete_idx <- stats::complete.cases(x, y) &
      is.finite(x) &
      is.finite(y)
    x_complete <- x[complete_idx]
    y_complete <- y[complete_idx]
    n_complete[i] <- sum(complete_idx)

    # Preserve the canonical status precedence exactly.
    if (n_complete[i] < min_n) {
      status[i] <- "insufficient_n"
      next
    }
    if (length(unique(x_complete)) < 2L) {
      status[i] <- "constant_feature"
      next
    }
    if (length(unique(y_complete)) < 2L) {
      status[i] <- "constant_target"
      next
    }

    test_result <- tryCatch(
      suppressWarnings(
        stats::cor.test(
          x = x_complete,
          y = y_complete,
          method = correlation_method,
          exact = FALSE
        )
      ),
      error = function(e) e
    )
    if (inherits(test_result, "error")) {
      status[i] <- "calculation_error"
      next
    }

    rho[i] <- unname(test_result$estimate)
    pvalue[i] <- test_result$p.value
    status[i] <- "OK"
  }

  result_df <- data.frame(
    feature = feature_name,
    target = target_col,
    group = group_label,
    rho = rho,
    abs_rho = abs(rho),
    pvalue = pvalue,
    padj = NA_real_,
    n = n_complete,
    direction = ifelse(
      is.na(rho),
      NA_character_,
      ifelse(rho > 0, "Positive", ifelse(rho < 0, "Negative", "Zero"))
    ),
    status = status,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  valid_idx <- which(result_df$status == "OK" & !is.na(result_df$pvalue))
  if (length(valid_idx) > 0L) {
    # Adjustment is intentionally local to this target x group result block.
    result_df$padj[valid_idx] <- stats::p.adjust(
      result_df$pvalue[valid_idx],
      method = "BH"
    )
  }

  result_df <- result_df[
    order(is.na(result_df$pvalue), result_df$pvalue),
    ,
    drop = FALSE
  ]
  rownames(result_df) <- NULL
  result_df
}

.spearman_core_validate_prepared <- function(prepared) {
  required_fields <- c(
    "feature_df_aligned",
    "target_df_aligned",
    "feature_mat",
    "features_used",
    "targets_used",
    "groups_used",
    "sample_report",
    "sample_col",
    "group_col"
  )
  if (!is.list(prepared)) {
    stop("prepared must be the list returned by prepare_spearman_data().", call. = FALSE)
  }
  missing_fields <- setdiff(required_fields, names(prepared))
  if (length(missing_fields) > 0L) {
    stop(
      "prepared is missing fields: ",
      paste(missing_fields, collapse = ", "),
      ". Re-run prepare_spearman_data().",
      call. = FALSE
    )
  }
  if (!is.matrix(prepared$feature_mat) || !is.numeric(prepared$feature_mat)) {
    stop("prepared$feature_mat must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(prepared$feature_mat) == 0L) {
    stop("prepared$feature_mat contains zero features.", call. = FALSE)
  }
  if (!is.data.frame(prepared$target_df_aligned)) {
    stop("prepared$target_df_aligned must be a data.frame.", call. = FALSE)
  }
  if (!is.character(prepared$sample_col) || length(prepared$sample_col) != 1L ||
      is.na(prepared$sample_col) || !nzchar(prepared$sample_col)) {
    stop("prepared$sample_col is invalid; re-run data preparation.", call. = FALSE)
  }
  if (!prepared$sample_col %in% colnames(prepared$target_df_aligned)) {
    stop(
      "prepared$target_df_aligned is missing sample column '",
      prepared$sample_col,
      "'.",
      call. = FALSE
    )
  }
  target_samples <- as.character(
    prepared$target_df_aligned[[prepared$sample_col]]
  )
  if (!identical(colnames(prepared$feature_mat), target_samples)) {
    stop(
      "Prepared expression and metadata sample order no longer matches; ",
      "re-run data preparation.",
      call. = FALSE
    )
  }
  if (!identical(rownames(prepared$feature_mat), as.character(prepared$features_used))) {
    stop(
      "prepared$features_used does not match feature_mat row order; ",
      "re-run data preparation.",
      call. = FALSE
    )
  }
  if (length(prepared$targets_used) == 0L) {
    stop("prepared$targets_used is empty.", call. = FALSE)
  }
  missing_targets <- setdiff(
    prepared$targets_used,
    colnames(prepared$target_df_aligned)
  )
  if (length(missing_targets) > 0L) {
    stop(
      "Prepared metadata is missing target columns: ",
      .spearman_core_preview(missing_targets),
      call. = FALSE
    )
  }
  target_numeric <- vapply(
    prepared$target_df_aligned[, prepared$targets_used, drop = FALSE],
    is.numeric,
    logical(1)
  )
  if (!all(target_numeric)) {
    stop(
      "Prepared target columns are not numeric: ",
      .spearman_core_preview(names(target_numeric)[!target_numeric]),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Run canonical Spearman correlations on prepared inputs.
#'
#' @return A data.frame with the ten canonical raw result columns.
run_spearman_correlations <- function(
    prepared,
    run_all_samples,
    run_group_analysis,
    min_n,
    correlation_method = "spearman",
    logger
) {
  .spearman_core_validate_logger(logger)
  .spearman_core_validate_flag(run_all_samples, "run_all_samples")
  .spearman_core_validate_flag(run_group_analysis, "run_group_analysis")
  min_n <- .spearman_core_validate_min_n(min_n)
  .spearman_core_validate_prepared(prepared)

  if (!run_all_samples && !run_group_analysis) {
    stop(
      "run_all_samples and run_group_analysis cannot both be FALSE.",
      call. = FALSE
    )
  }
  if (!is.character(correlation_method) ||
      length(correlation_method) != 1L ||
      is.na(correlation_method) ||
      !nzchar(trimws(correlation_method))) {
    stop("correlation_method must be one non-empty string.", call. = FALSE)
  }
  correlation_method <- tolower(trimws(correlation_method))
  if (!identical(correlation_method, "spearman")) {
    stop(
      "correlation_method must be 'spearman' to preserve this workflow's ",
      "statistical definition; received '",
      correlation_method,
      "'.",
      call. = FALSE
    )
  }

  feature_mat <- prepared$feature_mat
  target_df_aligned <- prepared$target_df_aligned
  target_cols <- as.character(prepared$targets_used)
  groups_used <- as.character(prepared$groups_used)
  sample_col <- prepared$sample_col
  group_col <- prepared$group_col

  if (ncol(feature_mat) < min_n) {
    stop(
      "Prepared data has only ",
      ncol(feature_mat),
      " matched sample(s), fewer than min_n = ",
      min_n,
      ".",
      call. = FALSE
    )
  }
  if (run_group_analysis) {
    if (is.null(group_col) ||
        length(group_col) != 1L ||
        is.na(group_col) ||
        !group_col %in% colnames(target_df_aligned)) {
      stop(
        "run_group_analysis is TRUE, but prepared data has no usable group_col.",
        call. = FALSE
      )
    }
    if (length(groups_used) == 0L) {
      stop(
        "run_group_analysis is TRUE, but prepared$groups_used is empty.",
        call. = FALSE
      )
    }
  }

  if (
    run_all_samples &&
      run_group_analysis &&
      "All" %in% groups_used
  ) {
    .spearman_core_warn(
      logger,
      "Metadata contains a real group named 'All'. For backward ",
      "compatibility it will be analyzed, but its result label collides ",
      "with the all-sample label; interpret and export those rows carefully."
    )
  }

  result_list <- list()
  result_index <- 1L

  if (run_all_samples) {
    .spearman_core_log(logger, "INFO", "Running all-sample correlations.")
    for (target_name in target_cols) {
      result_list[[result_index]] <- .spearman_core_one_target_group(
        feature_mat = feature_mat,
        target_sub = target_df_aligned,
        target_col = target_name,
        sample_col = sample_col,
        group_label = "All",
        min_n = min_n,
        correlation_method = correlation_method
      )
      result_index <- result_index + 1L
    }
  }

  if (run_group_analysis) {
    .spearman_core_log(logger, "INFO", "Running within-group correlations.")
    for (group_name in groups_used) {
      group_vector <- target_df_aligned[[group_col]]
      group_index <- !is.na(group_vector) & group_vector == group_name
      target_group <- target_df_aligned[group_index, , drop = FALSE]

      .spearman_core_log(
        logger,
        "INFO",
        "Group '",
        group_name,
        "' has ",
        nrow(target_group),
        " matched sample(s)."
      )
      if (nrow(target_group) < min_n) {
        .spearman_core_warn(
          logger,
          "Group '",
          group_name,
          "' has ",
          nrow(target_group),
          " matched sample(s), fewer than min_n = ",
          min_n,
          "; skipping this group."
        )
        next
      }

      for (target_name in target_cols) {
        result_list[[result_index]] <- .spearman_core_one_target_group(
          feature_mat = feature_mat,
          target_sub = target_group,
          target_col = target_name,
          sample_col = sample_col,
          group_label = group_name,
          min_n = min_n,
          correlation_method = correlation_method
        )
        result_index <- result_index + 1L
      }
    }
  }

  if (length(result_list) == 0L) {
    stop(
      "No correlation result blocks were produced. Check run switches, ",
      "selected_groups, and min_n.",
      call. = FALSE
    )
  }

  cor_results <- do.call(rbind, result_list)
  rownames(cor_results) <- NULL
  canonical_columns <- c(
    "feature",
    "target",
    "group",
    "rho",
    "abs_rho",
    "pvalue",
    "padj",
    "n",
    "direction",
    "status"
  )
  cor_results <- cor_results[, canonical_columns, drop = FALSE]

  .spearman_core_log(
    logger,
    "INFO",
    "Completed ",
    length(result_list),
    " target-by-group result block(s) and ",
    nrow(cor_results),
    " feature correlation row(s)."
  )
  cor_results
}
