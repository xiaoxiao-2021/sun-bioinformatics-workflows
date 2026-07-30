# Runtime helpers shared by the command-line workflow.

.spearman_create_logger <- function(log_file, verbose = TRUE) {
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("logging.verbose must be a single TRUE/FALSE value.", call. = FALSE)
  }

  log_dir <- dirname(log_file)
  if (!dir.exists(log_dir) &&
      !dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop(
      "Unable to create log directory: ",
      log_dir,
      call. = FALSE
    )
  }

  function(level = "INFO", ...) {
    level <- toupper(as.character(level)[1L])
    text <- paste0(..., collapse = "")
    line <- sprintf(
      "[%s] [%s] %s",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
      level,
      text
    )

    cat(line, "\n", file = log_file, append = TRUE, sep = "")

    if (isTRUE(verbose) || level %in% c("WARN", "ERROR")) {
      message(line)
    }

    invisible(line)
  }
}

.spearman_safe_filename <- function(x) {
  value <- gsub(
    pattern = "[^A-Za-z0-9._-]+",
    replacement = "_",
    x = as.character(x)
  )
  value[value == ""] <- "unnamed"
  value
}

.spearman_assert_unique_filenames <- function(labels, context) {
  safe <- .spearman_safe_filename(labels)
  duplicated_safe <- unique(safe[duplicated(safe)])

  if (length(duplicated_safe) > 0L) {
    collisions <- vapply(
      duplicated_safe,
      function(one_safe) {
        paste(labels[safe == one_safe], collapse = " / ")
      },
      character(1)
    )
    stop(
      "Filename collision in ",
      context,
      " after sanitizing labels: ",
      paste(collisions, collapse = "; "),
      ". Rename the configured target/group labels.",
      call. = FALSE
    )
  }

  invisible(safe)
}

.spearman_write_csv <- function(x, path) {
  parent <- dirname(path)
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create output directory: ", parent, call. = FALSE)
  }

  utils::write.csv(
    x = x,
    file = path,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )

  invisible(path)
}

.spearman_format_threshold <- function(x) {
  as.character(x)
}

.spearman_config_signature <- function(config, input_md5) {
  list(
    schema_version = 3L,
    workflow_code_md5 = .spearman_workflow_code_signatures,
    input_md5 = input_md5,
    paths = list(
      expression_matrix = config$inputs$expression_matrix,
      score_metadata = config$inputs$score_metadata
    ),
    input_formats = list(
      expression = config$inputs$expression_format,
      metadata = config$inputs$metadata_format
    ),
    columns = config$columns,
    selection = list(
      selected_features = config$analysis$selected_features,
      target_columns = config$analysis$target_columns,
      selected_groups = config$analysis$selected_groups
    ),
    analysis = list(
      run_all_samples = config$analysis$run_all_samples,
      run_group_analysis = config$analysis$run_group_analysis,
      correlation_method = config$analysis$correlation_method,
      min_n = config$analysis$min_n,
      strict_sample_match = config$analysis$strict_sample_match
    )
  )
}

.spearman_stage_signature <- function(signature, stage) {
  code_signatures <- signature$workflow_code_md5
  if (
    is.null(names(code_signatures)) ||
      !stage %in% names(code_signatures) ||
      length(code_signatures[[stage]]) != 1L ||
      is.na(code_signatures[[stage]]) ||
      !nzchar(code_signatures[[stage]])
  ) {
    stop(
      "Missing workflow code signature for intermediate stage: ",
      stage,
      call. = FALSE
    )
  }

  fields <- if (identical(stage, "inputs")) {
    c(
      "schema_version",
      "input_md5",
      "paths",
      "input_formats"
    )
  } else if (identical(stage, "prepared")) {
    c(
      "schema_version",
      "input_md5",
      "paths",
      "input_formats",
      "columns",
      "selection",
      "analysis"
    )
  } else if (identical(stage, "correlations")) {
    setdiff(names(signature), "workflow_code_md5")
  } else {
    stop("Unknown intermediate stage: ", stage, call. = FALSE)
  }

  append(
    signature[fields],
    list(
      workflow_code_md5 = unname(code_signatures[[stage]])
    ),
    after = 1L
  )
}

.spearman_load_checkpoint <- function(
    path,
    expected_signature,
    stage,
    reuse,
    logger
) {
  if (!isTRUE(reuse) || !file.exists(path)) {
    return(NULL)
  }

  checkpoint <- tryCatch(
    readRDS(path),
    error = function(e) e
  )

  if (inherits(checkpoint, "error")) {
    logger(
      "WARN",
      "Ignoring unreadable ",
      stage,
      " checkpoint ",
      path,
      ": ",
      conditionMessage(checkpoint)
    )
    return(NULL)
  }

  if (!is.list(checkpoint) ||
      is.null(checkpoint$signature) ||
      is.null(checkpoint$value) ||
      !identical(checkpoint$signature, expected_signature)) {
    logger(
      "INFO",
      "Checkpoint is stale for stage '",
      stage,
      "'; recomputing."
    )
    return(NULL)
  }

  logger("INFO", "Reused checkpoint: ", path)
  checkpoint$value
}

.spearman_save_checkpoint <- function(
    value,
    signature,
    path,
    stage,
    enabled,
    logger
) {
  if (!isTRUE(enabled)) {
    return(invisible(NULL))
  }

  parent <- dirname(path)
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop(
      "Unable to create intermediate directory: ",
      parent,
      call. = FALSE
    )
  }

  tmp <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = parent
  )
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(
    list(
      stage = stage,
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      signature = signature,
      value = value
    ),
    file = tmp,
    version = 3L
  )

  if (file.exists(path) && !file.remove(path)) {
    stop("Unable to replace checkpoint: ", path, call. = FALSE)
  }
  if (!file.rename(tmp, path)) {
    stop("Unable to finalize checkpoint: ", path, call. = FALSE)
  }

  logger("INFO", "Saved checkpoint: ", path)
  invisible(path)
}
