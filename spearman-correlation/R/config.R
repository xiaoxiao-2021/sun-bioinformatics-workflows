# Configuration loading and validation for the Spearman workflow.

.spearman_config_error <- function(field, message) {
  stop(
    "Invalid configuration field `",
    field,
    "`: ",
    message,
    call. = FALSE
  )
}

.spearman_config_schema <- function() {
  list(
    project_id = TRUE,
    inputs = list(
      expression_matrix = TRUE,
      score_metadata = TRUE,
      expression_format = TRUE,
      metadata_format = TRUE
    ),
    output_dir = TRUE,
    columns = list(
      feature = TRUE,
      sample = TRUE,
      group = TRUE
    ),
    analysis = list(
      selected_features = TRUE,
      target_columns = TRUE,
      run_all_samples = TRUE,
      run_group_analysis = TRUE,
      selected_groups = TRUE,
      correlation_method = TRUE,
      min_n = TRUE,
      strict_sample_match = TRUE
    ),
    thresholds = list(
      r_cutoff = TRUE,
      p_cutoff = TRUE,
      padj_cutoff = TRUE
    ),
    filters = list(
      include_positive = TRUE,
      include_negative = TRUE
    ),
    quadrant = list(
      enabled = TRUE,
      x_target = TRUE,
      x_group = TRUE,
      y_target = TRUE,
      y_group = TRUE,
      significance = TRUE,
      require_significant_both = TRUE
    ),
    outputs = list(
      aligned_inputs = TRUE,
      combined_tables = TRUE,
      per_analysis_tables = TRUE,
      filtered_tables = TRUE,
      result_rds = TRUE
    ),
    intermediate = list(
      enabled = TRUE,
      reuse_existing = TRUE
    ),
    plots = list(
      overview = list(
        enabled = TRUE,
        label_top_n_each = TRUE,
        width = TRUE,
        height = TRUE,
        dpi = TRUE,
        formats = TRUE
      ),
      single_feature = list(
        enabled = TRUE,
        width = TRUE,
        height = TRUE,
        dpi = TRUE,
        formats = TRUE,
        items = list(
          .item = list(
            feature = TRUE,
            target = TRUE,
            group = TRUE,
            type = TRUE,
            add_lm = TRUE
          )
        )
      ),
      quadrant = list(
        enabled = TRUE,
        width = TRUE,
        height = TRUE,
        dpi = TRUE,
        formats = TRUE
      )
    ),
    validation = list(
      manual = list(
        enabled = TRUE,
        items = list(
          .item = list(
            feature = TRUE,
            target = TRUE,
            group = TRUE
          )
        )
      )
    ),
    logging = list(
      verbose = TRUE
    )
  )
}

.spearman_config_defaults <- function() {
  list(
    project_id = NULL,
    inputs = list(
      expression_matrix = NULL,
      score_metadata = NULL,
      expression_format = "csv",
      metadata_format = "csv"
    ),
    output_dir = NULL,
    columns = list(
      feature = "feature",
      sample = "sample",
      group = NULL
    ),
    analysis = list(
      selected_features = NULL,
      target_columns = NULL,
      run_all_samples = TRUE,
      run_group_analysis = FALSE,
      selected_groups = NULL,
      correlation_method = "spearman",
      min_n = 5L,
      strict_sample_match = FALSE
    ),
    thresholds = list(
      r_cutoff = 0.5,
      p_cutoff = 0.05,
      padj_cutoff = 0.05
    ),
    filters = list(
      include_positive = TRUE,
      include_negative = TRUE
    ),
    quadrant = list(
      enabled = FALSE,
      x_target = NULL,
      x_group = NULL,
      y_target = NULL,
      y_group = NULL,
      significance = "padj",
      require_significant_both = TRUE
    ),
    outputs = list(
      aligned_inputs = TRUE,
      combined_tables = TRUE,
      per_analysis_tables = TRUE,
      filtered_tables = TRUE,
      result_rds = TRUE
    ),
    intermediate = list(
      enabled = TRUE,
      reuse_existing = TRUE
    ),
    plots = list(
      overview = list(
        enabled = TRUE,
        label_top_n_each = 0L,
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
      manual = list(
        enabled = FALSE,
        items = list()
      )
    ),
    logging = list(
      verbose = TRUE
    )
  )
}

.spearman_check_unknown_keys <- function(value, schema, path) {
  if (!is.list(schema)) {
    return(invisible(NULL))
  }

  if (identical(names(schema), ".item")) {
    if (is.list(value)) {
      for (i in seq_along(value)) {
        .spearman_check_unknown_keys(
          value[[i]],
          schema$.item,
          paste0(path, "[", i, "]")
        )
      }
    }
    return(invisible(NULL))
  }

  if (!is.list(value)) {
    return(invisible(NULL))
  }

  value_names <- names(value)
  if (is.null(value_names) || any(!nzchar(value_names))) {
    return(invisible(NULL))
  }

  unknown <- setdiff(value_names, names(schema))
  if (length(unknown) > 0L) {
    stop(
      "Unknown configuration field",
      if (length(unknown) > 1L) "s" else "",
      " under `",
      path,
      "`: ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  nested <- intersect(value_names, names(schema))
  for (name in nested) {
    if (is.list(schema[[name]])) {
      .spearman_check_unknown_keys(
        value[[name]],
        schema[[name]],
        paste0(path, ".", name)
      )
    }
  }
  invisible(NULL)
}

.spearman_merge_config <- function(defaults, supplied) {
  result <- defaults
  for (name in names(supplied)) {
    default_value <- defaults[[name]]
    supplied_value <- supplied[[name]]

    if (
      is.list(default_value) &&
        !is.null(names(default_value)) &&
        is.list(supplied_value) &&
        !is.null(names(supplied_value))
    ) {
      result[[name]] <- .spearman_merge_config(
        default_value,
        supplied_value
      )
    } else {
      # Single-bracket assignment preserves an explicit YAML null as NULL.
      result[name] <- list(supplied_value)
    }
  }
  result
}

.spearman_assert_mapping <- function(x, field) {
  if (
    !is.list(x) ||
      is.null(names(x)) ||
      any(!nzchar(names(x)))
  ) {
    .spearman_config_error(field, "must be a YAML mapping.")
  }
  invisible(x)
}

.spearman_assert_string <- function(x, field) {
  if (
    !is.character(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(trimws(x))
  ) {
    .spearman_config_error(
      field,
      "must be one non-empty character string."
    )
  }
  x
}

.spearman_assert_optional_string <- function(x, field) {
  if (is.null(x)) {
    return(NULL)
  }
  .spearman_assert_string(x, field)
}

.spearman_assert_logical <- function(x, field) {
  if (
    !is.logical(x) ||
      length(x) != 1L ||
      is.na(x)
  ) {
    .spearman_config_error(
      field,
      "must be exactly `true` or `false`."
    )
  }
  x
}

.spearman_assert_number <- function(
    x,
    field,
    minimum = -Inf,
    maximum = Inf,
    minimum_inclusive = TRUE,
    maximum_inclusive = TRUE
) {
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x)
  ) {
    .spearman_config_error(
      field,
      "must be one finite numeric value."
    )
  }

  valid_minimum <- if (minimum_inclusive) {
    x >= minimum
  } else {
    x > minimum
  }
  valid_maximum <- if (maximum_inclusive) {
    x <= maximum
  } else {
    x < maximum
  }

  if (
    !isTRUE(valid_minimum) ||
      !isTRUE(valid_maximum)
  ) {
    left <- if (minimum_inclusive) "[" else "("
    right <- if (maximum_inclusive) "]" else ")"
    .spearman_config_error(
      field,
      paste0(
        "must be one finite number in ",
        left,
        minimum,
        ", ",
        maximum,
        right,
        "."
      )
    )
  }
  as.numeric(x)
}

.spearman_assert_integer <- function(
    x,
    field,
    minimum = -.Machine$integer.max,
    maximum = .Machine$integer.max
) {
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x != floor(x) ||
      x < minimum ||
      x > maximum
  ) {
    .spearman_config_error(
      field,
      paste0(
        "must be one integer in [",
        minimum,
        ", ",
        maximum,
        "]."
      )
    )
  }
  as.integer(x)
}

.spearman_assert_enum <- function(x, field, choices) {
  x <- .spearman_assert_string(x, field)
  if (!x %in% choices) {
    .spearman_config_error(
      field,
      paste0(
        "must be one of: ",
        paste(sprintf("`%s`", choices), collapse = ", "),
        "."
      )
    )
  }
  x
}

.spearman_assert_character_vector <- function(
    x,
    field,
    nullable = FALSE,
    allow_empty = FALSE
) {
  if (is.null(x)) {
    if (nullable) {
      return(NULL)
    }
    .spearman_config_error(field, "must not be null.")
  }

  if (is.list(x) && is.null(names(x))) {
    scalar_strings <- vapply(
      x,
      function(item) {
        is.character(item) &&
          length(item) == 1L &&
          !is.na(item)
      },
      logical(1)
    )
    if (!all(scalar_strings)) {
      .spearman_config_error(
        field,
        "must be a YAML sequence of character strings."
      )
    }
    x <- unlist(x, use.names = FALSE)
  }

  if (!is.character(x)) {
    .spearman_config_error(
      field,
      "must be a character string or YAML sequence of strings."
    )
  }
  if (!allow_empty && length(x) == 0L) {
    .spearman_config_error(field, "must contain at least one value.")
  }
  if (anyNA(x) || any(!nzchar(trimws(x)))) {
    .spearman_config_error(
      field,
      "must not contain null, missing, or empty values."
    )
  }
  if (anyDuplicated(x) > 0L) {
    duplicates <- unique(x[duplicated(x)])
    .spearman_config_error(
      field,
      paste0(
        "contains duplicate value(s): ",
        paste(duplicates, collapse = ", "),
        "."
      )
    )
  }
  unname(x)
}

.spearman_normalize_item_defaults <- function(items, defaults) {
  if (!is.list(items)) {
    return(items)
  }
  lapply(
    items,
    function(item) {
      if (
        is.list(item) &&
          !is.null(names(item)) &&
          all(nzchar(names(item)))
      ) {
        .spearman_merge_config(defaults, item)
      } else {
        item
      }
    }
  )
}

.spearman_assert_item_sequence <- function(
    items,
    field,
    require_items
) {
  if (!is.list(items)) {
    .spearman_config_error(
      field,
      "must be a YAML sequence of mappings."
    )
  }
  if (
    length(items) > 0L &&
      !is.null(names(items)) &&
      any(nzchar(names(items)))
  ) {
    .spearman_config_error(
      field,
      "must be a YAML sequence (each item begins with `-`)."
    )
  }
  if (isTRUE(require_items) && length(items) == 0L) {
    .spearman_config_error(
      field,
      "must contain at least one item when its module is enabled."
    )
  }
  invisible(items)
}

.spearman_validate_config <- function(config) {
  mapping_fields <- c(
    "inputs",
    "columns",
    "analysis",
    "thresholds",
    "filters",
    "quadrant",
    "outputs",
    "intermediate",
    "plots",
    "validation",
    "logging"
  )
  mapping_values <- list(
    config$inputs,
    config$columns,
    config$analysis,
    config$thresholds,
    config$filters,
    config$quadrant,
    config$outputs,
    config$intermediate,
    config$plots,
    config$validation,
    config$logging
  )
  for (i in seq_along(mapping_fields)) {
    .spearman_assert_mapping(
      mapping_values[[i]],
      mapping_fields[[i]]
    )
  }
  .spearman_assert_mapping(
    config$plots$overview,
    "plots.overview"
  )
  .spearman_assert_mapping(
    config$plots$single_feature,
    "plots.single_feature"
  )
  .spearman_assert_mapping(
    config$plots$quadrant,
    "plots.quadrant"
  )
  .spearman_assert_mapping(
    config$validation$manual,
    "validation.manual"
  )

  config$project_id <- .spearman_assert_string(
    config$project_id,
    "project_id"
  )
  config$inputs$expression_matrix <- .spearman_assert_string(
    config$inputs$expression_matrix,
    "inputs.expression_matrix"
  )
  config$inputs$score_metadata <- .spearman_assert_string(
    config$inputs$score_metadata,
    "inputs.score_metadata"
  )
  input_formats <- c("csv", "tsv", "rds")
  config$inputs$expression_format <- .spearman_assert_enum(
    config$inputs$expression_format,
    "inputs.expression_format",
    input_formats
  )
  config$inputs$metadata_format <- .spearman_assert_enum(
    config$inputs$metadata_format,
    "inputs.metadata_format",
    input_formats
  )
  config$output_dir <- .spearman_assert_string(
    config$output_dir,
    "output_dir"
  )

  config$columns$feature <- .spearman_assert_string(
    config$columns$feature,
    "columns.feature"
  )
  config$columns$sample <- .spearman_assert_string(
    config$columns$sample,
    "columns.sample"
  )
  config$columns["group"] <- list(
    .spearman_assert_optional_string(
      config$columns$group,
      "columns.group"
    )
  )

  config$analysis["selected_features"] <- list(
    .spearman_assert_character_vector(
      config$analysis$selected_features,
      "analysis.selected_features",
      nullable = TRUE
    )
  )
  config$analysis["target_columns"] <- list(
    .spearman_assert_character_vector(
      config$analysis$target_columns,
      "analysis.target_columns",
      nullable = TRUE
    )
  )
  config$analysis$run_all_samples <- .spearman_assert_logical(
    config$analysis$run_all_samples,
    "analysis.run_all_samples"
  )
  config$analysis$run_group_analysis <- .spearman_assert_logical(
    config$analysis$run_group_analysis,
    "analysis.run_group_analysis"
  )
  if (
    !config$analysis$run_all_samples &&
      !config$analysis$run_group_analysis
  ) {
    .spearman_config_error(
      "analysis",
      paste0(
        "at least one of `run_all_samples` or ",
        "`run_group_analysis` must be true."
      )
    )
  }
  if (
    config$analysis$run_group_analysis &&
      is.null(config$columns$group)
  ) {
    .spearman_config_error(
      "columns.group",
      "must be set when analysis.run_group_analysis is true."
    )
  }
  config$analysis["selected_groups"] <- list(
    .spearman_assert_character_vector(
      config$analysis$selected_groups,
      "analysis.selected_groups",
      nullable = TRUE
    )
  )
  config$analysis$correlation_method <- .spearman_assert_enum(
    config$analysis$correlation_method,
    "analysis.correlation_method",
    "spearman"
  )
  config$analysis$min_n <- .spearman_assert_integer(
    config$analysis$min_n,
    "analysis.min_n",
    minimum = 3L
  )
  config$analysis$strict_sample_match <- .spearman_assert_logical(
    config$analysis$strict_sample_match,
    "analysis.strict_sample_match"
  )

  config$thresholds$r_cutoff <- .spearman_assert_number(
    config$thresholds$r_cutoff,
    "thresholds.r_cutoff",
    minimum = 0,
    maximum = 1
  )
  config$thresholds$p_cutoff <- .spearman_assert_number(
    config$thresholds$p_cutoff,
    "thresholds.p_cutoff",
    minimum = 0,
    maximum = 1,
    minimum_inclusive = FALSE
  )
  config$thresholds$padj_cutoff <- .spearman_assert_number(
    config$thresholds$padj_cutoff,
    "thresholds.padj_cutoff",
    minimum = 0,
    maximum = 1,
    minimum_inclusive = FALSE
  )

  config$filters$include_positive <- .spearman_assert_logical(
    config$filters$include_positive,
    "filters.include_positive"
  )
  config$filters$include_negative <- .spearman_assert_logical(
    config$filters$include_negative,
    "filters.include_negative"
  )

  config$quadrant$enabled <- .spearman_assert_logical(
    config$quadrant$enabled,
    "quadrant.enabled"
  )
  axis_fields <- c("x_target", "x_group", "y_target", "y_group")
  for (field in axis_fields) {
    config$quadrant[field] <- list(
      .spearman_assert_optional_string(
        config$quadrant[[field]],
        paste0("quadrant.", field)
      )
    )
    if (
      config$quadrant$enabled &&
        is.null(config$quadrant[[field]])
    ) {
      .spearman_config_error(
        paste0("quadrant.", field),
        "must be set when quadrant analysis is enabled."
      )
    }
  }
  config$quadrant$significance <- .spearman_assert_enum(
    config$quadrant$significance,
    "quadrant.significance",
    c("p", "padj", "none")
  )
  config$quadrant$require_significant_both <-
    .spearman_assert_logical(
      config$quadrant$require_significant_both,
      "quadrant.require_significant_both"
    )

  output_fields <- c(
    "aligned_inputs",
    "combined_tables",
    "per_analysis_tables",
    "filtered_tables",
    "result_rds"
  )
  for (field in output_fields) {
    config$outputs[[field]] <- .spearman_assert_logical(
      config$outputs[[field]],
      paste0("outputs.", field)
    )
  }

  config$intermediate$enabled <- .spearman_assert_logical(
    config$intermediate$enabled,
    "intermediate.enabled"
  )
  config$intermediate$reuse_existing <- .spearman_assert_logical(
    config$intermediate$reuse_existing,
    "intermediate.reuse_existing"
  )

  plot_formats <- c("png", "pdf")
  config$plots$overview$enabled <- .spearman_assert_logical(
    config$plots$overview$enabled,
    "plots.overview.enabled"
  )
  config$plots$overview$label_top_n_each <-
    .spearman_assert_integer(
      config$plots$overview$label_top_n_each,
      "plots.overview.label_top_n_each",
      minimum = 0L
    )
  config$plots$overview$width <- .spearman_assert_number(
    config$plots$overview$width,
    "plots.overview.width",
    minimum = 0,
    minimum_inclusive = FALSE
  )
  config$plots$overview$height <- .spearman_assert_number(
    config$plots$overview$height,
    "plots.overview.height",
    minimum = 0,
    minimum_inclusive = FALSE
  )
  config$plots$overview$dpi <- .spearman_assert_integer(
    config$plots$overview$dpi,
    "plots.overview.dpi",
    minimum = 1L
  )
  config$plots$overview$formats <-
    .spearman_assert_character_vector(
      config$plots$overview$formats,
      "plots.overview.formats"
    )
  invalid_overview_formats <- setdiff(
    config$plots$overview$formats,
    plot_formats
  )
  if (length(invalid_overview_formats) > 0L) {
    .spearman_config_error(
      "plots.overview.formats",
      paste0(
        "unsupported format(s): ",
        paste(invalid_overview_formats, collapse = ", "),
        ". Allowed formats: ",
        paste(plot_formats, collapse = ", "),
        "."
      )
    )
  }

  config$plots$single_feature$enabled <- .spearman_assert_logical(
    config$plots$single_feature$enabled,
    "plots.single_feature.enabled"
  )
  config$plots$single_feature$width <- .spearman_assert_number(
    config$plots$single_feature$width,
    "plots.single_feature.width",
    minimum = 0,
    minimum_inclusive = FALSE
  )
  config$plots$single_feature$height <- .spearman_assert_number(
    config$plots$single_feature$height,
    "plots.single_feature.height",
    minimum = 0,
    minimum_inclusive = FALSE
  )
  config$plots$single_feature$dpi <- .spearman_assert_integer(
    config$plots$single_feature$dpi,
    "plots.single_feature.dpi",
    minimum = 1L
  )
  config$plots$single_feature$formats <-
    .spearman_assert_character_vector(
      config$plots$single_feature$formats,
      "plots.single_feature.formats"
    )
  invalid_single_formats <- setdiff(
    config$plots$single_feature$formats,
    plot_formats
  )
  if (length(invalid_single_formats) > 0L) {
    .spearman_config_error(
      "plots.single_feature.formats",
      paste0(
        "unsupported format(s): ",
        paste(invalid_single_formats, collapse = ", "),
        ". Allowed formats: ",
        paste(plot_formats, collapse = ", "),
        "."
      )
    )
  }

  single_items <- config$plots$single_feature$items
  .spearman_assert_item_sequence(
    single_items,
    "plots.single_feature.items",
    config$plots$single_feature$enabled
  )
  single_items <- .spearman_normalize_item_defaults(
    single_items,
    list(
      feature = NULL,
      target = NULL,
      group = "All",
      type = "scatter",
      add_lm = TRUE
    )
  )
  for (i in seq_along(single_items)) {
    prefix <- paste0("plots.single_feature.items[", i, "]")
    .spearman_assert_mapping(single_items[[i]], prefix)
    single_items[[i]]$feature <- .spearman_assert_string(
      single_items[[i]]$feature,
      paste0(prefix, ".feature")
    )
    single_items[[i]]$target <- .spearman_assert_string(
      single_items[[i]]$target,
      paste0(prefix, ".target")
    )
    single_items[[i]]$group <- .spearman_assert_string(
      single_items[[i]]$group,
      paste0(prefix, ".group")
    )
    single_items[[i]]$type <- .spearman_assert_enum(
      single_items[[i]]$type,
      paste0(prefix, ".type"),
      c("scatter", "ordinal")
    )
    single_items[[i]]$add_lm <- .spearman_assert_logical(
      single_items[[i]]$add_lm,
      paste0(prefix, ".add_lm")
    )
  }
  config$plots$single_feature$items <- single_items

  config$plots$quadrant$enabled <- .spearman_assert_logical(
    config$plots$quadrant$enabled,
    "plots.quadrant.enabled"
  )
  config$plots$quadrant$width <- .spearman_assert_number(
    config$plots$quadrant$width,
    "plots.quadrant.width",
    minimum = 0,
    minimum_inclusive = FALSE
  )
  config$plots$quadrant$height <- .spearman_assert_number(
    config$plots$quadrant$height,
    "plots.quadrant.height",
    minimum = 0,
    minimum_inclusive = FALSE
  )
  config$plots$quadrant$dpi <- .spearman_assert_integer(
    config$plots$quadrant$dpi,
    "plots.quadrant.dpi",
    minimum = 1L
  )
  config$plots$quadrant$formats <-
    .spearman_assert_character_vector(
      config$plots$quadrant$formats,
      "plots.quadrant.formats"
    )
  invalid_quadrant_formats <- setdiff(
    config$plots$quadrant$formats,
    plot_formats
  )
  if (length(invalid_quadrant_formats) > 0L) {
    .spearman_config_error(
      "plots.quadrant.formats",
      paste0(
        "unsupported format(s): ",
        paste(invalid_quadrant_formats, collapse = ", "),
        ". Allowed formats: ",
        paste(plot_formats, collapse = ", "),
        "."
      )
    )
  }

  config$validation$manual$enabled <- .spearman_assert_logical(
    config$validation$manual$enabled,
    "validation.manual.enabled"
  )
  manual_items <- config$validation$manual$items
  .spearman_assert_item_sequence(
    manual_items,
    "validation.manual.items",
    config$validation$manual$enabled
  )
  manual_items <- .spearman_normalize_item_defaults(
    manual_items,
    list(
      feature = NULL,
      target = NULL,
      group = "All"
    )
  )
  for (i in seq_along(manual_items)) {
    prefix <- paste0("validation.manual.items[", i, "]")
    .spearman_assert_mapping(manual_items[[i]], prefix)
    manual_items[[i]]$feature <- .spearman_assert_string(
      manual_items[[i]]$feature,
      paste0(prefix, ".feature")
    )
    manual_items[[i]]$target <- .spearman_assert_string(
      manual_items[[i]]$target,
      paste0(prefix, ".target")
    )
    manual_items[[i]]$group <- .spearman_assert_string(
      manual_items[[i]]$group,
      paste0(prefix, ".group")
    )
  }
  config$validation$manual$items <- manual_items

  config$logging$verbose <- .spearman_assert_logical(
    config$logging$verbose,
    "logging.verbose"
  )
  config
}

.spearman_is_absolute_path <- function(path) {
  startsWith(path, "/") ||
    startsWith(path, "\\\\") ||
    grepl("^[A-Za-z]:[/\\\\]", path) ||
    startsWith(path, "~")
}

.spearman_resolve_config_path <- function(path, config_dir) {
  path <- path.expand(path)
  candidate <- if (.spearman_is_absolute_path(path)) {
    path
  } else {
    file.path(config_dir, path)
  }
  normalizePath(
    candidate,
    winslash = "/",
    mustWork = FALSE
  )
}

#' Read and validate a Spearman workflow YAML configuration.
#'
#' Every relative data or output path is resolved against the directory that
#' contains `path`, never against the process working directory.
#'
#' @param path Path to a YAML configuration file.
#' @return A normalized configuration list with absolute `config_path`,
#'   `config_dir`, input paths, and output path.
read_spearman_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop(
      paste0(
        "Package `yaml` is required to read the configuration. ",
        "Activate the `spearman-r` environment or install `r-yaml` there."
      ),
      call. = FALSE
    )
  }
  path <- .spearman_assert_string(path, "config path")
  if (!file.exists(path)) {
    stop(
      "Configuration file does not exist: ",
      path,
      call. = FALSE
    )
  }
  if (dir.exists(path)) {
    stop(
      "Configuration path is a directory, not a file: ",
      path,
      call. = FALSE
    )
  }

  config_path <- normalizePath(
    path,
    winslash = "/",
    mustWork = TRUE
  )
  config_dir <- dirname(config_path)
  supplied <- tryCatch(
    yaml::read_yaml(config_path),
    error = function(error) {
      stop(
        "Unable to parse YAML configuration `",
        config_path,
        "`: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
  if (
    !is.list(supplied) ||
      is.null(names(supplied)) ||
      any(!nzchar(names(supplied)))
  ) {
    stop(
      "Configuration root must be a non-empty YAML mapping: ",
      config_path,
      call. = FALSE
    )
  }

  .spearman_check_unknown_keys(
    supplied,
    .spearman_config_schema(),
    "config"
  )
  config <- .spearman_merge_config(
    .spearman_config_defaults(),
    supplied
  )
  config <- .spearman_validate_config(config)

  config$inputs$expression_matrix <-
    .spearman_resolve_config_path(
      config$inputs$expression_matrix,
      config_dir
    )
  config$inputs$score_metadata <-
    .spearman_resolve_config_path(
      config$inputs$score_metadata,
      config_dir
    )
  config$output_dir <- .spearman_resolve_config_path(
    config$output_dir,
    config_dir
  )

  input_fields <- c("expression_matrix", "score_metadata")
  for (field in input_fields) {
    input_path <- config$inputs[[field]]
    if (!file.exists(input_path)) {
      .spearman_config_error(
        paste0("inputs.", field),
        paste0(
          "resolved to a file that does not exist: ",
          input_path,
          ". Relative paths are resolved from ",
          config_dir,
          "."
        )
      )
    }
    if (dir.exists(input_path)) {
      .spearman_config_error(
        paste0("inputs.", field),
        paste0("must refer to a file, not a directory: ", input_path)
      )
    }
  }
  if (
    file.exists(config$output_dir) &&
      !dir.exists(config$output_dir)
  ) {
    .spearman_config_error(
      "output_dir",
      paste0(
        "resolves to an existing file, not a directory: ",
        config$output_dir
      )
    )
  }

  config$config_path <- config_path
  config$config_dir <- config_dir
  config
}
