# ============================================================
# Reusable Spearman correlation workflow functions
# ============================================================
#
# This file defines reusable functions only.
# Running source("R/run_spearman_workflow.R") loads the functions
# but does not immediately start an analysis.
#
# Main function:
#   run_spearman_workflow()
#
# Plotting functions:
#   plot_spearman_overview()
#   plot_spearman_feature()
#
# Core statistical analysis uses base R.
# Plotting requires ggplot2; labels additionally require ggrepel.
# ============================================================


# ------------------------------------------------------------
# Internal helper: safe filename
# ------------------------------------------------------------

.spearman_safe_filename <- function(x) {
  gsub(
    pattern = "[^A-Za-z0-9._-]+",
    replacement = "_",
    x = as.character(x)
  )
}


# ------------------------------------------------------------
# Internal helper: write CSV without row names
# ------------------------------------------------------------

.spearman_write_csv <- function(x, path) {
  utils::write.csv(
    x = x,
    file = path,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}


# ------------------------------------------------------------
# Internal helper: check optional plotting packages
# ------------------------------------------------------------

.spearman_check_plot_packages <- function(label_top_n_each = 0L) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "绘图需要安装 ggplot2：install.packages('ggplot2')",
      call. = FALSE
    )
  }

  if (
    label_top_n_each > 0L &&
    !requireNamespace("ggrepel", quietly = TRUE)
  ) {
    stop(
      "标注 feature 需要安装 ggrepel：",
      "install.packages('ggrepel')",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------
# Internal helper: one target × one analysis group
# ------------------------------------------------------------

.spearman_run_one_target_group <- function(
    feature_mat,
    target_sub,
    target_col,
    sample_col,
    group_label,
    min_n
) {

  samples_use <- as.character(
    target_sub[[sample_col]]
  )

  missing_matrix_samples <- setdiff(
    samples_use,
    colnames(feature_mat)
  )

  if (length(missing_matrix_samples) > 0L) {
    stop(
      "以下样本不存在于 feature_mat：",
      paste(missing_matrix_samples, collapse = ", "),
      call. = FALSE
    )
  }

  feature_sub <- feature_mat[
    ,
    samples_use,
    drop = FALSE
  ]

  if (!identical(colnames(feature_sub), samples_use)) {
    stop(
      "feature_sub 与 target_sub 的样本顺序不一致。",
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

    x <- as.numeric(
      feature_sub[i, ]
    )

    complete_idx <- (
      stats::complete.cases(x, y) &
        is.finite(x) &
        is.finite(y)
    )

    x_complete <- x[complete_idx]
    y_complete <- y[complete_idx]

    n_complete[i] <- sum(complete_idx)

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
          method = "spearman",
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
      ifelse(
        rho > 0,
        "Positive",
        ifelse(rho < 0, "Negative", "Zero")
      )
    ),
    status = status,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  valid_idx <- which(
    result_df$status == "OK" &
      !is.na(result_df$pvalue)
  )

  if (length(valid_idx) > 0L) {
    result_df$padj[valid_idx] <- stats::p.adjust(
      result_df$pvalue[valid_idx],
      method = "BH"
    )
  }

  result_df <- result_df[
    order(
      is.na(result_df$pvalue),
      result_df$pvalue
    ),
    ,
    drop = FALSE
  ]

  rownames(result_df) <- NULL

  result_df
}


#' Plot a Spearman correlation overview
#'
#' Each point represents one feature.
#'
#' @param result_df Complete result table returned in
#'   `workflow_result$results`.
#' @param target_name One target name.
#' @param group_name One group label, such as `"All"` or `"IA"`.
#' @param r_cutoff Absolute rho threshold.
#' @param p_cutoff Raw p-value threshold.
#' @param label_top_n_each Number of positive and negative features
#'   to label separately.
#'
#' @return A ggplot object.
plot_spearman_overview <- function(
    result_df,
    target_name,
    group_name,
    r_cutoff = 0.5,
    p_cutoff = 0.05,
    label_top_n_each = 0L
) {

  label_top_n_each <- as.integer(label_top_n_each)

  .spearman_check_plot_packages(
    label_top_n_each = label_top_n_each
  )

  required_cols <- c(
    "feature",
    "target",
    "group",
    "rho",
    "pvalue",
    "status"
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(result_df)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "result_df 缺少以下列：",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  plot_df <- result_df[
    result_df$target == target_name &
      result_df$group == group_name &
      result_df$status == "OK" &
      !is.na(result_df$rho) &
      !is.na(result_df$pvalue),
    ,
    drop = FALSE
  ]

  if (nrow(plot_df) == 0L) {
    stop(
      "未找到指定 target × group 的有效结果。",
      call. = FALSE
    )
  }

  plot_df$neg_log10_p <- -log10(
    pmax(
      plot_df$pvalue,
      .Machine$double.xmin
    )
  )

  plot_df$correlation_class <- "Not significant"

  plot_df$correlation_class[
    plot_df$rho >= r_cutoff &
      plot_df$pvalue < p_cutoff
  ] <- "Positive"

  plot_df$correlation_class[
    plot_df$rho <= -r_cutoff &
      plot_df$pvalue < p_cutoff
  ] <- "Negative"

  plot_df$correlation_class <- factor(
    plot_df$correlation_class,
    levels = c(
      "Negative",
      "Not significant",
      "Positive"
    )
  )

  title_text <- if (identical(group_name, "All")) {
    "All-sample correlation overview"
  } else {
    paste0(group_name, " within-group correlation overview")
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = rho,
      y = neg_log10_p,
      color = correlation_class
    )
  ) +
    ggplot2::geom_point(
      size = 1.7,
      alpha = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = c(-r_cutoff, r_cutoff),
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Negative" = "#3B6FB6",
        "Not significant" = "grey75",
        "Positive" = "#D94B3D"
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = title_text,
      subtitle = paste0("Target: ", target_name),
      x = "Spearman correlation coefficient",
      y = "-log10(p value)",
      color = "Correlation"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5
      ),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (label_top_n_each > 0L) {

    positive_df <- plot_df[
      plot_df$correlation_class == "Positive",
      ,
      drop = FALSE
    ]

    negative_df <- plot_df[
      plot_df$correlation_class == "Negative",
      ,
      drop = FALSE
    ]

    if (nrow(positive_df) > 0L) {
      positive_df <- positive_df[
        order(
          positive_df$pvalue,
          -positive_df$rho
        ),
        ,
        drop = FALSE
      ]

      positive_df <- utils::head(
        positive_df,
        label_top_n_each
      )
    }

    if (nrow(negative_df) > 0L) {
      negative_df <- negative_df[
        order(
          negative_df$pvalue,
          negative_df$rho
        ),
        ,
        drop = FALSE
      ]

      negative_df <- utils::head(
        negative_df,
        label_top_n_each
      )
    }

    label_df <- rbind(
      positive_df,
      negative_df
    )

    if (nrow(label_df) > 0L) {
      p <- p +
        ggrepel::geom_text_repel(
          data = label_df,
          ggplot2::aes(label = feature),
          size = 3,
          max.overlaps = Inf,
          box.padding = 0.3,
          point.padding = 0.2,
          show.legend = FALSE
        )
    }
  }

  p
}


#' Run a reusable Spearman correlation workflow
#'
#' @param feature_df A data frame whose rows are features and whose
#'   sample values are stored in columns.
#' @param target_df A sample-level data frame containing sample names,
#'   target variables and optionally group information.
#' @param feature_col Column in `feature_df` containing feature names.
#' @param sample_col Column in `target_df` containing sample names.
#' @param selected_features `NULL` for all features, or a character
#'   vector of selected feature names.
#' @param target_cols `NULL` for all eligible numeric target columns,
#'   or a character vector of selected target columns.
#' @param run_all_samples Whether to run one analysis using all matched
#'   samples.
#' @param run_group_analysis Whether to run separate within-group
#'   analyses.
#' @param group_col Group column in `target_df`; required when
#'   `run_group_analysis = TRUE`.
#' @param selected_groups `NULL` for all non-missing groups, or a
#'   character vector of selected groups.
#' @param min_n Minimum complete sample size for each correlation.
#' @param r_cutoff Absolute rho threshold used for result screening.
#' @param p_cutoff Raw p-value threshold used for result screening.
#' @param padj_cutoff BH-adjusted p-value threshold.
#' @param strict_sample_match Stop when the two inputs contain any
#'   unmatched sample.
#' @param outdir Output directory.
#' @param save_results Whether to save aligned inputs, tables, plots
#'   and an RDS result object.
#' @param make_overview_plot Whether to create overview plots.
#' @param label_top_n_each Number of positive and negative features
#'   labelled on each overview plot.
#' @param verbose Whether to print progress messages.
#'
#' @return An object of class `spearman_workflow_result`, implemented
#'   as a list containing aligned data, result tables and plots.
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

  # ----------------------------------------------------------
  # 1. Check scalar parameters
  # ----------------------------------------------------------

  feature_col <- as.character(feature_col)
  sample_col <- as.character(sample_col)

  if (
    length(feature_col) != 1L ||
    is.na(feature_col) ||
    feature_col == ""
  ) {
    stop(
      "feature_col 必须是一个非空列名。",
      call. = FALSE
    )
  }

  if (
    length(sample_col) != 1L ||
    is.na(sample_col) ||
    sample_col == ""
  ) {
    stop(
      "sample_col 必须是一个非空列名。",
      call. = FALSE
    )
  }

  logical_parameters <- list(
    run_all_samples = run_all_samples,
    run_group_analysis = run_group_analysis,
    strict_sample_match = strict_sample_match,
    save_results = save_results,
    make_overview_plot = make_overview_plot,
    verbose = verbose
  )

  invalid_logical <- names(logical_parameters)[
    !vapply(
      logical_parameters,
      function(x) {
        is.logical(x) &&
          length(x) == 1L &&
          !is.na(x)
      },
      logical(1)
    )
  ]

  if (length(invalid_logical) > 0L) {
    stop(
      "以下参数必须是单个 TRUE/FALSE：",
      paste(invalid_logical, collapse = ", "),
      call. = FALSE
    )
  }

  if (!run_all_samples && !run_group_analysis) {
    stop(
      "run_all_samples 和 run_group_analysis ",
      "不能同时为 FALSE。",
      call. = FALSE
    )
  }

  min_n <- as.integer(min_n)
  label_top_n_each <- as.integer(label_top_n_each)

  if (
    length(min_n) != 1L ||
    is.na(min_n) ||
    min_n < 3L
  ) {
    stop(
      "min_n 必须是大于或等于 3 的整数。",
      call. = FALSE
    )
  }

  if (
    length(label_top_n_each) != 1L ||
    is.na(label_top_n_each) ||
    label_top_n_each < 0L
  ) {
    stop(
      "label_top_n_each 必须是非负整数。",
      call. = FALSE
    )
  }

  if (
    length(r_cutoff) != 1L ||
    is.na(r_cutoff) ||
    r_cutoff < 0 ||
    r_cutoff > 1
  ) {
    stop(
      "r_cutoff 必须位于 0 到 1 之间。",
      call. = FALSE
    )
  }

  for (
    parameter_name in c("p_cutoff", "padj_cutoff")
  ) {
    parameter_value <- get(parameter_name)

    if (
      length(parameter_value) != 1L ||
      is.na(parameter_value) ||
      parameter_value <= 0 ||
      parameter_value > 1
    ) {
      stop(
        parameter_name,
        " 必须位于 0 到 1 之间，且不能等于 0。",
        call. = FALSE
      )
    }
  }

  if (save_results) {
    if (
      length(outdir) != 1L ||
      is.na(outdir) ||
      outdir == ""
    ) {
      stop(
        "save_results = TRUE 时必须提供有效 outdir。",
        call. = FALSE
      )
    }
  }

  if (make_overview_plot) {
    .spearman_check_plot_packages(
      label_top_n_each = label_top_n_each
    )
  }

  # ----------------------------------------------------------
  # 2. Copy inputs and standardize names
  # ----------------------------------------------------------

  if (!is.data.frame(feature_df)) {
    stop(
      "feature_df 必须是 data.frame。",
      call. = FALSE
    )
  }

  if (!is.data.frame(target_df)) {
    stop(
      "target_df 必须是 data.frame。",
      call. = FALSE
    )
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

  colnames(feature_df) <- trimws(
    colnames(feature_df)
  )

  colnames(target_df) <- trimws(
    colnames(target_df)
  )

  if (!feature_col %in% colnames(feature_df)) {
    stop(
      "feature_df 中不存在 feature_col：",
      feature_col,
      call. = FALSE
    )
  }

  if (!sample_col %in% colnames(target_df)) {
    stop(
      "target_df 中不存在 sample_col：",
      sample_col,
      call. = FALSE
    )
  }

  if (run_group_analysis) {

    if (
      is.null(group_col) ||
      length(group_col) != 1L ||
      is.na(group_col) ||
      group_col == ""
    ) {
      stop(
        "组内分析时必须提供有效 group_col。",
        call. = FALSE
      )
    }

    group_col <- as.character(group_col)

    if (!group_col %in% colnames(target_df)) {
      stop(
        "target_df 中不存在 group_col：",
        group_col,
        call. = FALSE
      )
    }
  } else if (!is.null(group_col)) {

    group_col <- as.character(group_col)

    if (length(group_col) != 1L) {
      stop(
        "group_col 必须为 NULL 或单个列名。",
        call. = FALSE
      )
    }
  }

  feature_df[[feature_col]] <- trimws(
    as.character(feature_df[[feature_col]])
  )

  target_df[[sample_col]] <- trimws(
    as.character(target_df[[sample_col]])
  )

  if (
    !is.null(group_col) &&
    group_col %in% colnames(target_df)
  ) {
    target_df[[group_col]] <- trimws(
      as.character(target_df[[group_col]])
    )
  }

  # ----------------------------------------------------------
  # 3. Check names and duplicated identifiers
  # ----------------------------------------------------------

  if (
    any(is.na(feature_df[[feature_col]])) ||
    any(feature_df[[feature_col]] == "")
  ) {
    stop(
      "feature_df 的 feature 名称中存在 NA 或空字符串。",
      call. = FALSE
    )
  }

  if (
    any(is.na(target_df[[sample_col]])) ||
    any(target_df[[sample_col]] == "")
  ) {
    stop(
      "target_df 的 sample 名称中存在 NA 或空字符串。",
      call. = FALSE
    )
  }

  duplicated_features <- unique(
    feature_df[[feature_col]][
      duplicated(feature_df[[feature_col]])
    ]
  )

  if (length(duplicated_features) > 0L) {
    stop(
      "feature_df 中存在重复 feature：",
      paste(
        utils::head(duplicated_features, 20L),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  duplicated_target_samples <- unique(
    target_df[[sample_col]][
      duplicated(target_df[[sample_col]])
    ]
  )

  if (length(duplicated_target_samples) > 0L) {
    stop(
      "target_df 中存在重复 sample：",
      paste(
        utils::head(duplicated_target_samples, 20L),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  feature_samples <- setdiff(
    colnames(feature_df),
    feature_col
  )

  if (length(feature_samples) == 0L) {
    stop(
      "feature_df 中没有检测到样本列。",
      call. = FALSE
    )
  }

  duplicated_feature_samples <- unique(
    feature_samples[
      duplicated(feature_samples)
    ]
  )

  if (length(duplicated_feature_samples) > 0L) {
    stop(
      "feature_df 中存在重复样本列名：",
      paste(
        duplicated_feature_samples,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # 4. Resolve targets
  # ----------------------------------------------------------

  excluded_target_cols <- sample_col

  if (
    !is.null(group_col) &&
    group_col %in% colnames(target_df)
  ) {
    excluded_target_cols <- c(
      excluded_target_cols,
      group_col
    )
  }

  if (is.null(target_cols)) {

    candidate_target_cols <- setdiff(
      colnames(target_df),
      excluded_target_cols
    )

    candidate_numeric <- vapply(
      target_df[
        ,
        candidate_target_cols,
        drop = FALSE
      ],
      is.numeric,
      logical(1)
    )

    target_cols <- candidate_target_cols[
      candidate_numeric
    ]

    if (length(target_cols) == 0L) {
      stop(
        "target_cols = NULL，但未检测到可用的 numeric target。",
        call. = FALSE
      )
    }

    if (verbose) {
      message(
        "将分析全部 numeric target：",
        paste(target_cols, collapse = ", ")
      )
    }

  } else {

    target_cols <- trimws(
      as.character(target_cols)
    )

    if (
      length(target_cols) == 0L ||
      any(is.na(target_cols)) ||
      any(target_cols == "")
    ) {
      stop(
        "target_cols 不能为空。",
        call. = FALSE
      )
    }

    if (anyDuplicated(target_cols) > 0L) {
      stop(
        "target_cols 中存在重复 target。",
        call. = FALSE
      )
    }

    missing_target_cols <- setdiff(
      target_cols,
      colnames(target_df)
    )

    if (length(missing_target_cols) > 0L) {
      stop(
        "target_df 中缺少以下 target：",
        paste(missing_target_cols, collapse = ", "),
        call. = FALSE
      )
    }

    invalid_target_cols <- intersect(
      target_cols,
      excluded_target_cols
    )

    if (length(invalid_target_cols) > 0L) {
      stop(
        "以下列不能作为 target：",
        paste(invalid_target_cols, collapse = ", "),
        call. = FALSE
      )
    }
  }

  # ----------------------------------------------------------
  # 5. Check numeric data
  # ----------------------------------------------------------

  feature_numeric_check <- vapply(
    feature_df[
      ,
      feature_samples,
      drop = FALSE
    ],
    is.numeric,
    logical(1)
  )

  if (!all(feature_numeric_check)) {

    non_numeric_feature_cols <- names(
      feature_numeric_check[
        !feature_numeric_check
      ]
    )

    stop(
      "feature_df 中以下样本列不是 numeric：",
      paste(
        non_numeric_feature_cols,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  target_numeric_check <- vapply(
    target_df[
      ,
      target_cols,
      drop = FALSE
    ],
    is.numeric,
    logical(1)
  )

  if (!all(target_numeric_check)) {

    non_numeric_target_cols <- names(
      target_numeric_check[
        !target_numeric_check
      ]
    )

    stop(
      "target_df 中以下 target 不是 numeric：",
      paste(
        non_numeric_target_cols,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # 6. Resolve selected features
  # ----------------------------------------------------------

  if (is.null(selected_features)) {

    feature_df_selected <- feature_df

    if (verbose) {
      message(
        "将分析全部 feature；数量：",
        nrow(feature_df_selected)
      )
    }

  } else {

    selected_features <- trimws(
      as.character(selected_features)
    )

    if (
      length(selected_features) == 0L ||
      any(is.na(selected_features)) ||
      any(selected_features == "")
    ) {
      stop(
        "selected_features 不能为空。",
        call. = FALSE
      )
    }

    if (anyDuplicated(selected_features) > 0L) {
      stop(
        "selected_features 中存在重复 feature。",
        call. = FALSE
      )
    }

    missing_features <- setdiff(
      selected_features,
      feature_df[[feature_col]]
    )

    if (length(missing_features) > 0L) {
      stop(
        "以下 selected_features 不存在：",
        paste(missing_features, collapse = ", "),
        call. = FALSE
      )
    }

    feature_df_selected <- feature_df[
      match(
        selected_features,
        feature_df[[feature_col]]
      ),
      ,
      drop = FALSE
    ]

    if (verbose) {
      message(
        "仅分析指定 feature；数量：",
        nrow(feature_df_selected)
      )
    }
  }

  # ----------------------------------------------------------
  # 7. Resolve groups
  # ----------------------------------------------------------

  groups_used <- character(0)

  if (run_group_analysis) {

    group_values <- target_df[[group_col]]

    available_groups <- unique(
      group_values[
        !is.na(group_values) &
          group_values != ""
      ]
    )

    if (length(available_groups) == 0L) {
      stop(
        "group_col 中没有可用分组。",
        call. = FALSE
      )
    }

    if (is.null(selected_groups)) {

      selected_groups <- available_groups

      if (verbose) {
        message(
          "将分析全部组：",
          paste(selected_groups, collapse = ", ")
        )
      }

    } else {

      selected_groups <- trimws(
        as.character(selected_groups)
      )

      if (
        length(selected_groups) == 0L ||
        any(is.na(selected_groups)) ||
        any(selected_groups == "")
      ) {
        stop(
          "selected_groups 不能为空。",
          call. = FALSE
        )
      }

      if (anyDuplicated(selected_groups) > 0L) {
        stop(
          "selected_groups 中存在重复组。",
          call. = FALSE
        )
      }

      missing_groups <- setdiff(
        selected_groups,
        available_groups
      )

      if (length(missing_groups) > 0L) {
        stop(
          "以下 selected_groups 不存在：",
          paste(missing_groups, collapse = ", "),
          call. = FALSE
        )
      }
    }

    groups_used <- selected_groups

    if (verbose) {
      message("分组样本量：")
      print(
        table(
          target_df[[group_col]],
          useNA = "ifany"
        )
      )
    }
  }

  # ----------------------------------------------------------
  # 8. Match and align samples
  # ----------------------------------------------------------

  target_samples <- target_df[[sample_col]]

  only_in_target <- setdiff(
    target_samples,
    feature_samples
  )

  only_in_feature <- setdiff(
    feature_samples,
    target_samples
  )

  if (verbose) {
    message(
      "target_df 中有、feature_df 中没有的样本数：",
      length(only_in_target)
    )

    message(
      "feature_df 中有、target_df 中没有的样本数：",
      length(only_in_feature)
    )
  }

  if (
    strict_sample_match &&
    (
      length(only_in_target) > 0L ||
        length(only_in_feature) > 0L
    )
  ) {
    stop(
      "检测到不匹配样本；",
      "请检查样本名或关闭 strict_sample_match。",
      call. = FALSE
    )
  }

  common_samples <- intersect(
    feature_samples,
    target_samples
  )

  if (verbose) {
    message(
      "共同样本数：",
      length(common_samples)
    )
  }

  if (length(common_samples) < min_n) {
    stop(
      "共同样本数小于 min_n，无法继续分析。",
      call. = FALSE
    )
  }

  feature_df_aligned <- feature_df_selected[
    ,
    c(feature_col, common_samples),
    drop = FALSE
  ]

  target_row_index <- match(
    common_samples,
    target_df[[sample_col]]
  )

  if (anyNA(target_row_index)) {
    stop(
      "match() 产生 NA，样本匹配失败。",
      call. = FALSE
    )
  }

  target_df_aligned <- target_df[
    target_row_index,
    ,
    drop = FALSE
  ]

  sample_order_ok <- identical(
    common_samples,
    target_df_aligned[[sample_col]]
  )

  if (!sample_order_ok) {
    stop(
      "样本顺序统一失败。",
      call. = FALSE
    )
  }

  feature_mat <- as.matrix(
    feature_df_aligned[
      ,
      common_samples,
      drop = FALSE
    ]
  )

  storage.mode(feature_mat) <- "double"

  rownames(feature_mat) <- feature_df_aligned[[feature_col]]

  # ----------------------------------------------------------
  # 9. Run all-sample and group-wise analyses
  # ----------------------------------------------------------

  if (verbose) {
    message(
      "实际分析 feature 数量：",
      nrow(feature_mat)
    )

    message(
      "实际分析 target 数量：",
      length(target_cols)
    )
  }

  result_list <- list()
  result_index <- 1L

  if (run_all_samples) {

    if (verbose) {
      message("开始运行全样本相关性分析。")
    }

    for (target_name in target_cols) {

      result_list[[result_index]] <-
        .spearman_run_one_target_group(
          feature_mat = feature_mat,
          target_sub = target_df_aligned,
          target_col = target_name,
          sample_col = sample_col,
          group_label = "All",
          min_n = min_n
        )

      result_index <- result_index + 1L
    }
  }

  if (run_group_analysis) {

    if (verbose) {
      message("开始运行组内相关性分析。")
    }

    for (group_name in selected_groups) {

      group_vector <- target_df_aligned[[group_col]]

      group_index <- (
        !is.na(group_vector) &
          group_vector == group_name
      )

      target_group <- target_df_aligned[
        group_index,
        ,
        drop = FALSE
      ]

      if (verbose) {
        message(
          "当前组：",
          group_name,
          "；样本数：",
          nrow(target_group)
        )
      }

      if (nrow(target_group) < min_n) {
        warning(
          "组 ",
          group_name,
          " 的样本数小于 min_n，已跳过。",
          call. = FALSE
        )
        next
      }

      for (target_name in target_cols) {

        result_list[[result_index]] <-
          .spearman_run_one_target_group(
            feature_mat = feature_mat,
            target_sub = target_group,
            target_col = target_name,
            sample_col = sample_col,
            group_label = group_name,
            min_n = min_n
          )

        result_index <- result_index + 1L
      }
    }
  }

  if (length(result_list) == 0L) {
    stop(
      "没有产生任何分析结果。",
      call. = FALSE
    )
  }

  cor_results <- do.call(
    rbind,
    result_list
  )

  rownames(cor_results) <- NULL

  # ----------------------------------------------------------
  # 10. Add screening flags
  # ----------------------------------------------------------

  cor_results$significant_by_p <- (
    cor_results$status == "OK" &
      !is.na(cor_results$rho) &
      !is.na(cor_results$pvalue) &
      cor_results$abs_rho >= r_cutoff &
      cor_results$pvalue < p_cutoff
  )

  cor_results$significant_by_padj <- (
    cor_results$status == "OK" &
      !is.na(cor_results$rho) &
      !is.na(cor_results$padj) &
      cor_results$abs_rho >= r_cutoff &
      cor_results$padj < padj_cutoff
  )

  cor_sig_p <- cor_results[
    cor_results$significant_by_p,
    ,
    drop = FALSE
  ]

  cor_sig_p <- cor_sig_p[
    order(
      cor_sig_p$group,
      cor_sig_p$target,
      cor_sig_p$pvalue
    ),
    ,
    drop = FALSE
  ]

  cor_sig_padj <- cor_results[
    cor_results$significant_by_padj,
    ,
    drop = FALSE
  ]

  cor_sig_padj <- cor_sig_padj[
    order(
      cor_sig_padj$group,
      cor_sig_padj$target,
      cor_sig_padj$padj
    ),
    ,
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # 11. Build summary table
  # ----------------------------------------------------------

  result_combinations <- unique(
    cor_results[
      ,
      c("group", "target"),
      drop = FALSE
    ]
  )

  summary_list <- vector(
    "list",
    nrow(result_combinations)
  )

  for (i in seq_len(nrow(result_combinations))) {

    group_name <- result_combinations$group[i]
    target_name <- result_combinations$target[i]

    result_sub <- cor_results[
      cor_results$group == group_name &
        cor_results$target == target_name,
      ,
      drop = FALSE
    ]

    summary_list[[i]] <- data.frame(
      group = group_name,
      target = target_name,
      total_features = nrow(result_sub),
      valid_results = sum(
        result_sub$status == "OK",
        na.rm = TRUE
      ),
      significant_by_p = sum(
        result_sub$significant_by_p,
        na.rm = TRUE
      ),
      significant_by_padj = sum(
        result_sub$significant_by_padj,
        na.rm = TRUE
      ),
      positive_by_p = sum(
        result_sub$significant_by_p &
          result_sub$direction == "Positive",
        na.rm = TRUE
      ),
      negative_by_p = sum(
        result_sub$significant_by_p &
          result_sub$direction == "Negative",
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }

  result_summary <- do.call(
    rbind,
    summary_list
  )

  rownames(result_summary) <- NULL

  # ----------------------------------------------------------
  # 12. Create overview plots
  # ----------------------------------------------------------

  overview_plots <- list()

  if (make_overview_plot) {

    for (i in seq_len(nrow(result_combinations))) {

      group_name <- result_combinations$group[i]
      target_name <- result_combinations$target[i]

      plot_name <- paste0(
        .spearman_safe_filename(group_name),
        "__",
        .spearman_safe_filename(target_name)
      )

      overview_plots[[plot_name]] <-
        plot_spearman_overview(
          result_df = cor_results,
          target_name = target_name,
          group_name = group_name,
          r_cutoff = r_cutoff,
          p_cutoff = p_cutoff,
          label_top_n_each = label_top_n_each
        )
    }
  }

  # ----------------------------------------------------------
  # 13. Assemble returned object
  # ----------------------------------------------------------

  workflow_result <- list(
    call = match.call(),
    parameters = list(
      feature_col = feature_col,
      sample_col = sample_col,
      selected_features = selected_features,
      target_cols = target_cols,
      run_all_samples = run_all_samples,
      run_group_analysis = run_group_analysis,
      group_col = group_col,
      selected_groups = groups_used,
      min_n = min_n,
      r_cutoff = r_cutoff,
      p_cutoff = p_cutoff,
      padj_cutoff = padj_cutoff,
      strict_sample_match = strict_sample_match
    ),
    sample_report = list(
      only_in_target = only_in_target,
      only_in_feature = only_in_feature,
      common_samples = common_samples,
      sample_order_ok = sample_order_ok
    ),
    features_used = feature_df_aligned[[feature_col]],
    targets_used = target_cols,
    groups_used = groups_used,
    feature_df_aligned = feature_df_aligned,
    target_df_aligned = target_df_aligned,
    feature_mat = feature_mat,
    results = cor_results,
    significant_p = cor_sig_p,
    significant_padj = cor_sig_padj,
    summary = result_summary,
    overview_plots = overview_plots
  )

  class(workflow_result) <- c(
    "spearman_workflow_result",
    "list"
  )

  # ----------------------------------------------------------
  # 14. Save files
  # ----------------------------------------------------------

  if (save_results) {

    table_dir <- file.path(
      outdir,
      "tables"
    )

    plot_dir <- file.path(
      outdir,
      "plots"
    )

    input_dir <- file.path(
      outdir,
      "aligned_inputs"
    )

    dir.create(
      table_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    dir.create(
      plot_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    dir.create(
      input_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    .spearman_write_csv(
      feature_df_aligned,
      file.path(
        input_dir,
        "feature_df_aligned.csv"
      )
    )

    .spearman_write_csv(
      target_df_aligned,
      file.path(
        input_dir,
        "target_df_aligned.csv"
      )
    )

    .spearman_write_csv(
      data.frame(
        feature = feature_df_aligned[[feature_col]],
        stringsAsFactors = FALSE
      ),
      file.path(
        input_dir,
        "features_used.csv"
      )
    )

    .spearman_write_csv(
      data.frame(
        target = target_cols,
        stringsAsFactors = FALSE
      ),
      file.path(
        input_dir,
        "targets_used.csv"
      )
    )

    if (run_group_analysis) {
      .spearman_write_csv(
        data.frame(
          group = groups_used,
          stringsAsFactors = FALSE
        ),
        file.path(
          input_dir,
          "groups_used.csv"
        )
      )
    }

    .spearman_write_csv(
      cor_results,
      file.path(
        table_dir,
        "Spearman_all_results.csv"
      )
    )

    .spearman_write_csv(
      cor_sig_p,
      file.path(
        table_dir,
        paste0(
          "Spearman_sig_absR",
          r_cutoff,
          "_p",
          p_cutoff,
          ".csv"
        )
      )
    )

    .spearman_write_csv(
      cor_sig_padj,
      file.path(
        table_dir,
        paste0(
          "Spearman_sig_absR",
          r_cutoff,
          "_padj",
          padj_cutoff,
          ".csv"
        )
      )
    )

    .spearman_write_csv(
      result_summary,
      file.path(
        table_dir,
        "Spearman_result_summary.csv"
      )
    )

    for (i in seq_len(nrow(result_combinations))) {

      group_name <- result_combinations$group[i]
      target_name <- result_combinations$target[i]

      result_sub <- cor_results[
        cor_results$group == group_name &
          cor_results$target == target_name,
        ,
        drop = FALSE
      ]

      result_filename <- paste0(
        .spearman_safe_filename(group_name),
        "__",
        .spearman_safe_filename(target_name),
        "__all_results.csv"
      )

      .spearman_write_csv(
        result_sub,
        file.path(
          table_dir,
          result_filename
        )
      )
    }

    if (make_overview_plot) {

      for (plot_name in names(overview_plots)) {

        ggplot2::ggsave(
          filename = file.path(
            plot_dir,
            paste0(
              plot_name,
              "__correlation_overview.png"
            )
          ),
          plot = overview_plots[[plot_name]],
          width = 7,
          height = 5.5,
          dpi = 300
        )

        ggplot2::ggsave(
          filename = file.path(
            plot_dir,
            paste0(
              plot_name,
              "__correlation_overview.pdf"
            )
          ),
          plot = overview_plots[[plot_name]],
          width = 7,
          height = 5.5
        )
      }
    }

    saveRDS(
      workflow_result,
      file = file.path(
        outdir,
        "Spearman_workflow_results.rds"
      )
    )
  }

  if (verbose) {
    print(result_summary)
    message("Spearman workflow completed.")
  }

  workflow_result
}


#' Plot one feature-target Spearman relationship
#'
#' @param workflow_result Result returned by
#'   `run_spearman_workflow()`.
#' @param feature_name One feature name.
#' @param target_name One target name.
#' @param group_name `"All"` or one group value.
#' @param plot_type Plot type:
#'   `"scatter"` for a continuous target;
#'   `"ordinal"` for an ordered target such as 0, 1, 2, 3, 4.
#' @param add_lm Whether to add a linear trend line when
#'   `plot_type = "scatter"`.
#'
#' @return A ggplot object.
plot_spearman_feature <- function(
    workflow_result,
    feature_name,
    target_name,
    group_name = "All",
    plot_type = c("scatter", "ordinal"),
    add_lm = TRUE
) {

  plot_type <- match.arg(plot_type)

  .spearman_check_plot_packages(
    label_top_n_each = 0L
  )

  if (
    !inherits(
      workflow_result,
      "spearman_workflow_result"
    )
  ) {
    stop(
      "workflow_result 必须来自 run_spearman_workflow()。",
      call. = FALSE
    )
  }

  feature_mat <- workflow_result$feature_mat
  target_df_aligned <- workflow_result$target_df_aligned

  sample_col <- workflow_result$parameters$sample_col
  group_col <- workflow_result$parameters$group_col
  min_n <- workflow_result$parameters$min_n

  if (!feature_name %in% rownames(feature_mat)) {
    stop(
      "feature_mat 中不存在 feature：",
      feature_name,
      call. = FALSE
    )
  }

  if (!target_name %in% colnames(target_df_aligned)) {
    stop(
      "target_df_aligned 中不存在 target：",
      target_name,
      call. = FALSE
    )
  }

  if (identical(group_name, "All")) {

    target_sub <- target_df_aligned

  } else {

    if (
      is.null(group_col) ||
      !group_col %in% colnames(target_df_aligned)
    ) {
      stop(
        "该结果对象没有可用的 group_col。",
        call. = FALSE
      )
    }

    group_vector <- target_df_aligned[[group_col]]

    target_sub <- target_df_aligned[
      !is.na(group_vector) &
        group_vector == group_name,
      ,
      drop = FALSE
    ]
  }

  if (nrow(target_sub) == 0L) {
    stop(
      "指定 group 中没有样本。",
      call. = FALSE
    )
  }

  samples_use <- target_sub[[sample_col]]

  scatter_df <- data.frame(
    sample = samples_use,
    feature_value = as.numeric(
      feature_mat[
        feature_name,
        samples_use
      ]
    ),
    target_value = target_sub[[target_name]],
    stringsAsFactors = FALSE
  )

  complete_idx <- (
    stats::complete.cases(
      scatter_df$feature_value,
      scatter_df$target_value
    ) &
      is.finite(scatter_df$feature_value) &
      is.finite(scatter_df$target_value)
  )

  scatter_df <- scatter_df[
    complete_idx,
    ,
    drop = FALSE
  ]

  if (nrow(scatter_df) < min_n) {
    stop(
      "完整样本数小于 min_n，无法绘图。",
      call. = FALSE
    )
  }

  if (
    length(unique(scatter_df$feature_value)) < 2L ||
    length(unique(scatter_df$target_value)) < 2L
  ) {
    stop(
      "feature 或 target 为常数，无法计算相关性。",
      call. = FALSE
    )
  }

  cor_result <- suppressWarnings(
    stats::cor.test(
      x = scatter_df$feature_value,
      y = scatter_df$target_value,
      method = "spearman",
      exact = FALSE
    )
  )

  rho_value <- unname(cor_result$estimate)
  p_value <- cor_result$p.value

  annotation_text <- paste0(
    "Spearman rho = ",
    round(rho_value, 3),
    "\np = ",
    signif(p_value, 3),
    "\nn = ",
    nrow(scatter_df)
  )

  if (identical(plot_type, "scatter")) {

    p <- ggplot2::ggplot(
      scatter_df,
      ggplot2::aes(
        x = target_value,
        y = feature_value
      )
    ) +
      ggplot2::geom_point(
        size = 2.5,
        alpha = 0.85
      ) +
      ggplot2::annotate(
        geom = "text",
        x = Inf,
        y = Inf,
        hjust = 1.1,
        vjust = 1.3,
        label = annotation_text
      ) +
      ggplot2::labs(
        title = paste0(
          feature_name,
          " vs ",
          target_name
        ),
        subtitle = paste0(
          "Group: ",
          group_name
        ),
        x = target_name,
        y = paste0(
          feature_name,
          " value"
        )
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold"
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5
        )
      )

    if (isTRUE(add_lm)) {
      p <- p +
        ggplot2::geom_smooth(
          method = "lm",
          se = TRUE
        )
    }

  } else {

    scatter_df$target_level <- factor(
      scatter_df$target_value,
      levels = sort(
        unique(scatter_df$target_value)
      ),
      ordered = TRUE
    )

    p <- ggplot2::ggplot(
      scatter_df,
      ggplot2::aes(
        x = target_level,
        y = feature_value
      )
    ) +
      ggplot2::geom_boxplot(
        width = 0.6,
        outlier.shape = NA
      ) +
      ggplot2::geom_jitter(
        width = 0.12,
        height = 0,
        size = 2,
        alpha = 0.65
      ) +
      ggplot2::annotate(
        geom = "text",
        x = Inf,
        y = Inf,
        hjust = 1.1,
        vjust = 1.3,
        label = annotation_text
      ) +
      ggplot2::labs(
        title = paste0(
          feature_name,
          " vs ",
          target_name
        ),
        subtitle = paste0(
          "Group: ",
          group_name
        ),
        x = target_name,
        y = paste0(
          feature_name,
          " value"
        )
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold"
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5
        )
      )
  }

  p
}


# ------------------------------------------------------------
# Print method for returned workflow objects
# ------------------------------------------------------------

print.spearman_workflow_result <- function(
    x,
    ...
) {

  cat(
    "<spearman_workflow_result>\n",
    sep = ""
  )

  cat(
    "Features: ",
    length(x$features_used),
    "\n",
    sep = ""
  )

  cat(
    "Targets: ",
    paste(x$targets_used, collapse = ", "),
    "\n",
    sep = ""
  )

  analysis_groups <- unique(
    x$results$group
  )

  cat(
    "Analyses: ",
    paste(analysis_groups, collapse = ", "),
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


# ============================================================
# Minimal usage example
# ============================================================
#
# source("spearman-correlation/R/run_spearman_workflow.R")
#
# result <- run_spearman_workflow(
#   feature_df = feature_df,
#   target_df = target_df,
#   feature_col = "feature",
#   sample_col = "Sample",
#   selected_features = NULL,
#   target_cols = c("Score_A"),
#   run_all_samples = TRUE,
#   run_group_analysis = TRUE,
#   group_col = "Group",
#   selected_groups = c("IA", "ENEG"),
#   min_n = 5,
#   outdir = "results/example_spearman"
# )
#
# result$results
# result$summary
#
# p_gene <- plot_spearman_feature(
#   workflow_result = result,
#   feature_name = "Gene_Pos",
#   target_name = "Score_A",
#   group_name = "IA",
#   plot_type = "scatter",
#   add_lm = TRUE
# )
#
# print(p_gene)
# ============================================================
