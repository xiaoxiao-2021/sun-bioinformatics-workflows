# Optional plots. Statistical calculations remain in analysis_core.R.

.spearman_check_plot_packages <- function(label_top_n_each = 0L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Plotting is enabled but package 'ggplot2' is not installed.",
      call. = FALSE
    )
  }

  if (label_top_n_each > 0L &&
      !requireNamespace("ggrepel", quietly = TRUE)) {
    stop(
      "Overview labels are enabled but package 'ggrepel' is not installed.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

plot_spearman_overview <- function(
    result_df,
    target_name,
    group_name,
    r_cutoff = 0.5,
    p_cutoff = 0.05,
    label_top_n_each = 0L
) {
  label_top_n_each <- as.integer(label_top_n_each)
  .spearman_check_plot_packages(label_top_n_each)

  required <- c(
    "feature",
    "target",
    "group",
    "rho",
    "pvalue",
    "status"
  )
  missing <- setdiff(required, colnames(result_df))
  if (length(missing) > 0L) {
    stop(
      "result_df is missing columns required for overview plotting: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  plot_df <- result_df[
    result_df$target == target_name &
      result_df$group == group_name &
      result_df$status == "OK" &
      is.finite(result_df$rho) &
      is.finite(result_df$pvalue),
    ,
    drop = FALSE
  ]
  if (nrow(plot_df) == 0L) {
    stop(
      "No valid results are available for overview plot target '",
      target_name,
      "', group '",
      group_name,
      "'.",
      call. = FALSE
    )
  }

  plot_df$neg_log10_p <- -log10(
    pmax(plot_df$pvalue, .Machine$double.xmin)
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
    levels = c("Negative", "Not significant", "Positive")
  )

  title <- if (identical(group_name, "All")) {
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
    ggplot2::geom_point(size = 1.7, alpha = 0.8) +
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
      title = title,
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
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (label_top_n_each > 0L) {
    positive <- plot_df[
      plot_df$correlation_class == "Positive",
      ,
      drop = FALSE
    ]
    negative <- plot_df[
      plot_df$correlation_class == "Negative",
      ,
      drop = FALSE
    ]

    positive <- utils::head(
      positive[
        order(positive$pvalue, -positive$rho),
        ,
        drop = FALSE
      ],
      label_top_n_each
    )
    negative <- utils::head(
      negative[
        order(negative$pvalue, negative$rho),
        ,
        drop = FALSE
      ],
      label_top_n_each
    )
    labels <- rbind(positive, negative)

    if (nrow(labels) > 0L) {
      p <- p + ggrepel::geom_text_repel(
        data = labels,
        ggplot2::aes(label = feature),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.3,
        point.padding = 0.2,
        show.legend = FALSE,
        seed = 1
      )
    }
  }

  p
}

plot_spearman_feature <- function(
    workflow_result,
    feature_name,
    target_name,
    group_name = "All",
    plot_type = c("scatter", "ordinal"),
    add_lm = TRUE
) {
  plot_type <- match.arg(plot_type)
  .spearman_check_plot_packages(0L)

  if (!inherits(workflow_result, "spearman_workflow_result")) {
    stop(
      "workflow_result must come from run_spearman_workflow().",
      call. = FALSE
    )
  }

  feature_mat <- workflow_result$feature_mat
  target_df <- workflow_result$target_df_aligned
  sample_col <- workflow_result$parameters$sample_col
  group_col <- workflow_result$parameters$group_col
  min_n <- workflow_result$parameters$min_n

  if (!feature_name %in% rownames(feature_mat)) {
    stop("Feature not found for plotting: ", feature_name, call. = FALSE)
  }
  if (!target_name %in% colnames(target_df)) {
    stop(
      "Target column not found for plotting: ",
      target_name,
      call. = FALSE
    )
  }

  target_sub <- target_df
  if (!identical(group_name, "All")) {
    if (is.null(group_col) || !group_col %in% colnames(target_df)) {
      stop(
        "A group was requested for plotting but no group column is available.",
        call. = FALSE
      )
    }
    groups <- target_df[[group_col]]
    target_sub <- target_df[
      !is.na(groups) & groups == group_name,
      ,
      drop = FALSE
    ]
  }
  if (nrow(target_sub) == 0L) {
    stop(
      "No samples found for plot group: ",
      group_name,
      call. = FALSE
    )
  }

  samples <- target_sub[[sample_col]]
  plot_df <- data.frame(
    sample = samples,
    feature_value = as.numeric(
      feature_mat[feature_name, samples]
    ),
    target_value = target_sub[[target_name]],
    stringsAsFactors = FALSE
  )
  complete <- (
    stats::complete.cases(
      plot_df$feature_value,
      plot_df$target_value
    ) &
      is.finite(plot_df$feature_value) &
      is.finite(plot_df$target_value)
  )
  plot_df <- plot_df[complete, , drop = FALSE]

  if (nrow(plot_df) < min_n) {
    stop(
      "Plot '",
      feature_name,
      " / ",
      target_name,
      " / ",
      group_name,
      "' has ",
      nrow(plot_df),
      " complete samples, fewer than min_n = ",
      min_n,
      ".",
      call. = FALSE
    )
  }
  if (length(unique(plot_df$feature_value)) < 2L ||
      length(unique(plot_df$target_value)) < 2L) {
    stop(
      "Feature or target is constant for plot '",
      feature_name,
      " / ",
      target_name,
      " / ",
      group_name,
      "'.",
      call. = FALSE
    )
  }

  test <- suppressWarnings(
    stats::cor.test(
      plot_df$feature_value,
      plot_df$target_value,
      method = "spearman",
      exact = FALSE
    )
  )
  annotation <- paste0(
    "Spearman rho = ",
    round(unname(test$estimate), 3),
    "\np = ",
    signif(test$p.value, 3),
    "\nn = ",
    nrow(plot_df)
  )

  if (identical(plot_type, "scatter")) {
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = target_value, y = feature_value)
    ) +
      ggplot2::geom_point(size = 2.5, alpha = 0.85) +
      ggplot2::annotate(
        "text",
        x = Inf,
        y = Inf,
        hjust = 1.1,
        vjust = 1.3,
        label = annotation
      )

    if (isTRUE(add_lm)) {
      p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE)
    }
  } else {
    plot_df$target_level <- factor(
      plot_df$target_value,
      levels = sort(unique(plot_df$target_value)),
      ordered = TRUE
    )
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = target_level, y = feature_value)
    ) +
      ggplot2::geom_boxplot(width = 0.6, outlier.shape = NA) +
      ggplot2::geom_point(
        position = ggplot2::position_jitter(
          width = 0.12,
          height = 0,
          seed = 1
        ),
        size = 2,
        alpha = 0.65
      ) +
      ggplot2::annotate(
        "text",
        x = Inf,
        y = Inf,
        hjust = 1.1,
        vjust = 1.3,
        label = annotation
      )
  }

  p +
    ggplot2::labs(
      title = paste0(feature_name, " vs ", target_name),
      subtitle = paste0("Group: ", group_name),
      x = target_name,
      y = paste0(feature_name, " value")
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
}

plot_spearman_quadrant <- function(quadrant_df) {
  .spearman_check_plot_packages(0L)

  ggplot2::ggplot(
    quadrant_df,
    ggplot2::aes(
      x = x_rho,
      y = y_rho,
      color = quadrant,
      alpha = selected
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_alpha_manual(
      values = c("TRUE" = 0.9, "FALSE" = 0.25),
      guide = "none"
    ) +
    ggplot2::coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
    ggplot2::labs(
      title = "Two-analysis Spearman quadrant",
      subtitle = paste0(
        unique(quadrant_df$x_target),
        " @ ",
        unique(quadrant_df$x_group),
        " versus ",
        unique(quadrant_df$y_target),
        " @ ",
        unique(quadrant_df$y_group)
      ),
      x = "X-axis Spearman rho",
      y = "Y-axis Spearman rho",
      color = "Quadrant"
    ) +
    ggplot2::theme_bw(base_size = 13)
}

create_spearman_plots <- function(workflow_result, config, logger) {
  overview <- list()
  single <- list()
  quadrant_plot <- NULL

  if (isTRUE(config$plots$overview$enabled)) {
    .spearman_check_plot_packages(
      config$plots$overview$label_top_n_each
    )
    combinations <- unique(
      workflow_result$results[
        ,
        c("group", "target"),
        drop = FALSE
      ]
    )

    labels <- paste(
      combinations$group,
      combinations$target,
      sep = "__"
    )
    .spearman_assert_unique_filenames(
      labels,
      "overview plots"
    )

    for (i in seq_len(nrow(combinations))) {
      group_name <- combinations$group[i]
      target_name <- combinations$target[i]
      valid <- workflow_result$results[
        workflow_result$results$group == group_name &
          workflow_result$results$target == target_name &
          workflow_result$results$status == "OK" &
          is.finite(workflow_result$results$rho) &
          is.finite(workflow_result$results$pvalue),
        ,
        drop = FALSE
      ]

      if (nrow(valid) == 0L) {
        logger(
          "WARN",
          "Skipping overview plot with no valid results: ",
          target_name,
          " @ ",
          group_name
        )
        next
      }

      name <- paste0(
        .spearman_safe_filename(group_name),
        "__",
        .spearman_safe_filename(target_name)
      )
      overview[[name]] <- plot_spearman_overview(
        workflow_result$results,
        target_name,
        group_name,
        workflow_result$parameters$r_cutoff,
        workflow_result$parameters$p_cutoff,
        config$plots$overview$label_top_n_each
      )
    }
  }

  if (isTRUE(config$plots$single_feature$enabled)) {
    .spearman_check_plot_packages(0L)
    items <- config$plots$single_feature$items

    for (i in seq_along(items)) {
      item <- items[[i]]
      name <- paste0(
        .spearman_safe_filename(item$group),
        "__",
        .spearman_safe_filename(item$feature),
        "__",
        .spearman_safe_filename(item$target),
        "__",
        item$type
      )
      if (name %in% names(single)) {
        stop(
          "Duplicate single-feature plot output name: ",
          name,
          call. = FALSE
        )
      }
      single[[name]] <- plot_spearman_feature(
        workflow_result,
        item$feature,
        item$target,
        item$group,
        item$type,
        item$add_lm
      )
    }
  }

  if (
    !is.null(workflow_result$quadrant) &&
      isTRUE(config$plots$quadrant$enabled)
  ) {
    quadrant_plot <- plot_spearman_quadrant(
      workflow_result$quadrant
    )
  }

  logger(
    "INFO",
    "Created ",
    length(overview),
    " overview plot(s), ",
    length(single),
    " single-feature plot(s)."
  )

  list(
    overview = overview,
    single_feature = single,
    quadrant = quadrant_plot
  )
}

.spearman_save_plot_set <- function(
    plots,
    directory,
    suffix,
    settings
) {
  if (length(plots) == 0L) {
    return(invisible(character(0)))
  }

  paths <- character(0)
  for (name in names(plots)) {
    for (format in settings$formats) {
      path <- file.path(
        directory,
        paste0(name, suffix, ".", format)
      )
      ggplot2::ggsave(
        filename = path,
        plot = plots[[name]],
        width = settings$width,
        height = settings$height,
        dpi = settings$dpi,
        units = "in"
      )
      paths <- c(paths, path)
    }
  }
  invisible(paths)
}

save_spearman_plots <- function(plots, config, logger) {
  has_plots <- (
    length(plots$overview) > 0L ||
      length(plots$single_feature) > 0L ||
      !is.null(plots$quadrant)
  )
  if (!has_plots) {
    logger("INFO", "Saved 0 plot file(s).")
    return(invisible(character(0)))
  }

  plot_dir <- file.path(config$output_dir, "plots")
  if (!dir.exists(plot_dir) &&
      !dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create plot directory: ", plot_dir, call. = FALSE)
  }

  saved <- character(0)
  saved <- c(
    saved,
    .spearman_save_plot_set(
      plots$overview,
      plot_dir,
      "__correlation_overview",
      config$plots$overview
    ),
    .spearman_save_plot_set(
      plots$single_feature,
      plot_dir,
      "",
      config$plots$single_feature
    )
  )

  if (!is.null(plots$quadrant)) {
    settings <- config$plots$quadrant
    for (format in settings$formats) {
      path <- file.path(
        plot_dir,
        paste0("Spearman_quadrant.", format)
      )
      ggplot2::ggsave(
        filename = path,
        plot = plots$quadrant,
        width = settings$width,
        height = settings$height,
        dpi = settings$dpi,
        units = "in"
      )
      saved <- c(saved, path)
    }
  }

  logger("INFO", "Saved ", length(saved), " plot file(s).")
  invisible(saved)
}
