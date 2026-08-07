# Result screening, summaries, optional quadrant analysis and validation.

apply_spearman_thresholds <- function(
    cor_results,
    r_cutoff,
    p_cutoff,
    padj_cutoff,
    include_positive = TRUE,
    include_negative = TRUE
) {
  results <- cor_results

  results$significant_by_p <- (
    results$status == "OK" &
      !is.na(results$rho) &
      !is.na(results$pvalue) &
      results$abs_rho >= r_cutoff &
      results$pvalue < p_cutoff
  )

  results$significant_by_padj <- (
    results$status == "OK" &
      !is.na(results$rho) &
      !is.na(results$padj) &
      results$abs_rho >= r_cutoff &
      results$padj < padj_cutoff
  )

  significant_p <- results[
    results$significant_by_p,
    ,
    drop = FALSE
  ]
  significant_p <- significant_p[
    order(
      significant_p$group,
      significant_p$target,
      significant_p$pvalue
    ),
    ,
    drop = FALSE
  ]

  significant_padj <- results[
    results$significant_by_padj,
    ,
    drop = FALSE
  ]
  significant_padj <- significant_padj[
    order(
      significant_padj$group,
      significant_padj$target,
      significant_padj$padj
    ),
    ,
    drop = FALSE
  ]

  directional <- list()

  if (isTRUE(include_positive)) {
    directional$positive_p <- significant_p[
      significant_p$direction == "Positive",
      ,
      drop = FALSE
    ]
    directional$positive_padj <- significant_padj[
      significant_padj$direction == "Positive",
      ,
      drop = FALSE
    ]
  }

  if (isTRUE(include_negative)) {
    directional$negative_p <- significant_p[
      significant_p$direction == "Negative",
      ,
      drop = FALSE
    ]
    directional$negative_padj <- significant_padj[
      significant_padj$direction == "Negative",
      ,
      drop = FALSE
    ]
  }

  list(
    results = results,
    significant_p = significant_p,
    significant_padj = significant_padj,
    directional = directional
  )
}

build_spearman_summary <- function(results) {
  combinations <- unique(
    results[, c("group", "target"), drop = FALSE]
  )

  rows <- vector("list", nrow(combinations))

  for (i in seq_len(nrow(combinations))) {
    group_name <- combinations$group[i]
    target_name <- combinations$target[i]
    one <- results[
      results$group == group_name &
        results$target == target_name,
      ,
      drop = FALSE
    ]

    rows[[i]] <- data.frame(
      group = group_name,
      target = target_name,
      total_features = nrow(one),
      valid_results = sum(one$status == "OK", na.rm = TRUE),
      significant_by_p = sum(
        one$significant_by_p,
        na.rm = TRUE
      ),
      significant_by_padj = sum(
        one$significant_by_padj,
        na.rm = TRUE
      ),
      positive_by_p = sum(
        one$significant_by_p &
          one$direction == "Positive",
        na.rm = TRUE
      ),
      negative_by_p = sum(
        one$significant_by_p &
          one$direction == "Negative",
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }

  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  summary
}

run_spearman_quadrant <- function(results, quadrant, logger) {
  if (!isTRUE(quadrant$enabled)) {
    return(NULL)
  }

  select_axis <- function(target_name, group_name, prefix) {
    one <- results[
      results$target == target_name &
        results$group == group_name,
      ,
      drop = FALSE
    ]

    if (nrow(one) == 0L) {
      available <- unique(
        paste(results$target, results$group, sep = " @ ")
      )
      stop(
        "Quadrant axis ",
        prefix,
        " refers to missing analysis '",
        target_name,
        " @ ",
        group_name,
        "'. Available analyses: ",
        paste(available, collapse = ", "),
        call. = FALSE
      )
    }

    if (anyDuplicated(one$feature) > 0L) {
      stop(
        "Quadrant axis ",
        prefix,
        " contains duplicate features.",
        call. = FALSE
      )
    }

    one <- one[
      one$status == "OK" &
        is.finite(one$rho) &
        is.finite(one$pvalue),
      c(
        "feature",
        "rho",
        "pvalue",
        "padj",
        "n",
        "significant_by_p",
        "significant_by_padj"
      ),
      drop = FALSE
    ]

    colnames(one)[-1L] <- paste0(
      prefix,
      "_",
      colnames(one)[-1L]
    )
    one
  }

  x <- select_axis(
    quadrant$x_target,
    quadrant$x_group,
    "x"
  )
  y <- select_axis(
    quadrant$y_target,
    quadrant$y_group,
    "y"
  )

  joined <- merge(
    x,
    y,
    by = "feature",
    all = FALSE,
    sort = FALSE
  )

  if (nrow(joined) == 0L) {
    stop(
      "Quadrant axes have no common valid features.",
      call. = FALSE
    )
  }

  x_order <- match(joined$feature, x$feature)
  joined <- joined[order(x_order), , drop = FALSE]
  rownames(joined) <- NULL

  if (identical(quadrant$significance, "p")) {
    joined$x_significant <- joined$x_significant_by_p
    joined$y_significant <- joined$y_significant_by_p
  } else if (identical(quadrant$significance, "padj")) {
    joined$x_significant <- joined$x_significant_by_padj
    joined$y_significant <- joined$y_significant_by_padj
  } else {
    joined$x_significant <- TRUE
    joined$y_significant <- TRUE
  }

  joined$quadrant <- ifelse(
    joined$x_rho > 0 & joined$y_rho > 0,
    "Q1_positive_positive",
    ifelse(
      joined$x_rho < 0 & joined$y_rho > 0,
      "Q2_negative_positive",
      ifelse(
        joined$x_rho < 0 & joined$y_rho < 0,
        "Q3_negative_negative",
        ifelse(
          joined$x_rho > 0 & joined$y_rho < 0,
          "Q4_positive_negative",
          "Axis"
        )
      )
    )
  )

  joined$selected <- if (isTRUE(
    quadrant$require_significant_both
  )) {
    joined$x_significant & joined$y_significant
  } else {
    rep(TRUE, nrow(joined))
  }

  joined$x_target <- quadrant$x_target
  joined$x_group <- quadrant$x_group
  joined$y_target <- quadrant$y_target
  joined$y_group <- quadrant$y_group

  front <- c(
    "feature",
    "quadrant",
    "selected",
    "x_target",
    "x_group",
    "y_target",
    "y_group"
  )
  joined <- joined[
    ,
    c(front, setdiff(colnames(joined), front)),
    drop = FALSE
  ]

  logger(
    "INFO",
    "Quadrant analysis produced ",
    nrow(joined),
    " common valid features; selected ",
    sum(joined$selected),
    "."
  )

  joined
}

run_spearman_manual_validation <- function(
    workflow_result,
    validation,
    logger
) {
  if (!isTRUE(validation$enabled)) {
    return(NULL)
  }

  items <- validation$items
  if (length(items) == 0L) {
    logger(
      "WARN",
      "Manual validation is enabled but no items are configured."
    )
    return(NULL)
  }

  feature_mat <- workflow_result$feature_mat
  target_df <- workflow_result$target_df_aligned
  results <- workflow_result$results
  sample_col <- workflow_result$parameters$sample_col
  group_col <- workflow_result$parameters$group_col

  rows <- vector("list", length(items))

  for (i in seq_along(items)) {
    item <- items[[i]]

    if (!item$feature %in% rownames(feature_mat)) {
      stop(
        "Manual validation feature not found: ",
        item$feature,
        call. = FALSE
      )
    }
    if (!item$target %in% colnames(target_df)) {
      stop(
        "Manual validation target column not found: ",
        item$target,
        call. = FALSE
      )
    }

    workflow_row <- results[
      results$feature == item$feature &
        results$target == item$target &
        results$group == item$group,
      ,
      drop = FALSE
    ]

    if (nrow(workflow_row) != 1L) {
      stop(
        "Manual validation expected exactly one workflow row for ",
        item$feature,
        " / ",
        item$target,
        " / ",
        item$group,
        "; found ",
        nrow(workflow_row),
        ".",
        call. = FALSE
      )
    }
    if (!identical(workflow_row$status, "OK")) {
      stop(
        "Manual validation item is not calculable; workflow status is ",
        workflow_row$status,
        " for ",
        item$feature,
        " / ",
        item$target,
        " / ",
        item$group,
        ".",
        call. = FALSE
      )
    }

    target_sub <- target_df
    if (!identical(item$group, "All")) {
      if (is.null(group_col) || !group_col %in% colnames(target_df)) {
        stop(
          "Manual validation group requested but no group column is available: ",
          item$group,
          call. = FALSE
        )
      }
      group_value <- target_df[[group_col]]
      target_sub <- target_df[
        !is.na(group_value) &
          group_value == item$group,
        ,
        drop = FALSE
      ]
    }

    samples <- target_sub[[sample_col]]
    x <- as.numeric(feature_mat[item$feature, samples])
    y <- target_sub[[item$target]]
    complete <- (
      stats::complete.cases(x, y) &
        is.finite(x) &
        is.finite(y)
    )

    manual <- suppressWarnings(
      stats::cor.test(
        x = x[complete],
        y = y[complete],
        method = "spearman",
        exact = FALSE
      )
    )

    manual_rho <- unname(manual$estimate)
    manual_p <- manual$p.value
    rho_match <- isTRUE(all.equal(
      manual_rho,
      workflow_row$rho,
      tolerance = 1e-12
    ))
    p_match <- isTRUE(all.equal(
      manual_p,
      workflow_row$pvalue,
      tolerance = 1e-12
    ))
    n_match <- identical(
      as.integer(sum(complete)),
      as.integer(workflow_row$n)
    )

    rows[[i]] <- data.frame(
      feature = item$feature,
      target = item$target,
      group = item$group,
      workflow_rho = workflow_row$rho,
      manual_rho = manual_rho,
      workflow_pvalue = workflow_row$pvalue,
      manual_pvalue = manual_p,
      workflow_n = workflow_row$n,
      manual_n = sum(complete),
      rho_match = rho_match,
      pvalue_match = p_match,
      n_match = n_match,
      passed = rho_match && p_match && n_match,
      stringsAsFactors = FALSE
    )
  }

  validation_result <- do.call(rbind, rows)
  rownames(validation_result) <- NULL

  if (!all(validation_result$passed)) {
    stop(
      "Manual validation failed for one or more configured items.",
      call. = FALSE
    )
  }

  logger(
    "INFO",
    "Manual validation passed for ",
    nrow(validation_result),
    " item(s)."
  )
  validation_result
}
