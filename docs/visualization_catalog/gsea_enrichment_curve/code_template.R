# GSEA enrichment curve post-analysis visualization module

plot_gsea_enrichment_curve <- function(
    pathway_id,
    gsea_result,
    ranked_genes,
    term2gene,
    gene_id_column,
    rank_column,
    case_group,
    control_group,
    gsea_exponent = 1,
    output_dir,
    width = 7,
    height = 5,
    dpi = 300) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Missing required R package: ggplot2", call. = FALSE)
  }
  if (!is.data.frame(gsea_result) || !is.data.frame(ranked_genes) ||
      !is.data.frame(term2gene)) {
    stop("gsea_result, ranked_genes, and term2gene must be data.frames.", call. = FALSE)
  }
  if (length(pathway_id) != 1L || !is.character(pathway_id) ||
      !nzchar(pathway_id)) {
    stop("pathway_id must be one non-empty pathway ID.", call. = FALSE)
  }
  if (length(gene_id_column) != 1L || !gene_id_column %in% names(ranked_genes)) {
    stop("gene_id_column is absent from ranked_genes.", call. = FALSE)
  }
  if (length(rank_column) != 1L || !rank_column %in% names(ranked_genes)) {
    stop("rank_column is absent from ranked_genes.", call. = FALSE)
  }
  if (!all(c("ID", "Description", "NES", "pvalue") %in% names(gsea_result))) {
    stop("gsea_result must contain ID, Description, NES, and pvalue.", call. = FALSE)
  }
  if (!all(c("TERM", "GENE") %in% names(term2gene))) {
    stop("term2gene must contain TERM and GENE columns.", call. = FALSE)
  }
  if (length(gsea_exponent) != 1L || !is.numeric(gsea_exponent) ||
      !is.finite(gsea_exponent) || gsea_exponent < 0) {
    stop("gsea_exponent must be one non-negative finite number.", call. = FALSE)
  }
  if (length(width) != 1L || !is.numeric(width) || !is.finite(width) || width <= 0 ||
      length(height) != 1L || !is.numeric(height) || !is.finite(height) || height <= 0 ||
      length(dpi) != 1L || !is.numeric(dpi) || !is.finite(dpi) || dpi <= 0) {
    stop("width, height, and dpi must be positive finite numbers.", call. = FALSE)
  }
  if (length(output_dir) != 1L || !is.character(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be a non-empty directory path.", call. = FALSE)
  }

  pathway_rows <- which(as.character(gsea_result$ID) == pathway_id)
  if (!length(pathway_rows)) {
    stop("Pathway is absent from gsea_result: ", pathway_id, call. = FALSE)
  }
  pathway_result <- gsea_result[pathway_rows[[1]], , drop = FALSE]

  # The original description is retained for metadata; only the display label
  # is transformed, so TERM/ID identity is never changed.
  original_description <- as.character(pathway_result$Description[[1]])
  display_label <- original_description
  if (is.na(display_label) || !nzchar(display_label) ||
      display_label == pathway_id || startsWith(display_label, "HALLMARK_")) {
    display_label <- pathway_id
  }
  display_label <- sub("^HALLMARK_", "", display_label)
  display_label <- gsub("_", " ", display_label, fixed = TRUE)

  ranked <- ranked_genes
  gene_ids <- as.character(ranked[[gene_id_column]])
  rank_metric <- suppressWarnings(as.numeric(as.character(ranked[[rank_column]])))
  valid_rows <- !is.na(gene_ids) & nzchar(gene_ids) & is.finite(rank_metric)
  if (any(!valid_rows)) {
    warning(
      sum(!valid_rows),
      " ranked gene row(s) with missing gene ID or non-finite rank were removed.",
      call. = FALSE
    )
    ranked <- ranked[valid_rows, , drop = FALSE]
    gene_ids <- gene_ids[valid_rows]
    rank_metric <- rank_metric[valid_rows]
  }
  if (!length(gene_ids)) {
    stop("No valid ranked genes remain after filtering.", call. = FALSE)
  }
  if (anyDuplicated(gene_ids)) {
    stop(
      "ranked_gene_file contains duplicate gene IDs; provide the gene-level ranking snapshot used by GSEA.",
      call. = FALSE
    )
  }
  ordered_rows <- order(rank_metric, decreasing = TRUE)
  gene_ids <- gene_ids[ordered_rows]
  rank_metric <- rank_metric[ordered_rows]

  term_genes <- unique(as.character(term2gene$GENE[
    as.character(term2gene$TERM) == pathway_id
  ]))
  term_genes <- term_genes[!is.na(term_genes) & nzchar(term_genes)]
  matched_gene_ids <- intersect(term_genes, gene_ids)
  if (!length(matched_gene_ids)) {
    warning(
      "No TERM2GENE genes overlap the ranked gene list for pathway: ", pathway_id,
      call. = FALSE
    )
    return(NULL)
  }

  n_genes <- length(gene_ids)
  is_hit <- gene_ids %in% matched_gene_ids
  n_hits <- sum(is_hit)
  if (n_hits >= n_genes) {
    stop("The pathway contains every ranked gene; miss penalty is undefined.", call. = FALSE)
  }

  # Reconstruct the weighted GSEA running sum for visualization only.
  # A hit contributes |rank statistic|^p normalized over all hits; a miss
  # contributes the uniform penalty 1 / (N - Nh). Thus the trajectory returns
  # to zero at the end up to floating-point error.
  hit_weights <- ifelse(is_hit, abs(rank_metric)^gsea_exponent, 0)
  total_hit_weight <- sum(hit_weights)
  if (!is.finite(total_hit_weight) || total_hit_weight <= 0) {
    stop("Cannot normalize weighted hit contributions for pathway: ", pathway_id, call. = FALSE)
  }
  running_increment <- ifelse(
    is_hit,
    hit_weights / total_hit_weight,
    -1 / (n_genes - n_hits)
  )
  running_es <- cumsum(running_increment)

  max_rank <- which.max(running_es)
  min_rank <- which.min(running_es)
  max_es <- running_es[[max_rank]]
  min_es <- running_es[[min_rank]]
  if (abs(max_es) >= abs(min_es)) {
    visual_es <- max_es
    visual_es_rank <- max_rank
  } else {
    visual_es <- min_es
    visual_es_rank <- min_rank
  }

  original_enrichment_score <- if ("enrichmentScore" %in% names(pathway_result)) {
    suppressWarnings(as.numeric(as.character(pathway_result$enrichmentScore[[1]])))
  } else {
    NA_real_
  }
  visual_es_difference <- visual_es - original_enrichment_score
  if (is.finite(original_enrichment_score)) {
    comparison_tolerance <- max(0.05, 0.1 * abs(original_enrichment_score))
    if (abs(visual_es_difference) > comparison_tolerance) {
      warning(
        "Reconstructed ES differs from the original GSEA enrichmentScore. ",
        "Check rank metric / exponent / TERM2GENE consistency.",
        call. = FALSE
      )
    }
  }

  curve_data <- data.frame(
    rank = seq_len(n_genes),
    gene_id = gene_ids,
    rank_metric = rank_metric,
    is_hit = is_hit,
    running_ES = running_es,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # The curve occupies the upper area; hit ticks and the ranked metric strip
  # occupy compact bands below it while sharing the same rank x-axis.
  es_scale <- max(abs(running_es), 1e-8)
  curve_top <- max(max(running_es), 0) + 0.20 * es_scale
  curve_bottom <- min(min(running_es), 0) - 0.12 * es_scale
  hit_base <- curve_bottom - 0.18 * es_scale
  hit_top <- hit_base + 0.12 * es_scale
  metric_zero <- hit_base - 0.22 * es_scale
  metric_height <- 0.16 * es_scale
  metric_scale <- max(abs(rank_metric), 1e-8)
  metric_end <- metric_zero + metric_height * rank_metric / metric_scale
  y_min <- metric_zero - metric_height - 0.25 * es_scale
  y_max <- curve_top

  curve_data$metric_zero <- metric_zero
  curve_data$metric_end <- metric_end
  # Add an explicit rank-zero origin to the plotted trajectory. The exported
  # curve data below still contains one row for each ranked gene position.
  plot_curve_data <- rbind(
    data.frame(
      rank = 0,
      gene_id = NA_character_,
      rank_metric = NA_real_,
      is_hit = FALSE,
      running_ES = 0,
      metric_zero = metric_zero,
      metric_end = metric_zero,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    curve_data
  )
  hit_data <- plot_curve_data[plot_curve_data$is_hit, , drop = FALSE]
  stats_label <- paste0(
    "NES = ", formatC(as.numeric(pathway_result$NES[[1]]), format = "g", digits = 4),
    "; P = ", format.pval(as.numeric(pathway_result$pvalue[[1]]), digits = 3, eps = 1e-04)
  )

  plot <- ggplot2::ggplot(plot_curve_data, ggplot2::aes(x = .data$rank, y = .data$running_ES)) +
    ggplot2::geom_line(linewidth = 0.7, color = "steelblue4") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    ggplot2::geom_vline(
      xintercept = visual_es_rank,
      linetype = "dotted",
      color = "grey40"
    ) +
    ggplot2::geom_point(
      data = plot_curve_data[plot_curve_data$rank == visual_es_rank, , drop = FALSE],
      ggplot2::aes(x = .data$rank, y = .data$running_ES),
      color = "firebrick3",
      size = 2
    ) +
    ggplot2::geom_segment(
      data = hit_data,
      ggplot2::aes(
        x = .data$rank,
        xend = .data$rank,
        y = hit_base,
        yend = hit_top
      ),
      inherit.aes = FALSE,
      linewidth = 0.35,
      color = "grey20"
    ) +
    ggplot2::geom_hline(yintercept = metric_zero, color = "grey75", linewidth = 0.3) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$rank,
        xend = .data$rank,
        y = metric_zero,
        yend = .data$metric_end
      ),
      linewidth = 0.25,
      alpha = 0.35,
      color = "grey35"
    ) +
    ggplot2::annotate(
      "text", x = 1, y = y_min + 0.08 * es_scale,
      label = paste0(case_group, " direction (+)"), hjust = 0, size = 3
    ) +
    ggplot2::annotate(
      "text", x = n_genes, y = y_min + 0.08 * es_scale,
      label = paste0(control_group, " direction (-)"), hjust = 1, size = 3
    ) +
    ggplot2::annotate(
      "text", x = 1, y = y_max, label = stats_label,
      hjust = 0, vjust = 1, size = 3.2
    ) +
    ggplot2::labs(
      title = display_label,
      subtitle = "Running ES, pathway hit positions, and ranked metric",
      x = "Rank position",
      y = "Running ES",
      caption = paste0("Hit genes: ", n_hits, " / ", length(term_genes),
                       " TERM2GENE genes; exponent p = ", gsea_exponent)
    ) +
    ggplot2::coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Unable to create output directory: ", output_dir, call. = FALSE)
  }

  file_stem <- sub("^HALLMARK_", "", pathway_id)
  file_stem <- gsub("[^A-Za-z0-9._-]+", "_", file_stem)
  png_file <- file.path(output_dir, paste0(file_stem, ".png"))
  pdf_file <- file.path(output_dir, paste0(file_stem, ".pdf"))
  data_file <- file.path(output_dir, paste0(file_stem, "_curve_data.tsv"))

  ggplot2::ggsave(png_file, plot = plot, width = width, height = height, units = "in", dpi = dpi)
  ggplot2::ggsave(pdf_file, plot = plot, width = width, height = height, units = "in")
  utils::write.table(
    curve_data[, c("rank", "gene_id", "rank_metric", "is_hit", "running_ES")],
    file = data_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )

  cat("Saved GSEA enrichment curve:\n", pathway_id, "\n", sep = "")
  cat("PNG: ", png_file, "\n", sep = "")
  cat("PDF: ", pdf_file, "\n", sep = "")
  cat("Data: ", data_file, "\n", sep = "")

  optional_result_value <- function(column) {
    if (column %in% names(pathway_result)) {
      suppressWarnings(as.numeric(as.character(pathway_result[[column]][[1]])))
    } else {
      NA_real_
    }
  }
  summary_row <- data.frame(
    ID = pathway_id,
    Description = original_description,
    NES = optional_result_value("NES"),
    pvalue = optional_result_value("pvalue"),
    p.adjust = optional_result_value("p.adjust"),
    qvalue = optional_result_value("qvalue"),
    setSize = optional_result_value("setSize"),
    matched_gene_count = n_hits,
    visual_ES = visual_es,
    visual_ES_rank = visual_es_rank,
    enrichmentScore = original_enrichment_score,
    visual_ES_difference = visual_es_difference,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  list(summary = summary_row, pathway_id = pathway_id)
}
