# Hallmark NES summary post-GSEA visualization module

plot_hallmark_nes_summary <- function(
    gsea_result,
    pathway_file,
    sort_by = "NES",
    color_by = "NES",
    size_by = "setSize",
    output_prefix = "hallmark_NES_summary",
    width = 7,
    height = 5,
    dpi = 300) {
  required_packages <- c("yaml", "dplyr", "ggplot2")
  missing_packages <- required_packages[!vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )]
  if (length(missing_packages) > 0) {
    stop(
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.data.frame(gsea_result)) {
    stop("gsea_result must be a data.frame.", call. = FALSE)
  }
  if (length(pathway_file) != 1L || !is.character(pathway_file) ||
      !file.exists(pathway_file)) {
    stop("pathway_file must be an existing YAML file path.", call. = FALSE)
  }
  if (length(output_prefix) != 1L || !is.character(output_prefix) ||
      !nzchar(output_prefix)) {
    stop("output_prefix must be a non-empty file path prefix.", call. = FALSE)
  }
  if (length(width) != 1L || !is.numeric(width) || !is.finite(width) ||
      width <= 0) {
    stop("width must be one positive finite number.", call. = FALSE)
  }
  if (length(height) != 1L || !is.numeric(height) || !is.finite(height) ||
      height <= 0) {
    stop("height must be one positive finite number.", call. = FALSE)
  }
  if (length(dpi) != 1L || !is.numeric(dpi) || !is.finite(dpi) || dpi <= 0) {
    stop("dpi must be one positive finite number.", call. = FALSE)
  }
  if (length(sort_by) != 1L || !sort_by %in% c("NES", "p.adjust")) {
    stop('sort_by must be either "NES" or "p.adjust".', call. = FALSE)
  }
  if (length(color_by) != 1L || !color_by %in% c("NES", "p.adjust", "none")) {
    stop('color_by must be "NES", "p.adjust", or "none".', call. = FALSE)
  }
  if (length(size_by) != 1L || !size_by %in% c("setSize", "none")) {
    stop('size_by must be either "setSize" or "none".', call. = FALSE)
  }

  required_columns <- c("Description", "NES")
  missing_columns <- setdiff(required_columns, names(gsea_result))
  if (length(missing_columns) > 0L) {
    stop(
      "gsea_result is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (sort_by == "p.adjust" && !"p.adjust" %in% names(gsea_result)) {
    stop(
      'sort_by = "p.adjust" requires a p.adjust column in gsea_result.',
      call. = FALSE
    )
  }
  if (color_by == "p.adjust" && !"p.adjust" %in% names(gsea_result)) {
    stop(
      'color_by = "p.adjust" requires a p.adjust column in gsea_result.',
      call. = FALSE
    )
  }

  pathway_config <- yaml::read_yaml(pathway_file)
  if (is.null(pathway_config$hallmark_set) ||
      !is.character(pathway_config$hallmark_set) ||
      length(pathway_config$hallmark_set) == 0L) {
    stop(
      'The YAML file must define a non-empty character vector named "hallmark_set".',
      call. = FALSE
    )
  }
  hallmark_set <- unique(as.character(pathway_config$hallmark_set))

  # Prefer the stable ID column. Description is used only as an exact-
  # match compatibility fallback when an input has no ID column at all.
  if ("ID" %in% names(gsea_result)) {
    match_values <- as.character(gsea_result$ID)
    match_column <- "ID"
  } else {
    warning(
      "Column `ID` not found; using exact Description matching as a compatibility fallback.",
      call. = FALSE
    )
    match_values <- as.character(gsea_result$Description)
    match_column <- "Description"
  }

  matched_pathways <- intersect(hallmark_set, unique(match_values))
  missing_pathways <- setdiff(hallmark_set, matched_pathways)
  cat("requested pathway count: ", length(hallmark_set), "\n", sep = "")
  cat("matched pathway count: ", length(matched_pathways), "\n", sep = "")
  cat("missing pathway count: ", length(missing_pathways), "\n", sep = "")

  if (length(missing_pathways) > 0L) {
    warning(
      "Requested Hallmark pathways missing from GSEA result (matched by ",
      match_column,
      "): ",
      paste(missing_pathways, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(matched_pathways) == 0L) {
    stop(
      "No requested Hallmark pathways were found in the GSEA result.",
      call. = FALSE
    )
  }

  plot_data <- gsea_result[match_values %in% hallmark_set, , drop = FALSE]
  plot_data$NES <- suppressWarnings(as.numeric(as.character(plot_data$NES)))
  valid_nes <- is.finite(plot_data$NES)
  if (any(!valid_nes)) {
    warning(
      sum(!valid_nes),
      " record(s) with NA or non-finite NES were removed before plotting.",
      call. = FALSE
    )
    plot_data <- plot_data[valid_nes, , drop = FALSE]
  }
  if (nrow(plot_data) == 0L) {
    stop("No finite NES values remain after filtering.", call. = FALSE)
  }

  if (sort_by == "NES") {
    # The first row is the highest NES and is intended to appear at the top.
    plot_data <- dplyr::arrange(plot_data, dplyr::desc(.data$NES))
  } else {
    plot_data <- dplyr::arrange(
      plot_data,
      suppressWarnings(as.numeric(as.character(.data$p.adjust)))
    )
  }

  original_description <- as.character(plot_data$Description)
  if ("ID" %in% names(plot_data)) {
    pathway_ids <- as.character(plot_data$ID)
  } else {
    pathway_ids <- original_description
    plot_data$ID <- pathway_ids
  }

  # Keep readable descriptions when supplied; only transform ID-like labels
  # for display. Matching identity remains in the ID column.
  display_description <- original_description
  use_id_label <- is.na(display_description) |
    !nzchar(display_description) |
    display_description == pathway_ids |
    startsWith(display_description, "HALLMARK_")
  display_description[use_id_label] <- pathway_ids[use_id_label]
  display_description <- sub("^HALLMARK_", "", display_description)
  display_description <- gsub("_", " ", display_description, fixed = TRUE)
  plot_data$Description_original <- original_description
  plot_data$Description <- display_description

  # Reverse levels so the first row after sorting is displayed at the top.
  plot_data$Description <- factor(
    plot_data$Description,
    levels = rev(unique(display_description))
  )

  if (size_by == "setSize" && !"setSize" %in% names(plot_data)) {
    warning(
      "size_by = \"setSize\" requested but setSize is absent; using a fixed point size.",
      call. = FALSE
    )
    size_by <- "none"
  }
  if (size_by == "setSize") {
    plot_data$.size_value <- suppressWarnings(
      as.numeric(as.character(plot_data$setSize))
    )
  }

  if (color_by == "NES") {
    plot_data$.color_value <- plot_data$NES
  } else if (color_by == "p.adjust") {
    p_adjust_value <- suppressWarnings(
      as.numeric(as.character(plot_data$p.adjust))
    )
    # Clamp zero p.adjust values to the smallest positive double before log10.
    safe_p_adjust <- ifelse(
      is.na(p_adjust_value) | !is.finite(p_adjust_value),
      NA_real_,
      pmax(p_adjust_value, .Machine$double.xmin)
    )
    plot_data$.color_value <- -log10(safe_p_adjust)
  }

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$NES, y = .data$Description))
  if (color_by != "none") {
    if (size_by == "setSize") {
      plot <- plot + ggplot2::geom_point(
        ggplot2::aes(color = .data$.color_value, size = .data$.size_value)
      )
    } else {
      plot <- plot + ggplot2::geom_point(ggplot2::aes(color = .data$.color_value))
    }
  } else if (size_by == "setSize") {
    plot <- plot + ggplot2::geom_point(ggplot2::aes(size = .data$.size_value))
  } else {
    plot <- plot + ggplot2::geom_point(size = 3)
  }

  plot <- plot +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      x = "Normalized enrichment score (NES)",
      y = NULL,
      color = switch(
        color_by,
        NES = "NES",
        p.adjust = "-log10(p.adjust)",
        NULL
      ),
      size = "setSize"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (color_by == "NES") {
    plot <- plot + ggplot2::scale_color_gradient2(midpoint = 0)
  } else if (color_by == "p.adjust") {
    plot <- plot + ggplot2::scale_color_viridis_c()
  }
  if (size_by == "setSize") {
    plot <- plot + ggplot2::scale_size_continuous()
  }

  output_dir <- dirname(output_prefix)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Unable to create output directory: ", output_dir, call. = FALSE)
  }

  ggplot2::ggsave(
    filename = paste0(output_prefix, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
  ggplot2::ggsave(
    filename = paste0(output_prefix, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )

  # Save the sorted, finite data actually used by the plot.
  plot_data_for_tsv <- plot_data
  plot_data_for_tsv$Description <- as.character(plot_data_for_tsv$Description)
  output_columns <- intersect(
    c("ID", "Description", "Description_original", "NES", "pvalue",
      "p.adjust", "qvalue", "setSize"),
    names(plot_data_for_tsv)
  )
  utils::write.table(
    plot_data_for_tsv[, output_columns, drop = FALSE],
    file = paste0(output_prefix, "_plot_data.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )

  invisible(plot)
}
