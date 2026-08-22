# Hallmark NES summary post-GSEA visualization module

plot_hallmark_nes_summary <- function(
    gsea_result,
    pathway_file,
    sort_by = "NES",
    x_variable = "NES",
    size_variable = "NES",
    size_transform = "abs",
    color_variable = "pvalue",
    color_transform = "-log10",
    color_low = "#e07e65",
    color_high = "#9b1f2a",
    label_wrap_width = 0,
    output_prefix = "hallmark_NES_summary",
    width = 7,
    height = 5,
    dpi = 300) {
  required_packages <- c("yaml", "ggplot2")
  missing_packages <- required_packages[!vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )]
  if (length(missing_packages) > 0L) {
    stop(
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.data.frame(gsea_result)) {
    stop("gsea_result must be a data.frame.", call. = FALSE)
  }
  if (!"pvalue" %in% colnames(gsea_result)) {
    stop("Required column pvalue missing", call. = FALSE)
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
  if (length(color_low) != 1L || !is.character(color_low) ||
      !nzchar(color_low) || length(color_high) != 1L ||
      !is.character(color_high) || !nzchar(color_high)) {
    stop("color_low and color_high must be non-empty color strings.", call. = FALSE)
  }
  if (length(label_wrap_width) != 1L || !is.numeric(label_wrap_width) ||
      !is.finite(label_wrap_width) || label_wrap_width < 0) {
    stop("label_wrap_width must be one non-negative finite number.", call. = FALSE)
  }

  mapping_config <- c(
    sort_by = sort_by,
    x_variable = x_variable,
    size_variable = size_variable,
    size_transform = size_transform,
    color_variable = color_variable,
    color_transform = color_transform
  )
  required_mapping <- c(
    sort_by = "NES",
    x_variable = "NES",
    size_variable = "NES",
    size_transform = "abs",
    color_variable = "pvalue",
    color_transform = "-log10"
  )
  if (!identical(unname(mapping_config), unname(required_mapping))) {
    stop(
      paste0(
        "This module requires sort_by=NES, x_variable=NES, ",
        "size_variable=NES, size_transform=abs, ",
        "color_variable=pvalue, and color_transform=-log10."
      ),
      call. = FALSE
    )
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

  # Prefer the stable ID column. Description is an exact-match compatibility
  # fallback only when an input has no ID column.
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

  plot_data$pvalue <- suppressWarnings(
    as.numeric(as.character(plot_data$pvalue))
  )
  invalid_pvalue <- !is.finite(plot_data$pvalue) |
    plot_data$pvalue < 0 |
    plot_data$pvalue > 1
  if (any(invalid_pvalue)) {
    stop(
      "pvalue must contain finite numeric values between 0 and 1.",
      call. = FALSE
    )
  }

  # Pathway display order is determined only by NES (highest at the top).
  plot_data <- plot_data[order(plot_data$NES, decreasing = TRUE), , drop = FALSE]

  original_description <- as.character(plot_data$Description)
  if ("ID" %in% names(plot_data)) {
    pathway_ids <- as.character(plot_data$ID)
  } else {
    pathway_ids <- original_description
    plot_data$ID <- pathway_ids
  }

  display_description <- original_description
  use_id_label <- is.na(display_description) |
    !nzchar(display_description) |
    display_description == pathway_ids |
    startsWith(display_description, "HALLMARK_")
  display_description[use_id_label] <- pathway_ids[use_id_label]
  display_description <- sub("^HALLMARK_", "", display_description)
  display_description <- gsub("_", " ", display_description, fixed = TRUE)
  unwrapped_description <- display_description
  if (label_wrap_width > 0) {
    display_description <- vapply(
      display_description,
      function(label) {
        paste(strwrap(label, width = label_wrap_width), collapse = "\n")
      },
      character(1)
    )
  }
  plot_data$Description_original <- original_description
  plot_data$Description <- factor(
    display_description,
    levels = rev(unique(display_description))
  )

  plot_data$abs_NES <- abs(plot_data$NES)
  # Zero nominal p-values are clamped only to keep the display transform finite.
  safe_pvalue <- pmax(plot_data$pvalue, .Machine$double.xmin)
  plot_data$neg_log10_pvalue <- -log10(safe_pvalue)

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = NES,
      y = Description,
      size = abs_NES,
      color = neg_log10_pvalue
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    ggplot2::scale_size_continuous() +
    ggplot2::scale_color_gradient(low = color_low, high = color_high) +
    ggplot2::labs(
      x = "Normalized enrichment score (NES)",
      y = NULL,
      size = "|NES|",
      color = "-log10(p-value)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

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
    dpi = dpi,
    bg = "white"
  )
  ggplot2::ggsave(
    filename = paste0(output_prefix, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )

  plot_data_for_tsv <- plot_data
  plot_data_for_tsv$Description <- unwrapped_description
  output_columns <- intersect(
    c(
      "ID", "Description", "Description_original", "NES", "pvalue",
      "abs_NES", "neg_log10_pvalue", "p.adjust", "qvalue", "setSize"
    ),
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
