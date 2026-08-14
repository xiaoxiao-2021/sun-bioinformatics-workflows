args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Expected one config.yml argument")
cfg <- yaml::read_yaml(normalizePath(args[1], mustWork = TRUE))
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Missing required GSEA visualization package: ggplot2")
}

project_dir <- normalizePath(cfg$project_dir, mustWork = TRUE)
result_dir <- file.path(project_dir, "results", cfg$dataset_id, "GSEA")
figure_dir <- file.path(project_dir, "figures", cfg$dataset_id, "GSEA")
plot_data_dir <- file.path(result_dir, "plot_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_data_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- paste(cfg$dataset_id, cfg$case_group, "vs", cfg$control_group, sep = "_")
gsea_p_tag <- format(cfg$gsea_pvalue_cutoff, trim = TRUE, scientific = FALSE)
gsea_fdr_tag <- format(cfg$gsea_FDR_cutoff, trim = TRUE, scientific = FALSE)

read_gsea_table <- function(database, evidence) {
  result_tag <- if (evidence == "FDR") paste0("FDR", gsea_fdr_tag) else paste0("P", gsea_p_tag)
  input_file <- file.path(
    result_dir, paste0(prefix, "_GSEA_", database, "_", result_tag, ".tsv")
  )
  if (!file.exists(input_file) || file.info(input_file)$size == 0) return(data.frame())
  read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
}

# Keep both NES directions when possible, then fill unused slots deterministically.
select_balanced <- function(data, n, metric_column) {
  required <- c("ID", "Description", "setSize", "NES", metric_column)
  if (!nrow(data) || !all(required %in% names(data))) return(data[0, , drop = FALSE])
  data$NES <- suppressWarnings(as.numeric(data$NES))
  data[[metric_column]] <- suppressWarnings(as.numeric(data[[metric_column]]))
  data <- data[
    is.finite(data$NES) & is.finite(data[[metric_column]]),
    , drop = FALSE
  ]
  if (!nrow(data)) return(data)
  order_rows <- function(x) {
    x[order(x[[metric_column]], -abs(x$NES), x$ID), , drop = FALSE]
  }
  positive <- order_rows(data[data$NES > 0, , drop = FALSE])
  negative <- order_rows(data[data$NES < 0, , drop = FALSE])
  if (!nrow(positive) || !nrow(negative)) return(head(order_rows(data), n))

  positive_n <- min(nrow(positive), ceiling(n / 2))
  negative_n <- min(nrow(negative), floor(n / 2))
  selected <- rbind(head(positive, positive_n), head(negative, negative_n))
  remaining_n <- n - nrow(selected)
  if (remaining_n > 0) {
    remaining <- rbind(
      positive[-seq_len(positive_n), , drop = FALSE],
      negative[-seq_len(negative_n), , drop = FALSE]
    )
    if (nrow(remaining)) selected <- rbind(selected, head(order_rows(remaining), remaining_n))
  }
  selected[!duplicated(selected$ID), , drop = FALSE]
}

wrap_pathway <- function(x, width = 48) {
  vapply(x, function(label) paste(strwrap(label, width = width), collapse = "\n"), character(1))
}

save_plot_pair <- function(plot, output_base, width, height) {
  ggplot2::ggsave(paste0(output_base, ".png"), plot, width = width, height = height, dpi = 300)
  ggplot2::ggsave(paste0(output_base, ".pdf"), plot, width = width, height = height)
}

make_summary_plot <- function(database, data, evidence, plot_type) {
  if (!nrow(data)) return(invisible(FALSE))
  database_label <- if (database == "GO_BP") "GO Biological Process" else "KEGG"
  metric_column <- if (evidence == "FDR") "p.adjust" else "pvalue"
  metric_label <- if (evidence == "FDR") "-log10(FDR)" else "-log10(P.Value)"
  evidence_label <- if (evidence == "FDR") {
    paste("Formal: FDR <", cfg$gsea_FDR_cutoff)
  } else {
    paste("Exploratory: nominal P <", cfg$gsea_pvalue_cutoff)
  }
  file_tag <- if (evidence == "FDR") "FDR" else "P_exploratory"
  plot_data <- select_balanced(data, cfg$gsea_show_category_n, metric_column)
  if (!nrow(plot_data)) return(invisible(FALSE))
  plot_data <- plot_data[order(plot_data$NES), , drop = FALSE]
  plot_data$plot_label <- make.unique(wrap_pathway(plot_data$Description))
  plot_data$plot_label <- factor(plot_data$plot_label, levels = plot_data$plot_label)
  plot_data$minus_log10_metric <- -log10(
    pmax(as.numeric(plot_data[[metric_column]]), .Machine$double.xmin)
  )
  plot_data$enriched_group <- ifelse(plot_data$NES > 0, cfg$case_group, cfg$control_group)
  direction_subtitle <- paste0(
    evidence_label, "\nNES > 0 -> ", cfg$case_group,
    "; NES < 0 -> ", cfg$control_group
  )
  title <- paste(
    cfg$dataset_id, cfg$case_group, "vs", cfg$control_group,
    "\nGSEA", database_label
  )

  if (plot_type == "dotplot") {
    plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(NES, plot_label, size = setSize, color = minus_log10_metric)
    ) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
      ggplot2::geom_point(alpha = .9) +
      ggplot2::scale_color_gradient(low = "#4DBBD5", high = "#E64B35") +
      ggplot2::labs(
        title = title, subtitle = direction_subtitle,
        x = "Normalized enrichment score (NES)", y = NULL,
        size = "Gene set size", color = metric_label
      ) +
      ggplot2::theme_classic(base_size = 12)
    suffix <- paste0(file_tag, "_dotplot")
  } else {
    plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(NES, plot_label, color = enriched_group)
    ) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = NES, y = plot_label, yend = plot_label),
        linewidth = .8
      ) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
      ggplot2::scale_color_manual(
        values = stats::setNames(c("#E64B35", "#4DBBD5"), c(cfg$case_group, cfg$control_group))
      ) +
      ggplot2::labs(
        title = title, subtitle = direction_subtitle,
        x = "Normalized enrichment score (NES)", y = NULL,
        color = "Enriched group"
      ) +
      ggplot2::theme_classic(base_size = 12)
    suffix <- paste0(file_tag, "_NES_lollipop")
  }

  output_base <- file.path(
    figure_dir, paste0(prefix, "_GSEA_", database, "_", suffix)
  )
  plot_data_file <- file.path(
    plot_data_dir, paste0(prefix, "_GSEA_", database, "_", suffix, "_pathways.tsv")
  )
  write.table(
    plot_data[, setdiff(names(plot_data), "plot_label"), drop = FALSE],
    plot_data_file, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
  plot_height <- max(4.5, 2.7 + .38 * nrow(plot_data))
  save_plot_pair(plot, output_base, width = 10.5, height = plot_height)
  cat(database, evidence, plot_type, ": plotted", nrow(plot_data), "pathways\n")
  invisible(TRUE)
}

select_curves <- function(data, metric_column) {
  if (!nrow(data)) return(data)
  data$NES <- suppressWarnings(as.numeric(data$NES))
  data[[metric_column]] <- suppressWarnings(as.numeric(data[[metric_column]]))
  data <- data[is.finite(data$NES) & is.finite(data[[metric_column]]), , drop = FALSE]
  order_rows <- function(x) {
    x[order(x[[metric_column]], -abs(x$NES), x$ID), , drop = FALSE]
  }
  positive <- head(order_rows(data[data$NES > 0, , drop = FALSE]), cfg$gsea_curve_n)
  negative <- head(order_rows(data[data$NES < 0, , drop = FALSE]), cfg$gsea_curve_n)
  rbind(positive, negative)
}

make_curves <- function(database, data, evidence) {
  if (!nrow(data)) {
    cat("No", database, evidence, "pathway available for GSEA curves. Curves skipped.\n")
    return(invisible(NULL))
  }
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    cat("Warning: enrichplot is unavailable;", database, "GSEA curves skipped.\n")
    return(invisible(NULL))
  }
  object_file <- file.path(result_dir, paste0(prefix, "_GSEA_", database, "_object.rds"))
  if (!file.exists(object_file)) {
    cat("Warning: GSEA object is missing for", database, "; curves skipped.\n")
    return(invisible(NULL))
  }
  gsea_object <- readRDS(object_file)
  if (is.null(gsea_object)) {
    cat("Warning: GSEA object is empty for", database, "; curves skipped.\n")
    return(invisible(NULL))
  }

  metric_column <- if (evidence == "FDR") "p.adjust" else "pvalue"
  evidence_label <- if (evidence == "FDR") {
    paste("Formal FDR <", cfg$gsea_FDR_cutoff)
  } else {
    paste("Exploratory nominal P <", cfg$gsea_pvalue_cutoff)
  }
  file_tag <- if (evidence == "FDR") "FDR" else "P_exploratory"
  curve_data <- select_curves(data, metric_column)
  curve_table <- file.path(
    plot_data_dir,
    paste0(prefix, "_GSEA_", database, "_", file_tag, "_curve_pathways.tsv")
  )
  write.table(curve_data, curve_table, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  if (!nrow(curve_data)) {
    cat("No", database, evidence, "pathway available for GSEA curves. Curves skipped.\n")
    return(invisible(NULL))
  }

  for (i in seq_len(nrow(curve_data))) {
    pathway_id <- as.character(curve_data$ID[i])
    safe_id <- gsub("[^A-Za-z0-9._-]", "_", pathway_id)
    curve_title <- paste0(
      cfg$dataset_id, " ", cfg$case_group, " vs ", cfg$control_group,
      "\n", curve_data$Description[i],
      "\n", evidence_label,
      "; NES > 0 -> ", cfg$case_group,
      "; NES < 0 -> ", cfg$control_group
    )
    output_base <- file.path(
      figure_dir,
      paste0(prefix, "_GSEA_", database, "_", file_tag, "_curve_", safe_id)
    )
    tryCatch(
      {
        curve <- enrichplot::gseaplot2(
          gsea_object, geneSetID = pathway_id,
          title = curve_title, base_size = 11
        )
        save_plot_pair(curve, output_base, width = 9.5, height = 7)
        cat(database, evidence, "curve:", pathway_id, "-", curve_data$Description[i], "\n")
      },
      error = function(e) {
        cat("Warning:", database, "curve skipped for", pathway_id, ":", conditionMessage(e), "\n")
      }
    )
  }
}

for (database in c("GO_BP", "KEGG")) {
  # Remove stale figures for this dataset/database before recreating the selected set.
  stale <- list.files(
    figure_dir,
    pattern = paste0("^", prefix, "_GSEA_", database, "_.*\\.(png|pdf)$"),
    full.names = TRUE
  )
  unlink(stale)

  formal <- read_gsea_table(database, "FDR")
  nominal <- read_gsea_table(database, "P")
  if (nrow(formal)) {
    make_summary_plot(database, formal, "FDR", "dotplot")
  } else {
    cat("No formal FDR", database, "GSEA result. Formal dotplot skipped.\n")
  }
  # With fewer than three formal pathways, also expose explicitly exploratory evidence.
  if (nrow(formal) < 3 && nrow(nominal)) {
    make_summary_plot(database, nominal, "P", "dotplot")
  } else if (nrow(formal) < 3) {
    cat("No nominal", database, "GSEA result. Exploratory dotplot skipped.\n")
  }

  if (nrow(formal)) {
    make_summary_plot(database, formal, "FDR", "lollipop")
    make_curves(database, formal, "FDR")
  } else if (nrow(nominal)) {
    make_summary_plot(database, nominal, "P", "lollipop")
    make_curves(database, nominal, "P")
  } else {
    cat("No", database, "GSEA pathway available for lollipop or curves. Plots skipped.\n")
  }
}
