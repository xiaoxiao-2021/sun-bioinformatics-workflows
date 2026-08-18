args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("07_GSEA_visualization.R requires one config path")
for (package in c("yaml", "ggplot2")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing R package: ", package)
}

config <- yaml::read_yaml(normalizePath(args[[1]], mustWork = TRUE))
if (!is.logical(config$gsea_draw_curves) || length(config$gsea_draw_curves) != 1L || is.na(config$gsea_draw_curves)) {
  stop("gsea_draw_curves must be true or false")
}
if (isTRUE(config$gsea_draw_curves) && !requireNamespace("enrichplot", quietly = TRUE)) {
  stop("Missing R package required for GSEA curves: enrichplot")
}
project_dir <- normalizePath(config$project_dir, mustWork = TRUE)
result_dir <- file.path(project_dir, "results", "GSEA")
figure_dir <- file.path(project_dir, "figures", "GSEA")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

balanced_top <- function(result, metric, n) {
  result <- result[order(result[[metric]], -abs(result$NES)), , drop = FALSE]
  if (nrow(result) <= n) return(result)
  positive <- result[result$NES > 0, , drop = FALSE]
  negative <- result[result$NES < 0, , drop = FALSE]
  half <- floor(n / 2)
  selected <- rbind(head(positive, half), head(negative, half))
  remaining <- result[!result$ID %in% selected$ID, , drop = FALSE]
  if (nrow(selected) < n) selected <- rbind(selected, head(remaining, n - nrow(selected)))
  selected
}

for (database in c("Hallmark", "GO_BP", "KEGG")) {
  stem <- paste0(config$dataset_id, "_GSEA_", database)
  table_file <- file.path(result_dir, paste0(stem, "_all.tsv"))
  object_file <- file.path(result_dir, paste0(stem, "_object.rds"))
  if (!file.exists(table_file)) {
    cat("GSEA plots skipped; table absent:", table_file, "\n")
    next
  }
  result <- utils::read.delim(table_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(result)) {
    cat("GSEA plots skipped; no pathways:", database, "\n")
    next
  }
  result$NES <- suppressWarnings(as.numeric(result$NES))
  result$setSize <- suppressWarnings(as.numeric(result$setSize))
  result$pvalue <- suppressWarnings(as.numeric(result$pvalue))
  result$p.adjust <- suppressWarnings(as.numeric(result$p.adjust))

  evidence_tables <- list(
    formal = result[is.finite(result$p.adjust) & result$p.adjust < config$gsea_FDR_cutoff, , drop = FALSE],
    exploratory = result[is.finite(result$pvalue) & result$pvalue < config$gsea_pvalue_cutoff, , drop = FALSE]
  )
  for (evidence in names(evidence_tables)) {
    metric <- if (evidence == "formal") "p.adjust" else "pvalue"
    cutoff <- if (evidence == "formal") config$gsea_FDR_cutoff else config$gsea_pvalue_cutoff
    selected <- evidence_tables[[evidence]]
    if (!nrow(selected)) {
      cat("GSEA", evidence, "plots skipped; no pathways:", database, "\n")
      next
    }
    selected <- balanced_top(selected, metric, as.integer(config$gsea_show_category_n))
    labels <- make.unique(selected$Description)
    selected$plot_description <- factor(labels, levels = rev(labels))
    selected$minus_log10 <- -log10(pmax(selected[[metric]], .Machine$double.xmin))
    evidence_title <- if (evidence == "formal") {
      paste0("formal adjusted evidence (BH FDR < ", cutoff, ")")
    } else {
      paste0("exploratory nominal evidence (P < ", cutoff, ")")
    }
    title <- paste(config$dataset_id, database, evidence_title, sep = " | ")
    dot <- ggplot2::ggplot(selected, ggplot2::aes(x = NES, y = plot_description, size = setSize, color = minus_log10)) +
      ggplot2::geom_vline(xintercept = 0, color = "grey60") + ggplot2::geom_point() +
      ggplot2::scale_color_gradient(low = "#4DBBD5", high = "#E64B35") +
      ggplot2::labs(title = title, x = paste0("NES (+ ", config$case_group, "; - ", config$control_group, ")"), y = NULL, size = "Set size", color = paste0("-log10(", metric, ")")) +
      ggplot2::theme_bw(base_size = 10) + ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
    lollipop <- ggplot2::ggplot(selected, ggplot2::aes(y = plot_description)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = plot_description), color = "grey65") +
      ggplot2::geom_point(ggplot2::aes(x = NES, size = setSize, color = minus_log10)) +
      ggplot2::geom_vline(xintercept = 0, color = "grey45") +
      ggplot2::scale_color_gradient(low = "#4DBBD5", high = "#E64B35") +
      ggplot2::labs(title = title, x = paste0("NES (+ ", config$case_group, "; - ", config$control_group, ")"), y = NULL, size = "Set size", color = paste0("-log10(", metric, ")")) +
      ggplot2::theme_bw(base_size = 10) + ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
    for (plot_type in c("NES_dotplot", "NES_lollipop")) {
      plot <- if (plot_type == "NES_dotplot") dot else lollipop
      output_stem <- file.path(figure_dir, paste0(stem, "_", evidence, "_", plot_type))
      ggplot2::ggsave(paste0(output_stem, ".png"), plot, width = 9, height = 6, dpi = 300)
      ggplot2::ggsave(paste0(output_stem, ".pdf"), plot, width = 9, height = 6)
    }
  }

  if (!isTRUE(config$gsea_draw_curves)) {
    cat("GSEA enrichment curves disabled by config for:", database, "\n")
    next
  }
  object <- if (file.exists(object_file)) readRDS(object_file) else NULL
  if (is.null(object)) {
    cat("GSEA curves skipped; enrichment object absent:", database, "\n")
    next
  }
  formal <- evidence_tables$formal
  exploratory <- evidence_tables$exploratory
  use_formal <- nrow(formal) >= 2L || !nrow(exploratory)
  curve_table <- if (use_formal) formal else exploratory
  curve_evidence <- if (use_formal) "formal_FDR" else "exploratory_nominal_P"
  curve_metric <- if (use_formal) "p.adjust" else "pvalue"
  if (!nrow(curve_table)) {
    cat("GSEA curves skipped; neither formal nor exploratory pathways available:", database, "\n")
    next
  }
  curve_table <- curve_table[order(curve_table[[curve_metric]], -abs(curve_table$NES)), , drop = FALSE]
  curve_table <- rbind(
    head(curve_table[curve_table$NES > 0, , drop = FALSE], as.integer(config$gsea_curve_n)),
    head(curve_table[curve_table$NES < 0, , drop = FALSE], as.integer(config$gsea_curve_n))
  )
  for (i in seq_len(nrow(curve_table))) {
    pathway_id <- curve_table$ID[[i]]
    safe_id <- gsub("[^A-Za-z0-9._-]+", "_", pathway_id)
    curve_title <- paste(
      database, curve_table$Description[[i]], curve_evidence,
      paste0("NES=", signif(curve_table$NES[[i]], 3)), sep = " | "
    )
    curve <- tryCatch(
      enrichplot::gseaplot2(object, geneSetID = pathway_id, title = curve_title, pvalue_table = TRUE),
      error = function(error) {
        warning("Curve skipped for ", pathway_id, ": ", conditionMessage(error))
        NULL
      }
    )
    if (is.null(curve)) next
    output_stem <- file.path(figure_dir, paste0(stem, "_curve_", curve_evidence, "_", safe_id))
    ggplot2::ggsave(paste0(output_stem, ".png"), curve, width = 9, height = 6, dpi = 300)
    ggplot2::ggsave(paste0(output_stem, ".pdf"), curve, width = 9, height = 6)
  }
}
