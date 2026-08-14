args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Expected one config.yml argument")
cfg <- yaml::read_yaml(normalizePath(args[1], mustWork = TRUE))

required_packages <- c("clusterProfiler", "org.Hs.eg.db")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop("Missing required GSEA package(s): ", paste(missing_packages, collapse = ", "))
}

# KEGG REST can be slow; use the same timeout policy as the ORA module.
options(timeout = max(300, getOption("timeout")))
project_dir <- normalizePath(cfg$project_dir, mustWork = TRUE)
de_dir <- file.path(project_dir, "results", cfg$dataset_id, "DE")
result_dir <- file.path(project_dir, "results", cfg$dataset_id, "GSEA")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- paste(cfg$dataset_id, cfg$case_group, "vs", cfg$control_group, sep = "_")
gsea_p_tag <- format(cfg$gsea_pvalue_cutoff, trim = TRUE, scientific = FALSE)
gsea_fdr_tag <- format(cfg$gsea_FDR_cutoff, trim = TRUE, scientific = FALSE)

# GSEA is always ranked from the standard limma all-gene result.
input_file <- file.path(de_dir, paste0(prefix, "_limma_all_genes_annotated.tsv"))
all_genes <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c(
  "ENSEMBL", "ENTREZID", "SYMBOL", "DISPLAY_NAME",
  "logFC", "t", "P.Value", "adj.P.Val"
)
if (!all(required_columns %in% names(all_genes))) {
  stop("Missing GSEA input column(s): ",
       paste(setdiff(required_columns, names(all_genes)), collapse = ", "))
}

all_genes$ENTREZID <- trimws(as.character(all_genes$ENTREZID))
all_genes$t <- suppressWarnings(as.numeric(all_genes$t))
all_genes$P.Value <- suppressWarnings(as.numeric(all_genes$P.Value))
valid_entrez <- !is.na(all_genes$ENTREZID) & all_genes$ENTREZID != ""
finite_t <- is.finite(all_genes$t)
ranking_data <- all_genes[valid_entrez & finite_t, , drop = FALSE]
ranking_data$.input_order <- seq_len(nrow(ranking_data))
ranking_data$.abs_t <- abs(ranking_data$t)
ranking_data$.p_tie <- ifelse(is.finite(ranking_data$P.Value), ranking_data$P.Value, Inf)

# One score per ENTREZID: largest |t|, then smaller P.Value, then first input row.
ranking_data <- ranking_data[
  order(
    ranking_data$ENTREZID, -ranking_data$.abs_t,
    ranking_data$.p_tie, ranking_data$.input_order
  ),
  , drop = FALSE
]
eligible_before_dedup <- nrow(ranking_data)
ranking_data <- ranking_data[!duplicated(ranking_data$ENTREZID), , drop = FALSE]
duplicate_removed <- eligible_before_dedup - nrow(ranking_data)
ranking_data <- ranking_data[
  order(-ranking_data$t, ranking_data$ENTREZID),
  , drop = FALSE
]

ranked_columns <- intersect(
  c(
    "ENTREZID", "ENSEMBL", "SYMBOL", "DISPLAY_NAME", "GENE_BIOTYPE",
    "logFC", "t", "P.Value", "adj.P.Val"
  ),
  names(ranking_data)
)
ranked_table <- ranking_data[, ranked_columns, drop = FALSE]
gene_list <- stats::setNames(ranking_data$t, ranking_data$ENTREZID)

if (!length(gene_list)) stop("GSEA ranking is empty after ENTREZID/t filtering")
if (anyDuplicated(names(gene_list))) stop("GSEA ranking contains duplicated ENTREZID")
if (any(!is.finite(gene_list))) stop("GSEA ranking contains non-finite t statistics")
if (length(gene_list) > 1 && any(diff(gene_list) > 0)) {
  stop("GSEA ranking is not sorted from largest positive t to largest negative t")
}

ranked_file <- file.path(result_dir, paste0(prefix, "_GSEA_ranked_genes.tsv"))
write.table(
  ranked_table, ranked_file, sep = "\t", quote = FALSE,
  row.names = FALSE, na = ""
)

cat("GSEA input summary\n")
cat("--------------------------------\n")
cat(sprintf("%-32s %d\n", "All limma-tested genes:", nrow(all_genes)))
cat(sprintf("%-32s %d\n", "Genes with valid ENTREZID:", sum(valid_entrez)))
cat(sprintf("%-32s %d\n", "Genes with finite t:", sum(finite_t)))
cat(sprintf("%-32s %d\n", "Valid ENTREZID + finite t:", eligible_before_dedup))
cat(sprintf("%-32s %d\n", "Duplicate ENTREZID removed:", duplicate_removed))
cat(sprintf("%-32s %d\n", "Final ranked genes:", length(gene_list)))
cat("\nRanking:\n")
cat("case group:", cfg$case_group, "\n")
cat("control group:", cfg$control_group, "\n\n")
cat("GSEA ranking direction:\n")
cat("positive t / NES ->", cfg$case_group, "\n")
cat("negative t / NES ->", cfg$control_group, "\n")
cat("Ranked gene table:", ranked_file, "\n")

adcy7 <- ranked_table[
  (!is.na(ranked_table$SYMBOL) & ranked_table$SYMBOL == "ADCY7") |
    (!is.na(ranked_table$DISPLAY_NAME) & ranked_table$DISPLAY_NAME == "ADCY7"),
  , drop = FALSE
]
if (nrow(adcy7)) {
  cat("ADCY7 ranked t:", adcy7$t[1], "\n")
  if (cfg$dataset_id == "LM3" && adcy7$t[1] <= 0) {
    stop("LM3 ADCY7 t is not positive; check comparison orientation before GSEA")
  }
} else {
  cat("ADCY7 is not present in the final ranked list.\n")
}

empty_gsea_result <- function() {
  data.frame(
    ID = character(), Description = character(), setSize = integer(),
    enrichmentScore = numeric(), NES = numeric(), pvalue = numeric(),
    p.adjust = numeric(), qvalue = numeric(), rank = integer(),
    leading_edge = character(), core_enrichment = character(),
    stringsAsFactors = FALSE
  )
}

run_gsea <- function(database) {
  database_label <- if (database == "GO_BP") "GO Biological Process" else "KEGG"
  set.seed(123)
  gsea_object <- tryCatch(
    {
      if (database == "GO_BP") {
        clusterProfiler::gseGO(
          geneList = gene_list, ont = "BP",
          OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
          exponent = 1, minGSSize = cfg$gsea_min_GS_size,
          maxGSSize = cfg$gsea_max_GS_size, eps = 0,
          pvalueCutoff = 1, pAdjustMethod = "BH",
          verbose = FALSE, seed = TRUE
        )
      } else {
        clusterProfiler::gseKEGG(
          geneList = gene_list, organism = "hsa", keyType = "kegg",
          exponent = 1, minGSSize = cfg$gsea_min_GS_size,
          maxGSSize = cfg$gsea_max_GS_size, eps = 0,
          pvalueCutoff = 1, pAdjustMethod = "BH",
          verbose = FALSE, seed = TRUE
        )
      }
    },
    error = function(e) {
      cat(database_label, "GSEA skipped:", conditionMessage(e), "\n")
      NULL
    }
  )

  all_result <- if (is.null(gsea_object)) {
    empty_gsea_result()
  } else {
    result <- as.data.frame(gsea_object)
    if (nrow(result)) result else empty_gsea_result()
  }
  nominal_result <- if (nrow(all_result)) {
    all_result[
      !is.na(all_result$pvalue) & all_result$pvalue < cfg$gsea_pvalue_cutoff,
      , drop = FALSE
    ]
  } else {
    all_result
  }
  formal_result <- if (nrow(all_result)) {
    all_result[
      !is.na(all_result$p.adjust) & all_result$p.adjust < cfg$gsea_FDR_cutoff,
      , drop = FALSE
    ]
  } else {
    all_result
  }

  output_base <- paste0(prefix, "_GSEA_", database)
  write.table(
    all_result, file.path(result_dir, paste0(output_base, "_all.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
  write.table(
    nominal_result,
    file.path(result_dir, paste0(output_base, "_P", gsea_p_tag, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
  write.table(
    formal_result,
    file.path(result_dir, paste0(output_base, "_FDR", gsea_fdr_tag, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
  saveRDS(gsea_object, file.path(result_dir, paste0(output_base, "_object.rds")))

  cat(database_label, "GSEA\n")
  cat("All pathways:", nrow(all_result), "\n")
  cat("Exploratory nominal P <", cfg$gsea_pvalue_cutoff, ":", nrow(nominal_result), "\n")
  cat("Formal FDR <", cfg$gsea_FDR_cutoff, ":", nrow(formal_result), "\n")
  if (!nrow(all_result)) cat("No", database_label, "GSEA result returned.\n")
  if (!nrow(nominal_result)) cat("No nominally significant", database_label, "GSEA result.\n")
  if (!nrow(formal_result)) cat("No FDR-significant", database_label, "GSEA result.\n")

  list(object = gsea_object, result = all_result)
}

make_leading_edge <- function(database, all_result) {
  output_file <- file.path(
    result_dir, paste0(prefix, "_GSEA_", database, "_leading_edge_genes.tsv")
  )
  annotation_columns <- intersect(
    c(
      "ENSEMBL", "SYMBOL", "DISPLAY_NAME", "GENE_BIOTYPE",
      "logFC", "t", "P.Value", "adj.P.Val"
    ),
    names(ranked_table)
  )
  empty_leading_edge <- data.frame(
    Pathway_ID = character(), Pathway = character(), NES = numeric(),
    pathway_pvalue = numeric(), pathway_p_adjust = numeric(),
    ENTREZID = character(), stringsAsFactors = FALSE
  )
  for (column in annotation_columns) empty_leading_edge[[column]] <- character()

  if (!nrow(all_result) || !"core_enrichment" %in% names(all_result)) {
    leading_edge <- empty_leading_edge
  } else {
    pieces <- lapply(seq_len(nrow(all_result)), function(i) {
      entrez <- strsplit(as.character(all_result$core_enrichment[i]), "/", fixed = TRUE)[[1]]
      entrez <- entrez[!is.na(entrez) & entrez != ""]
      if (!length(entrez)) return(NULL)
      data.frame(
        Pathway_ID = rep(as.character(all_result$ID[i]), length(entrez)),
        Pathway = rep(as.character(all_result$Description[i]), length(entrez)),
        NES = rep(as.numeric(all_result$NES[i]), length(entrez)),
        pathway_pvalue = rep(as.numeric(all_result$pvalue[i]), length(entrez)),
        pathway_p_adjust = rep(as.numeric(all_result$p.adjust[i]), length(entrez)),
        ENTREZID = entrez,
        stringsAsFactors = FALSE
      )
    })
    pieces <- Filter(Negate(is.null), pieces)
    if (!length(pieces)) {
      leading_edge <- empty_leading_edge
    } else {
      leading_edge <- do.call(rbind, pieces)
      match_index <- match(leading_edge$ENTREZID, ranked_table$ENTREZID)
      for (column in annotation_columns) {
        leading_edge[[column]] <- ranked_table[[column]][match_index]
      }
    }
  }
  write.table(
    leading_edge, output_file, sep = "\t", quote = FALSE,
    row.names = FALSE, na = ""
  )
  cat(database, "leading-edge rows:", nrow(leading_edge), "\n")
  cat("Leading-edge table:", output_file, "\n")
}

go_gsea <- run_gsea("GO_BP")
make_leading_edge("GO_BP", go_gsea$result)
kegg_gsea <- run_gsea("KEGG")
make_leading_edge("KEGG", kegg_gsea$result)
