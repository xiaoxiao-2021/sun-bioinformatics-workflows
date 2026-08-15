args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Expected one config.yml argument")
cfg <- yaml::read_yaml(normalizePath(args[1], mustWork = TRUE))
# KEGG REST can be slow; keep the same ORA method but allow a longer download.
options(timeout = max(300, getOption("timeout")))
project_dir <- normalizePath(cfg$project_dir, mustWork = TRUE)
de_dir <- file.path(project_dir, "results", cfg$dataset_id, "DE")
result_dir <- file.path(project_dir, "results", cfg$dataset_id, "enrichment")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- paste(cfg$dataset_id, cfg$case_group, "vs", cfg$control_group, sep = "_")
lfc_tag <- format(cfg$TREAT_lfc_cutoff, trim = TRUE, scientific = FALSE)
fdr_tag <- format(cfg$TREAT_FDR_cutoff, trim = TRUE, scientific = FALSE)
deg_fdr_tag <- format(cfg$deg_FDR_cutoff, trim = TRUE, scientific = FALSE)
deg_lfc_tag <- format(cfg$deg_logFC_cutoff, trim = TRUE, scientific = FALSE)
pvalue_tag <- format(cfg$pvalue_cutoff, trim = TRUE, scientific = FALSE)
enrichment_p_tag <- format(cfg$enrichment_pvalue_cutoff, trim = TRUE, scientific = FALSE)
enrichment_fdr_tag <- format(cfg$enrichment_FDR_cutoff, trim = TRUE, scientific = FALSE)
enrichment_gene_logFC_cutoff <- cfg$enrichment_gene_logFC_cutoff
if (is.null(enrichment_gene_logFC_cutoff)) {
  enrichment_gene_logFC_cutoff <- cfg$deg_logFC_cutoff
}
if (is.null(enrichment_gene_logFC_cutoff)) {
  stop("Missing both enrichment_gene_logFC_cutoff and fallback deg_logFC_cutoff")
}
if (length(enrichment_gene_logFC_cutoff) != 1 ||
    !is.finite(enrichment_gene_logFC_cutoff) ||
    enrichment_gene_logFC_cutoff < 0) {
  stop("enrichment_gene_logFC_cutoff must be one finite, non-negative number")
}
enrichment_gene_lfc_tag <- format(
  enrichment_gene_logFC_cutoff, trim = TRUE, scientific = FALSE
)

# Background: all standard limma-tested genes with mapped ENTREZID
background_data <- read.delim(
  file.path(de_dir, paste0(prefix, "_limma_all_genes_annotated.tsv")),
  check.names = FALSE, stringsAsFactors = FALSE
)
background_valid_entrez <- !is.na(background_data$ENTREZID) &
  background_data$ENTREZID != ""
background_entrez <- unique(as.character(background_data$ENTREZID[background_valid_entrez]))

# Build ORA foreground independently from visualization/DEG-table effect-size filters.
if (cfg$downstream_deg_method == "limma") {
  foreground_candidates <- background_data[
    !is.na(background_data$adj.P.Val) &
      background_data$adj.P.Val < cfg$deg_FDR_cutoff &
      is.finite(background_data$logFC) &
      abs(background_data$logFC) >= enrichment_gene_logFC_cutoff &
      background_data$logFC != 0,
    , drop = FALSE
  ]
  foreground_significance_label <- paste("Gene FDR cutoff:", cfg$deg_FDR_cutoff)
  ora_base <- paste0(
    prefix, "_ORA_geneFDR", deg_fdr_tag,
    "_logFC", enrichment_gene_lfc_tag
  )
} else if (cfg$downstream_deg_method == "limma_pvalue") {
  foreground_candidates <- background_data[
    !is.na(background_data$P.Value) &
      background_data$P.Value < cfg$pvalue_cutoff &
      is.finite(background_data$logFC) &
      abs(background_data$logFC) >= enrichment_gene_logFC_cutoff &
      background_data$logFC != 0,
    , drop = FALSE
  ]
  foreground_significance_label <- paste("Gene P-value cutoff:", cfg$pvalue_cutoff)
  ora_base <- paste0(
    prefix, "_ORA_geneP", pvalue_tag,
    "_logFC", enrichment_gene_lfc_tag
  )
} else {
  treat_file <- file.path(
    de_dir,
    paste0(prefix, "_TREAT_lfc", lfc_tag, "_FDR", fdr_tag, "_DEG_annotated.tsv")
  )
  foreground_candidates <- read.delim(
    treat_file, check.names = FALSE, stringsAsFactors = FALSE
  )
  foreground_candidates <- foreground_candidates[
    is.finite(foreground_candidates$logFC) & foreground_candidates$logFC != 0,
    , drop = FALSE
  ]
  foreground_significance_label <- paste0(
    "TREAT hypothesis: |logFC| > ", cfg$TREAT_lfc_cutoff,
    "; FDR < ", cfg$TREAT_FDR_cutoff
  )
  ora_base <- paste0(
    prefix, "_ORA_TREAT_lfc", lfc_tag, "_FDR", fdr_tag
  )
}

# logFC defines foreground eligibility and direction only; enrichGO/enrichKEGG
# receive gene IDs, the unchanged universe, and pathway membership—not logFC.
foreground_candidates$direction <- ifelse(
  foreground_candidates$logFC > 0, "UP", "DOWN"
)
foreground_valid_entrez <- !is.na(foreground_candidates$ENTREZID) &
  foreground_candidates$ENTREZID != ""
foreground_mapped <- foreground_candidates[foreground_valid_entrez, , drop = FALSE]
foreground_mapped$ENTREZID <- as.character(foreground_mapped$ENTREZID)
# ORA consumes unique ENTREZID values; retain one annotated row per actual input ID.
foreground_mapped <- foreground_mapped[
  !duplicated(foreground_mapped$ENTREZID), , drop = FALSE
]
up_entrez <- foreground_mapped$ENTREZID[foreground_mapped$direction == "UP"]
down_entrez <- foreground_mapped$ENTREZID[foreground_mapped$direction == "DOWN"]

foreground_columns <- intersect(
  c(
    "ENSEMBL", "SYMBOL", "ENTREZID", "DISPLAY_NAME", "GENE_BIOTYPE",
    "logFC", "P.Value", "adj.P.Val", "direction"
  ),
  names(foreground_mapped)
)
foreground_file <- file.path(
  result_dir, paste0(ora_base, "_foreground_genes.tsv")
)
write.table(
  foreground_mapped[, foreground_columns, drop = FALSE], foreground_file,
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

cat("ORA background definition:\n")
cat("all limma-tested genes with valid ENTREZID\n")
cat("Background genes:", length(background_entrez), "\n\n")
cat("============================================================\n")
cat("ORA foreground summary\n")
cat("============================================================\n")
cat("Method:", cfg$downstream_deg_method, "\n")
cat(foreground_significance_label, "\n")
if (cfg$downstream_deg_method == "treat") {
  cat("Gene |logFC| cutoff: not additionally applied for TREAT\n")
} else {
  cat("Gene |logFC| cutoff:", enrichment_gene_logFC_cutoff, "\n")
  if (enrichment_gene_logFC_cutoff == 0) {
    cat("Gene |logFC| cutoff = 0 -> no effect-size filtering\n")
  }
}
cat("\nBackground genes with valid ENTREZID:", length(background_entrez), "\n\n")
cat("Foreground before ENTREZID mapping:\n")
cat("Total:", nrow(foreground_candidates), "\n")
cat("UP:", sum(foreground_candidates$direction == "UP"), "\n")
cat("DOWN:", sum(foreground_candidates$direction == "DOWN"), "\n\n")
cat("Foreground with valid ENTREZID:\n")
cat("Total:", nrow(foreground_mapped), "\n")
cat("UP:", length(up_entrez), "\n")
cat("DOWN:", length(down_entrez), "\n")
cat("Foreground audit table:", foreground_file, "\n")
if (!length(up_entrez)) cat("No UP genes available for ORA. GO/KEGG UP skipped.\n")
if (!length(down_entrez)) cat("No DOWN genes available for ORA. GO/KEGG DOWN skipped.\n")

# Run the four requested ORA analyses with identical output handling
run_ora <- function(genes, database, direction) {
  analysis_name <- paste(database, direction, sep = "_")
  output_base <- paste(ora_base, analysis_name, sep = "_")
  empty_result <- data.frame(
    ID = character(), Description = character(), GeneRatio = character(),
    Count = integer(), pvalue = numeric(), p.adjust = numeric(),
    stringsAsFactors = FALSE
  )
  if (!length(genes)) {
    all_result <- empty_result
  } else if (database == "GO_BP") {
    enrichment <- tryCatch(
      clusterProfiler::enrichGO(
        gene = genes, universe = background_entrez, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1,
        minGSSize = cfg$min_GS_size, maxGSSize = cfg$max_GS_size,
        readable = TRUE
      ),
      error = function(e) {
        cat("GO-BP ORA skipped for", direction, ":", conditionMessage(e), "\n")
        NULL
      }
    )
    all_result <- if (is.null(enrichment)) empty_result else as.data.frame(enrichment)
  } else {
    enrichment <- tryCatch(
      clusterProfiler::enrichKEGG(
        gene = genes, universe = background_entrez, organism = "hsa",
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
        minGSSize = cfg$min_GS_size, maxGSSize = cfg$max_GS_size
      ),
      error = function(e) {
        cat("KEGG ORA skipped for", direction, ":", conditionMessage(e), "\n")
        NULL
      }
    )
    all_result <- if (is.null(enrichment)) empty_result else as.data.frame(enrichment)
  }
  # Exploratory nominal pathway result
  nominal_result <- if (nrow(all_result)) {
    all_result[
      !is.na(all_result$pvalue) &
        all_result$pvalue < cfg$enrichment_pvalue_cutoff,
      , drop = FALSE
    ]
  } else {
    all_result
  }

  # Formal FDR-adjusted pathway result
  significant_result <- if (nrow(all_result)) {
    all_result[
      !is.na(all_result$p.adjust) &
        all_result$p.adjust < cfg$enrichment_FDR_cutoff,
      , drop = FALSE
    ]
  } else {
    all_result
  }

  write.table(all_result, file.path(result_dir, paste0(output_base, "_all.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  # Save nominal enrichment table
  write.table(nominal_result, file.path(result_dir, paste0(output_base, "_P", enrichment_p_tag, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  write.table(significant_result, file.path(result_dir, paste0(output_base, "_FDR", enrichment_fdr_tag, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  cat(analysis_name, "\n")
  cat("All terms:", nrow(all_result), "\n")
  cat("Nominal P <", cfg$enrichment_pvalue_cutoff, ":", nrow(nominal_result), "\n")
  cat("FDR <", cfg$enrichment_FDR_cutoff, ":", nrow(significant_result), "\n")
  if (!nrow(nominal_result)) {
    cat("No nominally significant enrichment result for", analysis_name, ".\n")
  }
  if (!nrow(significant_result)) {
    cat("No FDR-significant enrichment result for", analysis_name, ".\n")
  }
}

run_ora(up_entrez, "GO_BP", "UP")
run_ora(down_entrez, "GO_BP", "DOWN")
run_ora(up_entrez, "KEGG", "UP")
run_ora(down_entrez, "KEGG", "DOWN")
