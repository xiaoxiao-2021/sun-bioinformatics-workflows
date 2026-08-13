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

# Background: all standard limma-tested genes with mapped ENTREZID
background_data <- read.delim(file.path(de_dir, paste0(prefix, "_limma_all_genes_annotated.tsv")), check.names = FALSE)
background_entrez <- unique(as.character(background_data$ENTREZID[!is.na(background_data$ENTREZID) & background_data$ENTREZID != ""]))

# Select downstream DEG source
if (cfg$downstream_deg_method == "limma") {
  deg_file <- file.path(
    de_dir,
    paste0(prefix, "_limma_FDR", deg_fdr_tag, "_logFC", deg_lfc_tag, "_DEG_annotated.tsv")
  )
} else if (cfg$downstream_deg_method == "limma_pvalue") {
  deg_file <- file.path(
    de_dir,
    paste0(prefix, "_limma_P", pvalue_tag, "_logFC", deg_lfc_tag, "_DEG_annotated.tsv")
  )
} else {
  deg_file <- file.path(
    de_dir,
    paste0(prefix, "_TREAT_lfc", lfc_tag, "_FDR", fdr_tag, "_DEG_annotated.tsv")
  )
}
deg <- read.delim(deg_file, check.names = FALSE)
cat("ORA DEG source:", cfg$downstream_deg_method, "\n")

# Build ORA foreground
up_entrez <- unique(as.character(deg$ENTREZID[deg$logFC > 0 & !is.na(deg$ENTREZID) & deg$ENTREZID != ""]))
down_entrez <- unique(as.character(deg$ENTREZID[deg$logFC < 0 & !is.na(deg$ENTREZID) & deg$ENTREZID != ""]))
cat("ORA background:", length(background_entrez), "genes\n")
cat("ORA foreground: Up", length(up_entrez), "Down", length(down_entrez), "\n")
if (!length(up_entrez)) cat("No UP genes available for ORA. GO/KEGG UP skipped.\n")
if (!length(down_entrez)) cat("No DOWN genes available for ORA. GO/KEGG DOWN skipped.\n")

# Run the four requested ORA analyses with identical output handling
run_ora <- function(genes, database, direction) {
  analysis_name <- paste(database, direction, sep = "_")
  output_base <- paste(prefix, analysis_name, sep = "_")
  empty_result <- data.frame(
    ID = character(), Description = character(), GeneRatio = character(),
    Count = integer(), p.adjust = numeric(), stringsAsFactors = FALSE
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
  significant <- if (nrow(all_result)) all_result[!is.na(all_result$p.adjust) & all_result$p.adjust < cfg$enrichment_FDR_cutoff, , drop = FALSE] else all_result
  write.table(all_result, file.path(result_dir, paste0(output_base, "_all.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  write.table(significant, file.path(result_dir, paste0(output_base, "_FDR", cfg$enrichment_FDR_cutoff, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  if (nrow(significant)) {
    cat(analysis_name, ":", nrow(significant), "significant terms\n")
  } else {
    cat("No significant enrichment result for", analysis_name, ".\n")
  }
}

run_ora(up_entrez, "GO_BP", "UP")
run_ora(down_entrez, "GO_BP", "DOWN")
run_ora(up_entrez, "KEGG", "UP")
run_ora(down_entrez, "KEGG", "DOWN")
