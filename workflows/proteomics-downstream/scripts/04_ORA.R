args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("04_ORA.R requires one config path")
for (package in c("yaml", "clusterProfiler", "org.Hs.eg.db")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing R package: ", package)
}

config <- yaml::read_yaml(normalizePath(args[[1]], mustWork = TRUE))
project_dir <- normalizePath(config$project_dir, mustWork = TRUE)
processed_dir <- file.path(project_dir, "datasets", "processed")
result_dir <- file.path(project_dir, "results", "ORA")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
data <- utils::read.delim(
  file.path(processed_dir, paste0(config$dataset_id, "_all_proteins.tsv")),
  check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA")
)
flags <- readRDS(file.path(processed_dir, paste0(config$dataset_id, "_QC_flags.rds")))
if (!isTRUE(flags$entrez_valid_for_enrichment)) {
  stop("ORA stopped: GeneID values failed the Entrez-format/mapping QC. Volcano output remains valid.")
}

data$logFC <- suppressWarnings(as.numeric(data$logFC))
data$P.Value <- suppressWarnings(as.numeric(data$P.Value))
data$Q.Value <- suppressWarnings(as.numeric(data$Q.Value))
valid_entrez <- !is.na(data$ENTREZID) & grepl("^[0-9]+$", data$ENTREZID)

# The ORA universe is the mapped part of this complete unscreened protein table.
background_rows <- data[valid_entrez, c("ENTREZID", "PROTEIN_ID", "SYMBOL", "DESCRIPTION"), drop = FALSE]
background_rows <- background_rows[!duplicated(background_rows$ENTREZID), , drop = FALSE]
background <- background_rows$ENTREZID
utils::write.table(background_rows, file.path(result_dir, "ORA_background.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

make_foreground <- function(evidence, direction) {
  evidence_pass <- if (evidence == "formal") {
    is.finite(data$Q.Value) & data$Q.Value < config$enrichment_qvalue_cutoff
  } else {
    is.finite(data$P.Value) & data$P.Value < config$enrichment_pvalue_cutoff
  }
  direction_pass <- if (direction == "UP") data$logFC >= config$enrichment_gene_logFC_cutoff else data$logFC <= -config$enrichment_gene_logFC_cutoff
  rows <- which(valid_entrez & evidence_pass & is.finite(data$logFC) & direction_pass)
  if (!length(rows)) return(data[FALSE, c("ENTREZID", "PROTEIN_ID", "SYMBOL", "DESCRIPTION", "logFC", "P.Value", "Q.Value"), drop = FALSE])
  metric <- if (evidence == "formal") data$Q.Value[rows] else data$P.Value[rows]
  rows <- rows[order(metric, -abs(data$logFC[rows]), na.last = TRUE)]
  selected <- data[rows, c("ENTREZID", "PROTEIN_ID", "SYMBOL", "DESCRIPTION", "logFC", "P.Value", "Q.Value"), drop = FALSE]
  selected[!duplicated(selected$ENTREZID), , drop = FALSE]
}

foregrounds <- list()
for (evidence in c("formal", "exploratory")) {
  for (direction in c("UP", "DOWN")) {
    key <- paste(evidence, direction, sep = "_")
    foregrounds[[key]] <- make_foreground(evidence, direction)
    utils::write.table(
      foregrounds[[key]], file.path(result_dir, paste0("ORA_foreground_", key, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE, na = ""
    )
  }
}

empty_result <- data.frame(
  ID = character(), Description = character(), GeneRatio = character(), BgRatio = character(),
  pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(), geneID = character(),
  Count = integer(), stringsAsFactors = FALSE, check.names = FALSE
)
run_ora <- function(gene, database) {
  if (length(gene) < 2L) return(empty_result)
  object <- tryCatch({
    if (database == "GO_BP") {
      clusterProfiler::enrichGO(
        gene = gene, universe = background, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1,
        minGSSize = config$min_GS_size, maxGSSize = config$max_GS_size,
        readable = FALSE
      )
    } else {
      clusterProfiler::enrichKEGG(
        gene = gene, universe = background, organism = "hsa", keyType = "kegg",
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
        minGSSize = config$min_GS_size, maxGSSize = config$max_GS_size,
        use_internal_data = FALSE
      )
    }
  }, error = function(error) {
    warning(database, " ORA returned an error: ", conditionMessage(error))
    NULL
  })
  if (is.null(object)) return(empty_result)
  result <- as.data.frame(object)
  if (!nrow(result)) empty_result else result
}
cutoff_tag <- function(prefix, value) paste0(prefix, gsub("\\.", "p", format(value, scientific = FALSE, trim = TRUE)))
p_tag <- cutoff_tag("P", config$enrichment_pvalue_cutoff)
fdr_tag <- cutoff_tag("FDR", config$enrichment_qvalue_cutoff)

mapping_rows <- list(data.frame(
  evidence = "background", direction = "ALL", input_proteins = sum(valid_entrez),
  unique_mapped_entrez = length(background), stringsAsFactors = FALSE
))
for (evidence in c("formal", "exploratory")) {
  for (direction in c("UP", "DOWN")) {
    foreground <- foregrounds[[paste(evidence, direction, sep = "_")]]
    mapping_rows[[length(mapping_rows) + 1L]] <- data.frame(
      evidence = evidence, direction = direction,
      input_proteins = nrow(foreground), unique_mapped_entrez = length(unique(foreground$ENTREZID)),
      stringsAsFactors = FALSE
    )
    for (database in c("GO_BP", "KEGG")) {
      result <- run_ora(unique(foreground$ENTREZID), database)
      stem <- paste(config$dataset_id, "ORA", evidence, database, direction, sep = "_")
      utils::write.table(result, file.path(result_dir, paste0(stem, "_all.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      p_result <- result[is.finite(result$pvalue) & result$pvalue < config$enrichment_pvalue_cutoff, , drop = FALSE]
      fdr_result <- result[is.finite(result$p.adjust) & result$p.adjust < config$enrichment_qvalue_cutoff, , drop = FALSE]
      utils::write.table(p_result, file.path(result_dir, paste0(stem, "_", p_tag, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      utils::write.table(fdr_result, file.path(result_dir, paste0(stem, "_", fdr_tag, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      cat(evidence, database, direction, ": foreground", nrow(foreground), "; pathways", nrow(result), "\n")
    }
  }
}
mapping_summary <- do.call(rbind, mapping_rows)
utils::write.table(
  mapping_summary, file.path(result_dir, paste0(config$dataset_id, "_ORA_mapping_summary.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
