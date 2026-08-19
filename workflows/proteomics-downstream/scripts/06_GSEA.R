args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("06_GSEA.R requires one config path")
for (package in c("yaml", "limma", "clusterProfiler", "org.Hs.eg.db", "msigdbr")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing R package: ", package)
}

config <- yaml::read_yaml(normalizePath(args[[1]], mustWork = TRUE))
project_dir <- normalizePath(config$project_dir, mustWork = TRUE)
processed_dir <- file.path(project_dir, "datasets", "processed")
result_dir <- file.path(project_dir, "results", "GSEA")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
data <- utils::read.delim(
  file.path(processed_dir, paste0(config$dataset_id, "_all_proteins.tsv")),
  check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA")
)
metadata <- utils::read.delim(file.path(project_dir, config$sample_metadata), check.names = FALSE, stringsAsFactors = FALSE)
flags <- readRDS(file.path(processed_dir, paste0(config$dataset_id, "_QC_flags.rds")))
if (!isTRUE(flags$abundance_scale_ok)) {
  stop("GSEA stopped: abundance values failed the log2-scale/vendor-logFC agreement check; values were not transformed.")
}
if (!isTRUE(flags$entrez_valid_for_enrichment)) {
  stop("GSEA stopped: GeneID values failed the Entrez-format/mapping QC.")
}

sample_ids <- metadata$sample_id
abundance <- as.matrix(data[sample_ids])
storage.mode(abundance) <- "numeric"
case_columns <- which(metadata$group == config$case_group)
control_columns <- which(metadata$group == config$control_group)
keep <- rowSums(is.finite(abundance[, case_columns, drop = FALSE])) >= 2L &
  rowSums(is.finite(abundance[, control_columns, drop = FALSE])) >= 2L
if (sum(keep) < 10L) stop("GSEA stopped: fewer than 10 proteins have sufficient abundance replicates for limma")

# limma is used only to construct a stable signed GSEA ranking statistic.
group <- factor(metadata$group, levels = c(config$control_group, config$case_group))
design <- stats::model.matrix(~0 + group)
colnames(design) <- make.names(levels(group))
contrast_text <- paste0(make.names(config$case_group), "-", make.names(config$control_group))
contrast <- limma::makeContrasts(contrasts = contrast_text, levels = design)
fit <- limma::lmFit(abundance[keep, , drop = FALSE], design)
fit <- limma::contrasts.fit(fit, contrast)
fit <- limma::eBayes(fit)
limma_table <- limma::topTable(fit, number = Inf, sort.by = "none")
source_rows <- which(keep)

protein_stats <- data.frame(
  PROTEIN_ID = data$PROTEIN_ID[source_rows],
  SYMBOL = data$SYMBOL[source_rows],
  ENTREZID = data$ENTREZID[source_rows],
  limma_logFC = limma_table$logFC,
  limma_t = limma_table$t,
  limma_P.Value = limma_table$P.Value,
  limma_adj.P.Val = limma_table$adj.P.Val,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
protein_stats_file <- file.path(result_dir, paste0(config$dataset_id, "_GSEA_limma_protein_statistics.tsv"))
utils::write.table(protein_stats, protein_stats_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

vendor_logfc <- suppressWarnings(as.numeric(data$logFC[source_rows]))
comparable <- is.finite(vendor_logfc) & is.finite(protein_stats$limma_logFC)
direction_rows <- comparable & vendor_logfc != 0 & protein_stats$limma_logFC != 0
comparison <- data.frame(
  metric = c("compared_proteins", "pearson_correlation", "direction_concordance"),
  value = c(
    sum(comparable),
    if (sum(comparable) >= 3L) stats::cor(vendor_logfc[comparable], protein_stats$limma_logFC[comparable]) else NA_real_,
    if (any(direction_rows)) mean(sign(vendor_logfc[direction_rows]) == sign(protein_stats$limma_logFC[direction_rows])) else NA_real_
  ),
  stringsAsFactors = FALSE
)
utils::write.table(
  comparison, file.path(result_dir, paste0(config$dataset_id, "_vendor_vs_limma_logFC.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

valid_rank <- !is.na(protein_stats$ENTREZID) & grepl("^[0-9]+$", protein_stats$ENTREZID) &
  is.finite(protein_stats$limma_t) & is.finite(protein_stats$limma_P.Value)
rank_rows <- which(valid_rank)
# One protein represents each Entrez ID: largest |moderated t|, then smallest P.
rank_rows <- rank_rows[order(-abs(protein_stats$limma_t[rank_rows]), protein_stats$limma_P.Value[rank_rows])]
rank_rows <- rank_rows[!duplicated(protein_stats$ENTREZID[rank_rows])]
ranked_table <- protein_stats[rank_rows, , drop = FALSE]
ranked_table <- ranked_table[order(ranked_table$limma_t, decreasing = TRUE), , drop = FALSE]
ranked_file <- file.path(result_dir, paste0(config$dataset_id, "_GSEA_ranked_genes.tsv"))
utils::write.table(ranked_table, ranked_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
if (nrow(ranked_table) < config$gsea_min_GS_size) {
  stop("GSEA stopped: ranked Entrez list is smaller than gsea_min_GS_size")
}
gene_list <- ranked_table$limma_t
names(gene_list) <- ranked_table$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
# clusterProfiler 4.8/fgsea reads .Random.seed when seed = TRUE. Initialize it
# explicitly so a fresh Rscript session is reproducible instead of erroring.
set.seed(1L)

empty_result <- data.frame(
  ID = character(), Description = character(), setSize = integer(),
  enrichmentScore = numeric(), NES = numeric(), pvalue = numeric(),
  p.adjust = numeric(), qvalue = numeric(), rank = integer(),
  leading_edge = character(), core_enrichment = character(),
  stringsAsFactors = FALSE, check.names = FALSE
)
run_gsea <- function(database) {
  object <- tryCatch({
    if (database == "Hallmark") {
      msig_args <- list(species = "Homo sapiens")
      msig_formals <- names(formals(msigdbr::msigdbr))
      if ("collection" %in% msig_formals) msig_args$collection <- "H" else msig_args$category <- "H"
      hallmark <- do.call(msigdbr::msigdbr, msig_args)
      gene_column <- if ("ncbi_gene" %in% names(hallmark)) "ncbi_gene" else "entrez_gene"
      # Save the exact deduplicated TERM2GENE object used by this Hallmark GSEA.
      # This snapshot keeps downstream enrichment curves/GSEA visualization
      # aligned with the Hallmark gene-set membership used for this run.
      hallmark_t2g <- unique(data.frame(
        TERM = hallmark$gs_name,
        GENE = as.character(hallmark[[gene_column]]),
        stringsAsFactors = FALSE,
        check.names = FALSE
      ))
      term2gene_output <- file.path(
        result_dir,
        paste0(config$dataset_id, "_GSEA_Hallmark_TERM2GENE.tsv")
      )
      utils::write.table(
        hallmark_t2g,
        file = term2gene_output,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
      cat("Saved Hallmark TERM2GENE snapshot:\n", term2gene_output, "\n", sep = "")
      cat("Hallmark TERM2GENE rows: ", nrow(hallmark_t2g), "\n", sep = "")
      cat("Hallmark gene sets: ", length(unique(hallmark_t2g$TERM)), "\n", sep = "")
      cat("Hallmark unique genes: ", length(unique(hallmark_t2g$GENE)), "\n", sep = "")
      clusterProfiler::GSEA(
        geneList = gene_list, TERM2GENE = hallmark_t2g, pvalueCutoff = 1,
        pAdjustMethod = "BH", minGSSize = config$gsea_min_GS_size,
        maxGSSize = config$gsea_max_GS_size, eps = 0, verbose = FALSE,
        seed = TRUE, by = "fgsea"
      )
    } else if (database == "GO_BP") {
      clusterProfiler::gseGO(
        geneList = gene_list, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID", ont = "BP", pvalueCutoff = 1,
        pAdjustMethod = "BH", minGSSize = config$gsea_min_GS_size,
        maxGSSize = config$gsea_max_GS_size, eps = 0, verbose = FALSE,
        seed = TRUE, by = "fgsea"
      )
    } else {
      clusterProfiler::gseKEGG(
        geneList = gene_list, organism = "hsa", keyType = "ncbi-geneid",
        pvalueCutoff = 1, pAdjustMethod = "BH",
        minGSSize = config$gsea_min_GS_size, maxGSSize = config$gsea_max_GS_size,
        eps = 0, verbose = FALSE, seed = TRUE, by = "fgsea",
        use_internal_data = FALSE
      )
    }
  }, error = function(error) {
    warning(database, " GSEA returned an error: ", conditionMessage(error))
    NULL
  })
  object
}
cutoff_tag <- function(prefix, value) paste0(prefix, gsub("\\.", "p", format(value, scientific = FALSE, trim = TRUE)))
p_tag <- cutoff_tag("P", config$gsea_pvalue_cutoff)
fdr_tag <- cutoff_tag("FDR", config$gsea_FDR_cutoff)

for (database in c("Hallmark", "GO_BP", "KEGG")) {
  object <- run_gsea(database)
  saveRDS(object, file.path(result_dir, paste0(config$dataset_id, "_GSEA_", database, "_object.rds")))
  result <- if (is.null(object)) empty_result else as.data.frame(object)
  if (!nrow(result)) result <- empty_result
  if (nrow(result)) {
    result$enriched_in <- ifelse(result$NES > 0, config$case_group, config$control_group)
  } else {
    result$enriched_in <- character()
  }
  stem <- paste0(config$dataset_id, "_GSEA_", database)
  utils::write.table(result, file.path(result_dir, paste0(stem, "_all.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  p_result <- result[is.finite(result$pvalue) & result$pvalue < config$gsea_pvalue_cutoff, , drop = FALSE]
  fdr_result <- result[is.finite(result$p.adjust) & result$p.adjust < config$gsea_FDR_cutoff, , drop = FALSE]
  utils::write.table(p_result, file.path(result_dir, paste0(stem, "_", p_tag, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  utils::write.table(fdr_result, file.path(result_dir, paste0(stem, "_", fdr_tag, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  leading_rows <- list()
  if (nrow(result) && "core_enrichment" %in% names(result)) {
    for (i in seq_len(nrow(result))) {
      genes <- unique(strsplit(as.character(result$core_enrichment[[i]]), "/", fixed = TRUE)[[1]])
      genes <- genes[nzchar(genes)]
      if (!length(genes)) next
      annotation <- ranked_table[match(genes, ranked_table$ENTREZID), , drop = FALSE]
      annotation$database <- database
      annotation$pathway_id <- result$ID[[i]]
      annotation$pathway_description <- result$Description[[i]]
      annotation$pathway_NES <- result$NES[[i]]
      annotation$pathway_pvalue <- result$pvalue[[i]]
      annotation$pathway_p_adjust <- result$p.adjust[[i]]
      annotation$enriched_in <- if (result$NES[[i]] > 0) config$case_group else config$control_group
      leading_rows[[length(leading_rows) + 1L]] <- annotation
    }
  }
  leading <- if (length(leading_rows)) do.call(rbind, leading_rows) else data.frame(
    PROTEIN_ID = character(), SYMBOL = character(), ENTREZID = character(),
    limma_logFC = numeric(), limma_t = numeric(), limma_P.Value = numeric(),
    limma_adj.P.Val = numeric(), database = character(), pathway_id = character(),
    pathway_description = character(), pathway_NES = numeric(), pathway_pvalue = numeric(),
    pathway_p_adjust = numeric(), enriched_in = character(), stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.table(leading, file.path(result_dir, paste0(stem, "_leading_edge.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  cat(database, "GSEA pathways:", nrow(result), "; formal:", nrow(fdr_result), "; exploratory:", nrow(p_result), "\n")
}
