args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Expected one config.yml argument")
cfg <- yaml::read_yaml(normalizePath(args[1], mustWork = TRUE))
project_dir <- normalizePath(cfg$project_dir, mustWork = TRUE)
expr_file <- file.path(project_dir, "datasets", "processed", paste0(cfg$dataset_id, "_bulk_workflow_expression_log2_filtered.rds"))
gene_ids <- rownames(readRDS(expr_file))

# Map Ensembl IDs while preserving ENSEMBL as the primary key
annotation <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = gene_ids,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID", "GENENAME")
)
annotation <- annotation[!duplicated(annotation$ENSEMBL), c("ENSEMBL", "SYMBOL", "ENTREZID", "GENENAME")]
annotation <- merge(data.frame(ENSEMBL = gene_ids), annotation, by = "ENSEMBL", all.x = TRUE, sort = FALSE)
annotation <- annotation[match(gene_ids, annotation$ENSEMBL), ]
annotation_file <- file.path(project_dir, "datasets", "processed", paste0(cfg$dataset_id, "_bulk_workflow_gene_annotation.tsv"))
write.table(annotation, annotation_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("Annotation mapping rates: SYMBOL", mean(!is.na(annotation$SYMBOL) & annotation$SYMBOL != ""), "ENTREZID", mean(!is.na(annotation$ENTREZID) & annotation$ENTREZID != ""), "GENENAME", mean(!is.na(annotation$GENENAME) & annotation$GENENAME != ""), "\n")
cat("Saved:", annotation_file, "\n")
