args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Expected one config.yml argument")
cfg <- yaml::read_yaml(normalizePath(args[1], mustWork = TRUE))
project_dir <- normalizePath(cfg$project_dir, mustWork = TRUE)
expr_file <- file.path(
  project_dir, "datasets", "processed",
  paste0(cfg$dataset_id, "_bulk_workflow_expression_log2_filtered.rds")
)
gene_ids <- rownames(readRDS(expr_file))

# Preserve original ENSEMBL as primary key; use base IDs only for mapping
ensembl_base <- sub("\\.[0-9]+$", "", gene_ids)
mapping_keys <- unique(ensembl_base)

# org.Hs.eg.db annotation
symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = mapping_keys, keytype = "ENSEMBL", column = "SYMBOL",
  multiVals = "first"
)
entrez_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = mapping_keys, keytype = "ENSEMBL", column = "ENTREZID",
  multiVals = "first"
)
genename_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = mapping_keys, keytype = "ENSEMBL", column = "GENENAME",
  multiVals = "first"
)

# Ensembl biotype annotation at gene level
ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
ensdb_annotation <- as.data.frame(
  ensembldb::genes(
    ensdb,
    columns = c("gene_id", "gene_name", "gene_biotype"),
    return.type = "DataFrame"
  )
)
ensdb_annotation <- ensdb_annotation[
  !duplicated(ensdb_annotation$gene_id),
  c("gene_id", "gene_name", "gene_biotype"),
  drop = FALSE
]
ensdb_match <- match(ensembl_base, ensdb_annotation$gene_id)

annotation <- data.frame(
  ENSEMBL = gene_ids,
  SYMBOL = unname(symbol_map[ensembl_base]),
  ENTREZID = unname(entrez_map[ensembl_base]),
  GENENAME = unname(genename_map[ensembl_base]),
  ENSEMBL_GENE_NAME = ensdb_annotation$gene_name[ensdb_match],
  GENE_BIOTYPE = ensdb_annotation$gene_biotype[ensdb_match],
  # EnsDb.Hsapiens.v86 does not contain a gene description column
  ENSEMBL_DESCRIPTION = rep(NA_character_, length(gene_ids)),
  stringsAsFactors = FALSE
)

# Build display name with fallback without overwriting either annotation source
has_symbol <- !is.na(annotation$SYMBOL) & annotation$SYMBOL != ""
has_ensembl_name <- !is.na(annotation$ENSEMBL_GENE_NAME) &
  annotation$ENSEMBL_GENE_NAME != ""
loc_symbol <- has_symbol & grepl("^LOC[0-9]+$", annotation$SYMBOL, ignore.case = TRUE)
readable_ensembl_name <- has_ensembl_name &
  !grepl("^LOC[0-9]+$", annotation$ENSEMBL_GENE_NAME, ignore.case = TRUE) &
  !grepl("^ENSG[0-9]+(\\.[0-9]+)?$", annotation$ENSEMBL_GENE_NAME, ignore.case = TRUE)
prefer_ensembl_name <- loc_symbol & readable_ensembl_name
annotation$DISPLAY_NAME <- annotation$ENSEMBL
annotation$DISPLAY_NAME[has_ensembl_name] <-
  annotation$ENSEMBL_GENE_NAME[has_ensembl_name]
annotation$DISPLAY_NAME[has_symbol] <- annotation$SYMBOL[has_symbol]
annotation$DISPLAY_NAME[prefer_ensembl_name] <-
  annotation$ENSEMBL_GENE_NAME[prefer_ensembl_name]

if (nrow(annotation) != length(gene_ids) || anyDuplicated(annotation$ENSEMBL)) {
  stop("Master annotation must contain exactly one row per original ENSEMBL")
}

annotation_file <- file.path(
  project_dir, "datasets", "processed",
  paste0(cfg$dataset_id, "_bulk_workflow_gene_annotation.tsv")
)
write.table(
  annotation, annotation_file, sep = "\t", quote = FALSE,
  row.names = FALSE, na = ""
)

cat("Total genes:", nrow(annotation), "\n")
cat("SYMBOL mapping rate:", mean(has_symbol), "\n")
cat("ENTREZID mapping rate:", mean(!is.na(annotation$ENTREZID) & annotation$ENTREZID != ""), "\n")
cat("GENENAME mapping rate:", mean(!is.na(annotation$GENENAME) & annotation$GENENAME != ""), "\n")
cat("ENSEMBL_GENE_NAME mapping rate:", mean(has_ensembl_name), "\n")
cat("GENE_BIOTYPE mapping rate:", mean(!is.na(annotation$GENE_BIOTYPE) & annotation$GENE_BIOTYPE != ""), "\n")
cat("Gene biotype distribution:\n")
print(sort(table(annotation$GENE_BIOTYPE, useNA = "ifany"), decreasing = TRUE))
cat("Number of genes using SYMBOL as DISPLAY_NAME:", sum(has_symbol & !prefer_ensembl_name), "\n")
cat("Number using ENSEMBL_GENE_NAME as DISPLAY_NAME:",
    sum((!has_symbol & has_ensembl_name) | prefer_ensembl_name), "\n")
cat("Number of LOC symbols replaced by readable ENSEMBL_GENE_NAME:",
    sum(prefer_ensembl_name), "\n")
cat("Number still using raw ENSEMBL fallback:", sum(!has_symbol & !has_ensembl_name), "\n")
cat("Saved:", annotation_file, "\n")
