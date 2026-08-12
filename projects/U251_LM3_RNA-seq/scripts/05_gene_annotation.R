# ============================================================
# 05_gene_annotation.R
#
# Purpose:
#   Annotate Ensembl gene IDs with:
#   - Gene Symbol
#   - Entrez ID
#   - Gene Name
#
# Input:
#   datasets/processed/U251_expression_log2_filtered.rds
#
# Output:
#   datasets/processed/U251_gene_annotation.tsv
# ============================================================


library(AnnotationDbi)
library(org.Hs.eg.db)


# ============================================================
# 1. Load all genes that entered differential analysis
# ============================================================

expr <- readRDS(
  "datasets/processed/U251_expression_log2_filtered.rds"
)

gene_ids <- rownames(expr)

cat("\nTotal genes for annotation:\n")
print(length(gene_ids))

cat("\nFirst Ensembl IDs:\n")
print(head(gene_ids))

# ============================================================
# 2. Ensembl -> Gene Symbol
# ============================================================

gene_symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)


# ============================================================
# 3. Ensembl -> Entrez ID
# ============================================================

entrez_id <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)


# ============================================================
# 4. Ensembl -> Gene Name
# ============================================================

gene_name <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  keytype = "ENSEMBL",
  column = "GENENAME",
  multiVals = "first"
)
# ============================================================
# 5. Build master annotation table
# ============================================================

annotation <- data.frame(
  gene_id = gene_ids,

  SYMBOL = unname(
    gene_symbol[gene_ids]
  ),

  ENTREZID = unname(
    entrez_id[gene_ids]
  ),

  GENENAME = unname(
    gene_name[gene_ids]
  ),

  stringsAsFactors = FALSE
)

cat("\nAnnotation preview:\n")
print(head(annotation))

# ============================================================
# 6. Check annotation success
# ============================================================

total_genes <- nrow(annotation)

symbol_mapped <- sum(
  !is.na(annotation$SYMBOL)
)

entrez_mapped <- sum(
  !is.na(annotation$ENTREZID)
)

genename_mapped <- sum(
  !is.na(annotation$GENENAME)
)


cat("\n========================================\n")
cat("Annotation summary\n")
cat("========================================\n")

cat("\nTotal genes:\n")
print(total_genes)

cat("\nSYMBOL mapped:\n")
print(symbol_mapped)

cat("\nSYMBOL mapping rate:\n")
print(symbol_mapped / total_genes)

cat("\nENTREZID mapped:\n")
print(entrez_mapped)

cat("\nENTREZID mapping rate:\n")
print(entrez_mapped / total_genes)

cat("\nGENENAME mapped:\n")
print(genename_mapped)

cat("\nGENENAME mapping rate:\n")
print(genename_mapped / total_genes)
unmapped <- annotation[
  is.na(annotation$SYMBOL),
]

cat("\nUnmapped genes:\n")
print(nrow(unmapped))

cat("\nFirst unmapped Ensembl IDs:\n")
print(head(unmapped, 20))

# ============================================================
# 7. Save project-level annotation table
# ============================================================

write.table(
  annotation,
  "datasets/processed/U251_gene_annotation.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nGene annotation completed successfully.\n")