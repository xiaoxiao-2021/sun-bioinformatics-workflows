args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Expected one config.yml argument")
cfg <- yaml::read_yaml(normalizePath(args[1], mustWork = TRUE))
project_dir <- normalizePath(cfg$project_dir, mustWork = TRUE)
de_dir <- file.path(project_dir, "results", cfg$dataset_id, "DE")
annotation_file <- file.path(project_dir, "datasets", "processed", paste0(cfg$dataset_id, "_bulk_workflow_gene_annotation.tsv"))
annotation <- read.delim(annotation_file, check.names = FALSE, stringsAsFactors = FALSE)
if (anyDuplicated(annotation$ENSEMBL)) stop("Master annotation contains duplicated ENSEMBL keys")
prefix <- paste(cfg$dataset_id, cfg$case_group, "vs", cfg$control_group, sep = "_")
lfc_tag <- format(cfg$TREAT_lfc_cutoff, trim = TRUE, scientific = FALSE)
fdr_tag <- format(cfg$TREAT_FDR_cutoff, trim = TRUE, scientific = FALSE)
deg_fdr_tag <- format(cfg$deg_FDR_cutoff, trim = TRUE, scientific = FALSE)
deg_lfc_tag <- format(cfg$deg_logFC_cutoff, trim = TRUE, scientific = FALSE)
pvalue_tag <- format(cfg$pvalue_cutoff, trim = TRUE, scientific = FALSE)
deg_base <- paste0(prefix, "_limma_FDR", deg_fdr_tag, "_logFC", deg_lfc_tag)
pvalue_base <- paste0(prefix, "_limma_P", pvalue_tag, "_logFC", deg_lfc_tag)
files <- c(
  paste0(prefix, "_limma_all_genes.tsv"),
  paste0(deg_base, "_DEG.tsv"),
  paste0(deg_base, "_UP.tsv"),
  paste0(deg_base, "_DOWN.tsv"),
  paste0(pvalue_base, "_DEG.tsv"),
  paste0(pvalue_base, "_UP.tsv"),
  paste0(pvalue_base, "_DOWN.tsv"),
  paste0(prefix, "_TREAT_lfc", lfc_tag, "_all_genes.tsv"),
  paste0(prefix, "_TREAT_lfc", lfc_tag, "_FDR", fdr_tag, "_DEG.tsv")
)

# Add the same master annotation to all DE outputs
for (filename in files) {
  input <- file.path(de_dir, filename)
  result <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)
  original_order <- result$ENSEMBL
  result <- merge(result, annotation, by = "ENSEMBL", all.x = TRUE, sort = FALSE)
  result <- result[match(original_order, result$ENSEMBL), ]
  output <- file.path(de_dir, sub("\\.tsv$", "_annotated.tsv", filename))
  write.table(result, output, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  cat("Saved:", output, "\n")
}
