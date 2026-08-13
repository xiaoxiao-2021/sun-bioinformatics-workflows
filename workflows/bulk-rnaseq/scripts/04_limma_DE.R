args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Expected one config.yml argument")
cfg <- yaml::read_yaml(normalizePath(args[1], mustWork = TRUE))
project_dir <- normalizePath(cfg$project_dir, mustWork = TRUE)
metadata <- read.delim(file.path(project_dir, cfg$sample_metadata), check.names = FALSE, stringsAsFactors = FALSE)
expr_file <- file.path(project_dir, "datasets", "processed", paste0(cfg$dataset_id, "_bulk_workflow_expression_log2_filtered.rds"))
expr <- readRDS(expr_file)
metadata <- metadata[match(colnames(expr), metadata$sample), , drop = FALSE]
group <- factor(metadata$group, levels = c(cfg$control_group, cfg$case_group))

# Build case vs control design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
contrast <- matrix(0, nrow = ncol(design), ncol = 1, dimnames = list(colnames(design), paste0(cfg$case_group, "_vs_", cfg$control_group)))
contrast[cfg$case_group, 1] <- 1
contrast[cfg$control_group, 1] <- -1
cat("Comparison:", cfg$case_group, "-", cfg$control_group, "\n")

# Fit case - control contrast
fit_raw <- limma::lmFit(expr, design)
fit_contrast <- limma::contrasts.fit(fit_raw, contrast)

# Standard limma test
fit_ebayes <- limma::eBayes(fit_contrast)
result <- limma::topTable(fit_ebayes, coef = 1, number = Inf, adjust.method = "BH", sort.by = "P")
result$ENSEMBL <- rownames(result)
result <- result[, c("ENSEMBL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

# Standard limma DEG
deg_limma <- result[
  result$adj.P.Val < cfg$deg_FDR_cutoff &
    abs(result$logFC) >= cfg$deg_logFC_cutoff,
  , drop = FALSE
]
deg_limma$change <- ifelse(deg_limma$logFC > 0, "Up", "Down")

# Split DEG by direction
deg_limma_up <- deg_limma[deg_limma$change == "Up", , drop = FALSE]
deg_limma_down <- deg_limma[deg_limma$change == "Down", , drop = FALSE]

# Nominal P-value exploratory DEG
deg_pvalue <- result[
  result$P.Value < cfg$pvalue_cutoff &
    abs(result$logFC) >= cfg$deg_logFC_cutoff,
  , drop = FALSE
]
deg_pvalue$change <- ifelse(deg_pvalue$logFC > 0, "Up", "Down")
deg_pvalue_up <- deg_pvalue[deg_pvalue$change == "Up", , drop = FALSE]
deg_pvalue_down <- deg_pvalue[deg_pvalue$change == "Down", , drop = FALSE]

# Strict TREAT test; DEG uses FDR only because lfc is in the hypothesis
fit_treat <- limma::treat(fit_contrast, lfc = cfg$TREAT_lfc_cutoff)
result_treat <- limma::topTreat(fit_treat, coef = 1, number = Inf, adjust.method = "BH", sort.by = "P")
result_treat$ENSEMBL <- rownames(result_treat)
result_treat <- result_treat[, c("ENSEMBL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
deg_treat <- result_treat[result_treat$adj.P.Val < cfg$TREAT_FDR_cutoff, , drop = FALSE]
deg_treat$change <- ifelse(deg_treat$logFC > 0, "Up", "Down")

de_dir <- file.path(project_dir, "results", cfg$dataset_id, "DE")
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- paste(cfg$dataset_id, cfg$case_group, "vs", cfg$control_group, sep = "_")
lfc_tag <- format(cfg$TREAT_lfc_cutoff, trim = TRUE, scientific = FALSE)
fdr_tag <- format(cfg$TREAT_FDR_cutoff, trim = TRUE, scientific = FALSE)
deg_fdr_tag <- format(cfg$deg_FDR_cutoff, trim = TRUE, scientific = FALSE)
deg_lfc_tag <- format(cfg$deg_logFC_cutoff, trim = TRUE, scientific = FALSE)
pvalue_tag <- format(cfg$pvalue_cutoff, trim = TRUE, scientific = FALSE)
deg_base <- paste0(prefix, "_limma_FDR", deg_fdr_tag, "_logFC", deg_lfc_tag)
pvalue_base <- paste0(prefix, "_limma_P", pvalue_tag, "_logFC", deg_lfc_tag)
write.table(result, file.path(de_dir, paste0(prefix, "_limma_all_genes.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_limma, file.path(de_dir, paste0(deg_base, "_DEG.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_limma_up, file.path(de_dir, paste0(deg_base, "_UP.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_limma_down, file.path(de_dir, paste0(deg_base, "_DOWN.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_pvalue, file.path(de_dir, paste0(pvalue_base, "_DEG.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_pvalue_up, file.path(de_dir, paste0(pvalue_base, "_UP.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_pvalue_down, file.path(de_dir, paste0(pvalue_base, "_DOWN.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(result_treat, file.path(de_dir, paste0(prefix, "_TREAT_lfc", lfc_tag, "_all_genes.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_treat, file.path(de_dir, paste0(prefix, "_TREAT_lfc", lfc_tag, "_FDR", fdr_tag, "_DEG.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Standard limma FDR hits:", sum(result$adj.P.Val < cfg$limma_FDR_cutoff), "\n")
cat("Standard limma DEG summary:\n")
cat("Up:", nrow(deg_limma_up), "\n")
cat("Down:", nrow(deg_limma_down), "\n")
cat("Total:", nrow(deg_limma), "\n")
cat("Nominal P-value exploratory DEG summary:\n")
cat("Up:", nrow(deg_pvalue_up), "\n")
cat("Down:", nrow(deg_pvalue_down), "\n")
cat("Total:", nrow(deg_pvalue), "\n")
cat("Strict TREAT DEG summary:\n")
cat("Up:", sum(deg_treat$change == "Up"), "\n")
cat("Down:", sum(deg_treat$change == "Down"), "\n")
cat("Total:", nrow(deg_treat), "\n")
