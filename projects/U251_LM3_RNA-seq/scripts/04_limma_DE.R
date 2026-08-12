library(limma)

# ============================================================
# 04_limma_DE.R
#
# Purpose:
#   Differential expression analysis: OE vs NC
#
# Input:
#   log2(expression + 1) cleaned expression matrix
#
# Strategy:
#   1. lmFit()       -> linear model
#   2. eBayes()      -> standard differential-expression ranking
#   3. treat(lfc=1)  -> formal test for effect size exceeding threshold
#
# Comparison:
#   OE - NC
#
# Interpretation:
#   logFC > 0 -> higher in OE
#   logFC < 0 -> lower in OE
# ============================================================


# ============================================================
# 1. Read cleaned expression matrix
# ============================================================

expr <- readRDS(
  "datasets/processed/U251_expression_log2_filtered.rds"
)

cat("\nExpression matrix:\n")
print(dim(expr))

cat("\nSamples:\n")
print(colnames(expr))


# ============================================================
# 2. Define experimental groups
# ============================================================

group <- factor(
  c("OE", "OE", "OE", "NC", "NC", "NC"),
  levels = c("NC", "OE")
)

cat("\nGroup information:\n")
print(group)


# ============================================================
# 3. Build design matrix
#
# NC is the reference group.
# Therefore:
#   Intercept = mean expression of NC
#   groupOE   = OE - NC
# ============================================================

design <- model.matrix(~ group)

cat("\nDesign matrix:\n")
print(design)


# ============================================================
# 4. Fit linear model
# ============================================================

fit_raw <- lmFit(
  expr,
  design
)


# ============================================================
# 5. Standard empirical Bayes analysis
#
# This analysis tests:
#
#   H0: logFC = 0
#
# It is used for:
#   - complete gene ranking
#   - GSEA
#   - exploratory enrichment
#   - general inspection of differential expression
# ============================================================

fit_ebayes <- eBayes(
  fit_raw
)

result <- topTable(
  fit_ebayes,
  coef = "groupOE",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)


# ============================================================
# 6. Add gene IDs
# ============================================================

result$gene_id <- rownames(result)

result <- result[
  ,
  c(
    "gene_id",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )
]


# ============================================================
# 7. Create output directory
# ============================================================

dir.create(
  "results/DE",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 8. Save complete standard limma result
#
# Contains all tested genes.
#
# Recommended uses:
#   - GSEA
#   - ranked gene analysis
#   - checking individual genes
# ============================================================

write.table(
  result,
  "results/DE/U251_OE_vs_NC_limma_all_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 9. Exploratory enrichment candidate set
#
# FDR < 0.05
# |logFC| >= 0.5
#
# This is a relatively permissive candidate set intended
# to retain moderate expression changes for exploratory ORA
# such as GO / KEGG enrichment.
#
# IMPORTANT:
# GSEA should use the complete ranked gene list instead.
# ============================================================

enrichment_candidates <- result[
  result$adj.P.Val < 0.05 &
  abs(result$logFC) >= 0.5,
]

cat("\nEnrichment candidate genes:\n")
print(nrow(enrichment_candidates))


# Split by direction

enrichment_up <- enrichment_candidates[
  enrichment_candidates$logFC >= 0.5,
]

enrichment_down <- enrichment_candidates[
  enrichment_candidates$logFC <= -0.5,
]

cat("\nEnrichment UP genes:\n")
print(nrow(enrichment_up))

cat("\nEnrichment DOWN genes:\n")
print(nrow(enrichment_down))


# Save enrichment candidate sets

write.table(
  enrichment_candidates,
  "results/DE/U251_OE_vs_NC_enrichment_FDR0.05_logFC0.5.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  enrichment_up,
  "results/DE/U251_OE_vs_NC_enrichment_UP_FDR0.05_logFC0.5.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  enrichment_down,
  "results/DE/U251_OE_vs_NC_enrichment_DOWN_FDR0.05_logFC0.5.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 10. TREAT analysis
#
# Instead of testing:
#
#   H0: logFC = 0
#
# TREAT tests whether the true effect exceeds
# the specified minimum effect size.
#
# With lfc = 1:
#
#   H0: |logFC| <= 1
#   H1: |logFC| > 1
#
# NOTE:
# Because our matrix is log2(expression + 1),
# lfc = 1 means a one-unit difference on this transformed
# scale. It approximates a two-fold difference for genes
# with sufficiently high expression, but is not exactly
# equivalent to raw-expression two-fold change at low levels.
# ============================================================

TREAT_lfc_cutoff <- 1
TREAT_FDR_cutoff <- 0.05

fit_treat <- treat(
  fit_raw,
  lfc = TREAT_lfc_cutoff
)


# ============================================================
# 11. Extract complete TREAT results
#
# topTreat() ranks genes according to evidence that their
# true effect exceeds the specified lfc threshold.
#
# TREAT does not produce a B statistic.
# ============================================================

result_treat <- topTreat(
  fit_treat,
  coef = "groupOE",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

result_treat$gene_id <- rownames(result_treat)

result_treat <- result_treat[
  ,
  c(
    "gene_id",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val"
  )
]


# ============================================================
# 12. Save complete TREAT results
# ============================================================

write.table(
  result_treat,
  "results/DE/U251_OE_vs_NC_TREAT_lfc1_all_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 13. Define formal TREAT DEG set
#
# No additional |logFC| cutoff is needed here.
#
# The fold-change threshold has already been incorporated
# into the statistical hypothesis by treat().
# ============================================================

deg_treat <- result_treat[
  result_treat$adj.P.Val < TREAT_FDR_cutoff,
]


# ============================================================
# 14. Classify TREAT DEGs by direction
# ============================================================

deg_treat$change <- ifelse(
  deg_treat$logFC > 0,
  "Up",
  "Down"
)

cat("\nFormal TREAT DEG summary:\n")
print(table(deg_treat$change))

cat("\nTotal formal TREAT DEGs:\n")
print(nrow(deg_treat))


# ============================================================
# 15. Save formal TREAT DEG set
# ============================================================

write.table(
  deg_treat,
  "results/DE/U251_OE_vs_NC_TREAT_lfc1_FDR0.05_DEG.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 16. Analysis summary
# ============================================================

cat("\n========================================\n")
cat("Differential expression analysis done\n")
cat("========================================\n")

cat("\nTotal tested genes:\n")
print(nrow(result))

cat("\nStandard limma FDR < 0.05:\n")
print(
  sum(result$adj.P.Val < 0.05)
)

cat("\nExploratory enrichment candidates:\n")
print(nrow(enrichment_candidates))

cat("\nFormal TREAT DEG:\n")
print(nrow(deg_treat))

cat("\nTREAT DEG direction:\n")
print(table(deg_treat$change))