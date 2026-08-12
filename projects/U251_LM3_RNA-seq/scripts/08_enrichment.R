# ============================================================
# 08_enrichment.R
#
# Over-Representation Analysis (ORA)
#
# 1. GO Biological Process (GO-BP)
# 2. KEGG pathway enrichment
#
# Up- and down-regulated genes are analyzed separately.
#
#
# Foreground:
#
#   Strict TREAT DEG:
#   results/DE/
#   U251_OE_vs_NC_TREAT_lfc1_FDR0.05_DEG_annotated.tsv
#
#
# Background / Universe:
#
#   All genes tested by limma with successful ENTREZID mapping:
#   results/DE/
#   U251_OE_vs_NC_limma_all_genes_annotated.tsv
#
#
# Output:
#
#   results/enrichment/
#       GO_BP_UP_all.tsv
#       GO_BP_UP_FDR0.05.tsv
#       GO_BP_DOWN_all.tsv
#       GO_BP_DOWN_FDR0.05.tsv
#
#       KEGG_UP_all.tsv
#       KEGG_UP_FDR0.05.tsv
#       KEGG_DOWN_all.tsv
#       KEGG_DOWN_FDR0.05.tsv
#
#
#   figures/enrichment/
#       GO_BP_UP_dotplot.png/pdf
#       GO_BP_DOWN_dotplot.png/pdf
#       KEGG_UP_dotplot.png/pdf
#       KEGG_DOWN_dotplot.png/pdf
#
# ============================================================



# ============================================================
# 1. Environment
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)


cat("\n")
cat("============================================================\n")
cat("08 Enrichment Analysis Started\n")
cat("============================================================\n")



# ============================================================
# 2. Parameters
# ============================================================

FDR_cutoff <- 0.05

logFC_cutoff <- 1

enrichment_FDR_cutoff <- 0.05

show_category_n <- 15


# Gene-set size limits
min_GS_size <- 10
max_GS_size <- 500



# ============================================================
# 3. Paths
# ============================================================

background_file <-
  "results/DE/U251_OE_vs_NC_limma_all_genes_annotated.tsv"


deg_file <-
  "results/DE/U251_OE_vs_NC_TREAT_lfc1_FDR0.05_DEG_annotated.tsv"


result_dir <-
  "results/enrichment"


figure_dir <-
  "figures/enrichment"



# ============================================================
# 4. Create output directories
# ============================================================

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat("\nOutput directories checked:\n")
cat(result_dir, "\n")
cat(figure_dir, "\n")



# ============================================================
# 5. Check input files
# ============================================================

input_files <- c(
  background_file,
  deg_file
)


missing_files <- input_files[
  !file.exists(input_files)
]


if(length(missing_files) > 0){

  stop(
    paste(
      "Missing input file(s):\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )

}


cat("\nAll input files detected.\n")



# ============================================================
# 6. Load data
# ============================================================

background_data <- read.table(

  background_file,

  header = TRUE,

  sep = "\t",

  quote = "",

  stringsAsFactors = FALSE

)



deg <- read.table(

  deg_file,

  header = TRUE,

  sep = "\t",

  quote = "",

  stringsAsFactors = FALSE

)



cat(
  "\nAll limma genes:",
  nrow(background_data),
  "\n"
)


cat(
  "Strict TREAT DEG:",
  nrow(deg),
  "\n"
)



# ============================================================
# 7. Check required columns
# ============================================================

required_background_columns <- c(
  "ENTREZID"
)


required_deg_columns <- c(
  "ENTREZID",
  "SYMBOL",
  "logFC",
  "adj.P.Val"
)



missing_background_columns <- setdiff(

  required_background_columns,

  colnames(background_data)

)


missing_deg_columns <- setdiff(

  required_deg_columns,

  colnames(deg)

)



if(length(missing_background_columns) > 0){

  stop(
    paste(
      "Missing background column(s):",
      paste(
        missing_background_columns,
        collapse = ", "
      )
    )
  )

}



if(length(missing_deg_columns) > 0){

  stop(
    paste(
      "Missing DEG column(s):",
      paste(
        missing_deg_columns,
        collapse = ", "
      )
    )
  )

}



# ============================================================
# 8. Build enrichment background
#
# Universe =
# all genes tested by limma
# AND successfully mapped to ENTREZID
# ============================================================

background_entrez <- background_data %>%

  filter(
    !is.na(ENTREZID),
    ENTREZID != ""
  ) %>%

  pull(
    ENTREZID
  ) %>%

  as.character() %>%

  unique()



cat(
  "\nBackground ENTREZ genes:",
  length(background_entrez),
  "\n"
)



# ============================================================
# 9. Build strict Up / Down gene lists
# ============================================================

up_entrez <- deg %>%

  filter(

    adj.P.Val < FDR_cutoff,

    logFC >= logFC_cutoff,

    !is.na(ENTREZID),

    ENTREZID != ""

  ) %>%

  pull(
    ENTREZID
  ) %>%

  as.character() %>%

  unique()



down_entrez <- deg %>%

  filter(

    adj.P.Val < FDR_cutoff,

    logFC <= -logFC_cutoff,

    !is.na(ENTREZID),

    ENTREZID != ""

  ) %>%

  pull(
    ENTREZID
  ) %>%

  as.character() %>%

  unique()



cat("\nGene lists for ORA:\n")

cat(
  "Up genes with ENTREZID:",
  length(up_entrez),
  "\n"
)

cat(
  "Down genes with ENTREZID:",
  length(down_entrez),
  "\n"
)



# ============================================================
# 10. Helper: convert GeneRatio to numeric
# ============================================================

gene_ratio_numeric <- function(x){

  sapply(

    strsplit(
      x,
      "/",
      fixed = TRUE
    ),

    function(z){

      as.numeric(z[1]) /
        as.numeric(z[2])

    }

  )

}



# ============================================================
# 11. Helper: save enrichment results
# ============================================================

save_enrichment_results <- function(

    enrichment_object,

    analysis_name,

    plot_title

){

  # ----------------------------------------------------------
  # Convert enrichment object to data.frame
  # ----------------------------------------------------------

  result_all <- as.data.frame(
    enrichment_object
  )



  # ----------------------------------------------------------
  # File paths
  # ----------------------------------------------------------

  all_file <- file.path(

    result_dir,

    paste0(
      analysis_name,
      "_all.tsv"
    )

  )


  sig_file <- file.path(

    result_dir,

    paste0(
      analysis_name,
      "_FDR0.05.tsv"
    )

  )



  # ----------------------------------------------------------
  # Save complete result
  # ----------------------------------------------------------

  write.table(

    result_all,

    all_file,

    sep = "\t",

    quote = FALSE,

    row.names = FALSE

  )



  # ----------------------------------------------------------
  # Significant enrichment
  # ----------------------------------------------------------

  if(nrow(result_all) == 0){

    cat(
      "\n",
      analysis_name,
      ": no enrichment terms returned.\n",
      sep = ""
    )

    return(
      invisible(NULL)
    )

  }



  result_sig <- result_all %>%

    filter(

      !is.na(p.adjust),

      p.adjust <
        enrichment_FDR_cutoff

    )



  write.table(

    result_sig,

    sig_file,

    sep = "\t",

    quote = FALSE,

    row.names = FALSE

  )



  cat(
    "\n",
    analysis_name,
    "\n",
    sep = ""
  )

  cat(
    "All terms:",
    nrow(result_all),
    "\n"
  )

  cat(
    "FDR <",
    enrichment_FDR_cutoff,
    ":",
    nrow(result_sig),
    "\n"
  )



  # ----------------------------------------------------------
  # No significant result -> skip plot
  # ----------------------------------------------------------

  if(nrow(result_sig) == 0){

    cat(
      "No significant terms. Dotplot skipped.\n"
    )

    return(
      invisible(NULL)
    )

  }



  # ==========================================================
  # Dotplot data
  # ==========================================================

  plot_data <- result_sig %>%

    arrange(
      p.adjust
    ) %>%

    slice_head(
      n = show_category_n
    ) %>%

    mutate(

      GeneRatio_numeric =
        gene_ratio_numeric(
          GeneRatio
        ),

      minus_log10_FDR =
        -log10(
          pmax(
            p.adjust,
            .Machine$double.xmin
          )
        )

    )



  # Preserve display order
  plot_data$Description <- factor(

    plot_data$Description,

    levels =
      rev(
        plot_data$Description
      )

  )



  # ==========================================================
  # Dotplot
  # ==========================================================

  enrichment_plot <- ggplot(

    plot_data,

    aes(

      x =
        GeneRatio_numeric,

      y =
        Description,

      size =
        Count,

      color =
        minus_log10_FDR

    )

  ) +

    geom_point(
      alpha = 0.85
    ) +

    scale_color_gradient(

      low =
        "#4DBBD5",

      high =
        "#E64B35"

    ) +

    theme_classic(
      base_size = 12
    ) +

    labs(

      title =
        plot_title,

      x =
        "Gene Ratio",

      y =
        NULL,

      size =
        "Gene Count",

      color =
        "-log10(FDR)"

    ) +

    theme(

      plot.title = element_text(

        face = "bold",

        size = 15

      ),

      axis.text.y = element_text(
        size = 10
      ),

      legend.position =
        "right"

    )



  # ==========================================================
  # Save dotplot
  # ==========================================================

  plot_png <- file.path(

    figure_dir,

    paste0(
      analysis_name,
      "_dotplot.png"
    )

  )


  plot_pdf <- file.path(

    figure_dir,

    paste0(
      analysis_name,
      "_dotplot.pdf"
    )

  )



  ggsave(

    filename =
      plot_png,

    plot =
      enrichment_plot,

    width = 8,

    height = 6,

    dpi = 300,

    bg = "white"

  )



  ggsave(

    filename =
      plot_pdf,

    plot =
      enrichment_plot,

    width = 8,

    height = 6,

    device = "pdf",

    bg = "white"

  )



  cat(
    "Dotplot saved.\n"
  )



  invisible(
    result_sig
  )

}



# ============================================================
# 12. GO Biological Process ORA
# ============================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("GO Biological Process ORA\n")
cat("------------------------------------------------------------\n")



# ============================================================
# 12.1 GO-BP UP
# ============================================================

cat("\nRunning GO-BP UP...\n")


go_up <- enrichGO(

  gene =
    up_entrez,

  universe =
    background_entrez,

  OrgDb =
    org.Hs.eg.db,

  keyType =
    "ENTREZID",

  ont =
    "BP",

  pAdjustMethod =
    "BH",

  # Return broad result first;
  # explicit FDR filtering is performed later.
  pvalueCutoff =
    1,

  qvalueCutoff =
    1,

  minGSSize =
    min_GS_size,

  maxGSSize =
    max_GS_size,

  readable =
    TRUE

)



save_enrichment_results(

  go_up,

  analysis_name =
    "GO_BP_UP",

  plot_title =
    "GO Biological Process - Upregulated Genes"

)



# ============================================================
# 12.2 GO-BP DOWN
# ============================================================

cat("\nRunning GO-BP DOWN...\n")


go_down <- enrichGO(

  gene =
    down_entrez,

  universe =
    background_entrez,

  OrgDb =
    org.Hs.eg.db,

  keyType =
    "ENTREZID",

  ont =
    "BP",

  pAdjustMethod =
    "BH",

  pvalueCutoff =
    1,

  qvalueCutoff =
    1,

  minGSSize =
    min_GS_size,

  maxGSSize =
    max_GS_size,

  readable =
    TRUE

)



save_enrichment_results(

  go_down,

  analysis_name =
    "GO_BP_DOWN",

  plot_title =
    "GO Biological Process - Downregulated Genes"

)



# ============================================================
# 13. KEGG ORA
# ============================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("KEGG Pathway ORA\n")
cat("------------------------------------------------------------\n")



# ============================================================
# 13.1 KEGG UP
# ============================================================

cat("\nRunning KEGG UP...\n")


kegg_up <- tryCatch(

  enrichKEGG(

    gene =
      up_entrez,

    organism =
      "hsa",

    keyType =
      "ncbi-geneid",

    universe =
      background_entrez,

    pAdjustMethod =
      "BH",

    pvalueCutoff =
      1,

    qvalueCutoff =
      1,

    minGSSize =
      min_GS_size,

    maxGSSize =
      max_GS_size,

    use_internal_data =
      FALSE

  ),

  error = function(e){

    cat(
      "\nKEGG UP failed:\n",
      conditionMessage(e),
      "\n"
    )

    NULL

  }

)



if(!is.null(kegg_up)){

  save_enrichment_results(

    kegg_up,

    analysis_name =
      "KEGG_UP",

    plot_title =
      "KEGG Pathways - Upregulated Genes"

  )

}



# ============================================================
# 13.2 KEGG DOWN
# ============================================================

cat("\nRunning KEGG DOWN...\n")


kegg_down <- tryCatch(

  enrichKEGG(

    gene =
      down_entrez,

    organism =
      "hsa",

    keyType =
      "ncbi-geneid",

    universe =
      background_entrez,

    pAdjustMethod =
      "BH",

    pvalueCutoff =
      1,

    qvalueCutoff =
      1,

    minGSSize =
      min_GS_size,

    maxGSSize =
      max_GS_size,

    use_internal_data =
      FALSE

  ),

  error = function(e){

    cat(
      "\nKEGG DOWN failed:\n",
      conditionMessage(e),
      "\n"
    )

    NULL

  }

)



if(!is.null(kegg_down)){

  save_enrichment_results(

    kegg_down,

    analysis_name =
      "KEGG_DOWN",

    plot_title =
      "KEGG Pathways - Downregulated Genes"

  )

}



# ============================================================
# 14. Summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Enrichment Summary\n")
cat("============================================================\n")


cat(
  "\nBackground genes:",
  length(background_entrez),
  "\n"
)


cat(
  "Up genes:",
  length(up_entrez),
  "\n"
)


cat(
  "Down genes:",
  length(down_entrez),
  "\n"
)



# ============================================================
# 15. Finished
# ============================================================

cat("\n")
cat("============================================================\n")
cat("08 Enrichment Analysis Finished\n")
cat("============================================================\n")