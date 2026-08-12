# ============================================================
# 09_enrichment_visualization.R
#
# Visualization of ORA enrichment results
#
# Input:
#
#   results/enrichment/
#   ├── GO_BP_UP_FDR0.05.tsv
#   ├── GO_BP_DOWN_FDR0.05.tsv
#   ├── KEGG_UP_FDR0.05.tsv
#   └── KEGG_DOWN_FDR0.05.tsv
#
#
# Output:
#
#   figures/enrichment/
#   ├── GO_BP_UP_dotplot.png
#   ├── GO_BP_UP_dotplot.pdf
#   ├── GO_BP_DOWN_dotplot.png
#   ├── GO_BP_DOWN_dotplot.pdf
#   ├── KEGG_UP_dotplot.png
#   ├── KEGG_UP_dotplot.pdf
#   ├── KEGG_DOWN_dotplot.png
#   └── KEGG_DOWN_dotplot.pdf
#
#
#   results/enrichment/plot_data/
#   ├── GO_BP_UP_plot_terms.tsv
#   ├── GO_BP_DOWN_plot_terms.tsv
#   ├── KEGG_UP_plot_terms.tsv
#   └── KEGG_DOWN_plot_terms.tsv
#
# ============================================================



# ============================================================
# 1. Environment
# ============================================================

library(ggplot2)
library(dplyr)


cat("\n")
cat("============================================================\n")
cat("09 Enrichment Visualization Started\n")
cat("============================================================\n")



# ============================================================
# 2. Parameters
# ============================================================

show_category_n <- 15

FDR_cutoff <- 0.05



# ============================================================
# 3. Input / output paths
# ============================================================

enrichment_dir <-
  "results/enrichment"


figure_dir <-
  "figures/enrichment"


plot_data_dir <-
  "results/enrichment/plot_data"



go_up_file <- file.path(
  enrichment_dir,
  "GO_BP_UP_FDR0.05.tsv"
)


go_down_file <- file.path(
  enrichment_dir,
  "GO_BP_DOWN_FDR0.05.tsv"
)


kegg_up_file <- file.path(
  enrichment_dir,
  "KEGG_UP_FDR0.05.tsv"
)


kegg_down_file <- file.path(
  enrichment_dir,
  "KEGG_DOWN_FDR0.05.tsv"
)



# ============================================================
# 4. Create output directories automatically
# ============================================================

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  plot_data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



cat("\nOutput directories checked:\n")

cat(
  "Figures:",
  figure_dir,
  "\n"
)

cat(
  "Plot data:",
  plot_data_dir,
  "\n"
)



# ============================================================
# 5. Check input files
# ============================================================

input_files <- c(
  go_up_file,
  go_down_file,
  kegg_up_file,
  kegg_down_file
)


missing_files <- input_files[
  !file.exists(input_files)
]


if(length(missing_files) > 0){

  stop(
    paste(
      "Missing enrichment file(s):\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )

}


cat("\nAll enrichment files detected.\n")



# ============================================================
# 6. Helper function:
#    GeneRatio string -> numeric
#
# Example:
#
# 13/68
#   ↓
# 0.191
# ============================================================

gene_ratio_to_numeric <- function(x){

  sapply(

    strsplit(
      x,
      "/",
      fixed = TRUE
    ),

    function(z){

      if(length(z) != 2){
        return(NA_real_)
      }

      numerator <-
        suppressWarnings(
          as.numeric(z[1])
        )

      denominator <-
        suppressWarnings(
          as.numeric(z[2])
        )


      if(
        is.na(numerator) |
        is.na(denominator) |
        denominator == 0
      ){
        return(NA_real_)
      }


      numerator / denominator

    }

  )

}



# ============================================================
# 7. Helper function:
#    detect FDR column
#
# Different software/file viewers may display:
#
# p.adjust
# P.Adjust
# padj
#
# This function makes the workflow more robust.
# ============================================================

find_fdr_column <- function(df){

  possible_names <- c(
    "p.adjust",
    "P.Adjust",
    "p_adjust",
    "padj",
    "FDR"
  )


  found <- possible_names[
    possible_names %in% colnames(df)
  ]


  if(length(found) == 0){

    stop(
      paste(
        "Cannot find FDR column.\n",
        "Available columns:\n",
        paste(
          colnames(df),
          collapse = ", "
        )
      )
    )

  }


  found[1]

}



# ============================================================
# 8. Main plotting function
# ============================================================

make_enrichment_dotplot <- function(

    input_file,

    analysis_name,

    plot_title

){


  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("Processing:", analysis_name, "\n")
  cat("------------------------------------------------------------\n")



  # ----------------------------------------------------------
  # 8.1 Load enrichment table
  # ----------------------------------------------------------

  enrichment <- read.delim(

    input_file,

    header = TRUE,

    sep = "\t",

    quote = "",

    stringsAsFactors = FALSE,

    check.names = FALSE

  )



  cat(
    "Significant terms loaded:",
    nrow(enrichment),
    "\n"
  )



  # ----------------------------------------------------------
  # 8.2 Empty result protection
  # ----------------------------------------------------------

  if(nrow(enrichment) == 0){

    cat(
      "No significant enrichment terms. Plot skipped.\n"
    )

    return(
      invisible(NULL)
    )

  }



  # ----------------------------------------------------------
  # 8.3 Check required columns
  # ----------------------------------------------------------

  required_columns <- c(
    "Description",
    "GeneRatio",
    "Count"
  )


  missing_columns <- setdiff(
    required_columns,
    colnames(enrichment)
  )


  if(length(missing_columns) > 0){

    stop(
      paste(
        analysis_name,
        "is missing column(s):",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )

  }



  # ----------------------------------------------------------
  # 8.4 Detect FDR column
  # ----------------------------------------------------------

  fdr_col <- find_fdr_column(
    enrichment
  )


  cat(
    "FDR column detected:",
    fdr_col,
    "\n"
  )



  # ----------------------------------------------------------
  # 8.5 Prepare plotting data
  # ----------------------------------------------------------

  enrichment <- enrichment %>%

    mutate(

      FDR =
        as.numeric(
          .data[[fdr_col]]
        ),

      GeneRatio_numeric =
        gene_ratio_to_numeric(
          GeneRatio
        ),

      minus_log10_FDR =
        -log10(
          pmax(
            FDR,
            .Machine$double.xmin
          )
        )

    ) %>%

    filter(

      !is.na(FDR),

      !is.na(GeneRatio_numeric),

      !is.na(Count),

      FDR < FDR_cutoff

    )



  # ----------------------------------------------------------
  # 8.6 Select top enrichment terms
  #
  # If there are fewer than 15 terms,
  # all significant terms are displayed.
  # ----------------------------------------------------------

  plot_n <- min(
    show_category_n,
    nrow(enrichment)
  )


  plot_data <- enrichment %>%

    arrange(
      FDR
    ) %>%

    slice_head(
      n = plot_n
    )



  cat(
    "Terms displayed:",
    nrow(plot_data),
    "\n"
  )



  # ----------------------------------------------------------
  # 8.7 Save plotting data
  # ----------------------------------------------------------

  plot_data_file <- file.path(

    plot_data_dir,

    paste0(
      analysis_name,
      "_plot_terms.tsv"
    )

  )


  write.table(

    plot_data,

    plot_data_file,

    sep = "\t",

    quote = FALSE,

    row.names = FALSE

  )



  cat(
    "Plot data saved:",
    plot_data_file,
    "\n"
  )



  # ----------------------------------------------------------
  # 8.8 Set pathway display order
  #
  # Most significant terms appear at the top.
  # ----------------------------------------------------------

  plot_data$Description <- factor(

    plot_data$Description,

    levels =
      rev(
        plot_data$Description
      )

  )



  # ----------------------------------------------------------
  # 8.9 Dynamic figure height
  #
  # GO-UP may have 15 terms,
  # while KEGG-DOWN may have only 2.
  #
  # Avoid making a 2-term plot excessively tall.
  # ----------------------------------------------------------

  plot_height <- max(
    3.5,
    2.5 + 0.30 * nrow(plot_data)
  )



  # ----------------------------------------------------------
  # 8.10 Dotplot
  # ----------------------------------------------------------

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

    scale_x_continuous(

      expand =
        expansion(
          mult = c(
            0.02,
            0.08
          )
        )

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

      plot.title =
        element_text(
          size = 15,
          face = "bold",
          hjust = 0
        ),

      axis.title.x =
        element_text(
          size = 12
        ),

      axis.text.x =
        element_text(
          size = 10
        ),

      axis.text.y =
        element_text(
          size = 10,
          color = "black"
        ),

      legend.title =
        element_text(
          size = 10
        ),

      legend.text =
        element_text(
          size = 9
        ),

      legend.position =
        "right",

      plot.margin =
        margin(
          t = 10,
          r = 15,
          b = 10,
          l = 10
        )

    )



  # ----------------------------------------------------------
  # 8.11 Output paths
  # ----------------------------------------------------------

  png_file <- file.path(

    figure_dir,

    paste0(
      analysis_name,
      "_dotplot.png"
    )

  )


  pdf_file <- file.path(

    figure_dir,

    paste0(
      analysis_name,
      "_dotplot.pdf"
    )

  )



  # ----------------------------------------------------------
  # 8.12 Save PNG
  # ----------------------------------------------------------

  ggsave(

    filename =
      png_file,

    plot =
      enrichment_plot,

    width = 9,

    height = plot_height,

    dpi = 300,

    bg = "white"

  )



  # ----------------------------------------------------------
  # 8.13 Save PDF
  # ----------------------------------------------------------

  ggsave(

    filename =
      pdf_file,

    plot =
      enrichment_plot,

    width = 9,

    height = plot_height,

    device = "pdf",

    bg = "white"

  )



  cat(
    "PNG saved:",
    png_file,
    "\n"
  )


  cat(
    "PDF saved:",
    pdf_file,
    "\n"
  )



  invisible(
    enrichment_plot
  )

}



# ============================================================
# 9. GO Biological Process - UP
# ============================================================

make_enrichment_dotplot(

  input_file =
    go_up_file,

  analysis_name =
    "GO_BP_UP",

  plot_title =
    "GO Biological Process - Upregulated Genes"

)



# ============================================================
# 10. GO Biological Process - DOWN
# ============================================================

make_enrichment_dotplot(

  input_file =
    go_down_file,

  analysis_name =
    "GO_BP_DOWN",

  plot_title =
    "GO Biological Process - Downregulated Genes"

)



# ============================================================
# 11. KEGG - UP
# ============================================================

make_enrichment_dotplot(

  input_file =
    kegg_up_file,

  analysis_name =
    "KEGG_UP",

  plot_title =
    "KEGG Pathways - Upregulated Genes"

)



# ============================================================
# 12. KEGG - DOWN
# ============================================================

make_enrichment_dotplot(

  input_file =
    kegg_down_file,

  analysis_name =
    "KEGG_DOWN",

  plot_title =
    "KEGG Pathways - Downregulated Genes"

)



# ============================================================
# 13. Finished
# ============================================================

cat("\n")
cat("============================================================\n")
cat("09 Enrichment Visualization Finished\n")
cat("============================================================\n")


cat("\nFigures saved in:\n")

cat(
  figure_dir,
  "\n"
)


cat("\nPlot data saved in:\n")

cat(
  plot_data_dir,
  "\n"
)