# ============================================================
# 07_visualization.R
#
# Visualization of differential expression results
#
# 1. Volcano plot
#    - Based on standard limma results
#    - DEG definition:
#         FDR < 0.05
#         |log2FC| >= 1
#
# 2. Top30 DEG heatmap
#    - Based on strict TREAT DEG results
#    - Top30 selected by |logFC|
#    - Row Z-score
#
#
# Input:
#
#   Volcano:
#   results/DE/
#   U251_OE_vs_NC_limma_all_genes_annotated.tsv
#
#   Heatmap:
#   results/DE/
#   U251_OE_vs_NC_TREAT_lfc1_FDR0.05_DEG_annotated.tsv
#
#   Expression matrix:
#   datasets/processed/
#   U251_expression_log2_filtered.rds
#
#
# Output:
#
#   figures/
#   U251_OE_vs_NC_volcano.png
#   U251_OE_vs_NC_volcano.pdf
#
#   U251_OE_vs_NC_top30_DEG_heatmap.png
#   U251_OE_vs_NC_top30_DEG_heatmap.pdf
#
# ============================================================



# ============================================================
# 1. Environment
# ============================================================

library(ggplot2)
library(dplyr)
library(ggrepel)
library(pheatmap)


cat("\n")
cat("============================================================\n")
cat("07 Visualization Started\n")
cat("============================================================\n")



# ============================================================
# 2. Parameters
# ============================================================

FDR_cutoff <- 0.05

logFC_cutoff <- 1

top_gene_n <- 30

label_gene_n <- 5



# ============================================================
# 3. Paths
# ============================================================

volcano_file <-
  "results/DE/U251_OE_vs_NC_limma_all_genes_annotated.tsv"


deg_file <-
  "results/DE/U251_OE_vs_NC_TREAT_lfc1_FDR0.05_DEG_annotated.tsv"


expr_file <-
  "datasets/processed/U251_expression_log2_filtered.rds"


figure_dir <-
  "figures"



# ============================================================
# 4. Create output directory
# ============================================================

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat("\nOutput directory:\n")
cat(figure_dir, "\n")



# ============================================================
# 5. Check input files
# ============================================================

input_files <- c(
  volcano_file,
  deg_file,
  expr_file
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
# 6. Volcano plot
# ============================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("Generating volcano plot\n")
cat("------------------------------------------------------------\n")



# ------------------------------------------------------------
# 6.1 Load standard limma results
# ------------------------------------------------------------

volcano_data <- read.table(
  volcano_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE
)


cat(
  "\nTotal limma genes:",
  nrow(volcano_data),
  "\n"
)



# ------------------------------------------------------------
# 6.2 Check required columns
# ------------------------------------------------------------

required_volcano_columns <- c(
  "logFC",
  "adj.P.Val",
  "SYMBOL"
)


missing_columns <- setdiff(
  required_volcano_columns,
  colnames(volcano_data)
)


if(length(missing_columns) > 0){

  stop(
    paste(
      "Missing volcano column(s):",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )

}



# ------------------------------------------------------------
# 6.3 Remove invalid values
# ------------------------------------------------------------

volcano_data <- volcano_data %>%

  filter(
    !is.na(logFC),
    !is.na(adj.P.Val)
  )



# ------------------------------------------------------------
# 6.4 Define DEG groups
# ------------------------------------------------------------

volcano_data <- volcano_data %>%

  mutate(

    change = case_when(

      adj.P.Val < FDR_cutoff &
        logFC >= logFC_cutoff
      ~ "Up",

      adj.P.Val < FDR_cutoff &
        logFC <= -logFC_cutoff
      ~ "Down",

      TRUE
      ~ "NS"

    ),

    minus_log10_FDR =
      -log10(
        pmax(
          adj.P.Val,
          .Machine$double.xmin
        )
      )

  )



cat("\nVolcano classification:\n")

print(
  table(volcano_data$change)
)



# ------------------------------------------------------------
# 6.5 Remove FDR = 1 only for visualization
#
# FDR = 1:
#
# -log10(1) = 0
#
# These points form a horizontal pile at y = 0.
# Removing them here does NOT alter statistical results.
# ------------------------------------------------------------

volcano_plot_data <- volcano_data %>%

  filter(
    adj.P.Val < 1
  )


cat(
  "\nGenes shown in volcano:",
  nrow(volcano_plot_data),
  "\n"
)



# ------------------------------------------------------------
# 6.6 Select labels
#
# Top 5 Up + Top 5 Down by FDR
# ------------------------------------------------------------

label_up <- volcano_plot_data %>%

  filter(
    change == "Up",
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%

  arrange(
    adj.P.Val
  ) %>%

  slice_head(
    n = label_gene_n
  )



label_down <- volcano_plot_data %>%

  filter(
    change == "Down",
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%

  arrange(
    adj.P.Val
  ) %>%

  slice_head(
    n = label_gene_n
  )



label_genes <- bind_rows(
  label_up,
  label_down
)



cat("\nGenes labeled on volcano:\n")

print(
  label_genes[
    ,
    c(
      "SYMBOL",
      "logFC",
      "adj.P.Val",
      "change"
    )
  ]
)



# ------------------------------------------------------------
# 6.7 Build volcano plot
# ------------------------------------------------------------

volcano_plot <- ggplot()



# NS genes
volcano_plot <- volcano_plot +

  geom_point(

    data =
      volcano_plot_data %>%
      filter(change == "NS"),

    aes(
      x = logFC,
      y = minus_log10_FDR
    ),

    color = "grey75",

    size = 0.9,

    alpha = 0.45

  )



# Down-regulated genes
volcano_plot <- volcano_plot +

  geom_point(

    data =
      volcano_plot_data %>%
      filter(change == "Down"),

    aes(
      x = logFC,
      y = minus_log10_FDR,
      color = change
    ),

    size = 1.6,

    alpha = 0.80

  )



# Up-regulated genes
volcano_plot <- volcano_plot +

  geom_point(

    data =
      volcano_plot_data %>%
      filter(change == "Up"),

    aes(
      x = logFC,
      y = minus_log10_FDR,
      color = change
    ),

    size = 1.6,

    alpha = 0.80

  )



# Threshold lines
volcano_plot <- volcano_plot +

  geom_vline(

    xintercept = c(
      -logFC_cutoff,
      logFC_cutoff
    ),

    linetype = "dashed",

    linewidth = 0.5,

    color = "grey40"

  ) +

  geom_hline(

    yintercept =
      -log10(FDR_cutoff),

    linetype = "dashed",

    linewidth = 0.5,

    color = "grey40"

  )



# Gene labels
volcano_plot <- volcano_plot +

  geom_text_repel(

    data = label_genes,

    aes(
      x = logFC,
      y = minus_log10_FDR,
      label = SYMBOL
    ),

    size = 3.5,

    color = "black",

    box.padding = 0.5,

    point.padding = 0.3,

    min.segment.length = 0,

    segment.color = "grey50",

    max.overlaps = Inf,

    seed = 123

  )



# Colors
volcano_plot <- volcano_plot +

  scale_color_manual(

    values = c(

      "Down" = "#377EB8",

      "Up" = "#E41A1C"

    ),

    breaks = c(
      "Down",
      "Up"
    )

  )



# Theme
volcano_plot <- volcano_plot +

  theme_classic(
    base_size = 14
  ) +

  labs(

    title =
      "U251 OE vs NC",

    subtitle =
      "limma | FDR < 0.05 | |log2FC| >= 1",

    x =
      "log2 Fold Change",

    y =
      "-log10(FDR)",

    color =
      NULL

  ) +

  theme(

    plot.title = element_text(
      size = 18,
      face = "bold"
    ),

    plot.subtitle = element_text(
      size = 12
    ),

    axis.title = element_text(
      size = 14
    ),

    axis.text = element_text(
      size = 11
    ),

    legend.position = "right"

  )



# ------------------------------------------------------------
# 6.8 Save volcano
#
# PNG = quick preview
# PDF = publication / vector master
# ------------------------------------------------------------

volcano_png <- file.path(
  figure_dir,
  "U251_OE_vs_NC_volcano.png"
)


volcano_pdf <- file.path(
  figure_dir,
  "U251_OE_vs_NC_volcano.pdf"
)



ggsave(

  filename = volcano_png,

  plot = volcano_plot,

  width = 7,

  height = 6,

  dpi = 300,

  bg = "white"

)



ggsave(

  filename = volcano_pdf,

  plot = volcano_plot,

  width = 7,

  height = 6,

  device = "pdf",

  bg = "white"

)



cat(
  "\nVolcano saved:\n",
  volcano_png,
  "\n",
  volcano_pdf,
  "\n"
)



# ============================================================
# 7. Top30 TREAT DEG heatmap
# ============================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("Generating Top30 DEG heatmap\n")
cat("------------------------------------------------------------\n")



# ------------------------------------------------------------
# 7.1 Load strict TREAT DEG results
# ------------------------------------------------------------

deg <- read.table(

  deg_file,

  header = TRUE,

  sep = "\t",

  quote = "",

  stringsAsFactors = FALSE

)


cat(
  "\nStrict TREAT DEG number:",
  nrow(deg),
  "\n"
)



# ------------------------------------------------------------
# 7.2 Check required columns
# ------------------------------------------------------------

required_deg_columns <- c(
  "ENSEMBL",
  "SYMBOL",
  "logFC",
  "adj.P.Val"
)


missing_columns <- setdiff(
  required_deg_columns,
  colnames(deg)
)


if(length(missing_columns) > 0){

  stop(
    paste(
      "Missing DEG column(s):",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )

}



# ------------------------------------------------------------
# 7.3 Load expression matrix
# ------------------------------------------------------------

expr <- readRDS(
  expr_file
)


cat(
  "\nExpression matrix:",
  nrow(expr),
  "genes x",
  ncol(expr),
  "samples\n"
)



if(is.null(rownames(expr))){

  stop(
    "Expression matrix has no gene IDs in rownames."
  )

}



# ------------------------------------------------------------
# 7.4 Keep DEG present in expression matrix
# ------------------------------------------------------------

deg_heatmap <- deg %>%

  filter(
    !is.na(ENSEMBL),
    ENSEMBL %in% rownames(expr)
  )


cat(
  "\nDEGs present in expression matrix:",
  nrow(deg_heatmap),
  "\n"
)



# ------------------------------------------------------------
# 7.5 Select Top30 by absolute logFC
# ------------------------------------------------------------

deg_heatmap <- deg_heatmap %>%

  arrange(
    desc(abs(logFC))
  )


top_gene_n_actual <- min(
  top_gene_n,
  nrow(deg_heatmap)
)


top30 <- deg_heatmap %>%

  slice_head(
    n = top_gene_n_actual
  )


cat(
  "\nTop genes selected:",
  top_gene_n_actual,
  "\n"
)



print(

  top30[
    ,
    c(
      "ENSEMBL",
      "SYMBOL",
      "logFC",
      "adj.P.Val"
    )
  ]

)



# ------------------------------------------------------------
# 7.6 Extract expression
# ------------------------------------------------------------

heat_expr <- expr[
  top30$ENSEMBL,
  ,
  drop = FALSE
]



# ------------------------------------------------------------
# 7.7 Gene labels
#
# SYMBOL is used for display.
# ENSEMBL is used if SYMBOL is unavailable.
# ------------------------------------------------------------

gene_labels <- ifelse(

  is.na(top30$SYMBOL) |
    top30$SYMBOL == "",

  top30$ENSEMBL,

  top30$SYMBOL

)


gene_labels <- make.unique(
  gene_labels
)


rownames(heat_expr) <-
  gene_labels



# ------------------------------------------------------------
# 7.8 Remove zero-variance genes
# ------------------------------------------------------------

gene_sd <- apply(

  heat_expr,

  1,

  sd,

  na.rm = TRUE

)


keep_gene <- is.finite(gene_sd) &
  gene_sd > 0


heat_expr <- heat_expr[
  keep_gene,
  ,
  drop = FALSE
]


cat(
  "\nGenes retained after SD check:",
  nrow(heat_expr),
  "\n"
)



# ------------------------------------------------------------
# 7.9 Row Z-score
# ------------------------------------------------------------

heat_z <- t(
  scale(
    t(heat_expr)
  )
)



# ------------------------------------------------------------
# 7.10 Sample annotation
# ------------------------------------------------------------

sample_group <- case_when(

  grepl(
    "^OE",
    colnames(heat_z)
  ) ~ "OE",

  grepl(
    "^NC",
    colnames(heat_z)
  ) ~ "NC",

  TRUE ~ NA_character_

)



if(any(is.na(sample_group))){

  stop(
    paste(
      "Cannot determine OE/NC group for sample(s):",
      paste(
        colnames(heat_z)[
          is.na(sample_group)
        ],
        collapse = ", "
      )
    )
  )

}



annotation_col <- data.frame(

  Group = factor(
    sample_group,
    levels = c(
      "NC",
      "OE"
    )
  )

)


rownames(annotation_col) <-
  colnames(heat_z)



annotation_colors <- list(

  Group = c(

    "NC" = "#4DBBD5",

    "OE" = "#F28E85"

  )

)



# ------------------------------------------------------------
# 7.11 Symmetric Z-score color scale
# ------------------------------------------------------------

z_limit <- max(
  abs(heat_z),
  na.rm = TRUE
)


heat_breaks <- seq(
  -z_limit,
  z_limit,
  length.out = 101
)


heat_colors <- colorRampPalette(
  c(
    "#3B6FB6",
    "white",
    "#D73027"
  )
)(100)



# ------------------------------------------------------------
# 7.12 Heatmap output paths
# ------------------------------------------------------------

heatmap_png <- file.path(
  figure_dir,
  "U251_OE_vs_NC_top30_DEG_heatmap.png"
)


heatmap_pdf <- file.path(
  figure_dir,
  "U251_OE_vs_NC_top30_DEG_heatmap.pdf"
)



# ------------------------------------------------------------
# 7.13 Heatmap plotting function
#
# Same visual parameters are used for PNG and PDF.
# ------------------------------------------------------------

save_heatmap <- function(output_file){

  pheatmap(

    heat_z,

    # clustering
    cluster_rows = TRUE,
    cluster_cols = TRUE,

    # sample annotation
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,

    # color scale
    color = heat_colors,
    breaks = heat_breaks,

    # labels
    show_rownames = TRUE,
    show_colnames = TRUE,

    fontsize = 10,
    fontsize_row = 9,
    fontsize_col = 10,

    angle_col = 45,

    # appearance
    border_color = NA,

    treeheight_row = 90,
    treeheight_col = 50,

    cellwidth = 45,
    cellheight = 18,

    # title
    main =
      "Top 30 Differentially Expressed Genes",

    # output
    filename =
      output_file,

    width = 8,
    height = 10

  )

}



# ------------------------------------------------------------
# 7.14 Save heatmap
# ------------------------------------------------------------

save_heatmap(
  heatmap_png
)


save_heatmap(
  heatmap_pdf
)



cat(
  "\nHeatmap saved:\n",
  heatmap_png,
  "\n",
  heatmap_pdf,
  "\n"
)



# ============================================================
# 8. Finished
# ============================================================

cat("\n")
cat("============================================================\n")
cat("07 Visualization Finished\n")
cat("============================================================\n")