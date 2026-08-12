# ======================================
# 06_add_annotation.R
# Add gene annotation to DE results
# ======================================


library(dplyr)


# -------------------------
# 1. Paths
# -------------------------

annotation_file <- 
  "datasets/processed/U251_gene_annotation.tsv"


result_dir <- 
  "results/DE"


# -------------------------
# 2. Load annotation
# -------------------------

annotation <- read.table(
  annotation_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE
)
annotation <- annotation %>%
  rename(
    ENSEMBL = gene_id
  )

cat("Annotation loaded:",
    nrow(annotation),
    "genes\n")


# -------------------------
# Function
# -------------------------

add_annotation <- function(
    input_file,
    output_file
){

  df <- read.table(
    input_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )


  # check
 if(!"ENSEMBL" %in% colnames(df)){

  if("gene_id" %in% colnames(df)){

    df <- df %>%
      rename(
        ENSEMBL = gene_id
      )

  } else {

    stop(
      "Cannot find gene ID column in ",
      input_file
    )

  }

}

  df_anno <- df %>%
    left_join(
      annotation,
      by = "ENSEMBL"
    )

cat(
  "Annotation rate:",
  mean(!is.na(df_anno$SYMBOL)),
  "\n"
)

  write.table(
    df_anno,
    output_file,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )


  cat(
    basename(output_file),
    "completed:",
    nrow(df_anno),
    "rows\n"
  )

}



# -------------------------
# 3. Add annotation
# -------------------------


add_annotation(
  paste0(
    result_dir,
    "/U251_OE_vs_NC_limma_all_genes.tsv"
  ),

  paste0(
    result_dir,
    "/U251_OE_vs_NC_limma_all_genes_annotated.tsv"
  )
)



add_annotation(
  paste0(
    result_dir,
    "/U251_OE_vs_NC_TREAT_lfc1_all_genes.tsv"
  ),

  paste0(
    result_dir,
    "/U251_OE_vs_NC_TREAT_lfc1_all_genes_annotated.tsv"
  )
)



add_annotation(
  paste0(
    result_dir,
    "/U251_OE_vs_NC_TREAT_lfc1_FDR0.05_DEG.tsv"
  ),

  paste0(
    result_dir,
    "/U251_OE_vs_NC_TREAT_lfc1_FDR0.05_DEG_annotated.tsv"
  )
)


cat("\nAll annotation finished!\n")