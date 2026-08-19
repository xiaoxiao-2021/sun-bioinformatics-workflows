#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[[1]])) {
  stop(
    "Provide exactly one plot config YAML file as the command-line argument.",
    call. = FALSE
  )
}
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Missing required R package: yaml", call. = FALSE)
}

config_file <- normalizePath(args[[1]], mustWork = TRUE)
cfg <- yaml::read_yaml(config_file)
if (!is.list(cfg)) {
  stop("The plot config YAML must contain a named mapping.", call. = FALSE)
}

required_fields <- c(
  "gsea_result_file", "ranked_gene_file", "term2gene_file", "pathway_file",
  "gene_id_column", "rank_column", "case_group", "control_group", "output_dir"
)
missing_field <- vapply(required_fields, function(field) {
  value <- cfg[[field]]
  is.null(value) || length(value) != 1L || is.na(value) ||
    !nzchar(as.character(value))
}, logical(1))
if (any(missing_field)) {
  stop(
    "Missing required plot config field(s): ",
    paste(required_fields[missing_field], collapse = ", "),
    call. = FALSE
  )
}

resolve_path <- function(path, must_work = FALSE) {
  path <- as.character(path)
  candidate <- if (startsWith(path, "/")) path else file.path(getwd(), path)
  normalizePath(candidate, mustWork = must_work)
}

gsea_result_file <- resolve_path(cfg$gsea_result_file, must_work = TRUE)
ranked_gene_file <- resolve_path(cfg$ranked_gene_file, must_work = TRUE)
term2gene_file <- resolve_path(cfg$term2gene_file, must_work = TRUE)
pathway_file <- resolve_path(cfg$pathway_file, must_work = TRUE)
output_dir <- resolve_path(cfg$output_dir, must_work = FALSE)

optional_config <- function(name, default) {
  value <- cfg[[name]]
  if (is.null(value) || length(value) == 0L) default else value
}

gene_id_column <- as.character(cfg$gene_id_column)
rank_column <- as.character(cfg$rank_column)
case_group <- as.character(cfg$case_group)
control_group <- as.character(cfg$control_group)
gsea_exponent <- as.numeric(optional_config("gsea_exponent", 1))
width <- as.numeric(optional_config("width", 7))
height <- as.numeric(optional_config("height", 5))
dpi <- as.numeric(optional_config("dpi", 300))

gsea_result <- read.delim(
  gsea_result_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
ranked_genes <- read.delim(
  ranked_gene_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
term2gene <- read.delim(
  term2gene_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!all(c("ID", "Description", "NES", "pvalue") %in% names(gsea_result))) {
  stop("gsea_result_file must contain ID, Description, NES, and pvalue.", call. = FALSE)
}
if (!all(c("TERM", "GENE") %in% names(term2gene))) {
  stop("term2gene_file must contain TERM and GENE.", call. = FALSE)
}
if (!gene_id_column %in% names(ranked_genes) || !rank_column %in% names(ranked_genes)) {
  stop("ranked_gene_file is missing gene_id_column or rank_column.", call. = FALSE)
}

# The ranked file is a snapshot of the gene-level ranking used by GSEA.  Do
# not let the per-pathway error handler turn duplicated IDs into a partial
# visualization; duplicated IDs invalidate the ranking for every pathway.
ranked_ids_for_check <- as.character(ranked_genes[[gene_id_column]])
ranked_values_for_check <- suppressWarnings(
  as.numeric(as.character(ranked_genes[[rank_column]]))
)
valid_ranked_rows <- !is.na(ranked_ids_for_check) &
  nzchar(ranked_ids_for_check) & is.finite(ranked_values_for_check)
if (anyDuplicated(ranked_ids_for_check[valid_ranked_rows])) {
  stop(
    "ranked_gene_file contains duplicate gene IDs; provide the gene-level ranking snapshot used by GSEA.",
    call. = FALSE
  )
}

pathway_config <- yaml::read_yaml(pathway_file)
if (is.null(pathway_config$pathway_set) ||
    !is.character(pathway_config$pathway_set) ||
    length(pathway_config$pathway_set) == 0L) {
  stop('The pathway YAML must define a non-empty character vector named "pathway_set".', call. = FALSE)
}
requested_pathways <- unique(as.character(pathway_config$pathway_set))
gsea_ids <- unique(as.character(gsea_result$ID))
term_ids <- unique(as.character(term2gene$TERM))
matched_gsea <- intersect(requested_pathways, gsea_ids)
matched_term2gene <- intersect(requested_pathways, term_ids)
plot_pathways <- requested_pathways[
  requested_pathways %in% matched_gsea & requested_pathways %in% matched_term2gene
]
missing_gsea <- setdiff(requested_pathways, matched_gsea)
missing_term2gene <- setdiff(requested_pathways, matched_term2gene)

cat("requested pathway count: ", length(requested_pathways), "\n", sep = "")
cat("matched in GSEA result: ", length(matched_gsea), "\n", sep = "")
cat("matched in TERM2GENE: ", length(matched_term2gene), "\n", sep = "")
cat("missing pathways: ", length(union(missing_gsea, missing_term2gene)), "\n", sep = "")
if (length(missing_gsea)) {
  warning(
    "Pathways missing from GSEA result: ", paste(missing_gsea, collapse = ", "),
    call. = FALSE
  )
}
if (length(missing_term2gene)) {
  warning(
    "Pathways missing from TERM2GENE: ", paste(missing_term2gene, collapse = ", "),
    call. = FALSE
  )
}
if (!length(plot_pathways)) {
  stop("No selected pathways are present in both the GSEA result and TERM2GENE snapshot.", call. = FALSE)
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("Unable to determine the path of run_plot.R.", call. = FALSE)
}
module_dir <- dirname(normalizePath(
  sub("^--file=", "", script_argument[[1]]),
  mustWork = TRUE
))
source(file.path(module_dir, "code_template.R"))

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(output_dir)) {
  stop("Unable to create output directory: ", output_dir, call. = FALSE)
}

summary_rows <- list()
skipped_pathways <- character()
skipped_reasons <- character()
for (pathway_id in plot_pathways) {
  result <- tryCatch(
    plot_gsea_enrichment_curve(
      pathway_id = pathway_id,
      gsea_result = gsea_result,
      ranked_genes = ranked_genes,
      term2gene = term2gene,
      gene_id_column = gene_id_column,
      rank_column = rank_column,
      case_group = case_group,
      control_group = control_group,
      gsea_exponent = gsea_exponent,
      output_dir = output_dir,
      width = width,
      height = height,
      dpi = dpi
    ),
    error = function(error) {
      skipped_reasons[pathway_id] <<- conditionMessage(error)
      warning(
        "Skipping pathway ", pathway_id, ": ", conditionMessage(error),
        call. = FALSE
      )
      NULL
    }
  )
  if (is.null(result)) {
    skipped_pathways <- c(skipped_pathways, pathway_id)
  } else {
    summary_row <- result$summary
    summary_row$status <- "plotted"
    summary_row$skip_reason <- NA_character_
    summary_rows[[length(summary_rows) + 1L]] <- summary_row
  }
}

if (length(summary_rows)) {
  summary_table <- do.call(rbind, summary_rows)
} else {
  stop("No selected pathway could be plotted after ranked-gene matching.", call. = FALSE)
}

# Keep one audit row for every requested pathway.  Missing or skipped
# pathways retain an explicit status while curve quantities remain NA.
for (pathway_id in setdiff(requested_pathways, as.character(summary_table$ID))) {
  reason <- skipped_reasons[[pathway_id]]
  if (is.null(reason) || is.na(reason) || !nzchar(reason)) {
    reason <- if (pathway_id %in% missing_gsea) {
      "Missing from GSEA result"
    } else if (pathway_id %in% missing_term2gene) {
      "Missing from TERM2GENE snapshot"
    } else {
      "Not plotted"
    }
  }
  summary_table <- rbind(
    summary_table,
    data.frame(
      ID = pathway_id,
      Description = NA_character_,
      NES = NA_real_,
      pvalue = NA_real_,
      p.adjust = NA_real_,
      qvalue = NA_real_,
      setSize = NA_real_,
      matched_gene_count = NA_integer_,
      visual_ES = NA_real_,
      visual_ES_rank = NA_integer_,
      enrichmentScore = NA_real_,
      visual_ES_difference = NA_real_,
      status = "skipped",
      skip_reason = reason,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
}
utils::write.table(
  summary_table,
  file = file.path(output_dir, "gsea_enrichment_curve_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

cat("Requested pathways: ", length(requested_pathways), "\n", sep = "")
cat("Successfully plotted: ", sum(summary_table$status == "plotted"), "\n", sep = "")
cat("Skipped: ", length(skipped_pathways), "\n", sep = "")
if (length(skipped_pathways)) {
  cat(paste(skipped_pathways, collapse = ", "), "\n", sep = "")
}
