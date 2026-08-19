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

required_fields <- c("input_file", "pathway_file", "output_prefix")
is_missing_field <- vapply(required_fields, function(field) {
  value <- cfg[[field]]
  is.null(value) || length(value) != 1L || is.na(value) ||
    !nzchar(as.character(value))
}, logical(1))
if (any(is_missing_field)) {
  stop(
    "Missing required plot config field(s): ",
    paste(required_fields[is_missing_field], collapse = ", "),
    call. = FALSE
  )
}

resolve_path <- function(path, must_work = FALSE) {
  path <- as.character(path)
  if (startsWith(path, "/")) {
    candidate <- path
  } else {
    # Config paths are interpreted from the directory where Rscript is run;
    # this supports the documented repository-root invocation without baking
    # any project-specific directory into the module.
    candidate <- file.path(getwd(), path)
  }
  normalizePath(candidate, mustWork = must_work)
}

input_file <- resolve_path(cfg$input_file, must_work = TRUE)
pathway_file <- resolve_path(cfg$pathway_file, must_work = TRUE)
output_prefix <- resolve_path(cfg$output_prefix, must_work = FALSE)

optional_config <- function(name, default) {
  value <- cfg[[name]]
  if (is.null(value) || length(value) == 0L) {
    return(default)
  }
  value
}

sort_by <- optional_config("sort_by", "NES")
color_by <- optional_config("color_by", "NES")
size_by <- optional_config("size_by", "setSize")
width <- as.numeric(optional_config("width", 7))
height <- as.numeric(optional_config("height", 5))
dpi <- as.numeric(optional_config("dpi", 300))

gsea_result <- read.delim(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_argument) != 1L) {
  stop("Unable to determine the path of run_plot.R.", call. = FALSE)
}
module_dir <- dirname(normalizePath(
  sub("^--file=", "", script_argument[[1]]),
  mustWork = TRUE
))
source(file.path(module_dir, "code_template.R"))

plot_hallmark_nes_summary(
  gsea_result = gsea_result,
  pathway_file = pathway_file,
  sort_by = sort_by,
  color_by = color_by,
  size_by = size_by,
  output_prefix = output_prefix,
  width = width,
  height = height,
  dpi = dpi
)
