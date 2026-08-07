#!/usr/bin/env Rscript

# The only supported command-line entry point for this workflow.

.script_path <- function() {
  file_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (length(file_arg) != 1L) {
    stop(
      "Unable to determine the path of run_spearman_analysis.R.",
      call. = FALSE
    )
  }
  normalizePath(
    sub("^--file=", "", file_arg),
    winslash = "/",
    mustWork = TRUE
  )
}

.exit_with_error <- function(message_text, status = 1L) {
  message(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    "] [ERROR] ",
    message_text
  )
  quit(save = "no", status = status, runLast = FALSE)
}

if (!requireNamespace("optparse", quietly = TRUE)) {
  .exit_with_error(
    paste0(
      "Required package 'optparse' is not installed. ",
      "Activate the spearman-r environment defined by environment.yml."
    ),
    status = 2L
  )
}

script_file <- tryCatch(
  .script_path(),
  error = function(e) .exit_with_error(conditionMessage(e), 2L)
)
workflow_dir <- dirname(script_file)
workflow_file <- file.path(
  workflow_dir,
  "R",
  "run_spearman_workflow.R"
)
workflow_env <- new.env(
  parent = asNamespace("stats")
)

source_error <- tryCatch(
  {
    source(
      workflow_file,
      local = workflow_env,
      encoding = "UTF-8"
    )
    NULL
  },
  error = function(e) e
)
if (inherits(source_error, "error")) {
  .exit_with_error(
    paste0(
      "Unable to load workflow modules: ",
      conditionMessage(source_error)
    ),
    2L
  )
}

parser <- optparse::OptionParser(
  usage = "%prog --config PATH",
  description = paste(
    "Run a configuration-driven Spearman correlation analysis.",
    "Relative paths in the YAML file are resolved from the YAML directory."
  ),
  option_list = list(
    optparse::make_option(
      c("-c", "--config"),
      type = "character",
      metavar = "PATH",
      help = "Path to a YAML configuration file."
    )
  )
)

options <- tryCatch(
  optparse::parse_args(parser),
  error = function(e) .exit_with_error(conditionMessage(e), 2L)
)

if (is.null(options$config) ||
    length(options$config) != 1L ||
    is.na(options$config) ||
    !nzchar(options$config)) {
  optparse::print_help(parser)
  .exit_with_error("--config PATH is required.", 2L)
}

config <- tryCatch(
  workflow_env$read_spearman_config(options$config),
  error = function(e) .exit_with_error(
    paste0("Configuration error: ", conditionMessage(e)),
    2L
  )
)

log_file <- file.path(
  config$output_dir,
  "logs",
  paste0(
    workflow_env$.spearman_safe_filename(config$project_id),
    ".log"
  )
)
logger <- tryCatch(
  workflow_env$.spearman_create_logger(
    log_file,
    config$logging$verbose
  ),
  error = function(e) .exit_with_error(conditionMessage(e), 1L)
)

logger("INFO", "Configuration: ", config$config_path)
logger("INFO", "Output directory: ", config$output_dir)
logger(
  "INFO",
  "R version: ",
  paste(R.version$major, R.version$minor, sep = ".")
)

status <- tryCatch(
  {
    result <- withCallingHandlers(
      workflow_env$run_spearman_project(config, logger),
      warning = function(w) {
        logger("WARN", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    logger(
      "INFO",
      "Completed with ",
      nrow(result$results),
      " result row(s)."
    )
    0L
  },
  error = function(e) {
    logger("ERROR", conditionMessage(e))
    1L
  }
)

quit(save = "no", status = status, runLast = FALSE)
