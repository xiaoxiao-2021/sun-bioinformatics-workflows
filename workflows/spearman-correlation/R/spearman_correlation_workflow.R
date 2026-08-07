# Compatibility loader for the former monolithic script.
#
# This file no longer starts an analysis or reads variables from the Global
# Environment. Source R/run_spearman_workflow.R directly for the function API,
# or use run_spearman_analysis.R --config PATH for command-line analyses.

.legacy_loader_file <- NULL
.legacy_loader_frames <- sys.frames()
.legacy_loader_i <- NULL
.legacy_loader_candidate <- NULL
if (length(.legacy_loader_frames) > 0L) {
  for (.legacy_loader_i in rev(seq_along(.legacy_loader_frames))) {
    .legacy_loader_candidate <- .legacy_loader_frames[[.legacy_loader_i]]$ofile
    if (is.null(.legacy_loader_candidate) ||
        !is.character(.legacy_loader_candidate) ||
        length(.legacy_loader_candidate) != 1L ||
        is.na(.legacy_loader_candidate) ||
        !nzchar(.legacy_loader_candidate)) {
      next
    }
    .legacy_loader_candidate <- tryCatch(
      normalizePath(
        .legacy_loader_candidate,
        winslash = "/",
        mustWork = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(.legacy_loader_candidate) &&
        identical(
          basename(.legacy_loader_candidate),
          "spearman_correlation_workflow.R"
        )) {
      .legacy_loader_file <- .legacy_loader_candidate
      break
    }
  }
}

if (is.null(.legacy_loader_file)) {
  stop(
    "Unable to locate R/spearman_correlation_workflow.R.",
    call. = FALSE
  )
}

source(
  file.path(
    dirname(.legacy_loader_file),
    "run_spearman_workflow.R"
  ),
  local = environment(),
  encoding = "UTF-8"
)

message(
  "Loaded compatibility functions. For formal runs use: ",
  "Rscript run_spearman_analysis.R --config PATH"
)

rm(
  .legacy_loader_file,
  .legacy_loader_frames,
  .legacy_loader_i,
  .legacy_loader_candidate
)
