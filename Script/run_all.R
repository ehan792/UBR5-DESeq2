################################################################################
# run_all.R — Strict runner for refactored UBR5 pipeline
################################################################################

# Locate and source config from the same directory as this script when possible.
this_file_candidates <- c(
  Sys.getenv("UBR5_CONFIG", unset = ""),
  file.path(getwd(), "00_config.R"),
  file.path(dirname(getwd()), "00_config.R")
)
this_file_candidates <- this_file_candidates[file.exists(this_file_candidates)]
if (!length(this_file_candidates)) {
  hits <- list.files(getwd(), pattern = "^00_config\\.R$", recursive = TRUE, full.names = TRUE)
  hits <- hits[!grepl("/\\.git/|/results/|/renv/|/\\.Rproj\\.user/", hits)]
  if (!length(hits)) stop("Could not find 00_config.R. Run from the script folder or set UBR5_CONFIG.")
  config_file <- normalizePath(hits[1], winslash = "/", mustWork = TRUE)
} else {
  config_file <- normalizePath(this_file_candidates[1], winslash = "/", mustWork = TRUE)
}
script_dir <- dirname(config_file)
source(config_file, local = FALSE)

scripts <- c(
  "01_import_deseq2.R",
  "02_cell_state_qc.R",
  "03_gsea_and_plots.R",
  "04_biology_panel_concordance.R"
)

run_log <- tibble::tibble(
  script = character(),
  status = character(),
  started = as.POSIXct(character()),
  finished = as.POSIXct(character()),
  message = character()
)

pipeline_failed <- FALSE
for (s in scripts) {
  script_path <- file.path(script_dir, s)
  if (!file.exists(script_path)) stop("Missing script: ", script_path)
  started <- Sys.time()
  message("\n============================================================")
  message("Running ", s)
  message("============================================================")
  status <- "success"
  msg <- ""
  warning_messages <- character()
  tryCatch(
    {
      withCallingHandlers(
        {
          source(script_path, local = FALSE)
        },
        warning = function(w) {
          warning_messages <<- c(warning_messages, conditionMessage(w))
          message("WARNING in ", s, ": ", conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    },
    error = function(e) {
      status <<- "error"
      msg <<- conditionMessage(e)
      pipeline_failed <<- TRUE
      message("ERROR in ", s, ": ", msg)
    }
  )
  if (status == "success" && length(warning_messages) > 0) {
    status <- "success_with_warnings"
    msg <- paste(unique(warning_messages), collapse = " | ")
  }
  finished <- Sys.time()
  run_log <- dplyr::bind_rows(run_log, tibble::tibble(
    script = s, status = status, started = started, finished = finished, message = msg
  ))
  utils::write.csv(run_log, p_log(paste0("pipeline_run_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")), row.names = FALSE)
  if (pipeline_failed) break
}

if (pipeline_failed) {
  stop("Pipeline stopped because at least one script failed. See results/logs for the run log.")
} else {
  message("\nPipeline completed successfully.")
}
