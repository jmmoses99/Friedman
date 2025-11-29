# Function for parallel processing

# parallel_processing_backend.R
#' @keywords internal
#' @import BiocParallel

# parallel_backend.R
# Internal helper to choose a parallel backend

parallel_backend <- function(cores = parallel::detectCores() - 1) {
  cores <- max(1, cores)
  if (.detect_macOs_internal(check_systems = TRUE)) {
    BiocParallel::MulticoreParam(workers = cores)
  } else {
    # Use SOCK/Snow cluster for Windows / non-macOS

    BiocParallel::SnowParam(
      workers       = cores,
      type          = "SOCK",
      exportglobals = TRUE,
      exportvariables = TRUE,
      stop.on.error = FALSE
    )
  }
}
