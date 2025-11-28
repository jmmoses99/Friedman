# Function for parallel processing; Need to install Bioconductor
# Utility Functions for parallel back end selection
 #' @keywords internal  (helper)

.detect_macOs_internal <- function(systems = NULL, cores = 2, check_systems = FALSE) {

  # Simple, safe OS check
  if (check_systems) {
    return(Sys.info()[["sysname"]] == "Darwin")
  }

  # If using the vector-check mode:
  if (is.null(systems)) {
    stop("Please provide vector of system strings")
  }

  # This DOES NOT use parallel processing (safe on macOS)

  results <- vapply(
    systems,
    function(x) grepl("Macintosh|Darwin", x, ignore.case = TRUE),
    logical(1)
  )

  return(results)
}
