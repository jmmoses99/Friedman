# Function for parallel processing; Need to install Bioconductor
# Utility Functions for parallel back end selection
#' Internal function to detect macOS systems
#' This is an internal helper function used by IDPsBio to detect macOS.
#' Detect macOS Systems (Internal)
#' It returns `TRUE` if the system is macOS, otherwise `FALSE`.
#'
#' @param cores Integer; number of cores (currently unused, reserved for future parallel use)
#' @param check_systems Logical; if `TRUE`, perform a simple system check
#'
#' @return Logical indicating if the system is macOS
#' @keywords internal
#' @name detect_macOs_internal
NULL

.detect_macOs_internal <- function(check_systems = FALSE, cores = 2) {
  # cores argument is kept for backward compatibility, but not used
  # check_systems must be a single logical value
  if (!is.logical(check_systems) || length(check_systems) != 1) {
    stop("'check_systems' must be a single logical value (TRUE/FALSE)")
  }

  # Detect current OS
  current_system <- Sys.info()[["sysname"]]

  # If simple check requested, return TRUE if Darwin/Macintosh
  if (check_systems) {
    return(current_system %in% c("Darwin", "Macintosh"))
  }

  # Otherwise, return a vectorized check (still just one element here)
  results <- grepl("Macintosh|Darwin", current_system, ignore.case = TRUE)
  return(results)
}
