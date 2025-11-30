#' Choose a parallel backend for IDPR processing
#'
#' Automatically selects the appropriate BiocParallel backend depending on the operating system:
#' * On macOS, uses `MulticoreParam`.
#' * On Windows or other OS, uses `SnowParam` (SOCK cluster).
#'
#' This ensures consistent parallel processing behavior across platforms for functions like
#' `idpr_parallel_processing()`.
#'
#' @param cores Number of cores to use (default = detected cores - 1). Minimum is 1.
#'
#' @return A `BiocParallelParam` object suitable for use in parallel Bioconductor workflows.
#'
#' @examples
#' \dontrun{
#' BPPARAM <- IDPsBio::parallel_backend(cores = 4)
#' results <- idpr_parallel_processing(
#'   idpr_functions = c("iupred", "scaledHydropathyLocal"),
#'   input = P53_Human_Protein,
#'   BPPARAM = BPPARAM
#' )
#' }
#'
#' @export
#' @import BiocParallel
parallel_backend <- function(cores = parallel::detectCores() - 1) {

  # Ensure at least 1 core
  cores <- max(1, cores)

  # Platform check helper
  is_mac <- .detect_macOs_internal(check_systems = TRUE)

  if (is_mac) {
    # macOS: use multicore
    BiocParallel::MulticoreParam(workers = cores)
  } else {
    # Windows / other: use Snow/SOCK cluster
    BiocParallel::SnowParam(
      workers       = cores,
      type          = "SOCK",
      exportglobals = TRUE,
      exportvariables = TRUE,
      stop.on.error = FALSE
    )
  }
}
