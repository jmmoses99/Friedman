# Function for parallel processing

# parallel_processing_backend.R
#' @keywords internal
#' @import BiocParallel

parallel_backend <- function() {

  # Detect number of cores, leave 1 free
  cores <- max(parallel::detectCores() - 1, 1)

  # Choose backend depending on OS
  if (.detect_macOs_internal(check_systems = TRUE)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = cores)
  } else {
    BPPARAM <- BiocParallel::SnowParam(workers = cores)

    # Automatically load required packages and functions on each worker
    parallel::clusterEvalQ(BPPARAM$workers, {
      # Load required packages
      require(IDPsBio)
      require(Biostrings)

      # Export utility functions to workers

      assign("cleaned_protein_seq_aa", IDPsBio::cleaned_protein_seq_aa, envir = .GlobalEnv)
      assign("idprofile", IDPsBio::idprofile, envir = .GlobalEnv)
      assign("iupred", IDPsBio::iupred, envir = .GlobalEnv)
      assign("iupredRedox", IDPsBio::iupredRedox, envir = .GlobalEnv)
      assign("chargeCalculationLocal", IDPsBio::chargeCalculationLocal, envir = .GlobalEnv)
      assign("chargeCalculationGlobal", IDPsBio::chargeCalculationGlobal, envir = .GlobalEnv)
      assign("foldIndexR", IDPsBio::foldIndexR, envir = .GlobalEnv)
      assign("scaledHydropathyLocal", IDPsBio::scaledHydropathyLocal, envir = .GlobalEnv)
      assign("meanScaledHydropathy", IDPsBio::meanScaledHydropathy, envir = .GlobalEnv)
    })
  }

  return(BPPARAM)
}
