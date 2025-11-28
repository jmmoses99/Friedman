# Function for parallel processing

#' @keywords internal

parallel_backend <- function(){

  # Detect number of cores, leaving 1 free
  # Very critical because it prevents all cores from being used on parallel processing;
  #ensures at least 1 core is used even on single-core systems

  cores <- max(parallel::detectCores() - 1, 1)

  # Calls internal Utility Function to detect macOs and then decides between
  # MulticoreParam for Mac/Linux or SnowParam for Windows
  # Had to add (check_systems = TRUE) otherwise .detect_macOs_internal passes NULL and error

  if (.detect_macOs_internal(check_systems = TRUE)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = cores)
  } else {
    BPPARAM <- BiocParallel::SnowParam(workers = cores)
  }

  return(BPPARAM) # Explicitly returns the BiocParallelParam object
}



