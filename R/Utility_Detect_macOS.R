# Function for parallel processing; Need to install Bioconductor
# Utility Functions for parallel back end selection
 #' @keywords internal  (helper)
 #' @import BiocParallel

#BiocParallel
BiocParallel::bplapply


# Detection of System/Correct BiocParallel Back end
# Set up Cluster/Auto-Detect OS function
# Parallel processing requires you to understand the processing power of the computer you are using

# Internal Utility for parallel macOS detection
.detect_macOs_internal <- function(systems = NULL, cores = 2, check_systems = FALSE){

  if(check_systems){
    return(Sys.info()["sysname"]=="Darwin")
  }

  if(is.null(systems)) stop("Please provide vector of system strings")

  # Only use available CPU cores
  max_cores <- parallel::detectCores()
  cores <- min(cores, max_cores)
  if(cores < 1) cores <- 1  # Ensure at least 1 core

  # Choose backend: Multicore for macOs/Linux and Snow for Windows
  cluster <- if(.Platform$OS.type == "windows") {
    SnowParam(workers = cores)
  } else {
    MulticoreParam(workers = cores)
  }

  # Detection function for macOs/Parallel execution
  results <- bplapply(
    systems,
    function(x){
      # explicitly making sure stats package is loaded in each worker

      base::grepl("Macintosh|Darwin", x, ignore.case = TRUE)
    },
    BPPARAM = cluster
  )

  return(unlist(results))
}
