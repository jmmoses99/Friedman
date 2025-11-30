#' Parallel processing for IDPR functions
#'Safely parallelizes numeric-returning IDPR functions while running plotting functions sequentially.
#'
#' @param sequences Character vector of protein sequences
#' @param idpr_functions Character vector of IDPR function names to run
#' @param accessions Optional character vector of UniProt accessions
#' @param BPPARAM BiocParallel backend (default uses internal `parallel_backend()`)
#' @return List of results for each sequence
#' @import idpr
#'
#' @export
idpr_parallel_processing <- function(
    sequences,
    idpr_functions,
    accessions = NULL,
    BPPARAM = BiocParallel::bpparam()
) {

  # --- Load idpr namespace for safety ---
  if (!requireNamespace("idpr", quietly = TRUE)) {
    stop("The 'idpr' package must be installed.", call. = FALSE)
  }

  # --- Input handling -------------------------------------------------------

  # Convert AAStringSet → character
  if (inherits(sequences, "AAStringSet")) {
    accessions <- names(sequences)
    sequences  <- as.character(sequences)
  }

  # Data.frame input
  if (is.data.frame(sequences)) {
    stop("Pass df$sequence, not the whole data.frame, to `sequences`.", call. = FALSE)
  }

  # Named character vector input
  if (is.null(accessions)) {
    if (!is.null(names(sequences))) {
      accessions <- names(sequences)
    } else {
      stop("Accessions must be supplied or sequences must be named.")
    }
  }

  # --- Clean sequences + accessions ----------------------------------------

  cleaned <- clean_sequences_and_accessions(
    sequences = sequences,
    accessions = accessions,
    warn = TRUE
  )

  sequences  <- cleaned$sequences
  accessions <- cleaned$accessions

  # --- Classify idpr functions ---------------------------------------------

  numeric_functions <- c(
    "idprofile",
    "iupred", "iupredAnchor", "iupredRedox",
    "chargeCalculationLocal", "chargeCalculationGlobal",
    "netCharge",
    "scaledHydropathyLocal", "scaledHydropathyGlobal",
    "meanScaledHydropathy",
    "foldIndexR"
  )

  plot_functions <- c(
    "chargeHydropathyPlot", "sequencePlot",
    "sequenceMap", "sequenceMapCoordinates",
    "structuralTendencyPlot"
  )

  # Validate names
  valid <- c(numeric_functions, plot_functions)
  invalid <- setdiff(idpr_functions, valid)
  if (length(invalid) > 0) {
    stop(sprintf(
      "Invalid idpr functions: %s",
      paste(invalid, collapse = ", ")
    ))
  }

  # --- Safe wrapper to call idpr functions ----------------------------------

  run_idpr_func_safe <- function(func_name, sequence, uniprotAccession) {

    func <- tryCatch(
      get(func_name, asNamespace("idpr")),
      error = function(e) NULL
    )

    if (is.null(func)) {
      warning("Function not found in idpr: ", func_name)
      return(NULL)
    }

    # Try all argument combinations
    attempts <- list(
      list(sequence = sequence, uniprotAccession = uniprotAccession),
      list(sequence = sequence),
      list(uniprotAccession = uniprotAccession)
    )

    for (a in attempts) {
      out <- tryCatch(
        do.call(func, a),
        error = function(e) e
      )
      if (!inherits(out, "error"))
        return(out)
    }

    warning(sprintf(
      "Skipping '%s': no valid argument combination worked.",
      func_name
    ))
    return(NULL)
  }

  # --- Parallel execution ---------------------------------------------------

  results <- BiocParallel::bplapply(seq_along(sequences), function(i) {

    seq_i <- sequences[i]
    acc_i <- accessions[i]

    # Numeric
    numeric_out <- lapply(
      intersect(idpr_functions, numeric_functions),
      run_idpr_func_safe,
      sequence = seq_i,
      uniprotAccession = acc_i
    )
    numeric_out <- numeric_out[!sapply(numeric_out, is.null)]
    names(numeric_out) <- names(numeric_out)[!sapply(numeric_out, is.null)]

    # Plots
    plot_out <- lapply(
      intersect(idpr_functions, plot_functions),
      run_idpr_func_safe,
      sequence = seq_i,
      uniprotAccession = acc_i
    )
    plot_out <- plot_out[!sapply(plot_out, is.null)]
    names(plot_out) <- names(plot_out)[!sapply(plot_out, is.null)]

    list(
      numeric = numeric_out,
      plots = plot_out
    )
  }, BPPARAM = BPPARAM)

  names(results) <- accessions
  return(results)
}
