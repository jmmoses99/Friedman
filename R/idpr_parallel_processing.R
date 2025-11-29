#' Parallel execution of IDPsBio functions on a protein set
#'
#' @description
#' This function runs selected IDPsBio functions in parallel over a named vector of
#' protein sequences. It supports Windows, macOS, and Linux backends through BiocParallel.
#' Avoids relying on `library()` in workers; functions are explicitly passed.
#'Run multiple IDPR computations in parallel
#'
#' This function runs a set of IDPR-related functions on a collection of
#' protein sequences in parallel using BiocParallel.
#'
#' @param idpr_functions Character vector of functions to run (must be valid functions)
#' @param valid_functions Character vector of valid function names (used for validation)
#' @param input Named character vector of protein sequences
#' @param BPPARAM A BiocParallelParam object (default: `parallel_backend()`)
#'
#' @return A list of results, one element per input sequence. Each element
#'   is a named list with entries for each requested function.
#'
#' @examples
#'
#' \dontrun{
#' sequences <- c("P12345" = "MKTFFVAGA", "Q9XYZ1" = "ACDEFGHIKLMNPQRSTVWY")
#' idpr_functions <- c("idprofile", "chargeCalculationLocal")
#' valid_functions <- idpr_functions
#' if (requireNamespace("BiocParallel", quietly = TRUE)) {
#'   BPPARAM <- BiocParallel::SnowParam(1)
#'   results <- idpr_parallel_processing(
#'       idpr_functions = idpr_functions,
#'       valid_functions = valid_functions,
#'       input = sequences,
#'       BPPARAM = BPPARAM
#'   )
#'   str(results)
#' }
#' }
#' @export



idpr_parallel_processing <- function(
    idpr_functions,
    valid_functions,
    input,
    BPPARAM = parallel_backend()
) {
  # Convert input list to data.frame
  input_df <- data.frame(
    accession = names(input),
    sequence  = unlist(input),
    stringsAsFactors = FALSE
  )

  # Clean sequences
  input_df$sequence <- cleaned_protein_seq_aa(input_df$sequence, warn = TRUE)

  # Validate requested functions
  validate_idpr_functions(
    idpr_functions = idpr_functions,
    valid_functions = valid_functions,
    accessions = input_df$accession,
    sequences = input_df$sequence,
    input = input_df
  )

  # Map function names to actual idpr functions
  func_map <- list(
    idprofile              = idpr::idprofile,
    iupred                 = idpr::iupred,
    iupredAnchor           = idpr::iupredAnchor,
    iupredRedox            = idpr::iupredRedox,
    chargeCalculationLocal = idpr::chargeCalculationLocal,
    chargeCalculationGlobal= idpr::chargeCalculationGlobal,
    foldIndexR             = idpr::foldIndexR,
    scaledHydropathyLocal  = idpr::scaledHydropathyLocal,
    meanScaledHydropathy   = idpr::meanScaledHydropathy
  )

  # Worker function for each sequence
  worker_fun <- function(i, df, funcs, fmap) {
    # Ensure needed namespaces inside worker
    requireNamespace("stats", quietly = TRUE)
    requireNamespace("Biostrings", quietly = TRUE)
    requireNamespace("IDPsBio", quietly = TRUE)

    row <- df[i, , drop = FALSE]
    out <- list(accession = row$accession)

    # Apply each requested function to the sequence
    for (fn in funcs) {
      if (!fn %in% names(fmap)) {
        out[[fn]] <- NA
      } else {
        out[[fn]] <- tryCatch(
          fmap[[fn]](row$sequence),
          error = function(e) {
            warning("Function ", fn, " failed for accession ", row$accession, ": ", e$message)
            NA
          }
        )
      }
    }
    out
  }

  # Run in parallel
  results_raw <- tryCatch({
    BiocParallel::bplapply(
      seq_len(nrow(input_df)),
      worker_fun,
      df    = input_df,
      funcs = idpr_functions,
      fmap  = func_map,
      BPPARAM = BPPARAM
    )
  }, error = function(e) {
    # Warn and fallback to serial if parallel fails
    warning("Parallel run failed: ", conditionMessage(e),
            "; falling back to serial processing.")
    NULL
  })

  # Fallback to serial processing if parallel failed
  if (is.null(results_raw)) {
    results_raw <- lapply(
      seq_len(nrow(input_df)),
      function(i) worker_fun(i, input_df, idpr_functions, func_map)
    )
  }

  # Ensure every element is a list with all requested functions
  results <- lapply(results_raw, function(el) {
    if (!is.list(el)) {
      el <- list(accession = NA)
    }
    for (fn in idpr_functions) {
      if (!fn %in% names(el)) el[[fn]] <- NA
    }
    el
  })

  # Warn if any workers failed
  failed <- sapply(results_raw, function(el) !is.list(el))
  if (any(failed)) {
    warning(sum(failed), " of ", length(results),
            " parallel tasks failed or returned invalid results.")
  }

  # Convert results to data.frame
  results_df <- do.call(rbind, lapply(results, function(x) {
    # Ensure order: accession first, then functions
    c(accession = x$accession, unlist(x[idpr_functions]))
  }))

  rownames(results_df) <- NULL
  as.data.frame(results_df, stringsAsFactors = FALSE)
}
