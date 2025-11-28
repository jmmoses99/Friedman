#' Parallel execution of IDPsBio functions on a protein set
#'
#' @description
#' This function runs selected IDPsBio functions in parallel over a named vector of
#' protein sequences. It supports Windows, macOS, and Linux backends through BiocParallel.
#' Avoids relying on `library()` in workers; functions are explicitly passed.
#'
#' @param idpr_functions Character vector of all requested functions (including potentially invalid)
#' @param valid_functions Character vector of valid functions to actually execute
#' @param input Named character vector of protein sequences (names = accessions)
#' @param BPPARAM BiocParallelParam object, e.g., from `parallel_backend()`
#'
#' @return List of results for each protein, each element containing `accession` and function outputs
#'
#' @import BiocParallel
#' @importFrom stats median sd
#' @importFrom utils write.csv
#' @export
idpr_parallel_processing <- function(idpr_functions, valid_functions, input, BPPARAM) {


  # 1. Prepare input Data.frame

  input_df <- data.frame(
    accession = names(input),
    sequence  = unlist(input),
    stringsAsFactors = FALSE
  )

  # Clean sequences before parallel execution
  # Removes invalid amino acids
  input_df$sequence <- cleaned_protein_seq_aa(input_df$sequence, warn = TRUE)

  # Validate requested functions
  validate_idpr_functions(
    idpr_functions = idpr_functions,
    valid_functions = valid_functions,
    accessions = input_df$accession,
    sequences = input_df$sequence,
    input = input_df
  )


  # 2. Define mapping of function names to actual R functions

  # Fix: Explicitly map names to functions to avoid relying on library() in workers
  func_map <- list(
    idprofile              = IDPsBio::idprofile,
    iupred                 = IDPsBio::iupred,
    iupredAnchor           = IDPsBio::iupredAnchor,
    iupredRedox            = IDPsBio::iupredRedox,
    chargeCalculationLocal = IDPsBio::chargeCalculationLocal,
    chargeCalculationGlobal= IDPsBio::chargeCalculationGlobal,
    foldIndexR             = IDPsBio::foldIndexR,
    scaledHydropathyLocal  = IDPsBio::scaledHydropathyLocal,
    meanScaledHydropathy   = IDPsBio::meanScaledHydropathy
  )


  # 3. Parallel processing

  results <- BiocParallel::bplapply(
    seq_len(nrow(input_df)),

    function(i, df, funcs, fmap) {


      # Keep row as a data.frame to avoid vectorization issues

      row <- df[i, , drop = FALSE]
      out <- list(accession = row$accession)

      # Loop over requested functions
      for (func in funcs) {
        # Only run valid functions
        if (!func %in% names(fmap)) {
          out[[func]] <- NA
          next
        }

        # Execute function safely
        out[[func]] <- tryCatch(
          fmap[[func]](row$sequence),
          error = function(e) {
            warning(
              paste("Function", func, "failed for accession",
                    row$accession, ":", e$message)
            )
            return(NA)
          }
        )
      }

      out
    },

    df = input_df,
    funcs = idpr_functions,
    fmap = func_map,
    BPPARAM = BPPARAM
  )


  # 4. Return results

  return(results)
}
