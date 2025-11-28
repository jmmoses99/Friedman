#' Parallel Processing Wrapper for intrinsically disordered proteins (IDPs) and intrinsically disordered regions (IDRs)
#'
#' This function runs multiple IDPR analyses in parallel over
#' a set of protein sequences. It performs sequence cleaning, input validation,
#' and dispatches tasks using a BiocParallel backend.
#'
#' @param idpr_functions Character vector of IDPR function names to run.
#' @param valid_functions Character vector listing all valid allowed functions.
#' @param input Named character vector of protein sequences, where
#'   names correspond to accession IDs from UniProt.
#' @param BPPARAM A BiocParallel parameter object (e.g. from
#'   \code{parallel_processing_backend()}).
#'
#' @return A list where each element contains:
#'   \describe{
#'     \item{accession}{Protein accession ID}
#'     \item{<function>}{Output of each IDPR function requested}
#'   }
#'
#' @details This function automatically:
#'   \itemize{
#'     \item Cleans protein sequences using \code{cleaned_protein_seq_aa()}
#'     \item Validates input using \code{validate_idpr_functions()}
#'     \item Runs selected IDPR computations in parallel
#'     \item Handles errors safely using \code{tryCatch()}
#'   }
#' The input must be a *named* vector, where each entry corresponds to one protein sequence.
#'
#' @examples
#' # Fake UniProt-style accessions and example sequences
#' sequences <- c(
#'   "P12345" = "MKTFFVAGA",
#'   "Q9XYZ1" = "ACDEFGHIKLMNPQRSTVWY"
#' )
#'
#' idpr_functions <- c("idprofile", "chargeCalculationLocal")
#' valid_functions <- idpr_functions
#'
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
#'
#' @import BiocParallel
#' @import idpr
#' @importFrom stats median sd
#' @importFrom utils write.csv
#' @export




idpr_parallel_processing <- function(idpr_functions, valid_functions, input, BPPARAM) {

  # Prepare input/Data.frame creation

  input_df <- data.frame(
    accession = names(input),
    sequence  = unlist(input),
    stringsAsFactors = FALSE
  )

  # Need to clean raw data/clean sequences to remove invalid residues

  # Cleaned_protein_seq_aa() already runs BEFORE parallel workers

  input_df$sequence <- cleaned_protein_seq_aa(input_df$sequence, warn = TRUE)

  # Running Utility_Validation

  validate_idpr_functions(
    idpr_functions = idpr_functions,
    valid_functions = valid_functions,
    accessions = input_df$accession,
    sequences = input_df$sequence,
    input = input_df
  )

  # Parallel processing

  results <- BiocParallel::bplapply(
    seq_len(nrow(input_df)),

    function(i, df, funcs) {

      # Keep row as a data.frame to avoid dropping to vector

      row <- df[i, , drop = FALSE]
      out <- list(accession = row$accession)

      # Loop over requested functions

      for (func in funcs) {

        out[[func]] <- tryCatch(
          switch(func,
                 idprofile              = idprofile(row$sequence),
                 iupred                 = iupred(row$sequence),
                 iupredAnchor           = iupredAnchor(row$sequence),
                 iupredRedox            = iupredRedox(row$sequence),
                 chargeCalculationLocal = chargeCalculationLocal(row$sequence),
                 chargeCalculationGlobal= chargeCalculationGlobal(row$sequence),
                 foldIndexR             = foldIndexR(row$sequence),
                 scaledHydropathyLocal  = scaledHydropathyLocal(row$sequence),
                 meanScaledHydropathy   = meanScaledHydropathy(row$sequence),
                 NA
          ),
          error = function(e) {
            warning(paste("Function", func, "failed for accession",
                          row$accession, ":", e$message))
            return(NA)
          }
        )
      }

      out
    },

    # Pass data and utility paths into worker

    df = input_df,
    funcs = idpr_functions,
    BPPARAM = BPPARAM
  )

  # Return results

  return(results)
}
