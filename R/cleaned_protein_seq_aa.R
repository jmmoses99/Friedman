#' Clean protein sequences to keep only valid amino acids
#'
#' Created because using `idpr` on real-world data sometimes produces errors
#' when sequences contain invalid residues. This function takes a character
#' vector of protein sequences and removes all characters that are not valid
#' one-letter amino acid codes. It converts sequences to uppercase.
#'
#' Example of potential error avoided:
#' Error in sequenceCheck(sequence = sequence, method = "stop", outputType = "vector", :
#' Protein contains the following invalid residues
#'
#' Valid amino acids are:
#' \code{A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V}
#'
#' @param sequences A character vector of protein sequences.
#' @param warn Logical, if \code{TRUE}, prints a warning for sequences with removed residues.
#'
#' @return A character vector of cleaned protein sequences with only valid uppercase amino acids.
#'
#' @examples
#' sequences <- c("P12345" = "MKQTALV123$", "Q9XYZ1" = "ACD-ZEFG!")
#' cleaned <- cleaned_protein_seq_aa(sequences)
#' cleaned
#' @export
cleaned_protein_seq_aa <- function(sequences, warn = TRUE) {

  # sequences: a character vector of protein sequences
  # warn: whether to print warnings if residues are removed

  # Ensure input is a character vector to prevent errors if user passes numeric or factor

  if (!is.character(sequences)) {
    stop("'sequences' must be a character vector", call. = FALSE)
  }

  # Keep only valid Amino acids; toupper() is a base R function that converts all characters in a string to uppercase.

  # Clarified regex comment
  cleaned <- toupper(gsub("[^ARNDCEQGHILKMFPSTWYV]", "", sequences))  # gsub() replaces all invalid residues (anything not in [ARNDCEQGHILKMFPSTWYV]) with ""

  if (warn) {
    for (i in seq_along(sequences)) {
      if (nchar(sequences[i]) != nchar(cleaned[i])) {  # Compares the original sequence length to the cleaned sequence length. If they differ → some characters were removed and invalid residues existed
        warning(paste("Sequence", i, "had invalid residues removed"))
      }
    }
  }

  return(cleaned) # Explicitly return cleaned sequences
}
