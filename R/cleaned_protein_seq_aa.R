#' Clean protein sequences and accessions for IDPR processing
#'
#' This function prepares real-world protein data for the `idpr` package.
#' It ensures that sequences contain only valid amino acids (A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V)
#' and converts them to uppercase. It also normalizes UniProt-style accessions from
#' "sp|P04637|P53_HUMAN" → "P04637".
#'
#' @param sequences A character vector of protein sequences.
#' @param accessions A character vector of protein accessions (names of sequences or explicit vector).
#' @param warn Logical; if TRUE, prints warnings for sequences where invalid residues were removed.
#'
#' @return A list with cleaned sequences and accessions:
#' \item{sequences}{Character vector of cleaned protein sequences.}
#' \item{accessions}{Character vector of cleaned UniProt accessions.}
#'
#' @examples
#' raw_seq <- c("sp|P04637|P53_HUMAN" = "MKQTALV123$")
#' cleaned <- clean_sequences_and_accessions(names(raw_seq), raw_seq)
#' cleaned$sequences
#' cleaned$accessions
#' @export
clean_sequences_and_accessions <- function(sequences, accessions, warn = TRUE) {

  # --- validate inputs ---
  if (!is.character(sequences))
    stop("'sequences' must be a character vector.", call. = FALSE)

  if (!is.character(accessions))
    stop("'accessions' must be a character vector.", call. = FALSE)

  if (length(sequences) != length(accessions))
    stop("Lengths of 'sequences' and 'accessions' must match", call. = FALSE)

  # --- valid amino acids ---
  validAA <- "ACDEFGHIKLMNPQRSTVWY"

  cleaned_sequences <- character(length(sequences))
  cleaned_accessions <- character(length(accessions))

  for (i in seq_along(sequences)) {

    seq_raw <- sequences[i]

    # uppercase first
    seq_up <- toupper(seq_raw)

    # remove invalid characters
    seq_clean <- gsub(paste0("[^", validAA, "]"), "", seq_up)

    # warn if changed
    if (warn && seq_clean != seq_up) {
      warning(
        sprintf(
          "Sequence %d: invalid residues removed (%s → %s)",
          i, seq_raw, seq_clean
        ),
        call. = FALSE
      )
    }

    cleaned_sequences[i] <- seq_clean

    # ---- accession cleaning ----

    acc_raw <- accessions[i]

    # Try UniProt-style split
    if (grepl("\\|", acc_raw)) {
      # sp|P04637|NAME → take middle
      parts <- strsplit(acc_raw, "\\|")[[1]]
      if (length(parts) >= 2) {
        cleaned_accessions[i] <- parts[2]
      } else {
        cleaned_accessions[i] <- acc_raw
      }
    } else {
      cleaned_accessions[i] <- acc_raw
    }
  }

  # return
  list(
    sequences  = cleaned_sequences,
    accessions = cleaned_accessions
  )
}
