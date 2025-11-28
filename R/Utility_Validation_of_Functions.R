#Utility Function to ensure the validation of correct idpr functions using either a data.frame or FASTA

#' @keywords internal

#Argument validation Function for idpr Function names

validate_idpr_functions <- function(idpr_functions, valid_functions, input, accessions, sequences) {

  # Validate function names
  invalid_functions <- setdiff(idpr_functions, valid_functions)
  if (length(invalid_functions) > 0) {
    stop(
      paste0(
        "Invalid IDPR functions specified: ",
        paste(invalid_functions, collapse = ", "),
        ". Valid functions are: ",
        paste(valid_functions, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  # Functions requiring UniProt accessions
  required_accession <- c("idprofile", "iupred", "iupredAnchor", "iupredRedox")
  needs_accession <- any(idpr_functions %in% required_accession)

  # Validate accessions if needed
  if (needs_accession) {
    if (is.null(accessions)) {
      stop(
        "The selected idpr functions require UniProt accessions, but 'accessions' was not supplied.",
        call. = FALSE
      )
    }
    if (!is.character(accessions)) {
      stop("'accessions' must be a character vector.", call. = FALSE)
    }
  }

  # Validate sequences
  if (!is.null(sequences)) {
    if (!is.character(sequences)) {
      stop("'sequences' must be a character vector.", call. = FALSE)
    }

    if (!is.null(accessions) && length(sequences) != length(accessions)) {
      stop(
        paste0(
          "Lengths of 'accessions' (", length(accessions),
          ") and 'sequences' (", length(sequences), ") do not match."
        ),
        call. = FALSE
      )
    }

    if (any(is.na(sequences))) {
      stop("Sequences cannot contain missing values (NA).", call. = FALSE)
    }
  }

  # Validate 'input'
  if (!is.null(input)) {

    # FASTA file path
    if (is.character(input) && length(input) == 1 && file.exists(input)) {
      # valid case → nothing to do
    }

    # data.frame with accession + sequence
    else if (is.data.frame(input)) {

      required_columns <- c("accession", "sequence")

      if (!all(required_columns %in% names(input))) {
        stop(
          paste0(
            "Input data.frame must contain columns: ",
            paste(required_columns, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      if (any(is.na(input$accession)) || any(is.na(input$sequence))) {
        stop(
          "Columns 'accession' and 'sequence' cannot contain NA values.",
          call. = FALSE
        )
      }
    }

    else {
      stop(
        "'input' must be a FASTA file path or a data.frame with 'accession' and 'sequence' columns.",
        call. = FALSE
      )
    }
  }

}






