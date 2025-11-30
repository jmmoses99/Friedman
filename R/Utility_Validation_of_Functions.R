#'@keywords internal
validate_idpr_functions <- function(idpr_functions,
                                    input = NULL,
                                    accessions = NULL,
                                    sequences = NULL) {

  # -------------------------------------------------
  # 1. Valid idpr function sets
  # -------------------------------------------------
  numeric_functions <- c(
    "idprofile", "iupred", "iupredAnchor", "iupredRedox",
    "chargeCalculationLocal", "chargeCalculationGlobal",
    "netCharge",
    "scaledHydropathyGlobal", "scaledHydropathyLocal",
    "meanScaledHydropathy",
    "foldIndexR"
  )

  plot_functions <- c(
    "chargeHydropathyPlot", "sequencePlot", "structuralTendencyPlot",
    "sequenceMap", "sequenceMapCoordinates"
  )

  valid_functions <- c(numeric_functions, plot_functions)

  # -------------------------------------------------
  # 2. Validate requested function names
  # -------------------------------------------------
  invalid <- setdiff(idpr_functions, valid_functions)
  if (length(invalid) > 0) {
    stop(
      sprintf(
        "Invalid IDPR functions specified: %s. Valid functions are: %s.",
        paste(invalid, collapse = ", "),
        paste(valid_functions, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # -------------------------------------------------
  # 3. Functions requiring UniProt accessions
  # -------------------------------------------------
  must_use_accession <- c("idprofile", "iupred", "iupredAnchor", "iupredRedox")

  if (any(idpr_functions %in% must_use_accession)) {
    if (is.null(accessions)) {
      stop(
        paste(
          "One or more selected idpr functions require UniProt accessions.",
          "Provide 'accessions' when using:",
          paste(must_use_accession, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  # -------------------------------------------------
  # 4. Validate sequences (after cleaning happens outside)
  # -------------------------------------------------
  if (!is.null(sequences)) {

    if (!is.character(sequences)) {
      stop("'sequences' must be a character vector.", call. = FALSE)
    }

    if (!is.null(accessions) && length(sequences) != length(accessions)) {
      stop("'sequences' and 'accessions' must have equal length.", call. = FALSE)
    }

    if (any(is.na(sequences))) {
      stop("'sequences' cannot contain NA values.", call. = FALSE)
    }
  }

  # -------------------------------------------------
  # 5. Validate input argument
  # -------------------------------------------------
  if (!is.null(input)) {

    # FASTA file
    if (is.character(input) && length(input) == 1 && file.exists(input)) {
      # valid
    }

    # Data.frame
    else if (is.data.frame(input)) {
      req_cols <- c("accession", "sequence")

      if (!all(req_cols %in% names(input))) {
        stop("Input data.frame must contain columns: accession, sequence.", call. = FALSE)
      }

      if (any(is.na(input$accession)) || any(is.na(input$sequence))) {
        stop("Columns 'accession' and 'sequence' cannot contain NA.", call. = FALSE)
      }
    }

    # Named character vector
    else if (is.character(input)) {
      if (is.null(names(input))) {
        stop("Character vector input must be named: names = UniProt accessions.", call. = FALSE)
      }
    }

    # AAStringSet
    else if (inherits(input, "AAStringSet")) {
      # valid
    }

    else {
      stop(
        "'input' must be: FASTA file, AAStringSet, named character vector, or data.frame.",
        call. = FALSE
      )
    }
  }

  # -------------------------------------------------
  # 6. Return validated lists
  # -------------------------------------------------
  list(
    numeric_functions = intersect(idpr_functions, numeric_functions),
    plot_functions    = intersect(idpr_functions, plot_functions)
  )
}
