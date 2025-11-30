library(testthat)

# --- Function to test ---
clean_sequences_and_accessions <- function(sequences, accessions, warn = TRUE) {
  # Validate inputs
  if (!is.character(sequences)) stop("'sequences' must be a character vector.", call. = FALSE)
  if (!is.character(accessions)) stop("'accessions' must be a character vector.", call. = FALSE)
  if (length(sequences) != length(accessions)) stop("Lengths of 'sequences' and 'accessions' must match", call. = FALSE)

  validAA <- "ACDEFGHIKLMNPQRSTVWY"
  cleaned_sequences <- character(length(sequences))
  cleaned_accessions <- character(length(accessions))

  for (i in seq_along(sequences)) {
    # Clean sequence
    seq_up <- toupper(sequences[i])
    seq_clean <- gsub(paste0("[^", validAA, "]"), "", seq_up)
    if (warn && seq_clean != seq_up) {
      warning(sprintf("Sequence %d: invalid residues removed (%s → %s)", i, sequences[i], seq_clean), call. = FALSE)
    }
    cleaned_sequences[i] <- seq_clean

    # Clean accession
    acc_raw <- accessions[i]
    if (grepl("\\|", acc_raw)) {
      parts <- strsplit(acc_raw, "\\|")[[1]]
      # Only take middle field if it exists and is non-empty
      if(length(parts) >= 2 && nzchar(parts[2])) {
        cleaned_accessions[i] <- parts[2]
      } else {
        cleaned_accessions[i] <- acc_raw
      }
    } else {
      cleaned_accessions[i] <- acc_raw
    }
  }

  list(sequences = cleaned_sequences, accessions = cleaned_accessions)
}

# --- Tests ---

test_that("Basic cleaning works", {
  sequences  <- c("MKQTALV123$", "acgtxyzMKQ")
  accessions <- c("sp|P04637|P53_HUMAN", "tr|Q9XYZ1|SOME_PROT")

  expect_warning(
    result <- clean_sequences_and_accessions(sequences, accessions, warn = TRUE),
    "invalid residues removed"
  )

  # Y is valid, so it should remain
  expect_equal(result$sequences, c("MKQTALV", "ACGTYMKQ"))
  expect_equal(result$accessions, c("P04637", "Q9XYZ1"))
})

test_that("No warning if warn = FALSE", {
  sequences  <- c("MKQ123")
  accessions <- c("sp|A0A000|TEST")

  expect_silent(
    result <- clean_sequences_and_accessions(sequences, accessions, warn = FALSE)
  )

  expect_equal(result$sequences, "MKQ")
  expect_equal(result$accessions, "A0A000")
})

test_that("Uppercasing works", {
  sequences  <- c("mkqa")
  accessions <- c("ABC123")

  result <- clean_sequences_and_accessions(sequences, accessions, warn = FALSE)
  expect_equal(result$sequences, "MKQA")
})

test_that("Accessions without pipes remain unchanged", {
  sequences  <- c("MKQA")
  accessions <- c("P12345")

  result <- clean_sequences_and_accessions(sequences, accessions, warn = FALSE)
  expect_equal(result$accessions, "P12345")
})

test_that("Input validation works", {
  expect_error(clean_sequences_and_accessions(123, c("A")), "'sequences' must be a character vector.")
  expect_error(clean_sequences_and_accessions("ABC", 5), "'accessions' must be a character vector.")
  expect_error(clean_sequences_and_accessions(c("ABC", "XYZ"), "P04637"), "Lengths of 'sequences' and 'accessions' must match")
})

test_that("Edge cases: empty sequences", {
  sequences <- c("", "ACD")
  accessions <- c("P1", "P2")

  result <- clean_sequences_and_accessions(sequences, accessions, warn = TRUE)
  expect_equal(result$sequences, c("", "ACD"))
  expect_equal(result$accessions, c("P1", "P2"))
})

test_that("Edge cases: sequences with only invalid characters", {
  sequences <- c("123$#@", "!!@@")
  accessions <- c("P1", "P2")

  expect_warning(
    result <- clean_sequences_and_accessions(sequences, accessions, warn = TRUE),
    "invalid residues removed"
  )

  expect_equal(result$sequences, c("", ""))
})

test_that("Edge cases: accessions with missing middle field", {
  sequences <- c("ACDE")
  accessions <- c("sp||PROT")  # missing middle

  result <- clean_sequences_and_accessions(sequences, accessions, warn = FALSE)
  expect_equal(result$accessions, "sp||PROT")
})

test_that("Edge cases: sequences with all valid amino acids", {
  sequences <- c("ACDEFGHIKLMNPQRSTVWY")
  accessions <- c("P12345")

  result <- clean_sequences_and_accessions(sequences, accessions, warn = TRUE)
  expect_equal(result$sequences, "ACDEFGHIKLMNPQRSTVWY")
})

test_that("Mixed valid/invalid sequences produce correct cleaning", {
  sequences <- c("ACDXYZEFG")  # X,Z invalid, Y valid
  accessions <- c("P67890")

  expect_warning(
    result <- clean_sequences_and_accessions(sequences, accessions, warn = TRUE),
    "invalid residues removed"
  )

  expect_equal(result$sequences, "ACDYEFG")
})
