testthat::context("Testing idpr_parallel_processing on human P53")

library(IDPsBio)
library(Biostrings)
library(testthat)

# ------------------------------------------------
# Load once (for all tests)
# ------------------------------------------------
data("P53_Human_Protein", package = "IDPsBio")

# ------------------------------------------------
# IDPR function set & backend
# ------------------------------------------------
idpr_funcs <- c(
  "idprofile",
  "iupred",
  "iupredAnchor",
  "iupredRedox",
  "chargeCalculationLocal",
  "chargeCalculationGlobal",
  "foldIndexR",
  "scaledHydropathyLocal",
  "meanScaledHydropathy"
)

backend <- IDPsBio::parallel_backend(cores = 1)

# ------------------------------------------------
# Test 0 — Dataset validity
# ------------------------------------------------
test_that("P53_Human_Protein loads correctly", {

  expect_true(exists("P53_Human_Protein"))
  expect_s4_class(P53_Human_Protein, "AAStringSet")
  expect_equal(length(P53_Human_Protein), 1)

  expect_true(nzchar(as.character(P53_Human_Protein)))
  expect_true(nzchar(names(P53_Human_Protein)))
})

# ------------------------------------------------
# Test 1 — AAStringSet input
# ------------------------------------------------
test_that("idpr_parallel_processing works on P53 with AAStringSet", {

  raw_seq <- as.character(P53_Human_Protein)
  raw_acc <- names(P53_Human_Protein)

  cleaned <- IDPsBio::clean_sequences_and_accessions(
    sequences  = raw_seq,
    accessions = raw_acc,
    warn = FALSE
  )

  results <- IDPsBio::idpr_parallel_processing(
    sequences       = cleaned$sequences,
    idpr_functions  = idpr_funcs,
    accessions      = cleaned$accessions,
    BPPARAM         = backend
  )

  expect_type(results, "list")
  expect_equal(names(results), cleaned$accessions)

  for (res in results) {
    expect_true(all(c("numeric", "plots") %in% names(res)))
    expect_type(res$numeric, "list")
    expect_type(res$plots, "list")
  }
})

# ------------------------------------------------
# Test 2 — Named character vector input
# ------------------------------------------------
test_that("idpr_parallel_processing works with named character vector", {

  sequences <- as.character(P53_Human_Protein)
  names(sequences) <- names(P53_Human_Protein)
  accessions <- names(sequences)

  cleaned <- IDPsBio::clean_sequences_and_accessions(
    sequences  = sequences,
    accessions = accessions,
    warn = FALSE
  )

  results <- IDPsBio::idpr_parallel_processing(
    sequences      = cleaned$sequences,
    idpr_functions = idpr_funcs,
    accessions     = cleaned$accessions,
    BPPARAM        = backend
  )

  expect_type(results, "list")
  expect_equal(names(results), cleaned$accessions)

  for (res in results) {
    expect_true(all(c("numeric", "plots") %in% names(res)))
  }
})

# ------------------------------------------------
# Test 3 — data.frame input
# ------------------------------------------------
test_that("idpr_parallel_processing works with data.frame input", {

  df_input <- data.frame(
    accession = names(P53_Human_Protein),
    sequence  = as.character(P53_Human_Protein),
    stringsAsFactors = FALSE
  )

  cleaned <- IDPsBio::clean_sequences_and_accessions(
    sequences  = df_input$sequence,
    accessions = df_input$accession,
    warn = FALSE
  )

  results <- IDPsBio::idpr_parallel_processing(
    sequences      = cleaned$sequences,
    idpr_functions = idpr_funcs,
    accessions     = cleaned$accessions,
    BPPARAM        = backend
  )
  expect_type(results, "list")
  expect_equal(names(results), cleaned$accessions)

  for (res in results) {
    expect_true(all(c("numeric", "plots") %in% names(res)))
  }
})
