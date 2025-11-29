test_that("idpr_parallel_processing and summarize_idpr_results work on synthetic proteome", {

  # --- Step 1: Create a small synthetic proteome with edge cases ---
  sequences <- c(
    "P00001" = "MKTFFVAGAXZ",   # X and Z included for cleaning test
    "P00002" = "ACDEFGHIKLMNPQRSTVWY",
    "P00003" = "GGGGGGGGGG",    # repeated residues
    "P00004" = "MZZTCC"         # invalid residues to test cleaning
  )

  # --- Step 2: Define functions to run ---
  idpr_functions <- c("idprofile", "chargeCalculationLocal")
  valid_functions <- idpr_functions

  # --- Step 3: Use minimal parallel backend for testing ---
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BPPARAM <- BiocParallel::SnowParam(1)

    # --- Step 4: Run the parallel processing function ---
    results <- idpr_parallel_processing(
      idpr_functions = idpr_functions,
      valid_functions = valid_functions,
      input = sequences,
      BPPARAM = BPPARAM
    )

    # --- Step 5: Run the summarization function ---
    summary_results <- summarize_idpr_results(
      results,
      protein_outfile = NULL,   # do not write files during testing
      dataset_outfile = NULL
    )

    # --- Step 6: Check outputs ---
    expect_type(results, "list")
    expect_true(all(sapply(results, function(x) all(c("accession", idpr_functions) %in% names(x)))))
    expect_type(summary_results, "list")
    expect_true(all(c("protein_summary", "dataset_summary") %in% names(summary_results)))

    # --- Step 7: Check cleaning worked ---
    cleaned_seqs <- sapply(results, function(x) x$idprofile$sequence %||% NA)
    expect_false(any(grepl("[XZ]", cleaned_seqs)))  # X and Z should be removed
  }

})
