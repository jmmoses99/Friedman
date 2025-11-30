# Test: Writing and summarizing idpr_parallel_processing results

library(testthat)
library(openxlsx)

test_that("IDPR results were written correctly", {

  # TEMP DIRECTORY — created earlier during the pipeline test
  temp_dir <- tempdir()

  combined_path  <- file.path(temp_dir, "combined_idpr_results.xlsx")
  summary_path   <- file.path(temp_dir, "summarized_idpr_results.xlsx")

  # --- Create a minimal results object for testing ---
  results <- list(
    list(
      accession = "P12345",
      idprofile = NA,
      iupred = NA,
      iupredAnchor = NA,
      iupredRedox = NA,
      chargeCalculationLocal = data.frame(
        Position = 1,
        CenterResidue = "A",
        Window = "ACDEFGHIK",
        windowCharge = 0.1
      ),
      chargeCalculationGlobal = NA,
      foldIndexR = NA,
      scaledHydropathyLocal = NA,
      meanScaledHydropathy = NA
    ),
    list(
      accession = "Q67890",
      idprofile = NA
      # other columns can be NA
    )
  )

  expected_cols <- c(
    "accession", "idprofile", "iupred", "iupredAnchor", "iupredRedox",
    "chargeCalculationLocal", "chargeCalculationGlobal",
    "foldIndexR", "scaledHydropathyLocal", "meanScaledHydropathy"
  )

  # --- Write combined results ---
  write_combined_idpr_results(
    results = results,
    out_file = combined_path,
    file_type = "xlsx",
    expected_functions = expected_cols[-1]  # exclude "accession"
  )

  # --- Check that output file exists ---
  expect_true(file.exists(combined_path))

  # --- Read and test basic structure ---
  combined_df <- openxlsx::read.xlsx(combined_path)
  expect_true("accession" %in% names(combined_df))
  expect_true(all(expected_cols %in% names(combined_df)))

})
