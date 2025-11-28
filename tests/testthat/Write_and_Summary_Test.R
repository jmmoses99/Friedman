# Test: Writing and summarizing idpr_parallel_processing results

# Assumes the earlier test created `results` object and wrote:
#   combined_idpr_results.xlsx
#   summarized_idpr_results.xlsx

# These are written inside a temporary directory created by testthat

library(testthat)
library(openxlsx)

test_that("IDPR results were written correctly", {


  # TEMP DIRECTORY — created earlier during the pipeline test
  # testthat::with_tempdir() ensures all files stay inside temp space

  temp_dir <- tempdir()

  combined_path  <- file.path(temp_dir, "combined_idpr_results.xlsx")
  summary_path   <- file.path(temp_dir, "summarized_idpr_results.xlsx")


  # 1. Check that output files exist

  expect_true(file.exists(combined_path))
  expect_true(file.exists(summary_path))


  # 2. Load written data

  combined_df <- openxlsx::read.xlsx(combined_path)
  summary_df  <- openxlsx::read.xlsx(summary_path)


  # 3. Basic structure sanity checks

  # MUST contain accession column
  expect_true("accession" %in% names(combined_df))

  # All IDPR functions must appear in combined output
  expected_cols <- c(
    "accession",
    "idprofile", "iupred", "iupredAnchor", "iupredRedox",
    "chargeCalculationLocal", "chargeCalculationGlobal",
    "foldIndexR", "scaledHydropathyLocal", "meanScaledHydropathy"
  )
  expect_true(all(expected_cols %in% names(combined_df)))

  # Summary should contain fewer columns

  expect_true("accession" %in% names(summary_df))
  expect_gt(ncol(summary_df), 1)   # at least one summary metric


  # Row count consistency
  # Must match the number of proteins processed

  expect_equal(nrow(combined_df), length(results))
  expect_equal(nrow(summary_df), length(results))

  # Check for absence of missing accessions

  expect_false(any(is.na(combined_df$accession)))
  expect_false(any(is.na(summary_df$accession)))


 # Sample-based sanity check: first 5 rows must contain valid results

  sample_rows <- seq_len(min(5, nrow(combined_df)))

  for (i in sample_rows) {
    expect_true(!is.na(combined_df$accession[i]))
    expect_true(!is.na(summary_df$accession[i]))
  }
})

