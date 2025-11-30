library(testthat)
library(IDPsBio) 

library(testthat)
library(IDPsBio)

test_that("write_combined_idpr_results rejects built-in S4 data sets", {
  # load the dataset that is known to be S4
  data("Tau_Human_Protein", package = "IDPsBio")
  
  # confirm it's S4
  expect_true(isS4(Tau_Human_Protein))
  
  results <- list(
    good1 = 1.23,
    bad = Tau_Human_Protein
  )
  
  expect_error(
    write_combined_idpr_results(
      results = results,
      out_file = tempfile(fileext = ".xlsx"),
      file_type = "xlsx",
      expected_functions = NULL
    ),
    "Invalid value detected"
  )
})
