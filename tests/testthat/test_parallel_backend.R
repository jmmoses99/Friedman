# Test parallel_backend function
library(testthat)
library(IDPsBio)

test_that("parallel_backend returns a valid BiocParallelParam", {
  # Function for parallel processing

  BPPARAM <- parallel_backend()

  # Debug: print the backend returned

  print(BPPARAM)

  # Check that the object is the right BiocParallelParam class

  expect_s4_class(BPPARAM, "BiocParallelParam") # Explicitly tests correct object type
})
