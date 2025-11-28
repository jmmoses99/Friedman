# Test 1 code for function (Based on Strings); Not a good test because of the warning message

library(testthat)

# Example strings to test detection
Test_Sys <- c(
  "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)", # For macOs
  "Mozilla/5.0 (Windows NT 10.0; Win64; x64)",       # For Windows
  "Linux x86_64",                                    # For Linux
  "Darwin Kernel Version 20.3.0"                     # For Darwin Operating Systems
)

# Test 1: Detect MacOS from strings
test_that("Detect MacOS works correctly", {
  # Call the internal function directly (loaded by devtools::load_all())
  results <- .detect_macOs_internal(Test_Sys, cores = 2) # Warning message expected on Windows; BioConductor notes cause of SnowParam
  expect_equal(results, c(TRUE, FALSE, FALSE, TRUE))
})

# Test 2: Detect current OS
test_that("Detect current OS returns logical", {
  # Call the internal function directly
  Running_macOs <- .detect_macOs_internal(check_systems = TRUE)
  print(Running_macOs)
  # Results should be TRUE if running on a macOs, FALSE if Windows/Linux
  expect_true(is.logical(Running_macOs))
})
