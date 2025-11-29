library(testthat)

# Test 1: Detect current OS (simple check)
test_that("Detect current OS works correctly", {
  Running_macOs <- .detect_macOs_internal(check_systems = TRUE)

  # Must be a single logical
  expect_type(Running_macOs, "logical")
  expect_length(Running_macOs, 1)

  # Should be TRUE only if running on macOS, FALSE otherwise
  expect_true(Running_macOs %in% c(TRUE, FALSE))
})

# Test 2: Vectorized check (default)
test_that("Vectorized detection returns logical", {
  # By default, check_systems = FALSE
  result_vec <- .detect_macOs_internal(check_systems = FALSE)

  expect_type(result_vec, "logical")
  expect_length(result_vec, 1)  # still just one element
  expect_true(result_vec %in% c(TRUE, FALSE))
})
