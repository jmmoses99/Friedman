# Run_all_tests.R
# Testing all functions with a small subset of the UniProt yeast proteome
# Required packages: BiocParallel, Biostrings, idpr, testthat
# Run_all_tests.R
# Self-contained testing of IDPsBio functions
# Works on any machine without external FASTA


library(testthat)
library(Biostrings)
library(BiocParallel)
library(IDPsBio)

# -------------------------------
# 1. Load required utility scripts
# -------------------------------
R_path <- file.path(getwd(), "R")
required_scripts <- c(
  "Utility_Detect_macOS.R",
  "Utility_Validation_of_Functions.R",
  "parallel_processing_backend.R"
)

for (f in required_scripts) {
  f_path <- file.path(R_path, f)
  if (!file.exists(f_path)) stop("Required R script not found: ", f_path)
  source(f_path)
}

# -------------------------------
# 2. Create or load proteome subset
# -------------------------------
subset_file <- file.path("fasta_test", "UP000002311_proteome_subset.fasta")

if (!file.exists(subset_file)) {
  message("Subset FASTA not found — generating automatically.")

  # Load full proteome from your package dataset
  data("Saccharomyces_cerevisiae_proteome", package = "IDPsBio")

  # Select a random subset of 10 sequences
  set.seed(123)
  subset_idx <- sample(seq_along(Saccharomyces_cerevisiae_proteome), 10)
  subset_obj <- Saccharomyces_cerevisiae_proteome[subset_idx]

  # Ensure output directory exists
  if (!dir.exists("fasta_test")) dir.create("fasta_test")

  # Save subset as FASTA
  writeXStringSet(subset_obj, filepath = subset_file)
  message("Subset FASTA created at: ", subset_file)
}

# Load subset
proteome_subset <- readAAStringSet(subset_file)
names(proteome_subset) <- sub(" .*", "", names(proteome_subset))

# Convert to named character vector
input_vector <- as.character(proteome_subset)
names(input_vector) <- names(proteome_subset)

# -------------------------------
# 3. Clean sequences
# -------------------------------
input_vector <- cleaned_protein_seq_aa(input_vector, warn = TRUE)
input_vector <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", input_vector)

# -------------------------------
# 4. Set up parallel backend
# -------------------------------
BPPARAM <- parallel_backend()  # SnowParam workers load packages automatically

# -------------------------------
# 5. Define IDPR functions to test
# -------------------------------
idpr_functions <- c(
  "idprofile", "iupred", "iupredRedox", "chargeCalculationLocal",
  "chargeCalculationGlobal", "foldIndexR", "scaledHydropathyLocal",
  "meanScaledHydropathy"
)
valid_functions <- idpr_functions  # excludes iupredAnchor for safety

# -------------------------------
# 6. Run tests
# -------------------------------
tmp_dir <- tempdir()
message("Temporary directory for results: ", tmp_dir)

test_that("idpr_parallel_processing completes safely", {
  results <- NULL
  result_try <- try({
    results <<- idpr_parallel_processing(
      idpr_functions,
      valid_functions,
      input_vector,
      BPPARAM = BPPARAM
    )
  }, silent = TRUE)

  expect_false(inherits(result_try, "try-error"))
  expect_type(results, "list")
  expect_true(length(results) > 0)

  first <- results[[1]]
  expect_true("accession" %in% names(first))
  expect_true(all(valid_functions %in% names(first)))
})

test_that("write_combined_idpr_results writes output safely", {
  results <- idpr_parallel_processing(
    idpr_functions,
    valid_functions,
    input_vector,
    BPPARAM = BPPARAM
  )

  write_combined_idpr_results(results)

  out_file <- file.path(getwd(), "IDPsBio_Combined_Results.xlsx")
  expect_true(file.exists(out_file))
})

test_that("summarize_idpr_results writes summary safely", {
  results <- idpr_parallel_processing(
    idpr_functions,
    valid_functions,
    input_vector,
    BPPARAM = BPPARAM
  )

  summarize_idpr_results(results)

  summary_file <- file.path(getwd(), "IDPsBio_Summary.csv")
  expect_true(file.exists(summary_file))
})
