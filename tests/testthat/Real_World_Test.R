# Real_World_Test.R
# Testing all Functions with a Real World Data Set from UniProt
# Required packages: BiocParallel, Biostrings, idpr, openxlsx

library(testthat)
library(Biostrings)
library(BiocParallel)
library(idpr)
library(openxlsx)


# IMPORTANT — Directory (must exist) otherwise workers may stop early

R_path <- file.path(getwd(), "R")


# Paths to utility scripts and other functions

cleaned_path   <- file.path(R_path, "Utility_Clean_AA.R")
Detect_path    <- file.path(R_path, "Utility_Detect_macOS.R")
validate_path  <- file.path(R_path, "Utility_Validation_of_Functions.R")
parallel_path  <- file.path(R_path, "parallel_processing_backend.R")
write_path     <- file.path(R_path,"write_combined_idpr_results.R")
summarize_path <- file.path(R_path, "summarize_idpr_results.R")


# Load utility scripts

#These source() calls now check existence first

for (p in c(cleaned_path, Detect_path, validate_path, parallel_path, write_path, summarize_path)) {
  if (!file.exists(p)) stop("Required R script not found: ", p)
  source(p)  # Restored proper sourcing order
}

# Retrieve Proteome from UniProt

proteome_id <- "UP000002311"
url <- paste0(
  "https://rest.uniprot.org/uniprotkb/stream?",
  "query=proteome:", proteome_id,
  "&format=fasta&compressed=false&includeIsoform=true"
)
fasta_file <- paste0(proteome_id, "_proteome.fasta")

# download.file now runs only if file not already present

if (!file.exists(fasta_file)) {
  download.file(url, fasta_file, quiet = TRUE)
}

UP000002311_proteome_fasta <- readAAStringSet(fasta_file)


# Prepare input for parallel processing

input_vector <- setNames(
  as.character(UP000002311_proteome_fasta),
  names(UP000002311_proteome_fasta)
)


# Clean sequences to remove invalid residues

# Using the package version of cleaned_protein_seq_aa()

input_vector <- cleaned_protein_seq_aa(input_vector, warn = TRUE)


# Define IDPR functions

idpr_functions <- c(
  "idprofile", "iupred", "iupredAnchor", "iupredRedox",
  "chargeCalculationLocal", "chargeCalculationGlobal",
  "foldIndexR", "scaledHydropathyLocal", "meanScaledHydropathy"
)
valid_functions <- idpr_functions


# Set up parallel backend

BPPARAM <- parallel_backend()  # FIX: Corrected to use your backend function
print(BPPARAM)


# Use tempdir() for output tests

tmp_dir <- tempdir()
message("TEMP DIRECTORY USED FOR TESTS: ", tmp_dir)  # FIX: debug print

-
# Run tests

test_that("idpr_parallel_processing completes without error", {

  results <- NULL

  result_try <- try({
    results <<- idpr_parallel_processing(
      idpr_functions,
      valid_functions,
      input_vector,
      BPPARAM = BPPARAM
    )
  }, silent = TRUE)

  expect_false(inherits(result_try, "try-error"),
               info = paste("Parallel processing failed:", attr(result_try, "condition")$message))

  expect_type(results, "list")
  expect_true(length(results) > 0)

  first <- results[[1]]

  expect_true("accession" %in% names(first))
  expect_true(all(idpr_functions %in% names(first)))

})

# Test writing combined results

test_that("write_combined_idpr_results writes output", {

  results <- idpr_parallel_processing(
    idpr_functions,
    valid_functions,
    input_vector,
    BPPARAM = BPPARAM
  )

  # Write into temp dir instead of package directory

  write_combined_idpr_results(results, dir = tmp_dir)

  out_file <- file.path(tmp_dir, "IDPsBio_Combined_Results.xlsx")
  expect_true(file.exists(out_file))

})


# Test summarize results

test_that("summarize_idpr_results writes summary files", {

  results <- idpr_parallel_processing(
    idpr_functions,
    valid_functions,
    input_vector,
    BPPARAM = BPPARAM
  )

  # Summary  written into temp dir

  summarize_idpr_results(results, dir = tmp_dir)

  summary_file <- file.path(tmp_dir, "IDPsBio_Summary.csv")
  expect_true(file.exists(summary_file))

})

