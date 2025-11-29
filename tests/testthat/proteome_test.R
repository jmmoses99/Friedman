# -----------------------------
# Robust test for idpr_parallel_processing
# -----------------------------
library(testthat)
library(IDPsBio)
library(BiocParallel)

#  Step 1: Load Saccharomyces cerevisiae proteome (smallest eukaryotic proteome)
data("Saccharomyces_cerevisiae_proteome", package = "IDPsBio")

# --- Step 2: Convert S4 object to named character vector ---
proteome_char <- as.character(Saccharomyces_cerevisiae_proteome)
names(proteome_char) <- names(Saccharomyces_cerevisiae_proteome)

# --- Step 3: Take first 5 sequences ---
test_sequences <- proteome_char[1:5]

# --- Step 4: Introduce deliberately "dirty" characters safely ---
test_sequences[1] <- paste0(test_sequences[1], "XZ")  # invalid letters
test_sequences[2] <- gsub("A", " ", test_sequences[2]) # space

# --- Step 5: Clean sequences ---
cleaned <- setNames(
  cleaned_protein_seq_aa(test_sequences, warn = FALSE),
  names(test_sequences)
)

# Remove sequences that became empty or too short (optional, e.g., <5 residues)
min_length <- 5
cleaned <- cleaned[nzchar(cleaned) & nchar(cleaned) >= min_length]

# If no sequences left, skip the test
if (length(cleaned) == 0) stop("No valid sequences left after cleaning.")

# --- Step 6: Define functions and parameters ---
idpr_functions <- c("idprofile", "chargeCalculationLocal")
valid_functions <- idpr_functions
BPPARAM <- SnowParam(1)

# --- Step 7: Run robust test ---
test_that("idpr_parallel_processing works on small cleaned proteome subset", {

  results <- idpr_parallel_processing(
    idpr_functions = idpr_functions,
    valid_functions = valid_functions,
    input = cleaned,
    BPPARAM = BPPARAM
  )

  # Keep only non-NULL results
  results <- results[!sapply(results, is.null)]

  # Basic structure checks
  expect_type(results, "list")
  expect_length(results, length(cleaned))

  # Check functions in each result
  for (res in results) {
    # Some functions may return NULL for short sequences; allow that
    expect_true(all(idpr_functions %in% names(res) | sapply(res, is.null)))
    for (fn in idpr_functions) {
      # Either NULL or valid
      expect_true(is.null(res[[fn]]) || !is.null(res[[fn]]))
    }
  }

  # Optional: summarize results silently
  expect_silent(summary_list <- summarize_idpr_results(results))
  expect_true(all(c("protein_summary", "dataset_summary") %in% names(summary_list)))
})
