#' Function to combine all protein computation results in one file (.csv or .xlsx) organized by each protein (row)
#' and each function result becomes a column
#'
#' Writes the aggregated results from multiple IDPR functions to disk.
#'
#' @param results List of results from `idpr_parallel_processing()`.
#' @param out_file Character, path to output file.
#' @param file_type Character, format to save ("csv" or "xlsx").
#' @param expected_functions Character vector of IDPR functions expected in the results.
#'
#' @return Invisibly returns the path to the written file.
#'
#' @examples
#' \donttest{
#' # Minimal fake results mimicking actual output
#' results <- list(
#'   list(
#'     accession = "P12345",
#'     idprofile = NA,
#'     chargeCalculationLocal = data.frame(
#'       Position = 5,
#'       CenterResidue = "F",
#'       Window = "MKTFFVAGA",
#'       windowCharge = 0.11
#'     )
#'   )
#' )
#'
#' # Use a temporary file to avoid creating files in the check directory
#' tmp_file <- tempfile(fileext = ".csv")
#' write_combined_idpr_results(
#'   results = results,
#'   out_file = tmp_file,
#'   file_type = "csv",
#'   expected_functions = c("idprofile", "chargeCalculationLocal")
#' )
#' }
#'
#' @export
write_combined_idpr_results <- function(results,
                                        out_file = "combined_idpr_results.csv",
                                        file_type = "csv",
                                        expected_functions = c(
                                          "idprofile", "iupred", "iupredAnchor", "iupredRedox",
                                          "chargeCalculationLocal", "chargeCalculationGlobal",
                                          "foldIndexR", "scaledHydropathyLocal", "meanScaledHydropathy"
                                        )) {

  #  Convert list of lists to data frame
  results_df <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))

  #  Ensure all expected columns exist (fill missing columns with NA)
  expected_cols <- c(intersect(c("accession","name"), colnames(results_df)), expected_functions)
  missing_cols <- setdiff(expected_cols, colnames(results_df))
  if (length(missing_cols) > 0) {
    results_df[missing_cols] <- NA
  }

  #  Dynamic row sorting
  # Sorts rows by accession if it exists; then by name; or the first column
  if ("accession" %in% colnames(results_df)) {
    results_df <- results_df[order(results_df$accession), ]
  } else if ("name" %in% colnames(results_df)) {
    results_df <- results_df[order(results_df$name), ]
  } else {
    results_df <- results_df[order(results_df[[1]]), ]
  }

  # Column ordering
  # Ensures only columns present in the data frame are selected (avoids errors)
  results_df <- results_df[, intersect(expected_cols, colnames(results_df))]

  # Write file type into either a csv or an xlsx
  if (file_type == "csv") {
    utils::write.csv(results_df, out_file, row.names = FALSE)
  } else if (file_type == "xlsx") {
    openxlsx::write.xlsx(results_df, out_file, rowNames = FALSE)
  } else {
    stop("Unsupported file_type. Use 'csv' or 'xlsx'.")
  }

  # Feedback & return
  # Provides a message that combined results were saved and the path they were saved
  message("Combined results saved to: ", normalizePath(out_file))
  invisible(results_df)  # return the combined data frame invisibly
}
