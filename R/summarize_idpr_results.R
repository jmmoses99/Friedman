#'Function to help summarize all protein computation results to include a protein-level summary and a Data set Level summary
#' writes it into Summary files

#' This function takes the output list from \code{idpr_parallel_processing()}
#' and produces a summarized table containing key metrics for a protein-level summary of each protein
#' and a data-set level summary, where the entire data set is taken into account.

#' Aggregates multiple IDPR computation outputs into summary tables for further analysis.
#'
#' @param results A list of results from `idpr_parallel_processing()`.
#' @param protein_outfile Optional file path to save per-protein summary results as CSV.
#' @param dataset_outfile Optional file path to save dataset-level summary results as CSV.
#' @param file_type Character, format to write files, e.g., "csv".
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{protein_summary}{Data frame with per-protein summary statistics.}
#'     \item{dataset_summary}{Data frame with dataset-level summary statistics.}
#'   }
#'

#' @examples
#' \donttest{
#' sequences <- c(
#'   "P12345" = "MKTFFVAGA",
#'   "Q9XYZ1" = "ACDEFGHIKLMNPQRSTVWY"
#' )
#'
#' idpr_functions <- c("idprofile", "chargeCalculationLocal")
#' valid_functions <- idpr_functions
#'
#' if (requireNamespace("BiocParallel", quietly = TRUE)) {
#'   BPPARAM <- BiocParallel::SnowParam(1)
#'   results <- idpr_parallel_processing(
#'       idpr_functions = idpr_functions,
#'       valid_functions = valid_functions,
#'       input = sequences,
#'       BPPARAM = BPPARAM
#'   )
#'   str(results)
#' }
#' }

#' @export




summarize_idpr_results <- function(results,
                                 protein_outfile = "protein-level_summary.xlsx",
                                 dataset_outfile = "dataset-level_summary.xlsx",
                                 file_type = "xlsx") {

  # Step 1: Convert results to data frame; Converts your parallel results (a list of lists, one per protein)

  results_df <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))

  # Step 2: Identify numeric columns (IDPR function scores)
  # Finds which columns in results_df are numeric (these are the IDPR function outputs).
  # Creates numeric_data containing only the numeric values

  num_cols <- sapply(results_df, is.numeric)
  numeric_data <- results_df[, num_cols, drop = FALSE]

  # Step 3: Protein-level summary:
  # Quick Summary of each protein by either the accession # from UniProt or the name;
  # Calculates mean of all IDPR function scores for that protein
  # Calculates maximum IDPR function score
  # number of missing values for that protein

  protein_summary <- data.frame(
    accession = results_df$accession,
    name = results_df$name,
    mean_score = rowMeans(numeric_data, na.rm = TRUE),
    max_score = apply(numeric_data, 1, max, na.rm = TRUE),
    n_missing = apply(numeric_data, 1, function(x) sum(is.na(x))),
    stringsAsFactors = FALSE
  )

  # Sort by accession if available, else by the name of the protein

  if ("accession" %in% colnames(protein_summary)) {
    protein_summary <- protein_summary[order(protein_summary$accession), ]
  } else {
    protein_summary <- protein_summary[order(protein_summary$name), ]
  }

  # Step 4: Data set-level summary for each idpr functions across all proteins
  # 2 is the logic that allows the statistic to be calculate column wise


  dataset_summary <- data.frame(
    function_name = colnames(numeric_data),
    mean = colMeans(numeric_data, na.rm = TRUE),
    median = apply(numeric_data, 2, median, na.rm = TRUE),
    sd = apply(numeric_data, 2, sd, na.rm = TRUE),
    min = apply(numeric_data, 2, min, na.rm = TRUE),
    max = apply(numeric_data, 2, max, na.rm = TRUE),
    n_missing = apply(numeric_data, 2, function(x) sum(is.na(x))),
    stringsAsFactors = FALSE
  )

  # Step 5: Write summaries into two files for protein-level summary and a Data set Level summary

  write_combined_idpr_results(protein_summary, protein_outfile, file_type, colnames(protein_summary))
  write_combined_idpr_results(dataset_summary, dataset_outfile, file_type, colnames(dataset_summary))

  invisible(list(protein_summary = protein_summary, dataset_summary = dataset_summary))
}
