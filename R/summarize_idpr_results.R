#' Summarize IDPR Results
#'
#' Takes the output list from `idpr_parallel_processing()` and produces summarized
#' tables containing key metrics for each protein and for the entire dataset.
#'
#' @param results List. Output from `idpr_parallel_processing()`.
#' @param protein_outfile Character. Path to save per-protein summary (optional).
#' @param dataset_outfile Character. Path to save dataset summary (optional).
#' @param file_type Character. File type to save ("csv" or "tsv"). Default is "csv".
#'
#' @return A list containing:
#' \describe{
#'   \item{protein_summary}{Data frame with per-protein metrics.}
#'   \item{dataset_summary}{Data frame summarizing the entire dataset.}
#' }
#'
#' @examples
#' sequences <- c(
#'   "P12345" = "MKTFFVAGA",
#'   "Q9XYZ1" = "ACDEFGHIKLMNPQRSTVWY"
#' )
#' idpr_functions <- c("idprofile", "chargeCalculationLocal")
#' valid_functions <- idpr_functions
#' \dontrun{
#'   if (requireNamespace("BiocParallel", quietly = TRUE)) {
#'     BPPARAM <- BiocParallel::SnowParam(1)
#'     results <- idpr_parallel_processing(
#'       idpr_functions = idpr_functions,
#'       valid_functions = valid_functions,
#'       input = sequences,
#'       BPPARAM = BPPARAM
#'     )
#'     summarize_idpr_results(
#'       results,
#'       protein_outfile = "protein_summary.csv",
#'       dataset_outfile = "dataset_summary.csv"
#'     )
#'   }
#' }
#'
#' @importFrom stats median sd
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
