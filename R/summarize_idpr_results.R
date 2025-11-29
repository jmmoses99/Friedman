#' Summarize IDPR Results (Safe Version)
#'
#' Safely summarizes outputs from `idpr_parallel_processing()`, ignoring non-numeric or non-coercible elements.
#'
#' @inheritParams summarize_idpr_results
#' @importFrom stats median sd
#' @export
summarize_idpr_results <- function(results,
                                   protein_outfile = "protein-level_summary.xlsx",
                                   dataset_outfile = "dataset-level_summary.xlsx",
                                   file_type = "xlsx") {

  # Step 0: Helper to safely convert elements to data.frame
  safe_as_df <- function(x) {
    # Keep only numeric, character, or logical elements
    x <- lapply(x, function(el) {
      if (is.numeric(el) || is.character(el) || is.logical(el)) {
        return(el)
      } else {
        return(NA) # replace plots/S4 objects with NA
      }
    })
    as.data.frame(x, stringsAsFactors = FALSE)
  }

  # Step 1: Convert results safely
  results_df_list <- lapply(results, safe_as_df)

  # Step 2: Bind rows; warn if fails
  results_df <- tryCatch(
    do.call(rbind, results_df_list),
    error = function(e) {
      warning("Some elements could not be coerced to data.frame and are replaced with NA.")
      do.call(rbind, lapply(results_df_list, function(x) {
        # replace completely non-coercible rows with NA row
        if (is.null(ncol(x))) x <- as.data.frame(t(rep(NA, length(x))), stringsAsFactors = FALSE)
        x
      }))
    }
  )

  # Step 3: Identify numeric columns
  numeric_data <- results_df[, sapply(results_df, is.numeric), drop = FALSE]

  # Step 4: Protein-level summary
  protein_summary <- data.frame(
    accession = if ("accession" %in% colnames(results_df)) results_df$accession else NA,
    name = if ("name" %in% colnames(results_df)) results_df$name else NA,
    mean_score = rowMeans(numeric_data, na.rm = TRUE),
    max_score = apply(numeric_data, 1, max, na.rm = TRUE),
    n_missing = apply(numeric_data, 1, function(x) sum(is.na(x))),
    stringsAsFactors = FALSE
  )

  # Step 5: Dataset-level summary
  dataset_summary <- data.frame(
    function_name = colnames(numeric_data),
    mean   = colMeans(numeric_data, na.rm = TRUE),
    median = apply(numeric_data, 2, median, na.rm = TRUE),
    sd     = apply(numeric_data, 2, sd, na.rm = TRUE),
    min    = apply(numeric_data, 2, min, na.rm = TRUE),
    max    = apply(numeric_data, 2, max, na.rm = TRUE),
    n_missing = apply(numeric_data, 2, function(x) sum(is.na(x))),
    stringsAsFactors = FALSE
  )

  # Step 6: Save summaries
  write_combined_idpr_results(protein_summary, protein_outfile, file_type, colnames(protein_summary))
  write_combined_idpr_results(dataset_summary, dataset_outfile, file_type, colnames(dataset_summary))

  invisible(list(protein_summary = protein_summary, dataset_summary = dataset_summary))
}
