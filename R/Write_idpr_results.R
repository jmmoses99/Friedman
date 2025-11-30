#' Write combined IDPR results
#'
#' This function writes processed IDPR results to a file. It accepts only
#' numeric (single values), NA, or S3 lists. S4, S7, and ggplot objects are
#' rejected. Missing columns or rows are filled with NA to ensure consistent output.
#'
#' @param results A named list of processed results (numeric(1), NA, or S3 list).
#' @param out_file Output file path (xlsx or csv supported).
#' @param file_type File type: "xlsx" or "csv".
#' @param expected_functions Optional character vector of expected function names.
#' @importFrom utils write.csv
#' @export
write_combined_idpr_results <- function(results,
                                        out_file,
                                        file_type = c("xlsx", "csv"),
                                        expected_functions = NULL) {

  file_type <- match.arg(file_type)

  # --- internal validator ---
  validate_value <- function(val, acc = NULL, fname = NULL) {
    if ((is.numeric(val) && length(val) == 1) || (is.na(val) && length(val) == 1)) return(TRUE)
    if (is.list(val) && !isS4(val)) return(TRUE)
    stop(sprintf(
      "Invalid value detected%s%s: Only numeric(1), NA, or S3 list objects are allowed.",
      if (!is.null(acc)) paste0(" for accession '", acc, "'") else "",
      if (!is.null(fname)) paste0(" in function '", fname, "'") else ""
    ))
  }

  # --- validate each element ---
  lapply(results, validate_value)

  # --- convert each result into a row (data.frame) ---
  df_list <- lapply(names(results), function(nm) {
    val <- results[[nm]]

    if (is.null(val)) val <- NA  # replace NULL with NA

    if (is.numeric(val) || is.na(val)) {
      data.frame(accession = nm, value = val, stringsAsFactors = FALSE)
    } else if (is.list(val)) {
      # flatten simple named lists into a single-row data.frame
      df <- as.data.frame(val, stringsAsFactors = FALSE)
      df$accession <- nm
      df
    } else {
      # should never get here due to validation
      stop("Unexpected value type")
    }
  })

  # --- align columns ---
  all_cols <- unique(unlist(lapply(df_list, names)))
  df_list <- lapply(df_list, function(df) {
    missing <- setdiff(all_cols, names(df))
    for (col in missing) df[[col]] <- NA
    df[all_cols]  # reorder columns
  })

  combined_df <- do.call(rbind, df_list)

  # --- write output ---
  if (file_type == "csv") {
    write.csv(combined_df, out_file, row.names = FALSE)
  } else if (file_type == "xlsx") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' required for writing xlsx files")
    }
    openxlsx::write.xlsx(combined_df, out_file, rowNames = FALSE)
  }

  invisible(combined_df)
}
