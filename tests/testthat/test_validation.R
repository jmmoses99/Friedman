library(testthat)
library(Biostrings)

#Test Code Real Example; Important because this utility is critical to idpr_parallel_processing

valid_functions <- c("idprofile", "iupred", "iupredAnchor", "iupredRedox",
                     "chargeCalculationLocal", "chargeCalculationGlobal",
                     "foldIndexR", "scaledHydropathyLocal", "meanScaledHydropathy")

# Need to pull accession # from UniProt
accessions <- c(
  "P10636", #Tau Human Protein - Known IDP
  "P68871", #Hemoglobin subunit beta
  "P04637"  #P53 Human Protein - Known IDP
)

#Downloading real FASTA Files
dir.create("fasta_test", showWarnings = FALSE) # FIXED: avoid warning if folder exists

for(access in accessions){  # each accession temporarily stored in access
  url <- paste0("https://www.uniprot.org/uniprot/", access, ".fasta") # builds download link for fasta of each protein
  destination_file <- file.path("fasta_test", paste0(access, ".fasta")) # builds file path
  download.file(url,destination_file, quiet = TRUE ) #downloads file from the url and saves it to the file path
}

#Reads FASTA Files into a data.frame and returns
read_fasta <- function(fasta_files){
  rows <- lapply(fasta_files, function(file){
    sequ <- readAAStringSet(file)
    #Build a small data.frame for file
    data.frame(
      accession = sub("^>", "", names(sequ)),# sub() finds the first occurrence of a pattern in a string and replaces it
      sequence = as.character(sequ),
      file = basename(file),
      stringsAsFactors = FALSE
    )
  })
  #Combine into a singular data.frame
  do.call(rbind,rows)
}

fasta_files <- list.files("fasta_test", full.names = TRUE) #storing fasta_test in fasta_files
input_df <- read_fasta(fasta_files) #Creates data.frame from fasta files

# Test 1: Validate FASTA input
test_that("FASTA file validation works", {
  for (file in fasta_files) {
    expect_silent( # FIXED: use expect_silent() to assert no errors occur
      validate_idpr_functions(
        idpr_functions = valid_functions,
        valid_functions = valid_functions,
        input = input_df,                #  provide data.frame input
        accessions = input_df$accession, #  provide accessions
        sequences = input_df$sequence    #  provide sequences
      )
    )
  }
  print("FASTA file validation passed") # print inside test_that()
})

# Test 2: Validate data.frame input
test_that("Data.frame validation works", {
  expect_silent(
    validate_idpr_functions(
      idpr_functions = valid_functions,
      valid_functions = valid_functions,
      accessions = input_df$accession,
      sequences = input_df$sequence,
      input = input_df
    )
  )
  print("Data.frame validation passed") # print inside test_that()
})

# FIXED: Clean up folder after test run
unlink("fasta_test", recursive = TRUE)
