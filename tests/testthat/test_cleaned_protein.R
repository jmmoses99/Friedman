#Test Code for cleaning protein sequences; Important to ensure only valid amino acids are kept

library(testthat)
library(Biostrings)

# Use FASTA files from previous example

fasta_files <- list.files("fasta_test", full.names = TRUE) # storing fasta_test in fasta_files

#Reads FASTA Files into a data.frame and returns

read_fasta <- function(fasta_files){

  rows <- lapply(fasta_files, function(file){
    sequ <- readAAStringSet(file)
    #Build a small data.frame for file
    data.frame(
      accession = sub("^>", "", names(sequ)),# sub() finds the first occurrence of a pattern in a string and replaces it
      sequence = as.character(sequ),          # convert AAStringSet to character
      file = basename(file),
      stringsAsFactors = FALSE
    )
  })

  #Combine into a singular data.frame

  do.call(rbind,rows)
}

# Create data.frame from FASTA files

fasta_df <- read_fasta(fasta_files)

#Test 1: Ensure sequences are cleaned without errors

expect_silent(
  IDPsBio:::cleaned_protein_seq_aa(as.character(fasta_df$sequence)) #convert to character before passing
)
print("FASTA sequences cleaned correctly")

#Test 2: Ensure sequences with invalid residues trigger warning

seqs <- c("MKT*Y$V", "GHIJKL", "ARNDCEQ") # deliberately add invalid residues
expect_warning(
  IDPsBio:::cleaned_protein_seq_aa(seqs)  # warning expected
)
print("Sequences with invalid residues are cleaned with warning")
