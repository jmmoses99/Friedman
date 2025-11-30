# Script: prepare_datasets.R
# Purpose: Convert raw FASTA datasets into package-ready R objects

library(Biostrings)  # for AAStringSet
library(usethis)     # for use_data()

# Set path to your folder with FASTA files
fasta_dir <- "fasta_test"


# 1. P53 Human Protein (UniProt P04637)

p53_file <- file.path(fasta_dir, "P04637.fasta")
P53_Human_Protein <- readAAStringSet(p53_file)
names(P53_Human_Protein) <- sub(" .*", "", names(P53_Human_Protein))
use_data_raw(P53_Human_Protein, overwrite = TRUE)


# 2. Hemoglobin Subunit Beta (UniProt P68871)

hbb_file <- file.path(fasta_dir, "P68871.fasta")
Hemoglobin_subunit_beta <- readAAStringSet(hbb_file)
names(Hemoglobin_subunit_beta) <- sub(" .*", "", names(Hemoglobin_subunit_beta))
use_data(Hemoglobin_subunit_beta, overwrite = TRUE)

#
# 3. Tau Human Protein (UniProt P10636)

tau_file <- file.path(fasta_dir, "P10636.fasta")
Tau_Human_Protein <- readAAStringSet(tau_file)
names(Tau_Human_Protein) <- sub(" .*", "", names(Tau_Human_Protein))
use_data(Tau_Human_Protein, overwrite = TRUE)


# 4. Saccharomyces cerevisiae Proteome (UP000002311)

yeast_file <- file.path(fasta_dir, "UP000002311_proteome.fasta")
Saccharomyces_cerevisiae_proteome <- readAAStringSet(yeast_file)
names(Saccharomyces_cerevisiae_proteome) <- sub(" .*", "", names(Saccharomyces_cerevisiae_proteome))
use_data(Saccharomyces_cerevisiae_proteome, overwrite = TRUE)


# 5. Small subset for testing (first 3 proteins)

Saccharomyces_cerevisiae_proteome_small <- Saccharomyces_cerevisiae_proteome[1:3]
use_data(Saccharomyces_cerevisiae_proteome_small, overwrite = TRUE)
