getwd()
file.exists("fasta_test/P10636.fasta")

library(Biostrings)

# Example for one FASTA
Tau_Human_Protein <- readAAStringSet("fasta_test/P10636.fasta")
save(Tau_Human_Protein, file = "data/Tau_Human_Protein.rda")

# Repeat for each FASTA
Hemoglobin_subunit_beta <- readAAStringSet("fasta_test/P68871.fasta")
save(Hemoglobin_subunit_beta, file = "data/Hemoglobin_subunit_beta.rda")

P53_Human_Protein <- readAAStringSet("fasta_test/P04637.fasta")
save(P53_Human_Protein, file = "data/P04637.rda")

Saccharomyces_cerevisiae_proteome <- readAAStringSet("fasta_test/UP000002311_proteome.fasta")
save(Saccharomyces_cerevisiae_proteome, file = "data/Saccharomyces_cerevisiae_proteome.rda", compress = "xz")
