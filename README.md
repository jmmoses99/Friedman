# Friedman
Blog IDPsBio R package: https://jmmoses99-drjac.wordpress.com/2025/11/30/week-14-idpsbio-r-package/
Need to install Biostrings and idpr. I import idpr, BiocParallel, and openxlsx.
This package provides additional functions for data analysis of intrinsically disordered proteins (IDPs)
         and intrinsically disordered regions (IDPRs). It works on Windows and macOs systems.
         It includes utilities for batch processing, data validation, 
         sequence cleaning, and integration with the 'idpr' package.
         Parallel processing is supported via the BiocParallel framework to
         efficiently handle large proteome-scale data sets.

 Tests can be found in tests/testthat folders. I used the P53_Human_Protein for tests and it should be easy to load the data after installing all packages. I also included other protein data sets for the user to use if they would like.
 
The IDPsBio R package uses both S3 and S4 objects. The protein data sets that I included are S4 objects. They can be read using Biostrings.
