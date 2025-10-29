## code to prepare `sub_mat` dataset goes here

submat = as.matrix( read.table("inst/extdata/TCRdist_matrix_mega.tsv", sep = "\t") )
usethis::use_data(submat, overwrite = TRUE)
