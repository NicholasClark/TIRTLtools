## code to prepare `params` dataset goes here

#params = readr::read_tsv(fs::path_package("TIRTLtools", "inst", "extdata", package = "params_v2.tsv"))
#params = read.table("inst/extdata/params_v2.tsv", sep = "\t")

params = readr::read_tsv("inst/extdata/params_v2.tsv", col_names = FALSE) %>%
  magrittr::set_colnames(c("feature", "value"))

usethis::use_data(params, overwrite = TRUE)
