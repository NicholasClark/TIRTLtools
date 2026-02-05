## code to prepare `vdj_db` dataset goes here

## downloaded from https://github.com/antigenomics/vdjdb-db/releases/tag/pyvdjdb-2025-09-25
# library(dplyr)
vdj_db = data.table::fread("inst/extdata/vdjdb/vdjdb_full_2025_autumn.txt") %>%
  dplyr::filter(
    species == "HomoSapiens", cdr3.alpha != "", cdr3.beta != "",
    v.alpha != "", j.alpha != "", j.beta != "", v.beta != "",
    antigen.species != "HomoSapiens" ## TCRs are too similar, makes networks huge
    ) %>%
  as_tibble() %>%
  dplyr::rename(cdr3a = cdr3.alpha, va = v.alpha, ja = j.alpha, cdr3b = cdr3.beta, vb = v.beta, db = d.beta,
                jb = j.beta) %>%
  mutate(tcr = paste(va,vb,cdr3a,cdr3b)) %>%
  filter(!duplicated(tcr)) %>%
  select(-tcr)

usethis::use_data(vdj_db, overwrite = TRUE)
