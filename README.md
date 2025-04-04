# TIRTLtools R package

A suite of tools for analyzing paired T-cell receptor (TCR) datasets created using TIRTL-seq.

Package may change frequently as it is in active development.

## Installation

You can install the package from GitHub:
```R
remotes::install_github("NicholasClark/TIRTLtools")
```
or you can clone the repository to a local folder
```
git clone https://github.com/NicholasClark/TIRTLtools.git
```
and load it using devtools:
```R
devtools::load_all("./TIRTLtools") ### location of cloned repo on your machine
```

## Quick start

### TCRdist GPU implementation

This function computes the similarity of two T-cell receptors (TCRs) using TCRdist ([Dash et al., Nature 2017](https://doi.org/10.1038/nature22383)). The function can efficiently calculate pairwise similarities of thousands or millions of TCRs using a GPU if it is available. It relies on python code using 'cupy' (Nvidia GPUs) or 'mlx' (Apple Silicon GPUs), which is called through the 'reticulate' R package.

The following example computes TCRdist values for all pairs of annotated TCRs in the [VDJdb](https://vdjdb.cdr3.net/) database. It outputs all values with TCRdist <= 90.

```R
#library(TIRTLtools)
devtools::load_all("./TIRTLtools") ### location of cloned repo on your machine
library(dplyr)

### Load an example dataset (annotated TCRs from VDJ-db)
tcr1 = TIRTLtools::vdj_db
### note: You may replace this with a file of your choice - needs to have columns "va", "vb", "cdr3a", and "cdr3b"

as_tibble(tcr1)
# # A tibble: 26,213 × 15
#       complex.id cdr3a         va       ja    cdr3b vb    jb    Species MHC.A MHC.B MHC.class Epitope Epitope.gene
#         <int>   <chr>         <chr>    <chr> <chr> <chr> <chr> <chr>   <chr> <chr> <chr>     <chr>   <chr>
# 1          5 CAVRDGGTGFQKLVF   TRAV3*01 TRAJ… CASR… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      HPVGEA… EBNA1
# 2          6 CAARGIGSGTYKYIF   TRAV13-… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      HPVGEA… EBNA1
# 3          8 CAVRDGTANNLFF     TRAV3*01 TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 4          9 CAFMKDAGGTSYGKLTF TRAV38-… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 5         10 CAGGGSQGNLIF      TRAV27*… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 6         11 CALAGSQGNLIF      TRAV6*01 TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 7         12 CAGADGGSQGNLIF    TRAV27*… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 8         13 CAAGGSQGNLIF      TRAV27*… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 9         14 CAAPPSNTGKLIF     TRAV3*01 TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 10         15 CAVRDLLTNAGKSTF   TRAV3*01 TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# # ℹ 26,203 more rows
# # ℹ 2 more variables: Epitope.species <chr>, Reference <chr>
# # ℹ Use `print(n = ...)` to see more rows
result = TCRdist(tcr1 = tcr1, tcrdist_cutoff = 90, chunk_size = 5000)
# Checking for available GPU...
#
# Apple Silicon GPU detected:
#   Apple Silicon GPU (M1/M2/M3)
# Checking for GPU-related Python modules...
#
# 'mlx' is installed (for Apple Silicon GPUs).
# Loading mlx to perform TCRdist
# total number of chunks (rows): 5
# total number of chunks (cols): 5
# Processing chunk (rows) 0
# Processing chunk (rows) 5000
# Processing chunk (rows) 10000
# Processing chunk (rows) 15000
# Processing chunk (rows) 20000
# Processing chunk (rows) 25000
# Time taken: 2.970693 seconds

edge_df = result[['TCRdist_df']] ### table of TCRdist values <= cutoff
meta_df = result[['tcr1']] ### table of input data with indices

as_tibble(edge_df)
# # A tibble: 101,797 × 3
#       edge1_0index edge2_0index TCRdist
#           <int>        <int>   <int>
# 1            5            4      87
# 2            7            4       3
# 3            7            5      90
# 4           12            4      15
# 5           12            5      78
# 6           12            7      12
# 7           13            4      78
# 8           13            6      72
# 9           13            7      81
# 10           14            6      45
# # ℹ 101,787 more rows
# # ℹ Use `print(n = ...)` to see more rows

as_tibble(meta_df)
# # A tibble: 26,213 × 16
#   complex.id cdr3a      va    ja    cdr3b vb    jb    Species MHC.A MHC.B MHC.class Epitope Epitope.gene
#       <int>   <chr>     <chr> <chr> <chr> <chr> <chr> <chr>   <chr> <chr> <chr>     <chr>   <chr>
# 1          5 CAVRDGGTG… TRAV… TRAJ… CASR… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      HPVGEA… EBNA1
# 2          6 CAARGIGSG… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      HPVGEA… EBNA1
# 3          8 CAVRDGTAN… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 4          9 CAFMKDAGG… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 5         10 CAGGGSQGN… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 6         11 CALAGSQGN… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 7         12 CAGADGGSQ… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 8         13 CAAGGSQGN… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 9         14 CAAPPSNTG… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 10         15 CAVRDLLTN… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# # ℹ 26,203 more rows
# # ℹ 3 more variables: Epitope.species <chr>, Reference <chr>, tcr_index <dbl>
# # ℹ Use `print(n = ...)` to see more rows

#' TCRdist function output:
#' A list with two items:
#'    'TCRdist_df' - A dataframe containing all TCRdist values less than the cutoff (default = 90).
#'                  Contains the indices of the two TCRs and the TCRdist value.
#'    'tcr1' - A dataframe with the input data, along with a column "tcr_index" with the index of each TCR
```

## References

Please see the [TIRTLseq preprint (Pogorelyy and Kirk et al.)](https://www.biorxiv.org/content/10.1101/2024.09.16.613345v1) for details.
