# Parallelized C++ implementation of TCRdist (no GPU required)

**\[experimental\]** This is an alternative to the GPU version of
TCRdist that is still very fast for large datasets (tens of thousands of
TCRs). It is written in C++ and will run in parallel across available
CPU cores.

## Usage

``` r
TCRdist_cpp(tcr1, tcr2 = NULL)
```

## Arguments

- tcr1:

  a data frame with one TCR per row. It must have the columns "va",
  "vb", "cdr3a", and "cdr3b"

- tcr2:

  (optional) another data frame of TCRs. If supplied, TCRdist will be
  calculated for every combination of one TCR from tcr1 and one TCR from
  tcr2. Otherwise, it will calculate TCRdist for every pair of TCRs in
  tcr1.

## Value

a list with two objects (or three if tcr2 is not null):

- matrix - a matrix of TCRdist values

- tcr1 - the input matrix tcr1, after pre-processing and removing
  unacceptable TCRs

- tcr2 (if supplied) - the input matrix tcr2, after pre-processing and
  removing unacceptable TCRs

## Details

This version of TCRdist is currently less feature-rich than
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
and returns a dense matrix as output. It does not yet allow for sparse
output or writing output directly to a file.

## See also

Other tcr_similarity:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md)

## Examples

``` r
folder = system.file("extdata/SJTRC_TIRTLseq_minimal",
  package = "TIRTLtools")
sjtrc = load_tirtlseq(folder,
  meta_columns = c("marker", "timepoint", "version"), sep = "_",
  chain = "paired", verbose = FALSE)
#> Loading files from: /Users/nclark52/git/TIRTLtools/inst/extdata/SJTRC_TIRTLseq_minimal...
#> Found 2 beta chain pseudo-bulk files.
#> Found 2 paired chain files.
#> Loaded 2 files from 2 samples.
#> 0.2 seconds
df = get_all_tcrs(sjtrc, chain="paired", remove_duplicates = TRUE)

result = TCRdist_cpp(df)
#> Removed 384 TCRs with unknown V-segments (1.2%) from a total of 32,164 TCRs.
#> Removed 12 TCRs with short CDR3 segments (0.038%) from a total of 31,780 TCRs.
#> Removed 13,324 TCRs with non-functional CDR3 amino acid sequences (42%) from a total of 31,768 TCRs.
#> Removed 812 MAIT TCRs (4.4%) from a total of 18,444 TCRs.
#> Filtered data frame contains 17,632 TCRs (55%) of original 32,164 TCRs.
#> Removed 384 TCRs with unknown V-segments (1.2%) from a total of 32,164 TCRs.
#> Removed 12 TCRs with short CDR3 segments (0.038%) from a total of 31,780 TCRs.
#> Removed 13,324 TCRs with non-functional CDR3 amino acid sequences (42%) from a total of 31,768 TCRs.
#> Removed 812 MAIT TCRs (4.4%) from a total of 18,444 TCRs.
#> Filtered data frame contains 17,632 TCRs (55%) of original 32,164 TCRs.

mat = result$matrix
node_df = result$tcr1

mat[1:5,1:5]
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   NA   NA   NA   NA   NA
#> [2,]  338   NA   NA   NA   NA
#> [3,]  287  270   NA   NA   NA
#> [4,]  374  320  293   NA   NA
#> [5,]  382  372  386  387   NA
## note: TCRdist is symmetric. Returned matrix contains only lower-triangular values

node_df %>%
  mutate(alpha_nuc = paste(substr(alpha_nuc, 0, 20), "...", sep = ""),
         beta_nuc = paste(substr(beta_nuc, 0, 20), "...", sep = "")) %>%
  data.table::as.data.table()
#> Error in mutate(., alpha_nuc = paste(substr(alpha_nuc, 0, 20), "...",     sep = ""), beta_nuc = paste(substr(beta_nuc, 0, 20), "...",     sep = "")): could not find function "mutate"

```
