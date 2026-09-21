# Parallelized C++ implementation of TCRdist (no GPU required)

**\[experimental\]**

This is an alternative to the GPU version of TCRdist that is still very
fast for large datasets (tens of thousands of TCRs). It is written in
C++ and will run in parallel across available CPU cores.

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
[`TCRdist_old()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_old.md),
[`TCRdist_to_igraph()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_to_igraph.md),
[`TCRdist_to_sparse_matrix()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_to_sparse_matrix.md),
[`add_to_tcr_network()`](https://nicholasclark.github.io/TIRTLtools/reference/add_to_tcr_network.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md)

## Examples

``` r
load_example_data(dataset = "SJTRC_minimal")
#> Example data already loaded into object: 'SJTRC_minimal'
df = get_all_tcrs(SJTRC_minimal, chain="paired", remove_duplicates = TRUE)

result = TCRdist_cpp(df)
#> Error in TCRdist_cpp(df): could not find function "TCRdist_cpp"

mat = result$matrix
#> Error: object 'result' not found
node_df = result$tcr1
#> Error: object 'result' not found

mat[1:5,1:5]
#> Error: object 'mat' not found
## note: TCRdist is symmetric. Returned matrix contains only lower-triangular values

node_df %>%
  mutate(alpha_nuc = paste(substr(alpha_nuc, 0, 20), "...", sep = ""),
         beta_nuc = paste(substr(beta_nuc, 0, 20), "...", sep = "")) %>%
  data.table::as.data.table()
#> Error in mutate(., alpha_nuc = paste(substr(alpha_nuc, 0, 20), "...",     sep = ""), beta_nuc = paste(substr(beta_nuc, 0, 20), "...",     sep = "")): could not find function "mutate"
```
