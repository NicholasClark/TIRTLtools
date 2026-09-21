# Build a sparse adjacency matrix from TCRdist() results

**\[experimental\]**

Builds a sparse adjacency matrix from the `edges_df`/`nodes_df` produced
by
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
(or
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md)
using
[`sparseMatrix()`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html).

## Usage

``` r
TCRdist_to_sparse_matrix(edges_df, nodes_df, binary = FALSE)
```

## Arguments

- edges_df:

  a data frame with columns "node1_idx" and "node2_idx", such as
  [`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)'s
  `$edges_df`.

- nodes_df:

  a data frame with one row per node, such as
  [`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)'s
  `$nodes_df`.

- binary:

  if TRUE, the nonzero values in the matrix will all be 1, otherwise
  they will be equal to TCRdist between the two TCRs (default if FALSE).

## Value

an `n x n` symmetric sparse matrix (class `dsCMatrix`), where
`n = nrow(nodes_df)`. Entry `[i, j]` is 1 if the two TCRs are connected
by an edge (TCRdist \<= cutoff) and 0 otherwise. Row/column names are
`nodes_df$tcr_index`.

## See also

[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`TCRdist_to_igraph()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_to_igraph.md)

Other tcr_similarity:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md),
[`TCRdist_old()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_old.md),
[`TCRdist_to_igraph()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_to_igraph.md),
[`add_to_tcr_network()`](https://nicholasclark.github.io/TIRTLtools/reference/add_to_tcr_network.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md)

## Examples

``` r
load_example_data(dataset = "SJTRC_minimal")
#> Example data already loaded into object: 'SJTRC_minimal'
df = get_all_tcrs(SJTRC_minimal, chain="paired", remove_duplicates = TRUE)
result = TCRdist(df, tcrdist_cutoff = 90)
#> Removed 384 TCRs with unknown V-segments (1.2%) from a total of 32,164 TCRs.
#> Removed 10 TCRs with short CDR3 segments (0.031%) from a total of 31,780 TCRs.
#> Removed 13,326 TCRs with non-functional CDR3 amino acid sequences (42%) from a total of 31,770 TCRs.
#> Filtered data frame contains 18,444 TCRs (57%) of original 32,164 TCRs.
#> ℹ Number of chunks: 190
#> ℹ 10% done — time taken so far: 0.26 seconds
#> ℹ 20% done — time taken so far: 0.36 seconds
#> ℹ 30% done — time taken so far: 0.45 seconds
#> ℹ 40% done — time taken so far: 0.7 seconds
#> ℹ 50% done — time taken so far: 0.79 seconds
#> ℹ 60% done — time taken so far: 1.04 seconds
#> ℹ 70% done — time taken so far: 1.14 seconds
#> ℹ 80% done — time taken so far: 1.39 seconds
#> ℹ 90% done — time taken so far: 1.49 seconds
#> ℹ 100% done — time taken so far: 1.53 seconds
#> ✔ Total time taken: 1.53 seconds
adj_mat = TCRdist_to_sparse_matrix(result$edges_df, result$nodes_df)
```
