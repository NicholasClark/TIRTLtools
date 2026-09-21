# Build an undirected igraph graph from TCRdist() results

**\[experimental\]**

Builds an undirected `igraph` graph from the `edges_df`/`nodes_df`
produced by
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
(or
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
which uses the same names). Every edge gets weight 1 regardless of its
TCRdist value – i.e. this is a binary adjacency graph, not one weighted
by TCRdist.

## Usage

``` r
TCRdist_to_igraph(edges_df, nodes_df)
```

## Arguments

- edges_df:

  a data frame with columns "node1_idx" and "node2_idx" (1-indexed pairs
  of connected nodes), such as
  [`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)'s
  `$edges_df`. Any additional columns (e.g. "TCRdist") are kept as
  igraph edge attributes.

- nodes_df:

  a data frame with one row per node and a "tcr_index" column giving
  each node's 1-indexed id (matching `node1_idx`/`node2_idx`), such as
  [`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)'s
  `$nodes_df`. Its other columns become igraph vertex attributes. Nodes
  with no edges are included in the graph as isolated vertices.

## Value

an undirected `igraph` object with `vcount() == nrow(nodes_df)` and an
edge attribute `weight` equal to 1 for every edge.

## See also

[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`TCRdist_to_sparse_matrix()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_to_sparse_matrix.md)

Other tcr_similarity:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md),
[`TCRdist_old()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_old.md),
[`TCRdist_to_sparse_matrix()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_to_sparse_matrix.md),
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
#> ℹ 10% done — time taken so far: 0.25 seconds
#> ℹ 20% done — time taken so far: 0.5 seconds
#> ℹ 30% done — time taken so far: 0.74 seconds
#> ℹ 40% done — time taken so far: 0.98 seconds
#> ℹ 50% done — time taken so far: 1.21 seconds
#> ℹ 60% done — time taken so far: 1.31 seconds
#> ℹ 70% done — time taken so far: 1.87 seconds
#> ℹ 80% done — time taken so far: 1.97 seconds
#> ℹ 90% done — time taken so far: 2.22 seconds
#> ℹ 100% done — time taken so far: 2.26 seconds
#> ✔ Total time taken: 2.26 seconds
gr = TCRdist_to_igraph(result$edges_df, result$nodes_df)
```
