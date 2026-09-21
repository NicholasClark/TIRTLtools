# Add new TCRs to an existing TCR similarity network

**\[experimental\]**

Given a network previously built with
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
(or
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
this function calculates TCRdist between a new set of TCRs
(`nodes_new_df`) and all of the existing TCRs, then appends the new
nodes and edges to the network and returns the updated
`edges_df`/`nodes_df`.

## Usage

``` r
add_to_tcr_network(
  edges_df,
  nodes_df,
  nodes_new_df,
  remove_MAIT = FALSE,
  params = NULL,
  submat = NULL,
  tcrdist_cutoff = NULL,
  chunk_size = 1000,
  backend = c("auto", "cupy", "mlx", "cpp")
)
```

## Arguments

- edges_df:

  the `edges_df` from a previous
  [`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
  (or
  [`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md))
  result: a data frame with columns "node1_idx", "node2_idx", and
  "TCRdist".

- nodes_df:

  the `nodes_df` from a previous
  [`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
  (or
  [`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md))
  result: one row per existing TCR.

- nodes_new_df:

  a data frame of new TCRs to add to the network. Must have the columns
  "va", "vb", "cdr3a", and "cdr3b" (the same requirements as
  `tcr1`/`tcr2` in
  [`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)).

- remove_MAIT:

  whether to remove TCRs from MAIT cells (default is FALSE)

- params:

  (optional) a table of valid parameters for amino acids and va/vb
  segments. (default is NULL, which uses TIRTLtools::params)

- submat:

  (optional) a substitution matrix with mismatch penalties for each
  combination of amino acids or va/vb segments (default is NULL, which
  uses TIRTLtools::submat).

- tcrdist_cutoff:

  (optional) discard all TCRdist values above this cutoff. If not
  supplied by the user, this will default to 90 for dual-chain TCRdist
  or 45 for single-chain TCRdist.

- chunk_size:

  (optional) the chunk size to use in calculation of TCRdist (default
  1000). See
  [`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
  for details.

- backend:

  (optional) the backend to use for the chunk computation (default
  "auto"). See
  [`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
  for details.

## Value

a list with the same shape as
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)'s
output:

`$edges_df` - `edges_df` with the new TCRdist edges appended: new nodes
vs. existing nodes, and new nodes vs. each other.

`$nodes_df` - `nodes_df` with the new TCRs appended. New rows get
`source = "new"` and `tcr_index` values starting right after the highest
`tcr_index` already present in `nodes_df`.

## Details

`nodes_new_df` must contain the columns "va", "vb", "cdr3a", and
"cdr3b". These columns must contain the V-alpha segment, V-beta segment,
the CDR3-alpha amino acid sequence, and the CDR3-beta amino acid
sequence, respectively.

`nodes_df` must contain the columns "va", "vb", "cdr3a", "cdr3b", and
"tcr_index", where "tcr_index" is an integer and the other columns are
as above.

`edges_df` must contain the columns "node1_idx", "node2_idx", and
"TCRdist", where the first two columns are integers that map to the
"tcr_index" column in `nodes_df` and "TCRdist" is the value of TCRdist
between the two TCRs.

New nodes are numbered starting right after the highest `tcr_index`
already present in `nodes_df`.

## See also

[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md)

Other tcr_similarity:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md),
[`TCRdist_old()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_old.md),
[`TCRdist_to_igraph()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_to_igraph.md),
[`TCRdist_to_sparse_matrix()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_to_sparse_matrix.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md)

## Examples

``` r
load_example_data(dataset = "SJTRC_minimal")
#> Example data already loaded into object: 'SJTRC_minimal'
df = get_all_tcrs(SJTRC_minimal, chain = "paired", remove_duplicates = TRUE)
result = TCRdist(df[1:500, ], tcrdist_cutoff = 90)
#> Removed 10 TCRs with unknown V-segments (2%) from a total of 500 TCRs.
#> Removed 0 TCRs with short CDR3 segments (0%) from a total of 490 TCRs.
#> Removed 196 TCRs with non-functional CDR3 amino acid sequences (40%) from a total of 490 TCRs.
#> Filtered data frame contains 294 TCRs (59%) of original 500 TCRs.
#> ℹ Number of chunks: 1
#> ℹ 100% done — time taken so far: 0 seconds
#> ✔ Total time taken: 0 seconds
result2 = add_to_tcr_network(result$edges_df, result$nodes_df, df[501:600, ])
#> ℹ Both cdr3a and cdr3b found — using "tcrdist_cutoff = 90"
#> Removed 1 TCRs with unknown V-segments (1%) from a total of 100 TCRs.
#> Removed 0 TCRs with short CDR3 segments (0%) from a total of 99 TCRs.
#> Removed 42 TCRs with non-functional CDR3 amino acid sequences (42%) from a total of 99 TCRs.
#> Filtered data frame contains 57 TCRs (57%) of original 100 TCRs.
#> Removed 0 TCRs with unknown V-segments (0%) from a total of 294 TCRs.
#> Removed 0 TCRs with short CDR3 segments (0%) from a total of 294 TCRs.
#> Removed 0 TCRs with non-functional CDR3 amino acid sequences (0%) from a total of 294 TCRs.
#> Filtered data frame contains 294 TCRs (100%) of original 294 TCRs.
#> ℹ Number of chunks: 1
#> ℹ 100% done — time taken so far: 0 seconds
#> ✔ Total time taken: 0 seconds
#> ℹ Both cdr3a and cdr3b found — using "tcrdist_cutoff = 90"
#> Removed 1 TCRs with unknown V-segments (1%) from a total of 100 TCRs.
#> Removed 0 TCRs with short CDR3 segments (0%) from a total of 99 TCRs.
#> Removed 42 TCRs with non-functional CDR3 amino acid sequences (42%) from a total of 99 TCRs.
#> Filtered data frame contains 57 TCRs (57%) of original 100 TCRs.
#> ℹ Number of chunks: 1
#> ℹ 100% done — time taken so far: 0 seconds
#> ✔ Total time taken: 0 seconds

## compare to a single TCRdist() run on the combined set:
result_full = TCRdist(df[1:600, ], tcrdist_cutoff = 90)
#> Removed 11 TCRs with unknown V-segments (1.8%) from a total of 600 TCRs.
#> Removed 0 TCRs with short CDR3 segments (0%) from a total of 589 TCRs.
#> Removed 238 TCRs with non-functional CDR3 amino acid sequences (40%) from a total of 589 TCRs.
#> Filtered data frame contains 351 TCRs (58%) of original 600 TCRs.
#> ℹ Number of chunks: 1
#> ℹ 100% done — time taken so far: 0 seconds
#> ✔ Total time taken: 0 seconds
nrow(result2$edges_df) == nrow(result_full$edges_df)
#> [1] TRUE
```
