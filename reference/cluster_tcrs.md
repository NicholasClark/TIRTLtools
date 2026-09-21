# Cluster TCRs (using the Leiden algorithm) based on their pairwise TCRdist values

**\[experimental\]**

The `cluster_tcrs()` function aggregates all of the paired TCRs found in
the data, calculates pairwise similarity using the va, vb, cdr3a, and
cdr3b regions (via TCRdist), and clusters the results using the Leiden
algorithm.

## Usage

``` r
cluster_tcrs(
  data,
  tcrdist_cutoff = NULL,
  resolution = 0.1,
  with_db = TRUE,
  db = TIRTLtools::vdj_db,
  allow_self_edges = TRUE,
  remove_MAIT = TRUE
)
```

## Arguments

- data:

  a list of TIRTLseq TCR data for samples created with
  [`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)

- tcrdist_cutoff:

  discard all TCRdist values above this cutoff. If not supplied by the
  user, this will default to 90 for dual-chain TCRdist or 45 for
  single-chain TCRdist. Note: Higher cutoffs will return more data, at
  most NxN where N is the number of unique TCRs.

- resolution:

  the "resolution" parameter for the Leiden algorithm. A lower value
  will produce larger clusters and a higher value will produce smaller
  clusters. Typical values are in the 0.1 - 2.0 range. A higher value
  may be better for densely connected data while a lower value may be
  better for moderately connected data. Default is 0.1.

- with_db:

  if TRUE, observed clones will be compared and clustered with a
  dataframe of annotated clones. By default, a dataframe with VDJ-db
  annotations is used.

- db:

  a data frame with annotated TCRs. The default is the VDJ-db database.

- allow_self_edges:

  if FALSE, only calculate TCRdist between input data TCRs and the TCR
  annotation database (db). If TRUE, calculate pairwise TCRdist for all
  of the data including the input and the annotated TCRs.

- remove_MAIT:

  remove MAIT TCRs before clustering (default is TRUE)

## Value

Returns a list with the following elements:

`$nodes_df` - a data frame with all unique TCRs along with cluster
annotations

`$edges_df` - a data frame with distances (TCRdist) between TCR pairs in
long format

`$sparse_adj_mat` - an adjacency matrix (in sparse format) marking TCR
pairs with TCRdist \<= tcrdist_cutoff

`$graph_adj` - an igraph object created from the adjacency matrix

`$tcrdist_cutoff` - the cutoff used for TCRdist

`$resolution` - the resolution parameter used for the Leiden algorithm

## Details

The function also filters the dataset to TCRs that are valid for
TCRdist.

The following TCRs are removed:

- TCRs that contain stop codons (\*) or frame shifts (\_) in their cdr3a
  or cdr3b regions

- TCRs that contain a cdr3 region with 5 or less amino acids

- TCRs that contain a v segment allele not found in our parameter table

V-segments that do not specify an allele (e.g. "TRAV1-2" instead of
"TRAV1-2\*01") will be assigned to the "\*01" allele.

## See also

[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)

Other tcr_similarity:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md),
[`TCRdist_old()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_old.md),
[`TCRdist_to_igraph()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_to_igraph.md),
[`TCRdist_to_sparse_matrix()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_to_sparse_matrix.md),
[`add_to_tcr_network()`](https://nicholasclark.github.io/TIRTLtools/reference/add_to_tcr_network.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md)

## Examples

``` r
load_example_data(dataset = "SJTRC_minimal")
#> Example data already loaded into object: 'SJTRC_minimal'
df = get_all_tcrs(SJTRC_minimal, chain="paired", remove_duplicates = TRUE)

result = cluster_tcrs(df)
#> ℹ Both cdr3a and cdr3b found — using "tcrdist_cutoff = 90"
#> Removed 1,583 MAIT TCRs (2.2%) from a total of 71,206 TCRs.
#> Removed 999 MAIT TCRs (3.1%) from a total of 32,164 TCRs.
#> Removed 452 TCRs with unknown V-segments (0.65%) from a total of 69,623 TCRs.
#> Removed 13 TCRs with short CDR3 segments (0.019%) from a total of 69,171 TCRs.
#> Removed 13,139 TCRs with non-functional CDR3 amino acid sequences (19%) from a total of 69,158 TCRs.
#> Filtered data frame contains 56,019 TCRs (80%) of original 69,623 TCRs.
#> ℹ Number of chunks: 1653
#> ℹ 10% done — time taken so far: 1.38 seconds
#> ℹ 20% done — time taken so far: 4.52 seconds
#> ℹ 30% done — time taken so far: 5.32 seconds
#> ℹ 40% done — time taken so far: 6.66 seconds
#> ℹ 50% done — time taken so far: 7.99 seconds
#> ℹ 60% done — time taken so far: 9.34 seconds
#> ℹ 70% done — time taken so far: 10.68 seconds
#> ℹ 80% done — time taken so far: 11.99 seconds
#> ℹ 90% done — time taken so far: 13.31 seconds
#> ℹ 100% done — time taken so far: 13.86 seconds
#> ✔ Total time taken: 13.86 seconds
#> Out of 56019 valid TCRs, 5724 clusters detected and 37998 singleton TCRs.
#> 128 clusters of size >= 10, 13 clusters of size >= 50, 5 clusters of size >=100.
```
