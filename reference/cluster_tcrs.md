# Cluster TCRs (using the Leiden algorithm) based on their pairwise TCRdist values

**\[experimental\]** The `cluster_tcrs()` function aggregates all of the
paired TCRs found in the data, calculates pairwise similarity using the
va, vb, cdr3a, and cdr3b regions (via TCRdist), and clusters the results
using the Leiden algorithm.

## Usage

``` r
cluster_tcrs(
  data,
  tcrdist_cutoff = 90,
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

  the
  [`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
  function will only record TCRdist values less than or equal to the
  cutoff. Default is 90. Note: Higher cutoffs will return more data, at
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

`$df` - a data frame with all unique TCRs along with cluster annotations

`$dist_df` - a data frame with distances (TCRdist) between TCR pairs in
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

Other repertoire_analysis:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md),
[`get_all_div_metrics()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_div_metrics.md),
[`get_all_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_tcrs.md),
[`get_pair_stats()`](https://nicholasclark.github.io/TIRTLtools/reference/get_pair_stats.md),
[`get_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/get_paired_by_read_fraction_range.md),
[`summarize_data()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_data.md)

## Examples

``` r
# example code
# paired = load_tirtlseq("your_directory/")
# obj = cluster_tcrs(paired)
```
