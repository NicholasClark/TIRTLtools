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

Other tcr_similarity:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md),
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
#> 0.1 seconds
df = get_all_tcrs(sjtrc, chain="paired", remove_duplicates = TRUE)

result = cluster_tcrs(df)
#> Removed 1,470 MAIT TCRs (2.5%) from a total of 58,377 TCRs.
#> Removed 999 MAIT TCRs (3.1%) from a total of 32,164 TCRs.
#> Removed 384 TCRs with unknown V-segments (0.67%) from a total of 56,907 TCRs.
#> Removed 12 TCRs with short CDR3 segments (0.021%) from a total of 56,523 TCRs.
#> Removed 13,137 TCRs with non-functional CDR3 amino acid sequences (23%) from a total of 56,511 TCRs.
#> Filtered data frame contains 43,374 TCRs (76%) of original 56,907 TCRs.
#> Out of 43374 valid TCRs, 4914 clusters detected and 28631 singleton TCRs.
#> 94 clusters of size >= 10, 7 clusters of size >= 50, 3 clusters of size >=100.
```
