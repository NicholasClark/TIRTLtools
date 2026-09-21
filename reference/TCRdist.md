# A fast implementation of TCRdist, a distance/similarity metric for TCRs

**\[experimental\]**

An efficient GPU-enabled version of TCRdist with an almost-as-fast CPU
version backup. If a GPU is available it will run a version of TCRdist
using the `cupy` (NVIDIA) or `mlx` (Apple Silicon) Python package via
`reticulate`. If no GPU is available, it will run a CPU-only C++
version. The C++ version is slower than the GPU, but still relatively
fast.

## Usage

``` r
TCRdist(
  tcr1,
  tcr2 = NULL,
  remove_MAIT = FALSE,
  params = NULL,
  submat = NULL,
  tcrdist_cutoff = NULL,
  chunk_size = 1000,
  write_to_tsv = FALSE,
  output_folder = ".",
  backend = c("auto", "cupy", "mlx", "cpp")
)
```

## Arguments

- tcr1:

  a data frame with one TCR per row. It must have the columns "va",
  "vb", "cdr3a", and "cdr3b". These columns must contain the V-alpha
  segment, V-beta segment, the CDR3-alpha amino acid sequence, and the
  CDR3-beta amino acid sequence, respectively.

- tcr2:

  (optional) another data frame of TCRs. If supplied, TCRdist will be
  calculated for every combination of one TCR from tcr1 and one TCR from
  tcr2. Otherwise, it will calculate TCRdist for each pair of TCRs in
  tcr1.

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

  (optional) The chunk size to use in calculation of TCRdist (default
  1000). If set at n, it will calculate pairwise TCRdist for n x n TCRs
  at once.

- write_to_tsv:

  (optional) write the results to a tab-separated file ".tsv" (default
  is FALSE, does not write .tsv file)

- output_folder:

  (optional) folder to write output ".tsv" files to, if `write_to_tsv`
  is TRUE (default is the current directory).

- backend:

  (optional) the backend to use for the chunk computation (default
  "auto"). One of "auto" (a GPU backend if available, otherwise C++),
  "cpp" (a parallel C++ implementation via RcppParallel – fast,
  CPU-only, and does not require python at all), "cupy" (NVIDIA GPU, via
  python), or "mlx" (Apple Silicon GPU, via python).

## Value

If write_to_tsv is TRUE (default is FALSE), the function will write
output .tsv files and return NULL. Otherwise, it will return a list with
entries:

`$edges_df` - a data frame with three columns: "node1_idx", "node2_idx",
and "TCRdist". The first two columns contain the indices of the TCRs for
each pair, matching `nodes_df$tcr_index`. The last column contains the
TCRdist if it is below the cutoff. The output is sparse in that it only
contains pairs that have TCRdist \<= cutoff.

`$nodes_df` - a data frame of the TCRs supplied to the function. It
contains an additional column "tcr_index" with the index of each TCR. If
`tcr2` was supplied, this is `bind_rows(tcr1, tcr2)`: tcr1's rows are
numbered first (starting at 1), and tcr2's `tcr_index` values continue
on immediately after the highest `tcr_index` in tcr1. Note that any TCRs
with invalid V-segments or frameshifts/stop-codons in their amino acid
sequence will be removed.

## Details

This function calculates pairwise TCRdist (Dash et al., Nature 2017) for
a set of TCRs (or between two sets of TCRs) and returns a sparse output
with the TCRdist and indices of all pairs that have TCRdist less than or
equal to a desired cutoff (default cutoff is 90).

## See also

[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md),
and
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md)

Other tcr_similarity:
[`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md),
[`TCRdist_old()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_old.md),
[`TCRdist_to_igraph()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_to_igraph.md),
[`TCRdist_to_sparse_matrix()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_to_sparse_matrix.md),
[`add_to_tcr_network()`](https://nicholasclark.github.io/TIRTLtools/reference/add_to_tcr_network.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md)

## Examples

``` r
load_example_data(dataset = "SJTRC_minimal")
#> Loading file: SJTRC_minimal.qs2...
#> 8.971 sec elapsed
df = get_all_tcrs(SJTRC_minimal, chain="paired", remove_duplicates = TRUE)
result = TCRdist(df, tcrdist_cutoff = 90)
#> Removed 384 TCRs with unknown V-segments (1.2%) from a total of 32,164 TCRs.
#> Removed 10 TCRs with short CDR3 segments (0.031%) from a total of 31,780 TCRs.
#> Removed 13,326 TCRs with non-functional CDR3 amino acid sequences (42%) from a total of 31,770 TCRs.
#> Filtered data frame contains 18,444 TCRs (57%) of original 32,164 TCRs.
#> ℹ Number of chunks: 190
#> ℹ 10% done — time taken so far: 0.75 seconds
#> ℹ 20% done — time taken so far: 0.99 seconds
#> ℹ 30% done — time taken so far: 1.09 seconds
#> ℹ 40% done — time taken so far: 1.34 seconds
#> ℹ 50% done — time taken so far: 1.59 seconds
#> ℹ 60% done — time taken so far: 1.68 seconds
#> ℹ 70% done — time taken so far: 1.92 seconds
#> ℹ 80% done — time taken so far: 2.16 seconds
#> ℹ 90% done — time taken so far: 2.4 seconds
#> ℹ 100% done — time taken so far: 2.44 seconds
#> ✔ Total time taken: 2.45 seconds
edge_df = result[['edges_df']] ### table of TCRdist values <= cutoff
node_df = result[['nodes_df']] ### table of input metadata with indices
```
