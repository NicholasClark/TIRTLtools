# GPU implementation of TCRdist, a distance/similarity metric for pairs of TCRs

**\[experimental\]** An efficient, batched version of TCRdist that is
compatible with both NVIDIA and Apple Silicon GPUs.

## Usage

``` r
TCRdist(
  tcr1,
  tcr2 = NULL,
  remove_MAIT = FALSE,
  params = NULL,
  submat = NULL,
  tcrdist_cutoff = 90,
  chunk_size = 1000,
  print_chunk_size = 10,
  print_res = TRUE,
  only_lower_tri = TRUE,
  return_data = TRUE,
  write_to_tsv = FALSE,
  backend = c("auto", "cpu", "cupy", "mlx"),
  fork = NULL,
  shared = NULL
)
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

  (optional) discard all TCRdist values above this cutoff (default is
  90).

- chunk_size:

  (optional) The chunk size to use in calculation of TCRdist (default
  1000). If set at n, it will calculate pairwise TCRdist for n x n TCRs
  at once. This may be as high as allowable by GPU memory (in our
  testing, a chunk_size of 1000 to 5000 provided the fastest runtime and
  chunk_size of over 7500 resulted in memory errors on some GPUs).

- print_chunk_size:

  (optional) print a line of output for every n TCRs processed (default
  1000)

- print_res:

  (optional) print summary of results (default is TRUE)

- only_lower_tri:

  (optional) return one TCRdist value for each pair (like the lower
  triangle of a symmetric matrix). Default is TRUE.

- return_data:

  (optional) whether to return the output result from the function. With
  large data it may be desirable to write the result to disk instead.
  (default is TRUE, returns output)

- write_to_tsv:

  (optional) write the results to a tab-separated file ".tsv" (default
  is FALSE, does not write .tsv file)

- backend:

  (optional) the CPU or GPU backend to use (default "auto")

- fork:

  (optional) a TRUE/FALSE value for whether to "fork" a new Python
  process for running TCRdist via the "basilisk" package. Default is
  NULL, which should use choose a safe value based on how the package is
  loaded.

- shared:

  (optional) a TRUE/FALSE value for whether to "share" the Python
  process for running TCRdist via the "basilisk" package. Default is
  NULL, which should use choose a safe value based on how the package is
  loaded.

## Value

A list with entries:

`$TCRdist_df` - a data frame with three columns: "node1_0index",
"node2_0index", and "TCRdist". The first two columns contain the indices
(0-indexed) of the TCRs for each pair. The last column contains the
TCRdist if it is below the cutoff. The output is sparse in that it only
contains pairs that have TCRdist \<= cutoff.

`$tcr1` - a data frame of the TCRs supplied to the function. It contains
an additional column "tcr_index" with the (0-indexed) index of each TCR.

`$tcr2` - a similar data frame for tcr2, if it was supplied.

## Details

This function calculates pairwise TCRdist (Dash et al., Nature 2017) for
a set of TCRs (or between two sets of TCRs) and returns a sparse output
with the TCRdist and indices of all pairs that have TCRdist less than or
equal to a desired cutoff (default cutoff is 90).

The function uses the `reticulate` package to call a python script that
uses `cupy` (NVIDIA GPUs), `mlx` (Apple Silicon GPUs), or `numpy` (no
GPU) to calculate TCRdist efficiently.

## See also

[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md),
and
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md)

Other tcr_similarity:
[`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md)

## Examples

``` r
folder = system.file("extdata/SJTRC_TIRTLseq_minimal",
  package = "TIRTLtools")
sjtrc = load_tirtlseq(folder,
  meta_columns = c("marker", "timepoint", "version"), sep = "_",
  chain = "paired", verbose = FALSE)
#> Loading files from: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/TIRTLtools/extdata/SJTRC_TIRTLseq_minimal...
#> Found 2 beta chain pseudo-bulk files.
#> Found 2 paired chain files.
#> Loaded 2 files from 2 samples.
#> 0.2 seconds
df = get_all_tcrs(sjtrc, chain="paired", remove_duplicates = TRUE)
result = TCRdist(df, tcrdist_cutoff = 90)
#> Removed 384 TCRs with unknown V-segments (1.2%) from a total of 32,164 TCRs.
#> Removed 12 TCRs with short CDR3 segments (0.038%) from a total of 31,780 TCRs.
#> Removed 13,324 TCRs with non-functional CDR3 amino acid sequences (42%) from a total of 31,768 TCRs.
#> Filtered data frame contains 18,444 TCRs (57%) of original 32,164 TCRs.
edge_df = result[['TCRdist_df']] %>%
  data.table::as.data.table() ### table of TCRdist values <= cutoff
node_df = result[['tcr1']] %>%
  data.table::as.data.table() ### table of input data with indices

edge_df ## sparse 3-column output: node1, node2, TCRdist
#>        node1_0index node2_0index TCRdist
#>               <int>        <int>   <int>
#>     1:          158          157       0
#>     2:          286          284      84
#>     3:          287          283      89
#>     4:          288          284      89
#>     5:          289          285      89
#>    ---                                  
#> 18422:        18434        18073       0
#> 18423:        18437        18406      72
#> 18424:        18440        18241      90
#> 18425:        18440        18324      66
#> 18426:        18440        18422      84
## note that indices start at 0 and are found in node_df$tcr_index

node_df %>%
  select(tcr_index, everything()) %>%
  mutate(alpha_nuc = paste(substr(alpha_nuc, 0, 20), "...", sep = ""),
         beta_nuc = paste(substr(beta_nuc, 0, 20), "...", sep = ""))
#> Error in mutate(., alpha_nuc = paste(substr(alpha_nuc, 0, 20), "...",     sep = ""), beta_nuc = paste(substr(beta_nuc, 0, 20), "...",     sep = "")): could not find function "mutate"

```
