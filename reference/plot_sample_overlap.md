# Plot the overlap/agreement between samples (in terms of most frequent clones)

**\[experimental\]** This function returns a heatmap showing the overlap
among the most frequent `n_seq` TCRs (default is n_seq=200) among pairs
of samples in a dataset.

## Usage

``` r
plot_sample_overlap(
  data,
  chain = c("paired", "alpha", "beta"),
  n_seq = 200,
  show_row_names = FALSE,
  show_column_names = FALSE,
  label_col = "Sample",
  title = "",
  return_data = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)

sample_overlap(...)
```

## Arguments

- data:

  the dataset, an object loaded using the
  [`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)
  function

- chain:

  which chain to plot: either paired or alpha-/beta-pseudobulk. (default
  "paired")

- n_seq:

  the number of most frequent TCR sequences to use (default 200)

- show_row_names:

  whether to show row names for the heatmap (default FALSE)

- show_column_names:

  whether to show column names for the heatmap (default FALSE)

- label_col:

  a column of the metadata to use as labels for rows and columns
  (default "Sample", uses the sample_id)

- title:

  a title for the heatmap

- return_data:

  whether to return the data used for plotting (a matrix with the
  overlap values) instead of a heatmap (default is FALSE)

- cluster_rows:

  whether to cluster the rows of the heatmap (default TRUE)

- cluster_cols:

  whether to cluster the columns of the heatmap (default TRUE)

## Value

A heatmap with hierarchically clustered rows and columns showing the
number of TCRs shared between each pair of samples, among their top N
most frequent TCRs.

If `return_data` is TRUE, a matrix of overlap values will be returned
instead.

## See also

[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)

Other qc:
[`get_all_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_tcrs.md),
[`get_pair_stats()`](https://nicholasclark.github.io/TIRTLtools/reference/get_pair_stats.md),
[`get_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/get_paired_by_read_fraction_range.md),
[`plot_n_reads()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_n_reads.md),
[`plot_num_partners()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_num_partners.md),
[`plot_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired.md),
[`plot_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_by_read_fraction_range.md),
[`plot_paired_vs_rank()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_vs_rank.md),
[`plot_pairs_with_eachother()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_pairs_with_eachother.md),
[`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md),
[`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md),
[`summarize_data()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_data.md)

## Examples

``` r
folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)
#> Loading files from: /Users/nclark52/git/TIRTLtools/inst/extdata/SJTRC_TIRTL_seq_longitudinal...
#> Found 6 beta chain pseudo-bulk files.
#> Found 6 paired chain files.
#> Loaded 18 files from 6 samples.
#> 12.6 seconds

plot_sample_overlap(ts_data, chain = "beta")

plot_sample_overlap(ts_data, chain = "alpha")

plot_sample_overlap(ts_data, chain = "paired")

```
