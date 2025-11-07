# Plot the overlap/agreement between samples (in terms of most frequent clones)

**\[experimental\]** `plot_sample_overlap()` returns a heatmap showing
the overlap in most frequent TCRs among pairs of samples in a dataset

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

## Value

A heatmap with hierarchically clustered rows and columns showing the
number of TCRs shared between each pair of samples, among their top N
most frequent TCRs.

## Details

The function takes the top N most frequent TCRs found in each dataset
(default 200) and compares their overlap between samples.

## See also

[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)

Other plotting:
[`plot_clone_size_across_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clone_size_across_samples.md),
[`plot_clonotype_indices()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clonotype_indices.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md),
[`plot_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_diversity.md),
[`plot_n_reads()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_n_reads.md),
[`plot_num_partners()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_num_partners.md),
[`plot_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired.md),
[`plot_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_by_read_fraction_range.md),
[`plot_paired_vs_rank()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_vs_rank.md),
[`plot_pairs_with_eachother()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_pairs_with_eachother.md),
[`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md),
[`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md),
[`plot_sample_vs_sample()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_vs_sample.md)

## Examples

``` r
# example code
# data = load_tirtlseq("your_directory/")
# plot_sample_overlap(data, chain = "paired", n_seq=200)
```
