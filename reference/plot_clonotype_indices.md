# Stacked bar chart with fractions of reads attributed to the most frequent clonotypes

**\[experimental\]** `plot_clonotype_indices()` creates a stacked bar
chart containing the fraction of reads for the top 10 most frequent
clonotypes, the 11th to 100th most frequent clonotypes, the 101st to
1000th most frequent, and so on, for each sample in the dataset.

## Usage

``` r
plot_clonotype_indices(
  data,
  chain = c("beta", "alpha"),
  cutoffs = 10^(1:5),
  group_col = NULL,
  label_col = "Sample",
  flip = FALSE,
  return_data = FALSE,
  color_scheme = NULL
)
```

## Arguments

- data:

  a TIRTLseqData object

- chain:

  the TCR chain to use (default is "beta")

- cutoffs:

  the indices used for the end of each group in the bar chart. The
  default is 10^(1:5), i.e. the 10th most-frequent clone, the 100th
  most-frequent clone, the 1,000th, the 10,000th, and the 100,000th.

- group_col:

  (optional) if supplied, a column of the metadata that will be used to
  group samples

- label_col:

  (optional) labels for the samples

- flip:

  (optional) if TRUE, flip the x and y-axes (default is FALSE)

- return_data:

  (optional) if TRUE, return the data used to make the plot (default is
  FALSE)

- color_scheme:

  (optional) a color scheme to use

## Value

(default) Returns a stacked bar chart of relative frequencies of
most-frequent clonotypes. If return_data is TRUE, the data frame used to
create the plot is returned instead.

## See also

`diversity()`,
[`plot_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_diversity.md)

Other plotting:
[`plot_clone_size_across_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clone_size_across_samples.md),
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
[`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md),
[`plot_sample_vs_sample()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_vs_sample.md)

## Examples

``` r
# example code

```
