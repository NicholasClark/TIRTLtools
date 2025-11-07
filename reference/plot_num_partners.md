# Stacked bar plot of the fraction of single-chains with different numbers of partners for each sample

**\[experimental\]** `plot_num_partners()` creates bar plots for alpha
and beta chains showing how many partners they are paired with by the
MAD-HYPE and/or T-shell algorithms.

## Usage

``` r
plot_num_partners(
  data,
  group_col = NULL,
  fraction = TRUE,
  include_non_functional = FALSE,
  max_partners = 5,
  samples = NULL,
  return_data = FALSE,
  color_scheme = NULL
)
```

## Arguments

- data:

  a TIRTLseqData object

- group_col:

  a column of the metadata to use to group multiple samples into one bar
  plot

- fraction:

  whether to plot the fraction of chains or the total number of chains
  (default is TRUE, i.e. plot fractions)

- include_non_functional:

  whether to include chains with non-functional cdr3 sequences when
  tabulating the output.

- max_partners:

  the maximum number of partners, N, to include in the plots. All chains
  with more than N partners will be grouped together under the "\>N"
  category.

- samples:

  (optional) which samples to plot

- return_data:

  if TRUE, return the data used to make the plots

- color_scheme:

  the color scheme to use for the plot

## Value

Either a bar chart (ggplot object) with facets (sub-plots) for each
sample or a list with two objects:

`$plot` the plot referenced above

`$data` the data used to create the plot

## Details

For each sample, the function creates stacked bar charts for alpha and
beta chains, showing the proportion of them (among all called pairs)
that are paired with 1 chain, 2 chains, 3 chains, etc.

## See also

[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md)

Other plotting:
[`plot_clone_size_across_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clone_size_across_samples.md),
[`plot_clonotype_indices()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clonotype_indices.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md),
[`plot_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_diversity.md),
[`plot_n_reads()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_n_reads.md),
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
