# Bar plot of the fraction of single chains that were paired within different read fraction ranges for each sample.

**\[experimental\]**

## Usage

``` r
plot_paired_by_read_fraction_range(
  data,
  chain = c("beta", "alpha"),
  cutoffs = 10^(-6:-1),
  freq = FALSE,
  samples = NULL,
  color_scheme = NULL
)
```

## Arguments

- data:

  a TIRTLseqData object or a data frame created using
  [`get_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/get_paired_by_read_fraction_range.md)

- chain:

  the TCR chain to plot (default is "beta")

- cutoffs:

  a vector of cutoffs for the read fraction ranges

- freq:

  if TRUE, plot the number of pairs, if FALSE plot the fraction paired
  (default is FALSE, plot the fraction paired)

- samples:

  (optional) the samples to include in the plot

- color_scheme:

  (optional) the color scheme to use in the plot

## See also

Other plotting:
[`plot_clone_size_across_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clone_size_across_samples.md),
[`plot_clonotype_indices()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clonotype_indices.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md),
[`plot_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_diversity.md),
[`plot_n_reads()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_n_reads.md),
[`plot_num_partners()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_num_partners.md),
[`plot_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired.md),
[`plot_paired_vs_rank()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_vs_rank.md),
[`plot_pairs_with_eachother()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_pairs_with_eachother.md),
[`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md),
[`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md),
[`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md),
[`plot_sample_vs_sample()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_vs_sample.md)
