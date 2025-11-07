# A stepped plot of the cumulative number of paired (or unpaired) single-chains for the N most frequent single-chains

**\[experimental\]**

## Usage

``` r
plot_paired_vs_rank(
  data,
  sample = 1,
  y_axis = c("n_not_paired", "n_paired"),
  chain = c("both", "beta", "alpha"),
  n_max = 100,
  color_scheme = NULL
)
```

## Arguments

- data:

  a TIRTLseqData object

- sample:

  the sample to plot

- y_axis:

  whether to plot the cumulative number of paired or un-paired
  single-chains (default is "n_not_paired")

- chain:

  the TCR chain to plot (default is both alpha and beta)

- n_max:

  the number of most frequent single-chains to plot

- color_scheme:

  (optional) the color scheme for the plot

## See also

Other plotting:
[`plot_clone_size_across_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clone_size_across_samples.md),
[`plot_clonotype_indices()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clonotype_indices.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md),
[`plot_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_diversity.md),
[`plot_n_reads()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_n_reads.md),
[`plot_num_partners()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_num_partners.md),
[`plot_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired.md),
[`plot_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_by_read_fraction_range.md),
[`plot_pairs_with_eachother()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_pairs_with_eachother.md),
[`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md),
[`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md),
[`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md),
[`plot_sample_vs_sample()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_vs_sample.md)
