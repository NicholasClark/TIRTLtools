# Line plot of clone read fraction across multiple samples

**\[experimental\]**

## Usage

``` r
plot_clone_size_across_samples(
  data,
  clones,
  chain = c("beta", "alpha"),
  pseudo = 1e-06,
  group_vec = NULL,
  sum_readFraction = TRUE,
  samples = NULL,
  return_data = FALSE,
  label_zero = FALSE,
  show_legend = TRUE,
  log_scale = TRUE,
  x_var = NULL
)
```

## Arguments

- data:

  a TIRTLseqData object

- clones:

  a list of nucleotide sequences of TCR clones

- chain:

  the TCR chain used (default is "beta")

- pseudo:

  the value of the pseudocount added to all read fractions (default is
  10^-6)

- group_vec:

  (optional) a vector of "groups" for the clones

- sum_readFraction:

  if TRUE, plot the sum of read fractions of clones in each group. If
  FALSE, plot each clone with a separate line.

- samples:

  (optional) which samples to use in the plot (default is all samples)

- return_data:

  whether to return the data used for plotting (default is FALSE)

- label_zero:

  whether to label zero on the y-axis (default is FALSE)

- show_legend:

  whether to show the legend (default is TRUE)

- log_scale:

  (optional) if TRUE, use log-scale for the y-axis (default is FALSE)

## See also

Other plotting:
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
[`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md),
[`plot_sample_vs_sample()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_vs_sample.md)
