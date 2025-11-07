# Alternative point plot of read fraction vs. rank for the N most frequent alpha and beta chains, with a cross indicating an unpaired single-chain

**\[experimental\]**

## Usage

``` r
plot_read_fraction_vs_pair_status(
  data,
  sample = 1,
  chain = c("both", "beta", "alpha"),
  n_max = 100,
  show_num_partners = F,
  color_scheme = NULL
)
```

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
[`plot_paired_vs_rank()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_vs_rank.md),
[`plot_pairs_with_eachother()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_pairs_with_eachother.md),
[`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md),
[`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md),
[`plot_sample_vs_sample()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_vs_sample.md)
