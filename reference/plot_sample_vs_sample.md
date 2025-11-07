# Scatterplot of TCR clone read fraction of clones between two samples

**\[experimental\]** `plot_sample_vs_sample()` returns a scatterplot of
read frequencies of TCRs between two samples

## Usage

``` r
plot_sample_vs_sample(
  data1,
  data2,
  chain = c("beta", "alpha"),
  log2fc_cutoff = 3,
  sem_cutoff = 2.5,
  smooth_sem = c("window", "none"),
  window_size = 30,
  end_window_size = 5,
  pseudo1 = 1e-06,
  pseudo2 = 1e-06,
  labelx = "Frequency on timepoint 1",
  labely = "Frequency on timepoint 2",
  return_data = FALSE
)
```

## Arguments

- data1:

  a list of three data frames (alpha, beta, and paired) for one sample

- data2:

  a list of three data frames (alpha, beta, and paired) for one sample

- log2_cutoff:

  the log2 fold-change cutoff to call a TCR up- or down-regulated
  (default 1.5)

## Value

A scatterplot (ggplot object) with read frequencies (proportions),
colored by whether each TCR is up-regulated, down-regulated, or neither,
given the log2 fold-change cutoff.

## Details

The function labels each TCR as up-regulated, down-regulated, or stable,
based on the log2 fold-change cutoff supplied (default 1.5).

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
[`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md),
[`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md)

## Examples

``` r
# example code

```
