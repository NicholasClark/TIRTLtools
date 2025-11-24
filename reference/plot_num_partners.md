# Stacked bar plot of the fraction of alpha/beta chains with different numbers of partners

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
  plot The default is NULL, which uses the sample id.

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

A bar chart (ggplot object) with facets (sub-plots) for each sample. If
`return_data` is TRUE, the data frame used to create the plot is
returned instead.

## Details

For each sample, the function creates stacked bar charts for alpha and
beta chains, showing the proportion of them (among all called pairs)
that are paired with 1 chain, 2 chains, 3 chains, etc.

## See also

[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md)

Other qc:
[`get_all_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_tcrs.md),
[`get_pair_stats()`](https://nicholasclark.github.io/TIRTLtools/reference/get_pair_stats.md),
[`get_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/get_paired_by_read_fraction_range.md),
[`plot_n_reads()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_n_reads.md),
[`plot_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired.md),
[`plot_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_by_read_fraction_range.md),
[`plot_paired_vs_rank()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_vs_rank.md),
[`plot_pairs_with_eachother()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_pairs_with_eachother.md),
[`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md),
[`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md),
[`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md),
[`summarize_data()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_data.md)

## Examples

``` r
folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal",
  package = "TIRTLtools")
ts_data = load_tirtlseq(folder,
  meta_columns = c("marker","timepoint", "version"),
  sep = "_", verbose = FALSE)
#> Loading files from: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/TIRTLtools/extdata/SJTRC_TIRTL_seq_longitudinal...
#> Found 6 beta chain pseudo-bulk files.
#> Found 6 paired chain files.
#> Loaded 18 files from 6 samples.
#> 14.6 seconds

plot_num_partners(ts_data)


```
