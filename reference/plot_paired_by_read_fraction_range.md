# Bar plot of the fraction of paired single chains by frequency

**\[experimental\]** This function returns a bar plot of the fraction of
single-chains (alpha or beta, default is beta) that were paired within
different frequency ranges (default is
`[10^-1, 10^-2], [10^-2, 10^-3], ... , [10^-5,10^-6], [10^-6, 0]`) for
each sample.

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

## Value

a ggplot object with the bar plot of fraction of paired chains by
frequency

## See also

Other qc:
[`get_all_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_tcrs.md),
[`get_pair_stats()`](https://nicholasclark.github.io/TIRTLtools/reference/get_pair_stats.md),
[`get_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/get_paired_by_read_fraction_range.md),
[`plot_n_reads()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_n_reads.md),
[`plot_num_partners()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_num_partners.md),
[`plot_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired.md),
[`plot_paired_vs_rank()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_vs_rank.md),
[`plot_pairs_with_eachother()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_pairs_with_eachother.md),
[`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md),
[`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md),
[`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md),
[`summarize_data()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_data.md)

## Examples

``` r
folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)
#> Loading files from: /Users/nclark52/git/TIRTLtools/inst/extdata/SJTRC_TIRTL_seq_longitudinal...
#> Found 6 beta chain pseudo-bulk files.
#> Found 6 paired chain files.
#> Loaded 18 files from 6 samples.
#> 12.7 seconds

plot_paired_by_read_fraction_range(ts_data, chain = "beta")

```
