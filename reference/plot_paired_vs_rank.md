# A step plot of the cumulative number of paired/unpaired alpha/beta chains among the most frequent chains

**\[experimental\]** The function creates a stepped line plot of the
cumulative number of paired/unpaired alpha or beta chains (default is
both) among the `n_max` most frequent chains.

## Usage

``` r
plot_paired_vs_rank(
  data,
  sample,
  y_axis = c("n_not_paired", "n_paired"),
  chain = c("both", "beta", "alpha"),
  n_max = 100,
  color_scheme = NULL
)
```

## Arguments

- data:

  a TIRTLseqDataSet object

- sample:

  the sample to plot (either by number or sample id)

- y_axis:

  whether to plot the cumulative number of paired or un-paired
  single-chains (default is "n_not_paired")

- chain:

  the TCR chain to plot (default is both alpha and beta)

- n_max:

  the number of most frequent single-chains to plot

- color_scheme:

  (optional) the color scheme for the plot

## Value

A ggplot object with the plot

## See also

Other qc:
[`get_all_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_tcrs.md),
[`get_pair_stats()`](https://nicholasclark.github.io/TIRTLtools/reference/get_pair_stats.md),
[`get_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/get_paired_by_read_fraction_range.md),
[`plot_n_reads()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_n_reads.md),
[`plot_num_partners()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_num_partners.md),
[`plot_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired.md),
[`plot_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_by_read_fraction_range.md),
[`plot_pairs_with_eachother()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_pairs_with_eachother.md),
[`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md),
[`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md),
[`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md),
[`summarize_data()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_data.md)

## Examples

``` r
folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)
#> Loading files from: /Users/nclark52/git/temp_build/TIRTLtools/inst/extdata/SJTRC_TIRTL_seq_longitudinal...
#> Found 6 beta chain pseudo-bulk files.
#> Found 6 paired chain files.
#> Loaded 18 files from 6 samples.
#> 15.5 seconds

plot_paired_vs_rank(ts_data, sample = 1, y_axis = "n_not_paired", chain = "both", n_max = 100)


plot_paired_vs_rank(ts_data, sample = 1, y_axis = "n_paired", chain = "both", n_max = 100)


```
