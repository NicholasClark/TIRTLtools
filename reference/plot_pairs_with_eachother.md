# A connected point plot of read fraction vs. rank for the most frequent alpha/beta chains

**\[experimental\]** This function creates a point plot of read fraction
vs. rank for the N most frequent alpha (left, mirrored) and beta (right)
chains with lines between alpha and beta chains indicating a pair and a
cross indicating an unpaired single-chain.

## Usage

``` r
plot_pairs_with_eachother(
  data,
  sample,
  n_max = 100,
  show_num_partners = FALSE,
  color_scheme = NULL
)
```

## Arguments

- data:

  a TIRTLseqData object

- sample:

  the sample to plot (either by number or sample id)

- n_max:

  the number of most frequent single-chains to plot

- show_num_partners:

  whether to show the number of partners for each single-chain (default
  is FALSE)

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
[`plot_paired_vs_rank()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_vs_rank.md),
[`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md),
[`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md),
[`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md),
[`summarize_data()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_data.md)

## Examples

``` r
folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)
#> Loading files from: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/TIRTLtools/extdata/SJTRC_TIRTL_seq_longitudinal...
#> Found 6 beta chain pseudo-bulk files.
#> Found 6 paired chain files.
#> Loaded 18 files from 6 samples.
#> 15 seconds

plot_pairs_with_eachother(ts_data, sample = 1, n_max = 100)


plot_pairs_with_eachother(ts_data, sample = 1, n_max = 100, show_num_partners = TRUE)


```
