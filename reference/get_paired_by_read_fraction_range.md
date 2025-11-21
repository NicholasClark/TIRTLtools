# Calculate the number and fraction of single chains that were paired by frequency

**\[experimental\]** This function returns the number and fraction of
single-chains (alpha or beta, default is beta) that were paired within
different frequency ranges (default is
`[10^-1, 10^-2], [10^-2, 10^-3], ... , [10^-5,10^-6], [10^-6, 0]`) for
each sample. You can use this function to get the data used in the plot
for
[`plot_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_by_read_fraction_range.md).

## Usage

``` r
get_paired_by_read_fraction_range(
  data,
  chain = c("beta", "alpha"),
  cutoffs = 10^(-6:-1)
)
```

## Arguments

- data:

  a TIRTLseqData object

- chain:

  the TCR chain to use for read fraction (default is "beta")

- cutoffs:

  a vector of cutoffs for the read fraction ranges

## Value

a data frame with the number and fraction of chains paired for each read
fraction range for each sample.

## See also

Other qc:
[`get_all_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_tcrs.md),
[`get_pair_stats()`](https://nicholasclark.github.io/TIRTLtools/reference/get_pair_stats.md),
[`plot_n_reads()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_n_reads.md),
[`plot_num_partners()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_num_partners.md),
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
folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)
#> Loading files from: /Users/nclark52/git/TIRTLtools/inst/extdata/SJTRC_TIRTL_seq_longitudinal...
#> Found 6 beta chain pseudo-bulk files.
#> Found 6 paired chain files.
#> Loaded 18 files from 6 samples.
#> 12.1 seconds

get_paired_by_read_fraction_range(ts_data, chain = "beta")
#> # A tibble: 33 × 13
#>    range          n_paired_sum n_paired_sum_tshell n_paired_sum_madhype n_total
#>    <fct>                 <int>               <int>                <int>   <int>
#>  1 [0,1e-06)              6325                2114                 5361 1584466
#>  2 [1e-06,1e-05)         16357               10807                15716  276468
#>  3 [1e-05,0.0001)         6038                4774                 6014    6725
#>  4 [0.0001,0.001)          224                 191                  221     231
#>  5 [0.001,0.01)              1                   1                    0       1
#>  6 [0,1e-06)              4258                1789                 3392 1482037
#>  7 [1e-06,1e-05)         11577                8386                10988  288341
#>  8 [1e-05,0.0001)         4846                3529                 4819    6227
#>  9 [0.0001,0.001)          167                 131                  155     176
#> 10 [0.001,0.01)              1                   1                    0       1
#> # ℹ 23 more rows
#> # ℹ 8 more variables: fraction_paired <dbl>, fraction_paired_tshell <dbl>,
#> #   fraction_paired_madhype <dbl>, sample_id <chr>, marker <chr>,
#> #   timepoint <chr>, version <chr>, label <chr>

```
