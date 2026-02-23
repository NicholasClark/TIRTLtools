# Count the number of pairs called by each algorithm

**\[experimental\]** This function tabulates the number of TCR pairs
called for each sample by the MAD-HYPE algorithm, the T-SHELL algorithm,
or both. It is used in the
[`plot_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired.md)
function.

## Usage

``` r
get_pair_stats(data, verbose = TRUE)
```

## Arguments

- data:

  a TIRTLseqDataSet object

- verbose:

  (optional) whether to print progress of the function (default is
  TRUE).

## Value

A data frame with the number of pairs called by each algorithm for each
sample.

## Details

The function can count the number of alpha chains paired, or beta
chains, or the number of pairs overall. By default, it counts the number
of pairs overall.

## See also

Other qc:
[`get_all_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_tcrs.md),
[`get_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/get_paired_by_read_fraction_range.md),
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
#> Loading files from: /Users/nclark52/git/temp_build/TIRTLtools/inst/extdata/SJTRC_TIRTL_seq_longitudinal...
#> Found 6 beta chain pseudo-bulk files.
#> Found 6 paired chain files.
#> Loaded 18 files from 6 samples.
#> 18.2 seconds

get_pair_stats(ts_data)
#> Calculating pairing stats for sample 1
#> Calculating pairing stats for sample 2
#> Calculating pairing stats for sample 3
#> Calculating pairing stats for sample 4
#> Calculating pairing stats for sample 5
#> Calculating pairing stats for sample 6
#> # A tibble: 66 × 8
#>    category         Freq chain  sample_id  marker timepoint version label       
#>    <fct>           <int> <chr>  <chr>      <chr>  <chr>     <chr>   <chr>       
#>  1 MAD-HYPE only   25327 paired cd4_tp1_v2 cd4    tp1       v2      marker: cd4…
#>  2 T-SHELL only     2498 paired cd4_tp1_v2 cd4    tp1       v2      marker: cd4…
#>  3 both            20003 paired cd4_tp1_v2 cd4    tp1       v2      marker: cd4…
#>  4 MAD-HYPE only   14343 alpha  cd4_tp1_v2 cd4    tp1       v2      marker: cd4…
#>  5 T-SHELL only     1991 alpha  cd4_tp1_v2 cd4    tp1       v2      marker: cd4…
#>  6 both            19137 alpha  cd4_tp1_v2 cd4    tp1       v2      marker: cd4…
#>  7 neither       1936905 alpha  cd4_tp1_v2 cd4    tp1       v2      marker: cd4…
#>  8 MAD-HYPE only   11058 beta   cd4_tp1_v2 cd4    tp1       v2      marker: cd4…
#>  9 T-SHELL only     1633 beta   cd4_tp1_v2 cd4    tp1       v2      marker: cd4…
#> 10 both            16254 beta   cd4_tp1_v2 cd4    tp1       v2      marker: cd4…
#> # ℹ 56 more rows
```
