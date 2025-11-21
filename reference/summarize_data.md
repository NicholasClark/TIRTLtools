# Create a summary table with number of reads and unique alpha/beta chains

**\[experimental\]** This function creates a summary table with the
number of reads and unique alpha/beta chains observed for each sample.

## Usage

``` r
summarize_data(data)
```

## Arguments

- data:

  a TIRTLseqData object

## Value

a data frame with the number of reads and unique chains observed for
each sample.

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
[`plot_pairs_with_eachother()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_pairs_with_eachother.md),
[`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md),
[`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md),
[`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md)

## Examples

``` r
folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)
#> Loading files from: /Users/nclark52/git/TIRTLtools/inst/extdata/SJTRC_TIRTL_seq_longitudinal...
#> Found 6 beta chain pseudo-bulk files.
#> Found 6 paired chain files.
#> Loaded 18 files from 6 samples.
#> 12.4 seconds

summarize_data(ts_data)
#> # A tibble: 6 × 6
#>   sample_id  n_alpha_clones n_beta_clones n_clone_pairs n_reads_alpha
#>   <chr>               <int>         <int>         <int>         <int>
#> 1 cd4_tp1_v2        1972376       1867891         47828      76225265
#> 2 cd4_tp2_v2        1966962       1776782         32672      74060977
#> 3 cd4_tp3_v2        2062013       1930208         39012      66235760
#> 4 cd8_tp1_v2         937875        928777         14273      68573045
#> 5 cd8_tp2_v2         980522        895688         17891      68560963
#> 6 cd8_tp3_v2         850578        830552         13995      58978542
#> # ℹ 1 more variable: n_reads_beta <int>
```
