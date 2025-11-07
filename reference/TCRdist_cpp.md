# Parallelized C++ implementation of TCRdist (no GPU required)

**\[experimental\]**

## Usage

``` r
TCRdist_cpp(tcr1, tcr2 = NULL)
```

## Arguments

- tcr1:

  a data frame with one TCR per row. It must have the columns "va",
  "vb", "cdr3a", and "cdr3b"

- tcr2:

  (optional) another data frame of TCRs. If supplied, TCRdist will be
  calculated for every combination of one TCR from tcr1 and one TCR from
  tcr2. Otherwise, it will calculate TCRdist for every pair of TCRs in
  tcr1.

## See also

Other repertoire_analysis:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`get_all_div_metrics()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_div_metrics.md),
[`get_all_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_tcrs.md),
[`get_pair_stats()`](https://nicholasclark.github.io/TIRTLtools/reference/get_pair_stats.md),
[`get_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/get_paired_by_read_fraction_range.md),
[`summarize_data()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_data.md)
