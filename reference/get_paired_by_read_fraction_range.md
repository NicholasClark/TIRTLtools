# Calculate the number and fraction of single chains that were paired within different read fraction ranges

**\[experimental\]**

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

## See also

Other repertoire_analysis:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`get_all_div_metrics()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_div_metrics.md),
[`get_all_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_tcrs.md),
[`get_pair_stats()`](https://nicholasclark.github.io/TIRTLtools/reference/get_pair_stats.md),
[`summarize_data()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_data.md)
