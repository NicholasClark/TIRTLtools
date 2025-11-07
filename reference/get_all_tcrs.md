# Returns all of the paired TCRs from all samples in a dataset

**\[experimental\]** The `get_all_tcrs()` function aggregates the TCRs
from all samples of a dataset and puts them into one table.

## Usage

``` r
get_all_tcrs(
  data,
  chain = c("paired", "alpha", "beta"),
  remove_duplicates = TRUE
)
```

## Arguments

- data:

  a TIRTLseqData object created by
  [`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)

- chain:

  the TCR chain, "alpha", "beta", or "paired" (default is paired)

- remove_duplicates:

  only return one TCR for TCRs paired by both the T-SHELL and MAD-HYPE
  algorithms (default is TRUE).

## Value

A dataframe including all of the TCRs in a dataset.

## Details

A pair of TCRs is included twice in the TIRTLseq data if it is
recognized by both the T-SHELL and MAD-HYPE algorithms. If
remove_duplicates is TRUE (default) the function will only return one of
these pairs of TCRs.

## See also

[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)

Other repertoire_analysis:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`get_all_div_metrics()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_div_metrics.md),
[`get_pair_stats()`](https://nicholasclark.github.io/TIRTLtools/reference/get_pair_stats.md),
[`get_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/get_paired_by_read_fraction_range.md),
[`summarize_data()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_data.md)

## Examples

``` r
# example code

```
