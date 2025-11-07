# Returns all diversity metric options for [`calculate_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/calculate_diversity.md)

**\[experimental\]**

## Usage

``` r
get_all_div_metrics()
```

## Value

A vector of all available diversity metrics for the
[`calculate_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/calculate_diversity.md)
function.

## See also

`diversity()`,
[`plot_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_diversity.md)

Other repertoire_analysis:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`get_all_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_tcrs.md),
[`get_pair_stats()`](https://nicholasclark.github.io/TIRTLtools/reference/get_pair_stats.md),
[`get_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/get_paired_by_read_fraction_range.md),
[`summarize_data()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_data.md)

## Examples

``` r
# example code
get_all_div_metrics()
#>  [1] "simpson"        "gini"           "gini.simpson"   "inv.simpson"   
#>  [5] "shannon"        "berger.parker"  "richness"       "d50"           
#>  [9] "dXX"            "renyi"          "hill"           "top10fraction" 
#> [13] "top100fraction" "topNfraction"  
```
