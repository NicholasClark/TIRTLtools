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

[`diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/calculate_diversity.md),
[`plot_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_diversity.md)

Other diversity:
[`calculate_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/calculate_diversity.md),
[`plot_clonotype_indices()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clonotype_indices.md),
[`plot_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_diversity.md)

## Examples

``` r
# example code
get_all_div_metrics()
#>  [1] "simpson"        "gini"           "gini.simpson"   "inv.simpson"   
#>  [5] "shannon"        "berger.parker"  "richness"       "d50"           
#>  [9] "dXX"            "renyi"          "hill"           "top10fraction" 
#> [13] "top100fraction" "topNfraction"  
```
