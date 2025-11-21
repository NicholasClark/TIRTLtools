# Substitution penalty matrix for TCRdist amino acids and V-segments

Substitution penalty matrix for TCRdist amino acids and V-segments

## Usage

``` r
submat
```

## Format

A symmetric matrix (161 x 161) with substitution penalties for pairs of
amino acids or V-segments for TCRdist calculation. The amino
acids/V-segments corresponding to each row/column are listed in the
"params" table.

## Examples

``` r
dim(submat)
#> [1] 161 161
submat[1:5,1:5]
#>      V1 V2 V3 V4 V5
#> [1,]  0 12 12 12 12
#> [2,] 12  0 12 12 12
#> [3,] 12 12  0 12 12
#> [4,] 12 12 12  0  6
#> [5,] 12 12 12  6  0
## to find the amino acid or v-segment(s) that map to a specific row/column
index = 106 ## 106th row/column of matrix
params[params$value == index-1,] ## need to subtract 1 because values in params are 0-indexed
#> # A tibble: 6 × 2
#>   feature     value
#>   <chr>       <dbl>
#> 1 TRBV20-1*01   105
#> 2 TRBV20-1*02   105
#> 3 TRBV20-1*04   105
#> 4 TRBV20-1*05   105
#> 5 TRBV20-1*06   105
#> 6 TRBV20-1*07   105
```
