# Table of permissible amino acids and V-segments for TCRdist

Table of permissible amino acids and V-segments for TCRdist

## Usage

``` r
params
```

## Format

A data frame with two columns. The data frame contains 266 features,
mapping to 161 row/column indices.

The first column is called "feature" and contains TCR V-segments (e.g.
"TRAV12-1\*02") and one-letter abbreviations for amino acids. The second
column is called "value" and contains the row/column index in the
substitution matrix corresponding to that amino acid or V-segment. Note
that the first row/column is 0 (instead of 1) because the TCRdist code
is in Python, which is 0-indexed.

## Examples

``` r
params
#> # A tibble: 266 × 2
#>    feature value
#>    <chr>   <dbl>
#>  1 _           0
#>  2 A           1
#>  3 C           2
#>  4 D           3
#>  5 E           4
#>  6 F           5
#>  7 G           6
#>  8 H           7
#>  9 I           8
#> 10 K           9
#> # ℹ 256 more rows
```
