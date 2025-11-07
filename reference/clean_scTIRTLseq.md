# Imputes missing alpha and beta chains where possible for single-cell TIRTLseq data

**\[experimental\]**

## Usage

``` r
clean_scTIRTLseq(df, verbose = TRUE, keep_all_columns = FALSE)
```

## Arguments

- df:

  a data frame with single-cell TIRTLseq data

- verbose:

  print a before and after summary of clonotypes (default TRUE)

- keep_all_columns:

  keep intermediate columns created by the function (default FALSE)

## See also

Other single-cell:
[`summarize_scTIRTLseq()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_scTIRTLseq.md)
