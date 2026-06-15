# Plot a plate of single-cell TIRTL-seq data

**\[experimental\]** This function plots which wells in a plate are
missing an alpha and/or beta chain for single-cell TIRTL-seq data.

## Usage

``` r
plot_scTIRTLseq_plate(df, title = "")
```

## Arguments

- df:

  a data frame, output from the
  [`process_scTIRTLseq()`](https://nicholasclark.github.io/TIRTLtools/reference/process_scTIRTLseq.md)
  function

- title:

  an optional title for the plot

## Value

A ggplot object.

## See also

Other single-cell:
[`process_scTIRTLseq()`](https://nicholasclark.github.io/TIRTLtools/reference/process_scTIRTLseq.md)
