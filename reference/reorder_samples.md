# Re-order samples in a TIRTLseqData object

**\[experimental\]** This function changes the order of the samples in a
TIRTLseqData object to the order that the user specifies.

## Usage

``` r
reorder_samples(data, samples)
```

## Arguments

- data:

  a TIRTLseqData object

- samples:

  the samples of 'data', in their desired order.

  Either:

  1.  A numeric vector of indices

  2.  A character vector of sample names

  3.  A character vector of metadata conditions where each entry is of
      the form `"column==value"`

## Value

A TIRTLseqData object with re-ordered samples.

## See also

Other data_wrangling:
[`add_metadata()`](https://nicholasclark.github.io/TIRTLtools/reference/add_metadata.md),
[`filter_dataset()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_dataset.md)
