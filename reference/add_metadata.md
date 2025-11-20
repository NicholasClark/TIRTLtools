# Add metadata to a TIRTLseqData object

**\[experimental\]** This function adds a metadata data frame,
`obj$meta`, to the input TIRTLseqData object based on the sample names.
This can be useful when the metadata was not specified when the initial
object was loaded.

## Usage

``` r
add_metadata(obj, meta_columns = NULL, sep = "_")
```

## Arguments

- obj:

  a TIRTLseqData object

- meta_columns:

  (optional) a vector of names for metadata for each field in the sample
  names.

- sep:

  (optional) the character separating fields in the sample names of the
  data For example `c("marker", "timepoint", "donor")` for samples named
  something like "cd8_timepoint2_donor1".

## Value

a TIRTLseqData object with added metadata

## See also

Other data_wrangling:
[`filter_dataset()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_dataset.md),
[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md),
[`reorder_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/reorder_samples.md)
