# Add metadata to a TIRTLseqData object

**\[experimental\]**

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

## See also

Other data_wrangling:
[`filter_dataset()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_dataset.md),
[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md),
[`reorder_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/reorder_samples.md)
