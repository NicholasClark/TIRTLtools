# Change the order of samples in a TIRTLseqData object

**\[experimental\]**

## Usage

``` r
reorder_samples(data, samples)
```

## Arguments

- data:

  a TIRTLseqData object

- samples:

  the samples of 'data', in their desired order. Either 1) a numeric
  vector of indices, 2) a character vector of sample names, or 3) a
  character vector of metadata conditions where each entry is of the
  form `"column==value"`.

## See also

Other data_wrangling:
[`add_metadata()`](https://nicholasclark.github.io/TIRTLtools/reference/add_metadata.md),
[`filter_dataset()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_dataset.md),
[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)
