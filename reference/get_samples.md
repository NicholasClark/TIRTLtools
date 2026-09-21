# Returns sample names that are available to be imported from a directory

The `get_samples()` takes a directory path and searches the directory
for TIRTL-seq files – files ending in `_TIRTLoutput.tsv(.gz)`,
`_pseudobulk_TRA.tsv(.gz)`, or `_pseudobulk_TRB.tsv(.gz)`. You can give
a regex pattern to return only samples that match this pattern.

## Usage

``` r
get_samples(data_dirs, pattern = NULL)
```

## Arguments

- data_dirs:

  a directory path (or a vector with multiple paths) with TIRTL-seq data
  files

- pattern:

  a filename pattern

## Value

A dataframe with a column for the sample name and other columns to
indicate if the sample has paired TCR data, TCRalpha pseudobulk data,
and/or TCRbeta pseudobulk data.

## See also

Other data_import:
[`load_example_data()`](https://nicholasclark.github.io/TIRTLtools/reference/load_example_data.md),
[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md),
[`read_external_bulk()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_bulk.md),
[`read_external_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_paired.md)
