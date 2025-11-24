# Subset a TIRTLseqData object

**\[experimental\]** The `filter_dataset()` function is used to select a
subset of samples from a loaded TIRTLseq dataset and create a new
dataset object.

## Usage

``` r
filter_dataset(data, samples)
```

## Arguments

- data:

  a TIRTLseqData object created by
  [`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)

- samples:

  the selected samples. Either 1) a numeric vector of indices, 2) a
  character vector of sample names, or 3) a character vector of metadata
  conditions where each entry is of the form `"column==value"`.

## Value

A dataset object similar to that created by
[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md),
but with only the selected samples.

## Details

The function accepts 1) a numeric vector of indices, 2) a character
vector of sample names, or 3) a character vector of metadata conditions
where each entry is of the form `"column==value"`.

In the third case, `c("cell_type==cd4", "timepoint==tp2")` would, for
example, select samples whose `cell_type` is `cd4` and whose `timepoint`
is `tp2` in the sample metadata.

## See also

Other data_wrangling:
[`add_metadata()`](https://nicholasclark.github.io/TIRTLtools/reference/add_metadata.md),
[`reorder_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/reorder_samples.md)

## Examples

``` r
folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal",
  package = "TIRTLtools")
sjtrc = load_tirtlseq(folder,
  meta_columns = c("marker", "timepoint", "version"), sep = "_",
  chain = "paired", verbose = FALSE)
#> Loading files from: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/TIRTLtools/extdata/SJTRC_TIRTL_seq_longitudinal...
#> Found 6 beta chain pseudo-bulk files.
#> Found 6 paired chain files.
#> Loaded 6 files from 6 samples.
#> 0.5 seconds
p2 = filter_dataset(sjtrc, 1:3) ### by indices
p3 = filter_dataset(sjtrc, c("cd8_tp1_v2", "cd8_tp2_v2", "cd8_tp3_v2")) ### by sample names
p4 = filter_dataset(sjtrc, "marker==cd4") ### by sample metadata condition
p5 = filter_dataset(sjtrc, c("marker==cd4", "timepoint==tp2")) ### by multiple sample metadata conditions
```
