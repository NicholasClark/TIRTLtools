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

## Examples

``` r
folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal",
  package = "TIRTLtools")
sjtrc = load_tirtlseq(folder,
  meta_columns = c("marker", "timepoint", "version"), sep = "_",
  chain = "paired", verbose = FALSE)
#> Loading files from: /Users/nclark52/git/TIRTLtools/inst/extdata/SJTRC_TIRTL_seq_longitudinal...
#> Found 6 beta chain pseudo-bulk files.
#> Found 6 paired chain files.
#> Loaded 6 files from 6 samples.
#> 0.5 seconds
print(sjtrc$meta)
#> # A tibble: 6 × 5
#>   sample_id  marker timepoint version label                                     
#>   <chr>      <chr>  <chr>     <chr>   <chr>                                     
#> 1 cd4_tp1_v2 cd4    tp1       v2      marker: cd4 | timepoint: tp1 | version: v2
#> 2 cd4_tp2_v2 cd4    tp2       v2      marker: cd4 | timepoint: tp2 | version: v2
#> 3 cd4_tp3_v2 cd4    tp3       v2      marker: cd4 | timepoint: tp3 | version: v2
#> 4 cd8_tp1_v2 cd8    tp1       v2      marker: cd8 | timepoint: tp1 | version: v2
#> 5 cd8_tp2_v2 cd8    tp2       v2      marker: cd8 | timepoint: tp2 | version: v2
#> 6 cd8_tp3_v2 cd8    tp3       v2      marker: cd8 | timepoint: tp3 | version: v2
new_order = names(sjtrc$data) %>% rev()
print(new_order)
#> [1] "cd8_tp3_v2" "cd8_tp2_v2" "cd8_tp1_v2" "cd4_tp3_v2" "cd4_tp2_v2"
#> [6] "cd4_tp1_v2"
sjtrc = reorder_samples(sjtrc, new_order)
print(sjtrc$meta)
#> # A tibble: 6 × 5
#>   sample_id  marker timepoint version label                                     
#>   <chr>      <chr>  <chr>     <chr>   <chr>                                     
#> 1 cd8_tp3_v2 cd8    tp3       v2      marker: cd8 | timepoint: tp3 | version: v2
#> 2 cd8_tp2_v2 cd8    tp2       v2      marker: cd8 | timepoint: tp2 | version: v2
#> 3 cd8_tp1_v2 cd8    tp1       v2      marker: cd8 | timepoint: tp1 | version: v2
#> 4 cd4_tp3_v2 cd4    tp3       v2      marker: cd4 | timepoint: tp3 | version: v2
#> 5 cd4_tp2_v2 cd4    tp2       v2      marker: cd4 | timepoint: tp2 | version: v2
#> 6 cd4_tp1_v2 cd4    tp1       v2      marker: cd4 | timepoint: tp1 | version: v2
```
