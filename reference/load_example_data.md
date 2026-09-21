# Load example TIRTL-seq data

**\[experimental\]**

This is a helper function to load example TIRTL-seq data. The function
will also download the data if needed from
https://github.com/NicholasClark/TIRTLtools/releases/tag/data-v1 and
save it in the package's R cache directory.

The function currently loads one of two available TIRTL-seq datasets,
either "SJTRC_minimal" or "SJTRC_longitudinal". The data is loaded to
the global environment and assigned with the same name as the dataset.
For example, if `dataset` is "SJTRC_minimal", the data will be set to
`SJTRC_minimal` in the global environment. If the dataset is already
loaded, it will not be re-loaded.

## Usage

``` r
load_example_data(
  dataset = c("SJTRC_minimal", "SJTRC_longitudinal"),
  verbose = TRUE
)
```

## Arguments

- dataset:

  either "SJTRC_minimal" or "SJTRC_longitudinal".

- verbose:

  whether to print messages to the user (default TRUE)

## Value

The function returns NULL invisibly, but loads a dataset and assigns it
to either `SJTRC_minimal` or `SJTRC_longitudinal` in the global
environment.

## Details

The `SJTRC_longitudinal` dataset contains 6 samples of TIRTL-seq data
from the St. Jude Tracking Study of Immune Responses Associated with
COVID-19 (SJTRC). The data is from one donor at three timepoints,
isolated for either CD4+ or CD8+ T-cells.

The `SJTRC_minimal` dataset is a subset of the above dataset, containing
only two samples: CD8+ isolated T-cells for timepoints 1 and 2.

## See also

Other data_import:
[`get_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/get_samples.md),
[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md),
[`read_external_bulk()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_bulk.md),
[`read_external_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_paired.md)

## Examples

``` r
load_example_data("SJTRC_minimal")
#> Example data already loaded into object: 'SJTRC_minimal'
print(SJTRC_minimal)
#> <TIRTLseqDataSet>
#>   Number of samples: 2 
#>   Samples: cd8_tp1_v2 cd8_tp2_v2
summary(SJTRC_minimal)
#> <TIRTLseqDataSet>
#>   Number of samples: 2 
#>   Per-sample TCR counts (alpha / beta / paired):
#> # A tibble: 2 × 5
#>   sample_id  `Paired TCRs` `Unique Paired TCRs` `Alpha Chains` `Beta Chains`
#>   <chr>              <int>                <int>          <int>         <int>
#> 1 cd8_tp1_v2         20690                14273         937875        928777
#> 2 cd8_tp2_v2         26260                17891         980522        895688
load_example_data("SJTRC_longitudinal")
#> Example data already loaded into object: 'SJTRC_longitudinal'
print(SJTRC_longitudinal)
#> <TIRTLseqDataSet>
#>   Number of samples: 6 
#>   Samples: cd4_tp1_v2 cd4_tp2_v2 cd4_tp3_v2 cd8_tp1_v2 cd8_tp2_v2 cd8_tp3_v2
summary(SJTRC_longitudinal)
#> <TIRTLseqDataSet>
#>   Number of samples: 6 
#>   Per-sample TCR counts (alpha / beta / paired):
#> # A tibble: 6 × 5
#>   sample_id  `Paired TCRs` `Unique Paired TCRs` `Alpha Chains` `Beta Chains`
#>   <chr>              <int>                <int>          <int>         <int>
#> 1 cd4_tp1_v2         67831                47828        1972376       1867891
#> 2 cd4_tp2_v2         48365                32672        1966962       1776782
#> 3 cd4_tp3_v2         57632                39012        2062013       1930208
#> 4 cd8_tp1_v2         20690                14273         937875        928777
#> 5 cd8_tp2_v2         26260                17891         980522        895688
#> 6 cd8_tp3_v2         20550                13995         850578        830552
```
