# Read and process bulk single-chain TCR-seq data

**\[experimental\]**

This function reads and processes bulk single-chain TCR-sequencing data
from a non-TIRTLseq assay. Currently, only MiXCR format is supported.

## Usage

``` r
read_external_bulk(path, format = "MiXCR")
```

## Arguments

- path:

  the path to the data file

- format:

  the data format. Currently only "MiXCR" is supported.

## Value

A list containing the following slots:

- df - a data frame with a few columns modified and renamed

- df_raw - the original data, un-modified

## Details

Supported data types:

- `"MiXCR"` - "\<sample_name\>\_TRA.tsv" and "\<sample_name\>\_TRB.tsv"

## See also

Other data_import:
[`get_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/get_samples.md),
[`load_example_data()`](https://nicholasclark.github.io/TIRTLtools/reference/load_example_data.md),
[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md),
[`read_external_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_paired.md)
