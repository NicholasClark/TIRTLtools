# Process single-cell TIRTL-seq data

**\[experimental\]** This function reads a directory of MiXCR output
files from single-cell TIRTL-seq and outputs a data frame with a single
TCR for each well.

## Usage

``` r
process_scTIRTLseq(
  folder,
  min_read_fraction = 0.1,
  min_read_count = 15,
  require_functional = TRUE
)
```

## Arguments

- folder:

  the path of the folder with .tsv (or .tsv.gz) output files from MiXCR

- min_read_fraction:

  the minimum read fraction required to select a chain for a well

- min_read_count:

  the minimum read count required to select a chain for a well

- require_functional:

  whether to require a functional chain (i.e. without frameshift or stop
  codons in the CDR3 region)

## Value

A data frame with one TCR pair selected per well.

## See also

Other single-cell:
[`plot_scTIRTLseq_plate()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_scTIRTLseq_plate.md)
