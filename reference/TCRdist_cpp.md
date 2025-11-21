# Parallelized C++ implementation of TCRdist (no GPU required)

**\[experimental\]** This is an alternative to the GPU version of
TCRdist that is still very fast for large datasets (tens of thousands of
TCRs). It is written in C++ and will run in parallel across available
CPU cores.

## Usage

``` r
TCRdist_cpp(tcr1, tcr2 = NULL)
```

## Arguments

- tcr1:

  a data frame with one TCR per row. It must have the columns "va",
  "vb", "cdr3a", and "cdr3b"

- tcr2:

  (optional) another data frame of TCRs. If supplied, TCRdist will be
  calculated for every combination of one TCR from tcr1 and one TCR from
  tcr2. Otherwise, it will calculate TCRdist for every pair of TCRs in
  tcr1.

## Value

a list with two objects (or three if tcr2 is not null):

- matrix - a matrix of TCRdist values

- tcr1 - the input matrix tcr1, after pre-processing and removing
  unacceptable TCRs

- tcr2 (if supplied) - the input matrix tcr2, after pre-processing and
  removing unacceptable TCRs

## Details

This version of TCRdist is currently less feature-rich than
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
and returns a dense matrix as output. It does not yet allow for sparse
output or writing output directly to a file.

## See also

Other tcr_similarity:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md)
