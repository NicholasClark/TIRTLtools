# Load well-level data written to binary format

**\[experimental\]**

This function loads well-level data written to binary files (.parquet
and .rds) by the
[`write_well_data_to_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/write_well_data_to_binary.md)
function. The data contains sparse matrices (well x clone) of read
counts for TCRalpha and TCRbeta along with metadata data frames for each
clone.

Data loaded this way can be used as input to the
[`plot_tshell()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_tshell.md)
function.

## Usage

``` r
load_well_counts_binary(folder, prefix, lazy = TRUE)
```

## Arguments

- folder:

  the directory with ".parquet" and ".rds" files written by
  [`write_well_data_to_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/write_well_data_to_binary.md)

- prefix:

  a prefix with the sample name that will be prepended to the output
  file names

- lazy:

  whether to "lazy load" the clone metadata files. This can speed up
  loading because these files are usually very large (default is TRUE).

## Value

This function returns a list with the following elements:

- `alpha` - a sparse matrix ("dgCMatrix") of TCRalpha read counts in
  column-oriented format. See
  [`` Matrix::`dgCMatrix-class`() ``](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)

- `beta` - a sparse matrix ("dgCMatrix") of TCRbeta read counts in
  column-oriented format. See
  [`` Matrix::`dgCMatrix-class`() ``](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)

- `alpha_meta` - a data frame of clone metadata for TCRalpha clones. Row
  i relates to column i of `alpha`.

- `beta_meta` - a data frame of clone metadata for TCRbeta clones. Row i
  relates to column i of `beta`.

- `col_meta` - a data frame of well metadata (i.e. well names)

- `alpha_meta_file` - the file for the TCRalpha clone metadata

- `beta_meta_file` - the file for the TCRbeta clone metadata

## See also

Other well:
[`choose_pair_manual()`](https://nicholasclark.github.io/TIRTLtools/reference/choose_pair_manual.md),
[`get_well_subset()`](https://nicholasclark.github.io/TIRTLtools/reference/get_well_subset.md),
[`get_wells_from_edges()`](https://nicholasclark.github.io/TIRTLtools/reference/get_wells_from_edges.md),
[`plot_tshell()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_tshell.md),
[`write_well_data_to_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/write_well_data_to_binary.md)
