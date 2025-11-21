# Load data from TIRTLseq experiments

**\[experimental\]**

`load_tirtlseq()` loads paired-TCR and pseudo-bulk data from TIRTLseq
experiments from a given directory. It can also automatically assemble
metadata from filenames.

## Usage

``` r
load_tirtlseq(
  directory,
  chain = c("all", "paired", "alpha", "beta"),
  sep = "_",
  meta_columns = NULL,
  samples = NULL,
  pseudobulk_columns = "auto",
  paired_columns = "auto",
  n_threads = data.table::getDTthreads(),
  verbose = TRUE,
  stringsAsFactors = FALSE,
  n_max = Inf
)
```

## Arguments

- directory:

  the directory to look in for ".tsv" or ".tsv.gz" files of TIRTLseq
  data

- chain:

  which TCR chain data to load data for – "all" chains (alpha, beta, and
  paired) by default.

- sep:

  (optional) separator in the filename for metadata information ("\_" by
  default)

- meta_columns:

  (optional) a vector of names identifying the metadata contained in the
  filenames, for example `c("marker", "timepoint", "donor")` for files
  named something like "cd8_timepoint2_donor1 ... .tsv".

- samples:

  (optional) specific sample ids (the part of the filename before
  "\_pseudobulk" or "\_TIRTLoutput") to load. Default is NULL (loads all
  samples in the directory).

- pseudobulk_columns:

  (optional) the columns of the pseudobulk .tsv(.gz) to read. Either a
  list of columns or one of "auto", "all", or "minimal". "auto"
  (default) loads all columns except for some redundant ones. "all"
  loads all columns. "minimal" loads a small number of the most
  important columns.

- paired_columns:

  (optional) the columns of the paired .tsv(.gz) to read. Either a list
  of columns or one of "auto", "all", or "minimal". "auto" (default)
  loads all columns except for some redundant ones. "all" loads all
  columns. "minimal" loads a small number of the most important columns.

- n_threads:

  (optional) number of CPU threads to use for reading .tsv(.gz) files

- verbose:

  (optional) whether to print the name of each file loaded (default is
  TRUE).

- stringsAsFactors:

  (optional) read character strings in as factors

- n_max:

  (optional) the maximum number of files to read in – used mostly for
  testing purposes (default is Inf, i.e. read all files).

## Value

The function returns a list with the following structure, containing the
data from each sample and the sample metadata.

    your_data_object (list)
    ├───meta (metadata dataframe)
    └───data (list)
        └───sample_1 (list)
            ├───alpha (alpha pseudobulk dataframe)
            ├───beta (beta pseudobulk dataframe)
            └───paired (paired pseudobulk dataframe)
        ...
        └───sample_n (list)
            ├───alpha (alpha pseudobulk dataframe)
            ├───beta (beta pseudobulk dataframe)
            └───paired (paired pseudobulk dataframe)

## Details

The function expects ".tsv" (or ".tsv.gz") files. It looks for files
ending in "\_pseudobulk_TRA.tsv" (alpha-chain pseudo-bulk),
"\_pseudobulk_TRB.tsv" (beta-chain pseudo-bulk), and "\_TIRTLoutput.tsv"
(paired alpha and beta chains).

By default, the function will construct a metadata table with a row for
each sample, based on unique strings at the beginning of filenames
(before "\_TIRTLoutput.tsv" or similar). If the filename contains sample
metadata, then it can add multiple columns to the metadata table with
this information. For example, if a typical file looks like
"cd8_timepoint2_donor1_TIRTLoutput.tsv" and the user supplies
`c("cell_type", "timepoint", "donor")` for `meta_columns` and `"_"` for
`sep`, then the metadata table will look like something like this:

       sample_id             cell_type   timepoint       donor     label
         <chr>               <chr>         <chr>         <chr>     <chr>
    1 cd8_timepoint2_donor1    cd8       timepoint2      donor1    cell_type: cd8 | timepoint: timepoint2 | donor: donor1
    2 ...
    3 cd4_timepoint1_donor3    cd4       timepoint1      donor3    cell_type: cd4 | timepoint: timepoint1 | donor: donor3

## Examples

``` r
folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_")
#> Loading files from: /Users/nclark52/git/TIRTLtools/inst/extdata/SJTRC_TIRTL_seq_longitudinal...
#> Found 6 alpha chain pseudo-bulk files.
#> Found 6 beta chain pseudo-bulk files.
#> Found 6 paired chain files.
#> -- Loading files for sample: cd4_tp1_v2...
#> ---- Loading file -- Pseudobulk (alpha chain) -- cd4_tp1_v2_pseudobulk_TRA.tsv...
#> ---- Loading file -- Pseudobulk (beta chain)  -- cd4_tp1_v2_pseudobulk_TRB.tsv...
#> ---- Loading file -- TIRTLseq (paired chain)  -- cd4_tp1_v2_TIRTLoutput.tsv...
#> -- Loading files for sample: cd4_tp2_v2...
#> ---- Loading file -- Pseudobulk (alpha chain) -- cd4_tp2_v2_pseudobulk_TRA.tsv...
#> ---- Loading file -- Pseudobulk (beta chain)  -- cd4_tp2_v2_pseudobulk_TRB.tsv...
#> ---- Loading file -- TIRTLseq (paired chain)  -- cd4_tp2_v2_TIRTLoutput.tsv...
#> -- Loading files for sample: cd4_tp3_v2...
#> ---- Loading file -- Pseudobulk (alpha chain) -- cd4_tp3_v2_pseudobulk_TRA.tsv...
#> ---- Loading file -- Pseudobulk (beta chain)  -- cd4_tp3_v2_pseudobulk_TRB.tsv...
#> ---- Loading file -- TIRTLseq (paired chain)  -- cd4_tp3_v2_TIRTLoutput.tsv...
#> -- Loading files for sample: cd8_tp1_v2...
#> ---- Loading file -- Pseudobulk (alpha chain) -- cd8_tp1_v2_pseudobulk_TRA.tsv...
#> ---- Loading file -- Pseudobulk (beta chain)  -- cd8_tp1_v2_pseudobulk_TRB.tsv...
#> ---- Loading file -- TIRTLseq (paired chain)  -- cd8_tp1_v2_TIRTLoutput.tsv...
#> -- Loading files for sample: cd8_tp2_v2...
#> ---- Loading file -- Pseudobulk (alpha chain) -- cd8_tp2_v2_pseudobulk_TRA.tsv...
#> ---- Loading file -- Pseudobulk (beta chain)  -- cd8_tp2_v2_pseudobulk_TRB.tsv...
#> ---- Loading file -- TIRTLseq (paired chain)  -- cd8_tp2_v2_TIRTLoutput.tsv...
#> -- Loading files for sample: cd8_tp3_v2...
#> ---- Loading file -- Pseudobulk (alpha chain) -- cd8_tp3_v2_pseudobulk_TRA.tsv...
#> ---- Loading file -- Pseudobulk (beta chain)  -- cd8_tp3_v2_pseudobulk_TRB.tsv...
#> ---- Loading file -- TIRTLseq (paired chain)  -- cd8_tp3_v2_TIRTLoutput.tsv...
#> Loaded 18 files from 6 samples.
#> 12.3 seconds
```
