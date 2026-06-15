# Find TCRalpha/beta pairs from individual well read counts

**\[experimental\]** This runs the MAD-HYPE and T-SHELL algorithms to
find TCRalpha-beta pairs originating from the same clone.

## Usage

``` r
run_pairing(
  folder_path,
  folder_out,
  prefix,
  well_filter_thres = 0.5,
  min_reads = 0,
  min_wells = 2,
  well_pos = 3,
  tshell_settings = get_tshell_settings(format = "auto"),
  wellset = get_well_subset(1:16, 1:24),
  wellset1 = lifecycle::deprecated(),
  compute = TRUE,
  backend = c("auto", "cpu", "cupy", "mlx"),
  pval_thres_tshell = lifecycle::deprecated(),
  wij_thres_tshell = lifecycle::deprecated(),
  verbose = TRUE,
  write_extra_files = FALSE,
  filter_before_top3 = FALSE,
  fork = lifecycle::deprecated(),
  shared = lifecycle::deprecated(),
  chunk_size = 500,
  exclude_nonfunctional = FALSE,
  select_best_madhype = FALSE,
  select_best_tshell = FALSE
)
```

## Arguments

- folder_path:

  the path of the folder with well-level data

- folder_out:

  the path of the folder to write results to. The function will create
  the folder if it does not exist.

- prefix:

  a prefix for the output file names

- well_filter_thres:

  wells are removed if they have fewer unique clones than:
  wellfilter_thres\*(Avg. \# of unique clones per well). The default
  value is 0.5

- min_reads:

  minimum number of reads a chain must have in a well to be considered
  observed (note: actual minimum is min_reads+1. Default value is 0,
  i.e. chain must have \>= 1 read in a well)

- min_wells:

  minimum number of wells a chain must be observed in to be paired.

- well_pos:

  the position of the well ID (e.g. "B5") in the file names. For
  example, files named "\<well_id\>\_TCRalpha.tsv" would use well_pos=3.
  (default is 3)

- tshell_settings:

  a named list with pairing settings for T-SHELL: `pval_thres_tshell`
  and `wij_thres_tshell`. You may pass your own list or call
  [`get_tshell_settings()`](https://nicholasclark.github.io/TIRTLtools/reference/get_tshell_settings.md)
  with setting "auto", "384_well", or "96_well". By default, "auto" is
  selected and "384_well" settings will be used if there are \>= 150
  wells passing QC, otherwise "96_well" settings will be used.

  - `pval_thres_tshell` is the adjusted p-value threshold for T-SHELL
    significance (default 1e-10 for 384-well plate, 1e-3 for 96-well
    plate)

  - `wij_thres_tshell` is the threshold for the number of wells
    containing both chains for T-SHELL significance (default \>2 wells
    for 384-well plate, \>3 wells for 96-well plate)

- wellset:

  a vector of wells to use for the pairing

- compute:

  whether or not to run the pairing algorithms after tabulating and
  writing pseudobulk data (default TRUE)

- backend:

  the computing backend to use. The function looks for a GPU and
  automatically chooses an appropriate backend by default.

- verbose:

  whether to print out messages (default TRUE)

- write_extra_files:

  whether to write un-necessary intermediate files (default FALSE)

- filter_before_top3:

  whether to filter by loss fraction before extracting top 3 correlation
  values for T-SHELL (default FALSE)

- fork:

  whether to "fork" the python process for basilisk (default is NULL,
  which automatically chooses an appropriate option)

- shared:

  whether to use a "shared" python process for basilisk (default is
  NULL, which automatically chooses an appropriate option)

- chunk_size:

  batch size for calculations in pairing scripts

- exclude_nonfunctional:

  whether to exclude non-functional chains before pairing (default is
  FALSE)

- select_best_madhype:

  whether to use a secondary algorithm on the pairs from the MAD-HYPE
  algorithm to select the best pairs for each clone (default is FALSE)

- select_best_tshell:

  whether to use a secondary algorithm on the pairs from the T-SHELL
  algorithm to select the best pairs for each clone (default is FALSE)

## Value

A data frame with the TCR-alpha/TCR-beta pairs.

The function also writes three files to "folder_out":

- A data frame ("\_pseudobulk_TRA.tsv") of pseudobulk counts for
  TCRalpha chains

- A data frame ("\_pseudobulk_TRB.tsv") of pseudobulk counts for TCRbeta
  chains

- A data frame ("\_TIRTLoutput.tsv") of TCR-alpha/TCR-beta pairs.

These files can be loaded using the
[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)
function.

If write_extra_files is TRUE, the function also writes sparse matrices
of per-well read counts (well x chain) for TCR-alpha and beta to
"\_alpha_mat.rds" and "\_beta_mat.rds". Metadata for the chains in these
matrices are written to "\_alpha_meta.parquet" and "\_beta_meta.parquet"
and metadata for the wells is written to "\_well_meta.parquet".

These files can be loaded using the
[`load_well_counts_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/load_well_counts_binary.md)
function.

## See also

Other pairing:
[`get_tshell_settings()`](https://nicholasclark.github.io/TIRTLtools/reference/get_tshell_settings.md)
