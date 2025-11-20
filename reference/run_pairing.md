# Find TCRalpha/beta pairs from individual well read counts

**\[experimental\]** This runs the MAD-HYPE and T-SHELL algorithms to
find TCRalpha-beta pairs originating from the same clone.

## Usage

``` r
run_pairing(
  folder_path,
  folder_out,
  prefix = "tmp",
  well_filter_thres = 0.75,
  min_reads = 0,
  min_wells = 2,
  well_pos = 3,
  wellset1 = get_well_subset(1:16, 1:24),
  compute = T,
  backend = c("auto", "cpu", "cupy", "mlx"),
  pval_thres_tshell = 1e-10,
  wij_thres_tshell = 2,
  verbose = TRUE,
  write_extra_files = FALSE,
  filter_before_top3 = FALSE,
  fork = NULL,
  shared = NULL
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

- wellset1:

  a vector of wells to use for the pairing

- compute:

  whether or not to run the pairing algorithms after tabulating and
  writing pseudobulk data (default TRUE)

- backend:

  the computing backend to use. The function looks for a GPU and
  automatically chooses an appropriate backend by default.

- pval_thres_tshell:

  the adjusted p-value threshold for T-SHELL significance (default
  1e-10)

- wij_thres_tshell:

  the threshold for the number of wells containing both chains for
  T-SHELL significance (default \>2 wells)

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

- wellfilter_thres:

  wells are removed if they have fewer unique clones than:
  wellfilter_thres\*(Avg. \# of unique clones per well). The default
  value is 0.75

## Value

A data frame with the TCR-alpha/TCR-beta pairs.

The function also writes three files to "folder_out":

- A data frame ("\_pseudobulk_TRA.tsv") of pseudobulk counts for
  TCRalpha chains

- A data frame ("\_pseudobulk_TRB.tsv") of pseudobulk counts for TCRbeta
  chains

- A data frame ("\_TIRTLoutput.tsv") of TCR-alpha/TCR-beta pairs.
