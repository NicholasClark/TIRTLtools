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
  clean = FALSE,
  remove_nonfunctional = FALSE,
  process = TRUE,
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

- clean:

  (optional) a TRUE/FALSE value, whether or not to "clean" the paired
  data by removing excess pairs for individual alpha and beta chains
  (default is FALSE).

- remove_nonfunctional:

  whether to remove non-functional TCR chains (default is FALSE)

- process:

  whether to process the raw data, adding columns for alpha and beta
  readcounts to the paired data frame, adding columns to the pseudobulk
  data to identify which chains are paired, etc.

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

## See also

Other data_import:
[`get_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/get_samples.md),
[`load_example_data()`](https://nicholasclark.github.io/TIRTLtools/reference/load_example_data.md),
[`read_external_bulk()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_bulk.md),
[`read_external_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_paired.md)

## Examples

``` r
load_example_data(dataset = "SJTRC_minimal")
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
SJTRC_minimal$meta ## metadata table
#> # A tibble: 2 × 5
#>   sample_id  marker timepoint version label                                     
#>   <chr>      <chr>  <chr>     <chr>   <chr>                                     
#> 1 cd8_tp1_v2 cd8    tp1       v2      marker: cd8 | timepoint: tp1 | version: v2
#> 2 cd8_tp2_v2 cd8    tp2       v2      marker: cd8 | timepoint: tp2 | version: v2
names(SJTRC_minimal$data) ## samples
#> [1] "cd8_tp1_v2" "cd8_tp2_v2"
SJTRC_minimal$data$cd8_tp1_v2$alpha %>% dplyr::as_tibble() ## alpha chain pseudobulk
#> # A tibble: 937,875 × 23
#>    id        rank aaSeqCDR3           v     j     readFraction readCount n_wells
#>    <chr>    <int> <chr>               <chr> <chr>        <dbl>     <int>   <int>
#>  1 alpha_1      1 CAMRGRSNRDDKIIF     TRAV… TRAJ…      0.0206    1411609     191
#>  2 alpha_2      2 CALKTSYDKVIF        TRAV… TRAJ…      0.0148    1016256     191
#>  3 alpha_3      3 CAASRLPDDMRF        TRAV… TRAJ…      0.00956    655369     191
#>  4 alpha_4      4 CAVRDSNYQLIW        TRAV… TRAJ…      0.00548    375648     191
#>  5 alpha_5      5 CLLAHLRSGPNSGGSNYK… TRAV… TRAJ…      0.00470    322182     191
#>  6 alpha_6      6 CAYRTPGETGANNLFF    TRAV… TRAJ…      0.00452    309912     191
#>  7 alpha_7      7 CALSEVPYRGNTPLVF    TRAV… TRAJ…      0.00384    263000     191
#>  8 alpha_8      8 CAVSDSNYQLIW        TRAV… TRAJ…      0.00367    251464     191
#>  9 alpha_9      9 CAVMDSNYQLIW        TRAV… TRAJ…      0.00361    247566     191
#> 10 alpha_10    10 CAMREGYNTDKLIF      TRAV… TRAJ…      0.00330    226187     191
#> # ℹ 937,865 more rows
#> # ℹ 15 more variables: paired_status <chr>, n_paired <dbl>,
#> #   targetSequences <chr>, readCount_median <dbl>, readCount_max <int>,
#> #   sem <dbl>, max_wells <int>, n_paired_madhype <dbl>, n_paired_tshell <dbl>,
#> #   has_stop_codon <lgl>, has_frameshift <lgl>, is_functional <lgl>,
#> #   is_paired <lgl>, is_paired_tshell <lgl>, is_paired_madhype <lgl>
SJTRC_minimal$data$cd8_tp1_v2$beta %>% dplyr::as_tibble() ## beta chain pseudobulk
#> # A tibble: 928,777 × 23
#>    id       rank aaSeqCDR3         v        j     readFraction readCount n_wells
#>    <chr>   <int> <chr>             <chr>    <chr>        <dbl>     <int>   <int>
#>  1 beta_1      1 CSVERVNYNEQFF     TRBV29-1 TRBJ…      0.0135    1325768     191
#>  2 beta_2      2 CASSTGLPRDTQYF    TRBV5-8  TRBJ…      0.0114    1124522     191
#>  3 beta_3      3 CSATAGRNYGYTF     TRBV20-1 TRBJ…      0.0107    1048157     191
#>  4 beta_4      4 CASSLTYEQYF       TRBV6-2  TRBJ…      0.00841    827268     191
#>  5 beta_5      5 CASSWPQGSGSLDEQFF TRBV5-1  TRBJ…      0.00518    509093     191
#>  6 beta_6      6 CASSTTWGTGELFF    TRBV7-8  TRBJ…      0.00491    483264     191
#>  7 beta_7      7 CSVEQNTEAFF       TRBV29-1 TRBJ…      0.00491    482789     191
#>  8 beta_8      8 CASSKEGASPLHF     TRBV5-1  TRBJ…      0.00456    448792     191
#>  9 beta_9      9 CAWSKGEANVLTF     TRBV30   TRBJ…      0.00443    435608     191
#> 10 beta_10    10 CASSYEGSWGDEQFF   TRBV6-2  TRBJ…      0.00396    389542     191
#> # ℹ 928,767 more rows
#> # ℹ 15 more variables: paired_status <chr>, n_paired <dbl>,
#> #   targetSequences <chr>, readCount_median <dbl>, readCount_max <int>,
#> #   sem <dbl>, max_wells <int>, n_paired_madhype <dbl>, n_paired_tshell <dbl>,
#> #   has_stop_codon <lgl>, has_frameshift <lgl>, is_functional <lgl>,
#> #   is_paired <lgl>, is_paired_tshell <lgl>, is_paired_madhype <lgl>
SJTRC_minimal$data$cd8_tp1_v2$paired %>% dplyr::as_tibble() ## all pairs (rows duplicated if called by both pairing algorithms)
#> # A tibble: 20,690 × 47
#>    method va         ja     cdr3a         cdr3b vb    jb    alpha_rank beta_rank
#>    <chr>  <chr>      <chr>  <chr>         <chr> <chr> <chr>      <int>     <int>
#>  1 tshell TRAV13-1   TRAJ32 CAAKGVYGGATN… CSVE… TRBV… TRBJ…         12         1
#>  2 tshell TRAV10     TRAJ10 CVVNHGR_GNKL… CSVE… TRBV… TRBJ…         79         1
#>  3 tshell TRAV14/DV4 TRAJ30 CAMRGRSNRDDK… CASS… TRBV… TRBJ…          1         2
#>  4 tshell TRAV12-2   TRAJ23 CAVKQ*P_GGKL… CASS… TRBV… TRBJ…         30         2
#>  5 tshell TRAV29/DV5 TRAJ43 CAASRLPDDMRF  CSAT… TRBV… TRBJ…          3         3
#>  6 tshell TRAV19     TRAJ29 CALSEVPYRGNT… CASS… TRBV… TRBJ…          7         4
#>  7 tshell TRAV2      TRAJ42 CAVEASLW_GSQ… CASS… TRBV… TRBJ…         99         4
#>  8 tshell TRAV8-1    TRAJ32 CAVDATNKLIF   CASS… TRBV… TRBJ…         28         5
#>  9 tshell TRAV30     TRAJ48 CGTDAFISNFGN… CASS… TRBV… TRBJ…         67         6
#> 10 tshell TRAV12-2   TRAJ24 CAVNTISDSWGK… CSVE… TRBV… TRBJ…         18         7
#> # ℹ 20,680 more rows
#> # ℹ 38 more variables: alpha_readFraction <dbl>, beta_readFraction <dbl>,
#> #   alpha_readCount <int>, beta_readCount <int>, is_functional <lgl>,
#> #   alpha_nuc <chr>, beta_nuc <chr>, alpha_beta <chr>, wi <int>, wj <int>,
#> #   wij <int>, wa <int>, wb <int>, score <dbl>, r <dbl>, ts <dbl>, pval <dbl>,
#> #   pval_adj <dbl>, loss_a_frac <dbl>, loss_b_frac <dbl>,
#> #   alpha_readCount_max <int>, alpha_readCount_median <dbl>, alpha_sem <dbl>, …
SJTRC_minimal$data$cd8_tp1_v2$paired_alt %>% dplyr::as_tibble() ## all pairs (each row unique)
#> # A tibble: 14,273 × 49
#>    paired_status va         ja     cdr3a  cdr3b vb    jb    alpha_rank beta_rank
#>    <chr>         <chr>      <chr>  <chr>  <chr> <chr> <chr>      <int>     <int>
#>  1 T-SHELL only  TRAV13-1   TRAJ32 CAAKG… CSVE… TRBV… TRBJ…         12         1
#>  2 T-SHELL only  TRAV10     TRAJ10 CVVNH… CSVE… TRBV… TRBJ…         79         1
#>  3 T-SHELL only  TRAV14/DV4 TRAJ30 CAMRG… CASS… TRBV… TRBJ…          1         2
#>  4 T-SHELL only  TRAV12-2   TRAJ23 CAVKQ… CASS… TRBV… TRBJ…         30         2
#>  5 T-SHELL only  TRAV29/DV5 TRAJ43 CAASR… CSAT… TRBV… TRBJ…          3         3
#>  6 T-SHELL only  TRAV19     TRAJ29 CALSE… CASS… TRBV… TRBJ…          7         4
#>  7 T-SHELL only  TRAV2      TRAJ42 CAVEA… CASS… TRBV… TRBJ…         99         4
#>  8 T-SHELL only  TRAV8-1    TRAJ32 CAVDA… CASS… TRBV… TRBJ…         28         5
#>  9 T-SHELL only  TRAV30     TRAJ48 CGTDA… CASS… TRBV… TRBJ…         67         6
#> 10 T-SHELL only  TRAV12-2   TRAJ24 CAVNT… CSVE… TRBV… TRBJ…         18         7
#> # ℹ 14,263 more rows
#> # ℹ 40 more variables: alpha_readFraction <dbl>, beta_readFraction <dbl>,
#> #   alpha_readCount <int>, beta_readCount <int>, is_functional <lgl>,
#> #   alpha_nuc <chr>, beta_nuc <chr>, alpha_beta <chr>, wi <int>, wj <int>,
#> #   wij <int>, wa <int>, wb <int>, score <dbl>, r <dbl>, ts <dbl>, pval <dbl>,
#> #   pval_adj <dbl>, loss_a_frac <dbl>, loss_b_frac <dbl>,
#> #   alpha_readCount_max <int>, alpha_readCount_median <dbl>, alpha_sem <dbl>, …
SJTRC_minimal$info ## information
#> $data_version
#> [1] "v1"
#> 
#> $processing_version
#> [1] "v1"
#> 
#> $package_version
#> [1] ‘0.2.1’
#> 
#> $call
#> $call$directory
#> [1] "/Users/nclark2/git/TIRTLtools/inst/extdata/SJTRC_TIRTLseq_minimal"
#> 
#> $call$chain
#> [1] "all"    "paired" "alpha"  "beta"  
#> 
#> $call$sep
#> [1] "_"
#> 
#> $call$meta_columns
#> [1] "marker"    "timepoint" "version"  
#> 
#> $call$samples
#> NULL
#> 
#> $call$clean
#> [1] FALSE
#> 
#> $call$remove_nonfunctional
#> [1] FALSE
#> 
#> $call$process
#> [1] TRUE
#> 
#> $call$pseudobulk_columns
#> [1] "auto"
#> 
#> $call$paired_columns
#> [1] "auto"
#> 
#> $call$n_threads
#> [1] 9
#> 
#> $call$verbose
#> [1] TRUE
#> 
#> $call$stringsAsFactors
#> [1] FALSE
#> 
#> $call$n_max
#> [1] Inf
#> 
#> 
```
