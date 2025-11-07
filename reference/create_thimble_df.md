# Convert paired TCRs to a "thimble" data frame for stitching together TCR nucleotide coding sequences with the Stitchr Python package

**\[experimental\]**

## Usage

``` r
create_thimble_df(
  df,
  preset = c("default", "none"),
  TCR_names = "TCR",
  exclude_non_functional = TRUE,
  remove_duplicates = TRUE,
  verbose = TRUE,
  Linker = NULL,
  Link_order = NULL,
  TRAC = NULL,
  TRBC = NULL,
  TRA_leader = NULL,
  TRB_leader = NULL,
  TRA_5_prime_seq = NULL,
  TRA_3_prime_seq = NULL,
  TRB_5_prime_seq = NULL,
  TRB_3_prime_seq = NULL
)
```

## See also

Other stitchr:
[`write_stitchr_tsv()`](https://nicholasclark.github.io/TIRTLtools/reference/write_stitchr_tsv.md)
