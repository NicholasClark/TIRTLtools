# Convert paired TCRs to a "thimble" data frame for Stitchr

**\[experimental\]** This function creates a "thimble" data frame for
use with the Stitchr Python package
(https://jamieheather.github.io/stitchr/) to stitch together TCR
nucleotide coding sequences. The resulting data frame can be written to
a .tsv with
[`write_stitchr_tsv()`](https://nicholasclark.github.io/TIRTLtools/reference/write_stitchr_tsv.md).

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

## Arguments

- df:

  a dataframe with paired TCRs.

- preset:

  a preset for linker, link order, and TRAC/TRBC settings. The default
  ("default") setting uses "T2A" for the linker, "BA" for link order,
  "mTRAC\\01" for TRAC, and "mTRBC2\\01" for TRBC.

- TCR_names:

  a vector of names for the TCRs. If a vector of length one, then TCR
  names will be "\<TCR_names\>\_1", "\<TCR_names\>\_2", ... ,
  "\<TCR_names\>\_N".

- exclude_non_functional:

  whether to exclude TCRs with non-functional CDR3 sequences (default is
  TRUE).

- remove_duplicates:

  whether to remove any duplicate TCRs (default is TRUE).

- verbose:

  whether to print messages

- Linker:

  the linker for TCRs (default "T2A")

- Link_order:

  the link order for TCRs (default "BA")

- TRAC:

  TCR alpha chain constant region (default "mTRAC\*01")

- TRBC:

  TCR beta chain constant region (default "mTRBC2\*01")

- TRA_leader:

  (optional) An alternative leader sequence for the alpha chain. This
  can be either a specified gene entry in the pre-programmed IMGT data
  or alternatively a simple DNA string can be entered (e.g. ‘ATG’ for a
  minimal start codon in the absence of a leader sequence).

- TRB_leader:

  (optional) An alternative leader sequence for the beta chain. This can
  be either a specified gene entry in the pre-programmed IMGT data or
  alternatively a simple DNA string can be entered (e.g. ‘ATG’ for a
  minimal start codon in the absence of a leader sequence).

- TRA_5_prime_seq:

  Optional arbitrary sequence to be appended to the 5’ of the alpha
  chain.

- TRA_3_prime_seq:

  Optional arbitrary sequence to be appended to the 3’ of the alpha
  chain.

- TRB_5_prime_seq:

  Optional arbitrary sequence to be appended to the 5’ of the beta
  chain.

- TRB_3_prime_seq:

  Optional arbitrary sequence to be appended to the 3’ of the beta
  chain.

## Value

A data frame with TCRs properly formatted for use with Stitchr.

## See also

Other stitchr:
[`write_stitchr_tsv()`](https://nicholasclark.github.io/TIRTLtools/reference/write_stitchr_tsv.md)
