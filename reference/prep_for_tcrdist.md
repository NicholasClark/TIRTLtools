# Prepare paired TCRs for TCRdist calculation

**\[experimental\]** This function prepares a paired TCR data frame for
running TCRdist, adding alleles to V-segments (i.e. "\*01") if they are
missing and removing TCRs with unrecognized V-segments or non-functional
CDR3 sequences.

## Usage

``` r
prep_for_tcrdist(df, params = NULL, remove_MAIT = TRUE, verbose = TRUE)
```

## Arguments

- df:

  a data frame of paired TCRs

- params:

  (optional) a data frame with permissible amino acids and V-segments.
  By default, the TIRTLtools::params data frame is used.

- remove_MAIT:

  (optional) whether or not to remove MAIT cells from the data frame
  (default is TRUE).

- verbose:

  (optional) whether to print progress of the function (default is
  TRUE).

## See also

Other data_processing:
[`TIRTL_process()`](https://nicholasclark.github.io/TIRTLtools/reference/TIRTL_process.md),
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md),
[`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)
