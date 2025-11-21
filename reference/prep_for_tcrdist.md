# Prepare paired TCRs for TCRdist calculation

**\[experimental\]** This function prepares a paired TCR data frame for
running TCRdist, adding alleles to V-segments (i.e. "\*01") if they are
missing and removing TCRs with unrecognized V-segments or non-functional
CDR3 sequences.

## Usage

``` r
prep_for_tcrdist(
  df,
  params = NULL,
  remove_duplicate_tcrs = FALSE,
  remove_MAIT = TRUE,
  verbose = TRUE
)
```

## Arguments

- df:

  a data frame of paired TCRs

- params:

  (optional) a data frame with permissible amino acids and V-segments.
  By default, the
  [`params`](https://nicholasclark.github.io/TIRTLtools/reference/params.md)
  data frame is used.

- remove_duplicate_tcrs:

  (optional) whether or not to remove duplicate TCRs (default is FALSE,
  do not remove duplicates).

- remove_MAIT:

  (optional) whether or not to remove MAIT cells from the data frame
  (default is TRUE).

- verbose:

  (optional) whether to print progress of the function (default is
  TRUE).

## Value

Returns a paired TCR data frame that is compatible with TCRdist python
scripts.

## See also

Other data_processing:
[`TIRTL_process()`](https://nicholasclark.github.io/TIRTLtools/reference/TIRTL_process.md),
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md),
[`filter_duplicate_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_duplicate_tcrs.md),
[`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md),
[`filter_nonfunctional_TCRs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_nonfunctional_TCRs.md),
[`filter_short_cdr3s()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_short_cdr3s.md),
[`filter_v_alleles()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_v_alleles.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)
