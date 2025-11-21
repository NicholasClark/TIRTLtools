# Remove TCRs with unknown V-segments

**\[experimental\]** This function removes TCR pairs with V-segments
that are not found in
[`params`](https://nicholasclark.github.io/TIRTLtools/reference/params.md).
This is needed for pre-processing for
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
and is part of
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md),
which is automatically run during
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md).

## Usage

``` r
filter_v_alleles(df, params = NULL, verbose = TRUE)
```

## Arguments

- df:

  data frame with paired TCRs

- params:

  a table with acceptable V-segments (default is
  [`TIRTLtools::params`](https://nicholasclark.github.io/TIRTLtools/reference/params.md))

- verbose:

  whether to print number of TCRs removed

## Value

a paired TCR data frame with TCRs with unknown V-segments removed.

## See also

Other data_processing:
[`TIRTL_process()`](https://nicholasclark.github.io/TIRTLtools/reference/TIRTL_process.md),
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md),
[`filter_duplicate_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_duplicate_tcrs.md),
[`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md),
[`filter_nonfunctional_TCRs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_nonfunctional_TCRs.md),
[`filter_short_cdr3s()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_short_cdr3s.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)
