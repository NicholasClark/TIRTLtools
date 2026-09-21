# Remove duplicate TCRs from a data frame

**\[experimental\]**

By default, the TIRTLseq pairing scripts write a "paired" data table
that has one row per pair per pairing algorithm (MAD-HYPE or T-SHELL).
Thus, any pairs identified by both methods will be found in two rows.
This function removes those duplicates.

## Usage

``` r
filter_duplicate_tcrs(df, verbose = TRUE)
```

## Arguments

- df:

  data frame with paired TCRs

- verbose:

  whether to print number of TCRs removed

## See also

Other data_processing:
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md),
[`combine_bulk_and_paired_data()`](https://nicholasclark.github.io/TIRTLtools/reference/combine_bulk_and_paired_data.md),
[`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md),
[`filter_nonfunctional_TCRs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_nonfunctional_TCRs.md),
[`filter_short_cdr3s()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_short_cdr3s.md),
[`filter_v_alleles()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_v_alleles.md),
[`make_tcr_schema()`](https://nicholasclark.github.io/TIRTLtools/reference/make_tcr_schema.md),
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)
