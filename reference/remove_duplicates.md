# Removes duplicate paired TCRs

**\[experimental\]**

By default, paired TCRs are listed in one row in the paired chain data
for each pairing algorithm they were called by (either T-SHELL or
MAD-HYPE). This means some TCRs are present in two rows, once for each
algorithm. This function removes duplicate TCRs and returns a data frame
where each TCR pair is found in only one row.

## Usage

``` r
remove_duplicates(data)
```

## Arguments

- data:

  either a single data frame (paired chain) or a list of data frames

## Value

Returns a data frame (or list of data frames, depending on input) of
paired TCRs where each TCR is listed only once.

## See also

Other data_processing:
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md),
[`combine_bulk_and_paired_data()`](https://nicholasclark.github.io/TIRTLtools/reference/combine_bulk_and_paired_data.md),
[`filter_duplicate_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_duplicate_tcrs.md),
[`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md),
[`filter_nonfunctional_TCRs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_nonfunctional_TCRs.md),
[`filter_short_cdr3s()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_short_cdr3s.md),
[`filter_v_alleles()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_v_alleles.md),
[`make_tcr_schema()`](https://nicholasclark.github.io/TIRTLtools/reference/make_tcr_schema.md),
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md)
