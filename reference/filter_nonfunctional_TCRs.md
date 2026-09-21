# Remove TCRs with nonfunctional CDR3 amino acid sequences

**\[experimental\]**

This function removes TCR pairs where the CDR3-alpha or beta segment
does not make a functional protein, i.e. amino acid sequence contains a
stop codon (\*) or a frameshift (\_).

## Usage

``` r
filter_nonfunctional_TCRs(data, verbose = TRUE)
```

## Arguments

- data:

  a TIRTLseqDataSet object

- verbose:

  whether to print number of TCRs removed

## Value

Returns a TIRTLseqDataSet object with nonfunctional TCRs removed.

## See also

Other data_processing:
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md),
[`combine_bulk_and_paired_data()`](https://nicholasclark.github.io/TIRTLtools/reference/combine_bulk_and_paired_data.md),
[`filter_duplicate_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_duplicate_tcrs.md),
[`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md),
[`filter_short_cdr3s()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_short_cdr3s.md),
[`filter_v_alleles()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_v_alleles.md),
[`make_tcr_schema()`](https://nicholasclark.github.io/TIRTLtools/reference/make_tcr_schema.md),
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)
