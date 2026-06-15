# Run data processing functions on a TIRTLseqDataSet object

**\[experimental\]** This function runs annotation and data cleaning
functions on a TIRTLseqDataSet object. Specifically, it calls the
functions
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
and (optionally)
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md).

## Usage

``` r
TIRTL_process(...)

process_TIRTLseq(data, clean = FALSE, remove_nonfunctional = FALSE)
```

## Arguments

- data:

  a TIRTLseqDataSet object

- clean:

  (optional) a TRUE/FALSE value, whether or not to "clean" the paired
  data by removing excess pairs for individual alpha and beta chains
  (default is FALSE).

- remove_nonfunctional:

  whether to remove non-functional TCR chains (default is FALSE)

## Value

a TIRTLseqDataSet object with annotated and (optionally) cleaned data

## See also

Other data_processing:
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md),
[`combine_bulk_and_paired_data()`](https://nicholasclark.github.io/TIRTLtools/reference/combine_bulk_and_paired_data.md),
[`filter_duplicate_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_duplicate_tcrs.md),
[`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md),
[`filter_nonfunctional_TCRs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_nonfunctional_TCRs.md),
[`filter_short_cdr3s()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_short_cdr3s.md),
[`filter_v_alleles()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_v_alleles.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`make_tcr_schema()`](https://nicholasclark.github.io/TIRTLtools/reference/make_tcr_schema.md),
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md),
[`read_external_bulk()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_bulk.md),
[`read_external_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_paired.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)
