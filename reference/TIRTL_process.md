# Run data processing functions on a TIRTLseqData object

**\[experimental\]** This function runs annotation and data cleaning
functions on a TIRTLseqData object. Specifically, it calls the functions
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
and (optionally)
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md).

## Usage

``` r
TIRTL_process(data, clean = FALSE)
```

## Arguments

- data:

  a TIRTLseqData object

- clean:

  (optional) a TRUE/FALSE value, whether or not to "clean" the paired
  data by removing excess pairs for individual alpha and beta chains
  (default is FALSE).

## See also

Other data_processing:
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md),
[`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)
