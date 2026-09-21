# Run data processing functions on a TIRTLseqDataSet object

**\[deprecated\]** This function has been deprecated because its
functionality has been added to the load_tirtlseq() function

This function runs annotation and data cleaning functions on a
TIRTLseqDataSet object. Specifically, it calls the functions
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
and (optionally)
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md).

## Usage

``` r
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

Other deprecated:
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md)
