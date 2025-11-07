# Removes duplicate paired TCRs

**\[experimental\]** By default, paired TCRs are listed twice in the
paired chain data if they were called by both the T-SHELL and MAD-HYPE
pairing algorithms. This function removes duplicate TCRs and returns a
data frame with only one of each pair.

## Usage

``` r
remove_duplicates(data)
```

## Arguments

- data:

  either a single data frame (paired chain) or a list of data frames

## See also

Other data_processing:
[`TIRTL_process()`](https://nicholasclark.github.io/TIRTLtools/reference/TIRTL_process.md),
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md),
[`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md)
