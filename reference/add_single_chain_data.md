# Add single-chain read counts/fractions to the paired TCR data

**\[experimental\]** The `add_single_chain_data()` function adds read
counts and proportions from the single-chain pseudobulk data to the
paired data frame for each sample of a dataset.

## Usage

``` r
add_single_chain_data(data, verbose = TRUE)
```

## Arguments

- data:

  a TIRTLseq dataset created by
  [`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)

- verbose:

  (optional) whether to print progress of the function (default is
  TRUE).

## Value

A TIRTLseq dataset object where the paired data frames for each sample
have added columns for read counts and proportions from the single-chain
pseudobulk data.

## See also

Other data_processing:
[`TIRTL_process()`](https://nicholasclark.github.io/TIRTLtools/reference/TIRTL_process.md),
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md),
[`filter_duplicate_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_duplicate_tcrs.md),
[`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md),
[`filter_nonfunctional_TCRs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_nonfunctional_TCRs.md),
[`filter_short_cdr3s()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_short_cdr3s.md),
[`filter_v_alleles()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_v_alleles.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)

## Examples

``` r
# example code

```
