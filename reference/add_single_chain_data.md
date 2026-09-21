# Add single-chain read counts/fractions to the paired TCR data

**\[deprecated\]** This function has been deprecated because its
functionality has been added to the load_tirtlseq() function

The `add_single_chain_data()` function adds read counts and proportions
from the single-chain pseudobulk data to the paired data frame for each
sample of a dataset.

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

Other deprecated:
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`process_TIRTLseq()`](https://nicholasclark.github.io/TIRTLtools/reference/process_TIRTLseq.md)
