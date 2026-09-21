# Identify which single chains were paired

**\[deprecated\]**

This function has been deprecated because its functionality has been
added to the load_tirtlseq() function

For each sample in the dataset, `identify_paired()` annotates the alpha
and beta pseudobulk data with the number of distinct pairs each chain is
a part of in the paired data as well as a TRUE/FALSE column indicating
whether the chain is paired with any partners.

## Usage

``` r
identify_paired(data, verbose = TRUE, by_method = TRUE)
```

## Arguments

- data:

  a TIRTLseq dataset created by
  [`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)

- verbose:

  (optional) whether to print progress of the function (default is
  TRUE).

- by_method:

  (optional) whether to get stats for each pairing method

## Value

A dataset similar to that created by
[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md),
but with added columns `is_paired` and `n_paired` in the alpha and beta
pseudobulk data frames.

`is_paired` is TRUE if the chain is found in the paired data. `n_paired`
is the number of distinct chains that the particular chain is paired
with.

## See also

[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)

Other deprecated:
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`process_TIRTLseq()`](https://nicholasclark.github.io/TIRTLtools/reference/process_TIRTLseq.md)

## Examples

``` r
# example code
# paired = load_tirtlseq("path_to/your_directory", sep = "_", meta_columns = c("cell_type", "timepoint"))
# paired = identify_paired(paired)
```
