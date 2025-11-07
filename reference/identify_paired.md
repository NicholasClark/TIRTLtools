# Identify which single chains were paired

**\[experimental\]** For each sample in the dataset, `identify_paired()`
annotates the alpha and beta pseudobulk data with the number of distinct
pairs each chain is a part of in the paired data as well as a TRUE/FALSE
column indicating whether the chain is paired with any partners.

## Usage

``` r
identify_paired(data, verbose = TRUE)
```

## Arguments

- data:

  a TIRTLseq dataset created by
  [`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)

- verbose:

  (optional) whether to print progress of the function (default is
  TRUE).

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

Other data_processing:
[`TIRTL_process()`](https://nicholasclark.github.io/TIRTLtools/reference/TIRTL_process.md),
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md),
[`filter_duplicate_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_duplicate_tcrs.md),
[`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md),
[`filter_nonfunctional_TCRs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_nonfunctional_TCRs.md),
[`filter_short_cdr3s()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_short_cdr3s.md),
[`filter_v_alleles()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_v_alleles.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)

## Examples

``` r
# example code
# paired = load_tirtlseq("path_to/your_directory", sep = "_", meta_columns = c("cell_type", "timepoint"))
# paired = identify_paired(paired)
```
