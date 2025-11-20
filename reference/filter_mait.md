# Remove MAIT (Mucosal-Associated Invariant T cells) TCRs

**\[experimental\]** This function uses a heuristic to remove TCRs
associated with MAIT cells, which may not be desired for some
applications. For example, when creating a network of related TCRs with
TCRdist, these TCRs often form a very highly dense sub-network that
inflates output file size.

## Usage

``` r
filter_mait(df, verbose = TRUE)
```

## Arguments

- df:

  data frame with paired TCRs

- verbose:

  whether to print number of MAIT TCRs removed

## Value

A data frame with MAIT TCRs removed

## Details

We use a heiristic where we remove TCRs with V-alpha TRAV1-2 and either
J-alpha segment TRAJ33, TRAJ12, or TRAJ20 (see Garner et al. 2023 -
https://www.nature.com/articles/s41590-023-01575-1).

## See also

Other data_processing:
[`TIRTL_process()`](https://nicholasclark.github.io/TIRTLtools/reference/TIRTL_process.md),
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md),
[`filter_duplicate_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_duplicate_tcrs.md),
[`filter_nonfunctional_TCRs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_nonfunctional_TCRs.md),
[`filter_short_cdr3s()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_short_cdr3s.md),
[`filter_v_alleles()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_v_alleles.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)
