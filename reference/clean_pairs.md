# Remove excess pairs for individual single chains

**\[experimental\]** This function filters the paired TCR data, keeping
at most one beta-chain partner for each individual alpha-chain and at
most two alpha-chain partners for each individual beta-chain.

## Usage

``` r
clean_pairs(data, n_max_alpha = 2, n_max_beta = 1, verbose = TRUE)
```

## Arguments

- data:

  a TIRTLseqData object

- n_max_alpha:

  (optional) the maximum number of alpha chains allowed paired with an
  individual beta chain (default 2)

- n_max_beta:

  (optional) the maximum number of beta chains allowed paired with an
  individual alpha chain (default 1)

- verbose:

  (optional) whether to print progress of the function (default is
  TRUE).

## Details

Excess partners for an individual chain are often sequencing errors of
the true partner. The sequencing error-derived chains often share the
same V/J segment as the true partner, but are found at much lower read
fractions.

As a heuristic, to mitigate this phenomenon, for each unique beta chain
we group the partner alpha chains by their V/J segments and keep only
the alpha chains with the highest read fraction for each group. Out of
the remaining alpha chain partners, we keep up to two chains (can be
changed with n_max_alpha) with the highest read fractions.

We go through a similar process with each unique alpha chain, grouping
partner beta chains by their V/J segments. However, we keep only one
beta chain (can be changed with n_max_beta).

## See also

Other data_processing:
[`TIRTL_process()`](https://nicholasclark.github.io/TIRTLtools/reference/TIRTL_process.md),
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`filter_duplicate_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_duplicate_tcrs.md),
[`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md),
[`filter_nonfunctional_TCRs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_nonfunctional_TCRs.md),
[`filter_short_cdr3s()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_short_cdr3s.md),
[`filter_v_alleles()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_v_alleles.md),
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)
