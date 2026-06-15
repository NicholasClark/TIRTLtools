# Make a "schema" for defining a T-Cell Receptor

**\[experimental\]** This is a convenience function that simply returns
a vector of column names that are used by
[`read_external_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_paired.md)
to define a unique T-Cell receptor.

For example, using `features = c("v", "cdr3_aa")` will group all chains
with the same V-alpha/beta and CDR3-alpha/beta amino acid sequence
together, while using `features = c("v", "cdr3_aa", "cdr3_nt")` will
require that their nucleotide sequence is also matching.

## Usage

``` r
make_tcr_schema(
  features = c("v", "j", "cdr3_aa", "cdr3_nt"),
  second_alpha = FALSE
)
```

## Arguments

- features:

  some subset of `c("v", "j", "cdr3_aa", "cdr3_nt")` that you would like
  to use to define a unique receptor.

- second_alpha:

  whether to include a second alpha chain (default is FALSE)

## Value

A vector containing

## See also

Other data_processing:
[`TIRTL_process()`](https://nicholasclark.github.io/TIRTLtools/reference/process_TIRTLseq.md),
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
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md),
[`read_external_bulk()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_bulk.md),
[`read_external_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_paired.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)

## Examples

``` r
make_tcr_schema() ## default: v, j, and cdr3 amino acid and nucleotide sequence
#> [1] "va"        "ja"        "cdr3a"     "cdr3b"     "vb"        "jb"       
#> [7] "alpha_nuc" "beta_nuc" 
make_tcr_schema(c("v", "j")) ## only v and j segments
#> [1] "va" "ja" "vb" "jb"
```
