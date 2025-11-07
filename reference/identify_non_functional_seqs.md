# Identify TCRs that contain non-functional CDR3 sequences

**\[experimental\]** `identify_non_functional_seqs()` identifies CDR3
nucleotide sequences in the data that contain either stop codons (\*) or
frame shifts (\_) that would indicate a non-functional protein product.

## Usage

``` r
identify_non_functional_seqs(data)
```

## Arguments

- data:

  a TIRTLseqData object

## Value

A TIRTLseqData object with modified pseudobulk and paired data frames
for each sample. Each dataframe in the ouptut object has added columns
that identify whether the CDR3 alpha and beta nucleotide sequences
contain any stop codons (\*) or frame shifts (\_) that would indicate a
non-functional chain.

If column `is_functional` is TRUE if neither chain has a stop codon or a
frame shift. The columns `has_stop_codon` and `has_frameshift` are
similar, but specific to each kind of coding error. Other columns
identify if the alpha chain or beta chain has a stop codon or
frameshift, and if it is functional.

## See also

[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`plot_num_partners()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_num_partners.md)

Other data_processing:
[`TIRTL_process()`](https://nicholasclark.github.io/TIRTLtools/reference/TIRTL_process.md),
[`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md),
[`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md),
[`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md),
[`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md),
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)

## Examples

``` r
# example code

```
