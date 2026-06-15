# Read and process single-cell paired-chain TCR-seq data

**\[experimental\]** This function reads and processes paired
TCR-sequencing data from a non-TIRTLseq assay. Currently 10X and Parse
Biosciences data are supported.

## Usage

``` r
read_external_paired(
  path,
  format = c("auto", "10X", "ParseBio"),
  id_cols = make_tcr_schema(features = c("v", "j", "cdr3_aa", "cdr3_nt"), second_alpha =
    FALSE),
  multi = FALSE,
  separate_rows = TRUE,
  productive_only = TRUE,
  verbose = TRUE
)
```

## Arguments

- path:

  the path to the data file

- format:

  the format of the data, either `"10X"` or `"ParseBio"` or `"auto"`. If
  `"auto"` (default), the function will try to decipher which technology
  the data file was created by.

- id_cols:

  a vector of column names to be used to define a clone (e.g.
  `c("va","vb","ja","jb","alpha_nuc", "beta_nuc", "cdr3a", "cdr3b")`).
  This can be produced automatically using the
  [`make_tcr_schema()`](https://nicholasclark.github.io/TIRTLtools/reference/make_tcr_schema.md)
  function.

- multi:

  If `FALSE` (default), select only the best two alpha chains for each
  beta chain when processing the data and creating a data frame with
  paired TCRs. If `TRUE`, keep all alphas.

- separate_rows:

  If `TRUE`, when there are multiple alpha chains paired with one beta
  chain, put each pair in a separate row in the output data frame. If
  `FALSE`, add second alpha chain in extra columns on the same row.

- productive_only:

  If `TRUE`, keep only "productive" chains

- verbose:

  If `TRUE`, print messages.

## Value

A list containing the following slots:

- df_pairs_complete - (data frame) paired receptors - one row for each
  receptor (excluding those missing an alpha or beta chain)

- df_pairs - (data frame) with paired receptors - one row for each
  receptor (including those missing an alpha or beta chain)

- df_pairs_long - (data frame) with paired receptors - one row for each
  cell (including those missing an alpha or beta chain)

- df_pairs_long_complete - (data frame) with paired receptors - one row
  for each cell (excluding those missing an alpha or beta chain)

- df_raw - (data frame) un-edited input data

- chain_df - (data frame) summary of total number of each chain in input
  data

- barcode_df - (data frame) summary of number of chains found in each
  cell

- id_cols - (character vector) columns used as IDs to uniquely define
  receptor pairs

- n_cells_total -(integer) total number of cells

- n_cells_complete - (integer) number of cells with both chains

## Details

Supported data types:

- `"10X"` - "filtered_contig_annotations.csv" or
  "all_contig_annotations.csv" outputs from Cell Ranger
  (https://www.10xgenomics.com/support/software/cell-ranger/7.2/analysis/outputs/cr-5p-outputs-overview-vdj)

- `"ParseBio"` - "tcr_annotation_airr.tsv" output from Trailmaker

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
[`make_tcr_schema()`](https://nicholasclark.github.io/TIRTLtools/reference/make_tcr_schema.md),
[`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md),
[`read_external_bulk()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_bulk.md),
[`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)
