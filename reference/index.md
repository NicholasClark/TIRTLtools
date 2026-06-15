# Package index

## TCR-alpha/beta pairing

- [`get_tshell_settings()`](https://nicholasclark.github.io/TIRTLtools/reference/get_tshell_settings.md)
  **\[experimental\]** : Get p-value and well threshold settings for
  T-SHELL
- [`run_pairing()`](https://nicholasclark.github.io/TIRTLtools/reference/run_pairing.md)
  **\[experimental\]** : Find TCRalpha/beta pairs from individual well
  read counts

## Data loading

- [`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)
  **\[experimental\]** : Load data from TIRTLseq experiments

## Data wrangling

- [`add_metadata()`](https://nicholasclark.github.io/TIRTLtools/reference/add_metadata.md)
  **\[experimental\]** : Add metadata to a TIRTLseqDataSet object
- [`filter_dataset()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_dataset.md)
  **\[experimental\]** : Subset a TIRTLseqDataSet object
- [`reorder_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/reorder_samples.md)
  **\[experimental\]** : Re-order samples in a TIRTLseqDataSet object

## Data processing and cleaning

- [`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md)
  **\[experimental\]** : Add single-chain read counts/fractions to the
  paired TCR data
- [`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md)
  **\[experimental\]** : Remove excess pairs for individual single
  chains
- [`combine_bulk_and_paired_data()`](https://nicholasclark.github.io/TIRTLtools/reference/combine_bulk_and_paired_data.md)
  **\[experimental\]** : Combine bulk single-chain TCR data with paired
  TCR data
- [`filter_duplicate_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_duplicate_tcrs.md)
  **\[experimental\]** : Remove duplicate TCRs from a data frame
- [`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md)
  **\[experimental\]** : Remove MAIT (Mucosal-Associated Invariant T
  cells) TCRs
- [`filter_nonfunctional_TCRs()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_nonfunctional_TCRs.md)
  **\[experimental\]** : Remove TCRs with nonfunctional CDR3 amino acid
  sequences
- [`filter_short_cdr3s()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_short_cdr3s.md)
  **\[experimental\]** : Remove TCRs with short CDR3 loops
- [`filter_v_alleles()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_v_alleles.md)
  **\[experimental\]** : Remove TCRs with unknown V-segments
- [`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md)
  **\[experimental\]** : Identify TCRs that contain non-functional CDR3
  sequences
- [`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md)
  **\[experimental\]** : Identify which single chains were paired
- [`make_tcr_schema()`](https://nicholasclark.github.io/TIRTLtools/reference/make_tcr_schema.md)
  **\[experimental\]** : Make a "schema" for defining a T-Cell Receptor
- [`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md)
  **\[experimental\]** : Prepare paired TCRs for TCRdist calculation
- [`TIRTL_process()`](https://nicholasclark.github.io/TIRTLtools/reference/process_TIRTLseq.md)
  [`process_TIRTLseq()`](https://nicholasclark.github.io/TIRTLtools/reference/process_TIRTLseq.md)
  **\[experimental\]** : Run data processing functions on a
  TIRTLseqDataSet object
- [`read_external_bulk()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_bulk.md)
  **\[experimental\]** : Read and process bulk single-chain TCR-seq data
- [`read_external_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/read_external_paired.md)
  **\[experimental\]** : Read and process single-cell paired-chain
  TCR-seq data
- [`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)
  **\[experimental\]** : Removes duplicate paired TCRs

## Quality control and data summary

- [`get_all_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_tcrs.md)
  **\[experimental\]** : Returns all of the paired TCRs from all samples
  in a dataset
- [`get_pair_stats()`](https://nicholasclark.github.io/TIRTLtools/reference/get_pair_stats.md)
  **\[experimental\]** : Count the number of pairs called by each
  algorithm
- [`get_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/get_paired_by_read_fraction_range.md)
  **\[experimental\]** : Calculate the number and fraction of single
  chains that were paired by frequency
- [`plot_n_reads()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_n_reads.md)
  **\[experimental\]** : Bar plot of the number of single-chain reads
  for each sample
- [`plot_num_partners()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_num_partners.md)
  **\[experimental\]** : Stacked bar plot of the fraction of alpha/beta
  chains with different numbers of partners
- [`plot_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired.md)
  **\[experimental\]** : Stacked bar plot of the number of TCRs paired
  by each algorithm
- [`plot_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_by_read_fraction_range.md)
  **\[experimental\]** : Bar plot of the fraction of paired single
  chains by frequency
- [`plot_paired_vs_rank()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_vs_rank.md)
  **\[experimental\]** : A step plot of the cumulative number of
  paired/unpaired alpha/beta chains among the most frequent chains
- [`plot_pairs_with_eachother()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_pairs_with_eachother.md)
  **\[experimental\]** : A connected point plot of read fraction vs.
  rank for the most frequent alpha/beta chains
- [`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md)
  [`rank_plot()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md)
  **\[experimental\]** : Line plot of clonotype rank vs. read fraction
  for each sample
- [`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md)
  **\[experimental\]** : A point plot of read fraction vs. rank for the
  most frequent alpha/beta chains
- [`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md)
  [`sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md)
  **\[experimental\]** : Plot the overlap/agreement between samples (in
  terms of most frequent clones)
- [`summarize_data()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_data.md)
  **\[experimental\]** : Create a summary table with number of reads and
  unique alpha/beta chains

## TCR similarity and clustering

- [`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
  **\[experimental\]** : GPU implementation of TCRdist, a
  distance/similarity metric for pairs of TCRs
- [`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md)
  **\[experimental\]** : Parallelized C++ implementation of TCRdist (no
  GPU required)
- [`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md)
  **\[experimental\]** : Cluster TCRs (using the Leiden algorithm) based
  on their pairwise TCRdist values
- [`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md)
  **\[experimental\]** : Plot clusters of similar TCRs

## Longitudinal analysis

- [`plot_clone_size_across_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clone_size_across_samples.md)
  **\[experimental\]** : Line plot of clone read fraction across
  multiple samples
- [`plot_sample_vs_sample()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_vs_sample.md)
  **\[experimental\]** : Scatterplot of TCR clone read fraction of
  clones between two samples

## Diversity

- [`calculate_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/calculate_diversity.md)
  [`diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/calculate_diversity.md)
  **\[experimental\]** : Calculate TCR repertoire diversity metrics

- [`get_all_div_metrics()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_div_metrics.md)
  **\[experimental\]** :

  Returns all diversity metric options for
  [`calculate_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/calculate_diversity.md)

- [`plot_clonotype_indices()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clonotype_indices.md)
  **\[experimental\]** : Stacked bar chart with fractions of reads
  attributed to the most frequent clonotypes

- [`plot_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_diversity.md)
  **\[experimental\]** : Bar plot of clonal diversity metrics

## Export for external tools

- [`create_thimble_df()`](https://nicholasclark.github.io/TIRTLtools/reference/create_thimble_df.md)
  **\[experimental\]** : Convert paired TCRs to a "thimble" data frame
  for Stitchr
- [`write_stitchr_tsv()`](https://nicholasclark.github.io/TIRTLtools/reference/write_stitchr_tsv.md)
  **\[experimental\]** : Write a tab-separated "thimble" text file for
  use with Stitchr

## Well-level data functions

- [`choose_pair_manual()`](https://nicholasclark.github.io/TIRTLtools/reference/choose_pair_manual.md)
  **\[experimental\]** : Choose a partner manually for an input chain
  using T-SHELL
- [`get_well_subset()`](https://nicholasclark.github.io/TIRTLtools/reference/get_well_subset.md)
  **\[experimental\]** : Get well names from numerical rows and columns
- [`get_wells_from_edges()`](https://nicholasclark.github.io/TIRTLtools/reference/get_wells_from_edges.md)
  **\[experimental\]** : Get well names from the top left and bottom
  right wells
- [`load_well_counts_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/load_well_counts_binary.md)
  **\[experimental\]** : Load well-level data written to binary format
- [`plot_tshell()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_tshell.md)
  **\[experimental\]** : Calculate and plot read fraction correlation
  across wells (T-SHELL)
- [`write_well_data_to_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/write_well_data_to_binary.md)
  **\[experimental\]** : Convert individual well .tsv files to faster
  loading binary format

## Single-cell TIRTLseq data

- [`plot_scTIRTLseq_plate()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_scTIRTLseq_plate.md)
  **\[experimental\]** : Plot a plate of single-cell TIRTL-seq data
- [`process_scTIRTLseq()`](https://nicholasclark.github.io/TIRTLtools/reference/process_scTIRTLseq.md)
  **\[experimental\]** : Process single-cell TIRTL-seq data

## Data

- [`params`](https://nicholasclark.github.io/TIRTLtools/reference/params.md)
  : Table of permissible amino acids and V-segments for TCRdist
- [`submat`](https://nicholasclark.github.io/TIRTLtools/reference/submat.md)
  : Substitution penalty matrix for TCRdist amino acids and V-segments
- [`vdj_db`](https://nicholasclark.github.io/TIRTLtools/reference/vdj_db.md)
  : Table of known TCRs from VDJ-db

## Other
