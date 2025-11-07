# Package index

## Bulk TIRTLseq data

### Data loading and wrangling

- [`add_metadata()`](https://nicholasclark.github.io/TIRTLtools/reference/add_metadata.md)
  **\[experimental\]** : Add metadata to a TIRTLseqData object
- [`filter_dataset()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_dataset.md)
  **\[experimental\]** : Subset a TIRTLseqData object
- [`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)
  **\[experimental\]** : Load data from TIRTLseq experiments
- [`reorder_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/reorder_samples.md)
  **\[experimental\]** : Change the order of samples in a TIRTLseqData
  object

### Data processing and cleaning

- [`TIRTL_process()`](https://nicholasclark.github.io/TIRTLtools/reference/TIRTL_process.md)
  **\[experimental\]** : Run data processing functions on a TIRTLseqData
  object
- [`add_single_chain_data()`](https://nicholasclark.github.io/TIRTLtools/reference/add_single_chain_data.md)
  **\[experimental\]** : Add single-chain read counts/fractions to the
  paired TCR data
- [`clean_pairs()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_pairs.md)
  **\[experimental\]** : Remove excess pairs for individual single
  chains
- [`filter_mait()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_mait.md)
  **\[experimental\]** : Remove MAIT (Mucosal-Associated Invariant T
  cells) TCRs
- [`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md)
  **\[experimental\]** : Identify TCRs that contain non-functional CDR3
  sequences
- [`identify_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_paired.md)
  **\[experimental\]** : Identify which single chains were paired
- [`prep_for_tcrdist()`](https://nicholasclark.github.io/TIRTLtools/reference/prep_for_tcrdist.md)
  **\[experimental\]** : Prepare paired TCRs for TCRdist calculation
- [`remove_duplicates()`](https://nicholasclark.github.io/TIRTLtools/reference/remove_duplicates.md)
  **\[experimental\]** : Removes duplicate paired TCRs

### TCR repertoire analysis

- [`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
  **\[experimental\]** : GPU implementation of TCRdist, a
  distance/similarity metric for pairs of TCRs

- [`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md)
  **\[experimental\]** : Parallelized C++ implementation of TCRdist (no
  GPU required)

- [`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md)
  **\[experimental\]** : Cluster TCRs (using the Leiden algorithm) based
  on their pairwise TCRdist values

- [`get_all_div_metrics()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_div_metrics.md)
  **\[experimental\]** :

  Returns all diversity metric options for
  [`calculate_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/calculate_diversity.md)

- [`get_all_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_tcrs.md)
  **\[experimental\]** : Returns all of the paired TCRs from all samples
  in a dataset

- [`get_pair_stats()`](https://nicholasclark.github.io/TIRTLtools/reference/get_pair_stats.md)
  **\[experimental\]** : Count the number of single-chains that were
  paired by each algorithm

- [`get_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/get_paired_by_read_fraction_range.md)
  **\[experimental\]** : Calculate the number and fraction of single
  chains that were paired within different read fraction ranges

- [`summarize_data()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_data.md)
  **\[experimental\]** : Create a summary table with number of reads and
  unique alpha/beta chains observed for each sample

### Diversity

- [`calculate_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/calculate_diversity.md)
  [`diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/calculate_diversity.md)
  **\[experimental\]** : Calculate TCR repertoire diversity metrics
- [`plot_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_diversity.md)
  **\[experimental\]** : Bar plot of clonal diversity metrics

### Plotting

- [`plot_clone_size_across_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clone_size_across_samples.md)
  **\[experimental\]** : Line plot of clone read fraction across
  multiple samples
- [`plot_clonotype_indices()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clonotype_indices.md)
  **\[experimental\]** : Stacked bar chart with fractions of reads
  attributed to the most frequent clonotypes
- [`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md)
  **\[experimental\]** : Plot clusters of similar TCRs
- [`plot_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_diversity.md)
  **\[experimental\]** : Bar plot of clonal diversity metrics
- [`plot_n_reads()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_n_reads.md)
  **\[experimental\]** : Bar plot of the number of single-chain reads
  for each sample
- [`plot_num_partners()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_num_partners.md)
  **\[experimental\]** : Stacked bar plot of the fraction of
  single-chains with different numbers of partners for each sample
- [`plot_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired.md)
  **\[experimental\]** : Stacked bar plot of the number alpha/beta
  chains paired by each pairing algorithm
- [`plot_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_by_read_fraction_range.md)
  **\[experimental\]** : Bar plot of the fraction of single chains that
  were paired within different read fraction ranges for each sample.
- [`plot_paired_vs_rank()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_vs_rank.md)
  **\[experimental\]** : A stepped plot of the cumulative number of
  paired (or unpaired) single-chains for the N most frequent
  single-chains
- [`plot_pairs_with_eachother()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_pairs_with_eachother.md)
  **\[experimental\]** : Point plot of read fraction vs. rank for the N
  most frequent alpha (left, mirrored) and beta (right) chains with
  lines between alpha and beta chains indicating a pair and a cross
  indicating an unpaired single-chain
- [`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md)
  [`rank_plot()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md)
  **\[experimental\]** : Line plot of clonotype rank vs. read fraction
  for each sample
- [`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md)
  **\[experimental\]** : Alternative point plot of read fraction vs.
  rank for the N most frequent alpha and beta chains, with a cross
  indicating an unpaired single-chain
- [`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md)
  [`sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md)
  **\[experimental\]** : Plot the overlap/agreement between samples (in
  terms of most frequent clones)
- [`plot_sample_vs_sample()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_vs_sample.md)
  **\[experimental\]** : Scatterplot of TCR clone read fraction of
  clones between two samples

### Export for external tools

- [`create_thimble_df()`](https://nicholasclark.github.io/TIRTLtools/reference/create_thimble_df.md)
  **\[experimental\]** : Convert paired TCRs to a "thimble" data frame
  for stitching together TCR nucleotide coding sequences with the
  Stitchr Python package
- [`write_stitchr_tsv()`](https://nicholasclark.github.io/TIRTLtools/reference/write_stitchr_tsv.md)
  **\[experimental\]** : Write a tab-separated "thimble" text file for
  use with Stitchr

### Well-level data functions

- [`filter_well_data()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_well_data.md)
  **\[experimental\]** : Remove rare clones from well-level TIRTLseq
  data
- [`get_well_subset()`](https://nicholasclark.github.io/TIRTLtools/reference/get_well_subset.md)
  **\[experimental\]** : Get well names from numerical rows and columns
- [`get_wells_from_edges()`](https://nicholasclark.github.io/TIRTLtools/reference/get_wells_from_edges.md)
  **\[experimental\]** : Get well names from top left and bottom right
  wells
- [`load_well_data()`](https://nicholasclark.github.io/TIRTLtools/reference/load_well_data.md)
  **\[experimental\]** : Load well-level TIRTLseq data
- [`load_well_data_to_sparse_multi()`](https://nicholasclark.github.io/TIRTLtools/reference/load_well_data_to_sparse_multi.md)
  **\[experimental\]** : Load well-level TIRTLseq data to sparse
  matrices

### Download datasets

- [`download_dataset()`](https://nicholasclark.github.io/TIRTLtools/reference/download_dataset.md)
  **\[experimental\]** : Download TIRTLseq datasets

### Data

- [`params`](https://nicholasclark.github.io/TIRTLtools/reference/params.md)
  : Table of permissible amino acids and V-segments for TCRdist
- [`submat`](https://nicholasclark.github.io/TIRTLtools/reference/submat.md)
  : Substitution penalty matrix for TCRdist amino acids and V-segments

## Single-cell TIRTLseq data

- [`clean_scTIRTLseq()`](https://nicholasclark.github.io/TIRTLtools/reference/clean_scTIRTLseq.md)
  **\[experimental\]** : Imputes missing alpha and beta chains where
  possible for single-cell TIRTLseq data
- [`summarize_scTIRTLseq()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_scTIRTLseq.md)
  **\[experimental\]** : Summarize single-cell TIRTLseq data

## Other

- [`calculate_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/calculate_diversity.md)
  [`diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/calculate_diversity.md)
  **\[experimental\]** : Calculate TCR repertoire diversity metrics
