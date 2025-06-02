# TIRTLtools R package

***Package is in active development and may change frequently.***

## Overview

`TIRTLtools` is a suite of tools for analyzing T-cell receptor (TCR) repertoires created using **TIRTL-seq** (**T**hroughput-**I**ntensive **R**apid **T**CR **L**ibrary **seq**uencing) [(Pogorelyy and Kirk et al., bioRxiv 2024)](https://www.biorxiv.org/content/10.1101/2024.09.16.613345v2).
We provide functions for analysis of **paired TCR repertoires** as well as **single-chain bulk data**.

In addition to various **analysis and plotting functions**, we provide an **efficient batched GPU implementation of TCRdist** ([Dash et al., Nature 2017](https://doi.org/10.1038/nature22383)) that works with both NVIDIA and Apple Silicon GPUs. In testing, we were able to calculate pairwise TCRdist for a repertoire of ~1 million TCRs in a few hours using a MacBook Pro (16-core GPU, M4 Pro).

### Input data for package

The package requires paired-chain and single-chain "pseudo-bulk" data files created by the [TIRTLseq pipeline](https://github.com/pogorely/TIRTL).

The pipeline reads [MiXCR](https://mixcr.com/) output files from TCR repertoire sequencing for individual wells in a multiwell plate (e.g. 96-well or 384-well plate). It uses two complementary algorithms,  [MAD-HYPE](https://doi.org/10.1093/bioinformatics/bty801) and [T-SHELL](https://www.biorxiv.org/content/10.1101/2024.09.16.613345v2), to pair alpha and beta chains. It outputs three files: a .tsv file containing alpha and beta chain pairs predicted by each algorithm (file ending in "TIRTLoutput.tsv"), and two .tsv files containing aggregated single-chain "pseudo-bulk" data for alpha or beta chains from all wells on a plate (files ending in "pseudobulk_TRA.tsv" and "pseudobulk_TRB.tsv").

**MAD-HYPE** (**M**ulticell **A**nalytical **D**econvolution for **H**igh **Y**ield **P**aired-chain **E**valuation, Holec and Berleant et al., 2019) is a Bayesian method that computes the likelihood of alpha and beta chain matches based on their co-occurence patterns in wells on a plate. We re-implemented the MAD-HYPE algorithm for GPU computation, increasing efficiency and speed by two orders of magnitude.

**T-SHELL** (**T**CR ɑβ **S**equence **H**ighly **E**fficient **L**inkage **L**earning, Pogorelyy and Kirk et al., bioRxiv 2024) is an algorithm we developed that uses the correlation between TCRalpha and TCRbeta clonotype relative frequencies across wells in a plate (rather than their presence or absence) to pair them. It works well in the case of the most abundant clones, which are found in all or almost all wells, making it impossible for the MAD-HYPE algorithm to pair them using their co-occurrence patterns.

Example input data for the TIRTLseq pipeline can be found on [Zenodo](https://zenodo.org/records/14010377). For example, the aph17pre_exp_clones.tar.gz (172 MB) and aph17post_exp_clones.tar.gz (1.1 GB) files contain well-level TCR clonotype data (MiXCR output) for pre- and post-antigen-independent T cell expansion experiments.

### Package functions

#### Data loading and filtering:
  * `load_tirtlseq()` -- load paired TCR and single-chain data from TIRTL-seq experiments
  * `filter_dataset()` -- subset a TIRTLseq dataset to only the requested samples
  
#### Repertoire analysis:
  * `TCRdist()` -- calculate TCRdist via batched GPU calculation
  * `cluster_tcrs()` -- cluster TCRs (with themselves and with TCRs from VDJ-db) based on pairwise TCRdist
  * `diversity()` -- calculate repertoire diversity metrics
  
#### Plotting:
  * `sample_overlap()` -- plot overlap/agreement between samples (in terms of most frequent clones)
  * `plot_sample_vs_sample()` -- scatterplot of readFraction between two samples
  * `rank_plot()` -- a line plot of clonotype rank vs. readFraction for different samples
  * `plot_clusters()` -- plot top "n" clusters with the most observed TCRs -- currently outputs a network graph, a umap, and a heatmap
  * `plot_gene_usage()` -- bar plots of V-gene segment usage (or J-gene segment usage)
  * `plot_num_partners()` -- plot number of partners for alpha and beta chains in paired data
  * `plot_clonotype_indices()` -- stacked bar chart with fractions of reads for most frequent clonotypes.
  * `plot_diversity()` -- bar plots of different diversity metrics

#### Helper functions:
  * `add_single_chain_data()` -- add alpha- and beta-chain read counts/fractions to the paired TCR data
  * `identify_non_functional_seqs()` -- mark cdr3 sequences that have stop codons (*) or frame shifts (_)
  * `identify_paired()` -- identify which alpha and beta chains in the pseudo-bulk data are found in pairs
  * `get_all_div_metrics()` -- returns all options for diversity metrics (currently: simpson, gini, gini.simpson, inv.simpson, shannon, berger.parker, richness, d50, dXX, renyi, hill)
  * `get_all_tcrs()` -- aggregate all of the TCRs from all samples in a dataset.

## Installation

You can install the package directly from GitHub:
```R
remotes::install_github("NicholasClark/TIRTLtools")
```
Alternatively, you may clone the repository to a local folder
```
git clone https://github.com/NicholasClark/TIRTLtools.git
```
and load it using devtools:
```R
devtools::load_all("./TIRTLtools") ### location of cloned repo on your machine
```

## Quick start

### TCRdist GPU implementation

This function computes the similarity of two T-cell receptors (TCRs) using TCRdist ([Dash et al., Nature 2017](https://doi.org/10.1038/nature22383)). The function can efficiently calculate pairwise similarities of thousands or millions of TCRs using a GPU if it is available. It relies on python code using [cupy](https://cupy.dev/) (Nvidia GPUs) or [mlx](https://opensource.apple.com/projects/mlx/) (Apple Silicon GPUs), which is called through the `reticulate` R package.

The following example computes TCRdist values for all pairs of annotated TCRs in the [VDJdb](https://vdjdb.cdr3.net/) database. It outputs all values with TCRdist <= 90.

```R
#library(TIRTLtools)
devtools::load_all("./TIRTLtools") ### location of cloned repo on your machine
library(dplyr)

### Load an example dataset (annotated TCRs from VDJ-db)
tcr1 = TIRTLtools::vdj_db
### note: You may replace this with a file of your choice - needs to have columns "va", "vb", "cdr3a", and "cdr3b"

as_tibble(tcr1)
# # A tibble: 26,213 × 15
#       complex.id cdr3a         va       ja    cdr3b vb    jb    Species MHC.A MHC.B MHC.class Epitope Epitope.gene
#         <int>   <chr>         <chr>    <chr> <chr> <chr> <chr> <chr>   <chr> <chr> <chr>     <chr>   <chr>
# 1          5 CAVRDGGTGFQKLVF   TRAV3*01 TRAJ… CASR… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      HPVGEA… EBNA1
# 2          6 CAARGIGSGTYKYIF   TRAV13-… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      HPVGEA… EBNA1
# 3          8 CAVRDGTANNLFF     TRAV3*01 TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 4          9 CAFMKDAGGTSYGKLTF TRAV38-… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 5         10 CAGGGSQGNLIF      TRAV27*… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 6         11 CALAGSQGNLIF      TRAV6*01 TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 7         12 CAGADGGSQGNLIF    TRAV27*… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 8         13 CAAGGSQGNLIF      TRAV27*… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 9         14 CAAPPSNTGKLIF     TRAV3*01 TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 10         15 CAVRDLLTNAGKSTF   TRAV3*01 TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# # ℹ 26,203 more rows
# # ℹ 2 more variables: Epitope.species <chr>, Reference <chr>
# # ℹ Use `print(n = ...)` to see more rows
result = TCRdist(tcr1 = tcr1, tcrdist_cutoff = 90, chunk_size = 5000)
# Checking for available GPU...
#
# Apple Silicon GPU detected:
#   Apple Silicon GPU (M1/M2/M3)
# Checking for GPU-related Python modules...
#
# 'mlx' is installed (for Apple Silicon GPUs).
# Loading mlx to perform TCRdist
# total number of chunks (rows): 5
# total number of chunks (cols): 5
# Processing chunk (rows) 0
# Processing chunk (rows) 5000
# Processing chunk (rows) 10000
# Processing chunk (rows) 15000
# Processing chunk (rows) 20000
# Processing chunk (rows) 25000
# Time taken: 2.970693 seconds

edge_df = result[['TCRdist_df']] ### table of TCRdist values <= cutoff
meta_df = result[['tcr1']] ### table of input data with indices

as_tibble(edge_df)
# # A tibble: 101,797 × 3
#       edge1_0index edge2_0index TCRdist
#           <int>        <int>   <int>
# 1            5            4      87
# 2            7            4       3
# 3            7            5      90
# 4           12            4      15
# 5           12            5      78
# 6           12            7      12
# 7           13            4      78
# 8           13            6      72
# 9           13            7      81
# 10           14            6      45
# # ℹ 101,787 more rows
# # ℹ Use `print(n = ...)` to see more rows

as_tibble(meta_df)
# # A tibble: 26,213 × 16
#   complex.id cdr3a      va    ja    cdr3b vb    jb    Species MHC.A MHC.B MHC.class Epitope Epitope.gene
#       <int>   <chr>     <chr> <chr> <chr> <chr> <chr> <chr>   <chr> <chr> <chr>     <chr>   <chr>
# 1          5 CAVRDGGTG… TRAV… TRAJ… CASR… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      HPVGEA… EBNA1
# 2          6 CAARGIGSG… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      HPVGEA… EBNA1
# 3          8 CAVRDGTAN… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 4          9 CAFMKDAGG… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 5         10 CAGGGSQGN… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 6         11 CALAGSQGN… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 7         12 CAGADGGSQ… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 8         13 CAAGGSQGN… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 9         14 CAAPPSNTG… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# 10         15 CAVRDLLTN… TRAV… TRAJ… CASS… TRBV… TRBJ… HomoSa… HLA-… B2M   MHCI      GILGFV… M
# # ℹ 26,203 more rows
# # ℹ 3 more variables: Epitope.species <chr>, Reference <chr>, tcr_index <dbl>
# # ℹ Use `print(n = ...)` to see more rows

#' TCRdist function output:
#' A list with two items:
#'    'TCRdist_df' - A dataframe containing all TCRdist values less than the cutoff (default = 90).
#'                  Contains the indices of the two TCRs and the TCRdist value.
#'    'tcr1' - A dataframe with the input data, along with a column "tcr_index" with the index of each TCR
```

## References

For details of our pairing pipeline, see the [TIRTLseq preprint (Pogorelyy and Kirk et al.)](https://www.biorxiv.org/content/10.1101/2024.09.16.613345v2) and our [github repository](https://github.com/pogorely/TIRTL).

For details of the MAD-HYPE algorithm, see [Holec and Berleant et al., 2019](https://academic.oup.com/bioinformatics/article/35/8/1318/5095649).

For details of MiXCR, see their [website](https://mixcr.com/) and publications in [Nature Methods (2015)](https://www.nature.com/articles/nmeth.3364) and [Nature Biotechnology (2017)](https://www.nature.com/articles/nbt.3979).

For details of TCRdist, see [Dash et al., Nature 2017](https://doi.org/10.1038/nature22383)
