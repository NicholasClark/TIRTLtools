# TIRTLtools R package

***Package is in active development and may change frequently.***

## Overview

`TIRTLtools` is a suite of tools for analyzing T-cell receptor (TCR) repertoires created using **TIRTL-seq** (**T**hroughput-**I**ntensive **R**apid **T**CR **L**ibrary **seq**uencing) [(Pogorelyy and Kirk et al., bioRxiv 2024)](https://www.biorxiv.org/content/10.1101/2024.09.16.613345v2).
We provide functions for analysis of **paired TCR repertoires** as well as **single-chain bulk data**.

In addition to various **analysis and plotting functions**, we provide an **efficient batched GPU implementation of TCRdist** ([Dash et al., Nature 2017](https://doi.org/10.1038/nature22383)) that works with both NVIDIA and Apple Silicon GPUs. In testing, we were able to calculate pairwise TCRdist for a repertoire of ~1 million TCRs in a few hours using a MacBook Pro (16-core GPU, M4 Pro).

## The TIRTLseq assay

[TIRTL-seq](https://doi.org/10.1101/2024.09.16.613345) is a method to obtain paired TCR-sequencing information from millions of human T cells.
Cells are split across wells in [384-well](https://www.protocols.io/view/tirtl-seq-384-well-eq2lyx37wgx9/) or [96-well](https://www.protocols.io/view/tirtl-seq-96-well-kqdg3k6xpv25/) plates with ~2,500 to ~25,000 cells per well and alpha and beta TCR chains are sequenced simultaneously in each well.

The initial output from the TIRTLseq assay consists of two tables of read counts for single-chain TCRs (one for alpha and one for beta) for each well.
Alpha and beta chains are paired using computational algorithms based on their occurrence patterns and read counts across wells.

## Installation

Instructions to [install TIRTLtools](articles/TIRTLtools.html).


## References

For details on our pairing pipeline, see the [TIRTLseq preprint (Pogorelyy and Kirk et al.)](https://www.biorxiv.org/content/10.1101/2024.09.16.613345v2) and our [github repository](https://github.com/pogorely/TIRTL).

For details on the MAD-HYPE algorithm, see [Holec and Berleant et al., 2019](https://academic.oup.com/bioinformatics/article/35/8/1318/5095649).

For details on MiXCR, see their [website](https://mixcr.com/) and publications in [Nature Methods (2015)](https://www.nature.com/articles/nmeth.3364) and [Nature Biotechnology (2017)](https://www.nature.com/articles/nbt.3979).

For details on TCRdist, see [Dash et al., Nature 2017](https://doi.org/10.1038/nature22383)
