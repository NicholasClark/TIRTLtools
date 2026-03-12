# TIRTLtools R package

***Package is in active development and may change frequently.***

## News

Dec. 8th 2025

- *Nature Methods* highlights TIRTL-seq in an *Immunology News & Views*
  article:  
  [**TIRTL-seq: heroes with
  T-SHELL**](https://doi.org/10.1038/s41592-025-02984-w).

Nov. 24th 2025

- The TIRTL-seq assay paper has been published in *Nature Methods*:  
  [**TIRTL-seq: deep, quantitative and affordable paired TCR repertoire
  sequencing**](http://doi.org/10.1038/s41592-025-02907-9).

## Installation

Run the following code in R or Rstudio to install TIRTLtools:

    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes") ## Run if you do not have the 'remotes' package installed
    }
    remotes::install_github("NicholasClark/TIRTLtools") ## install via the GitHub Repo

Load the package in R with:

    library(TIRTLtools)

## Documentation and usage instructions

For more information, see our website:
<https://nicholasclark.github.io/TIRTLtools/>

- [Installation
  instructions](https://nicholasclark.github.io/TIRTLtools/articles/TIRTLtools.html)
- [Reference
  manual](https://nicholasclark.github.io/TIRTLtools/reference/index.html)
- [Tutorials and
  vignettes](https://nicholasclark.github.io/TIRTLtools/articles/)

## Overview

`TIRTLtools` is a software package for analyzing T-cell receptor (TCR)
repertoires created using **TIRTL-seq** (**T**hroughput-**I**ntensive
**R**apid **T**CR **L**ibrary **seq**uencing) [(Pogorelyy and Kirk et
al., Nature Methods 2025)](http://doi.org/10.1038/s41592-025-02907-9).
We provide functions for analysis of **paired TCR repertoires** as well
as **single-chain bulk TCR-sequencing data**.

In addition to various **analysis and plotting functions**, we provide
an **efficient batched GPU implementation of TCRdist** ([Dash et al.,
Nature 2017](https://doi.org/10.1038/nature22383)) that works with both
NVIDIA and Apple Silicon GPUs. In testing, we were able to calculate
pairwise TCRdist for a repertoire of ~1 million TCRs in a few hours
using a MacBook Pro (16-core GPU, M4 Pro).

## The TIRTL-seq assay

The T-cell receptor (TCR), which allows for adaptive immune response, is
made up of two separate protein chains, TCR-alpha and TCR-beta. Bulk
single-chain TCR-sequencing allows for cost-effective in-depth
TCR-repertoire profiling, but does not provide chain pairings, which are
essential for determining T-cell specificity. Single-cell TCR sequencing
can produce paired chain data, but is much more expensive and thus
limited to thousands of cells in total and cost-prohibitive for
cohort-scale studies.

[TIRTL-seq](https://doi.org/10.1101/2024.09.16.613345) is a novel assay
(published Nov. 2025) that can rapidly and economically sequence
single-chain TCRs from **millions of human T-cells** and produce
**paired TCR-sequencing** information for tens of thousands of T cell
clones per sample at a **fraction of the cost of single-cell
sequencing**. Cells are split across wells in
[384-well](https://www.protocols.io/view/tirtl-seq-384-well-eq2lyx37wgx9/)
or
[96-well](https://www.protocols.io/view/tirtl-seq-96-well-kqdg3k6xpv25/)
plates with between 2,500 and 25,000 cells per well and alpha/beta TCR
chains are sequenced simultaneously in each well. Computational
algorithms based on chain co-occurrence and correlation of chain
frequency are used to identify alpha/beta chain pairs that derive from
the same T-cell clone.

We developed a **new algorithm (T-SHELL)** allows for **effective
pairing of the most abundant clones**, which had frustrated previous
algorithms. Along with paired TCR information, single-chain read
frequency can be used to identify expanded clones in one sample as well
as to track clone expansions and contractions longitudinally across
samples.

## TIRTL-seq assay protocols

See *protocols.io* for the 384-well and 96-well TIRTL-seq protocols.

- 384-well TIRTL-seq protocol:
  <https://www.protocols.io/view/tirtl-seq-384-well-eq2lyx37wgx9/>
- 96-well TIRTL-seq protocol:
  <https://www.protocols.io/view/tirtl-seq-96-well-kqdg3k6xpv25/>

## Bug reports

For any bug reports or new feature requests, please create a GitHub
Issue:

- <https://github.com/NicholasClark/TIRTLtools/issues>

## References

For details on our pairing pipeline, see the [**TIRTL-seq paper
(Pogorelyy and Kirk et al., *Nature Methods*
2025)**](http://doi.org/10.1038/s41592-025-02907-9) and our [github
repository](https://github.com/pogorely/TIRTL).

For a summary of the TIRTL-seq method and its importance, see the
*Immunology News & Views* article in *Nature Methods* highlighting
TIRTL-seq: [**TIRTL-seq: heroes with
T-SHELL**](https://doi.org/10.1038/s41592-025-02984-w).

For details on the MAD-HYPE algorithm, see [Holec and Berleant et al.,
2019](https://academic.oup.com/bioinformatics/article/35/8/1318/5095649).

For details on MiXCR, see their [website](https://mixcr.com/) and
publications in [*Nature Methods*
(2015)](https://www.nature.com/articles/nmeth.3364) and [*Nature
Biotechnology* (2017)](https://www.nature.com/articles/nbt.3979).

For details on TCRdist, see [Dash et al., *Nature*
2017](https://doi.org/10.1038/nature22383)

## Citing TIRTL-seq

Please cite the TIRTL-seq paper:

- Pogorelyy, M.V., Kirk, A.M., Adhikari, S. *et al.* TIRTL-seq: deep,
  quantitative and affordable paired TCR repertoire sequencing. *Nat
  Methods* 23, 56–64 (2026).
  <https://doi.org/10.1038/s41592-025-02907-9>

If using the 384-well or 96-well TIRTL-seq protocol, please cite the
*protocols.io* protocols:

- Mikhail Pogorelyy, Allison Kirk, Samir Adhikari, Balaji Sundararaman,
  Paul Thomas 2025. TIRTL-seq 384-well. *protocols.io*
  <https://dx.doi.org/10.17504/protocols.io.eq2lyx37wgx9/>

- Mikhail Pogorelyy, Allison Kirk, Samir Adhikari, Balaji Sundararaman,
  Paul Thomas 2025. TIRTL-seq 96-well. *protocols.io*
  <https://dx.doi.org/10.17504/protocols.io.kqdg3k6xpv25/>
