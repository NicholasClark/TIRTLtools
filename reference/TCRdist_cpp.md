# Parallelized C++ implementation of TCRdist (no GPU required)

**\[experimental\]** This is an alternative to the GPU version of
TCRdist that is still very fast for large datasets (tens of thousands of
TCRs). It is written in C++ and will run in parallel across available
CPU cores.

## Usage

``` r
TCRdist_cpp(tcr1, tcr2 = NULL)
```

## Arguments

- tcr1:

  a data frame with one TCR per row. It must have the columns "va",
  "vb", "cdr3a", and "cdr3b"

- tcr2:

  (optional) another data frame of TCRs. If supplied, TCRdist will be
  calculated for every combination of one TCR from tcr1 and one TCR from
  tcr2. Otherwise, it will calculate TCRdist for every pair of TCRs in
  tcr1.

## Value

a list with two objects (or three if tcr2 is not null):

- matrix - a matrix of TCRdist values

- tcr1 - the input matrix tcr1, after pre-processing and removing
  unacceptable TCRs

- tcr2 (if supplied) - the input matrix tcr2, after pre-processing and
  removing unacceptable TCRs

## Details

This version of TCRdist is currently less feature-rich than
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)
and returns a dense matrix as output. It does not yet allow for sparse
output or writing output directly to a file.

## See also

Other tcr_similarity:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md)

## Examples

``` r
folder = system.file("extdata/SJTRC_TIRTLseq_minimal",
  package = "TIRTLtools")
sjtrc = load_tirtlseq(folder,
  meta_columns = c("marker", "timepoint", "version"), sep = "_",
  chain = "paired", verbose = FALSE)
#> Loading files from: /Users/nclark52/git/TIRTLtools/inst/extdata/SJTRC_TIRTLseq_minimal...
#> Found 2 beta chain pseudo-bulk files.
#> Found 2 paired chain files.
#> Loaded 2 files from 2 samples.
#> 0.1 seconds
df = get_all_tcrs(sjtrc, chain="paired", remove_duplicates = TRUE)

result = TCRdist_cpp(df)
#> Removed 384 TCRs with unknown V-segments (1.2%) from a total of 32,164 TCRs.
#> Removed 12 TCRs with short CDR3 segments (0.038%) from a total of 31,780 TCRs.
#> Removed 13,324 TCRs with non-functional CDR3 amino acid sequences (42%) from a total of 31,768 TCRs.
#> Removed 812 MAIT TCRs (4.4%) from a total of 18,444 TCRs.
#> Filtered data frame contains 17,632 TCRs (55%) of original 32,164 TCRs.
#> Removed 384 TCRs with unknown V-segments (1.2%) from a total of 32,164 TCRs.
#> Removed 12 TCRs with short CDR3 segments (0.038%) from a total of 31,780 TCRs.
#> Removed 13,324 TCRs with non-functional CDR3 amino acid sequences (42%) from a total of 31,768 TCRs.
#> Removed 812 MAIT TCRs (4.4%) from a total of 18,444 TCRs.
#> Filtered data frame contains 17,632 TCRs (55%) of original 32,164 TCRs.

mat = result$matrix
node_df = result$tcr1

mat[1:5,1:5]
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   NA   NA   NA   NA   NA
#> [2,]  338   NA   NA   NA   NA
#> [3,]  287  270   NA   NA   NA
#> [4,]  374  320  293   NA   NA
#> [5,]  382  372  386  387   NA
## note: TCRdist is symmetric. Returned matrix contains only lower-triangular values

node_df %>%
  mutate(alpha_nuc = paste(substr(alpha_nuc, 0, 20), "...", sep = ""),
         beta_nuc = paste(substr(beta_nuc, 0, 20), "...", sep = "")) %>%
  data.table::as.data.table()
#>           wi    wj   wij               alpha_nuc                beta_nuc    wa
#>        <int> <int> <int>                  <char>                  <char> <int>
#>     1:     0     0     3 TGTGTGAGAGCCGGAGGCTT... TGCAGTGCTACATCTCGGAG...     3
#>     2:     0     0     3 TGTGCAATGAGAAATCGTAG... TGTGCCAGCAGTGAGACAGG...     3
#>     3:     0     0     3 TGTGCTCTGAGCCATACTGG... TGTGCCAGCAGCCATGATGG...     3
#>     4:     0     0     3 TGTGCTCTAATCCGCAATTC... TGTGCCAGCAGCTCCGGCGA...     3
#>     5:     0     0     3 TGTGCCGTGAAGGGCCTCGG... TGTGCCAGCAAACACATCCG...     3
#>    ---                                                                        
#> 17628:    74     0   117 TGCGGCACAGAAGAAGGAGC... TGCGCCAGCAGCTGGGACAA...   191
#> 17629:    75     5    48 TGTGCAGCAAGTGAGACCTC... TGTGCCAGCAGCTTTACAGG...   123
#> 17630:    76     0    41 TGTGCTGTGAGTGAGGGCTA... TGCGCCAGCAGCTTGGGACA...   117
#> 17631:    95     0    36 TGCATCCTGAGAGAAGGGGG... TGTGCCAGCAGCTTATGGGG...   131
#> 17632:    97     0    51 TGTGCAGAGAGGGGTGGAGG... TGCGCCAGCAGCTTGCGGGA...   148
#>           wb
#>        <int>
#>     1:     3
#>     2:     3
#>     3:     3
#>     4:     3
#>     5:     3
#>    ---      
#> 17628:   117
#> 17629:    53
#> 17630:    41
#> 17631:    36
#> 17632:    51
#>                                                                                                     alpha_beta
#>                                                                                                         <char>
#>     1:                            TGTGTGAGAGCCGGAGGCTTCAAAACTATCTTT_TGCAGTGCTACATCTCGGAGAGAGCCCTACGAGCAGTACTTC
#>     2:                   TGTGCAATGAGAAATCGTAGTTCCGGGTATGCACTCAACTTC_TGTGCCAGCAGTGAGACAGGGCTCTCCTACGAGCAGTACTTC
#>     3:                         TGTGCTCTGAGCCATACTGGAGGCTTCAAAACTATCTTT_TGTGCCAGCAGCCATGATGGCTCCTACGAGCAGTACTTC
#>     4:             TGTGCTCTAATCCGCAATTCAGGAAACACACCTCTTGTCTTT_TGTGCCAGCAGCTCCGGCGAAAGAAGGTCCTACAATGAGCAGTTCTTC
#>     5: TGTGCCGTGAAGGGCCTCGGAGGAAGCCAAGGAAATCTCATCTTT_TGTGCCAGCAAACACATCCGGGACAGGGGGGGGTGGAGCAATCAGCCCCAGCATTTT
#>    ---                                                                                                        
#> 17628:                      TGCGGCACAGAAGAAGGAGCTAGCAACTATAAACTGACATTT_TGCGCCAGCAGCTGGGACAAACCTTATGGCTACACCTTC
#> 17629:                         TGTGCAGCAAGTGAGACCTCCTACGACAAGGTGATATTT_TGTGCCAGCAGCTTTACAGGGGAAAACACCATATATTTT
#> 17630:                         TGTGCTGTGAGTGAGGGCTACAAGCTCAGCTTT_TGCGCCAGCAGCTTGGGACAGGGAGGCAACCAGCCCCAGCATTTT
#> 17631:                      TGCATCCTGAGAGAAGGGGGGAGTGACATGCGCTTT_TGTGCCAGCAGCTTATGGGGGGGGGGGGACACTGAAGCTTTCTTT
#> 17632:                      TGTGCAGAGAGGGGTGGAGGCTTCAAAACTATCTTT_TGCGCCAGCAGCTTGCGGGACCAGCCCTACAACGAGCAGTTCTTC
#>         method         r       ts         pval     pval_adj loss_a_frac
#>         <char>     <num>    <num>        <num>        <num>       <num>
#>     1: madhype        NA       NA           NA           NA   0.0000000
#>     2: madhype        NA       NA           NA           NA   0.0000000
#>     3: madhype        NA       NA           NA           NA   0.0000000
#>     4: madhype        NA       NA           NA           NA   0.0000000
#>     5: madhype        NA       NA           NA           NA   0.0000000
#>    ---                                                                 
#> 17628:  tshell 0.5679823 9.487338 1.044199e-17 6.459864e-11   0.0000000
#> 17629: madhype        NA       NA           NA           NA   0.0390625
#> 17630: madhype        NA       NA           NA           NA   0.0000000
#> 17631: madhype        NA       NA           NA           NA   0.0000000
#> 17632: madhype        NA       NA           NA           NA   0.0000000
#>        loss_b_frac    score           cdr3a            va     ja
#>              <num>    <num>          <char>        <char> <char>
#>     1:   0.0000000 0.618488     CVRAGGFKTIF   TRAV12-1*01  TRAJ9
#>     2:   0.0000000 0.618488  CAMRNRSSGYALNF TRAV14/DV4*01 TRAJ41
#>     3:   0.0000000 0.618488   CALSHTGGFKTIF    TRAV9-2*01  TRAJ9
#>     4:   0.0000000 0.618488  CALIRNSGNTPLVF     TRAV16*01 TRAJ29
#>     5:   0.0000000 0.618488 CAVKGLGGSQGNLIF   TRAV12-2*01 TRAJ42
#>    ---                                                          
#> 17628:   0.3874346     -Inf  CGTEEGASNYKLTF     TRAV30*01 TRAJ53
#> 17629:   0.5859375 0.605473   CAASETSYDKVIF   TRAV13-1*01 TRAJ50
#> 17630:   0.6495726 4.604051     CAVSEGYKLSF    TRAV8-6*01 TRAJ20
#> 17631:   0.7251908 1.073628    CILREGGSDMRF   TRAV26-2*01 TRAJ43
#> 17632:   0.6554054 1.113315    CAERGGGFKTIF      TRAV5*01  TRAJ9
#>                      cdr3b          vb      jb  sample_id marker timepoint
#>                     <char>      <char>  <char>     <char> <char>    <char>
#>     1:      CSATSRREPYEQYF TRBV20-1*01 TRBJ2-7 cd8_tp1_v2    cd8       tp1
#>     2:      CASSETGLSYEQYF TRBV25-1*01 TRBJ2-7 cd8_tp1_v2    cd8       tp1
#>     3:       CASSHDGSYEQYF  TRBV3-1*01 TRBJ2-7 cd8_tp1_v2    cd8       tp1
#>     4:    CASSSGERRSYNEQFF  TRBV5-4*01 TRBJ2-1 cd8_tp1_v2    cd8       tp1
#>     5: CASKHIRDRGGWSNQPQHF   TRBV28*01 TRBJ1-5 cd8_tp1_v2    cd8       tp1
#>    ---                                                                    
#> 17628:       CASSWDKPYGYTF  TRBV5-1*01 TRBJ1-2 cd8_tp2_v2    cd8       tp2
#> 17629:       CASSFTGENTIYF  TRBV7-9*01 TRBJ1-3 cd8_tp2_v2    cd8       tp2
#> 17630:     CASSLGQGGNQPQHF  TRBV5-1*01 TRBJ1-5 cd8_tp2_v2    cd8       tp2
#> 17631:     CASSLWGGGDTEAFF  TRBV7-2*01 TRBJ1-1 cd8_tp2_v2    cd8       tp2
#> 17632:     CASSLRDQPYNEQFF  TRBV5-1*01 TRBJ2-1 cd8_tp2_v2    cd8       tp2
#>        version                                      label    va_orig  vb_orig
#>         <char>                                     <char>     <char>   <char>
#>     1:      v2 marker: cd8 | timepoint: tp1 | version: v2   TRAV12-1 TRBV20-1
#>     2:      v2 marker: cd8 | timepoint: tp1 | version: v2 TRAV14/DV4 TRBV25-1
#>     3:      v2 marker: cd8 | timepoint: tp1 | version: v2    TRAV9-2  TRBV3-1
#>     4:      v2 marker: cd8 | timepoint: tp1 | version: v2     TRAV16  TRBV5-4
#>     5:      v2 marker: cd8 | timepoint: tp1 | version: v2   TRAV12-2   TRBV28
#>    ---                                                                       
#> 17628:      v2 marker: cd8 | timepoint: tp2 | version: v2     TRAV30  TRBV5-1
#> 17629:      v2 marker: cd8 | timepoint: tp2 | version: v2   TRAV13-1  TRBV7-9
#> 17630:      v2 marker: cd8 | timepoint: tp2 | version: v2    TRAV8-6  TRBV5-1
#> 17631:      v2 marker: cd8 | timepoint: tp2 | version: v2   TRAV26-2  TRBV7-2
#> 17632:      v2 marker: cd8 | timepoint: tp2 | version: v2      TRAV5  TRBV5-1
#>        is_functional
#>               <lgcl>
#>     1:          TRUE
#>     2:          TRUE
#>     3:          TRUE
#>     4:          TRUE
#>     5:          TRUE
#>    ---              
#> 17628:          TRUE
#> 17629:          TRUE
#> 17630:          TRUE
#> 17631:          TRUE
#> 17632:          TRUE

```
