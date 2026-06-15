# Returns all of the paired TCRs from all samples in a dataset

**\[experimental\]** The `get_all_tcrs()` function aggregates the TCRs
from all samples of a dataset and puts them into one table.

## Usage

``` r
get_all_tcrs(
  data,
  chain = c("paired", "alpha", "beta"),
  remove_duplicates = TRUE
)
```

## Arguments

- data:

  a TIRTLseqDataSet object created by
  [`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)

- chain:

  the TCR chain, "alpha", "beta", or "paired" (default is paired)

- remove_duplicates:

  only return one TCR for TCRs paired by both the T-SHELL and MAD-HYPE
  algorithms (default is TRUE).

## Value

A dataframe including all of the TCRs in a dataset.

## Details

A pair of TCRs is included twice in the TIRTLseq data if it is
recognized by both the T-SHELL and MAD-HYPE algorithms. If
remove_duplicates is TRUE (default) the function will only return one of
these pairs of TCRs.

## See also

[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)

Other qc:
[`get_pair_stats()`](https://nicholasclark.github.io/TIRTLtools/reference/get_pair_stats.md),
[`get_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/get_paired_by_read_fraction_range.md),
[`plot_n_reads()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_n_reads.md),
[`plot_num_partners()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_num_partners.md),
[`plot_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired.md),
[`plot_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_by_read_fraction_range.md),
[`plot_paired_vs_rank()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_vs_rank.md),
[`plot_pairs_with_eachother()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_pairs_with_eachother.md),
[`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md),
[`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md),
[`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md),
[`summarize_data()`](https://nicholasclark.github.io/TIRTLtools/reference/summarize_data.md)

## Examples

``` r
folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)

get_all_tcrs(ts_data, chain = "paired")
#> Key: <beta_nuc>
#>          method         va     ja            cdr3a               cdr3b       vb
#>          <char>     <char> <char>           <char>              <char>   <char>
#>      1:  tshell     TRAV21 TRAJ39      CAVNAGNMLTF        CSGTANTGELFF TRBV29-1
#>      2: madhype TRAV23/DV6 TRAJ40     CAAQSGTYKYIF         CSALAKNIQYF TRBV29-1
#>      3: madhype    TRAV9-2 TRAJ10 CALSDWGTGGGNKLTF        CSAGTSVETQYF TRBV29-1
#>      4: madhype    TRAV8-6 TRAJ23   CAVTHIT_GGKLIF        CSAGTSVETQYF TRBV29-1
#>      5: madhype     TRAV34 TRAJ54     CGPLQGAQKLVF       CASRVDGYNEQFF  TRBV5-1
#>     ---                                                                        
#> 165667: madhype     TRAV19 TRAJ58    CALSEASGSRLTF     CSAQLAA_FGDTQYF TRBV20-1
#> 165668: madhype    TRAV8-6 TRAJ48  CAVTPLDFRNEKLTF   CASSIFPG_RDPRGAVF  TRBV5-1
#> 165669: madhype   TRAV13-2 TRAJ26      CAENHGQNFVF    CASSATSR_PTLEQYF  TRBV7-3
#> 165670: madhype      TRAV6 TRAJ37    CALDTSNTGKLIF      CASSPDR_HRELFF  TRBV6-4
#> 165671:  tshell TRAV14/DV4 TRAJ24   CAMRESAYWGKLQF CASSLPI*REASGANVLTF   TRBV27
#>              jb alpha_rank beta_rank alpha_readFraction beta_readFraction
#>          <char>      <int>     <int>              <num>             <num>
#>      1: TRBJ2-2          8         1       3.812515e-04      1.379590e-03
#>      2: TRBJ2-4        214         2       1.026169e-04      8.040853e-04
#>      3: TRBJ2-5         94         3       1.406358e-04      5.757980e-04
#>      4: TRBJ2-5       2294         3       2.747121e-05      5.757980e-04
#>      5: TRBJ2-1         34         4       2.205174e-04      5.147174e-04
#>     ---                                                                  
#> 165667: TRBJ2-3       8434    492846       7.748581e-06      3.614594e-08
#> 165668: TRBJ2-2     322671    493012       2.034638e-07      3.614594e-08
#> 165669: TRBJ2-7     182443    493305       9.494979e-07      3.614594e-08
#> 165670: TRBJ2-2     206907    497592       7.629894e-07      3.614594e-08
#> 165671: TRBJ2-6     176039    497698       1.000364e-06      3.614594e-08
#>         alpha_readCount beta_readCount is_functional
#>                   <int>          <int>        <lgcl>
#>      1:           29061         145569          TRUE
#>      2:            7822          84844          TRUE
#>      3:           10720          60756          TRUE
#>      4:            2094          60756         FALSE
#>      5:           16809          54311          TRUE
#>     ---                                             
#> 165667:             457              3         FALSE
#> 165668:              12              3         FALSE
#> 165669:              56              3         FALSE
#> 165670:              45              3         FALSE
#> 165671:              59              3         FALSE
#>                                                alpha_nuc
#>                                                   <char>
#>      1:                TGTGCTGTGAACGCAGGCAACATGCTCACCTTT
#>      2:             TGTGCAGCTCAATCAGGAACCTACAAATACATCTTT
#>      3: TGTGCTCTGAGTGATTGGGGGACGGGAGGAGGAAACAAACTCACCTTT
#>      4:        TGTGCTGTGACCCATATAACCAGGGAGGAAAGCTTATCTTC
#>      5:             TGTGGACCTCTTCAGGGAGCCCAGAAGCTGGTATTT
#>     ---                                                 
#> 165667:          TGTGCTCTGAGTGAGGCCAGTGGCTCTAGGTTGACCTTT
#> 165668:    TGTGCTGTGACTCCCCTCGACTTTAGAAATGAGAAATTAACCTTT
#> 165669:                TGTGCAGAGAACCATGGTCAGAATTTTGTCTTT
#> 165670:          TGTGCTCTAGACACAAGCAACACAGGCAAACTAATCTTT
#> 165671:       TGTGCAATGAGAGAGTCCGCTTACTGGGGGAAATTGCAGTTT
#>                                                          beta_nuc
#>                                                            <char>
#>      1:                      TGCAGCGGGACAGCGAACACCGGGGAGCTGTTTTTT
#>      2:                         TGCAGCGCATTAGCCAAAAACATTCAGTACTTC
#>      3:                      TGCAGCGCCGGGACTAGCGTGGAGACCCAGTACTTC
#>      4:                      TGCAGCGCCGGGACTAGCGTGGAGACCCAGTACTTC
#>      5:                   TGCGCCAGCAGAGTCGATGGCTACAATGAGCAGTTCTTC
#>     ---                                                          
#> 165667:              TGCAGTGCTCAGCTAGCGGCTTTTTCGGGGATACGCAGTATTTT
#> 165668:         TGCGCCAGCAGCATCTTTCCGGGACAGAGATCCCCGGGGAGCTGTTTTT
#> 165669:           TGCGCCAGCAGCGCTACGTCCAGACACCTACTTTAGAGCAGTACTTC
#> 165670:                 TGTGCCAGCAGTCCCGACAGGGACACCGGGAGCTGTTTTTT
#> 165671: TGTGCCAGCAGTCTACCAATCTAGCGGGAGGCATCTGGGGCCAACGTCCTGACTTTC
#>                                                                                                   alpha_beta
#>                                                                                                       <char>
#>      1:                               TGTGCTGTGAACGCAGGCAACATGCTCACCTTT_TGCAGCGGGACAGCGAACACCGGGGAGCTGTTTTTT
#>      2:                               TGTGCAGCTCAATCAGGAACCTACAAATACATCTTT_TGCAGCGCATTAGCCAAAAACATTCAGTACTTC
#>      3:                TGTGCTCTGAGTGATTGGGGGACGGGAGGAGGAAACAAACTCACCTTT_TGCAGCGCCGGGACTAGCGTGGAGACCCAGTACTTC
#>      4:                       TGTGCTGTGACCCATATAACCAGGGAGGAAAGCTTATCTTC_TGCAGCGCCGGGACTAGCGTGGAGACCCAGTACTTC
#>      5:                         TGTGGACCTCTTCAGGGAGCCCAGAAGCTGGTATTT_TGCGCCAGCAGAGTCGATGGCTACAATGAGCAGTTCTTC
#>     ---                                                                                                     
#> 165667:                 TGTGCTCTGAGTGAGGCCAGTGGCTCTAGGTTGACCTTT_TGCAGTGCTCAGCTAGCGGCTTTTTCGGGGATACGCAGTATTTT
#> 165668:      TGTGCTGTGACTCCCCTCGACTTTAGAAATGAGAAATTAACCTTT_TGCGCCAGCAGCATCTTTCCGGGACAGAGATCCCCGGGGAGCTGTTTTT
#> 165669:                    TGTGCAGAGAACCATGGTCAGAATTTTGTCTTT_TGCGCCAGCAGCGCTACGTCCAGACACCTACTTTAGAGCAGTACTTC
#> 165670:                    TGTGCTCTAGACACAAGCAACACAGGCAAACTAATCTTT_TGTGCCAGCAGTCCCGACAGGGACACCGGGAGCTGTTTTTT
#> 165671: TGTGCAATGAGAGAGTCCGCTTACTGGGGGAAATTGCAGTTT_TGTGCCAGCAGTCTACCAATCTAGCGGGAGGCATCTGGGGCCAACGTCCTGACTTTC
#>            wi    wj   wij    wa    wb       score         r       ts
#>         <int> <int> <int> <int> <int>       <num>     <num>    <num>
#>      1:     0    24   166   166   190 -5.46801297 0.9433375 39.08193
#>      2:     1    58   110   111   168  2.22849497        NA       NA
#>      3:     1    65    76    77   141  6.17941220        NA       NA
#>      4:     3    73    68    71   141  2.22543597        NA       NA
#>      5:     4    16   151   155   167  7.40350947        NA       NA
#>     ---                                                             
#> 165667:     0     0     3     3     3  0.60578435        NA       NA
#> 165668:     0     0     3     3     3  0.60578435        NA       NA
#> 165669:     0     0     3     3     3  0.60578435        NA       NA
#> 165670:     0     0     3     3     3  0.60578435        NA       NA
#> 165671:     1     0     3     4     3  0.06720111 0.9371988 37.03713
#>                 pval     pval_adj loss_a_frac loss_b_frac alpha_readCount_max
#>                <num>        <num>       <num>       <num>               <int>
#>      1: 1.747597e-92 6.346881e-84  0.12631579 0.000000000                 467
#>      2:           NA           NA  0.34319527 0.005917160                 318
#>      3:           NA           NA  0.45774648 0.007042254                 399
#>      4:           NA           NA  0.50694444 0.020833333                 172
#>      5:           NA           NA  0.09356725 0.023391813                 380
#>     ---                                                                      
#> 165667:           NA           NA  0.00000000 0.000000000                 214
#> 165668:           NA           NA  0.00000000 0.000000000                   6
#> 165669:           NA           NA  0.00000000 0.000000000                  39
#> 165670:           NA           NA  0.00000000 0.000000000                  29
#> 165671: 7.557454e-89 1.735806e-53  0.00000000 0.250000000                  21
#>         alpha_readCount_median    alpha_sem alpha_max_wells beta_readCount_max
#>                          <num>        <num>           <int>              <int>
#>      1:                  167.5 2.012317e-05             191               2501
#>      2:                   55.0 1.011121e-05             191               2483
#>      3:                  117.0 1.598201e-05             191               2334
#>      4:                   24.0 3.776374e-06             191               2334
#>      5:                   78.0 1.505977e-05             191               1270
#>     ---                                                                       
#> 165667:                  135.0 3.696180e-06             192                  1
#> 165668:                    5.0 1.059311e-07             192                  1
#> 165669:                   10.0 9.904153e-07             192                  1
#> 165670:                   11.0 5.662021e-07             192                  1
#> 165671:                   13.0 5.283245e-07             192                  1
#>         beta_readCount_median     beta_sem beta_max_wells alpha_has_stop_codon
#>                         <num>        <num>          <int>               <lgcl>
#>      1:                 739.5 7.132553e-05            191                FALSE
#>      2:                 391.5 6.976097e-05            191                FALSE
#>      3:                 150.0 6.087061e-05            191                FALSE
#>      4:                 150.0 6.087061e-05            191                FALSE
#>      5:                 278.0 3.137988e-05            191                FALSE
#>     ---                                                                       
#> 165667:                   1.0 1.727190e-08            192                FALSE
#> 165668:                   1.0 1.698039e-08            192                FALSE
#> 165669:                   1.0 2.308962e-08            192                FALSE
#> 165670:                   1.0 2.047129e-08            192                FALSE
#> 165671:                   1.0 2.384106e-08            192                FALSE
#>         alpha_has_frameshift beta_has_stop_codon beta_has_frameshift
#>                       <lgcl>              <lgcl>              <lgcl>
#>      1:                FALSE               FALSE               FALSE
#>      2:                FALSE               FALSE               FALSE
#>      3:                FALSE               FALSE               FALSE
#>      4:                 TRUE               FALSE               FALSE
#>      5:                FALSE               FALSE               FALSE
#>     ---                                                             
#> 165667:                FALSE               FALSE                TRUE
#> 165668:                FALSE               FALSE                TRUE
#> 165669:                FALSE               FALSE                TRUE
#> 165670:                FALSE               FALSE                TRUE
#> 165671:                FALSE                TRUE               FALSE
#>         alpha_is_functional beta_is_functional     alpha_id     beta_id
#>                      <lgcl>             <lgcl>       <char>      <char>
#>      1:                TRUE               TRUE      alpha_8      beta_1
#>      2:                TRUE               TRUE    alpha_214      beta_2
#>      3:                TRUE               TRUE     alpha_94      beta_3
#>      4:               FALSE               TRUE   alpha_2294      beta_3
#>      5:                TRUE               TRUE     alpha_34      beta_4
#>     ---                                                                
#> 165667:                TRUE              FALSE   alpha_8434 beta_492846
#> 165668:                TRUE              FALSE alpha_322671 beta_493012
#> 165669:                TRUE              FALSE alpha_182443 beta_493305
#> 165670:                TRUE              FALSE alpha_206907 beta_497592
#> 165671:                TRUE              FALSE alpha_176039 beta_497698
#>         receptor_rank     receptor_id  sample_id marker timepoint version
#>                 <int>          <char>     <char> <char>    <char>  <char>
#>      1:             1      receptor_1 cd4_tp1_v2    cd4       tp1      v2
#>      2:             2      receptor_2 cd4_tp1_v2    cd4       tp1      v2
#>      3:             3      receptor_3 cd4_tp1_v2    cd4       tp1      v2
#>      4:             3      receptor_3 cd4_tp1_v2    cd4       tp1      v2
#>      5:             4      receptor_4 cd4_tp1_v2    cd4       tp1      v2
#>     ---                                                                  
#> 165667:        492846 receptor_492846 cd8_tp3_v2    cd8       tp3      v2
#> 165668:        493012 receptor_493012 cd8_tp3_v2    cd8       tp3      v2
#> 165669:        493305 receptor_493305 cd8_tp3_v2    cd8       tp3      v2
#> 165670:        497592 receptor_497592 cd8_tp3_v2    cd8       tp3      v2
#> 165671:        497698 receptor_497698 cd8_tp3_v2    cd8       tp3      v2
#>                                              label
#>                                             <char>
#>      1: marker: cd4 | timepoint: tp1 | version: v2
#>      2: marker: cd4 | timepoint: tp1 | version: v2
#>      3: marker: cd4 | timepoint: tp1 | version: v2
#>      4: marker: cd4 | timepoint: tp1 | version: v2
#>      5: marker: cd4 | timepoint: tp1 | version: v2
#>     ---                                           
#> 165667: marker: cd8 | timepoint: tp3 | version: v2
#> 165668: marker: cd8 | timepoint: tp3 | version: v2
#> 165669: marker: cd8 | timepoint: tp3 | version: v2
#> 165670: marker: cd8 | timepoint: tp3 | version: v2
#> 165671: marker: cd8 | timepoint: tp3 | version: v2

```
