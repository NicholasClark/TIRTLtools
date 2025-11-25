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

  a TIRTLseqData object created by
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
#> Loading files from: /Users/nclark52/git/TIRTLtools/inst/extdata/SJTRC_TIRTL_seq_longitudinal...
#> Found 6 beta chain pseudo-bulk files.
#> Found 6 paired chain files.
#> Loaded 18 files from 6 samples.
#> 15.3 seconds

get_all_tcrs(ts_data, chain = "paired")
#>            wi    wj   wij                                        alpha_nuc
#>         <int> <int> <int>                                           <char>
#>      1:     0     0     3          TGTGCCGTGAACGGAGGTAACGACTACAAGCTCAGCTTT
#>      2:     0     0     3       TGTGCTCTTCAAACTTCTGGTGGCTACAATAAGCTGATTTTT
#>      3:     0     0     3             TGTGCTCCCCGGTACACCGGTAACCAGTTCTATTTT
#>      4:     0     0     3 TGTGCAATGAGAGAGGTCGCTTCTGGTGGCTACAATAAGCTGATTTTT
#>      5:     0     0     3             TGTGCTGCCAGCGGGGGCTTTCAGAAACTTGTATTT
#>     ---                                                                   
#> 165667:    94     0    37    TGTGCCGTGGAGGGGGAGTCTGGTGGCTACAATAAGCTGATTTTT
#> 165668:    96     0    50             TGTGCCCCCATGGATAGCAACTATCAGTTAATCTGG
#> 165669:    97     0    36    TGCGCTGTGAGAGATAGGGACACCAATGCAGGCAAATCAACCTTT
#> 165670:   101     0    36                TGTGCAGACGGCCCCCAAGGAAATCTCATCTTT
#> 165671:   103     0    46       TGCGGCACAGAGAGTGGAGGTAGCAACTATAAACTGACATTT
#>                                                 beta_nuc    wa    wb
#>                                                   <char> <int> <int>
#>      1:          TGTGCCAGCTCACCAGAGCAGGGTCTGGCCCAGCATTTT     3     3
#>      2:                TGTGCCAGCAGCTATCAGGGCCAGCAGTTCTTC     3     3
#>      3:       TGTGCCATCAGTCGGGGACAGGGGGCAGGTGAGCAGTTCTTC     3     3
#>      4: TGTGCCAGCTCGACAGGGGGACAGCGACCTAATGAAAAACTGTTTTTT     3     3
#>      5:    TGCGCCAGCAGCCCAACAGGGGGGCAGGACTATGGCTACACCTTC     3     3
#>     ---                                                             
#> 165667:   TGTGCCAGCAGCGCCCACGGGGAAGGGAACGAAAAACTGTTTTTTT   131    37
#> 165668:      TGCAGTGCTAAATTAGCGGGCGACACCGGGGAGCTGTTTTTTT   146    50
#> 165669:      TGTGCCAGCAGTTTTTGCGGGGGCGGGCTACGAGCAGTACTTC   133    36
#> 165670:      TGTGCCAGCAGTTTTTGCGGGGGCGGGCTACGAGCAGTACTTC   137    36
#> 165671:          TGTGCCAGCAGCTGGGACAAAACCTACGAGCAGTACTTC   149    46
#>                                                                                                alpha_beta
#>                                                                                                    <char>
#>      1:                   TGTGCCGTGAACGGAGGTAACGACTACAAGCTCAGCTTT_TGTGCCAGCTCACCAGAGCAGGGTCTGGCCCAGCATTTT
#>      2:                      TGTGCTCTTCAAACTTCTGGTGGCTACAATAAGCTGATTTTT_TGTGCCAGCAGCTATCAGGGCCAGCAGTTCTTC
#>      3:                   TGTGCTCCCCGGTACACCGGTAACCAGTTCTATTTT_TGTGCCATCAGTCGGGGACAGGGGGCAGGTGAGCAGTTCTTC
#>      4: TGTGCAATGAGAGAGGTCGCTTCTGGTGGCTACAATAAGCTGATTTTT_TGTGCCAGCTCGACAGGGGGACAGCGACCTAATGAAAAACTGTTTTTT
#>      5:                TGTGCTGCCAGCGGGGGCTTTCAGAAACTTGTATTT_TGCGCCAGCAGCCCAACAGGGGGGCAGGACTATGGCTACACCTTC
#>     ---                                                                                                  
#> 165667:      TGTGCCGTGGAGGGGGAGTCTGGTGGCTACAATAAGCTGATTTTT_TGTGCCAGCAGCGCCCACGGGGAAGGGAACGAAAAACTGTTTTTTT
#> 165668:                  TGTGCCCCCATGGATAGCAACTATCAGTTAATCTGG_TGCAGTGCTAAATTAGCGGGCGACACCGGGGAGCTGTTTTTTT
#> 165669:         TGCGCTGTGAGAGATAGGGACACCAATGCAGGCAAATCAACCTTT_TGTGCCAGCAGTTTTTGCGGGGGCGGGCTACGAGCAGTACTTC
#> 165670:                     TGTGCAGACGGCCCCCAAGGAAATCTCATCTTT_TGTGCCAGCAGTTTTTGCGGGGGCGGGCTACGAGCAGTACTTC
#> 165671:                TGCGGCACAGAGAGTGGAGGTAGCAACTATAAACTGACATTT_TGTGCCAGCAGCTGGGACAAAACCTACGAGCAGTACTTC
#>          method     r    ts  pval pval_adj loss_a_frac loss_b_frac     score
#>          <char> <num> <num> <num>    <num>       <num>       <num>     <num>
#>      1: madhype    NA    NA    NA       NA           0   0.0000000 0.2011610
#>      2: madhype    NA    NA    NA       NA           0   0.0000000 0.2011610
#>      3: madhype    NA    NA    NA       NA           0   0.0000000 0.2011610
#>      4: madhype    NA    NA    NA       NA           0   0.0000000 0.2011610
#>      5: madhype    NA    NA    NA       NA           0   0.0000000 0.2011610
#>     ---                                                                     
#> 165667: madhype    NA    NA    NA       NA           0   0.7175573 1.5618088
#> 165668: madhype    NA    NA    NA       NA           0   0.6575342 1.6288014
#> 165669: madhype    NA    NA    NA       NA           0   0.7293233 1.0688280
#> 165670: madhype    NA    NA    NA       NA           0   0.7372263 0.5297975
#> 165671: madhype    NA    NA    NA       NA           0   0.6912752 0.4727439
#>                    cdr3a         va     ja            cdr3b       vb      jb
#>                   <char>     <char> <char>           <char>   <char>  <char>
#>      1:    CAVNGGNDYKLSF   TRAV12-2 TRAJ20    CASSPEQGLAQHF   TRBV18 TRBJ1-5
#>      2:   CALQTSGGYNKLIF    TRAV9-2  TRAJ4      CASSYQGQQFF  TRBV7-6 TRBJ2-1
#>      3:     CAPRYTGNQFYF     TRAV41 TRAJ49   CAISRGQGAGEQFF TRBV10-3 TRBJ2-1
#>      4: CAMREVASGGYNKLIF TRAV14/DV4  TRAJ4 CASSTGGQRPNEKLFF  TRBV6-1 TRBJ1-4
#>      5:     CAASGGFQKLVF      TRAV2  TRAJ8  CASSPTGGQDYGYTF  TRBV5-1 TRBJ1-2
#>     ---                                                                     
#> 165667:  CAVEGESGGYNKLIF   TRAV12-2  TRAJ4 CASSAHGE_GTKNCFF    TRBV9 TRBJ1-4
#> 165668:     CAPMDSNYQLIW    TRAV1-2 TRAJ33  CSAKLAG_TPGSCFF TRBV20-1 TRBJ2-2
#> 165669:  CAVRDRDTNAGKSTF    TRAV1-1 TRAJ27  CASSFCG_AGYEQYF   TRBV28 TRBJ2-7
#> 165670:      CADGPQGNLIF   TRAV13-2 TRAJ42  CASSFCG_AGYEQYF   TRBV28 TRBJ2-7
#> 165671:   CGTESGGSNYKLTF     TRAV30 TRAJ53    CASSWDKTYEQYF  TRBV5-1 TRBJ2-7
#>          sample_id marker timepoint version
#>             <char> <char>    <char>  <char>
#>      1: cd4_tp1_v2    cd4       tp1      v2
#>      2: cd4_tp1_v2    cd4       tp1      v2
#>      3: cd4_tp1_v2    cd4       tp1      v2
#>      4: cd4_tp1_v2    cd4       tp1      v2
#>      5: cd4_tp1_v2    cd4       tp1      v2
#>     ---                                    
#> 165667: cd8_tp3_v2    cd8       tp3      v2
#> 165668: cd8_tp3_v2    cd8       tp3      v2
#> 165669: cd8_tp3_v2    cd8       tp3      v2
#> 165670: cd8_tp3_v2    cd8       tp3      v2
#> 165671: cd8_tp3_v2    cd8       tp3      v2
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
