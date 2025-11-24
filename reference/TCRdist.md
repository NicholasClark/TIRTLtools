# GPU implementation of TCRdist, a distance/similarity metric for pairs of TCRs

**\[experimental\]** An efficient, batched version of TCRdist that is
compatible with both NVIDIA and Apple Silicon GPUs.

## Usage

``` r
TCRdist(
  tcr1,
  tcr2 = NULL,
  remove_MAIT = FALSE,
  params = NULL,
  submat = NULL,
  tcrdist_cutoff = 90,
  chunk_size = 1000,
  print_chunk_size = 10,
  print_res = TRUE,
  only_lower_tri = TRUE,
  return_data = TRUE,
  write_to_tsv = FALSE,
  backend = c("auto", "cpu", "cupy", "mlx"),
  fork = NULL,
  shared = NULL
)
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

- remove_MAIT:

  whether to remove TCRs from MAIT cells (default is FALSE)

- params:

  (optional) a table of valid parameters for amino acids and va/vb
  segments. (default is NULL, which uses TIRTLtools::params)

- submat:

  (optional) a substitution matrix with mismatch penalties for each
  combination of amino acids or va/vb segments (default is NULL, which
  uses TIRTLtools::submat).

- tcrdist_cutoff:

  (optional) discard all TCRdist values above this cutoff (default is
  90).

- chunk_size:

  (optional) The chunk size to use in calculation of TCRdist (default
  1000). If set at n, it will calculate pairwise TCRdist for n x n TCRs
  at once. This may be as high as allowable by GPU memory (in our
  testing, a chunk_size of 1000 to 5000 provided the fastest runtime and
  chunk_size of over 7500 resulted in memory errors on some GPUs).

- print_chunk_size:

  (optional) print a line of output for every n TCRs processed (default
  1000)

- print_res:

  (optional) print summary of results (default is TRUE)

- only_lower_tri:

  (optional) return one TCRdist value for each pair (like the lower
  triangle of a symmetric matrix). Default is TRUE.

- return_data:

  (optional) whether to return the output result from the function. With
  large data it may be desirable to write the result to disk instead.
  (default is TRUE, returns output)

- write_to_tsv:

  (optional) write the results to a tab-separated file ".tsv" (default
  is FALSE, does not write .tsv file)

- backend:

  (optional) the CPU or GPU backend to use (default "auto")

- fork:

  (optional) a TRUE/FALSE value for whether to "fork" a new Python
  process for running TCRdist via the "basilisk" package. Default is
  NULL, which should use choose a safe value based on how the package is
  loaded.

- shared:

  (optional) a TRUE/FALSE value for whether to "share" the Python
  process for running TCRdist via the "basilisk" package. Default is
  NULL, which should use choose a safe value based on how the package is
  loaded.

## Value

A list with entries:

`$TCRdist_df` - a data frame with three columns: "node1_0index",
"node2_0index", and "TCRdist". The first two columns contain the indices
(0-indexed) of the TCRs for each pair. The last column contains the
TCRdist if it is below the cutoff. The output is sparse in that it only
contains pairs that have TCRdist \<= cutoff.

`$tcr1` - a data frame of the TCRs supplied to the function. It contains
an additional column "tcr_index" with the (0-indexed) index of each TCR.

`$tcr2` - a similar data frame for tcr2, if it was supplied.

## Details

This function calculates pairwise TCRdist (Dash et al., Nature 2017) for
a set of TCRs (or between two sets of TCRs) and returns a sparse output
with the TCRdist and indices of all pairs that have TCRdist less than or
equal to a desired cutoff (default cutoff is 90).

The function uses the `reticulate` package to call a python script that
uses `cupy` (NVIDIA GPUs), `mlx` (Apple Silicon GPUs), or `numpy` (no
GPU) to calculate TCRdist efficiently.

## See also

[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md),
and
[`identify_non_functional_seqs()`](https://nicholasclark.github.io/TIRTLtools/reference/identify_non_functional_seqs.md)

Other tcr_similarity:
[`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md),
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
#> 1 seconds
df = get_all_tcrs(sjtrc, chain="paired", remove_duplicates = TRUE)
result = TCRdist(df, tcrdist_cutoff = 90)
#> Removed 384 TCRs with unknown V-segments (1.2%) from a total of 32,164 TCRs.
#> Removed 12 TCRs with short CDR3 segments (0.038%) from a total of 31,780 TCRs.
#> Removed 13,324 TCRs with non-functional CDR3 amino acid sequences (42%) from a total of 31,768 TCRs.
#> Filtered data frame contains 18,444 TCRs (57%) of original 32,164 TCRs.
edge_df = result[['TCRdist_df']] %>%
  data.table::as.data.table() ### table of TCRdist values <= cutoff
node_df = result[['tcr1']] %>%
  data.table::as.data.table() ### table of input data with indices

edge_df ## sparse 3-column output: node1, node2, TCRdist
#>        node1_0index node2_0index TCRdist
#>               <int>        <int>   <int>
#>     1:          158          157       0
#>     2:          286          284      84
#>     3:          287          283      89
#>     4:          288          284      89
#>     5:          289          285      89
#>    ---                                  
#> 18422:        18434        18073       0
#> 18423:        18437        18406      72
#> 18424:        18440        18241      90
#> 18425:        18440        18324      66
#> 18426:        18440        18422      84
## note that indices start at 0 and are found in node_df$tcr_index

node_df %>%
  select(tcr_index, everything()) %>%
  mutate(alpha_nuc = paste(substr(alpha_nuc, 0, 20), "...", sep = ""),
         beta_nuc = paste(substr(beta_nuc, 0, 20), "...", sep = ""))
#>        tcr_index    wi    wj   wij               alpha_nuc
#>            <num> <int> <int> <int>                  <char>
#>     1:         0     0     0     3 TGTGTGAGAGCCGGAGGCTT...
#>     2:         1     0     0     3 TGTGCAATGAGAAATCGTAG...
#>     3:         2     0     0     3 TGTGCTCTGAGCCATACTGG...
#>     4:         3     0     0     3 TGTGCTCTAATCCGCAATTC...
#>     5:         4     0     0     3 TGTGCCGTGAAGGGCCTCGG...
#>    ---                                                    
#> 18440:     18439    76     0    41 TGTGCTGTGAGTGAGGGCTA...
#> 18441:     18440    78     0    18 TGTGCTGTGTTGGATAGCAA...
#> 18442:     18441    85     1    29 TGTGCTGCTATGGATAGCAA...
#> 18443:     18442    95     0    36 TGCATCCTGAGAGAAGGGGG...
#> 18444:     18443    97     0    51 TGTGCAGAGAGGGGTGGAGG...
#>                       beta_nuc    wa    wb
#>                         <char> <int> <int>
#>     1: TGCAGTGCTACATCTCGGAG...     3     3
#>     2: TGTGCCAGCAGTGAGACAGG...     3     3
#>     3: TGTGCCAGCAGCCATGATGG...     3     3
#>     4: TGTGCCAGCAGCTCCGGCGA...     3     3
#>     5: TGTGCCAGCAAACACATCCG...     3     3
#>    ---                                    
#> 18440: TGCGCCAGCAGCTTGGGACA...   117    41
#> 18441: TGTGCCAGCAGCCAAGAAGG...    96    18
#> 18442: TGCAGTGCTAAACTAGCGGG...   114    30
#> 18443: TGTGCCAGCAGCTTATGGGG...   131    36
#> 18444: TGCGCCAGCAGCTTGCGGGA...   148    51
#>                                                                                                     alpha_beta
#>                                                                                                         <char>
#>     1:                            TGTGTGAGAGCCGGAGGCTTCAAAACTATCTTT_TGCAGTGCTACATCTCGGAGAGAGCCCTACGAGCAGTACTTC
#>     2:                   TGTGCAATGAGAAATCGTAGTTCCGGGTATGCACTCAACTTC_TGTGCCAGCAGTGAGACAGGGCTCTCCTACGAGCAGTACTTC
#>     3:                         TGTGCTCTGAGCCATACTGGAGGCTTCAAAACTATCTTT_TGTGCCAGCAGCCATGATGGCTCCTACGAGCAGTACTTC
#>     4:             TGTGCTCTAATCCGCAATTCAGGAAACACACCTCTTGTCTTT_TGTGCCAGCAGCTCCGGCGAAAGAAGGTCCTACAATGAGCAGTTCTTC
#>     5: TGTGCCGTGAAGGGCCTCGGAGGAAGCCAAGGAAATCTCATCTTT_TGTGCCAGCAAACACATCCGGGACAGGGGGGGGTGGAGCAATCAGCCCCAGCATTTT
#>    ---                                                                                                        
#> 18440:                         TGTGCTGTGAGTGAGGGCTACAAGCTCAGCTTT_TGCGCCAGCAGCTTGGGACAGGGAGGCAACCAGCCCCAGCATTTT
#> 18441:                   TGTGCTGTGTTGGATAGCAACTATCAGTTAATCTGG_TGTGCCAGCAGCCAAGAAGGCGGGGGGGACTACAATGAGCAGTTCTTC
#> 18442:                TGTGCTGCTATGGATAGCAACTATCAGTTAATCTGG_TGCAGTGCTAAACTAGCGGGGGGAATAAGGATCCAAGAGACCCAGTACTTC
#> 18443:                      TGCATCCTGAGAGAAGGGGGGAGTGACATGCGCTTT_TGTGCCAGCAGCTTATGGGGGGGGGGGGACACTGAAGCTTTCTTT
#> 18444:                      TGTGCAGAGAGGGGTGGAGGCTTCAAAACTATCTTT_TGCGCCAGCAGCTTGCGGGACCAGCCCTACAACGAGCAGTTCTTC
#>         method     r    ts  pval pval_adj loss_a_frac loss_b_frac     score
#>         <char> <num> <num> <num>    <num>       <num>       <num>     <num>
#>     1: madhype    NA    NA    NA       NA 0.000000000   0.0000000 0.6184880
#>     2: madhype    NA    NA    NA       NA 0.000000000   0.0000000 0.6184880
#>     3: madhype    NA    NA    NA       NA 0.000000000   0.0000000 0.6184880
#>     4: madhype    NA    NA    NA       NA 0.000000000   0.0000000 0.6184880
#>     5: madhype    NA    NA    NA       NA 0.000000000   0.0000000 0.6184880
#>    ---                                                                     
#> 18440: madhype    NA    NA    NA       NA 0.000000000   0.6495726 4.6040515
#> 18441: madhype    NA    NA    NA       NA 0.000000000   0.8125000 0.1574663
#> 18442: madhype    NA    NA    NA       NA 0.008695652   0.7391304 0.5224540
#> 18443: madhype    NA    NA    NA       NA 0.000000000   0.7251908 1.0736280
#> 18444: madhype    NA    NA    NA       NA 0.000000000   0.6554054 1.1133147
#>                  cdr3a            va     ja               cdr3b          vb
#>                 <char>        <char> <char>              <char>      <char>
#>     1:     CVRAGGFKTIF   TRAV12-1*01  TRAJ9      CSATSRREPYEQYF TRBV20-1*01
#>     2:  CAMRNRSSGYALNF TRAV14/DV4*01 TRAJ41      CASSETGLSYEQYF TRBV25-1*01
#>     3:   CALSHTGGFKTIF    TRAV9-2*01  TRAJ9       CASSHDGSYEQYF  TRBV3-1*01
#>     4:  CALIRNSGNTPLVF     TRAV16*01 TRAJ29    CASSSGERRSYNEQFF  TRBV5-4*01
#>     5: CAVKGLGGSQGNLIF   TRAV12-2*01 TRAJ42 CASKHIRDRGGWSNQPQHF   TRBV28*01
#>    ---                                                                     
#> 18440:     CAVSEGYKLSF    TRAV8-6*01 TRAJ20     CASSLGQGGNQPQHF  TRBV5-1*01
#> 18441:    CAVLDSNYQLIW    TRAV1-2*01 TRAJ33    CASSQEGGGDYNEQFF  TRBV4-2*01
#> 18442:    CAAMDSNYQLIW    TRAV1-2*01 TRAJ33   CSAKLAGGIRIQETQYF TRBV20-1*01
#> 18443:    CILREGGSDMRF   TRAV26-2*01 TRAJ43     CASSLWGGGDTEAFF  TRBV7-2*01
#> 18444:    CAERGGGFKTIF      TRAV5*01  TRAJ9     CASSLRDQPYNEQFF  TRBV5-1*01
#>             jb  sample_id marker timepoint version
#>         <char>     <char> <char>    <char>  <char>
#>     1: TRBJ2-7 cd8_tp1_v2    cd8       tp1      v2
#>     2: TRBJ2-7 cd8_tp1_v2    cd8       tp1      v2
#>     3: TRBJ2-7 cd8_tp1_v2    cd8       tp1      v2
#>     4: TRBJ2-1 cd8_tp1_v2    cd8       tp1      v2
#>     5: TRBJ1-5 cd8_tp1_v2    cd8       tp1      v2
#>    ---                                            
#> 18440: TRBJ1-5 cd8_tp2_v2    cd8       tp2      v2
#> 18441: TRBJ2-1 cd8_tp2_v2    cd8       tp2      v2
#> 18442: TRBJ2-5 cd8_tp2_v2    cd8       tp2      v2
#> 18443: TRBJ1-1 cd8_tp2_v2    cd8       tp2      v2
#> 18444: TRBJ2-1 cd8_tp2_v2    cd8       tp2      v2
#>                                             label    va_orig  vb_orig
#>                                            <char>     <char>   <char>
#>     1: marker: cd8 | timepoint: tp1 | version: v2   TRAV12-1 TRBV20-1
#>     2: marker: cd8 | timepoint: tp1 | version: v2 TRAV14/DV4 TRBV25-1
#>     3: marker: cd8 | timepoint: tp1 | version: v2    TRAV9-2  TRBV3-1
#>     4: marker: cd8 | timepoint: tp1 | version: v2     TRAV16  TRBV5-4
#>     5: marker: cd8 | timepoint: tp1 | version: v2   TRAV12-2   TRBV28
#>    ---                                                               
#> 18440: marker: cd8 | timepoint: tp2 | version: v2    TRAV8-6  TRBV5-1
#> 18441: marker: cd8 | timepoint: tp2 | version: v2    TRAV1-2  TRBV4-2
#> 18442: marker: cd8 | timepoint: tp2 | version: v2    TRAV1-2 TRBV20-1
#> 18443: marker: cd8 | timepoint: tp2 | version: v2   TRAV26-2  TRBV7-2
#> 18444: marker: cd8 | timepoint: tp2 | version: v2      TRAV5  TRBV5-1
#>        is_functional
#>               <lgcl>
#>     1:          TRUE
#>     2:          TRUE
#>     3:          TRUE
#>     4:          TRUE
#>     5:          TRUE
#>    ---              
#> 18440:          TRUE
#> 18441:          TRUE
#> 18442:          TRUE
#> 18443:          TRUE
#> 18444:          TRUE

```
