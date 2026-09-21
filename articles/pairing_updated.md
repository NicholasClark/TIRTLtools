# Pairing of alpha/beta TCR chains

## Tutorial overview

This vignette shows how to generate paired TCR and single-chain
“pseudo-bulk” data using TIRTL-seq TCR-alpha and TCR-beta read counts
from half of a 384-well plate.

## Study details

For this example, we use data from the St. Jude Tracking Study of Immune
Responses Associated with COVID-19 (SJTRC), a clinical trial launched in
2020 to study the T-cell receptor repertoires of adults with naturally
acquired COVID-19. This was a prospective, longitudinal cohort study
involving adult employees (18 years and older) at St. Jude Children’s
Research Hospital in Memphis, TN, USA.

Participants underwent weekly PCR screening for SARS-CoV-2 infection
while on the St. Jude campus. Blood samples were collected in 8-ml cell
preparation tubes, processed within 24 h into cellular and plasma
components, aliquoted and then frozen for future analysis.

Here, we use samples from a healthy 33-year-old female donor with
naturally acquired mild SARS-CoV-2 infection and no prior history of
SARS-CoV-2 infection or vaccination. Samples were collected 143 days
before this donor’s first positive SARS-CoV-2 PCR test (`baseline` or
`tp1` sample), 6 days after (`acute` or `tp2` sample) and 29 days after
(`convalescent` or `tp3` sample).

## Instructions

### Select plate rows and columns

CD4+ and CD8+ T-cells were isolated from peripheral blood mononuclear
cells (PBMCs) in each sample. For each timepoint, CD8+ cells were plated
on the left half of wells in a 384-well plate and CD4+ cells were plated
on the right half.

In this case, we will pair the CD8+ cells from the baseline timepoint
(tp1), which were in columns 1 through 12 (the left half of the plate)
and rows 1 through 16 (or A through P, all rows) of our 384-well plate.

``` r
library(TIRTLtools)
```

    ## Loading required package: ggplot2

``` r
TIRTLtools:::download_data("exp3_tp1_cd8.tar.gz") ## download well-level data
```

    ## ℹ Download skipped.
    ## File exp3_tp1_cd8.tar.gz already exists in /Users/nclark2/Library/Caches/org.R-project.R/R/TIRTLtools/data-v1.
    ## Use 'force = TRUE' to force download and overwrite existing file.

``` r
data_dir = file.path(TIRTLtools:::get_cache_dir(),"exp3_tp1_cd8") ## get the directory that data was downloaded and extracted to
print(data_dir)
```

    ## [1] "/Users/nclark2/Library/Caches/org.R-project.R/R/TIRTLtools/data-v1/exp3_tp1_cd8"

``` r
dir(data_dir) |> head()
```

    ## [1] "MVP113_TIRTLseq_A1_S72.TCRa_ID02.clones_TRA.tsv.gz"  
    ## [2] "MVP113_TIRTLseq_A1_S72.TCRb_ID02.clones_TRB.tsv.gz"  
    ## [3] "MVP113_TIRTLseq_A10_S216.TCRa_ID02.clones_TRA.tsv.gz"
    ## [4] "MVP113_TIRTLseq_A10_S216.TCRb_ID02.clones_TRB.tsv.gz"
    ## [5] "MVP113_TIRTLseq_A11_S232.TCRa_ID02.clones_TRA.tsv.gz"
    ## [6] "MVP113_TIRTLseq_A11_S232.TCRb_ID02.clones_TRB.tsv.gz"

``` r
write_folder = tools::R_user_dir("TIRTLtools", which = "cache") ### change this to the folder you would like to save paired output in (if the folder does not exist it will be created)
wells = get_well_subset(row_range = 1:16, col_range = 1:12) ## wells on left half of plate
```

### Run pairing algorithms

The
[`run_pairing()`](https://nicholasclark.github.io/TIRTLtools/reference/run_pairing.md)
function will run the pairing scripts, which implement the MAD-HYPE and
T-SHELL algorithms. These algorithms are implemented in Python using
GPU-accelerated packages (with numpy as a CPU-only backup) for fast
calculation. The package should automatically install a python
environment with the necessary packages and detect the GPU type (NVIDIA
or Apple Silicon) of your machine, if it has one available.

Each function call may take a few minutes to run.

The main arguments of the function are:

- `folder_path` - the path of the folder with well-level data
- `folder_out` - the path of the folder to write results to (the
  function will create the folder if it does not exist)
- `prefix` - a prefix (e.g. the sample name) for the output file names
- `well_filter_thres` - a threshold on the number of unique clones found
  in a well, for quality control
- `well_pos` - the position of the well ID (e.g. “B5”) in the file
  names.
- `wellset` - a vector of wells to use for the pairing.

For more information, check the documentation with
[`?run_pairing`](https://nicholasclark.github.io/TIRTLtools/reference/run_pairing.md).

The well-level data files in each folder are named something like
`MVP113_TIRTLseq_P12_S263.TCRa_ID02.clones_TRA.tsv`, so the well ID
(`P12` in this case) is in the third position (separate by underscores)
so we use `well_pos=3`. We use the
[`get_well_subset()`](https://nicholasclark.github.io/TIRTLtools/reference/get_well_subset.md)
function above to generate the vector of wells for the left side of the
plate. With this, we run the
[`run_pairing()`](https://nicholasclark.github.io/TIRTLtools/reference/run_pairing.md)
function, saving the output from CD8+ cells in each plate to the same
folder, with a different prefix for each:
`<experiment>_<timepoint>_<marker>`.

``` r
### This function will attempt to automatically install a Python environment with necessary packages for GPU computation
run_pairing(folder_path = data_dir, folder_out = write_folder, prefix = "exp3_tp1_cd8", wellset = wells)
```

    ## [1] "Folder already exists: /Users/nclark2/Library/Caches/org.R-project.R/R/TIRTLtools"

    ## Loading clone files (TCRalpha) for 192 wells...

    ## 25 of 192 files loaded

    ## 50 of 192 files loaded

    ## 75 of 192 files loaded

    ## 100 of 192 files loaded

    ## 125 of 192 files loaded

    ## 150 of 192 files loaded

    ## 175 of 192 files loaded

    ## 192 of 192 files loaded

    ## Loading clone files (TCRbeta) for 192 wells...

    ## 25 of 192 files loaded

    ## 50 of 192 files loaded

    ## 75 of 192 files loaded

    ## 100 of 192 files loaded

    ## 125 of 192 files loaded

    ## 150 of 192 files loaded

    ## 175 of 192 files loaded

    ## 192 of 192 files loaded

    ## Total number of TCRalpha reads: 6.86e+07

    ## Total number of TCRbeta reads: 9.83e+07

    ## Total number of unique TCRalpha chains: 937880

    ## Total number of unique TCRbeta chains: 928782

    ## Alpha clone threshold for QC: 3217

    ## Beta clone threshold for QC: 3217

    ## Alpha wells passing QC: 191

    ## Beta wells passing QC: 191

    ## Alpha wells failing QC: 1

    ## Beta wells failing QC: 1

    ## Wells passing alpha and beta QC: 191

    ## 1 out of 192 wells removed by QC: P3

    ## Using 'auto' T-SHELL settings, selecting 'pval_thres_tshell' and 'wij_thres_tshell' based on the number of wells passing QC

    ## 191 wells passing: Using '384_well' T-SHELL settings

    ## Tabulating TCRalpha pseudobulk counts

    ## Writing TCRalpha pseudobulk file... exp3_tp1_cd8_pseudobulk_TRA.tsv

    ## Tabulating TCRbeta pseudobulk counts

    ## Writing TCRbeta pseudobulk file... exp3_tp1_cd8_pseudobulk_TRB.tsv

    ## Pseudobulk done

    ## Merging alpha clonesets...

    ## Done! Unique alpha clones after filtering: 937875

    ## Unique alpha clones in more than 2 wells: 36679

    ## Done! Unique beta clones after filtering: 928777

    ## Unique beta clones and wells in more than 2 wells: 36080

    ## Pre-computing look-up table:

    ##   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%

    ## Running pairing algorithms...

    ## Checking for available GPU...
    ## 
    ## Apple Silicon GPU detected:
    ## Apple Silicon GPU (M1/M2/M3)
    ## Checking for GPU-related Python modules...
    ## 
    ## 'mlx' is installed (for Apple Silicon GPUs).
    ## Loading mlx
    ## total number of chunks 73
    ## start time for MAD-HYPE: 2026-09-21 15:02:52
    ## Progress: 3500 (10%)
    ## Progress: 7000 (20%)
    ## Progress: 14500 (40%)
    ## Progress: 18000 (50%)
    ## Progress: 22000 (60%)
    ## Progress: 25500 (70%)
    ## Progress: 33000 (90%)
    ## Progress: 36500 (100%)
    ## end time for MAD-HYPE: 2026-09-21 15:02:55
    ## Number of pairs: 13379
    ## total number of chunks 73
    ## start processing time for T-Shell: 2026-09-21 15:02:55
    ## Progress: 3500 (10%)
    ## Progress: 7000 (20%)
    ## Progress: 14500 (40%)
    ## Progress: 18000 (50%)
    ## Progress: 22000 (60%)
    ## Progress: 25500 (70%)
    ## Progress: 33000 (90%)
    ## Progress: 36500 (100%)
    ## end time for T-Shell: 2026-09-21 15:02:56

    ## Pairing is finished.

    ## [1] "Filtering results, adding amino acid and V segment information..."
    ## [1] "Scoring unique pairs..."

    ## Writing TCRalpha/beta pairs... exp3_tp1_cd8_TIRTLoutput.tsv

    ## All pairing is finished!

    ## Number of clones paired by MAD-HYPE algorithm: 13379

    ## Number of clones paired by T-SHELL algorithm: 7304

    ## Total number of unique TCRalpha/beta pairs: 14270

    ## 81.117 sec elapsed

    ## Key: <wi, wj, wij>
    ##           wi    wj   wij                                           alpha_nuc
    ##        <num> <num> <num>                                              <char>
    ##     1:     0     3    85          TGTGCAATGAGTTTTAACTTTGGAAATGAGAAATTAACCTTT
    ##     2:     1     0    56             TGTGCTGGACATACCGGCACTGCCAGTAAACTCACCTTT
    ##     3:     1     0    56 TGTGCAGGAGCGGAGGATGCTGGTGGTACTAGCTATGGAAAGCTGACATTT
    ##     4:     0     4    68             TGTGCCGTGGACGTGTACACCGGTAACCAGTTCTATTTT
    ##     5:     0     4    68          TGCGGCACAGAAAGCGGAGGTAGCAACTATAAACTGACATTT
    ##    ---                                                                      
    ## 20679:    40     0    86                TGTGCTCTGAGTGAGACCGGTAACCAGTTCTATTTT
    ## 20680:     0     1   190         TGTGCAGAGTTGAAGATCTTATAACACCGACAAGCTCATCTTT
    ## 20681:     0     0   191            TGTGCCGTGAAACAATAACCAGGGAGGAAAGCTTATCTTC
    ## 20682:     0    24   165           TGTGCAGGAGTTGGGGAGGAAGCCAAGGAAATCTCATCTTT
    ## 20683:     0     5   186        TGTGCTCTGAGTGACAAATAGGCTTTGGGAATGTGCTGCATTGC
    ##                                                      beta_nuc    wa    wb
    ##                                                        <char> <int> <int>
    ##     1:             TGTGCCAGCAGTTTATCGTCCCTCGAGGGGGGCTACACCTTC    85    88
    ##     2:          TGCGCCAGCAGCCGAGTCAGAAATCCCGACTACGAGCAGTACTTC    57    56
    ##     3:          TGCGCCAGCAGCCGAGTCAGAAATCCCGACTACGAGCAGTACTTC    57    56
    ##     4:                TGCGCCAGCAGCTTTGACAAGACCTATGGCTACACCTTC    68    72
    ##     5:                TGCGCCAGCAGCTTTGACAAGACCTATGGCTACACCTTC    68    72
    ##    ---                                                                   
    ## 20679:                 TGTGCCAGCAGCCAGGGATCTGGCGGAGACAGTTCTTC   126    86
    ## 20680:                TGTGCCTGGAGTAAGGGTGAGGCCAACGTCCTGACTTTC   190   191
    ## 20681:             TGTGCCAGCAGCACAGGACTGCCAAGAGATACGCAGTATTTT   191   191
    ## 20682:                TGCAGTGCTAGGGGGGACAGCTCCTACGAGCAGTACTTC   165   189
    ## 20683: TGTGCCAGCAGCTTTCAGGGCGGGGGAACGGTGGGCACAGATACGCAGTATTTT   186   191
    ##                                              alpha_nuc_seq
    ##                                                     <char>
    ##     1:          TGTGCAATGAGTTTTAACTTTGGAAATGAGAAATTAACCTTT
    ##     2:             TGTGCTGGACATACCGGCACTGCCAGTAAACTCACCTTT
    ##     3: TGTGCAGGAGCGGAGGATGCTGGTGGTACTAGCTATGGAAAGCTGACATTT
    ##     4:             TGTGCCGTGGACGTGTACACCGGTAACCAGTTCTATTTT
    ##     5:          TGCGGCACAGAAAGCGGAGGTAGCAACTATAAACTGACATTT
    ##    ---                                                    
    ## 20679:                TGTGCTCTGAGTGAGACCGGTAACCAGTTCTATTTT
    ## 20680:         TGTGCAGAGTTGAAGATCTTATAACACCGACAAGCTCATCTTT
    ## 20681:            TGTGCCGTGAAACAATAACCAGGGAGGAAAGCTTATCTTC
    ## 20682:           TGTGCAGGAGTTGGGGAGGAAGCCAAGGAAATCTCATCTTT
    ## 20683:        TGTGCTCTGAGTGACAAATAGGCTTTGGGAATGTGCTGCATTGC
    ##                                                  beta_nuc_seq
    ##                                                        <char>
    ##     1:             TGTGCCAGCAGTTTATCGTCCCTCGAGGGGGGCTACACCTTC
    ##     2:          TGCGCCAGCAGCCGAGTCAGAAATCCCGACTACGAGCAGTACTTC
    ##     3:          TGCGCCAGCAGCCGAGTCAGAAATCCCGACTACGAGCAGTACTTC
    ##     4:                TGCGCCAGCAGCTTTGACAAGACCTATGGCTACACCTTC
    ##     5:                TGCGCCAGCAGCTTTGACAAGACCTATGGCTACACCTTC
    ##    ---                                                       
    ## 20679:                 TGTGCCAGCAGCCAGGGATCTGGCGGAGACAGTTCTTC
    ## 20680:                TGTGCCTGGAGTAAGGGTGAGGCCAACGTCCTGACTTTC
    ## 20681:             TGTGCCAGCAGCACAGGACTGCCAAGAGATACGCAGTATTTT
    ## 20682:                TGCAGTGCTAGGGGGGACAGCTCCTACGAGCAGTACTTC
    ## 20683: TGTGCCAGCAGCTTTCAGGGCGGGGGAACGGTGGGCACAGATACGCAGTATTTT
    ##                                                                                                 alpha_beta
    ##                                                                                                     <char>
    ##     1:               TGTGCAATGAGTTTTAACTTTGGAAATGAGAAATTAACCTTT_TGTGCCAGCAGTTTATCGTCCCTCGAGGGGGGCTACACCTTC
    ##     2:               TGTGCTGGACATACCGGCACTGCCAGTAAACTCACCTTT_TGCGCCAGCAGCCGAGTCAGAAATCCCGACTACGAGCAGTACTTC
    ##     3:   TGTGCAGGAGCGGAGGATGCTGGTGGTACTAGCTATGGAAAGCTGACATTT_TGCGCCAGCAGCCGAGTCAGAAATCCCGACTACGAGCAGTACTTC
    ##     4:                     TGTGCCGTGGACGTGTACACCGGTAACCAGTTCTATTTT_TGCGCCAGCAGCTTTGACAAGACCTATGGCTACACCTTC
    ##     5:                  TGCGGCACAGAAAGCGGAGGTAGCAACTATAAACTGACATTT_TGCGCCAGCAGCTTTGACAAGACCTATGGCTACACCTTC
    ##    ---                                                                                                    
    ## 20679:                         TGTGCTCTGAGTGAGACCGGTAACCAGTTCTATTTT_TGTGCCAGCAGCCAGGGATCTGGCGGAGACAGTTCTTC
    ## 20680:                 TGTGCAGAGTTGAAGATCTTATAACACCGACAAGCTCATCTTT_TGTGCCTGGAGTAAGGGTGAGGCCAACGTCCTGACTTTC
    ## 20681:                 TGTGCCGTGAAACAATAACCAGGGAGGAAAGCTTATCTTC_TGTGCCAGCAGCACAGGACTGCCAAGAGATACGCAGTATTTT
    ## 20682:                   TGTGCAGGAGTTGGGGAGGAAGCCAAGGAAATCTCATCTTT_TGCAGTGCTAGGGGGGACAGCTCCTACGAGCAGTACTTC
    ## 20683: TGTGCTCTGAGTGACAAATAGGCTTTGGGAATGTGCTGCATTGC_TGTGCCAGCAGCTTTCAGGGCGGGGGAACGGTGGGCACAGATACGCAGTATTTT
    ##         method         r       ts         pval     pval_adj loss_a_frac
    ##         <char>     <num>    <num>        <num>        <num>       <num>
    ##     1: madhype        NA       NA           NA           NA 0.034090909
    ##     2: madhype        NA       NA           NA           NA 0.000000000
    ##     3: madhype        NA       NA           NA           NA 0.000000000
    ##     4: madhype        NA       NA           NA           NA 0.055555556
    ##     5: madhype        NA       NA           NA           NA 0.055555556
    ##    ---                                                                 
    ## 20679:  tshell 0.5611333 9.319887 3.088615e-17 8.477475e-11 0.000000000
    ## 20680:  tshell 0.5594879 9.280079 3.993057e-17 2.274191e-11 0.005235602
    ## 20681:  tshell 0.5515199 9.089538 1.357749e-16 9.696026e-13 0.000000000
    ## 20682:  tshell 0.5426992 8.882772 5.069335e-16 1.223313e-11 0.126984127
    ## 20683:  tshell 0.5270523 8.526113 4.781799e-15 3.515918e-11 0.026178010
    ##        loss_b_frac     score             cdr3a         va     ja
    ##              <num>     <num>            <char>     <char> <char>
    ##     1:  0.00000000 45.705465    CAMSFNFGNEKLTF TRAV14/DV4 TRAJ48
    ##     2:  0.01754386 42.204274     CAGHTGTASKLTF     TRAV21 TRAJ44
    ##     3:  0.01754386 42.204274 CAGAEDAGGTSYGKLTF     TRAV27 TRAJ52
    ##     4:  0.00000000 41.731629     CAVDVYTGNQFYF     TRAV39 TRAJ49
    ##     5:  0.00000000 41.731629    CGTESGGSNYKLTF     TRAV30 TRAJ53
    ##    ---                                                          
    ## 20679:  0.31746032 17.460590      CALSETGNQFYF     TRAV19 TRAJ49
    ## 20680:  0.00000000      -Inf   CAELKIL_NTDKLIF      TRAV5 TRAJ34
    ## 20681:  0.00000000 -4.560842    CAVKQ*P_GGKLIF   TRAV12-2 TRAJ23
    ## 20682:  0.00000000 -3.815652    CAGVGEE_QGNLIF     TRAV27 TRAJ42
    ## 20683:  0.00000000      -Inf   CALSDK*_FGNVLHC     TRAV19 TRAJ35
    ##                     cdr3b       vb      jb is_functional
    ##                    <char>   <char>  <char>        <lgcl>
    ##     1:     CASSLSSLEGGYTF   TRBV27 TRBJ1-2          TRUE
    ##     2:    CASSRVRNPDYEQYF  TRBV4-1 TRBJ2-7          TRUE
    ##     3:    CASSRVRNPDYEQYF  TRBV4-1 TRBJ2-7          TRUE
    ##     4:      CASSFDKTYGYTF  TRBV5-1 TRBJ1-2          TRUE
    ##     5:      CASSFDKTYGYTF  TRBV5-1 TRBJ1-2          TRUE
    ##    ---                                                  
    ## 20679:      CASSQG_WRRQFF  TRBV4-2 TRBJ2-1         FALSE
    ## 20680:      CAWSKGEANVLTF   TRBV30 TRBJ2-6         FALSE
    ## 20681:     CASSTGLPRDTQYF  TRBV5-8 TRBJ2-3         FALSE
    ## 20682:      CSARGDSSYEQYF TRBV20-1 TRBJ2-7         FALSE
    ## 20683: CASSFQGGGTVGTDTQYF  TRBV7-9 TRBJ2-3         FALSE

### Load paired TCR data

The output from each call to
[`run_pairing()`](https://nicholasclark.github.io/TIRTLtools/reference/run_pairing.md)
is 3 tab-separated text files (.tsv). If you would like to compress the
output files, you can use the option `gzip_output = TRUE`.

- `<sample-name>_TIRTLoutput.tsv(.gz)`– Paired TCRs – a table of
  computationally paired alpha/beta TCRs
- `<sample-name>_pseudobulk_TRA.tsv(.gz)` – Alpha chain “pseudo-bulk”
  data – a table of read count summary metrics for alpha-chains across
  all wells on a plate
- `<sample-name>_pseudobulk_TRB.tsv(.gz)` – Beta chain “pseudo-bulk”
  data – a table of read count summary metrics for beta-chains across
  all wells on a plate

The
[`load_tirtlseq()`](https://nicholasclark.github.io/TIRTLtools/reference/load_tirtlseq.md)
function can be used to load all of the data in a directory, whether it
is in `.tsv` or compressed `.tsv.gz` form.

This function loads the data into a “list”, which is an R data type that
can contain any type of unstructured data. Elements of a list can be
accessed by their name, e.g. `your_list$meta`, or their index,
e.g. `your_list[[2]]`. The list contains three slots: `meta`, `data`,
and `info`, and has the following structure:

    your_data_object (list, with class "TIRTLseqDataSet")
    ├───meta (a data frame of sample metadata)
    ├───info (a list of processing/pairing information)
    └───data (a list containing paired and pseudobulk data for each sample)
        └───sample_1 (list)
            ├───alpha (alpha pseudobulk dataframe)
            ├───beta (beta pseudobulk dataframe)
            └───paired (paired pseudobulk dataframe)
        ...
        └───sample_n (list)
            ├───alpha (alpha pseudobulk dataframe)
            ├───beta (beta pseudobulk dataframe)
            └───paired (paired pseudobulk dataframe)

Since we used `<experiment>_<timepoint>_<marker>` for the `prefix` for
each set of output files, we add this in the `meta_columns` argument, so
that the function creates a data frame of sample metadata with this
information.

``` r
ts = load_tirtlseq(write_folder, sep = "_", verbose = FALSE,
                     meta_columns = c("experiment", "timepoint", "marker"))
ts ## list with class "TIRTLseqDataSet"
```

    ## <TIRTLseqDataSet>
    ##   Number of samples: 1 
    ##   Samples: exp3_tp1_cd8

``` r
summary(ts)
```

    ## <TIRTLseqDataSet>
    ##   Number of samples: 1 
    ##   Per-sample TCR counts (alpha / beta / paired):
    ## # A tibble: 1 × 5
    ##   sample_id    `Paired TCRs` `Unique Paired TCRs` `Alpha Chains` `Beta Chains`
    ##   <chr>                <int>                <int>          <int>         <int>
    ## 1 exp3_tp1_cd8         20683                14270         937875        928777

You can inspect the loaded object and view the data frames with paired
TCRs and single-chain pseudobulk for each sample.

``` r
ts$meta ## sample metadata (data frame)
```

    ## # A tibble: 1 × 5
    ##   sample_id    experiment timepoint marker label                                
    ##   <chr>        <chr>      <chr>     <chr>  <chr>                                
    ## 1 exp3_tp1_cd8 exp3       tp1       cd8    experiment: exp3 | timepoint: tp1 | …

``` r
ts$info ## list of data pairing/processing information
```

    ## $data_version
    ## [1] "v1"
    ## 
    ## $processing_version
    ## [1] "v1"
    ## 
    ## $package_version
    ## [1] '0.2.6'
    ## 
    ## $call
    ## $call$directory
    ## [1] "/Users/nclark2/Library/Caches/org.R-project.R/R/TIRTLtools"
    ## 
    ## $call$chain
    ## [1] "all"    "paired" "alpha"  "beta"  
    ## 
    ## $call$sep
    ## [1] "_"
    ## 
    ## $call$meta_columns
    ## [1] "experiment" "timepoint"  "marker"    
    ## 
    ## $call$samples
    ## NULL
    ## 
    ## $call$clean
    ## [1] FALSE
    ## 
    ## $call$remove_nonfunctional
    ## [1] FALSE
    ## 
    ## $call$process
    ## [1] TRUE
    ## 
    ## $call$pseudobulk_columns
    ## [1] "auto"
    ## 
    ## $call$paired_columns
    ## [1] "auto"
    ## 
    ## $call$n_threads
    ## [1] 9
    ## 
    ## $call$verbose
    ## [1] FALSE
    ## 
    ## $call$stringsAsFactors
    ## [1] FALSE
    ## 
    ## $call$n_max
    ## [1] Inf

``` r
names(ts$data) ## list with data for each sample
```

    ## [1] "exp3_tp1_cd8"

``` r
names(ts$data$exp3_tp1_cd8) ## list with data for one sample (exp3_tp1_cd8)
```

    ## [1] "alpha"      "beta"       "paired"     "paired_alt"

``` r
ts$data$exp3_tp1_cd8$alpha ## TCR-alpha pseudobulk (data frame)
```

    ##                   id   rank             aaSeqCDR3            v      j
    ##               <char>  <int>                <char>       <char> <char>
    ##      1:      alpha_1      1       CAMRGRSNRDDKIIF   TRAV14/DV4 TRAJ30
    ##      2:      alpha_2      2          CALKTSYDKVIF TRAV38-2/DV8 TRAJ50
    ##      3:      alpha_3      3          CAASRLPDDMRF   TRAV29/DV5 TRAJ43
    ##      4:      alpha_4      4          CAVRDSNYQLIW      TRAV1-2 TRAJ33
    ##      5:      alpha_5      5 CLLAHLRSGPNSGGSNYKLTF       TRAV40 TRAJ53
    ##     ---                                                              
    ## 937871: alpha_937871 937871         FFLPNR_FQRLVF      TRAV1-2  TRAJ8
    ## 937872: alpha_937872 937872          FFQGPQ_NKLTF     TRAV13-1 TRAJ17
    ## 937873: alpha_937873 937873          FFQGLQ_NKLTF     TRAV13-1 TRAJ17
    ## 937874: alpha_937874 937874          FFQGLQ_NKLTF     TRAV13-1 TRAJ17
    ## 937875: alpha_937875 937875        FFLQRGW_SYQLTF      TRAV8-1 TRAJ28
    ##         readFraction readCount n_wells paired_status n_paired
    ##                <num>     <int>   <int>        <char>    <num>
    ##      1: 2.058548e-02   1411609     191  T-SHELL only        1
    ##      2: 1.482005e-02   1016256     191  T-SHELL only        1
    ##      3: 9.557239e-03    655369     191  T-SHELL only        1
    ##      4: 5.478071e-03    375648     191       neither        0
    ##      5: 4.698377e-03    322182     191  T-SHELL only        1
    ##     ---                                                      
    ## 937871: 1.458299e-08         1       1       neither        0
    ## 937872: 1.458299e-08         1       1       neither        0
    ## 937873: 1.458299e-08         1       1       neither        0
    ## 937874: 1.458299e-08         1       1       neither        0
    ## 937875: 1.458299e-08         1       1       neither        0
    ##                                                         targetSequences
    ##                                                                  <char>
    ##      1:                   TGTGCAATGAGGGGTCGATCGAACAGAGATGACAAGATCATCTTT
    ##      2:                            TGTGCTTTGAAAACCTCCTACGACAAGGTGATATTT
    ##      3:                            TGTGCAGCAAGCCGACTTCCGGATGACATGCGCTTT
    ##      4:                            TGTGCTGTGAGAGATAGCAACTATCAGTTAATCTGG
    ##      5: TGTCTTCTGGCACACCTCCGATCAGGACCCAATAGTGGAGGTAGCAACTATAAACTGACATTT
    ##     ---                                                                
    ## 937871:                           TTTTTCCTCCCGAACAGGCTTTCAGAGACTTGTATTT
    ## 937872:                             TTTTTTCAGGGGCCGCAGGCAACAAGCTAACTTTT
    ## 937873:                             TTTTTTCAGGGGCTGCAGGCAACAAGCTGACTTTT
    ## 937874:                             TTTTTTCAGGGGCTGCAGGCAACAAGTTAACTTTT
    ## 937875:                        TTTTTTTTACAACGGGGCTGGGAGTTACCAACTCACTTTC
    ##         readCount_median readCount_max          sem max_wells n_paired_madhype
    ##                    <num>         <int>        <num>     <int>            <num>
    ##      1:             7413         14066 2.138826e-04       191                0
    ##      2:             5305          9477 2.019982e-04       191                0
    ##      3:             3310          6673 1.399961e-04       191                0
    ##      4:             1909          3730 9.109789e-05       191                0
    ##      5:             1705          2812 7.981121e-05       191                0
    ##     ---                                                                       
    ## 937871:                1             1 1.202901e-08       191                0
    ## 937872:                1             1 1.449137e-08       191                0
    ## 937873:                1             1 1.518653e-08       191                0
    ## 937874:                1             1 1.518653e-08       191                0
    ## 937875:                1             1 1.105502e-08       191                0
    ##         n_paired_tshell has_stop_codon has_frameshift is_functional is_paired
    ##                   <num>         <lgcl>         <lgcl>        <lgcl>    <lgcl>
    ##      1:               1          FALSE          FALSE          TRUE      TRUE
    ##      2:               1          FALSE          FALSE          TRUE      TRUE
    ##      3:               1          FALSE          FALSE          TRUE      TRUE
    ##      4:               0          FALSE          FALSE          TRUE     FALSE
    ##      5:               1          FALSE          FALSE          TRUE      TRUE
    ##     ---                                                                      
    ## 937871:               0          FALSE           TRUE         FALSE     FALSE
    ## 937872:               0          FALSE           TRUE         FALSE     FALSE
    ## 937873:               0          FALSE           TRUE         FALSE     FALSE
    ## 937874:               0          FALSE           TRUE         FALSE     FALSE
    ## 937875:               0          FALSE           TRUE         FALSE     FALSE
    ##         is_paired_tshell is_paired_madhype
    ##                   <lgcl>            <lgcl>
    ##      1:             TRUE             FALSE
    ##      2:             TRUE             FALSE
    ##      3:             TRUE             FALSE
    ##      4:            FALSE             FALSE
    ##      5:             TRUE             FALSE
    ##     ---                                   
    ## 937871:            FALSE             FALSE
    ## 937872:            FALSE             FALSE
    ## 937873:            FALSE             FALSE
    ## 937874:            FALSE             FALSE
    ## 937875:            FALSE             FALSE

``` r
ts$data$exp3_tp1_cd8$beta ## TCR-beta pseudobulk (data frame)
```

    ##                  id   rank          aaSeqCDR3        v       j readFraction
    ##              <char>  <int>             <char>   <char>  <char>        <num>
    ##      1:      beta_1      1      CSVERVNYNEQFF TRBV29-1 TRBJ2-1 1.348328e-02
    ##      2:      beta_2      2     CASSTGLPRDTQYF  TRBV5-8 TRBJ2-3 1.143658e-02
    ##      3:      beta_3      3      CSATAGRNYGYTF TRBV20-1 TRBJ1-2 1.065993e-02
    ##      4:      beta_4      4        CASSLTYEQYF  TRBV6-2 TRBJ2-7 8.413453e-03
    ##      5:      beta_5      5  CASSWPQGSGSLDEQFF  TRBV5-1 TRBJ2-1 5.177560e-03
    ##     ---                                                                    
    ## 928773: beta_928773 928773     FASSLGV_GDEQFF TRBV12-3 TRBJ2-1 1.017017e-08
    ## 928774: beta_928774 928774 FASSIGDRA_LRNSPLHF   TRBV19 TRBJ1-6 1.017017e-08
    ## 928775: beta_928775 928775        FAISE_SEAVF TRBV10-3 TRBJ2-2 1.017017e-08
    ## 928776: beta_928776 928776     FARKQGN_DNEQFF   TRBV19 TRBJ2-1 1.017017e-08
    ## 928777: beta_928777 928777    FAADRGP_STDTQYF   TRBV19 TRBJ2-3 1.017017e-08
    ##         readCount n_wells paired_status n_paired
    ##             <int>   <int>        <char>    <num>
    ##      1:   1325768     191  T-SHELL only        2
    ##      2:   1124522     191  T-SHELL only        2
    ##      3:   1048157     191  T-SHELL only        1
    ##      4:    827268     191  T-SHELL only        2
    ##      5:    509093     191  T-SHELL only        1
    ##     ---                                         
    ## 928773:         1       1       neither        0
    ## 928774:         1       1       neither        0
    ## 928775:         1       1       neither        0
    ## 928776:         1       1       neither        0
    ## 928777:         1       1       neither        0
    ##                                              targetSequences readCount_median
    ##                                                       <char>            <num>
    ##      1:              TGCAGCGTTGAAAGGGTTAATTACAATGAGCAGTTCTTC             6728
    ##      2:           TGTGCCAGCAGCACAGGACTGCCAAGAGATACGCAGTATTTT             5933
    ##      3:              TGCAGTGCTACGGCAGGGAGGAACTATGGCTACACCTTC             5267
    ##      4:                    TGTGCCAGCAGTCTCACCTACGAGCAGTACTTC             4293
    ##      5:  TGCGCCAGCAGCTGGCCTCAGGGTAGCGGGAGTTTGGATGAGCAGTTCTTC             2516
    ##     ---                                                                      
    ## 928773:            TTTGCCAGCAGTTTGGGGGTATAGGGGATGAGCAGTTCTTC                1
    ## 928774: TTTGCCAGTAGTATCGGGGACAGGGCTCTTGAGAAATTCACCCCTCCACTTT                1
    ## 928775:                     TTTGCCATCAGTGAAAGAGCGAAGCTGTTTTT                1
    ## 928776:            TTTGCCCGCAAACAGGGAAATGGGACAATGAGCAGTTCTTC                1
    ## 928777:          TTTGCCGCCGACAGGGGGCCTAAGCACAGATACGCAGTATTTT                1
    ##         readCount_max          sem max_wells n_paired_madhype n_paired_tshell
    ##                 <int>        <num>     <int>            <num>           <num>
    ##      1:         14705 2.637763e-04       191                0               2
    ##      2:          9823 1.293336e-04       191                0               2
    ##      3:          9841 1.409173e-04       191                0               1
    ##      4:          8641 1.935800e-04       191                0               2
    ##      5:          6879 1.404382e-04       191                0               1
    ##     ---                                                                      
    ## 928773:             1 8.824574e-09       191                0               0
    ## 928774:             1 1.013799e-08       191                0               0
    ## 928775:             1 1.229442e-08       191                0               0
    ## 928776:             1 8.947667e-09       191                0               0
    ## 928777:             1 1.029493e-08       191                0               0
    ##         has_stop_codon has_frameshift is_functional is_paired is_paired_tshell
    ##                 <lgcl>         <lgcl>        <lgcl>    <lgcl>           <lgcl>
    ##      1:          FALSE          FALSE          TRUE      TRUE             TRUE
    ##      2:          FALSE          FALSE          TRUE      TRUE             TRUE
    ##      3:          FALSE          FALSE          TRUE      TRUE             TRUE
    ##      4:          FALSE          FALSE          TRUE      TRUE             TRUE
    ##      5:          FALSE          FALSE          TRUE      TRUE             TRUE
    ##     ---                                                                       
    ## 928773:          FALSE           TRUE         FALSE     FALSE            FALSE
    ## 928774:          FALSE           TRUE         FALSE     FALSE            FALSE
    ## 928775:          FALSE           TRUE         FALSE     FALSE            FALSE
    ## 928776:          FALSE           TRUE         FALSE     FALSE            FALSE
    ## 928777:          FALSE           TRUE         FALSE     FALSE            FALSE
    ##         is_paired_madhype
    ##                    <lgcl>
    ##      1:             FALSE
    ##      2:             FALSE
    ##      3:             FALSE
    ##      4:             FALSE
    ##      5:             FALSE
    ##     ---                  
    ## 928773:             FALSE
    ## 928774:             FALSE
    ## 928775:             FALSE
    ## 928776:             FALSE
    ## 928777:             FALSE

``` r
ts$data$exp3_tp1_cd8$paired ## paired TCRαβ (data frame)
```

    ## Key: <beta_nuc>
    ##         method           va     ja            cdr3a              cdr3b       vb
    ##         <char>       <char> <char>           <char>             <char>   <char>
    ##     1:  tshell     TRAV13-1 TRAJ32 CAAKGVYGGATNKLIF      CSVERVNYNEQFF TRBV29-1
    ##     2:  tshell       TRAV10 TRAJ10   CVVNHGR_GNKLTF      CSVERVNYNEQFF TRBV29-1
    ##     3:  tshell   TRAV14/DV4 TRAJ30  CAMRGRSNRDDKIIF     CASSTGLPRDTQYF  TRBV5-8
    ##     4:  tshell     TRAV12-2 TRAJ23   CAVKQ*P_GGKLIF     CASSTGLPRDTQYF  TRBV5-8
    ##     5:  tshell   TRAV29/DV5 TRAJ43     CAASRLPDDMRF      CSATAGRNYGYTF TRBV20-1
    ##    ---                                                                         
    ## 20679:  tshell       TRAV22  TRAJ6   CAVLASGGSYTPTF      CASSFG_DRGAVF  TRBV5-4
    ## 20680: madhype        TRAV3 TRAJ34    CAVRDMSTDELIF     CASSLVR_GKNCFF   TRBV27
    ## 20681:  tshell        TRAV3 TRAJ34    CAVRDMSTDELIF     CASSLVR_GKNCFF   TRBV27
    ## 20682:  tshell TRAV38-2/DV8 TRAJ45   CAYRSGGGADGLTF CASSLFM*G_TGNQPQHF   TRBV27
    ## 20683: madhype     TRAV12-3 TRAJ24    CAMNTDSWGKLQL   CASSFWGR_QMKKLFF   TRBV28
    ##             jb alpha_rank beta_rank alpha_readFraction beta_readFraction
    ##         <char>      <int>     <int>              <num>             <num>
    ##     1: TRBJ2-1         12         1       2.942424e-03      1.348328e-02
    ##     2: TRBJ2-1         79         1       7.625154e-04      1.348328e-02
    ##     3: TRBJ2-3          1         2       2.058548e-02      1.143658e-02
    ##     4: TRBJ2-3         30         2       1.645603e-03      1.143658e-02
    ##     5: TRBJ1-2          3         3       9.557239e-03      1.065993e-02
    ##    ---                                                                  
    ## 20679: TRBJ2-2     375555    572602       2.041619e-07      3.051050e-08
    ## 20680: TRBJ1-4     569015    572764       4.374897e-08      3.051050e-08
    ## 20681: TRBJ1-4     569015    572764       4.374897e-08      3.051050e-08
    ## 20682: TRBJ1-5       6296    574415       9.726854e-06      3.051050e-08
    ## 20683: TRBJ1-4     433237    574511       1.166639e-07      3.051050e-08
    ##        alpha_readCount beta_readCount is_functional
    ##                  <int>          <int>        <lgcl>
    ##     1:          201771        1325768          TRUE
    ##     2:           52288        1325768         FALSE
    ##     3:         1411609        1124522          TRUE
    ##     4:          112844        1124522         FALSE
    ##     5:          655369        1048157          TRUE
    ##    ---                                             
    ## 20679:              14              3         FALSE
    ## 20680:               3              3         FALSE
    ## 20681:               3              3         FALSE
    ## 20682:             667              3         FALSE
    ## 20683:               8              3         FALSE
    ##                                               alpha_nuc
    ##                                                  <char>
    ##     1: TGTGCAGCAAAGGGGGTTTATGGTGGTGCTACAAACAAGCTCATCTTT
    ##     2:         TGTGTGGTGAACCACGGGAGGAGGAAACAAACTCACCTTT
    ##     3:    TGTGCAATGAGGGGTCGATCGAACAGAGATGACAAGATCATCTTT
    ##     4:         TGTGCCGTGAAACAATAACCAGGGAGGAAAGCTTATCTTC
    ##     5:             TGTGCAGCAAGCCGACTTCCGGATGACATGCGCTTT
    ##    ---                                                 
    ## 20679:       TGTGCTGTACTCGCATCAGGAGGAAGCTACACACCTACATTT
    ## 20680:          TGTGCTGTGAGAGACATGAGCACCGACGAGCTCATCTTT
    ## 20681:          TGTGCTGTGAGAGACATGAGCACCGACGAGCTCATCTTT
    ## 20682:       TGTGCTTATAGGAGCGGAGGAGGTGCTGACGGACTCACCTTT
    ## 20683:          TGTGCAATGAATACGGACAGCTGGGGGAAATTGCAGCTT
    ##                                                    beta_nuc
    ##                                                      <char>
    ##     1:              TGCAGCGTTGAAAGGGTTAATTACAATGAGCAGTTCTTC
    ##     2:              TGCAGCGTTGAAAGGGTTAATTACAATGAGCAGTTCTTC
    ##     3:           TGTGCCAGCAGCACAGGACTGCCAAGAGATACGCAGTATTTT
    ##     4:           TGTGCCAGCAGCACAGGACTGCCAAGAGATACGCAGTATTTT
    ##     5:              TGCAGTGCTACGGCAGGGAGGAACTATGGCTACACCTTC
    ##    ---                                                     
    ## 20679:               TGTGCCAGCAGCTTCGGGGGGACCGGGGAGCTGTTTTT
    ## 20680:             TGTGCCAGCAGCTTGGTCCGAGGGGAAAAACTGTTTTTTT
    ## 20681:             TGTGCCAGCAGCTTGGTCCGAGGGGAAAAACTGTTTTTTT
    ## 20682: TGTGCCAGCAGTTTATTTATGTAGGGTAACGGGCAATCAGCCCCAGCATTTC
    ## 20683:       TGTGCCAGCAGTTTCTGGGGCAGGCCAGATGAAAAAACTGTTTTTT
    ##                                                                                             alpha_beta
    ##                                                                                                 <char>
    ##     1:        TGTGCAGCAAAGGGGGTTTATGGTGGTGCTACAAACAAGCTCATCTTT_TGCAGCGTTGAAAGGGTTAATTACAATGAGCAGTTCTTC
    ##     2:                TGTGTGGTGAACCACGGGAGGAGGAAACAAACTCACCTTT_TGCAGCGTTGAAAGGGTTAATTACAATGAGCAGTTCTTC
    ##     3:        TGTGCAATGAGGGGTCGATCGAACAGAGATGACAAGATCATCTTT_TGTGCCAGCAGCACAGGACTGCCAAGAGATACGCAGTATTTT
    ##     4:             TGTGCCGTGAAACAATAACCAGGGAGGAAAGCTTATCTTC_TGTGCCAGCAGCACAGGACTGCCAAGAGATACGCAGTATTTT
    ##     5:                    TGTGCAGCAAGCCGACTTCCGGATGACATGCGCTTT_TGCAGTGCTACGGCAGGGAGGAACTATGGCTACACCTTC
    ##    ---                                                                                                
    ## 20679:               TGTGCTGTACTCGCATCAGGAGGAAGCTACACACCTACATTT_TGTGCCAGCAGCTTCGGGGGGACCGGGGAGCTGTTTTT
    ## 20680:                TGTGCTGTGAGAGACATGAGCACCGACGAGCTCATCTTT_TGTGCCAGCAGCTTGGTCCGAGGGGAAAAACTGTTTTTTT
    ## 20681:                TGTGCTGTGAGAGACATGAGCACCGACGAGCTCATCTTT_TGTGCCAGCAGCTTGGTCCGAGGGGAAAAACTGTTTTTTT
    ## 20682: TGTGCTTATAGGAGCGGAGGAGGTGCTGACGGACTCACCTTT_TGTGCCAGCAGTTTATTTATGTAGGGTAACGGGCAATCAGCCCCAGCATTTC
    ## 20683:          TGTGCAATGAATACGGACAGCTGGGGGAAATTGCAGCTT_TGTGCCAGCAGTTTCTGGGGCAGGCCAGATGAAAAAACTGTTTTTT
    ##           wi    wj   wij    wa    wb       score         r         ts
    ##        <int> <int> <int> <int> <int>       <num>     <num>      <num>
    ##     1:     0     0   191   191   191 -4.56084201 0.8378690  21.101664
    ##     2:     0     0   191   191   191 -4.56084201 0.6553487  11.928026
    ##     3:     0     0   191   191   191 -4.56084201 0.7361605  14.953353
    ##     4:     0     0   191   191   191 -4.56084201 0.5515199   9.089538
    ##     5:     0     0   191   191   191 -4.56084201 0.8408442  21.356329
    ##    ---                                                               
    ## 20679:     1     0     3     4     3  0.07990462 0.9164956  31.495808
    ## 20680:     0     0     3     3     3  0.61848796        NA         NA
    ## 20681:     0     0     3     3     3  0.61848796 0.9989229 295.968667
    ## 20682:     1     0     3     4     3  0.07990462 0.8987059  28.172642
    ## 20683:     0     0     3     3     3  0.61848796        NA         NA
    ##                 pval      pval_adj loss_a_frac loss_b_frac alpha_readCount_max
    ##                <num>         <num>       <num>       <num>               <int>
    ##     1:  1.409514e-51  1.580244e-42           0        0.00                2192
    ##     2:  8.149960e-25  1.794423e-19           0        0.00                 622
    ##     3:  7.114002e-34  4.904156e-30           0        0.00               14066
    ##     4:  1.357749e-16  9.696026e-13           0        0.00                1002
    ##     5:  2.843320e-52  2.498831e-48           0        0.00                6673
    ##    ---                                                                        
    ## 20679:  3.972705e-77  4.643735e-30           0        0.25                   6
    ## 20680:            NA            NA           0        0.00                   1
    ## 20681: 5.436759e-254 6.601289e-221           0        0.00                   1
    ## 20682:  1.415573e-69  6.590745e-44           0        0.25                 186
    ## 20683:            NA            NA           0        0.00                   6
    ##        alpha_readCount_median    alpha_sem alpha_max_wells beta_readCount_max
    ##                         <num>        <num>           <int>              <int>
    ##     1:                   1024 6.329612e-05             191              14705
    ##     2:                    268 2.130842e-05             191              14705
    ##     3:                   7413 2.138826e-04             191               9823
    ##     4:                    580 2.581459e-05             191               9823
    ##     5:                   3310 1.399961e-04             191               9841
    ##    ---                                                                       
    ## 20679:                      3 1.509696e-07             191                  1
    ## 20680:                      1 2.784023e-08             191                  1
    ## 20681:                      1 2.784023e-08             191                  1
    ## 20682:                    166 5.002365e-06             191                  1
    ## 20683:                      1 8.189881e-08             191                  1
    ##        beta_readCount_median     beta_sem beta_max_wells alpha_has_stop_codon
    ##                        <num>        <num>          <int>               <lgcl>
    ##     1:                  6728 2.637763e-04            191                FALSE
    ##     2:                  6728 2.637763e-04            191                FALSE
    ##     3:                  5933 1.293336e-04            191                FALSE
    ##     4:                  5933 1.293336e-04            191                 TRUE
    ##     5:                  5267 1.409173e-04            191                FALSE
    ##    ---                                                                       
    ## 20679:                     1 3.114027e-08            191                FALSE
    ## 20680:                     1 1.957507e-08            191                FALSE
    ## 20681:                     1 1.957507e-08            191                FALSE
    ## 20682:                     1 1.859270e-08            191                FALSE
    ## 20683:                     1 1.601161e-08            191                FALSE
    ##        alpha_has_frameshift beta_has_stop_codon beta_has_frameshift
    ##                      <lgcl>              <lgcl>              <lgcl>
    ##     1:                FALSE               FALSE               FALSE
    ##     2:                 TRUE               FALSE               FALSE
    ##     3:                FALSE               FALSE               FALSE
    ##     4:                 TRUE               FALSE               FALSE
    ##     5:                FALSE               FALSE               FALSE
    ##    ---                                                             
    ## 20679:                FALSE               FALSE                TRUE
    ## 20680:                FALSE               FALSE                TRUE
    ## 20681:                FALSE               FALSE                TRUE
    ## 20682:                FALSE                TRUE                TRUE
    ## 20683:                FALSE               FALSE                TRUE
    ##        alpha_is_functional beta_is_functional     alpha_id     beta_id
    ##                     <lgcl>             <lgcl>       <char>      <char>
    ##     1:                TRUE               TRUE     alpha_12      beta_1
    ##     2:               FALSE               TRUE     alpha_79      beta_1
    ##     3:                TRUE               TRUE      alpha_1      beta_2
    ##     4:               FALSE               TRUE     alpha_30      beta_2
    ##     5:                TRUE               TRUE      alpha_3      beta_3
    ##    ---                                                                
    ## 20679:                TRUE              FALSE alpha_375555 beta_572602
    ## 20680:                TRUE              FALSE alpha_569015 beta_572764
    ## 20681:                TRUE              FALSE alpha_569015 beta_572764
    ## 20682:                TRUE              FALSE   alpha_6296 beta_574415
    ## 20683:                TRUE              FALSE alpha_433237 beta_574511
    ##        receptor_rank     receptor_id
    ##                <int>          <char>
    ##     1:             1      receptor_1
    ##     2:             1      receptor_1
    ##     3:             2      receptor_2
    ##     4:             2      receptor_2
    ##     5:             3      receptor_3
    ##    ---                              
    ## 20679:        572602 receptor_572602
    ## 20680:        572764 receptor_572764
    ## 20681:        572764 receptor_572764
    ## 20682:        574415 receptor_574415
    ## 20683:        574511 receptor_574511
