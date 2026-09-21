# Find expanded or contracted clones

**\[experimental\]**

This function returns clones that expand or contract between two
samples, based on their single-chain pseudo-bulk frequencies (beta chain
frequency is default). It returns data frames containing the most
expanded/contracted clones along with their direction and
log-fold-change. It returns both a single-chain data frame (α or β,
default is β) and a data frame containing all αβ TCR pairs corresponding
to these chains.

## Usage

``` r
get_expanded_clones(
  data1,
  data2,
  chain = c("beta", "alpha"),
  filter_pairs = FALSE,
  remove_nonfunctional = TRUE,
  log2fc_cutoff = 3,
  sem_cutoff = 2.5,
  pseudo1 = 1e-06,
  pseudo2 = 1e-06,
  smooth_sem = c("window", "none"),
  window_size = 30,
  end_window_size = 5
)
```

## Arguments

- data1:

  a sample from a TIRTLseqDataSet object (e.g.
  `<object>$data$<sample_tp1>`) - used for before frequencies

- data2:

  a sample from a TIRTLseqDataSet object (e.g.
  `<object>$data$<sample_tp2>`) - used for after frequencies

- chain:

  which chain to plot, alpha or beta (default is beta)

- filter_pairs:

  whether to keep only two alpha chains per unique beta chain. (default
  is FALSE, keep all pairs)

- remove_nonfunctional:

  whether to remove pairs with non-functional chains – e.g. that have
  stop codons or frameshifts. (default is FALSE, keep all pairs)

- log2fc_cutoff:

  the log2 fold-change cutoff to call a TCR expanded or contracted
  (default 1.5)

- sem_cutoff:

  the standard-error of the mean (SEM) to use as a cutoff in calling
  clones expanded or contracted (default is 2.5)

- pseudo1:

  the pseudocount to add to read frequency of the first sample (default
  is `10^-6`).

- pseudo2:

  the pseudocount to add to read frequency of the second sample (default
  is `10^-6`).

- smooth_sem:

  if "window", then SEM values for clones will be smoothed by comparing
  to other clones within a window of similar frequencies. Otherwise, no
  smoothing. (default is "window")

- window_size:

  the number of similar clones to include within a window.

- end_window_size:

  the number of clones to include in a window at the ends (most and
  least frequent)

## Value

A list with two slots (`expanded` and `contracted`). Each slot is a list
with two dataframes (`paired` and `single_chain`)

    (list)
    └───expanded (list)
        └───paired (dataframe of all expanded αβ TCR pairs)
        └───single_chain (dataframe of all expanded single-chains)
    └───contracted (list)
        └───paired (dataframe of all contracted αβ TCR pairs)
        └───single_chain (dataframe of all contracted single-chains)

## Details

If you would like a scatterplot of this data, you may call
[`plot_sample_vs_sample()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_vs_sample.md)
with the same arguments for `log2fc_cutoff` and `sem_cutoff`.

To call expanded and contracted clonotypes from TIRTL-seq data, we
calculated mean frequency and standard error of the mean (SEM) for each
TCRβ chain over all wells. We call clones significantly expanded or
contracted between time points if there is a log2 fold-change log2FC \>
3 between average frequencies and the difference between average
frequencies exceeds 5 SEM intervals. This matches the analysis in
Pogorelyy & Kirk et al. (2025).

Note: For each TCR, we actually calculate two SEMs, one for each
timepoint/sample. To calculate 5 SEM intervals, we multiply each SEM by
2.5 and sum them. The `sem_cutoff` argument controls this value, which
is why the default is 2.5.

## References

Pogorelyy, M, Kirk, A, Adhikari, S et al. (2025). "TIRTL-seq: deep,
quantitative and affordable paired TCR repertoire sequencing." *Nature
Methods*, 23, 56–64.
[doi:10.1038/s41592-025-02907-9](https://doi.org/10.1038/s41592-025-02907-9)

## See also

Other longitudinal:
[`plot_clone_size_across_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clone_size_across_samples.md),
[`plot_sample_vs_sample()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_vs_sample.md)

## Examples

``` r
load_example_data() ## loads minimal SJTRC dataset into "SJTRC_minimal" object
#> Example data already loaded into object: 'SJTRC_minimal'
clones = get_expanded_clones(
  data1 = SJTRC_minimal$data$cd8_tp1_v2,
  data2 = SJTRC_minimal$data$cd8_tp2_v2,
  chain = "beta",
  remove_nonfunctional = FALSE,
  filter_pairs = FALSE)
## expanded clones
clones$expanded$paired ## paired data frame
#> # A tibble: 1,362 × 32
#>    sign  log2FC va         ja     cdr3a     cdr3b vb    jb    alpha_nuc beta_nuc
#>    <chr>  <dbl> <chr>      <chr>  <chr>     <chr> <chr> <chr> <chr>     <chr>   
#>  1 up      14.6 TRAV29/DV5 TRAJ42 CAASGGSQ… CSSQ… TRBV… TRBJ… TGTGCAGC… TGCAGCT…
#>  2 up      14.6 TRAV41     TRAJ42 CAISRGGM… CSSQ… TRBV… TRBJ… TGTGCTAT… TGCAGCT…
#>  3 up      13.8 TRAV30     TRAJ49 CGTESGNQ… CASS… TRBV… TRBJ… TGCGGCAC… TGTGCCA…
#>  4 up      13.8 TRAV23/DV6 TRAJ32 CAAR_ILIF CASS… TRBV… TRBJ… TGTGCAGC… TGTGCCA…
#>  5 up      13.7 TRAV12-2   TRAJ16 CAVLRGQK… CASY… TRBV… TRBJ… TGTGCCGT… TGTGCCA…
#>  6 up      13.7 TRAV8-2    TRAJ28 CVVRLVSR… CASY… TRBV… TRBJ… TGTGTTGT… TGTGCCA…
#>  7 up      13.7 TRAV12-2   TRAJ16 CAGLSDGQ… CASA… TRBV… TRBJ… TGTGCCGG… TGTGCCA…
#>  8 up      13.7 TRAV13-2   TRAJ39 CAQARPSN… CASA… TRBV… TRBJ… TGTGCGCA… TGTGCCA…
#>  9 up      13.6 TRAV13-1   TRAJ50 CAAIETSY… CSSG… TRBV… TRBJ… TGTGCAGC… TGCAGCT…
#> 10 up      13.6 TRAV23/DV6 TRAJ48 CAAPHL_N… CSSG… TRBV… TRBJ… TGTGCAGC… TGCAGCT…
#> # ℹ 1,352 more rows
#> # ℹ 22 more variables: comparison_chain <chr>, readFraction.x <dbl>,
#> #   readFraction.y <dbl>, sem.x <dbl>, sem.y <dbl>, n_wells.x <int>,
#> #   n_wells.y <int>, max_wells.x <int>, max_wells.y <int>, rank.x <int>,
#> #   rank.y <int>, readCount.x <int>, readCount.y <int>,
#> #   readCount_median.x <int>, readCount_median.y <int>, readCount_max.x <int>,
#> #   readCount_max.y <int>, alpha_beta <chr>, is_functional <lgl>, …
clones$expanded$single_chain ## single_chain
#> # A tibble: 910 × 23
#>    sign  log2FC is_paired beta_nuc              cdr3b vb    jb    readFraction.x
#>    <chr>  <dbl> <lgl>     <chr>                 <chr> <chr> <chr>          <dbl>
#>  1 up      14.6 TRUE      TGCAGCTCCCAGGAGGGAGG… CSSQ… TRBV… TRBJ…   0.000000203 
#>  2 up      13.8 TRUE      TGTGCCAGCAGCTTATTAGC… CASS… TRBV… TRBJ…   0.000000203 
#>  3 up      13.7 TRUE      TGTGCCAGCTATTCCAGCTC… CASY… TRBV… TRBJ…   0.0000000407
#>  4 up      13.7 TRUE      TGTGCCAGCGCGGTAAGCTC… CASA… TRBV… TRBJ…   0.000000214 
#>  5 up      13.6 TRUE      TGCAGCTCGGGGGGTAACAC… CSSG… TRBV… TRBJ…   0.000000112 
#>  6 up      13.1 TRUE      TGTGCCAGCAGCCCTGTAGC… CASS… TRBV… TRBJ…   0.000000163 
#>  7 up      13.0 TRUE      TGTGCCAGTAGTGTTGGCAG… CASS… TRBV… TRBJ…   0.0000000102
#>  8 up      12.8 TRUE      TGTGCCAGCAGCCTAAGGGA… CASS… TRBV… TRBJ…   0           
#>  9 up      12.8 TRUE      TGTGCCAGCAGCTTAACGGG… CASS… TRBV… TRBJ…   0.0000000305
#> 10 up      12.7 TRUE      TGTGCCACGTCGGCCAGTCA… CATS… TRBV… TRBJ…   0.0000000814
#> # ℹ 900 more rows
#> # ℹ 15 more variables: readFraction.y <dbl>, sem.x <dbl>, sem.y <dbl>,
#> #   n_wells.x <int>, n_wells.y <int>, max_wells.x <int>, max_wells.y <int>,
#> #   rank.x <int>, rank.y <int>, readCount.x <int>, readCount.y <int>,
#> #   readCount_median.x <int>, readCount_median.y <int>, readCount_max.x <int>,
#> #   readCount_max.y <int>
## contracted clones
clones$contracted$paired ## paired data frame
#> # A tibble: 189 × 32
#>    sign  log2FC va           ja     cdr3a   cdr3b vb    jb    alpha_nuc beta_nuc
#>    <chr>  <dbl> <chr>        <chr>  <chr>   <chr> <chr> <chr> <chr>     <chr>   
#>  1 down   -7.56 TRAV9-2      TRAJ9  CARNTG… CASS… TRBV… TRBJ… TGTGCTCG… TGTGCCA…
#>  2 down   -7.56 TRAV4        TRAJ8  CLVETG… CASS… TRBV… TRBJ… TGCCTCGT… TGTGCCA…
#>  3 down   -7.06 TRAV19       TRAJ9  CALILT… CASS… TRBV9 TRBJ… TGTGCTCT… TGTGCCA…
#>  4 down   -7.06 TRAV19       TRAJ9  CALILT… CASS… TRBV9 TRBJ… TGTGCTCT… TGTGCCA…
#>  5 down   -7.06 TRAV14/DV4   TRAJ21 CAIGGS… CASS… TRBV9 TRBJ… TGTGCAAT… TGTGCCA…
#>  6 down   -6.87 TRAV29/DV5   TRAJ29 CAARLL… CASS… TRBV… TRBJ… TGTGCAGC… TGTGCCA…
#>  7 down   -6.87 TRDV1        TRAJ48 CALGEL… CASS… TRBV… TRBJ… TGTGCTCT… TGTGCCA…
#>  8 down   -6.75 TRAV1-2      TRAJ42 CAVRDY… CASS… TRBV… TRBJ… TGTGCTGT… TGTGCCA…
#>  9 down   -6.75 TRAV2        TRAJ45 CAVCSI… CASS… TRBV… TRBJ… TGTGCTGT… TGTGCCA…
#> 10 down   -6.67 TRAV38-2/DV8 TRAJ49 CAYLTG… CASS… TRBV… TRBJ… TGTGCTTA… TGTGCCA…
#> # ℹ 179 more rows
#> # ℹ 22 more variables: comparison_chain <chr>, readFraction.x <dbl>,
#> #   readFraction.y <dbl>, sem.x <dbl>, sem.y <dbl>, n_wells.x <int>,
#> #   n_wells.y <int>, max_wells.x <int>, max_wells.y <int>, rank.x <int>,
#> #   rank.y <int>, readCount.x <int>, readCount.y <int>,
#> #   readCount_median.x <int>, readCount_median.y <int>, readCount_max.x <int>,
#> #   readCount_max.y <int>, alpha_beta <chr>, is_functional <lgl>, …
clones$contracted$single_chain ## single_chain
#> # A tibble: 119 × 23
#>    sign  log2FC is_paired beta_nuc              cdr3b vb    jb    readFraction.x
#>    <chr>  <dbl> <lgl>     <chr>                 <chr> <chr> <chr>          <dbl>
#>  1 down   -7.56 TRUE      TGTGCCAGCAGCTACCCGGG… CASS… TRBV… TRBJ…      0.000187 
#>  2 down   -7.06 TRUE      TGTGCCAGCAGCGTAGATTC… CASS… TRBV9 TRBJ…      0.000134 
#>  3 down   -6.87 TRUE      TGTGCCAGCAGCCCCGGGAC… CASS… TRBV… TRBJ…      0.000211 
#>  4 down   -6.75 TRUE      TGTGCCAGCAGCTTAGTCGG… CASS… TRBV… TRBJ…      0.000113 
#>  5 down   -6.67 TRUE      TGTGCCAGCAGTTTAGAGGG… CASS… TRBV… TRBJ…      0.000101 
#>  6 down   -6.55 TRUE      TGCAGCGCACTTAGCATCCG… CSAL… TRBV… TRBJ…      0.0000926
#>  7 down   -6.13 TRUE      TGCAGCGTTAACAGGCCTCA… CSVN… TRBV… TRBJ…      0.0000691
#>  8 down   -5.91 TRUE      TGCAGTGCGAACCGGGGGGA… CSAN… TRBV… TRBJ…      0.0000592
#>  9 down   -5.78 TRUE      TGTGCCAGCAGCCAAGCCGG… CASS… TRBV… TRBJ…      0.0000539
#> 10 down   -5.66 TRUE      TGCGCCAGCAGCCAAGATAG… CASS… TRBV… TRBJ…      0.000109 
#> # ℹ 109 more rows
#> # ℹ 15 more variables: readFraction.y <dbl>, sem.x <dbl>, sem.y <dbl>,
#> #   n_wells.x <int>, n_wells.y <int>, max_wells.x <int>, max_wells.y <int>,
#> #   rank.x <int>, rank.y <int>, readCount.x <int>, readCount.y <int>,
#> #   readCount_median.x <int>, readCount_median.y <int>, readCount_max.x <int>,
#> #   readCount_max.y <int>
```
