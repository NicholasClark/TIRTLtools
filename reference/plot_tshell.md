# Calculate and plot read fraction correlation across wells (T-SHELL)

**\[experimental\]**

This function takes in an object with well-level read counts for a
sample, loaded by
[`load_well_counts_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/load_well_counts_binary.md)
along with a CDR3 sequence of interest for TCRalpha or TCRbeta. It
calculates the correlation of the input sequence with all of its
potential partner chains and returns the top N chains with the highest
correlation/T-SHELL value.

## Usage

``` r
plot_tshell(
  well_data,
  nuc_seq,
  chain,
  n_plot = 9,
  plot = TRUE,
  plot_log_scale = FALSE,
  interactive = FALSE,
  normalize_counts = TRUE,
  loss_frac_cutoff = 0.5,
  calc_fisher_pval_all = FALSE
)
```

## Arguments

- well_data:

  well-level read count data loaded by
  [`load_well_counts_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/load_well_counts_binary.md)

- nuc_seq:

  the CDR3 nucleotide sequence of the chain of interest

- chain:

  the chain of the input nucleotide sequence, either "alpha" or "beta"

- n_plot:

  the number of scatter plots to return (one for each of top N potential
  partners, default is 9)

- plot:

  whether to plot the data or only return data frames (default is TRUE)

- plot_log_scale:

  whether to plot log(readFraction) instead of readFraction (default is
  FALSE)

- interactive:

  whether to return interactive plots (default is FALSE)

- normalize_counts:

  whether to normalize the counts by the total count for each well
  (default is TRUE)

- loss_frac_cutoff:

  the cutoff to use for "loss fraction" of potential partners for plots.
  Must be between 0 and 1. 1 or higher will not filter any clones, lower
  values will apply more strict filtering. (default is 0.5).

- calc_fisher_pval_all:

  whether to calculate the Fisher exact test p-value of well overlap
  between the input chain and all of the filtered clones. Default is
  FALSE because this is computationally demanding.

## Value

A list of output with the following elements:

- data (list)

- df_top_n - a data frame with results for the top N potential partner
  chains

- df_all - a data frame with results for all potential partners passing
  the `loss_frac_cutoff`

- input_meta - a one-row data frame with the metadata for the input
  chain

- plots (list)

- scatter - a scatter plot of read fractions for the input chain vs. the
  top N potential partner chains

- r_vs_p - a point plot of correlation (r) vs. -log10(raw p-value)

- r_vs_p_adj - a point plot of correlation (r) vs. -log10(FDR adjusted
  p-value)

- rank_vs_p_adj - a Manhattan plot of chain rank vs. -log10(FDR adjusted
  p-value)

- call - a list with all arguments given to the function

## See also

Other well:
[`choose_pair_manual()`](https://nicholasclark.github.io/TIRTLtools/reference/choose_pair_manual.md),
[`get_well_subset()`](https://nicholasclark.github.io/TIRTLtools/reference/get_well_subset.md),
[`get_wells_from_edges()`](https://nicholasclark.github.io/TIRTLtools/reference/get_wells_from_edges.md),
[`load_well_counts_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/load_well_counts_binary.md),
[`write_well_data_to_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/write_well_data_to_binary.md)
