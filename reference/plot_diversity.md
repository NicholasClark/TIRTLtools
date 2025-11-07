# Bar plot of clonal diversity metrics

**\[experimental\]** The `plot_diversity()` plots the requested clonal
diversity metric

## Usage

``` r
plot_diversity(
  div,
  metric = get_all_div_metrics(),
  q = 2,
  percent = 90,
  group_col = NULL,
  label_col = "Sample",
  flip = FALSE,
  facet = FALSE,
  log_scale = FALSE,
  samples = NULL,
  return_data = FALSE,
  color_scheme = NULL,
  x_var = NULL
)
```

## Arguments

- div:

  a list created by the `diversity()` function with diversity metrics
  for each sample

- metric:

  the diversity metric to use (e.g. shannon, simpson, etc.)

- q:

  (optional) for 'renyi' and 'hill' metrics, the order q of the
  diversity index

- percent:

  (optional) for 'dXX' metric, the percentage 'XX' between 0 and 100

- group_col:

  (optional) if supplied, a column of the metadata that will be used to
  group samples

- label_col:

  (optional) labels for the samples

- flip:

  (optional) if TRUE, flip the x and y-axes (default is FALSE)

- facet:

  (optional) if TRUE, plot with separate facets, or sub-plots for each
  group (default is FALSE)

- log_scale:

  (optional) if TRUE, use log-scale for the y-axis (default is FALSE)

- samples:

  (optional) the samples to include in the plot (default is all samples)

- return_data:

  (optional) if TRUE, return the data used to make the plot (default is
  FALSE)

- color_scheme:

  (optional) a color scheme for the plot

- x_var:

  (optional) the variable to show on the x-axis, if other than
  "sample_id"

## Value

A list with two objects:

`$plot` - a ggplot object with the plot of the requested diversity
metric

`$data` - if return_data is TRUE, the data frame used to make the plot

## Details

This function can plot a variety of clonal diversity metrics for a
dataset (richness, Simpson diversity index, Shannon-Wiener index, etc.).
See
[`get_all_div_metrics()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_div_metrics.md)
for all available options. By default it return a barplot with one bar
for each sample in the dataset. If a grouping column (of the metadata)
is supplied, then samples will be grouped and bar heights will reflect
the average diversity metric across the group, with error bars showing
the standard deviation.

## See also

`diversity()`,
[`get_all_div_metrics()`](https://nicholasclark.github.io/TIRTLtools/reference/get_all_div_metrics.md)

Other plotting:
[`plot_clone_size_across_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clone_size_across_samples.md),
[`plot_clonotype_indices()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clonotype_indices.md),
[`plot_clusters()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clusters.md),
[`plot_n_reads()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_n_reads.md),
[`plot_num_partners()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_num_partners.md),
[`plot_paired()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired.md),
[`plot_paired_by_read_fraction_range()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_by_read_fraction_range.md),
[`plot_paired_vs_rank()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_paired_vs_rank.md),
[`plot_pairs_with_eachother()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_pairs_with_eachother.md),
[`plot_ranks()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_ranks.md),
[`plot_read_fraction_vs_pair_status()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_read_fraction_vs_pair_status.md),
[`plot_sample_overlap()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_overlap.md),
[`plot_sample_vs_sample()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_sample_vs_sample.md)

Other diversity:
[`calculate_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/calculate_diversity.md)

## Examples

``` r
# example code
# data = load_tirtlseq("your_directory/")
# div = diversity(data)
# plot_diversity(div, metric = "richness")
```
