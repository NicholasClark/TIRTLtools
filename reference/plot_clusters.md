# Plot clusters of similar TCRs

**\[experimental\]** `plot_clusters()` plots TCRs that were previously
clustered by the
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md)
function.

## Usage

``` r
plot_clusters(
  obj,
  n_clusters = 10,
  seed = 1234,
  annotation_cols = c("cluster", "source"),
  color_col = "cluster"
)
```

## Arguments

- obj:

  an object returned by the
  [`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md)
  function

- n_clusters:

  the number of clusters to plot (default 10)

- seed:

  a number to use as the seed. Using the same number across multiple
  runs allows for reproducible results.

- annotation_cols:

  columns to use for heatmap annotations (default is the cluster and the
  source, i.e. your data or vdj-db)

- color_col:

  a column to use to color nodes in the network visualization (default
  is "cluster")

## Value

Currently returns a list with three items: `$umap` - a umap
visualization of the TCRs from the top clusters

`$graph` - a graph/network visualization of TCRs from the top clusters

`$heatmap` - a heatmap visualization of TCRs from the top clusters

## Details

The function currently returns a list with a UMAP plot, a graph/network
plot, and a heatmap.

## See also

[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md),
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md)

Other plotting:
[`plot_clone_size_across_samples()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clone_size_across_samples.md),
[`plot_clonotype_indices()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_clonotype_indices.md),
[`plot_diversity()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_diversity.md),
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

## Examples

``` r
# example code
# paired = load_tirtlseq("your_directory/")
# obj = cluster_tcrs(paired)
# plot_clusters(obj)
```
