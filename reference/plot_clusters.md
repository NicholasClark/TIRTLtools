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

Other tcr_similarity:
[`TCRdist()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist.md),
[`TCRdist_cpp()`](https://nicholasclark.github.io/TIRTLtools/reference/TCRdist_cpp.md),
[`cluster_tcrs()`](https://nicholasclark.github.io/TIRTLtools/reference/cluster_tcrs.md)

## Examples

``` r
# example code
# paired = load_tirtlseq("your_directory/")
# obj = cluster_tcrs(paired)
# plot_clusters(obj)
```
