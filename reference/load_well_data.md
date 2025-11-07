# Load well-level TIRTLseq data

**\[experimental\]**

## Usage

``` r
load_well_data(
  folder_path,
  wells = get_well_subset(1:16, 1:24),
  well_pos = 3,
  chain = c("both", "alpha", "beta"),
  nproc = 1L,
  columns = NULL,
  max_files = Inf
)
```

## See also

Other well:
[`filter_well_data()`](https://nicholasclark.github.io/TIRTLtools/reference/filter_well_data.md),
[`get_well_subset()`](https://nicholasclark.github.io/TIRTLtools/reference/get_well_subset.md),
[`get_wells_from_edges()`](https://nicholasclark.github.io/TIRTLtools/reference/get_wells_from_edges.md),
[`load_well_data_to_sparse_multi()`](https://nicholasclark.github.io/TIRTLtools/reference/load_well_data_to_sparse_multi.md)
