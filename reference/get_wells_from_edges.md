# Get well names from the top left and bottom right wells

**\[experimental\]** This function returns a vector of well names in a
rectangle containing the top left well and the bottom right well.

## Usage

``` r
get_wells_from_edges(
  top_left,
  bottom_right,
  return_type = c("wells", "rows_and_columns")
)
```

## Arguments

- top_left:

  the well in the top left of the rectangle (e.g. "A1")

- bottom_right:

  the well in the top left of the rectangle (e.g. "H12")

## Value

If return_type = "wells" (default), returns a vector of well names
covering the specified rows and columns (e.g. "A1", "A2", ... "H12"). If
return_type = "rows_and_columns", returns a list with two vectors

- rows - a vector of the rows (e.g. "A", "B", ..., "H")

- columns - a vector of the columns (e.g. 1, 2, 3, ..., 24)

## See also

Other well:
[`choose_pair_manual()`](https://nicholasclark.github.io/TIRTLtools/reference/choose_pair_manual.md),
[`get_well_subset()`](https://nicholasclark.github.io/TIRTLtools/reference/get_well_subset.md),
[`load_well_counts_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/load_well_counts_binary.md),
[`plot_tshell()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_tshell.md),
[`write_well_data_to_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/write_well_data_to_binary.md)

## Examples

``` r
get_wells_from_edges("A1", "H12") ## top left half of plate
#> Error in get_wells_from_edges("A1", "H12"): could not find function "get_wells_from_edges"

get_wells_from_edges("A13", "P24") ## right half of plate
#> Error in get_wells_from_edges("A13", "P24"): could not find function "get_wells_from_edges"

get_wells_from_edges("A13", "P24", return_type = "rows_and_columns")
#> Error in get_wells_from_edges("A13", "P24", return_type = "rows_and_columns"): could not find function "get_wells_from_edges"
```
