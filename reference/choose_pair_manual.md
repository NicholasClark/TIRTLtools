# Choose a partner manually for an input chain using T-SHELL

**\[experimental\]**

This function takes in an object created by
[`plot_tshell()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_tshell.md)
and allows a user to choose the partner for the input chain used to
create that object.

## Usage

``` r
choose_pair_manual(res, rank)
```

## Arguments

- res:

  an object created by
  [`plot_tshell()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_tshell.md)

- rank:

  the rank of the potential partner chain to choose

## Value

A one-row data frame in the style of the "paired" data frame from
pairing output with the input chain and its chosen partner.

## See also

Other well:
[`get_well_subset()`](https://nicholasclark.github.io/TIRTLtools/reference/get_well_subset.md),
[`get_wells_from_edges()`](https://nicholasclark.github.io/TIRTLtools/reference/get_wells_from_edges.md),
[`load_well_counts_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/load_well_counts_binary.md),
[`plot_tshell()`](https://nicholasclark.github.io/TIRTLtools/reference/plot_tshell.md),
[`write_well_data_to_binary()`](https://nicholasclark.github.io/TIRTLtools/reference/write_well_data_to_binary.md)
