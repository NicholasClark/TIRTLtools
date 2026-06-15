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
get_wells_from_edges("A13", "P24") ## right half of plate
#>   [1] "A13" "A14" "A15" "A16" "A17" "A18" "A19" "A20" "A21" "A22" "A23" "A24"
#>  [13] "B13" "B14" "B15" "B16" "B17" "B18" "B19" "B20" "B21" "B22" "B23" "B24"
#>  [25] "C13" "C14" "C15" "C16" "C17" "C18" "C19" "C20" "C21" "C22" "C23" "C24"
#>  [37] "D13" "D14" "D15" "D16" "D17" "D18" "D19" "D20" "D21" "D22" "D23" "D24"
#>  [49] "E13" "E14" "E15" "E16" "E17" "E18" "E19" "E20" "E21" "E22" "E23" "E24"
#>  [61] "F13" "F14" "F15" "F16" "F17" "F18" "F19" "F20" "F21" "F22" "F23" "F24"
#>  [73] "G13" "G14" "G15" "G16" "G17" "G18" "G19" "G20" "G21" "G22" "G23" "G24"
#>  [85] "H13" "H14" "H15" "H16" "H17" "H18" "H19" "H20" "H21" "H22" "H23" "H24"
#>  [97] "I13" "I14" "I15" "I16" "I17" "I18" "I19" "I20" "I21" "I22" "I23" "I24"
#> [109] "J13" "J14" "J15" "J16" "J17" "J18" "J19" "J20" "J21" "J22" "J23" "J24"
#> [121] "K13" "K14" "K15" "K16" "K17" "K18" "K19" "K20" "K21" "K22" "K23" "K24"
#> [133] "L13" "L14" "L15" "L16" "L17" "L18" "L19" "L20" "L21" "L22" "L23" "L24"
#> [145] "M13" "M14" "M15" "M16" "M17" "M18" "M19" "M20" "M21" "M22" "M23" "M24"
#> [157] "N13" "N14" "N15" "N16" "N17" "N18" "N19" "N20" "N21" "N22" "N23" "N24"
#> [169] "O13" "O14" "O15" "O16" "O17" "O18" "O19" "O20" "O21" "O22" "O23" "O24"
#> [181] "P13" "P14" "P15" "P16" "P17" "P18" "P19" "P20" "P21" "P22" "P23" "P24"

get_wells_from_edges("A1", "H12") ## top left half of plate
#>  [1] "A1"  "A2"  "A3"  "A4"  "A5"  "A6"  "A7"  "A8"  "A9"  "A10" "A11" "A12"
#> [13] "B1"  "B2"  "B3"  "B4"  "B5"  "B6"  "B7"  "B8"  "B9"  "B10" "B11" "B12"
#> [25] "C1"  "C2"  "C3"  "C4"  "C5"  "C6"  "C7"  "C8"  "C9"  "C10" "C11" "C12"
#> [37] "D1"  "D2"  "D3"  "D4"  "D5"  "D6"  "D7"  "D8"  "D9"  "D10" "D11" "D12"
#> [49] "E1"  "E2"  "E3"  "E4"  "E5"  "E6"  "E7"  "E8"  "E9"  "E10" "E11" "E12"
#> [61] "F1"  "F2"  "F3"  "F4"  "F5"  "F6"  "F7"  "F8"  "F9"  "F10" "F11" "F12"
#> [73] "G1"  "G2"  "G3"  "G4"  "G5"  "G6"  "G7"  "G8"  "G9"  "G10" "G11" "G12"
#> [85] "H1"  "H2"  "H3"  "H4"  "H5"  "H6"  "H7"  "H8"  "H9"  "H10" "H11" "H12"

get_wells_from_edges("A13", "P24", return_type = "rows_and_columns")
#> $rows
#>  [1] "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "O" "P"
#> 
#> $columns
#>  [1] 13 14 15 16 17 18 19 20 21 22 23 24
#> 
```
