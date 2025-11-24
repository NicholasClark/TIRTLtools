# Get well names from numerical rows and columns

**\[experimental\]** This function returns a vector of well names from
the specified rows and columns.

## Usage

``` r
get_well_subset(row_range = 1:16, col_range = 1:24)
```

## Arguments

- row_range:

  a vector of row indices, from 1 (A) to 16 (P)

- col_range:

  a vector of column indices, from 1 to 24

## Value

A vector of well names covering the specified rows and columns (e.g.
"A1", "A2", ... "P24")

## Examples

``` r
get_well_subset(1:16, 1:12) ## left half of plate
#>   [1] "A1"  "A2"  "A3"  "A4"  "A5"  "A6"  "A7"  "A8"  "A9"  "A10" "A11" "A12"
#>  [13] "B1"  "B2"  "B3"  "B4"  "B5"  "B6"  "B7"  "B8"  "B9"  "B10" "B11" "B12"
#>  [25] "C1"  "C2"  "C3"  "C4"  "C5"  "C6"  "C7"  "C8"  "C9"  "C10" "C11" "C12"
#>  [37] "D1"  "D2"  "D3"  "D4"  "D5"  "D6"  "D7"  "D8"  "D9"  "D10" "D11" "D12"
#>  [49] "E1"  "E2"  "E3"  "E4"  "E5"  "E6"  "E7"  "E8"  "E9"  "E10" "E11" "E12"
#>  [61] "F1"  "F2"  "F3"  "F4"  "F5"  "F6"  "F7"  "F8"  "F9"  "F10" "F11" "F12"
#>  [73] "G1"  "G2"  "G3"  "G4"  "G5"  "G6"  "G7"  "G8"  "G9"  "G10" "G11" "G12"
#>  [85] "H1"  "H2"  "H3"  "H4"  "H5"  "H6"  "H7"  "H8"  "H9"  "H10" "H11" "H12"
#>  [97] "I1"  "I2"  "I3"  "I4"  "I5"  "I6"  "I7"  "I8"  "I9"  "I10" "I11" "I12"
#> [109] "J1"  "J2"  "J3"  "J4"  "J5"  "J6"  "J7"  "J8"  "J9"  "J10" "J11" "J12"
#> [121] "K1"  "K2"  "K3"  "K4"  "K5"  "K6"  "K7"  "K8"  "K9"  "K10" "K11" "K12"
#> [133] "L1"  "L2"  "L3"  "L4"  "L5"  "L6"  "L7"  "L8"  "L9"  "L10" "L11" "L12"
#> [145] "M1"  "M2"  "M3"  "M4"  "M5"  "M6"  "M7"  "M8"  "M9"  "M10" "M11" "M12"
#> [157] "N1"  "N2"  "N3"  "N4"  "N5"  "N6"  "N7"  "N8"  "N9"  "N10" "N11" "N12"
#> [169] "O1"  "O2"  "O3"  "O4"  "O5"  "O6"  "O7"  "O8"  "O9"  "O10" "O11" "O12"
#> [181] "P1"  "P2"  "P3"  "P4"  "P5"  "P6"  "P7"  "P8"  "P9"  "P10" "P11" "P12"

get_well_subset(1:16, 13:24) ## right half of plate
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
```
