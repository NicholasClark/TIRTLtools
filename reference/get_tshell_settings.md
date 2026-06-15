# Get p-value and well threshold settings for T-SHELL

**\[experimental\]** This function returns a list with p-value and well
threshold settings for 384-well or 96-well T-SHELL.

## Usage

``` r
get_tshell_settings(format = c("auto", "384_well", "96_well"))
```

## Arguments

- format:

  a selection from "auto", "384_well", and "96_well".

## Value

A named list with slots "pval_thres_tshell" and "wij_thres_tshell". If
format is "384_well", the default p-value threshold for a TCR pair is
\<1e-10 and threshold for number of wells with both chains is \>2 wells.
If format is "96_well", these are \<1e-3 and \>3 wells respectively. If
format is "auto", the selection will be made based on the number of
wells passing QC in the sample. If \>= 150 wells pass QC, the "384_well"
settings will be used. Otherwise the "96_well" settings will be used.

## See also

Other pairing:
[`run_pairing()`](https://nicholasclark.github.io/TIRTLtools/reference/run_pairing.md)
