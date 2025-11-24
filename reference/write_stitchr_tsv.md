# Write a tab-separated "thimble" text file for use with Stitchr

**\[experimental\]** This function takes a data frame created by
[`create_thimble_df()`](https://nicholasclark.github.io/TIRTLtools/reference/create_thimble_df.md)
and writes it to a .tsv file for use with Stitchr
(https://jamieheather.github.io/stitchr/).

## Usage

``` r
write_stitchr_tsv(df, output_name, output_folder = "")
```

## Arguments

- df:

  a data frame of properly formatted TCRs created by
  [`create_thimble_df()`](https://nicholasclark.github.io/TIRTLtools/reference/create_thimble_df.md)

- output_name:

  the output name for the .tsv file (".tsv" will be appended)

- output_folder:

  the folder to write the .tsv file to

## Value

this function returns NULL and writes a .tsv to the specified folder.

## See also

Other stitchr:
[`create_thimble_df()`](https://nicholasclark.github.io/TIRTLtools/reference/create_thimble_df.md)
