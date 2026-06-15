# Table of known TCRs from VDJ-db

Table of known TCRs from VDJ-db

## Usage

``` r
vdj_db
```

## Format

A data frame with metadata for known TCRs from VDJ-db
(<http://vdjdb.com/>)

The current table contains 39,042 TCRs annotated to 30 different antigen
species and 750 different epitopes. The current data is from the VDJ-db
2025 Autumn release
(<https://github.com/antigenomics/vdjdb-db/releases>). Only human TCRs
are included and self-antigens were excluded due to excessive
interconnection when clustering.

## Examples

``` r
vdj_db
#> # A tibble: 39,042 × 34
#>    cdr3a       va    ja    cdr3b vb    db    jb    species mhc.a mhc.b mhc.class
#>    <chr>       <chr> <chr> <chr> <chr> <chr> <chr> <chr>   <chr> <chr> <chr>    
#>  1 CILPLAGGTS… TRAV… TRAJ… CASS… TRBV… ""    TRBJ… HomoSa… HLA-… B2M   MHCI     
#>  2 CILPLAGGTS… TRAV… TRAJ… CASS… TRBV… ""    TRBJ… HomoSa… HLA-… B2M   MHCI     
#>  3 CAMEGNSGYS… TRAV… TRAJ… CASS… TRBV… ""    TRBJ… HomoSa… HLA-… B2M   MHCI     
#>  4 CALYNTDKLIF TRAV… TRAJ… CASS… TRBV… ""    TRBJ… HomoSa… HLA-… B2M   MHCI     
#>  5 CGTEEPNDYK… TRAV… TRAJ… CASS… TRBV… ""    TRBJ… HomoSa… HLA-… B2M   MHCI     
#>  6 CAPAGASGYS… TRAV… TRAJ… CASS… TRBV… ""    TRBJ… HomoSa… HLA-… B2M   MHCI     
#>  7 CAPAGPSGYS… TRAV… TRAJ… CASS… TRBV… ""    TRBJ… HomoSa… HLA-… B2M   MHCI     
#>  8 CAVDPWGNQF… TRAV… TRAJ… CATQ… TRBV… ""    TRBJ… HomoSa… HLA-… B2M   MHCI     
#>  9 CAVKDTDKLIF TRAV… TRAJ… CATS… TRBV… ""    TRBJ… HomoSa… HLA-… B2M   MHCI     
#> 10 CALIYYGGSQ… TRAV… TRAJ… CATF… TRBV… ""    TRBJ… HomoSa… HLA-… B2M   MHCI     
#> # ℹ 39,032 more rows
#> # ℹ 23 more variables: antigen.epitope <chr>, antigen.gene <chr>,
#> #   antigen.species <chr>, reference.id <chr>, method.identification <chr>,
#> #   method.frequency <chr>, method.singlecell <chr>, method.sequencing <chr>,
#> #   method.verification <chr>, meta.study.id <chr>, meta.cell.subset <chr>,
#> #   meta.subject.cohort <chr>, meta.subject.id <chr>, meta.replica.id <chr>,
#> #   meta.clone.id <chr>, meta.epitope.id <chr>, meta.tissue <chr>, …
```
