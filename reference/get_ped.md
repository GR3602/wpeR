# Organizes pedigree data

Offers an alternative to
[`get_colony()`](https://gr3602.github.io/wpeR/reference/get_colony.md)
function in cases where the pedigree was not reconstructed with
[COLONY](https://www.zsl.org/about-zsl/resources/software/colony)
software. It takes a pedigree dataframe and assigns sex to each
individual. The function also prepares data so that the output of the
function can be directly analyzed with
[`kinship2`](https://cran.r-project.org/package=kinship2),
[`pedtools`](https://cran.r-project.org/package=pedtools) or
[`FamAgg`](https://bioconductor.org/packages/FamAgg/) packages.

## Usage

``` r
get_ped(ped, sampledata, out = "FamAgg")
```

## Arguments

- ped:

  Data frame. Pedigree data frame with the most basic structure. Three
  columns corresponding to offspring, father and mother (see Details).
  Unknown parents should be represented by `NA` values.

- sampledata:

  Data frame. Metadata for all genetic samples that belong to the
  individuals included in pedigree reconstruction analysis. This data
  frame should adhere to the formatting and naming conventions outlined
  in the
  [`check_sampledata()`](https://gr3602.github.io/wpeR/reference/check_sampledata.md)
  documentation.

- out:

  Character string. For use with which package should the output be
  formatted? `kinship2` (out = "kinship2"), `pedtools` (out =
  "pedtools") or `FamAgg` (out = "FamAgg") or the created data.frame can
  be outputted as is (out = "table"). Defaults to "FamAgg"

## Value

A data frame describing a common pedigree structure. Each individual
included in pedigree represents one row. Columns describe individual
identifier code, identifier code for mother and father and sex of
individual. Column names and arrangement depends on selected output
(`out` parameter).

## Details

The custom pedigree specified through the `ped` parameter should mirror
the structure of a COLONY pedigree and share the same column names. It
should consist of three columns for each offspring: `OffspringID`,
`FatherID`, `MotherID`. When considering unknown parents they should be
represented by `NA` values.

## Examples

``` r
#example pedigree dataframe
ped <- data.frame(
  OffspringID = c(
    "M273P", "M20AM", "M2757", "M2ALK", "M2ETE", "M2EUJ", "MSV00E",
    "MSV018", "MSV05L", "MSV0M6", "MSV0T4", "MSV0T7", "MSV0TJ", "MSV0UL"
  ),
  FatherID = c(
    NA, NA, "M20AM", "M20AM", "M20AM", "M20AM", "M20AM",
    "M20AM", "M20AM", "M20AM", "M20AM", "M20AM", "M20AM", "M20AM"
  ),
  MotherID = c(
    NA, NA, "M273P", "M273P", "M273P", "M273P", "M273P",
    "M273P", "M273P", "M273P", "M273P", "M273P", "M273P", "M273P"
  )
)
#Get pedigree data in FamAgg format
get_ped(
    ped = ped,
    sampledata = wolf_samples
    )
#>        id father mother sex
#> 1   M273P   <NA>   <NA>   2
#> 2   M20AM   <NA>   <NA>   1
#> 3   M2757  M20AM  M273P   2
#> 4   M2ALK  M20AM  M273P   2
#> 5   M2ETE  M20AM  M273P   2
#> 6   M2EUJ  M20AM  M273P   2
#> 7  MSV00E  M20AM  M273P   1
#> 8  MSV018  M20AM  M273P   1
#> 9  MSV05L  M20AM  M273P   1
#> 10 MSV0M6  M20AM  M273P   1
#> 11 MSV0T4  M20AM  M273P   1
#> 12 MSV0T7  M20AM  M273P   1
#> 13 MSV0TJ  M20AM  M273P   2
#> 14 MSV0UL  M20AM  M273P   1
```
