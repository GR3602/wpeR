# Organizes COLONY output

Extends `BestConfig_Ordered` output from
[COLONY](https://www.zsl.org/about-zsl/resources/software/colony)
pedigree reconstruction software with additional data about individuals
included in pedigree. The function adds missing parents to
`OffspringID`, assigns sex to each individual included in `OffspringID`
and adds the computed probabilities of paternity and maternity
assignments (probability of assignments is visible only if the `out`
parameter is set to `"table"`). The function also prepares data so that
the output of the function can be directly analyzed with
[`kinship2`](https://cran.r-project.org/package=kinship2),
[`pedtools`](https://cran.r-project.org/package=pedtools) or
[`FamAgg`](https://bioconductor.org/packages/FamAgg/) packages.

## Usage

``` r
get_colony(
  colony_project_path,
  sampledata,
  rm_obsolete_parents = TRUE,
  out = "FamAgg"
)
```

## Arguments

- colony_project_path:

  Character string. Path to the folder where COLONY output files are
  saved. Has to include file path and project name (see Details).

- sampledata:

  Data frame. Metadata for all genetic samples that belong to the
  individuals included in pedigree reconstruction analysis. This data
  frame should adhere to the formatting and naming conventions outlined
  in the
  [`check_sampledata()`](https://gr3602.github.io/wpeR/reference/check_sampledata.md)
  documentation.

- rm_obsolete_parents:

  Logical. Should unknown parents be removed from output. Applies just
  to offspring for which both parents are unknown. Defaults to `TRUE`.

- out:

  Character string. For use with which package should the output be
  formatted? `kinship2` (out = "kinship2"), `pedtools` (out =
  "pedtools"), `FamAgg` (out = "FamAgg") or the created data.frame can
  be outputted as is (out = "table"). Defaults to "FamAgg"

## Value

A data frame describing a common pedigree structure. Each individual
included in pedigree represents one row. Columns describe individual
identifier code, identifier code for mother and father, sex and family
of individual. Column names and arrangement depends on selected output
(`out` parameter).

## Details

COLONY output tables needed for this function (`.BestConfig_Ordered`,
`.Maternity` and `.Paternity`) are read directly from the colony output
folder and do not need to be imported into R session. The path to the
outputs is defined with `colony_project_path` parameter. When defining
`colony_project_path` the user needs to define a complete path to the
directory where colony outputs are stored and also the file name (file
name of COLONY outputs equals the project name  
eg. /path/to/the/COLONY/output/folder/COLONY_project_name).

## Examples

``` r
# Define the path to COLONY output
path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")

# Get pedigree data in FamAgg format
get_colony(
    colony_project_path = path,
    sampledata = wolf_samples
    )
#>    ClusterIndex     id father mother sex
#> 1             1  M2AM8   <NA>   <NA>   1
#> 2             1  M273P   <NA>   <NA>   2
#> 3             1  M20AM   <NA>   <NA>   1
#> 4             1  M28TU   <NA>   <NA>   2
#> 5             1  M228J   <NA>   <NA>   1
#> 6             1  M200F   <NA>   <NA>   2
#> 7             1  M10XC  M228J  M200F   1
#> 8             1  M1YP0  M228J  M200F   2
#> 9             1  M220J  M228J  M200F   1
#> 10            1  M22AM  M228J  M200F   1
#> 11            1  M2772  M228J  M200F   1
#> 12            1  M28LU  M228J  M200F   2
#> 13            1  M2AXE  M228J  M200F   2
#> 14            1  M2C1T  M228J  M200F   2
#> 15            1 MSV10T MSV00E  M28LU   2
#> 16            1 MSV1C0 MSV00E  M28LU   1
#> 17            1 MSV1EX MSV00E  M28LU   1
#> 18            1 MSV1F5 MSV00E  M28LU   1
#> 19            1 MSV1F8 MSV00E  M28LU   2
#> 20            1 MSV1FJ MSV00E  M28LU   1
#> 21            1 MSV1FL MSV00E  M28LU   1
#> 22            1 MSV1FT MSV00E  M28LU   1
#> 23            1 MSV1TM MSV00E  M28LU   2
#> 24            1 MSV01X  M2772  M28TU   1
#> 25            1 MSV055  M2772  M28TU   2
#> 26            1 MSV0AL  M2772  M28TU   1
#> 27            1 MSV0TA  M2772  M28TU   1
#> 28            1 MSV0XT  M2772  M28TU   1
#> 29            1 MSV17U  M2772  M28TU   1
#> 30            1 MSV180  M2772  M28TU   2
#> 31            1 MSV18C  M2772  M28TU   2
#> 32            1 MSV1KT  M2772  M28TU   2
#> 33            1  M1J47  M2AM8  M200F   1
#> 34            1  M221C  M2AM8  M200F   2
#> 35            1  M2C8Y  M2AM8  M200F   2
#> 36            1  M2F1L  M2AM8  M200F   2
#> 37            1 MSV02F  M2AM8  M200F   1
#> 38            1 MSV02L  M2AM8  M200F   1
#> 39            1 MSV0CK  M2AM8  M200F   2
#> 40            1 MSV0FK  M2AM8  M200F   1
#> 41            1 MSV0H5  M2AM8  M200F   2
#> 42            1 MSV0P7  M2AM8  M200F   1
#> 43            1 MSV0UP  M2AM8  M200F   2
#> 44            1 MSV0UT  M2AM8  M200F   1
#> 45            1 MSV0UU  M2AM8  M200F   2
#> 46            1 MSV16T  M2AM8  M200F   2
#> 47            1 MSV16U  M2AM8  M200F   1
#> 48            1 MSV170  M2AM8  M200F   1
#> 49            1 MSV177  M2AM8  M200F   2
#> 50            1 MSV1FE  M2AM8  M200F   1
#> 51            1  M2757  M20AM  M273P   2
#> 52            1  M2ALK  M20AM  M273P   2
#> 53            1  M2ETE  M20AM  M273P   2
#> 54            1  M2EUJ  M20AM  M273P   2
#> 55            1 MSV00E  M20AM  M273P   1
#> 56            1 MSV018  M20AM  M273P   1
#> 57            1 MSV05L  M20AM  M273P   1
#> 58            1 MSV0M6  M20AM  M273P   1
#> 59            1 MSV0T4  M20AM  M273P   1
#> 60            1 MSV0T7  M20AM  M273P   1
#> 61            1 MSV0TJ  M20AM  M273P   2
#> 62            1 MSV0UL  M20AM  M273P   1
#> 63            1 MSV0X4  M20AM  M273P   1
#> 64            1 MSV17F  M20AM  M273P   2
#> 65            1 MSV1MH  M20AM  M273P   2


```
