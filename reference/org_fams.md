# Organize animals into families and expand pedigree data

Takes pedigree data from
[`get_colony()`](https://gr3602.github.io/wpeR/reference/get_colony.md)
or [`get_ped()`](https://gr3602.github.io/wpeR/reference/get_ped.md)
function and groups animals into families. It also expands the pedigree
data by adding information about the family that each individual was
born in and the family in which the individual is the reproductive
animal.

## Usage

``` r
org_fams(ped, sampledata, output = "both")
```

## Arguments

- ped:

  Data frame. `FamAgg` output of
  [`get_colony()`](https://gr3602.github.io/wpeR/reference/get_colony.md)
  or [`get_ped()`](https://gr3602.github.io/wpeR/reference/get_ped.md)
  function. With `rm_obsolete_parents` parameter set to `TRUE`.

- sampledata:

  Data frame. Metadata for all genetic samples that belong to the
  individuals included in pedigree reconstruction analysis. This data
  frame should adhere to the formatting and naming conventions outlined
  in the
  [`check_sampledata()`](https://gr3602.github.io/wpeR/reference/check_sampledata.md)
  documentation.

- output:

  Character string. Determines the format of the output. Options are:
  "ped": returns an extended pedigree data frame. "fams": returns a
  table of all families present in the pedigree. "both": returns a list
  with two data frames: "ped" and "fams". (Default)

## Value

Depending on the `output` parameter, the function returns either a data
frame (`ped` or `fams`) or a list containing both data frames (`ped` and
`fams`).

- `ped` data frame. An extended version of the pedigree data from
  [`get_colony()`](https://gr3602.github.io/wpeR/reference/get_colony.md)/[`get_ped()`](https://gr3602.github.io/wpeR/reference/get_ped.md).
  In addition to common pedigree information (individual, mother,
  father, sex, family), `ped` includes columns for:

  - `parents`: Identifier codes of both parents separated with `_`.

  - `FamID`: Numeric identifier for the family to which the individual
    belongs (see `fams` below).

  - `FirstSeen`: Date of first sample of individual.

  - `LastSeen`: Date of last sample of individual.

  - `IsDead`: Logical value (`TRUE/FALSE`) that identifies if the
    individual is dead.

  - `DadHSgroup`: Identifier of paternal half-sib group (see Details).

  - `MomHSgroup`: Identifier of maternal half-sib group (see Details).

  - `hsGroup`: Numeric value indicating if the individual is part of a
    half-sib group (see Details).

- `fams` data frame includes information on families that individuals in
  the pedigree belong to. The families are described by:

  - `parents`: Identifier codes of both parents separated with `_`.

  - `father`: Identifier code of the father.

  - `mother`: Identifier code of the mother.

  - `FamID`: Numeric identifier for the family.

  - `famStart`: Date when the first sample of one of the offspring from
    this family was collected (see Details).

  - `famEnd`: Date when the last sample of mother or father of this
    family was collected (see Details).

  - `FamDead`: Logical value (`TRUE/FALSE`) indicating if the family no
    longer exists.

  - `DadHSgroup`: Identifier connecting families that share the same
    father.

  - `MomHSgroup`: Identifier connecting families that share the same
    mother.

  - `hsGroup`: Numeric value connecting families that share one of the
    parents.

## Details

**Families and Half-sib Groups** The result of `org_fams()` function
introduces us to two important concepts within the context of this
package: family and half-sib group. A family in the output of this
function is defined as a group of animals where at least one parent and
at least one offspring is known. A half-sib group refers to a group of
half-siblings, either maternally or paternally related. In the function
output the `DadHSgroup` groups paternal half-siblings and `MomHSgroup`
maternal half-siblings.

**Lineage Isolation for Visualization** The `DadHSgroup` and
`MomHSgroup` identifiers are usefull for managing pedigree visualization
in long-lived or polygamous species. By using `DadHSgroup` or
`MomHSgroup`, users can identify all `FamID`s or individuals associated
with a specific parent. This vector of IDs or individuals of interest
can then be passed to the `plot_fams` or `plot_indivs` argument in
[`plot_table()`](https://gr3602.github.io/wpeR/reference/plot_table.md)
to isolate and visualize specific paternal or maternal lineages,
preventing visual clutter in
[`ped_satplot()`](https://gr3602.github.io/wpeR/reference/ped_satplot.md).

**Temporal Estimation (famStart and famEnd)** The `fams` output
dataframe contains `famStart` and `famEnd` columns, which estimate a
time window for the family based solely on sample collection dates
provided in `sampledata`. `famStart` marks the date of the earliest
sample collected from any offspring belonging to that family. `famEnd`
indicates the date of the latest sample collected from either the mother
or the father of that family. It is important to recognize that this
method relies on observation (sampling) times. Consequently, `famEnd`
(last parental sample date) can precede `famStart` (first offspring
sample date), creating a biologically impossible sequence and a negative
calculated family timespan. Users should interpret the interval between
`famStart` and `famEnd` with this understanding.

## Examples

``` r

# Prepare the data for usage with org_fams() function.
# Get animal timespan data using the anim_timespan() function.
animal_ts <- anim_timespan(
  wolf_samples$AnimalRef,
  wolf_samples$Date,
  wolf_samples$IsMortality
)
# Add animal timespan to the sampledata
sampledata <- merge(wolf_samples, animal_ts, by.x = "AnimalRef", by.y = "ID", all.x = TRUE)
# Define the path to the pedigree data file.
path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
# Retrieve the pedigree data from the get_colony function.
ped_colony <- get_colony(path, sampledata, rm_obsolete_parents = TRUE, out = "FamAgg")

# Run the function
# Organize families and expand pedigree data using the org_fams function.
org_fams(
    ped = ped_colony,
    sampledata = sampledata
    )
#> $ped
#>    ClusterIndex     id father mother sex      parents FamID  FirstSeen
#> 1             1  M2AM8   <NA>   <NA>   1         <NA>     0 2017-04-07
#> 2             1  M273P   <NA>   <NA>   2         <NA>     0 2018-01-02
#> 3             1  M20AM   <NA>   <NA>   1         <NA>     0 2016-08-29
#> 4             1  M28TU   <NA>   <NA>   2         <NA>     0 2017-12-18
#> 5             1  M228J   <NA>   <NA>   1         <NA>     0 2016-09-30
#> 6             1  M200F   <NA>   <NA>   2         <NA>     0 2015-07-27
#> 7             1  M10XC  M228J  M200F   1  M228J_M200F     1 2017-11-16
#> 8             1  M1YP0  M228J  M200F   2  M228J_M200F     1 2017-01-25
#> 9             1  M220J  M228J  M200F   1  M228J_M200F     1 2017-11-10
#> 10            1  M22AM  M228J  M200F   1  M228J_M200F     1 2017-01-26
#> 11            1  M2772  M228J  M200F   1  M228J_M200F     1 2017-11-12
#> 12            1  M28LU  M228J  M200F   2  M228J_M200F     1 2017-09-18
#> 13            1  M2AXE  M228J  M200F   2  M228J_M200F     1 2017-12-16
#> 14            1  M2C1T  M228J  M200F   2  M228J_M200F     1 2017-10-31
#> 15            1 MSV10T MSV00E  M28LU   2 MSV00E_M28LU     2 2020-08-16
#> 16            1 MSV1C0 MSV00E  M28LU   1 MSV00E_M28LU     2 2021-02-16
#> 17            1 MSV1EX MSV00E  M28LU   1 MSV00E_M28LU     2 2021-01-18
#> 18            1 MSV1F5 MSV00E  M28LU   1 MSV00E_M28LU     2 2021-01-23
#> 19            1 MSV1F8 MSV00E  M28LU   2 MSV00E_M28LU     2 2021-02-04
#> 20            1 MSV1FJ MSV00E  M28LU   1 MSV00E_M28LU     2 2021-02-16
#> 21            1 MSV1FL MSV00E  M28LU   1 MSV00E_M28LU     2 2021-02-17
#> 22            1 MSV1FT MSV00E  M28LU   1 MSV00E_M28LU     2 2021-01-18
#> 23            1 MSV1TM MSV00E  M28LU   2 MSV00E_M28LU     2 2020-11-15
#> 24            1 MSV01X  M2772  M28TU   1  M2772_M28TU     3 2019-08-28
#> 25            1 MSV055  M2772  M28TU   2  M2772_M28TU     3 2019-12-11
#> 26            1 MSV0AL  M2772  M28TU   1  M2772_M28TU     3 2019-12-11
#> 27            1 MSV0TA  M2772  M28TU   1  M2772_M28TU     3 2020-01-12
#> 28            1 MSV0XT  M2772  M28TU   1  M2772_M28TU     3 2020-10-08
#> 29            1 MSV17U  M2772  M28TU   1  M2772_M28TU     3 2020-09-29
#> 30            1 MSV180  M2772  M28TU   2  M2772_M28TU     3 2020-09-29
#> 31            1 MSV18C  M2772  M28TU   2  M2772_M28TU     3 2020-11-06
#> 32            1 MSV1KT  M2772  M28TU   2  M2772_M28TU     3 2020-10-06
#> 33            1  M1J47  M2AM8  M200F   1  M2AM8_M200F     4 2019-08-20
#> 34            1  M221C  M2AM8  M200F   2  M2AM8_M200F     4 2018-10-29
#> 35            1  M2C8Y  M2AM8  M200F   2  M2AM8_M200F     4 2018-10-30
#> 36            1  M2F1L  M2AM8  M200F   2  M2AM8_M200F     4 2019-04-07
#> 37            1 MSV02F  M2AM8  M200F   1  M2AM8_M200F     4 2019-10-07
#> 38            1 MSV02L  M2AM8  M200F   1  M2AM8_M200F     4 2019-09-15
#> 39            1 MSV0CK  M2AM8  M200F   2  M2AM8_M200F     4 2019-08-05
#> 40            1 MSV0FK  M2AM8  M200F   1  M2AM8_M200F     4 2020-01-18
#> 41            1 MSV0H5  M2AM8  M200F   2  M2AM8_M200F     4 2020-04-10
#> 42            1 MSV0P7  M2AM8  M200F   1  M2AM8_M200F     4 2019-11-23
#> 43            1 MSV0UP  M2AM8  M200F   2  M2AM8_M200F     4 2020-06-10
#> 44            1 MSV0UT  M2AM8  M200F   1  M2AM8_M200F     4 2020-06-10
#> 45            1 MSV0UU  M2AM8  M200F   2  M2AM8_M200F     4 2020-06-10
#> 46            1 MSV16T  M2AM8  M200F   2  M2AM8_M200F     4 2021-02-02
#> 47            1 MSV16U  M2AM8  M200F   1  M2AM8_M200F     4 2021-02-02
#> 48            1 MSV170  M2AM8  M200F   1  M2AM8_M200F     4 2021-02-02
#> 49            1 MSV177  M2AM8  M200F   2  M2AM8_M200F     4 2020-04-22
#> 50            1 MSV1FE  M2AM8  M200F   1  M2AM8_M200F     4 2021-02-10
#> 51            1  M2757  M20AM  M273P   2  M20AM_M273P     5 2018-01-05
#> 52            1  M2ALK  M20AM  M273P   2  M20AM_M273P     5 2019-01-03
#> 53            1  M2ETE  M20AM  M273P   2  M20AM_M273P     5 2019-03-20
#> 54            1  M2EUJ  M20AM  M273P   2  M20AM_M273P     5 2019-04-23
#> 55            1 MSV00E  M20AM  M273P   1  M20AM_M273P     5 2019-08-12
#> 56            1 MSV018  M20AM  M273P   1  M20AM_M273P     5 2019-09-12
#> 57            1 MSV05L  M20AM  M273P   1  M20AM_M273P     5 2020-01-27
#> 58            1 MSV0M6  M20AM  M273P   1  M20AM_M273P     5 2020-02-24
#> 59            1 MSV0T4  M20AM  M273P   1  M20AM_M273P     5 2020-02-07
#> 60            1 MSV0T7  M20AM  M273P   1  M20AM_M273P     5 2019-08-11
#> 61            1 MSV0TJ  M20AM  M273P   2  M20AM_M273P     5 2019-12-28
#> 62            1 MSV0UL  M20AM  M273P   1  M20AM_M273P     5 2020-07-15
#> 63            1 MSV0X4  M20AM  M273P   1  M20AM_M273P     5 2019-09-03
#> 64            1 MSV17F  M20AM  M273P   2  M20AM_M273P     5 2020-11-08
#> 65            1 MSV1MH  M20AM  M273P   2  M20AM_M273P     5 2021-02-25
#>      LastSeen IsDead DadHSgroup MomHSgroup hsGroup
#> 1  2021-03-23  FALSE       <NA>       <NA>       0
#> 2  2020-07-22  FALSE       <NA>       <NA>       0
#> 3  2020-08-02  FALSE       <NA>       <NA>       0
#> 4  2021-04-23  FALSE       <NA>       <NA>       0
#> 5  2018-02-09   TRUE       <NA>       <NA>       0
#> 6  2018-08-22  FALSE       <NA>       <NA>       0
#> 7  2017-12-22  FALSE       <NA>     MomP_1       1
#> 8  2017-01-25   TRUE       <NA>     MomP_1       1
#> 9  2018-02-17  FALSE       <NA>     MomP_1       1
#> 10 2017-08-07  FALSE       <NA>     MomP_1       1
#> 11 2020-09-29  FALSE       <NA>     MomP_1       1
#> 12 2021-01-19  FALSE       <NA>     MomP_1       1
#> 13 2019-11-06  FALSE       <NA>     MomP_1       1
#> 14 2018-01-26  FALSE       <NA>     MomP_1       1
#> 15 2020-08-16  FALSE       <NA>       <NA>       2
#> 16 2021-03-30  FALSE       <NA>       <NA>       2
#> 17 2021-01-18  FALSE       <NA>       <NA>       2
#> 18 2021-01-23  FALSE       <NA>       <NA>       2
#> 19 2021-02-04  FALSE       <NA>       <NA>       2
#> 20 2021-03-07  FALSE       <NA>       <NA>       2
#> 21 2021-02-17  FALSE       <NA>       <NA>       2
#> 22 2021-02-26  FALSE       <NA>       <NA>       2
#> 23 2020-11-15   TRUE       <NA>       <NA>       2
#> 24 2019-08-29  FALSE       <NA>       <NA>       3
#> 25 2020-03-24  FALSE       <NA>       <NA>       3
#> 26 2021-04-29  FALSE       <NA>       <NA>       3
#> 27 2020-01-12   TRUE       <NA>       <NA>       3
#> 28 2021-04-29  FALSE       <NA>       <NA>       3
#> 29 2021-04-29  FALSE       <NA>       <NA>       3
#> 30 2020-09-29  FALSE       <NA>       <NA>       3
#> 31 2020-11-06  FALSE       <NA>       <NA>       3
#> 32 2021-02-16  FALSE       <NA>       <NA>       3
#> 33 2021-01-07  FALSE       <NA>     MomP_1       1
#> 34 2020-05-22  FALSE       <NA>     MomP_1       1
#> 35 2018-10-30   TRUE       <NA>     MomP_1       1
#> 36 2021-03-13  FALSE       <NA>     MomP_1       1
#> 37 2019-10-07   TRUE       <NA>     MomP_1       1
#> 38 2019-09-15   TRUE       <NA>     MomP_1       1
#> 39 2019-12-09  FALSE       <NA>     MomP_1       1
#> 40 2021-04-15  FALSE       <NA>     MomP_1       1
#> 41 2021-03-17  FALSE       <NA>     MomP_1       1
#> 42 2020-07-17  FALSE       <NA>     MomP_1       1
#> 43 2021-03-17  FALSE       <NA>     MomP_1       1
#> 44 2020-06-10  FALSE       <NA>     MomP_1       1
#> 45 2020-06-10  FALSE       <NA>     MomP_1       1
#> 46 2021-02-02  FALSE       <NA>     MomP_1       1
#> 47 2021-02-02  FALSE       <NA>     MomP_1       1
#> 48 2021-02-02  FALSE       <NA>     MomP_1       1
#> 49 2020-07-11  FALSE       <NA>     MomP_1       1
#> 50 2021-02-10  FALSE       <NA>     MomP_1       1
#> 51 2018-02-09  FALSE       <NA>       <NA>       4
#> 52 2020-07-08  FALSE       <NA>       <NA>       4
#> 53 2019-08-13  FALSE       <NA>       <NA>       4
#> 54 2019-08-11  FALSE       <NA>       <NA>       4
#> 55 2021-03-05  FALSE       <NA>       <NA>       4
#> 56 2019-09-12  FALSE       <NA>       <NA>       4
#> 57 2020-02-15  FALSE       <NA>       <NA>       4
#> 58 2021-04-08  FALSE       <NA>       <NA>       4
#> 59 2020-02-07   TRUE       <NA>       <NA>       4
#> 60 2020-02-09   TRUE       <NA>       <NA>       4
#> 61 2019-12-28   TRUE       <NA>       <NA>       4
#> 62 2020-07-15  FALSE       <NA>       <NA>       4
#> 63 2019-10-23   TRUE       <NA>       <NA>       4
#> 64 2020-12-04  FALSE       <NA>       <NA>       4
#> 65 2021-07-15  FALSE       <NA>       <NA>       4
#> 
#> $fams
#>         parents   father   mother FamID   FamStart     FamEnd FamDead
#> 7   M228J_M200F    M228J    M200F     1 2017-01-25 2018-08-22    TRUE
#> 15 MSV00E_M28LU   MSV00E    M28LU     2 2020-08-16 2021-03-05   FALSE
#> 24  M2772_M28TU    M2772    M28TU     3 2019-08-28 2021-04-23   FALSE
#> 33  M2AM8_M200F    M2AM8    M200F     4 2018-10-29 2021-03-23   FALSE
#> 51  M20AM_M273P    M20AM    M273P     5 2018-01-05 2020-08-02   FALSE
#> NA      Unknown *Unknown #Unknown     0 2015-07-27 2021-04-23   FALSE
#>    DadHSgroup MomHSgroup hsGroup
#> 7        <NA>     MomP_1       1
#> 15       <NA>       <NA>       2
#> 24       <NA>       <NA>       3
#> 33       <NA>     MomP_1       1
#> 51       <NA>       <NA>       4
#> NA       <NA>       <NA>       0
#> 

```
