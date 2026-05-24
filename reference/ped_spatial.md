# Get files for spatial representation of pedigree

Creates georeferenced data for spatial pedigree representation form the
output of
[`plot_table()`](https://gr3602.github.io/wpeR/reference/plot_table.md)
function.

## Usage

``` r
ped_spatial(
  plottable,
  na.rm = TRUE,
  output = "list",
  fullsibdata = NULL,
  sibthreshold = 0,
  path = "",
  filename = "",
  out.format = "geopackage",
  time.limits = c(as.Date("1900-01-01"), as.Date("2100-01-01")),
  time.limit.rep = FALSE,
  time.limit.offspring = FALSE,
  time.limit.moves = FALSE
)
```

## Arguments

- plottable:

  Data frame. Output of
  [`plot_table()`](https://gr3602.github.io/wpeR/reference/plot_table.md)
  function.

- na.rm:

  Logical (`TRUE`/`FALSE`). Remove samples with missing coordinates
  and/or dates.

- output:

  Character vector specifying the desired output type ('list' - default
  or 'gis'). Available outputs: list: all spatial data returned as list,
  gis: all spatial data returned as georeferenced files.

- fullsibdata:

  Data frame with COLONY full-sibling data.

- sibthreshold:

  Numeric. P-value threshold for sibship assignment.

- path:

  System path for storing georeferenced files.

- filename:

  Common name for all georeferenced files.

- out.format:

  Character string. Type of georeferenced files to be generated. Can be
  ether `"geopackage"` or `"shapefile"`. Default is `"geopackage"`

- time.limits:

  Vector of two `Date` values as the time window.

- time.limit.rep:

  Logical (`TRUE`/`FALSE`). Apply time limits to reference samples of
  reproductive animals.

- time.limit.offspring:

  Logical (`TRUE`/`FALSE`). Apply time limits to reference samples of
  offspring.

- time.limit.moves:

  Logical (`TRUE`/`FALSE`). Apply time limits to movement data.

## Value

Depending on the `output` parameter the function can return a list of
[`sf`](https://r-spatial.github.io/sf/reference/sf.html) objects, a
georeferenced vector data files or both.

Most of the objects are created separately for mothers, fathers and
offspring, this include:

- Reference Points (`motherRpoints`, `fatherRpoints`, and
  `offspringRpoints`).

  - Each point corresponds to an animal included in the 'plot_table()'
    function output.

  - For reproductive animals (mothers and fathers), a reference point is
    the location of their last sample within the specified time window.

  - For offspring, the reference point is the location of their first
    sample within the time window.

- Movement Points (`motherMovePoints`, `fatherMovePoints`, and
  `offspringMovePoints`).

  - These points represent all the samples of the respective animals.

- Movement Lines (`motherMoveLines`, `fatherMoveLines` and
  `offspringMoveLines`).

  - Movement lines connect all '...MovePoints' of a specific animal in
    chronological order.

- Movement Polygons (`motherMovePolygons`, `fatherMovePolygons` and
  `offspringMovePolygons`):

  - Movement polygons represent a convex hull that encloses all the
    samples of an individual.

  - An individual must have more than two samples for this
    representation.

Besides that the function also produces lines that connect mothers and
their offspring (`maternityLines`), fathers and their offspring
(`paternityLines`), and if `fullsibdata` parameter is specified, full
siblings (`FullsibLines`).

## Details

The parameters `path`, `filename` and `out.format`, are used only when
`output` parameter is set to "gis", since they control which
georeferenced files should be created, where they will be saved and
which common file name will they have.

## Examples

``` r
# Prepare the data for usage with ped_spatial() function.
# Get animal timespan data using the anim_timespan() function.
animal_ts <- anim_timespan(wolf_samples$AnimalRef,
  wolf_samples$Date,
  wolf_samples$IsMortality
)
# Add animal timespan to the sampledata
sampledata <- merge(wolf_samples, animal_ts, by.x = "AnimalRef", by.y = "ID", all.x = TRUE)
# Define the path to the pedigree data file.
path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
# Retrieve the pedigree data from the get_colony function.
ped_colony <- get_colony(path, sampledata, rm_obsolete_parents = TRUE, out = "FamAgg")
# Organize families and expand pedigree data using the org_fams function.
org_tables <- org_fams(ped_colony, sampledata, output = "both")
# Prepare data for plotting.
pt <- plot_table(plot_fams = 1,
  all_fams = org_tables$fams,
  ped = org_tables$ped,
  sampledata = sampledata,
)

# Run the function
# Get files for spatial pedigree representation in list format.
ped_spatial(plottable = pt)
#> $motherRpoints
#> Simple feature collection with 1 feature and 15 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 14.14495 ymin: 45.71833 xmax: 14.14495 ymax: 45.71833
#> Geodetic CRS:  WGS 84
#>    Sample AnimalRef GeneticSex       Date SType  FirstSeen   LastSeen IsDead
#> 20  M2AMA     M200F          F 2018-08-22  Scat 2015-07-27 2018-08-22  FALSE
#>    IsMortality  rep later_rep isPolygamous  dead first_sample last_sample
#> 20       FALSE TRUE     FALSE         TRUE FALSE        FALSE        TRUE
#>                     geometry
#> 20 POINT (14.14495 45.71833)
#> 
#> $fatherRpoints
#> Simple feature collection with 1 feature and 15 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 14.15404 ymin: 45.67033 xmax: 14.15404 ymax: 45.67033
#> Geodetic CRS:  WGS 84
#>    Sample AnimalRef GeneticSex       Date  SType  FirstSeen   LastSeen IsDead
#> 58  M2C36     M228J          M 2018-02-09 Tissue 2016-09-30 2018-02-09   TRUE
#>    IsMortality  rep later_rep isPolygamous dead first_sample last_sample
#> 58        TRUE TRUE     FALSE        FALSE TRUE        FALSE        TRUE
#>                     geometry
#> 58 POINT (14.15404 45.67033)
#> 
#> $offspringRpoints
#> Simple feature collection with 8 features and 18 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 14.04242 ymin: 45.67397 xmax: 14.1563 ymax: 45.71845
#> Geodetic CRS:  WGS 84
#>     Sample AnimalRef GeneticSex       Date  SType  FirstSeen   LastSeen IsDead
#> 1    M10XC     M10XC          M 2017-11-16   Scat 2017-11-16 2017-12-22  FALSE
#> 8    M1XA7     M1YP0          F 2017-01-25   Scat 2017-01-25 2017-01-25   TRUE
#> 43   M2C2F     M220J          M 2017-11-10 Saliva 2017-11-10 2018-02-17  FALSE
#> 59   M22AM     M22AM          M 2017-01-26   Scat 2017-01-26 2017-08-07  FALSE
#> 71   M2772     M2772          M 2017-11-12   Scat 2017-11-12 2020-09-29  FALSE
#> 77   M28LU     M28LU          F 2017-09-18 Saliva 2017-09-18 2021-01-19  FALSE
#> 125  M275E     M2AXE          F 2017-12-16   Scat 2017-12-16 2019-11-06  FALSE
#> 131  M277F     M2C1T          F 2017-10-31   Scat 2017-10-31 2018-01-26  FALSE
#>     IsMortality plottingID FamID hsGroup   rep later_rep isPolygamous  dead
#> 1         FALSE          3     1       1 FALSE     FALSE        FALSE FALSE
#> 8         FALSE          4     1       1 FALSE     FALSE        FALSE FALSE
#> 43        FALSE          5     1       1 FALSE     FALSE        FALSE FALSE
#> 59        FALSE          6     1       1 FALSE     FALSE        FALSE FALSE
#> 71        FALSE          7     1       1 FALSE      TRUE        FALSE FALSE
#> 77        FALSE          8     1       1 FALSE      TRUE        FALSE FALSE
#> 125       FALSE          9     1       1 FALSE     FALSE        FALSE FALSE
#> 131       FALSE         10     1       1 FALSE     FALSE        FALSE FALSE
#>     first_sample last_sample                  geometry
#> 1           TRUE       FALSE POINT (14.12922 45.70766)
#> 8           TRUE       FALSE POINT (14.12954 45.71845)
#> 43          TRUE       FALSE  POINT (14.04242 45.7131)
#> 59          TRUE       FALSE POINT (14.12501 45.70004)
#> 71          TRUE       FALSE  POINT (14.1563 45.69639)
#> 77          TRUE       FALSE POINT (14.11112 45.67397)
#> 125         TRUE       FALSE POINT (14.14771 45.71293)
#> 131         TRUE       FALSE  POINT (14.1357 45.71323)
#> 
#> $motherMovePoints
#> Simple feature collection with 11 features and 15 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 14.03169 ymin: 45.68138 xmax: 14.14653 ymax: 45.7525
#> Geodetic CRS:  WGS 84
#> First 10 features:
#>     Sample AnimalRef GeneticSex       Date  SType  FirstSeen   LastSeen IsDead
#> 10 EX.1JH0     M200F          F 2015-07-27 Saliva 2015-07-27 2018-08-22  FALSE
#> 11 EX.1JJ5     M200F          F 2015-08-14 Saliva 2015-07-27 2018-08-22  FALSE
#> 12   M20A4     M200F          F 2016-10-11 Saliva 2015-07-27 2018-08-22  FALSE
#> 13   M200F     M200F          F 2017-01-11   Scat 2015-07-27 2018-08-22  FALSE
#> 14   M1TU6     M200F          F 2017-01-25   Scat 2015-07-27 2018-08-22  FALSE
#> 15   M1XPY     M200F          F 2017-01-25   Scat 2015-07-27 2018-08-22  FALSE
#> 16   M1YH8     M200F          F 2017-03-20   Scat 2015-07-27 2018-08-22  FALSE
#> 17   M1T7L     M200F          F 2017-04-21   Scat 2015-07-27 2018-08-22  FALSE
#> 18   M2ATF     M200F          F 2018-07-04   Scat 2015-07-27 2018-08-22  FALSE
#> 19   M2ATC     M200F          F 2018-07-19   Scat 2015-07-27 2018-08-22  FALSE
#>    IsMortality  rep later_rep isPolygamous  dead first_sample last_sample
#> 10       FALSE TRUE     FALSE         TRUE FALSE         TRUE       FALSE
#> 11       FALSE TRUE     FALSE         TRUE FALSE        FALSE       FALSE
#> 12       FALSE TRUE     FALSE         TRUE FALSE        FALSE       FALSE
#> 13       FALSE TRUE     FALSE         TRUE FALSE        FALSE       FALSE
#> 14       FALSE TRUE     FALSE         TRUE FALSE        FALSE       FALSE
#> 15       FALSE TRUE     FALSE         TRUE FALSE        FALSE       FALSE
#> 16       FALSE TRUE     FALSE         TRUE FALSE        FALSE       FALSE
#> 17       FALSE TRUE     FALSE         TRUE FALSE        FALSE       FALSE
#> 18       FALSE TRUE     FALSE         TRUE FALSE        FALSE       FALSE
#> 19       FALSE TRUE     FALSE         TRUE FALSE        FALSE       FALSE
#>                     geometry
#> 10  POINT (14.14653 45.7525)
#> 11 POINT (14.03169 45.69501)
#> 12 POINT (14.04434 45.71347)
#> 13 POINT (14.08802 45.71602)
#> 14 POINT (14.12843 45.69746)
#> 15 POINT (14.08937 45.71062)
#> 16 POINT (14.07037 45.70856)
#> 17 POINT (14.11843 45.68138)
#> 18 POINT (14.07282 45.70502)
#> 19 POINT (14.09111 45.70253)
#> 
#> $fatherMovePoints
#> Simple feature collection with 5 features and 15 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 14.01201 ymin: 45.67033 xmax: 14.15404 ymax: 45.7114
#> Geodetic CRS:  WGS 84
#>    Sample AnimalRef GeneticSex       Date  SType  FirstSeen   LastSeen IsDead
#> 54  M20AP     M228J          M 2016-09-30 Saliva 2016-09-30 2018-02-09   TRUE
#> 55  M228J     M228J          M 2017-01-26   Scat 2016-09-30 2018-02-09   TRUE
#> 56  M28ML     M228J          M 2017-08-18 Saliva 2016-09-30 2018-02-09   TRUE
#> 57  M28MM     M228J          M 2017-08-18 Saliva 2016-09-30 2018-02-09   TRUE
#> 58  M2C36     M228J          M 2018-02-09 Tissue 2016-09-30 2018-02-09   TRUE
#>    IsMortality  rep later_rep isPolygamous  dead first_sample last_sample
#> 54       FALSE TRUE     FALSE        FALSE FALSE         TRUE       FALSE
#> 55       FALSE TRUE     FALSE        FALSE FALSE        FALSE       FALSE
#> 56       FALSE TRUE     FALSE        FALSE FALSE        FALSE       FALSE
#> 57       FALSE TRUE     FALSE        FALSE FALSE        FALSE       FALSE
#> 58        TRUE TRUE     FALSE        FALSE  TRUE        FALSE        TRUE
#>                     geometry
#> 54  POINT (14.01201 45.7114)
#> 55 POINT (14.12798 45.70406)
#> 56  POINT (14.1115 45.67397)
#> 57  POINT (14.1115 45.67397)
#> 58 POINT (14.15404 45.67033)
#> 
#> $offspringMovePoints
#> Simple feature collection with 35 features and 18 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 13.80997 ymin: 45.52787 xmax: 14.39713 ymax: 46.28649
#> Geodetic CRS:  WGS 84
#> First 10 features:
#>    Sample AnimalRef GeneticSex       Date  SType  FirstSeen   LastSeen IsDead
#> 1   M10XC     M10XC          M 2017-11-16   Scat 2017-11-16 2017-12-22  FALSE
#> 2   M0PXH     M10XC          M 2017-11-22   Scat 2017-11-16 2017-12-22  FALSE
#> 3   M0PFL     M10XC          M 2017-12-22   Scat 2017-11-16 2017-12-22  FALSE
#> 8   M1XA7     M1YP0          F 2017-01-25   Scat 2017-01-25 2017-01-25   TRUE
#> 9   M1YP0     M1YP0          F 2017-01-25 Tissue 2017-01-25 2017-01-25   TRUE
#> 43  M2C2F     M220J          M 2017-11-10 Saliva 2017-11-10 2018-02-17  FALSE
#> 44  M220J     M220J          M 2018-02-17   Scat 2017-11-10 2018-02-17  FALSE
#> 59  M22AM     M22AM          M 2017-01-26   Scat 2017-01-26 2017-08-07  FALSE
#> 60  M28MF     M22AM          M 2017-08-07 Saliva 2017-01-26 2017-08-07  FALSE
#> 71  M2772     M2772          M 2017-11-12   Scat 2017-11-12 2020-09-29  FALSE
#>    IsMortality plottingID FamID hsGroup   rep later_rep isPolygamous  dead
#> 1        FALSE          3     1       1 FALSE     FALSE        FALSE FALSE
#> 2        FALSE          3     1       1 FALSE     FALSE        FALSE FALSE
#> 3        FALSE          3     1       1 FALSE     FALSE        FALSE FALSE
#> 8        FALSE          4     1       1 FALSE     FALSE        FALSE FALSE
#> 9         TRUE          4     1       1 FALSE     FALSE        FALSE  TRUE
#> 43       FALSE          5     1       1 FALSE     FALSE        FALSE FALSE
#> 44       FALSE          5     1       1 FALSE     FALSE        FALSE FALSE
#> 59       FALSE          6     1       1 FALSE     FALSE        FALSE FALSE
#> 60       FALSE          6     1       1 FALSE     FALSE        FALSE FALSE
#> 71       FALSE          7     1       1 FALSE      TRUE        FALSE FALSE
#>    first_sample last_sample                  geometry
#> 1          TRUE       FALSE POINT (14.12922 45.70766)
#> 2         FALSE       FALSE POINT (14.10497 45.71356)
#> 3         FALSE        TRUE POINT (14.07907 45.69898)
#> 8          TRUE       FALSE POINT (14.12954 45.71845)
#> 9         FALSE        TRUE POINT (14.11387 45.71246)
#> 43         TRUE       FALSE  POINT (14.04242 45.7131)
#> 44        FALSE        TRUE POINT (14.15949 45.71823)
#> 59         TRUE       FALSE POINT (14.12501 45.70004)
#> 60        FALSE        TRUE  POINT (13.80997 45.8081)
#> 71         TRUE       FALSE  POINT (14.1563 45.69639)
#> 
#> $maternityLines
#> Simple feature collection with 8 features and 7 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: 14.04242 ymin: 45.67397 xmax: 14.1563 ymax: 45.71845
#> Geodetic CRS:  WGS 84
#>   ID pair fam plyClust  relation child parent                       geometry
#> 1  1    1   1        1 maternity M10XC  M200F LINESTRING (14.14495 45.718...
#> 2  2    2   1        1 maternity M1YP0  M200F LINESTRING (14.14495 45.718...
#> 3  3    3   1        1 maternity M220J  M200F LINESTRING (14.14495 45.718...
#> 4  4    4   1        1 maternity M22AM  M200F LINESTRING (14.14495 45.718...
#> 5  5    5   1        1 maternity M2772  M200F LINESTRING (14.14495 45.718...
#> 6  6    6   1        1 maternity M28LU  M200F LINESTRING (14.14495 45.718...
#> 7  7    7   1        1 maternity M2AXE  M200F LINESTRING (14.14495 45.718...
#> 8  8    8   1        1 maternity M2C1T  M200F LINESTRING (14.14495 45.718...
#> 
#> $paternityLines
#> Simple feature collection with 8 features and 7 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: 14.04242 ymin: 45.67033 xmax: 14.1563 ymax: 45.71845
#> Geodetic CRS:  WGS 84
#>   ID pair fam plyClust  relation child parent                       geometry
#> 1  1    1   1        1 paternity M10XC  M228J LINESTRING (14.15404 45.670...
#> 2  2    2   1        1 paternity M1YP0  M228J LINESTRING (14.15404 45.670...
#> 3  3    3   1        1 paternity M220J  M228J LINESTRING (14.15404 45.670...
#> 4  4    4   1        1 paternity M22AM  M228J LINESTRING (14.15404 45.670...
#> 5  5    5   1        1 paternity M2772  M228J LINESTRING (14.15404 45.670...
#> 6  6    6   1        1 paternity M28LU  M228J LINESTRING (14.15404 45.670...
#> 7  7    7   1        1 paternity M2AXE  M228J LINESTRING (14.15404 45.670...
#> 8  8    8   1        1 paternity M2C1T  M228J LINESTRING (14.15404 45.670...
#> 
#> $motherMoveLines
#> Simple feature collection with 1 feature and 2 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: 14.03169 ymin: 45.68138 xmax: 14.14653 ymax: 45.7525
#> Geodetic CRS:  WGS 84
#> # A tibble: 1 × 3
#>   AnimalRef no_mvPoints                                                 geometry
#>   <chr>           <int>                                         <LINESTRING [°]>
#> 1 M200F              11 (14.14653 45.7525, 14.03169 45.69501, 14.04434 45.71347…
#> 
#> $fatherMoveLines
#> Simple feature collection with 1 feature and 2 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: 14.01201 ymin: 45.67033 xmax: 14.15404 ymax: 45.7114
#> Geodetic CRS:  WGS 84
#> # A tibble: 1 × 3
#>   AnimalRef no_mvPoints                                                 geometry
#>   <chr>           <int>                                         <LINESTRING [°]>
#> 1 M228J               5 (14.01201 45.7114, 14.12798 45.70406, 14.1115 45.67397,…
#> 
#> $offspringMoveLines
#> Simple feature collection with 8 features and 2 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: 13.80997 ymin: 45.52787 xmax: 14.39713 ymax: 46.28649
#> Geodetic CRS:  WGS 84
#> # A tibble: 8 × 3
#>   AnimalRef no_mvPoints                                                 geometry
#>   <chr>           <int>                                         <LINESTRING [°]>
#> 1 M10XC               3 (14.12922 45.70766, 14.10497 45.71356, 14.07907 45.6989…
#> 2 M1YP0               2                   (14.12954 45.71845, 14.11387 45.71246)
#> 3 M220J               2                    (14.04242 45.7131, 14.15949 45.71823)
#> 4 M22AM               2                    (14.12501 45.70004, 13.80997 45.8081)
#> 5 M2772               6 (14.1563 45.69639, 14.1013 46.2444, 14.0529 46.23509, 1…
#> 6 M28LU               9 (14.11112 45.67397, 14.07135 45.73484, 14.39713 45.6606…
#> 7 M2AXE               6 (14.14771 45.71293, 14.07355 45.7183, 14.06085 45.70965…
#> 8 M2C1T               5 (14.1357 45.71323, 14.22777 45.57041, 14.22777 45.57041…
#> 
#> $motherMovePolygons
#> Simple feature collection with 1 feature and 2 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 14.03169 ymin: 45.68138 xmax: 14.14653 ymax: 45.7525
#> Geodetic CRS:  WGS 84
#> # A tibble: 1 × 3
#>   AnimalRef no_mvPoints                                                 geometry
#>   <chr>           <int>                                            <POLYGON [°]>
#> 1 M200F              11 ((14.03169 45.69501, 14.11843 45.68138, 14.14495 45.718…
#> 
#> $fatherMovePolygons
#> Simple feature collection with 1 feature and 2 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 14.01201 ymin: 45.67033 xmax: 14.15404 ymax: 45.7114
#> Geodetic CRS:  WGS 84
#> # A tibble: 1 × 3
#>   AnimalRef no_mvPoints                                                 geometry
#>   <chr>           <int>                                            <POLYGON [°]>
#> 1 M228J               5 ((14.01201 45.7114, 14.1115 45.67397, 14.15404 45.67033…
#> 
#> $offspringMovePolygons
#> Simple feature collection with 8 features and 2 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 13.80997 ymin: 45.52787 xmax: 14.39713 ymax: 46.28649
#> Geodetic CRS:  WGS 84
#> # A tibble: 8 × 3
#>   AnimalRef no_mvPoints                                                 geometry
#>   <chr>           <int>                                            <POLYGON [°]>
#> 1 M10XC               3 ((14.07907 45.69898, 14.12922 45.70766, 14.10497 45.713…
#> 2 M1YP0               2 ((14.11387 45.71246, 14.12954 45.71845, 14.12171 45.715…
#> 3 M220J               2 ((14.10095 45.71568, 14.15949 45.71823, 14.04242 45.713…
#> 4 M22AM               2 ((13.80997 45.8081, 14.12501 45.70004, 13.96764 45.7541…
#> 5 M2772               6 ((13.97469 46.24017, 14.1563 45.69639, 14.1013 46.2444,…
#> 6 M28LU               9 ((14.07135 45.73484, 14.11112 45.67397, 14.27129 45.598…
#> 7 M2AXE               6 ((14.06085 45.70965, 14.02367 45.52787, 14.14771 45.712…
#> 8 M2C1T               5 ((14.18179 45.64183, 14.22777 45.57041, 14.1357 45.7132…
#> 

```
