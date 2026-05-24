# Check and prepare genetic sample metadata

Verifies the consistency of columns in the genetic sample metadata and
prepares it for use with other functions in the `wpeR` package. The
function ensures that the provided data is properly formatted and
conforms to the standards of functions that make up the `wpeR` package.

## Usage

``` r
check_sampledata(
  Sample,
  Date,
  AnimalRef,
  GeneticSex,
  lat,
  lng,
  SType,
  IsMortality = NULL,
  extraCols = NULL
)
```

## Arguments

- Sample:

  A vector of sample unique identifier codes.

- Date:

  A vector of sample collection dates in 'YYYY-MM-DD' format.

- AnimalRef:

  A vector of identifier codes of the particular individual that the
  sample belongs to.

- GeneticSex:

  A vector of genetic sex information ('F' for female, 'M' for male, NA
  for unknown).

- lat:

  A vector of latitude coordinates in the WGS84 coordinate system (EPSG:
  4326).

- lng:

  A vector of longitude coordinates in the WGS84 coordinate system
  (EPSG: 4326).

- SType:

  A vector of sample types eg.: scat, hair, tissue.

- IsMortality:

  Optional logical vector of the same length as the input vectors.
  `TRUE` marks samples that represent a mortality event (e.g. from a
  dead animal). Defaults to `NULL`.

- extraCols:

  A vector of extra column names that the user wants to include in
  sampledata data frame (see Details).

## Value

A data frame with a number of rows equal to the length of the input
vectors. Each column corresponds to one of the input parameters. If
`IsMortality` is provided, it is included as an additional column. If
the function executes without warnings or errors, the result from
`check_sampledata()` can be used as an input parameter for other
functions within this package:
[`get_colony()`](https://gr3602.github.io/wpeR/reference/get_colony.md),
[`get_ped()`](https://gr3602.github.io/wpeR/reference/get_ped.md),
[`org_fams()`](https://gr3602.github.io/wpeR/reference/org_fams.md) and
[`plot_table()`](https://gr3602.github.io/wpeR/reference/plot_table.md).

## Details

By specifying the `extraCols` parameter additional information can be
included in the sampledata dataframe. Such additional information is not
required for the functioning of the `wpeR` package functions, but can be
useful to the user when interpreting results. When including additional
columns the function inputs (Sample, Date, AnimalRef...) have to be
defined as a vector extracted from data frame column (eg. Sample =
dataframe\$column) and the `extraCols` parameter is defined, as a vector
of column names form the same data frame (eg. extraCols = c(column1,
column2, column3)).

The `IsMortality` column is a logical vector that flags samples
representing a mortality event (e.g. from a dead animal). It is used by
[`anim_timespan()`](https://gr3602.github.io/wpeR/reference/anim_timespan.md)
to determine the death status of individuals.

## Examples

``` r
sampledata <- check_sampledata(
  Sample = wolf_samples$Sample,
  Date = wolf_samples$Date,
  AnimalRef = wolf_samples$AnimalRef,
  GeneticSex = wolf_samples$GeneticSex,
  lat = wolf_samples$lat,
  lng = wolf_samples$lng,
  SType = wolf_samples$SType,
  IsMortality = wolf_samples$IsMortality
)
```
