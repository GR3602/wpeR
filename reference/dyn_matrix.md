# Get matrix of apparent survival

Creates a matrix that shows number of captured animals between multiple
seasons.

## Usage

``` r
dyn_matrix(animal_id, capture_date, start_dates, end_dates)
```

## Arguments

- animal_id:

  A column in the dataframe of all samples that stores individual animal
  identifier code.

- capture_date:

  A column in the dataframe of all samples that stores the date of
  sample collection. Must be in `Date` format.

- start_dates:

  Vector of dates in `Date` format that define the start of each season.

- end_dates:

  Vector of dates in `Date` format that define the end of each season.

## Value

A matrix with 1 + no. seasons rows and columns.

- diagonal: number of new captures in each session,

- above diagonal: number of recaptures from season x to season y,

- below diagonal: number of animals from season y that skipped season x.

Season x is defined in first row, season y in first column. Column
`Tot. Capts` gives all detected individuals in season y. Row
`Tot. Skipped` gives all individuals skipped in season x but detected
later.

## Examples

``` r
# Define start and end dates for sampling seasons.
seasons <- data.frame(
  start = c(
    as.Date("2017-01-01"),
    as.Date("2018-01-01"),
    as.Date("2019-01-01")
  ),
  end = c(
    as.Date("2017-12-31"),
    as.Date("2018-12-31"),
    as.Date("2019-12-31")
  )
)

# Create a dynamics matrix for animal captures.
dyn_matrix(
  animal_id = wolf_samples$AnimalRef,
  capture_date = wolf_samples$Date,
  start_dates = seasons$start,
  end_dates = seasons$end
)
#>                         2017-01-01 - 2017-12-31 2018-01-01 - 2018-12-31
#> 2017-01-01 - 2017-12-31                      13                       8
#> 2018-01-01 - 2018-12-31                       2                       4
#> 2019-01-01 - 2019-12-31                       0                       0
#> Tot. Skipped                                  0                       2
#>                         2019-01-01 - 2019-12-31 Tot. Capts
#> 2017-01-01 - 2017-12-31                       6         13
#> 2018-01-01 - 2018-12-31                       6         12
#> 2019-01-01 - 2019-12-31                      17         25
#> Tot. Skipped                                  0         NA

```
