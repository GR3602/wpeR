# Number of detected animals between two sampling seasons

Gives an numeric overview of individuals captured within the second
sampling season compared tho the first one.

## Usage

``` r
nbtw_seasons(
  animal_id,
  capture_date,
  season1_start,
  season1_end,
  season2_start,
  season2_end
)
```

## Arguments

- animal_id:

  A column in the dataframe of all samples that stores individual animal
  identifier code.

- capture_date:

  A column in the dataframe of all samples that stores the date of
  sample collection. Must be in `Date` format.

- season1_start:

  String in `Date` format. Start of fist capture season. Start and end
  date are included in the capture season.

- season1_end:

  String in `Date` format. End of fist capture season. Start and end
  date are included in the capture season.

- season2_start:

  String in `Date` format. Start of second capture season. Start and end
  date are included in the capture season.

- season2_end:

  String in `Date` format. End of second capture season. Start and end
  date are included in the capture season.

## Value

A data frame with one row and six columns corresponding to season 1 and
2 start and end dates, number of detected animals in season 2
(`total_cap`), number of new detentions in season 2 (`new_captures`),
umber of animals from season 1 detected within season 2 (`recaptured`)
and number of individuals skipped in season 2 but detected after the end
of that season (`skipped`).

## Examples

``` r
# Calculate the number of animals detected between two sampling seasons.
nbtw_seasons(
 animal_id = wolf_samples$AnimalRef,
 capture_date = wolf_samples$Date,
 season1_start = as.Date("2017-01-01"),
 season1_end = as.Date("2017-12-31"),
 season2_start = as.Date("2018-01-01"),
 season2_end = as.Date("2018-12-31")
)
#>                   season1                 season2 total_cap new_captures
#> 1 2017-01-01 - 2017-12-31 2018-01-01 - 2018-12-31        12            4
#>   recaptures skipped
#> 1          8       2


```
