# Get dates of individuals first and last sample

Takes data frame of all samples and returns the dates of individuals
first and last sample. Besides that the function determines if an
individual is dead based on the mortality_sample column.

## Usage

``` r
anim_timespan(individual_id, sample_date, mortality_sample)
```

## Arguments

- individual_id:

  Column in the dataframe of all samples containing individual animal
  identifier code. Defined as `dataframe$column`.

- sample_date:

  Column in the dataframe of all samples containing the date of sample
  collection. Must be in `Date` format. Defined as `dataframe$column`.

- mortality_sample:

  Logical vector or column in the dataframe of all samples that
  identifies samples that represent a mortality event (e.g. from a dead
  animal). `TRUE` values mark mortality samples. Defined as
  `dataframe$column` or a standalone logical vector the same length as
  `individual_id`.

## Value

A data frame with four columns and one row for each `individual_id`.
Returned data frame columns correspond to individual identification key
(`ID`), date of first (`FirstSeen`) and last (`LastSeen`) sample of
individual and logical (`TRUE/FALSE`) value that identifies if the
individual is dead (`IsDead`).

## Examples

``` r
anim_timespan(
  individual_id = wolf_samples$AnimalRef,
  sample_date = wolf_samples$Date,
  mortality_sample = wolf_samples$IsMortality
)
#>        ID  FirstSeen   LastSeen IsDead
#> 1   M10XC 2017-11-16 2017-12-22  FALSE
#> 2   M1J47 2019-08-20 2021-01-07  FALSE
#> 3   M1YP0 2017-01-25 2017-01-25   TRUE
#> 4   M200F 2015-07-27 2018-08-22  FALSE
#> 5   M20AM 2016-08-29 2020-08-02  FALSE
#> 6   M220J 2017-11-10 2018-02-17  FALSE
#> 7   M221C 2018-10-29 2020-05-22  FALSE
#> 8   M228J 2016-09-30 2018-02-09   TRUE
#> 9   M22AM 2017-01-26 2017-08-07  FALSE
#> 10  M273P 2018-01-02 2020-07-22  FALSE
#> 11  M2757 2018-01-05 2018-02-09  FALSE
#> 12  M2772 2017-11-12 2020-09-29  FALSE
#> 13  M28LU 2017-09-18 2021-01-19  FALSE
#> 14  M28TU 2017-12-18 2021-04-23  FALSE
#> 15  M2ALK 2019-01-03 2020-07-08  FALSE
#> 16  M2AM8 2017-04-07 2021-03-23  FALSE
#> 17  M2AXE 2017-12-16 2019-11-06  FALSE
#> 18  M2C1T 2017-10-31 2018-01-26  FALSE
#> 19  M2C8Y 2018-10-30 2018-10-30   TRUE
#> 20  M2ETE 2019-03-20 2019-08-13  FALSE
#> 21  M2EUJ 2019-04-23 2019-08-11  FALSE
#> 22  M2F1L 2019-04-07 2021-03-13  FALSE
#> 23 MSV00E 2019-08-12 2021-03-05  FALSE
#> 24 MSV018 2019-09-12 2019-09-12  FALSE
#> 25 MSV01X 2019-08-28 2019-08-29  FALSE
#> 26 MSV02F 2019-10-07 2019-10-07   TRUE
#> 27 MSV02L 2019-09-15 2019-09-15   TRUE
#> 28 MSV055 2019-12-11 2020-03-24  FALSE
#> 29 MSV05L 2020-01-27 2020-02-15  FALSE
#> 30 MSV0AL 2019-12-11 2021-04-29  FALSE
#> 31 MSV0CK 2019-08-05 2019-12-09  FALSE
#> 32 MSV0FK 2020-01-18 2021-04-15  FALSE
#> 33 MSV0H5 2020-04-10 2021-03-17  FALSE
#> 34 MSV0M6 2020-02-24 2021-04-08  FALSE
#> 35 MSV0P7 2019-11-23 2020-07-17  FALSE
#> 36 MSV0T4 2020-02-07 2020-02-07   TRUE
#> 37 MSV0T7 2019-08-11 2020-02-09   TRUE
#> 38 MSV0TA 2020-01-12 2020-01-12   TRUE
#> 39 MSV0TJ 2019-12-28 2019-12-28   TRUE
#> 40 MSV0UL 2020-07-15 2020-07-15  FALSE
#> 41 MSV0UP 2020-06-10 2021-03-17  FALSE
#> 42 MSV0UT 2020-06-10 2020-06-10  FALSE
#> 43 MSV0UU 2020-06-10 2020-06-10  FALSE
#> 44 MSV0X4 2019-09-03 2019-10-23   TRUE
#> 45 MSV0XT 2020-10-08 2021-04-29  FALSE
#> 46 MSV10T 2020-08-16 2020-08-16  FALSE
#> 47 MSV16T 2021-02-02 2021-02-02  FALSE
#> 48 MSV16U 2021-02-02 2021-02-02  FALSE
#> 49 MSV170 2021-02-02 2021-02-02  FALSE
#> 50 MSV177 2020-04-22 2020-07-11  FALSE
#> 51 MSV17F 2020-11-08 2020-12-04  FALSE
#> 52 MSV17U 2020-09-29 2021-04-29  FALSE
#> 53 MSV180 2020-09-29 2020-09-29  FALSE
#> 54 MSV18C 2020-11-06 2020-11-06  FALSE
#> 55 MSV1C0 2021-02-16 2021-03-30  FALSE
#> 56 MSV1EX 2021-01-18 2021-01-18  FALSE
#> 57 MSV1F5 2021-01-23 2021-01-23  FALSE
#> 58 MSV1F8 2021-02-04 2021-02-04  FALSE
#> 59 MSV1FE 2021-02-10 2021-02-10  FALSE
#> 60 MSV1FJ 2021-02-16 2021-03-07  FALSE
#> 61 MSV1FL 2021-02-17 2021-02-17  FALSE
#> 62 MSV1FT 2021-01-18 2021-02-26  FALSE
#> 63 MSV1KT 2020-10-06 2021-02-16  FALSE
#> 64 MSV1MH 2021-02-25 2021-07-15  FALSE
#> 65 MSV1TM 2020-11-15 2020-11-15   TRUE
```
