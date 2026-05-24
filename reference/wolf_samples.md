# Wolf monitoring genetic samples metadata

Metadata of selected genetic samples of wolves collected between 2015
and 2021, in the scope of Slovenian National Wolf Monitoring

## Usage

``` r
wolf_samples
```

## Format

A data frame with 407 rows and 8 columns:

- Sample:

  Sample unique identifier code

- Date:

  Date of sample collection (format: `YYYY-MM-DD`)

- AnimalRef:

  Identification string for particular animal

- GeneticSex:

  Sex of animal to which the sample belong (format: `M` = male, `F` =
  female)

- lat:

  latitude (N-S) of the sample (CRS: WGS84; EPSG: 4326)

- lng:

  longitude (W-E) of the sample (CRS: WGS84; EPSG: 4326)

- SType:

  Type of the sample.
  (`Direct Saliva, Scat, Urine, Saliva, Tissue, Decomposing Tissue, Blood`)

- IsMortality:

  A logical column (TRUE/FALSE) that marks samples collected from dead
  animals. TRUE = mortality event (e.g., a tissue sample from a
  carcass), FALSE = non-mortality sample (e.g., scat, hair, saliva from
  a live animal).

## Source

Slovenian National Wolf Monitoring
