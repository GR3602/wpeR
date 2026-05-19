test_that("check_sampledata returns a well-structured data frame", {
  result <- check_sampledata(
    Sample = wolf_samples$Sample,
    Date = wolf_samples$Date,
    AnimalRef = wolf_samples$AnimalRef,
    GeneticSex = wolf_samples$GeneticSex,
    lat = wolf_samples$lat,
    lng = wolf_samples$lng,
    SType = wolf_samples$SType
  )

  expect_true(inherits(result, "data.frame"))

  expected_columns <- c("Sample", "Date", "AnimalRef", "GeneticSex", "lat", "lng", "SType")
  expect_equal(colnames(result), expected_columns)
})

test_that("duplicated sample", {
  expect_error(
    check_sampledata(
      Sample = c("M10XC","M0PXH","M0PXH"),
      Date = c("2017-10-13","2017-11-22","2019-08-20"),
      AnimalRef = c("M10XC","M10XC","M1J47"),
      GeneticSex = c("M","M","F"),
      lat = c("45.707663","45.713559","45.698983"),
      lng = c("14.129219","14.104967","14.079067"),
      SType = c("Scat","Scat","Saliva")
    )
  )
})



test_that("date error", {
expect_error(
  check_sampledata(
    Sample = c("M10XC","M0PXH","M1J47"),
    Date = c("16.3.2023","2017-11-22","2019-08-20"),
    AnimalRef = c("M10XC","M10XC","M1J47"),
    GeneticSex = c("M","M","F"),
    lat = c("45.707663","45.713559","45.698983"),
    lng = c("14.129219","14.104967","14.079067"),
    SType = c("Scat","Scat","Saliva")
  )
)
})


test_that("sex error", {
  expect_error(
    check_sampledata(
      Sample = c("M10XC","M0PXH","M1J47"),
      Date = c("2017-10-13","2017-11-22","2019-08-20"),
      AnimalRef = c("M10XC","M10XC","M1J47"),
      GeneticSex = c("M","UNKNOWN","F"),
      lat = c("45.707663","45.713559","45.698983"),
      lng = c("14.129219","14.104967","14.079067"),
      SType = c("Scat","Scat","Saliva")
    )
  )
})




test_that("lat error", {
  expect_error(
    check_sampledata(
      Sample = c("M10XC","M0PXH","M1J47"),
      Date = c("2017-10-13","2017-11-22","2019-08-20"),
      AnimalRef = c("M10XC","M10XC","M1J47"),
      GeneticSex = c("M","M","F"),
      IsAnimalReference = c("1","0","1"),
      lat = c("348233","45.713559","45.698983"),
      lng = c("14.129219","14.104967","14.079067"),
      SType = c("Scat","Scat","Saliva")
    )
  )
})

test_that("lng error", {
  expect_error(
    check_sampledata(
      Sample = c("M10XC","M0PXH","M1J47"),
      Date = c("2017-10-13","2017-11-22","2019-08-20"),
      AnimalRef = c("M10XC","M10XC","M1J47"),
      GeneticSex = c("M","M","F"),
      IsAnimalReference = c("1","0","1"),
      lat = c("45.707663","45.713559","45.698983"),
      lng = c("1423122","14.104967","14.079067"),
      SType = c("Scat","Scat","Saliva")
    )
  )
})


test_that("IsMortality column is included when provided", {
  mortality <- wolf_samples$SType == "Tissue"
  result <- check_sampledata(
    Sample = wolf_samples$Sample,
    Date = wolf_samples$Date,
    AnimalRef = wolf_samples$AnimalRef,
    GeneticSex = wolf_samples$GeneticSex,
    lat = wolf_samples$lat,
    lng = wolf_samples$lng,
    SType = wolf_samples$SType,
    IsMortality = mortality
  )
  expect_true("IsMortality" %in% colnames(result))
  expect_type(result$IsMortality, "logical")
  expect_equal(sum(result$IsMortality), sum(mortality))
})

test_that("IsMortality length mismatch errors", {
  expect_error(
    check_sampledata(
      Sample = c("M10XC", "M0PXH", "M1J47"),
      Date = c("2017-10-13", "2017-11-22", "2019-08-20"),
      AnimalRef = c("M10XC", "M10XC", "M1J47"),
      GeneticSex = c("M", "M", "F"),
      lat = c("45.707663", "45.713559", "45.698983"),
      lng = c("14.129219", "14.104967", "14.079067"),
      SType = c("Scat", "Scat", "Saliva"),
      IsMortality = c(TRUE, FALSE)
    )
  )
})

test_that("IsMortality non-logical type errors", {
  expect_error(
    check_sampledata(
      Sample = c("M10XC", "M0PXH", "M1J47"),
      Date = c("2017-10-13", "2017-11-22", "2019-08-20"),
      AnimalRef = c("M10XC", "M10XC", "M1J47"),
      GeneticSex = c("M", "M", "F"),
      lat = c("45.707663", "45.713559", "45.698983"),
      lng = c("14.129219", "14.104967", "14.079067"),
      SType = c("Scat", "Scat", "Saliva"),
      IsMortality = c("TRUE", "FALSE", "TRUE")
    )
  )
})

#removed this limitation form the pacakge =)
# test_that("AnimalRef not in Sample", {
#   expect_error(
#     check_sampledata(
#       Sample = c("M10XC","M0PXH","M1J47"),
#       Date = c("2017-10-13","2017-11-22","2019-08-20"),
#       AnimalRef = c("M10XC","M10XC","M47"),
#       GeneticSex = c("M","M","F"),
#       lat = c("45.707663","45.713559","45.698983"),
#       lng = c("14.23122","14.104967","14.079067"),
#       SType = c("Scat","Scat","Saliva")
#     )
#   )
# })
