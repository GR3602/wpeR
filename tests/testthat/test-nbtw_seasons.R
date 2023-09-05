

dates <- as.Date(c("2017-01-01", "2017-12-31", "2018-01-01", "2018-12-31"))

test_that("wrong date sequence warning", {
  expect_warning(nbtw_seasons(wolf_samples$AnimalRef, wolf_samples$Date,
               season1_start = dates[1],
               season1_end = dates[3],
               season2_start = dates[2],
               season2_end = dates[4]))
})

test_that("wrong format error", {
  dates <- c("2017-01-01", "2017-12-31", "2018-01-01", "2018-12-31")
  expect_error(nbtw_seasons(wolf_samples$AnimalRef, wolf_samples$Date,
                            season1_start = dates[1],
                            season1_end = dates[2],
                            season2_start = dates[3],
                            season2_end = dates[4]))
})


test_that("simple result check", {
  result <- nbtw_seasons(wolf_samples$AnimalRef, wolf_samples$Date,
                         season1_start = dates[1],
                         season1_end = dates[2],
                         season2_start = dates[3],
                         season2_end = dates[4])
  expect_equal(result$total_cap, 12, label = "Total Captures")
  expect_equal(result$new_captures, 4, label = "Number of New Captures")
  expect_equal(result$recaptures, 8, label = "Number of Recaptures")
  expect_equal(result$skipped, 2, label = "Number of Skipped Animals")
})


