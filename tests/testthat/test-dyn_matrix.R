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

result <- dyn_matrix(wolf_samples$AnimalRef, wolf_samples$Date, seasons$start, seasons$end
)

test_that("dyn_matrix returns expected dimensions", {


  # Check if the dimensions of the matrix are correct
  expect_equal(dim(result)[1], length(seasons$start) + 1)
  expect_equal(dim(result)[2], length(seasons$start) + 1)
})



test_that("dyn_matrix returns expected values for specific inputs", {


  # Check the values in the matrix for specific seasons
  expect_equal(result[1, 1], 13)
  expect_equal(result[2, 2], 4)
  expect_equal(result[3, 3], 17)

  # Check total captures
  expect_equal(result[1, length(seasons$start) + 1], 13)
  expect_equal(result[2, length(seasons$start) + 1], 12)
  expect_equal(result[3, length(seasons$start) + 1], 25)

  # Check total skipped
  expect_equal(result[length(seasons$start) + 1, 1], 0)
  expect_equal(result[length(seasons$start) + 1, 2], 2)
  expect_equal(result[length(seasons$start) + 1, 3], 0)
})

