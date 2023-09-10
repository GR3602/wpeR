test_that("each row is one individual, all individuals are included", {
  expect_identical(nrow(anim_timespan(wolf_samples$AnimalRef, wolf_samples$Date, wolf_samples$SType, dead = "Tissue")),
                   length(unique(wolf_samples$AnimalRef)))
})

test_that("FirstSeen and LastSeen are min and max sample dates of particular individual", {
  subset <- wolf_samples[wolf_samples$AnimalRef == "M2ALK",]
  expect_identical(anim_timespan(subset$AnimalRef, subset$Date, subset$SType, dead = "Tissue")[c(2,3)],
                   data.frame(FirstSeen = min(subset$Date),LastSeen = max(subset$Date)))
})

test_that("No duplicated individuals", {
  result <- anim_timespan(wolf_samples$AnimalRef, wolf_samples$Date, wolf_samples$SType, dead = "Tissue")
  expect_true(anyDuplicated.default(result$ID) == 0)
})

test_that("FirstSeen is in Date format", {
  result <- anim_timespan(wolf_samples$AnimalRef, wolf_samples$Date, wolf_samples$SType, dead = "Tissue")
  expect_true(inherits(result$FirstSeen, "Date"))
})

test_that("LastSeen is in Date format", {
  result <- anim_timespan(wolf_samples$AnimalRef, wolf_samples$Date, wolf_samples$SType, dead = "Tissue")
  expect_true(inherits(result$LastSeen, "Date"))
})




