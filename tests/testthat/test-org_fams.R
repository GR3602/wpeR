animal_ts <- anim_timespan(wolf_samples$AnimalRef,
                           wolf_samples$Date,
                           wolf_samples$SType,
                           dead = c("Tissue"))
sampledata <- merge(wolf_samples, animal_ts, by.x = "AnimalRef", by.y = "ID", all.x = TRUE )
path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
ped_colony <- get_colony(path, sampledata, rm_obsolete_parents = TRUE, out = "FamAgg")

test_that("Expect error when output parameter uncorrectly specified", {
  expect_error(org_fams(ped_colony, sampledata, out = "something"))
})

test_that("ped data frame has correct columns", {
  result <- org_fams(ped_colony, sampledata, output = "ped")
  expect_s3_class(result, "data.frame")
  colnams <- c("ClusterIndex", "id", "father", "mother", "sex", "parents", "FamID",
               "FirstSeen", "LastSeen", "IsDead", "DadPclust", "MomPclust",
               "polyCluster")
  expect_true(all(colnams %in% names(result)))
})

test_that("fams data frame has correct columns", {
  result <- org_fams(ped_colony, sampledata, output = "fams")
  expect_s3_class(result, "data.frame")
  colnams <- c("parents", "father", "mother", "FamID", "FamStart", "FamEnd",
               "FamDead", "DadPclust", "MomPclust", "polyCluster")
  expect_true(all(colnams %in% names(result)))
})

test_that("number of families is the same in both data frames", {
  result <- org_fams(ped_colony, sampledata, output = "both")
  expect_type(result, "list")
  expect_equal(max(result$ped$FamID), max(result$fams$FamID))
})

test_that("number of families matches distinct parent pairs", {
  result <- org_fams(ped_colony, sampledata, output = "fams")
  distinct_parent_pairs <- length(unique(paste(ped_colony$father, ped_colony$mother, sep = "")))
  expect_equal(nrow(result), distinct_parent_pairs)
})

test_that("all individuals are included in ped data frame", {
  result <- org_fams(ped_colony, sampledata, output = "ped")
  expect_true(all(ped_colony$id %in% result$id))
})

