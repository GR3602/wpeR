path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")

test_that("Expect error when out parameter not correctly specified", {
  expect_error(get_colony(path, wolf_samples, out = "something"))
})

test_that("Output kinship2 structure", {
  ks2 <- c("id", "dadid", "momid", "sex", "famid")
  expect_identical(names(get_colony(path, wolf_samples, out = "kinship2")), ks2)
})

test_that("Output famagg structure", {
  famagg <- c("family", "id", "father", "mother", "sex")
  expect_identical(names(get_colony(path, wolf_samples, out = "FamAgg")), famagg)
})

test_that("Output pedtools structure", {
  pt <- c("id", "fid", "mid", "sex")
  expect_identical(names(get_colony(path, wolf_samples, out = "pedtools")), pt)
})

test_that("Pedigree stays the same", {
  colony <- utils::read.table(paste(path,".BestConfig_Ordered",sep=""),sep="",header=T,fill=T,comment.char = "?")
  colony$FatherID[grep("\\*", colony$FatherID)] <- NA
  colony$MotherID[grep("\\#", colony$MotherID)] <- NA
  result <- get_colony(path, wolf_samples,  out = "table")
  expect_identical(colony[c(1:3)], result[c(1:3)])
})

test_that("Correct sex", {
  result <- get_colony(path, wolf_samples,  out = "table")[c(1,7)]
  samples <- wolf_samples[wolf_samples$Sample %in% result$OffspringID,c(3,4)]
  expect_identical(samples$GeneticSex, result$sex[match(samples$AnimalRef, result$OffspringID)])
})

