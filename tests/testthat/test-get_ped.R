path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
colony <- utils::read.table(paste(path,".BestConfig_Ordered",sep=""),sep="",header=T,fill=T,comment.char = "?")
colony$FatherID[grep("\\*", colony$FatherID)] <- NA
colony$MotherID[grep("\\#", colony$MotherID)] <- NA

test_that("identcal output to get_colony", {
  expect_identical(get_ped(colony,wolf_samples, out = "FamAgg"), get_colony(path, wolf_samples, out = "FamAgg"))
  expect_identical(get_ped(colony,wolf_samples, out = "kinship2"), get_colony(path, wolf_samples, out = "kinship2"))
  expect_identical(get_ped(colony,wolf_samples, out = "pedtools"), get_colony(path, wolf_samples, out = "pedtools"))
})

test_that("error check", {
  expect_error(get_ped(colony,wolf_samples, out = "something"))
})
