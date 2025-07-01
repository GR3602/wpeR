path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
colony <- utils::read.table(paste(path,".BestConfig_Ordered",sep=""),sep="",header=T,fill=T,comment.char = "?")
colony$FatherID[grep("\\*", colony$FatherID)] <- NA
colony$MotherID[grep("\\#", colony$MotherID)] <- NA
#since ClusterIndex is not needed for get_ped() have to subset that is why colony[1:3]
colony <- colony[1:3]

#also have to remove cluster index column to compare for get_colony() that is why [-1][-4]
test_that("identcal output to get_colony", {
  expect_identical(get_ped(colony,wolf_samples, out = "FamAgg"), get_colony(path, wolf_samples, out = "FamAgg")[-1])
  expect_identical(get_ped(colony,wolf_samples, out = "kinship2"), get_colony(path, wolf_samples, out = "kinship2")[-5])
  expect_identical(get_ped(colony,wolf_samples, out = "pedtools"), get_colony(path, wolf_samples, out = "pedtools"))
})

test_that("error check", {
  expect_error(get_ped(colony,wolf_samples, out = "something"))
})
