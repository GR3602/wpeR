animal_ts <- anim_timespan(wolf_samples$AnimalRef,
                           wolf_samples$Date,
                           wolf_samples$SType,
                           dead = c("Tissue"))
sampledata <- merge(wolf_samples, animal_ts, by.x = "AnimalRef", by.y = "ID", all.x = TRUE )
path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
ped_colony <- get_colony(path, sampledata, rm_obsolete_parents = TRUE, out = "FamAgg")
org_tables <- org_fams(ped_colony, sampledata, output = "both")
pt <- plot_table(plot.fams = "all",
           org_tables$fams,
           org_tables$ped,
           sampledata,
           deadSample = c("Tissue", "Decomposing Tissue", "Blood"))


test_that("Number of individuals is the same as in pedigree",{
  expect_equal(length(unique(pt$AnimalRef)), nrow(ped_colony)  )
})

test_that("All samples are included",{
  expect_equal(length(unique(pt$Sample)), nrow(wolf_samples))
})

test_that("min date of animal has first sample flag", {
  min.sample <- aggregate(Date~AnimalRef,data = pt, FUN = min)
  first.flag <- distinct(pt[pt$first_sample == TRUE,c(2,4)])
  expect_true(all(min.sample[match(first.flag$AnimalRef, min.sample$AnimalRef),] == first.flag))
  expect_equal(nrow(first.flag), length(unique(pt$AnimalRef)))
})


test_that("max date of animal has last sample flag", {
  max.sample <- aggregate(Date~AnimalRef,data = pt, FUN = max)
  last.flag <- distinct(pt[pt$last_sample == TRUE,c(2,4)])
  expect_true(all(max.sample[match(last.flag$AnimalRef, max.sample$AnimalRef),] == last.flag))
  expect_equal(nrow(last.flag), length(unique(pt$AnimalRef)))
})

test_that("all reproductive animals are flagged", {
  expect_true(all(unique(pt$AnimalRef[pt$rep==TRUE]) %in% unlist(org_tables$fams[c(2,3)], use.names = FALSE)))
})



