animal_ts <- anim_timespan(wolf_samples$AnimalRef,
                           wolf_samples$Date,
                           wolf_samples$SType == c("Tissue"))
sampledata <- merge(wolf_samples, animal_ts, by.x = "AnimalRef", by.y = "ID", all.x = TRUE )
path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
ped_colony <- get_colony(path, sampledata, rm_obsolete_parents = TRUE, out = "FamAgg")
org_tables <- org_fams(ped_colony, sampledata, output = "both")
pt <- plot_table(plot_fams = "all",
           all_fams = org_tables$fams,
           ped = org_tables$ped,
           sampledata = sampledata)


test_that("Number of individuals is the same as in pedigree",{
  expect_equal(length(unique(pt$AnimalRef)), nrow(ped_colony)  )
})

test_that("All samples are included",{
  expect_equal(length(unique(pt$Sample)), nrow(wolf_samples))
})

test_that("min date of animal has first sample flag", {
  min.sample <- aggregate(Date~AnimalRef,data = pt, FUN = min)
  first.flag <- dplyr::distinct(pt[pt$first_sample == TRUE,c(2,4)])
  expect_true(all(min.sample[match(first.flag$AnimalRef, min.sample$AnimalRef),] == first.flag))
  expect_equal(nrow(first.flag), length(unique(pt$AnimalRef)))
})


test_that("max date of animal has last sample flag", {
  max.sample <- aggregate(Date~AnimalRef,data = pt, FUN = max)
  last.flag <- dplyr::distinct(pt[pt$last_sample == TRUE,c(2,4)])
  expect_true(all(max.sample[match(last.flag$AnimalRef, max.sample$AnimalRef),] == last.flag))
  expect_equal(nrow(last.flag), length(unique(pt$AnimalRef)))
})

test_that("all reproductive animals are flagged", {
  expect_true(all(unique(pt$AnimalRef[pt$rep==TRUE]) %in% unlist(org_tables$fams[c(2,3)], use.names = FALSE)))
})

test_that("plot_indivs filters families correctly", {
  # Pick an individual
  indiv <- org_tables$fams$mother[1]

  pt_indiv <- plot_table(all_fams = org_tables$fams,
                         plot_indivs = indiv,
                         ped = org_tables$ped,
                         sampledata = sampledata,
                         )

  # Find all families she is associated with (as mother, father, or offspring)
  all_associated_fams <- unique(c(
      org_tables$fams$FamID[org_tables$fams$father == indiv | org_tables$fams$mother == indiv],
      org_tables$ped$FamID[org_tables$ped$id == indiv]
  ))

  # Check if all families in pt_indiv are associated with this individual
  expect_true(all(unique(pt_indiv$FamID) %in% all_associated_fams))
})
