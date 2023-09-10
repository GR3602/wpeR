animal_ts <- anim_timespan(wolf_samples$AnimalRef,
                           wolf_samples$Date,
                           wolf_samples$SType,
                           dead = c("Tissue"))
sampledata <- merge(wolf_samples, animal_ts, by.x = "AnimalRef", by.y = "ID", all.x = TRUE )
path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
ped_colony <- get_colony(path, sampledata, rm_obsolete_parents = TRUE, out = "FamAgg")
org_tables <- org_fams(ped_colony, sampledata, output = "both")
pt <- plot_table(org_tables$fams,
                 org_tables$fams,
                 org_tables$ped,
                 sampledata,
                 deadSample = c("Tissue", "Decomposing Tissue", "Blood"))

test_that("Is plot produced", {
  test.plot <- ped_satplot(pt$Date, pt$AnimalRef, pt$plottingID, pt$GeneticSex, pt$FamID,
              pt$polyCluster, pt$isPolygamous, pt$rep, pt$later_rep, pt$dead)

  expect_true(!is.null(test.plot))
})
