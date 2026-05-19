animal_ts <- anim_timespan(wolf_samples$AnimalRef,
                           wolf_samples$Date,
                           wolf_samples$SType == c("Tissue")
)
sampledata <- merge(wolf_samples, animal_ts, by.x = "AnimalRef", by.y = "ID", all.x = TRUE)
path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
ped_colony <- get_colony(path, sampledata, rm_obsolete_parents = TRUE, out = "FamAgg")

org_tables <- org_fams(ped_colony, sampledata, output = "both")
pt <- plot_table(plot_fams = 1,
                 all_fams = org_tables$fams,
                 ped = org_tables$ped,
                 sampledata = sampledata
)



test_that("result is list with coorect length", {
result <- ped_spatial(plottable = pt,
             na.rm = TRUE,
             output = "list")

expect_type(result, "list")
expect_equal(length(result), 14)
})



