animal_ts <- anim_timespan(wolf_samples$AnimalRef,
                           wolf_samples$Date,
                           wolf_samples$SType,
                           dead = c("Tissue", "Decomposing Tissue", "Blood"))
sampledata <- merge(wolf_samples, animal_ts, by.x = "AnimalRef", by.y = "ID", all.x = TRUE )
path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
ped_colony <- get_colony(path, sampledata, remove_obsolete_parents = TRUE, out = "FamAgg")
org_tables <- org_fams(ped_colony, sampledata, output = "both")
pt <- plot_table(org_tables$fams,
                 org_tables$fams,
                 org_tables$ped,
                 sampledata,
                 deadSample = c("Tissue", "Decomposing Tissue", "Blood"))
pt$AnimalRef_Fam <- paste0(pt$AnimalRef, "_", pt$FamID)

####ppsList####
test_that("ppsList no limitations", {

  result <- ppsList(pt)

  expect_equal(nrow(result$motherRefs),
               length(unique(pt$AnimalRef_Fam[pt$rep == TRUE & pt$GeneticSex %in% "F"])),
               label = "mother refs"
  )
  expect_equal(nrow(result$motherAll),
               sum(pt$rep == TRUE & pt$GeneticSex %in% "F"),
               label = "mother all"
  )
  expect_equal(nrow(result$fatherRefs),
               length(unique(pt$AnimalRef_Fam[pt$rep == TRUE & pt$GeneticSex %in% "M"])),
               label = "father refs"
  )
  expect_equal(nrow(result$fatherAll),
               sum(pt$rep == TRUE & pt$GeneticSex %in% "M"),
               label = "father all"
  )
  expect_equal(nrow(result$offspringRefs),
               length(unique(pt$AnimalRef[pt$rep == FALSE])),
               label = "offspring refs"
  )
  expect_equal(nrow(result$offspring),
               sum(pt$rep == FALSE),
               label = "all offspring"
  )

})


tlmt <- as.Date(c("2020-01-01", "2020-12-31"))



test_that("ppsList time limits but no flags", {

  result <- ppsList(pt, time.limits =  tlmt)

  expect_equal(nrow(result$motherRefs),
               length(unique(pt$AnimalRef_Fam[pt$rep == TRUE & pt$GeneticSex %in% "F"])),
               label = "mother refs"
  )
  expect_equal(nrow(result$motherAll),
               sum(pt$rep == TRUE & pt$GeneticSex %in% "F"),
               label = "mother all"
  )
  expect_equal(nrow(result$fatherRefs),
               length(unique(pt$AnimalRef_Fam[pt$rep == TRUE & pt$GeneticSex %in% "M"])),
               label = "father refs"
  )
  expect_equal(nrow(result$fatherAll),
               sum(pt$rep == TRUE & pt$GeneticSex %in% "M"),
               label = "father all"
  )
  expect_equal(nrow(result$offspringRefs),
               length(unique(pt$AnimalRef[pt$rep == FALSE])),
               label = "offspring refs"
  )
  expect_equal(nrow(result$offspring),
               sum(pt$rep == FALSE),
               label = "all offspring"
  )

})


#ped_satplot(pt$Date, pt$AnimalRef, pt$plottingID, pt$GeneticSex, pt$FamID, pt$polyCluster, pt$isPolygamous, pt$rep, pt$later_rep, pt$dead)

as.Date("2020-01-01")
as.Date("2020-12-31")


#function (plottable,
#          time.limits = c(as.Date("1900-01-01"), as.Date("2100-01-01")),
#          na.rm = T,
#          time.limit.rep = F,
#          time.limit.offspring = F)


test_that("ppsList time limit for rep animals reference", {
  result <- ppsList(pt, time.limits =  tlmt,
                    time.limit.rep = T)

  reps <- pt[pt$rep == TRUE &
               pt$Date >= tlmt[1] &
               pt$Date <= tlmt[2],]

  expect_equal(nrow(result$motherRefs),
               length(unique(reps$AnimalRef_Fam[reps$rep == TRUE & reps$GeneticSex %in% "F"])),
               label = "mother refs"
  )
  expect_equal(nrow(result$motherAll),
               sum(pt$rep == TRUE & pt$GeneticSex %in% "F"),
               label = "mother all"
  )
  expect_equal(nrow(result$fatherRefs),
               length(unique(reps$AnimalRef_Fam[reps$rep == TRUE & reps$GeneticSex %in% "M"])),
               label = "father refs"
  )
  expect_equal(nrow(result$fatherAll),
               sum(pt$rep == TRUE & pt$GeneticSex %in% "M"),
               label = "father all"
  )
  expect_equal(nrow(result$offspringRefs),
               length(unique(pt$AnimalRef[pt$rep == FALSE])),
               label = "offspring refs"
  )
  expect_equal(nrow(result$offspring),
               sum(pt$rep == FALSE),
               label = "all offspring"
  )


})


test_that("ppsList time limit for offspring", {
  result <- ppsList(pt, time.limits =  tlmt,
                    time.limit.rep = F,
                    time.limit.offspring = T)

  offs <- pt[pt$rep == FALSE &
               pt$Date >= tlmt[1] &
               pt$Date <= tlmt[2],]

  expect_equal(nrow(result$motherRefs),
               length(unique(pt$AnimalRef_Fam[pt$rep == TRUE & pt$GeneticSex %in% "F"])),
               label = "mother refs"
  )
  expect_equal(nrow(result$motherAll),
               sum(pt$rep == TRUE & pt$GeneticSex %in% "F"),
               label = "mother all"
  )
  expect_equal(nrow(result$fatherRefs),
               length(unique(pt$AnimalRef_Fam[pt$rep == TRUE & pt$GeneticSex %in% "M"])),
               label = "father refs"
  )
  expect_equal(nrow(result$fatherAll),
               sum(pt$rep == TRUE & pt$GeneticSex %in% "M"),
               label = "father all"
  )
  expect_equal(nrow(result$offspringRefs),
               length(unique(offs$AnimalRef[offs$rep == FALSE])),
               label = "offspring refs"
  )
  expect_equal(nrow(result$offspring),
               sum(pt$rep == FALSE),
               label = "all offspring"
  )


})


test_that("ppsList moves timelimits", {
  result <- ppsList(pt, time.limits =  tlmt,
                    time.limit.rep = F,
                    time.limit.offspring = F,
                    time.limit.moves = T)

  moves <- pt[pt$Date >= tlmt[1] &
             pt$Date <= tlmt[2],]

  expect_equal(nrow(result$motherRefs),
               length(unique(pt$AnimalRef_Fam[pt$rep == TRUE & pt$GeneticSex %in% "F"])),
               label = "mother refs"
  )
  expect_equal(nrow(result$motherAll),
               sum(moves$rep == TRUE & moves$GeneticSex %in% "F"),
               label = "mother all"
  )
  expect_equal(nrow(result$fatherRefs),
               length(unique(pt$AnimalRef_Fam[pt$rep == TRUE & pt$GeneticSex %in% "M"])),
               label = "father refs"
  )
  expect_equal(nrow(result$fatherAll),
               sum(moves$rep == TRUE & moves$GeneticSex %in% "M"),
               label = "father all"
  )
  expect_equal(nrow(result$offspringRefs),
               length(unique(pt$AnimalRef[pt$rep == FALSE])),
               label = "offspring refs"
  )
  expect_equal(nrow(result$offspring),
               sum(moves$rep == FALSE),
               label = "all offspring"
  )

})

test_that("ppsList all flags", {
  result <- ppsList(pt, time.limits =  tlmt,
                    time.limit.rep = T,
                    time.limit.offspring = T,
                    time.limit.moves = T)

  pt <- pt[pt$Date >= tlmt[1] &
               pt$Date <= tlmt[2],]

  expect_equal(nrow(result$motherRefs),
               length(unique(pt$AnimalRef_Fam[pt$rep == TRUE & pt$GeneticSex %in% "F"])),
               label = "mother refs"
  )
  expect_equal(nrow(result$motherAll),
               sum(pt$rep == TRUE & pt$GeneticSex %in% "F"),
               label = "mother all"
  )
  expect_equal(nrow(result$fatherRefs),
               length(unique(pt$AnimalRef_Fam[pt$rep == TRUE & pt$GeneticSex %in% "M"])),
               label = "father refs"
  )
  expect_equal(nrow(result$fatherAll),
               sum(pt$rep == TRUE & pt$GeneticSex %in% "M"),
               label = "father all"
  )
  expect_equal(nrow(result$offspringRefs),
               length(unique(pt$AnimalRef[pt$rep == FALSE])),
               label = "offspring refs"
  )
  expect_equal(nrow(result$offspring),
               sum(pt$rep == FALSE),
               label = "all offspring"
  )

})

####ppsParLines####

list_data <- ppsList(pt)
test_that("ppsParLines no. paternity and maternity",{

  list_data <- ppsList(pt)

  result <- ppsParLines(list_data)

expect_equal(
  nrow(result$maternityLines),
  nrow(org_tables$ped[!is.na(org_tables$ped$mother),]),
  label =  "maternity"
)

expect_equal(
  nrow(result$paternityLines),
  nrow(org_tables$ped[!is.na(org_tables$ped$father),]),
  label = "paternity"
)
})


list_data <- ppsList(pt,
                     time.limits =  tlmt,
                     time.limit.rep = T,
                     time.limit.offspring = T,
                     time.limit.moves = T)

####ppsMvPoints####
test_that("ppsMvpoint equal to ppsList", {

  result <- ppsMvPoints(list_data)

  expect_equal(nrow(list_data$motherAll), nrow(result$motherMovePoints),
               label = "mother move")
  expect_equal(nrow(list_data$fatherAll), nrow(result$fatherMovePoints),
               label = "father move")
  expect_equal(nrow(list_data$offspring), nrow(result$offspringMovePoints),
               label = "offspring move")

})

####ppsRefPoints####
test_that("ppsRefPoints equal to ppsList", {

  result <- ppsRefPoints(list_data)

  expect_equal(nrow(list_data$motherRefs), nrow(result$motherRpoints),
               label = "mother refs")
  expect_equal(nrow(list_data$fatherRefs), nrow(result$fatherRpoints),
               label = "father refs")
  expect_equal(nrow(list_data$offspringRefs), nrow(result$offspringRpoints),
               label = "offspring refs")

})


####ppsMvLines####
test_that("ppsMvLines equal to ppsList", {

  result <- ppsMvLines(list_data)

  no_lines <- aggregate(Sample ~  AnimalRef, data = list_data$motherAll, function(x) length(x)-1)

  expect_equal(nrow(result$motherMoveLines), sum(no_lines$Sample != 0),
               label = "mother move lines")

  no_lines <- aggregate(Sample ~  AnimalRef, data = list_data$fatherAll, function(x) length(x)-1)

  expect_equal(nrow(result$fatherMoveLines), sum(no_lines$Sample != 0),
               label = "father move lines")

  no_lines <- aggregate(Sample ~  AnimalRef, data = list_data$offspring, function(x) length(x)-1)

  expect_equal(nrow(result$offspringMoveLines), sum(no_lines$Sample != 0),
               label = "offspring move lines")

})


