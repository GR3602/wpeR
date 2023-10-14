#### ppsList####
# Prepare Pedigree Spatial data preparation
#
# Takes data from fam_table function and arranges it so that
# georeferenced data for mothers, fathers and offspring can be created
#
# @param plottable output of fam_table function
# @param time.limits time window for movement and offspring reference samples data
# @param na.rm remove samples with missing coordinates
# @param time.limit.rep do time limits apply to reference samples of reproductive animals?
# @param time.limit.offspring do time limits apply to reference samples of offspring?
#

ppsList <- function(plottable,
                    time.limits = c(as.Date("1900-01-01"), as.Date("2100-01-01")),
                    na.rm = TRUE,
                    time.limit.rep = FALSE,
                    time.limit.offspring = FALSE,
                    time.limit.moves = FALSE) {
  # remove samples with missing coordinates
  if (sum(is.na(plottable$lat)) > 0 |
      sum(is.na(plottable$lng)) > 0 |
      sum(is.na(plottable$Date) > 0)) {
    if (na.rm == TRUE) {
      plottable <- plottable[!(is.na(plottable$lat) |
                                 is.na(plottable$lng) |
                                 is.na(plottable$Date)), ]
    } else {
      stop("Error, NA's in coordinates and/or dates. Use na.rm=TRUE to remove internally.")
    }
  }
  # sort plottable by individual & date
  plottable <- plottable[order(plottable$AnimalRef, plottable$Date), ]


  # find first&last sample within time.limits
  plottable$first_sample <- rep(FALSE, nrow(plottable))
  plottable$last_sample <- rep(FALSE, nrow(plottable))

  for (i in seq_len(nrow(plottable))) {
    curAnimal <- plottable[plottable$AnimalRef == plottable$AnimalRef[i] &
      plottable$Date <= time.limits[2] &
      plottable$Date >= time.limits[1], ]
    curAnimal <- curAnimal[order(curAnimal$Date), ]

    # of no samples are within the time window, use actual first and last sample
    if (nrow(curAnimal) == 0) {
      curAnimal <- plottable[plottable$AnimalRef == plottable$AnimalRef[i], ]
      curAnimal <- curAnimal[order(curAnimal$Date), ]
    }

    if (plottable$Sample[i] == curAnimal$Sample[1]) plottable$first_sample[i] <- TRUE
    if (plottable$Sample[i] == curAnimal$Sample[nrow(curAnimal)]) plottable$last_sample[i] <- TRUE
  } # for

  # First create pointdata
  # Use time window. Limit entire dataset, or only offspring (later)

  # all samples reproductive animal
  reps <- plottable[plottable$rep == TRUE, ]

  # reference samples reproductive animals = last sample
  # last sample of reproductive animals as spatial ref
  repRefs <- reps[reps$last_sample == TRUE, ]



  # use time limit window also for reproductive animals?
  if (time.limit.rep == TRUE) {
    repRefs <- repRefs[repRefs$Date <= time.limits[2] & repRefs$LastSeen >= time.limits[1], ]
  }

  fatherRefs <- repRefs[repRefs$GeneticSex == "M", ]
  motherRefs <- repRefs[repRefs$GeneticSex == "F", ]



  offspring <- plottable[plottable$rep == FALSE, ]
  if (time.limit.offspring == TRUE) {
    offspringRefs <- offspring[offspring$first_sample == TRUE &
                                 offspring$Date <= time.limits[2] &
                                 offspring$LastSeen >= time.limits[1], ]
  } else {
    offspringRefs <- offspring[offspring$first_sample == TRUE, ]
  }

  # limit

  # movement limit, limits all data

  if (time.limit.moves == TRUE) {
    reps <- reps[reps$Date >= time.limits[1] & reps$Date <= time.limits[2], ]
    offspring <- offspring[offspring$Date >= time.limits[1] &
                             offspring$Date <= time.limits[2], ]
  }


  fatherAll <- reps[reps$GeneticSex == "M", ]
  motherAll <- reps[reps$GeneticSex == "F", ]

  return(list(
    motherAll = motherAll,
    motherRefs = motherRefs,
    fatherAll = fatherAll,
    fatherRefs = fatherRefs,
    offspring = offspring,
    offspringRefs = offspringRefs,
    repRefs = repRefs
  ))
} # end ppsList



#### ppsParLines####
# Prepare Pedigree Spatial parent - offspring lines
#
# Creates sf dataframe with lines connecting reference samples of
# offspring to their parents
#
# @param ppsData output list of ppsList function
#
# @import sf

ppsParLines <- function(ppsData) {

  # Create Paternity / Maternity Lines (sf)
  maternityLines <- NULL
  paternityLines <- NULL

  dadcount <- 1
  momcount <- 1

  repRefs <- ppsData$repRefs
  offspringRefs <- ppsData$offspringRefs

  for (i in seq_len(nrow(offspringRefs))) {
    # offspring record
    if (nrow(offspringRefs) > 0) {
      child <- offspringRefs[i, ]
    } else {
      child <- NA
      stop("No offspring included, check data and/or time.limits parameter.")
    }

    father <- repRefs[repRefs$FamID == offspringRefs$FamID[i] &
                        repRefs$GeneticSex == "M", ]
    mother <- repRefs[repRefs$FamID == offspringRefs$FamID[i] &
                        repRefs$GeneticSex == "F", ]

    if (nrow(father) > 0 & nrow(child) > 0) {
      paternityLines[[dadcount]] <- rbind(father[, c("lat", "lng")], child[, c("lat", "lng")])
      paternityLines[[dadcount]] <- st_as_sf(paternityLines[[dadcount]],
                                             coords = c("lng", "lat"),
                                             crs = 4326)
      paternityLines[[dadcount]] <- st_combine(paternityLines[[dadcount]])
      paternityLines[[dadcount]] <- st_cast(paternityLines[[dadcount]], "LINESTRING")
      paternityLines[[dadcount]] <- data.frame(
        ID = dadcount,
        pair = i,
        fam = offspringRefs$FamID[i],
        plyClust = offspringRefs$polyCluster[i],
        relation = "paternity",
        child = child$AnimalRef,
        parent = father$AnimalRef,
        geometry = st_geometry(paternityLines[[dadcount]])
      )
      paternityLines[[dadcount]] <- st_as_sf(paternityLines[[dadcount]],
                                             sf_column_name = "geometry")
      dadcount <- dadcount + 1
    }

    if (nrow(mother) > 0 & nrow(child) > 0) {
      maternityLines[[momcount]] <- rbind(mother[, c("lat", "lng")], child[, c("lat", "lng")])
      maternityLines[[momcount]] <- st_as_sf(maternityLines[[momcount]],
                                             coords = c("lng", "lat"),
                                             crs = 4326)
      maternityLines[[momcount]] <- st_combine(maternityLines[[momcount]])
      maternityLines[[momcount]] <- st_cast(maternityLines[[momcount]], "LINESTRING")
      maternityLines[[momcount]] <- data.frame(
        ID = momcount,
        pair = i,
        fam = offspringRefs$FamID[i],
        plyClust = offspringRefs$polyCluster[i],
        relation = "maternity",
        child = child$AnimalRef,
        parent = mother$AnimalRef,
        geometry = st_geometry(maternityLines[[momcount]])
      )
      maternityLines[[momcount]] <- st_as_sf(maternityLines[[momcount]],
                                             sf_column_name = "geometry")
      momcount <- momcount + 1
    }
  } # for end

  #hashed 2023-10-14 it looks like that the if not needed
  #if (nrow(father) > 0 & nrow(child) > 0) {
    paternityLines <- dplyr::bind_rows(paternityLines)
  #}

  #if (nrow(father) > 0 & nrow(child) > 0) {
    maternityLines <- dplyr::bind_rows(maternityLines)
  #}

  return(list(
    maternityLines = maternityLines,
    paternityLines = paternityLines
  ))
} # end ppsParLines

#### ppsMvPoints####
# Prepare Pedigree Spatial - Animal Movement Points
#
# @param ppsData output list of ppsList function
# @param time.limits time window for movement and offspring reference samples data
# @param time.limit.moves time limit also movement data?
#
# @import sf
#


ppsMvPoints <- function(ppsData) {
  # columns to be removed for movelines, movepoints and refpoints to remove spatial data duplication
  dupCols <- c("plottingID", "FamID", "polyCluster")

  #### Offspring####

  offspringRefs <- ppsData$offspringRefs
  offspring <- ppsData$offspring

  offspringMoveData <- offspring[offspring$AnimalRef %in% offspringRefs$AnimalRef, ]

  if (nrow(offspringMoveData) > 0) {
    offspringMovePoints <- st_as_sf(offspringMoveData,
                                    coords = c("lng", "lat"),
                                    crs = 4326)
  } else {
    offspringMovePoints <- st_sf(1, st_sfc(st_point()))
    warning("No offspring samples within selected time period.
                  Creating empty sf data frame -> offspringMovePoints!")
  }

  #### Father####

  fatherRefs <- ppsData$fatherRefs
  fatherAll <- ppsData$fatherAll

  fatherMoveData <- fatherAll[fatherAll$AnimalRef %in% fatherRefs$AnimalRef, ]
  # drop duplication columns
  fatherMoveData <- fatherMoveData[, !names(fatherMoveData) %in% dupCols]
  fatherMoveData <- unique(fatherMoveData)

  if (nrow(fatherMoveData) > 0) {
    fatherMovePoints <- st_as_sf(fatherMoveData,
                                 coords = c("lng", "lat"),
                                 crs = 4326)
  } else {
    fatherMovePoints <- st_sf(1, st_sfc(st_point()))
    warning("No father samples within selected time period.
                  Creating empty sf data frame -> fatherMovePoints!")
  }

  #### Mother####

  motherRefs <- ppsData$motherRefs
  motherAll <- ppsData$motherAll

  motherMoveData <- motherAll[motherAll$AnimalRef %in% motherRefs$AnimalRef, ]
  # drop duplication columns
  motherMoveData <- motherMoveData[, !names(motherMoveData) %in% dupCols]
  motherMoveData <- unique(motherMoveData)

  if (nrow(motherMoveData) > 0) {
    motherMovePoints <- st_as_sf(motherMoveData,
                                 coords = c("lng", "lat"),
                                 crs = 4326)
  } else {
    motherMovePoints <- st_sf(1, st_sfc(st_point()))
    warning("No mother samples within selected time period.
                  Creating empty sf data frame -> motherMovePoints!")
  }

  return(list(
    offspringMovePoints = offspringMovePoints,
    fatherMovePoints = fatherMovePoints,
    motherMovePoints = motherMovePoints
  ))
} # end ppsMvPoints

#### ppsRefPoints####
# Prepare Pedigree Spatial - Animal Reference Points
#
# @param ppsData output list of ppsList function
#
# @import sf


ppsRefPoints <- function(ppsData) {
  # columns to be removed for movelines, movepoints and refpoints to remove spatial data duplication
  dupCols <- c("plottingID", "FamID", "polyCluster")

  #### Offspring####

  offspringRefs <- ppsData$offspringRefs

  if (nrow(offspringRefs) > 0) {
    offspringRpoints <- st_as_sf(offspringRefs,
                                 coords = c("lng", "lat"),
                                 crs = 4326)
  } else {
    offspringRpoints <- st_sf(1, st_sfc(st_point()))
    warning("No offspring samples within selected time period.
                  Creating empty sf data frame -> offspringRpoints!")
  }

  #### Father####

  fatherRefs <- ppsData$fatherRefs

  # drop duplication columns
  fatherRefs <- fatherRefs[, !names(fatherRefs) %in% dupCols]

  if (nrow(fatherRefs) > 0) {
    fatherRpoints <- st_as_sf(fatherRefs,
                              coords = c("lng", "lat"),
                              crs = 4326)
  } else {
    fatherRpoints <- st_sf(1, st_sfc(st_point()))
    warning("No father samples within selected time period.
                  Creating empty sf data frame -> fatherRpoints!")
  }

  #### Mother####

  motherRefs <- ppsData$motherRefs

  # drop duplication columns
  motherRefs <- motherRefs[, !names(motherRefs) %in% dupCols]
  motherRefs <- unique(motherRefs)

  if (nrow(motherRefs) > 0) {
    motherRpoints <- st_as_sf(motherRefs,
                              coords = c("lng", "lat"),
                              crs = 4326)
  } else {
    motherRpoints <- st_sf(1, st_sfc(st_point()))
    warning("No father samples within selected time period.
                  Creating empty sf data frame -> motherRpoints!")
  }

  return(list(
    offspringRpoints = offspringRpoints,
    fatherRpoints = fatherRpoints,
    motherRpoints = motherRpoints
  ))
} # end ppsRefPoints

#### ppsMvLines####
# Prepare Pedigree Spatial - Animal Movement Lines
#
# @param ppsData output list of ppsList function
# @param plottable output of the fam_table function
# @param time.limits time window for movement and offspring reference samples data
# @param time.limit.moves time limit also movement data?
#
# @import sf
#

ppsMvLines <- function(ppsData) {
  # columns to be removed for movelines, movepoints and refpoints to remove spatial data duplication
  dupCols <- c("plottingID", "FamID", "polyCluster")


  #### Offspring####

  indcount <- 1
  individualLines <- NULL
  offspringMoveLines <- NULL

  offspringRefs <- ppsData$offspringRefs
  offspring <- ppsData$offspring

  for (i in seq_len(nrow(offspringRefs))) {
    if (nrow(offspringRefs) == 0) break

    individual <- offspring[offspring$AnimalRef == offspringRefs$AnimalRef[i], !names(offspring) %in% dupCols]

    if (nrow(individual) == 0) next

    individual <- unique(individual)

    if (nrow(individual) >= 2) {
      individualLines[[indcount]] <- st_geometry(st_as_sf(individual,
                                                          coords = c("lng", "lat"),
                                                          crs = 4326))
      individualLines[[indcount]] <- st_cast(st_union(individualLines[[indcount]]),
                                             "LINESTRING")
      offspringMoveLines <- rbind(offspringMoveLines, data.frame(
        ID = indcount,
        AnimalID = offspringRefs$AnimalRef[i],
        fam = offspringRefs$FamID[i],
        plyClust = offspringRefs$polyCluster[i],
        geometry = st_geometry(individualLines[[indcount]])
      ))


      indcount <- indcount + 1
    }
  }

  # hashed are options to deal with examples where there is just one sample
  # option1
  # no df created -> not working properly 2023-05-20
  # if (nrow(individual) >= 2) {
  #  offspringMoveLines = st_as_sf(offspringMoveLines, sf_column_name = "geometry")
  # }

  # option2
  # empty file created
  if (length(individualLines) > 0) {
    offspringMoveLines <- st_as_sf(offspringMoveLines, sf_column_name = "geometry")
  } else {
    offspringMoveLines <- st_sf(1, st_sfc(st_linestring())) # dummy... empty object
    warning("None of the offspring has two samples needed to create a linestring.
            Creating empty sf data frame -> offspringMoveLines!")
  }



  #### Father####

  indcount <- 1
  individualLines <- NULL
  fatherMoveLines <- NULL

  fatherRefs <- ppsData$fatherRef
  fatherAll <- ppsData$fatherAll

  for (i in seq_len(nrow(fatherRefs))) {
    if (nrow(fatherRefs) == 0) break

    individual <- fatherAll[fatherAll$AnimalRef == fatherRefs$AnimalRef[i], !names(fatherAll) %in% dupCols]

    if (nrow(individual) == 0) next

    individual <- unique(individual)

    if (nrow(individual) >= 2) {
      individualLines[[indcount]] <- st_geometry(st_as_sf(individual,
                                                          coords = c("lng", "lat"),
                                                          crs = 4326))
      individualLines[[indcount]] <- st_cast(st_union(individualLines[[indcount]]),
                                             "LINESTRING")
      fatherMoveLines <- rbind(fatherMoveLines, data.frame(
        ID = indcount,
        AnimalID = fatherRefs$AnimalRef[i],
        fam = fatherRefs$FamID[i],
        plyClust = fatherRefs$polyCluster[i],
        geometry = st_geometry(individualLines[[indcount]])
      ))


      indcount <- indcount + 1
    }
  }

  # if (nrow(individual) >= 2) {
  #  fatherMoveLines= st_as_sf(fatherMoveLines, sf_column_name = "geometry")
  # }

  # option2
  # empty file created
  if (length(individualLines) > 0) {
    fatherMoveLines <- st_as_sf(fatherMoveLines, sf_column_name = "geometry")
  } else {
    fatherMoveLines <- st_sf(1, st_sfc(st_linestring())) # dummy... empty object
    warning("None of the fathers has two samples needed to create a linestring.
            Creating empty sf data frame -> fatherMoveLines!")
  }


  #### Mother####

  indcount <- 1
  individualLines <- NULL
  motherMoveLines <- NULL

  motherRefs <- ppsData$motherRefs
  motherAll <- ppsData$motherAll

  for (i in seq_len(nrow(motherRefs))) {
    if (nrow(motherRefs) == 0) break

    individual <- motherAll[motherAll$AnimalRef == motherRefs$AnimalRef[i], !names(motherAll) %in% dupCols]

    if (nrow(individual) == 0) next

    individual <- unique(individual)

    if (nrow(individual) >= 2) {
      individualLines[[indcount]] <- st_geometry(st_as_sf(individual,
                                                          coords = c("lng", "lat"),
                                                          crs = 4326))
      individualLines[[indcount]] <- st_cast(st_union(individualLines[[indcount]]),
                                             "LINESTRING")
      motherMoveLines <- rbind(motherMoveLines, data.frame(
        ID = indcount,
        AnimalID = motherRefs$AnimalRef[i],
        fam = motherRefs$FamID[i],
        plyClust = motherRefs$polyCluster[i],
        geometry = st_geometry(individualLines[[indcount]])
      ))


      indcount <- indcount + 1
    }
  }

  # if (nrow(individual) >= 2) {
  #  motherMoveLines= st_as_sf(motherMoveLines, sf_column_name = "geometry")
  # }

  # option2
  if (length(individualLines) > 0) {
    motherMoveLines <- st_as_sf(motherMoveLines, sf_column_name = "geometry")
  } else {
    motherMoveLines <- st_sf(1, st_sfc(st_linestring())) # dummy... empty object
    warning("None of the mothers has two samples needed to create a linestring.
            Creating empty sf data frame -> motherMoveLines!")
  }


  return(list(
    offspringMoveLines = offspringMoveLines,
    fatherMoveLines = fatherMoveLines,
    motherMoveLines = motherMoveLines
  ))
} # end ppsMvLines

#### ppsMvPolygons####
# Prepare Pedigree Spatial - Convex Hull For Animal Points
#
# @param ppsData output list of ppsList function
# @param MvPoints output of ppsMvPoints function
#
# @import sf
#
#
ppsMvPolygons <- function(ppsData, MvPoints) {
  #### Mother####
  # mother MovePolygons
  # mother MovePolygon
  indcount <- 1

  motherMovePoints <- MvPoints$motherMovePoints
  motherRefs <- ppsData$motherRefs

  motherMovePolygons <- NULL
  for (i in seq_len(nrow(motherRefs))) {
    if (nrow(motherRefs) == 0) break
    individualPoints <- motherMovePoints[motherMovePoints$AnimalRef == motherRefs$AnimalRef[i], ]
    if (nrow(individualPoints) >= 3) {
      animal <- unique(individualPoints$AnimalRef)
      individualPoints <- individualPoints[, 2]
      individualPoints <- st_combine(individualPoints)
      individualPolygon <- st_convex_hull(individualPoints)
      individualPolygon <- st_geometry(individualPolygon)
      individualPolygon <- data.frame(
        animal = animal,
        geometry = individualPolygon
      )
      individualPolygon <- st_as_sf(individualPolygon, sf_column_name = "geometry")
      motherMovePolygons <- rbind(motherMovePolygons, individualPolygon)
    } else {
      next
    }
  }

  if (is.null(motherMovePolygons)) {
    motherMovePolygons <- st_sf(1, st_sfc(st_polygon()))
    warning("Not enough mother samples to create at least one polygon.
            Creating empty sf data frame -> motherMovePolygons!")
  }

  #### Father####
  # father MovePolygons
  # father MovePolygon
  indcount <- 1

  fatherMovePoints <- MvPoints$fatherMovePoints
  fatherRefs <- ppsData$fatherRefs

  fatherMovePolygons <- NULL
  for (i in seq_len(nrow(fatherRefs))) {
    # offspring record
    if (nrow(fatherRefs) == 0) break
    individualPoints <- fatherMovePoints[fatherMovePoints$AnimalRef == fatherRefs$AnimalRef[i], ]
    if (nrow(individualPoints) >= 3) {
      animal <- unique(individualPoints$AnimalRef)
      individualPoints <- individualPoints[, 2]
      individualPoints <- st_combine(individualPoints)
      individualPolygon <- st_convex_hull(individualPoints)
      individualPolygon <- st_geometry(individualPolygon)
      individualPolygon <- data.frame(
        animal = animal,
        geometry = individualPolygon
      )
      individualPolygon <- st_as_sf(individualPolygon, sf_column_name = "geometry")
      fatherMovePolygons <- rbind(fatherMovePolygons, individualPolygon)
    } else {
      next
    }
  }

  if (is.null(fatherMovePolygons)) {
    fatherMovePolygons <- st_sf(1, st_sfc(st_polygon()))
    warning("Not enough father samples to create at least one polygon.
            Creating empty sf data frame -> fatherMovePolygons!")
  }



  #### Offspring####
  # offspring MovePolygons
  # offspring MovePolygon
  indcount <- 1

  offspringMovePoints <- MvPoints$offspringMovePoints
  offspringRefs <- ppsData$offspringRefs

  offspringMovePolygons <- NULL
  for (i in seq_len(nrow(offspringRefs))) {
    # offspring record
    if (nrow(offspringRefs) == 0) break
    individualPoints <- offspringMovePoints[offspringMovePoints$AnimalRef == offspringRefs$AnimalRef[i], ]
    if (nrow(individualPoints) >= 3) {
      animal <- unique(individualPoints$AnimalRef)
      individualPoints <- individualPoints[, 2]
      individualPoints <- st_combine(individualPoints)
      individualPolygon <- st_convex_hull(individualPoints)
      individualPolygon <- st_geometry(individualPolygon)
      individualPolygon <- data.frame(
        animal = animal,
        geometry = individualPolygon
      )
      individualPolygon <- st_as_sf(individualPolygon, sf_column_name = "geometry")
      offspringMovePolygons <- rbind(offspringMovePolygons, individualPolygon)
    } else {
      next
    }
  }

  if (is.null(offspringMovePolygons)) {
    offspringMovePolygons <- st_sf(1, st_sfc(st_polygon()))
    warning("Not enough offspring samples to create at least one polygon.
            Creating empty sf data frame -> offspringMovePolygons!")
  }


  return(list(
    motherMovePolygons = motherMovePolygons,
    fatherMovePolygons = fatherMovePolygons,
    offspringMovePolygons = offspringMovePolygons
  ))
} # end ppsMvPolygons



#### ppsFsLines####
# Prepare Pedigree Spatial - Connect Siblings
#
# @param ppsData output list of ppsList function
# @param fullsibdata COLONY2 fullsib data
# @param sibthreshold p value threshold for sibsihp assignment
#
# @import sf

ppsFsLines <- function(ppsData, fullsibdata, sibthreshold = 1) {
  # require(sf)
  # require(dplyr)

  offspringRefs <- ppsData$offspringRefs

  indcount <- 1
  individualLines <- NULL

  for (i in seq_len(nrow(offspringRefs))) {
    # offspring record
    if (nrow(offspringRefs) == 0) break

    allsibs <- fullsibdata[(fullsibdata$OffspringID1 == offspringRefs$AnimalRef[i] |
      fullsibdata$OffspringID2 == offspringRefs$AnimalRef[i]) &
      fullsibdata$Probability >= sibthreshold, ] # get full siblings

    # make a vector of all siblings
    sibgroup <- unique(c(as.character(allsibs[, 1]), as.character(allsibs[, 2])))

    # samples of sibgroup
    sibgroupsamples <- offspringRefs[offspringRefs$AnimalRef %in% sibgroup, ]

    # make lines between all reference samples within a siblings group
    if (nrow(sibgroupsamples) >= 2) {
      for (j in 1:(nrow(sibgroupsamples) - 1)) {
        for (k in (j + 1):nrow(sibgroupsamples)) {
          individualLines[[indcount]] <- sibgroupsamples[c(j, k), c("lng", "lat")]
          individualLines[[indcount]] <- st_as_sf(individualLines[[indcount]],
                                                  coords = c("lng", "lat"),
                                                  crs = 4326)
          individualLines[[indcount]] <- st_combine(individualLines[[indcount]])
          individualLines[[indcount]] <- st_cast(individualLines[[indcount]],
                                                 "LINESTRING")
          individualLines[[indcount]] <- data.frame(
            ID = indcount,
            AnimalID1 = sibgroupsamples$AnimalRef[j],
            AnimalID2 = sibgroupsamples$AnimalRef[k],
            geometry = st_geometry(individualLines[[indcount]])
          )
          indcount <- indcount + 1
        }
      }
    }
  }

  # FullsibLines <- do.call(rbind.data.frame, individualLines) #slow
  FullsibLines <- dplyr::bind_rows(individualLines)
  FullsibLines <- st_as_sf(FullsibLines, sf_column_name = "geometry")

  return(FullsibLines)
} # end pssFsLines
