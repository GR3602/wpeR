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
      paternityLines[[dadcount]] <- sf::st_as_sf(paternityLines[[dadcount]],
                                             coords = c("lng", "lat"),
                                             crs = 4326)
      paternityLines[[dadcount]] <- sf::st_combine(paternityLines[[dadcount]])
      paternityLines[[dadcount]] <- sf::st_cast(paternityLines[[dadcount]], "LINESTRING")
      paternityLines[[dadcount]] <- data.frame(
        ID = dadcount,
        pair = i,
        fam = offspringRefs$FamID[i],
        plyClust = offspringRefs$polyCluster[i],
        relation = "paternity",
        child = child$AnimalRef,
        parent = father$AnimalRef,
        geometry = sf::st_geometry(paternityLines[[dadcount]])
      )
      paternityLines[[dadcount]] <- sf::st_as_sf(paternityLines[[dadcount]],
                                             sf_column_name = "geometry")
      dadcount <- dadcount + 1
    }

    if (nrow(mother) > 0 & nrow(child) > 0) {
      maternityLines[[momcount]] <- rbind(mother[, c("lat", "lng")], child[, c("lat", "lng")])
      maternityLines[[momcount]] <- sf::st_as_sf(maternityLines[[momcount]],
                                             coords = c("lng", "lat"),
                                             crs = 4326)
      maternityLines[[momcount]] <- sf::st_combine(maternityLines[[momcount]])
      maternityLines[[momcount]] <- sf::st_cast(maternityLines[[momcount]], "LINESTRING")
      maternityLines[[momcount]] <- data.frame(
        ID = momcount,
        pair = i,
        fam = offspringRefs$FamID[i],
        plyClust = offspringRefs$polyCluster[i],
        relation = "maternity",
        child = child$AnimalRef,
        parent = mother$AnimalRef,
        geometry = sf::st_geometry(maternityLines[[momcount]])
      )
      maternityLines[[momcount]] <- sf::st_as_sf(maternityLines[[momcount]],
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
    offspringMovePoints <- sf::st_as_sf(offspringMoveData,
                                    coords = c("lng", "lat"),
                                    crs = 4326)
  } else {
    offspringMovePoints <- sf::st_sf(1, sf::st_sfc(sf::st_point()))
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
    fatherMovePoints <- sf::st_as_sf(fatherMoveData,
                                 coords = c("lng", "lat"),
                                 crs = 4326)
  } else {
    fatherMovePoints <- sf::st_sf(1, sf::st_sfc(sf::st_point()))
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
    motherMovePoints <- sf::st_as_sf(motherMoveData,
                                 coords = c("lng", "lat"),
                                 crs = 4326)
  } else {
    motherMovePoints <- sf::st_sf(1, sf::st_sfc(sf::st_point()))
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
    offspringRpoints <- sf::st_as_sf(offspringRefs,
                                 coords = c("lng", "lat"),
                                 crs = 4326)
  } else {
    offspringRpoints <- sf::st_sf(1, sf::st_sfc(sf::st_point()))
    warning("No offspring samples within selected time period.
                  Creating empty sf data frame -> offspringRpoints!")
  }

  #### Father####

  fatherRefs <- ppsData$fatherRefs

  # drop duplication columns
  fatherRefs <- fatherRefs[, !names(fatherRefs) %in% dupCols]
  fatherRefs <- unique(fatherRefs)

  if (nrow(fatherRefs) > 0) {
    fatherRpoints <- sf::st_as_sf(fatherRefs,
                              coords = c("lng", "lat"),
                              crs = 4326)
  } else {
    fatherRpoints <- sf::st_sf(1, sf::st_sfc(sf::st_point()))
    warning("No father samples within selected time period.
                  Creating empty sf data frame -> fatherRpoints!")
  }

  #### Mother####

  motherRefs <- ppsData$motherRefs

  # drop duplication columns
  motherRefs <- motherRefs[, !names(motherRefs) %in% dupCols]
  motherRefs <- unique(motherRefs)

  if (nrow(motherRefs) > 0) {
    motherRpoints <- sf::st_as_sf(motherRefs,
                              coords = c("lng", "lat"),
                              crs = 4326)
  } else {
    motherRpoints <- sf::st_sf(1, sf::st_sfc(sf::st_point()))
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

  # NEED THIS TO DEFINE VARIABLES NOT DEFINED BY FUNCTION
  # ELSE check() RETURNS NOTE no visible binding for global variable
  # solution found https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  AnimalRef <- NULL

  # columns to be removed for movelines, movepoints and refpoints to remove spatial data duplication
  dupCols <- c("plottingID", "FamID", "polyCluster")


  #### Offspring####
  if (nrow(ppsData$offspringRefs) == 0) {
    offspringMoveLines <- sf::st_sf(1, sf::st_sfc(sf::st_linestring())) # dummy... empty object
    warning("No offspring included in dataset.
          Creating empty sf data frame -> offspringMoveLines!")
  }

  offspring <- ppsData$offspring
  offspring <- unique(offspring[!names(offspring) %in% dupCols])

  offspring_sf <- sf::st_as_sf(offspring,
                               coords = c("lng", "lat"),
                               crs = 4326)

  #aggregate option just in case
  #offspring_sf <- offspring[2]
  #setNames(
    #aggregate(offspring_sf, by = list(offspring_sf$AnimalRef), FUN = length, do_union = FALSE),
    #c("animal", "no", "geometry"))

  offspringMoveLines <- offspring_sf |>
    dplyr::group_by(AnimalRef) |>
    dplyr::summarise(no_mvPoints = dplyr::n(),
                     do_union = FALSE) |>
    sf::st_cast("LINESTRING")


  if (length(which(offspringMoveLines$no_mvPoints > 1)) <= nrow(offspringMoveLines)){
    offspringMoveLines <- offspringMoveLines[which(offspringMoveLines$no_mvPoints > 1),]
  } else {
    offspringMoveLines <- sf::st_sf(1, sf::st_sfc(sf::st_linestring())) # dummy... empty object
    warning("None of the offspring has two samples needed to create a linestring.
            Creating empty sf data frame -> offspringMoveLines!")
  }

  #### Father####

  if (nrow(ppsData$fatherRefs) == 0) {
    fatherMoveLines <- sf::st_sf(1, sf::st_sfc(sf::st_linestring())) # dummy... empty object
    warning("No fathers included in dataset.
          Creating empty sf data frame -> fatherMoveLines!")
  }

  fathers <- ppsData$fatherAll
  fathers <- unique(fathers[!names(fathers) %in% dupCols])

  fathers_sf <- sf::st_as_sf(fathers,
                               coords = c("lng", "lat"),
                               crs = 4326)

  fatherMoveLines <- fathers_sf |>
    dplyr::group_by(AnimalRef) |>
    dplyr::summarise(no_mvPoints = dplyr::n(),
                     do_union = FALSE) |>
    sf::st_cast("LINESTRING")


  if (length(which(fatherMoveLines$no_mvPoints > 1)) <= nrow(fatherMoveLines)){
    fatherMoveLines <- fatherMoveLines[which(fatherMoveLines$no_mvPoints > 1),]
  } else {
    fatherMoveLines <- sf::st_sf(1, sf::st_sfc(sf::st_linestring())) # dummy... empty object
    warning("None of the fathers has two samples needed to create a linestring.
            Creating empty sf data frame -> fatherMoveLines!")
  }




  #### Mother####

  if (nrow(ppsData$motherRefs) == 0) {
    motherMoveLines <- sf::st_sf(1, sf::st_sfc(sf::st_linestring())) # dummy... empty object
    warning("No mothers included in dataset.
          Creating empty sf data frame -> motherMoveLines!")
  }

  mothers <- ppsData$motherAll
  mothers <- unique(mothers[!names(mothers) %in% dupCols])

  mothers_sf <- sf::st_as_sf(mothers,
                             coords = c("lng", "lat"),
                             crs = 4326)

  motherMoveLines <- mothers_sf |>
    dplyr::group_by(AnimalRef) |>
    dplyr::summarise(no_mvPoints = dplyr::n(),
                     do_union = FALSE) |>
    sf::st_cast("LINESTRING")


  if (length(which(motherMoveLines$no_mvPoints > 1)) <= nrow(motherMoveLines)){
    motherMoveLines <- motherMoveLines[which(motherMoveLines$no_mvPoints > 1),]
  } else {
    motherMoveLines <- sf::st_sf(1, sf::st_sfc(sf::st_linestring())) # dummy... empty object
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

  # NEED THIS TO DEFINE VARIABLES NOT DEFINED BY FUNCTION
  # ELSE check() RETURNS NOTE no visible binding for global variable
  # solution found https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  AnimalRef <- NULL

  excl_geom <- c("POINT", "LINESTRING")

  #### Mother####

  if (nrow(ppsData$motherRefs) == 0) {
    motherMovePolygons <- sf::st_sf(1, sf::st_sfc(sf::st_polygon())) # dummy... empty object
    warning("No mothers included in dataset.
          Creating empty sf data frame -> motherMovePolygons!")
  }

  motherMovePoints <- MvPoints$motherMovePoints

  motherMovePolygons <- motherMovePoints|>
    dplyr::group_by(AnimalRef) |>
    dplyr::summarise(no_mvPoints = dplyr::n()) |>
    sf::st_convex_hull()

  if(length(which(motherMovePolygons$no_mvPoints > 2)) <= nrow(motherMovePolygons)) {
    #motherMovePolygons <- motherMovePolygons[which(motherMovePolygons$no_mvPoints > 2),]
    corr_geom <- which(!(sf::st_geometry_type(motherMovePolygons) %in% excl_geom))
    motherMovePolygons <- motherMovePolygons[corr_geom,]

  } else {
    motherMovePolygons <- sf::st_sf(1, sf::st_sfc(sf::st_polygon())) # dummy... empty object
    warning("Not enough mother samples to create at least one polygon.
            Creating empty sf data frame -> motherMovePolygons!")
  }

  #### Father####

  if (nrow(ppsData$fatherRefs) == 0) {
    fatherMovePolygons <- sf::st_sf(1, sf::st_sfc(sf::st_polygon())) # dummy... empty object
    warning("No fathers included in dataset.
          Creating empty sf data frame -> fatherMovePolygons!")
  }

  fatherMovePoints <- MvPoints$fatherMovePoints

  fatherMovePolygons <- fatherMovePoints|>
    dplyr::group_by(AnimalRef) |>
    dplyr::summarise(no_mvPoints = dplyr::n()) |>
    sf::st_convex_hull()

  if(length(which(fatherMovePolygons$no_mvPoints > 2)) <= nrow(fatherMovePolygons)) {
    #fatherMovePolygons <- fatherMovePolygons[which(fatherMovePolygons$no_mvPoints > 2), ]
    corr_geom <- which(!(sf::st_geometry_type(fatherMovePolygons) %in% excl_geom))
    fatherMovePolygons <- fatherMovePolygons[corr_geom,]
  } else {
    fatherMovePolygons <- sf::st_sf(1, sf::st_sfc(sf::st_polygon())) # dummy... empty object
    warning("Not enough father samples to create at least one polygon.
            Creating empty sf data frame -> fatherMovePolygons!")
  }


  #### Offspring####

  if (nrow(ppsData$offspringRefs) == 0) {
    offspringMovePolygons <- sf::st_sf(1, sf::st_sfc(sf::st_polygon())) # dummy... empty object
    warning("No offspring included in dataset.
          Creating empty sf data frame -> offspringMovePolygons!")
  }

  offspringMovePoints <- MvPoints$offspringMovePoints

  offspringMovePolygons <- offspringMovePoints|>
    dplyr::group_by(AnimalRef) |>
    dplyr::summarise(no_mvPoints = dplyr::n()) |>
    sf::st_convex_hull()

  if(length(which(offspringMovePolygons$no_mvPoints > 2)) <= nrow(offspringMovePolygons)) {
    #offspringMovePolygons <- offspringMovePolygons[which(offspringMovePolygons$no_mvPoints > 2),]
    corr_geom <- which(!(sf::st_geometry_type(offspringMovePolygons) %in% excl_geom))
    offspringMovePolygons <- offspringMovePolygons[corr_geom,]
  } else {
    offspringMovePolygons <- sf::st_sf(1, sf::st_sfc(sf::st_polygon())) # dummy... empty object
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
          individualLines[[indcount]] <- sf::st_as_sf(individualLines[[indcount]],
                                                  coords = c("lng", "lat"),
                                                  crs = 4326)
          individualLines[[indcount]] <- sf::st_combine(individualLines[[indcount]])
          individualLines[[indcount]] <- sf::st_cast(individualLines[[indcount]],
                                                 "LINESTRING")
          individualLines[[indcount]] <- data.frame(
            ID = indcount,
            AnimalID1 = sibgroupsamples$AnimalRef[j],
            AnimalID2 = sibgroupsamples$AnimalRef[k],
            geometry = sf::st_geometry(individualLines[[indcount]])
          )
          indcount <- indcount + 1
        }
      }
    }
  }

  # FullsibLines <- do.call(rbind.data.frame, individualLines) #slow
  FullsibLines <- dplyr::bind_rows(individualLines)
  FullsibLines <- sf::st_as_sf(FullsibLines, sf_column_name = "geometry")

  return(FullsibLines)
} # end pssFsLines
