####INTRO####
#PREPARE PEDIGREE SPATIAL FUNCTION DISSECTION
#Function PreparePedigreeSpatial dissected into segments:
#
###ppsList - creates list of dataframes where sample data is separated into mother/father/offspring sample/reference sample
###ppsParLines - creates lines that connect mother/father and offspring reference samples
#
###ppsMvPoints - creates sf point df with mother/father/offspring move points
#
###ppsRefPoints - creates sf point df with mother/father/offspring reference sample points
#
###ppsMvLines - creates sf lines df that connects move points of single animal
####NEED TO CHECK AGAIN AND FIX SOMETHING!!!!
#
###ppsMvPolygons - creates sf polygon df that represents mcp polygon for all points of single animal
####STILL USES SP PACKAGE NEED TO CHANGE TO sf_convex_hull function of sf package!!!! -> CHANGED TO SF 20230227
#
###ppsFsLines - creates sf line df that connects samples of full siblings
#
###WRITE DESCRIPTION OF ALL FUNCTIONS!!!!!!

####FUNCTIONS:####

####ppsList####
# Prepare Pedigree Spatial data preparation
#
# Takes data from fam_table function and arranges it so that
# georeferenced data for mothers, fathers and offspring can be created
#
# @param pedplot output of PackTable function
# @param time.limits time window for movement and offspring reference samples data
# @param na.rm remove samples with missing coordinates
# @param time.limit.alpha do time limits apply to reference samples of alphas?
# @param time.limit.offspring do time limits apply to reference samples of offspring?
#

ppsList <- function (pedplot,
                     time.limits = c(as.Date("1900-01-01"), as.Date("2100-01-01")),
                     na.rm = T,
                     time.limit.alpha = F,
                     time.limit.offspring = F) {


  #remove samples with missing coordinates
  if(sum(is.na(pedplot$X)) > 0 | sum(is.na(pedplot$Y)) > 0 | sum(is.na(pedplot$Date) > 0)){
    if(na.rm == T) pedplot = pedplot[!(is.na(pedplot$X) | is.na(pedplot$Y) | is.na(pedplot$Date)),]
    else stop("Error, NA's in coordinates and/or dates. Use na.rm=TRUE to remove internally.")
  }
  #sort pedplot by individual & date
  pedplot = pedplot[order(pedplot$AnimalRef, pedplot$Date),]


  #find first&last sample within time.limits
  pedplot$first_sample = rep(FALSE, nrow(pedplot))
  pedplot$last_sample = rep(FALSE, nrow(pedplot))

  for(i in 1:nrow(pedplot)){
    curAnimal = pedplot[pedplot$AnimalRef == pedplot$AnimalRef[i] &
                          pedplot$Date <= time.limits[2] &
                          pedplot$Date >= time.limits[1],]
    curAnimal = curAnimal[order(curAnimal$Date),]

    #of no samples are within the time window, use actual first and last sample
    if (nrow(curAnimal) == 0) {
      curAnimal = pedplot[pedplot$AnimalRef == pedplot$AnimalRef[i],]
      curAnimal = curAnimal[order(curAnimal$Date),]
    }

    if (pedplot$Sample[i] == curAnimal$Sample[1]) pedplot$first_sample[i] = TRUE
    if (pedplot$Sample[i] == curAnimal$Sample[nrow(curAnimal)]) pedplot$last_sample[i] = TRUE
  }#for

  #First create pointdata
  #Use time window. Limit entire dataset, or only offspring (later)

  alphas = pedplot[pedplot$alpha == TRUE,]
  fatherAll = alphas[alphas$GeneticSex == "M",]
  motherAll = alphas[alphas$GeneticSex == "F",]

  alphaRefs = alphas[alphas$last_sample == TRUE,] #last sample of reproductive animals as spatial ref

  #use time limit window also for alphas?
  if (time.limit.alpha == T) {alphaRefs = alphaRefs[alphaRefs$Date <= time.limits[2] & alphaRefs$LastSeen >= time.limits[1], ]}

  fatherRefs = alphaRefs[alphaRefs$GeneticSex == "M",]
  motherRefs = alphaRefs[alphaRefs$GeneticSex == "F",]

  offspring = pedplot[pedplot$alpha == FALSE,]
  if (time.limit.offspring == T){
    offspringRefs = offspring[offspring$first_sample == TRUE & offspring$Date <= time.limits[2] & offspring$LastSeen >= time.limits[1], ]
  }  else { offspringRefs = offspring[offspring$first_sample == TRUE,]}

  return(list(motherAll = motherAll,
              motherRefs = motherRefs,
              fatherAll = fatherAll,
              fatherRefs = fatherRefs,
              offspring = offspring,
              offspringRefs = offspringRefs,
              alphaRefs = alphaRefs
  ))

}#end ppsList


####ppsParLines####
# Prepare Pedigree Spatial parent - offspring lines
#
# Creates sf dataframe with lines connectiong refrence samples of offspring to their parents
#
# @param ppsData output list of ppsList function
#
# @import sf

ppsParLines <- function (ppsData) {

  #require (dplyr)
  #require(sf)

  #Create Paternity / Maternity Lines (sf)
  maternityLines = NULL
  paternityLines = NULL

  dadcount = 1
  momcount = 1

  alphaRefs = ppsData$alphaRefs
  offspringRefs = ppsData$offspringRefs

  for (i in 1:nrow(offspringRefs)){
    #offspring record
    if (nrow(offspringRefs) > 0) {
      child = offspringRefs[i,]
    } else { child = NA }

    father =  alphaRefs[alphaRefs$PackID == offspringRefs$PackID[i] & alphaRefs$GeneticSex == "M",]
    mother =  alphaRefs[alphaRefs$PackID == offspringRefs$PackID[i] & alphaRefs$GeneticSex == "F",]

    if (nrow(father) > 0 & nrow(child) > 0) {
      paternityLines[[dadcount]] = rbind(father[,c("X", "Y")], child[,c("X", "Y")])
      paternityLines[[dadcount]] = st_as_sf(paternityLines[[dadcount]], coords = c("Y", "X"), crs = 4326)
      paternityLines[[dadcount]] = st_combine (paternityLines[[dadcount]])
      paternityLines[[dadcount]] = st_cast(paternityLines[[dadcount]], "LINESTRING")
      paternityLines[[dadcount]] = data.frame(ID=dadcount,
                                              pair = i,
                                              pack = offspringRefs$PackID[i],
                                              plyClust=offspringRefs$polyCluster[i],
                                              relation="paternity",
                                              child = child$AnimalRef,
                                              parent = father$AnimalRef,
                                              geometry = st_geometry(paternityLines[[dadcount]]))
      paternityLines[[dadcount]] = st_as_sf(paternityLines[[dadcount]], sf_column_name = "geometry")
      dadcount = dadcount + 1
    }

    if(nrow(mother) > 0 & nrow(child)>0){
      maternityLines[[momcount]] = rbind(mother[,c("X", "Y")], child[,c("X", "Y")])
      maternityLines[[momcount]] = st_as_sf(maternityLines[[momcount]], coords = c("Y", "X"), crs = 4326)
      maternityLines[[momcount]] = st_combine (maternityLines[[momcount]])
      maternityLines[[momcount]] = st_cast(maternityLines[[momcount]], "LINESTRING")
      maternityLines[[momcount]] = data.frame(ID=momcount,
                                              pair = i,
                                              pack = offspringRefs$PackID[i],
                                              plyClust=offspringRefs$polyCluster[i],
                                              relation="maternity",
                                              child = child$AnimalRef,
                                              parent = mother$AnimalRef,
                                              geometry = st_geometry(maternityLines[[momcount]]))
      maternityLines[[momcount]] = st_as_sf(maternityLines[[momcount]], sf_column_name = "geometry")
      momcount = momcount + 1
    }
  } #for end

  if (nrow(father) > 0 & nrow(child) > 0) {
    paternityLines = dplyr::bind_rows(paternityLines)
  }
  if (nrow(father) > 0 & nrow(child) > 0) {
    maternityLines = dplyr::bind_rows(maternityLines)
  }

  return(list(maternityLines = maternityLines,
              paternityLines = paternityLines))
}#end ppsParLines


####ppsMvPoints####
# Prepare Pedigree Spatial - Animal Movement Points
#
# @param ppsData output list of ppsList function
# @param time.limits time window for movement and offspring reference samples data
# @param time.limit.moves time limit also movement data?
#
# @import sf
#




ppsMvPoints <- function(ppsData,
                        time.limits = c(as.Date("1900-01-01"), as.Date("2100-01-01")),
                        time.limit.moves = F) {



  dupCols = c("plottingID", "PackID", "polyCluster") #columns to be removed for movelines, movepoints and refpoints to remove spatial data duplication


  ####Offspring####

  offspringRefs = ppsData$offspringRefs
  offspring = ppsData$offspring

  if (time.limit.moves == T) {
    offspringMoveData = offspring[offspring$Date <= time.limits[2] &
                                    offspring$Date >= time.limits[1] &
                                    offspring$AnimalRef %in% offspringRefs$AnimalRef,]
  } else { offspringMoveData = offspring[offspring$AnimalRef %in% offspringRefs$AnimalRef,]}

  # offspringMoveData = offspringMoveData[,!names(offspringMoveData) %in% dupCols] #drop duplication columns
  # offspringMoveData = unique(offspringMoveData)

  if(nrow(offspringMoveData) > 0) {
    offspringMovePoints = st_as_sf (offspringMoveData, coords = c("Y", "X"), crs = 4326)
  } else {offspringMovePoints = st_sf(1, st_sfc(st_point()))}

  ####Father####

  fatherRefs = ppsData$fatherRefs
  fatherAll = ppsData$fatherAll

  if (time.limit.moves == T) {
    fatherMoveData = fatherAll[fatherAll$Date <= time.limits[2] &
                                 fatherAll$Date >= time.limits[1] &
                                 fatherAll$AnimalRef %in% fatherRefs$AnimalRef,]}
  else fatherMoveData = fatherAll[fatherAll$AnimalRef %in% fatherRefs$AnimalRef,]

  fatherMoveData = fatherMoveData[,!names(fatherMoveData) %in% dupCols] #drop duplication columns
  fatherMoveData = unique(fatherMoveData)

  if(nrow(fatherMoveData) > 0) {
    fatherMovePoints  = st_as_sf (fatherMoveData, coords = c("Y", "X"), crs = 4326)
  } else {fatherMovePoints  = st_sf(1, st_sfc(st_point()))}

  ####Mother####

  motherRefs = ppsData$motherRefs
  motherAll = ppsData$motherAll

  if (time.limit.moves == T) {
    motherMoveData = motherAll[motherAll$Date <= time.limits[2] &
                                 motherAll$Date >= time.limits[1] &
                                 motherAll$AnimalRef %in% motherRefs$AnimalRef,]
  } else { motherMoveData = motherAll[motherAll$AnimalRef %in% motherRefs$AnimalRef,]}

  motherMoveData = motherMoveData[,!names(motherMoveData) %in% dupCols] #drop duplication columns
  motherMoveData = unique(motherMoveData)

  if(nrow(motherMoveData) > 0) {
    motherMovePoints <-st_as_sf (motherMoveData, coords = c("Y", "X"), crs = 4326)
  }else{motherMovePoints  = st_sf(1, st_sfc(st_point()))}

  return(list(offspringMovePoints = offspringMovePoints,
              fatherMovePoints = fatherMovePoints,
              motherMovePoints = motherMovePoints))


}#end ppsMvPoints



####ppsRefPoints####
# Prepare Pedigree Spatial - Animal Reference Points
#
# @param ppsData output list of ppsList function
#
# @import sf


ppsRefPoints <-function(ppsData) {


  dupCols = c("plottingID", "PackID", "polyCluster") #columns to be removed for movelines, movepoints and refpoints to remove spatial data duplication

  ####Offspring####

  offspringRefs = ppsData$offspringRefs

  # offspringRefs=offspringRefs[,!names(offspringRefs) %in% dupCols] #drop duplication columns
  # offspringRefs = unique(offspringRefs)

  if (nrow(offspringRefs) > 0) {
    offspringRpoints  = st_as_sf(offspringRefs, coords = c("Y", "X"), crs = 4326)
  } else {offspringRpoints = st_sf(1, st_sfc(st_point()))}

  ####Father####

  fatherRefs = ppsData$fatherRefs

  fatherRefs=fatherRefs[,!names(fatherRefs) %in% dupCols] #drop duplication columns
  fatherRefs = unique(fatherRefs)

  if(nrow(fatherRefs) > 0) {
    fatherRpoints = st_as_sf(fatherRefs, coords = c("Y", "X"), crs = 4326)
  } else {fatherRpoints = st_sf(1, st_sfc(st_point()))}

  ####Mother####

  motherRefs = ppsData$motherRefs

  motherRefs=motherRefs[,!names(motherRefs) %in% dupCols] #drop duplication columns
  motherRefs = unique(motherRefs)

  if(nrow(motherRefs) > 0) {
    motherRpoints = st_as_sf(motherRefs, coords = c("Y", "X"), crs = 4326)
  }else{motherRpoints = st_sf(1, st_sfc(st_point()))}

  return(list(offspringRpoints = offspringRpoints,
              fatherRpoints = fatherRpoints,
              motherRpoints = motherRpoints))

}#end ppsRefPoints


####ppsMvLines####
# Prepare Pedigree Spatial - Animal Movement Lines
#
# @param ppsData output list of ppsList function
# @param pedplot output of the PackTable function
# @param time.limits time window for movement and offspring reference samples data
# @param time.limit.moves time limit also movement data?
#
# @import sf
#

ppsMvLines <- function(ppsData,
                       pedplot,
                       time.limits = c(as.Date("1900-01-01"), as.Date("2100-01-01")),
                       time.limit.moves = F) {


  dupCols = c("plottingID", "PackID", "polyCluster") #columns to be removed for movelines, movepoints and refpoints to remove spatial data duplication


  ####Offspring####

  indcount=1
  individualLines = NULL
  offspringMoveLines = NULL

  offspringRefs = ppsData$offspringRefs

  for (i in 1:nrow(offspringRefs)){

    if(nrow(offspringRefs) == 0) break

    if(time.limit.moves == T){
      individual = pedplot[pedplot$AnimalRef == offspringRefs$AnimalRef[i] &
                             pedplot$Date >= time.limits[1] &
                             pedplot$Date <= time.limits[2],!names(pedplot) %in% dupCols]
    } else { individual = pedplot[pedplot$AnimalRef == offspringRefs$AnimalRef[i],!names(pedplot) %in% dupCols]}

    individual = unique(individual)

    if (nrow(individual) >= 2) {
      individualLines[[indcount]] = st_geometry(st_as_sf(individual, coords = c("Y", "X"), crs = 4326))
      individualLines[[indcount]] = st_cast(st_union(individualLines[[indcount]]), "LINESTRING")
      offspringMoveLines = rbind(offspringMoveLines,data.frame(ID=indcount,
                                                               AnimalID = offspringRefs$AnimalRef[i],
                                                               pack = offspringRefs$PackID[i],
                                                               plyClust=offspringRefs$polyCluster[i],
                                                               geometry = st_geometry(individualLines[[indcount]])))


      indcount = indcount + 1
    }
  }

  offspringMoveLines= st_as_sf(offspringMoveLines, sf_column_name = "geometry")

  #if(length(individualLines) > 0) {
  #  offspringMoveLines = SpatialLines(individualLines, CRS("+proj=longlat +datum=WGS84"))
  #  offspringMoveLines = SpatialLinesDataFrame(offspringMoveLines, data=individualDataFrame)
  #}
  #else offspringMoveLines = spLines(cbind(1,1), attr=data.frame(1))[-1,] #dummy... empty object


  ####Father####

  indcount=1
  individualLines = NULL
  fatherMoveLines = NULL

  fatherRefs = ppsData$fatherRefs

  for (i in 1:nrow(fatherRefs)){
    #offspring record
    if(nrow(fatherRefs) == 0) break

    if(time.limit.moves == T){
      individual = pedplot[pedplot$AnimalRef == fatherRefs$AnimalRef[i] &
                             pedplot$Date >= time.limits[1] &
                             pedplot$Date <= time.limits[2],!names(pedplot) %in% dupCols]
    } else { individual = pedplot[pedplot$AnimalRef == fatherRefs$AnimalRef[i],!names(pedplot) %in% dupCols]}

    individual = unique(individual)

    if (nrow(individual) >= 2) {
      individualLines[[indcount]] = st_geometry(st_as_sf(individual, coords = c("Y", "X"), crs = 4326))
      individualLines[[indcount]] = st_cast(st_union(individualLines[[indcount]]), "LINESTRING")
      fatherMoveLines = rbind(fatherMoveLines,data.frame(ID=indcount,
                                                         AnimalID = offspringRefs$AnimalRef[i],
                                                         pack = offspringRefs$PackID[i],
                                                         plyClust=offspringRefs$polyCluster[i],
                                                         geometry = st_geometry(individualLines[[indcount]])))


      indcount = indcount + 1
    }
  }

  if (nrow(individual) >= 2) {
    fatherMoveLines= st_as_sf(fatherMoveLines, sf_column_name = "geometry")
  }


  #if(length(individualLines) > 0) {
  # fatherMoveLines = SpatialLines(individualLines, CRS("+proj=longlat +datum=WGS84"))
  #  fatherMoveLines = SpatialLinesDataFrame(fatherMoveLines, data=individualDataFrame)
  #}  else {fatherMoveLines = spLines(cbind(1,1), attr=data.frame(1))[-1,]} #dummy... empty object



  ####Mother####

  indcount=1
  individualLines = NULL
  motherMoveLines = NULL

  motherRefs = ppsData$motherRefs

  for (i in 1:nrow(motherRefs)){
    #offspring record
    if(nrow(motherRefs) == 0) break

    if(time.limit.moves == T){
      individual = pedplot[pedplot$AnimalRef == motherRefs$AnimalRef[i] &
                             pedplot$Date >= time.limits[1] &
                             pedplot$Date <= time.limits[2],!names(pedplot) %in% dupCols]
    } else {individual = pedplot[pedplot$AnimalRef == motherRefs$AnimalRef[i],!names(pedplot) %in% dupCols]}

    individual = unique(individual)

    if (nrow(individual) >= 2) {
      individualLines[[indcount]] = st_geometry(st_as_sf(individual, coords = c("Y", "X"), crs = 4326))
      individualLines[[indcount]] = st_cast(st_union(individualLines[[indcount]]), "LINESTRING")
      motherMoveLines  = rbind(motherMoveLines ,data.frame(ID=indcount,
                                                           AnimalID = offspringRefs$AnimalRef[i],
                                                           pack = offspringRefs$PackID[i],
                                                           plyClust=offspringRefs$polyCluster[i],
                                                           geometry = st_geometry(individualLines[[indcount]])))


      indcount = indcount + 1
    }
  }

  motherMoveLines = st_as_sf(motherMoveLines, sf_column_name = "geometry")

  return(list(offspringMoveLines = offspringMoveLines,
              fatherMoveLines = fatherMoveLines,
              motherMoveLines = motherMoveLines))

}#end ppsMvLines


####ppsMvPolygons####
# Prepare Pedigree Spatial - Convex Hull For Animal Points
#
# @param ppsData output list of ppsList function
# @param MvPoints output of ppsMvPoints function
#
# @import sf
#
#
ppsMvPolygons <- function(ppsData, MvPoints) {

  ####Mother####
  #mother MovePolygons
  #mother MovePolygon
  indcount=1

  motherMovePoints = MvPoints$motherMovePoints
  motherRefs = ppsData$motherRefs

  motherMovePolygons = NULL
  for (i in 1:nrow(motherRefs)){
    if(nrow(motherRefs) == 0) break
    individualPoints = motherMovePoints[motherMovePoints$AnimalRef == motherRefs$AnimalRef[i],]
    if (nrow(individualPoints) >= 3) {
      animal =  unique(individualPoints$AnimalRef)
      individualPoints = individualPoints[,2]
      individualPoints = st_combine(individualPoints)
      individualPolygon = st_convex_hull(individualPoints)
      individualPolygon = st_geometry(individualPolygon)
      individualPolygon = data.frame(animal = animal,
                                     geometry = individualPolygon)
      individualPolygon = st_as_sf(individualPolygon, sf_column_name = "geometry")
      motherMovePolygons = rbind(motherMovePolygons, individualPolygon)
    } else {
      next
    }
  }


  ####Father####
  #father MovePolygons
  #father MovePolygon
  indcount=1

  fatherMovePoints = MvPoints$fatherMovePoints
  fatherRefs = ppsData$fatherRefs

  fatherMovePolygons = NULL
  for (i in 1:nrow(fatherRefs)){
    #offspring record
    if(nrow(fatherRefs) == 0) break
    individualPoints = fatherMovePoints[fatherMovePoints$AnimalRef == fatherRefs$AnimalRef[i],]
    if (nrow(individualPoints) >= 3) {
      animal =  unique(individualPoints$AnimalRef)
      individualPoints = individualPoints[,2]
      individualPoints = st_combine(individualPoints)
      individualPolygon = st_convex_hull(individualPoints)
      individualPolygon = st_geometry(individualPolygon)
      individualPolygon = data.frame(animal = animal,
                                     geometry = individualPolygon)
      individualPolygon = st_as_sf(individualPolygon, sf_column_name = "geometry")
      fatherMovePolygons = rbind(fatherMovePolygons, individualPolygon)
    } else {
      next
    }
  }





  ####Offspring####
  #offspring MovePolygons
  #offspring MovePolygon
  indcount=1

  offspringMovePoints = MvPoints$offspringMovePoints
  offspringRefs = ppsData$offspringRefs

  offspringMovePolygons = NULL
  for (i in 1:nrow(offspringRefs)){
    #offspring record
    if(nrow(offspringRefs) == 0) break
    individualPoints = offspringMovePoints[offspringMovePoints$AnimalRef == offspringRefs$AnimalRef[i],]
    if (nrow(individualPoints) >= 3) {
      animal =  unique(individualPoints$AnimalRef)
      individualPoints = individualPoints[,2]
      individualPoints = st_combine(individualPoints)
      individualPolygon = st_convex_hull(individualPoints)
      individualPolygon = st_geometry(individualPolygon)
      individualPolygon = data.frame(animal = animal,
                                     geometry = individualPolygon)
      individualPolygon = st_as_sf(individualPolygon, sf_column_name = "geometry")
      offspringMovePolygons = rbind(offspringMovePolygons, individualPolygon)
    } else {
      next
    }
  }


  return(list(motherMovePolygons = motherMovePolygons,
              fatherMovePolygons = fatherMovePolygons,
              offspringMovePolygons = offspringMovePolygons))

} #end ppsMvPolygons



####ppsFsLines####
# Prepare Pedigree Spatial - Connect Siblings
#
# @param ppsData output list of ppsList function
# @param fullsibdata COLONY2 fullsib data
# @param sibthreshold p value threshold for sibsihp assignment
#
# @import sf

ppsFsLines <- function(ppsData, fullsibdata, sibthreshold = 1) {

  #require(sf)
  #require(dplyr)

  offspringRefs = ppsData$offspringRefs

  indcount=1
  individualLines = NULL

  for (i in 1:nrow(offspringRefs)){
    #offspring record
    if(nrow(offspringRefs) == 0) break

    allsibs = fullsibdata[(fullsibdata$OffspringID1 == offspringRefs$AnimalRef[i] |
                             fullsibdata$OffspringID2 == offspringRefs$AnimalRef[i]) &
                            fullsibdata$Probability >= sibthreshold, ] #get full siblings

    #make a vector of all siblings
    sibgroup = unique(c(as.character(allsibs[,1]), as.character(allsibs[,2])))

    #samples of sibgroup
    sibgroupsamples = offspringRefs[offspringRefs$AnimalRef %in% sibgroup, ]

    #make lines between all reference samples within a siblings group
    if (nrow(sibgroupsamples) >= 2) {
      for(j in 1:(nrow(sibgroupsamples)-1)) {
        for(k in (j+1):nrow(sibgroupsamples)) {
          individualLines[[indcount]] = sibgroupsamples[c(j,k), c("Y", "X")]
          individualLines[[indcount]] = st_as_sf(individualLines[[indcount]], coords = c("Y", "X"), crs = 4326)
          individualLines[[indcount]] = st_combine(individualLines[[indcount]])
          individualLines[[indcount]] = st_cast(individualLines[[indcount]], "LINESTRING")
          individualLines[[indcount]] = data.frame(ID = indcount,
                                                   AnimalID1 = sibgroupsamples$AnimalRef[j],
                                                   AnimalID2 = sibgroupsamples$AnimalRef[k],
                                                   geometry = st_geometry(individualLines[[indcount]]))
          indcount = indcount + 1
        }
      }
    }
  }

  #FullsibLines <- do.call(rbind.data.frame, individualLines) #slow
  FullsibLines =  dplyr::bind_rows(individualLines)
  FullsibLines = st_as_sf( FullsibLines, sf_column_name = "geometry")

  return(FullsibLines)

}# end pssFsLines
