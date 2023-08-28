#' Organizes animals into families and expands pedigree data
#'
#' @description
#' `org_fams` takes pedigree data from [`get_colony`] function and groups animals into families.
#'  It also expands the pedigree data by adding information about the family that each individual was born in and the
#'  in which the individual is the reproductive animal. A family in this function is defined as a group of animals
#'  where at least one parent and at least one offspring is known.
#'
#' @param ped data frame. `FamAgg` output of [`get_colony`] function.
#'   With `RemoveObsoleteParents` parameter set to `TRUE`.
#' @param sampledata data frame. Metadata for all genetic samples that belong
#'   to the individuals included in pedigree reconstruction analysis.
#'   Must have $Sample with sample names and $GeneticSex coded as M/F/NA
#' @param output string. How the output data frame should be formatted.
#'   output = "ped" returns pedigree data (similar to [`org_fams`] output)
#'   with additional information on individuals. output = "fams" returns a table of all
#'   families present in the pedigree. output = "both" returns a list with two
#'   data frames: fams and ped. Defaults to "both"
#'
#' @return
#' Based on the ´output´ parameter the function can return a data frame (`ped` or `fams`)
#' or a list with two objects (`ped` and `fams`).
#'   * `ped` data frame. Extended output of [`get_colony`] function.
#'   Apart from common pedigree information (individual, mother, father, sex, family), `ped`
#'   also includes columns:
#'     - `parents`: identifier codes of both parents separated with `_`,
#'     - `FamID`: number of family that the individual belongs to (see `fams` below),
#'     - `FirstSeen`: date of first sample of individual,
#'     - `LastSeen`: date of last sample of individual,
#'     - `IsDead`: logical value (`TRUE/FALSE`) that identifies if the individual is dead,
#'     - `DadPclust`: identifier of fathers polygamy cluster,
#'     - `MomPclust`: identifier of mothers polygamy cuter,
#'     - `polyCluster`: polygamy cluster of the individual.
#'
#'  * `fams` data frame includes information on families that individuals in the pedigree
#'  belong to. The families are described by:
#'    - `parents`: identifier codes of both parents separated with `_`,
#'    - `father`: identifier code of the father,
#'    - `mother`: identifier code of the mother,
#'    - `FamID`: numeric value that identifies a particular family,
#'    - `famStart`: date when the first sample of any of the family members was collected,
#'    - `famEnd`: date when the last sample of any of the family members was collected,
#'    - `FamDead`: logical value (`TRUE/FALSE`) that identifies if the family does not exist any more,
#'    - `DadPclust`: identifier that connects families that share the same father,
#'    - `MomPclust`: identifier that connects families that share the same mother,
#'    - `polyCluster`: numeric value that connects families that share one of the parents.
#'
#'
#' @export
#'
#' @examples
#' animal_ts <- anim_timespan(wolf_samples$AnimalRef,
#'                                   wolf_samples$Date,
#'                                   wolf_samples$SType,
#'                                   dead = c("Tissue", "Decomposing Tissue", "Blood"))
#'
#' sampledata <- merge(wolf_samples, animal_ts, by.x = "AnimalRef", by.y = "ID", all.x = TRUE )
#'
#' path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
#'
#' ped_colony <- get_colony(path, sampledata, remove_obsolete_parents = TRUE, out = "FamAgg")
#'
#' org_fams(ped_colony, sampledata)
#'
#' @aliases org_fams organizePacks
#'
org_fams = function(ped, sampledata, output = "both") {

  #TODO check line 114 and 116 predispose same sample identification methoda as UL -> one of sample names is animal reference
  #maybe change


  #function creates two tables [FAMS] and [PED]

  parents = paste(ped$father, ped$mother, sep ="_")
  fams_extended = data.frame(parents, father=ped$father, mother=ped$mother)

  #[FAMILIES]
  fams = unique(fams_extended)
  fams = fams[!(is.na(fams$father) & is.na(fams$mother)),] #remove animals without both parents
  fams$FamID = 1:nrow(fams)


  #[PED]
  #write back into pedigree
  ##add family rep animals name to pedigree
  ped = cbind(ped,parents)
  ##add FamID to pedigree
  ped$FamID = fams$FamID[match(ped$parents, fams$parents)]
  ##FamID = NA for unknown animals
  ped$parents[is.na(ped$FamID)] = NA
  ##from sample data add firs/lastseena and is dead data
  ped$FirstSeen = sampledata$FirstSeen[match(ped$id, sampledata$Sample)]
  ped$LastSeen = sampledata$LastSeen[match(ped$id, sampledata$Sample)]
  ped$IsDead = sampledata$IsDead[match(ped$id, sampledata$Sample)]


  #[FAMS]
  ##empty columns
  fams$FamStart = rep(NA, nrow(fams))
  fams$FamEnd = rep(NA, nrow(fams))
  fams$FamDead = rep(NA, nrow(fams))


  #make family starts/family ends
  ##loop adds start and end dates for families if end date is present, than the function notes that the family is dead
  for(par in fams$parents) {

    fam = fams[fams$parents == par,]
    offspring = ped[ped$parents == par,]
    father = sampledata[sampledata$Sample == fam$father,]
    mother = sampledata[sampledata$Sample == fam$mother,]


    famstart = min(offspring$FirstSeen, na.rm=T) #family starts when first offspring seen
    famend = max(c(father$LastSeen, mother$LastSeen), na.rm=T) #family ends when last alpha seen

    famdead = FALSE
    if (!(length(mother$IsDead)) ==0) if (mother$IsDead) famdead = TRUE
    if (!(length(father$IsDead)) ==0) if (father$IsDead) famdead = TRUE

    fams$FamStart[fams$parents == par] = famstart
    fams$FamEnd[fams$parents == par] = famend
    if(!(length(famdead) == 0)) fams$FamDead[fams$parents == par] = famdead


  }

  ##dates to date format
  fams$FamStart = as.Date(fams$FamStart, origin = "1970-01-01")
  fams$FamEnd = as.Date(fams$FamEnd, origin = "1970-01-01")


  #[FAMS]/[PED]
  #sort out polygamy
  ##table that shows whic andimals are fathers in more than one family -> have 2+ mates
  ##unsampled fathers not included

  #DadPolyClusters = fams %>%
  #  group_by(father) %>%
  #  summarise(N = n()) %>%
  #  filter(N>1, !grepl("\\*", father))

  #Base R (stats) version
  DadPolyClusters = stats::aggregate(fams$father, by = list(fams$father),  FUN = length)
  DadPolyClusters = DadPolyClusters[DadPolyClusters$x > 1 &
                                      !grepl("\\*", DadPolyClusters$Group.1) ,]
  names(DadPolyClusters) = c("father", "N")


  ##polygamy for mothers (same as above)

  #MomPolyClusters = fams %>%
  #  group_by(mother) %>%
  #  summarise(N = n()) %>%
  #  filter(N>1, !grepl("#", mother))

  #Base R (stats) version
  MomPolyClusters = stats::aggregate(fams$mother, by = list(fams$mother),  FUN = length)
  MomPolyClusters = MomPolyClusters[MomPolyClusters$x > 1 &
                                      !grepl("\\*", MomPolyClusters$Group.1) ,]
  names(MomPolyClusters) = c("mother", "n")



  ##numbering polygamous fathers DadP_X
  if(nrow(DadPolyClusters) > 0)
    DadPolyClusters$PclustID = paste("DadP",1:nrow(DadPolyClusters),sep="_")
  else
    DadPolyClusters$PclustID = character()

  #numbering polygamous maother MomP_X
  if(nrow(MomPolyClusters) > 0)
    MomPolyClusters$PclustID = paste("MomP",1:nrow(MomPolyClusters),sep="_")
  else
    MomPolyClusters$PclustID = character()

  #[FAMS]
  #write to fams
  fams$DadPclust = DadPolyClusters$PclustID[match(fams$father, DadPolyClusters$father)]
  fams$MomPclust = MomPolyClusters$PclustID[match(fams$mother, MomPolyClusters$mother)]

  #[PED]
  #write polygamy clusters to ped
  ped$DadPclust = fams$DadPclust[match(ped$parents, fams$parents)]
  ped$MomPclust = fams$MomPclust[match(ped$parents, fams$parents)]

  #[FAMS]
  #join polygamy clusters
  fams$polyCluster = rep(NA, nrow(fams))

  counter = 1

  for (i in 1:nrow(fams)){

    if(is.na(fams$DadPclust[i]) & is.na(fams$MomPclust[i])) {
      fams$polyCluster[i] = counter
      counter = counter + 1
    }

    if(is.na(fams$polyCluster[i])) {
      if(!is.na(fams$DadPclust[i])) fams$polyCluster[fams$DadPclust == fams$DadPclust[i]] = counter
      if(!is.na(fams$MomPclust[i])) fams$polyCluster[fams$MomPclust == fams$MomPclust[i]] = counter
      counter = counter + 1
    }

  }#for

  #[PED]
  #write to ped
  ped$polyCluster = fams$polyCluster[match(ped$parents, fams$parents)]

  #sort out singleton animals (without polyclusters)
  #if they have family/polycluster NA, it's changed to 0

  ped$FamID[is.na(ped$FamID)] = 0
  ped$polyCluster[is.na(ped$polyCluster)] = 0

  #[FAMS]
  #make a family for singletons
  singletonFam = fams[Inf,]
  singletons = ped[ped$FamID == 0,]

  singletonFam$parents = "Unknown"
  singletonFam$father = "*Unknown"
  singletonFam$mother = "#Unknown"
  singletonFam$FamID = 0
  singletonFam$FamStart = min(singletons$FirstSeen, na.rm=T)
  singletonFam$FamEnd = max(singletons$LastSeen, na.rm=T)
  singletonFam$FamDead = FALSE
  singletonFam$DadPclust = NA
  singletonFam$MomPclust = NA
  singletonFam$polyCluster = 0

  fams = rbind(fams, singletonFam)



  # OUTPUT
  if (output == "fams"){
    return(fams)
  }

  if (output == "ped"){
    return(ped)
  }

  if (output == "both"){
    return(list(ped=ped, fams=fams))
  }

  stop('Unknown output seleced. Output shuld be defined as "fams", "ped" or "both"')

}
