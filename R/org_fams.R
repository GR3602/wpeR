#' Organize animals into families and expand pedigree data
#'
#' @description
#'  Takes pedigree data from [`get_colony()`] or [`get_ped()`] function and groups animals into families.
#'  It also expands the pedigree data by adding information about the family that each individual was born in and the
#'  family in which the individual is the reproductive animal.
#'
#' @details
#'  The result of `org_fams()` function introduces us to two important concepts
#'  within the context of this package: family and polygamy cluster. A family in the
#'  output of this function is defined as a group of animals where at least one
#'  parent and at least one offspring is known. A polygamy cluster refers to a
#'  group of half-siblings, either maternally or paternally related. In the
#'  function output the `DadPclust` groups paternal half-siblings and `MomPclust`
#'  maternal half-siblings.
#'
#' @param ped Data frame. `FamAgg` output of [`get_colony()`] or [`get_ped()`] function.
#'   With `rm_obsolete_parents` parameter set to `TRUE`.
#' @param sampledata Data frame. Metadata for all genetic samples that belong
#'   to the individuals included in pedigree reconstruction analysis.
#'   This data frame should adhere to the formatting and naming conventions
#'   outlined in the [`check_sampledata()`] documentation.
#' @param output Character string. Determines the format of the output. Options are:
#'   "ped": returns an extended pedigree data frame.
#'   "fams": returns a table of all families present in the pedigree.
#'   "both": returns a list with two data frames: "ped" and "fams". (Default)
#'
#' @return
#'  Depending on the `output` parameter, the function returns either a data frame
#'  (`ped` or `fams`) or a list containing both data frames (`ped` and `fams`).
#'
#'   * `ped` data frame. An extended version of the pedigree data from `get_colony()`.
#'   In addition to common pedigree information (individual, mother, father, sex,
#'   family), `ped` includes columns for:
#'     - `parents`: Identifier codes of both parents separated with `_`.
#'     - `FamID`: Numeric identifier for the family to which the individual belongs (see `fams` below).
#'     - `FirstSeen`: Date of first sample of individual.
#'     - `LastSeen`: Date of last sample of individual.
#'     - `IsDead`: Logical value (`TRUE/FALSE`) that identifies if the individual is dead.
#'     - `DadPclust`: Identifier of father's polygamy cluster (see Details).
#'     - `MomPclust`: Identifier of mother's polygamy cuter (see Details).
#'     - `polyCluster`: Numeric value indicating if the individual is part of
#'     a polygamy cluster (see Details).
#'
#'  * `fams` data frame includes information on families that individuals in the pedigree
#'  belong to. The families are described by:
#'    - `parents`: Identifier codes of both parents separated with `_`.
#'    - `father`: Identifier code of the father.
#'    - `mother`: Identifier code of the mother.
#'    - `FamID`: Numeric identifier for the family.
#'    - `famStart`: Date when the first sample of any family member was collected.
#'    - `famEnd`: Date when the last sample of any family member was collected.
#'    - `FamDead`: Logical value (`TRUE/FALSE`) indicating if the family no longer exists.
#'    - `DadPclust`: Identifier connecting families that share the same father.
#'    - `MomPclust`: Identifier connecting families that share the same mother.
#'    - `polyCluster`: Numeric value connecting families that share one of the parents.
#'
#'
#' @export
#'
#' @examples
#'
#' # Prepare the data for usage with org_fams() function.
#' # Get animal timespan data using the anim_timespan() function.
#' animal_ts <- anim_timespan(
#'   wolf_samples$AnimalRef,
#'   wolf_samples$Date,
#'   wolf_samples$SType,
#'   dead = c("Tissue")
#' )
#' # Add animal timespan to the sampledata
#' sampledata <- merge(wolf_samples, animal_ts, by.x = "AnimalRef", by.y = "ID", all.x = TRUE)
#' # Define the path to the pedigree data file.
#' path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
#' # Retrieve the pedigree data from the get_colony function.
#' ped_colony <- get_colony(path, sampledata, rm_obsolete_parents = TRUE, out = "FamAgg")
#'
#' # Run the function
#' # Organize families and expand pedigree data using the org_fams function.
#' org_fams(
#'     ped = ped_colony,
#'     sampledata = sampledata
#'     )
#'
#' @aliases org_fams organizePacks
#'
org_fams <- function(ped, sampledata, output = "both") {
  # TODO check line 114 and 116 predispose same sample identification methods as UL -> one of sample names is animal reference
  # maybe change


  # function creates two tables [FAMS] and [PED]

  parents <- paste(ped$father, ped$mother, sep = "_")
  fams_extended <- data.frame(parents, father = ped$father, mother = ped$mother)

  # [FAMILIES]
  fams <- unique(fams_extended)
  # remove animals without both parents
  fams <- fams[!(is.na(fams$father) & is.na(fams$mother)), ]
  fams$FamID <- seq_len(nrow(fams))


  # [PED]
  # add family rep animals name to pedigree
  ped <- cbind(ped, parents)
  ## add FamID to pedigree
  ped$FamID <- fams$FamID[match(ped$parents, fams$parents)]
  ## FamID = NA for unknown animals
  ped$parents[is.na(ped$FamID)] <- NA
  ## from sample data add firs/lastseen and is dead data
  ped$FirstSeen <- sampledata$FirstSeen[match(ped$id, sampledata$Sample)]
  ped$LastSeen <- sampledata$LastSeen[match(ped$id, sampledata$Sample)]
  ped$IsDead <- sampledata$IsDead[match(ped$id, sampledata$Sample)]


  # [FAMS]
  # empty columns
  fams$FamStart <- rep(NA, nrow(fams))
  fams$FamEnd <- rep(NA, nrow(fams))
  fams$FamDead <- rep(NA, nrow(fams))


  # make family starts/family ends
  # loop adds start and end dates for families if end date is present,
  # than the function notes that the family is dead
  for (par in fams$parents) {
    fam <- fams[fams$parents == par, ]
    offspring <- ped[ped$parents == par, ]
    father <- sampledata[sampledata$Sample == fam$father, ]
    mother <- sampledata[sampledata$Sample == fam$mother, ]

    # family starts when first offspring seen
    famstart <- min(offspring$FirstSeen, na.rm = TRUE)
    # family ends when last reproductive seen
    famend <- max(c(father$LastSeen, mother$LastSeen), na.rm = TRUE)

    famdead <- FALSE
    if (!(length(mother$IsDead)) == 0) if (mother$IsDead) famdead <- TRUE
    if (!(length(father$IsDead)) == 0) if (father$IsDead) famdead <- TRUE

    fams$FamStart[fams$parents == par] <- famstart
    fams$FamEnd[fams$parents == par] <- famend
    if (!(length(famdead) == 0)) fams$FamDead[fams$parents == par] <- famdead
  }

  # dates to date format
  fams$FamStart <- as.Date(fams$FamStart, origin = "1970-01-01")
  fams$FamEnd <- as.Date(fams$FamEnd, origin = "1970-01-01")


  # [FAMS]/[PED]
  # sort out polygamy
  # table that shows which animals are fathers in more than one family -> have 2+ mates
  # unsampled fathers not included

  # DadPolyClusters = fams %>%
  #  group_by(father) %>%
  #  summarise(N = n()) %>%
  #  filter(N>1, !grepl("\\*", father))

  # Base R (stats) version
  DadPolyClusters <- stats::aggregate(fams$father, by = list(fams$father), FUN = length)
  DadPolyClusters <- DadPolyClusters[DadPolyClusters$x > 1 &
    !grepl("\\*", DadPolyClusters$Group.1), ]
  names(DadPolyClusters) <- c("father", "N")


  ## polygamy for mothers (same as above)

  # MomPolyClusters = fams %>%
  #  group_by(mother) %>%
  #  summarise(N = n()) %>%
  #  filter(N>1, !grepl("#", mother))

  # Base R (stats) version
  MomPolyClusters <- stats::aggregate(fams$mother, by = list(fams$mother), FUN = length)
  MomPolyClusters <- MomPolyClusters[MomPolyClusters$x > 1 &
    !grepl("\\*", MomPolyClusters$Group.1), ]
  names(MomPolyClusters) <- c("mother", "n")



  # numbering polygamous fathers DadP_X
  if (nrow(DadPolyClusters) > 0) {
    DadPolyClusters$PclustID <- paste("DadP", seq_len(nrow(DadPolyClusters)), sep = "_")
  } else {
    DadPolyClusters$PclustID <- character()
  }

  # numbering polygamous maother MomP_X
  if (nrow(MomPolyClusters) > 0) {
    MomPolyClusters$PclustID <- paste("MomP", seq_len(nrow(MomPolyClusters)), sep = "_")
  } else {
    MomPolyClusters$PclustID <- character()
  }

  # [FAMS]
  # write to fams
  fams$DadPclust <- DadPolyClusters$PclustID[match(fams$father, DadPolyClusters$father)]
  fams$MomPclust <- MomPolyClusters$PclustID[match(fams$mother, MomPolyClusters$mother)]

  # [PED]
  # write polygamy clusters to ped
  ped$DadPclust <- fams$DadPclust[match(ped$parents, fams$parents)]
  ped$MomPclust <- fams$MomPclust[match(ped$parents, fams$parents)]

  # [FAMS]
  # join polygamy clusters
  fams$polyCluster <- rep(NA, nrow(fams))

  counter <- 1

  for (i in seq_len(nrow(fams))) {
    if (is.na(fams$DadPclust[i]) & is.na(fams$MomPclust[i])) {
      fams$polyCluster[i] <- counter
      counter <- counter + 1
    }

    if (is.na(fams$polyCluster[i])) {
      if (!is.na(fams$DadPclust[i])) fams$polyCluster[fams$DadPclust == fams$DadPclust[i]] <- counter
      if (!is.na(fams$MomPclust[i])) fams$polyCluster[fams$MomPclust == fams$MomPclust[i]] <- counter
      counter <- counter + 1
    }
  }

  # [PED]
  # write to ped
  ped$polyCluster <- fams$polyCluster[match(ped$parents, fams$parents)]

  # sort out singleton animals (without polyclusters)
  # if they have family/polycluster NA, it's changed to 0

  ped$FamID[is.na(ped$FamID)] <- 0
  ped$polyCluster[is.na(ped$polyCluster)] <- 0

  # [FAMS]
  # make a family for singletons
  singletonFam <- fams[Inf, ]
  singletons <- ped[ped$FamID == 0, ]

  singletonFam$parents <- "Unknown"
  singletonFam$father <- "*Unknown"
  singletonFam$mother <- "#Unknown"
  singletonFam$FamID <- 0
  singletonFam$FamStart <- min(singletons$FirstSeen, na.rm = TRUE)
  singletonFam$FamEnd <- max(singletons$LastSeen, na.rm = TRUE)
  singletonFam$FamDead <- FALSE
  singletonFam$DadPclust <- NA
  singletonFam$MomPclust <- NA
  singletonFam$polyCluster <- 0

  fams <- rbind(fams, singletonFam)



  # OUTPUT
  if (output == "fams") {
    return(fams)
  }

  if (output == "ped") {
    return(ped)
  }

  if (output == "both") {
    return(list(ped = ped, fams = fams))
  }

  stop('Unknown output selected. Output should be defined as "fams", "ped" or "both"')
}
