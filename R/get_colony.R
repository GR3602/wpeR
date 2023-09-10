#' Organizes COLONY2 output
#'
#' @description
#' `get_colony()` extends `BestConfig_Ordered` output from [COLONY2](https://www.zsl.org/about-zsl/resources/software/colony)
#'  pedigree reconstruction software with additional data about individuals included
#'  in pedigree. The function adds missing parents to `OffspringID`, assigns
#'  sex to each individual included in `OffspringID` and adds the computed
#'  probabilities of paternity and maternity assignments (probability of assignments is visible
#'  only if the `out` parameter is set to `"table"`).
#'  The function also prepares data so that the output of the function can be directly analysed with
#'  [`kinship2`](https://cran.r-project.org/web/packages/kinship2/index.html),
#'  [`pedtools`](https://cran.r-project.org/web/packages/pedtools/index.html) or
#'  [`FamAgg`](http://bioconductor.org/packages/release/bioc/html/FamAgg.html) packages.
#'
#' @details
#'  COLONY2 output tables needed for this function (`.BestConfig_Ordered`, `.Maternity`
#'  and `.Paternity`)
#'  are read directly from the colony output folder and do not need to be imported into R session. The path
#'  to the outputs is defined with `bestconf_path` parameter. When defining `bestconf_path` the user needs to
#'  define a complete path to the directory where colony outputs are stored and also the file name
#'  (file name of COLONY2 outputs equals the project name \cr
#'  eg. /path/to/the/COLONY2/output/folder/COLONY2_project_name).
#'
#'
#'
#' @param bestconf_path Character string. Path to the folder where COLONY2 output files are saved.
#'   Has to include file path and project name (see Details).
#' @param sampledata Data frame. Metadata for all genetic samples that belong
#'   to the individuals included in pedigree reconstruction analysis.
#'   Must have `$Sample` with sample names and `$GeneticSex` coded as `M/F/NA`
#' @param rm_obsolete_parents Logical. Should unknown parents be removed from output.
#'   Applies just to offspring for which both parents are unknown. Defaults to `TRUE`.
#' @param out Character string. For use with which package should the output be formatted?
#'   `kinship2` (out = "kinship2"), `pedtools` (out = "pedtools"),
#'   `FamAgg` (out = "FamAgg") or the created data.frame can be outputted as is
#'   (out = "table"). Defaults to "FamAgg"
#'
#' @return
#'  A data frame describing a common pedigree structure. Each individual included in
#'  pedigree represents one row. Columns describe individual identifier code, identifier code for
#'  mother and father, sex and family of individual. Column names and arrangement depends on selected
#'  output (`out` parameter).
#'
#'
#' @export
#'
#' @examples
#' # Define the path to COLONY2 output
#' path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
#'
#' # Get pedigree data in FamAgg format
#' get_colony(path, wolf_samples)
#'
#' @aliases get_colony GetDigestColony
#'


get_colony <- function(bestconf_path,
                       sampledata,
                       rm_obsolete_parents = TRUE,
                       out = "FamAgg") {

  # reads COLONY outputs needed for the function
  bestconfig <- utils::read.table(
    paste(bestconf_path, ".BestConfig_Ordered", sep = ""),
    sep = "",
    header = TRUE,
    fill = TRUE,
    comment.char = "?"
  )

  paternity <- utils::read.table(
    paste(bestconf_path, ".Paternity", sep = ""),
    sep = ",",
    header = TRUE,
    fill = TRUE,
    comment.char = "?"
  )

  maternity <- utils::read.table(
    paste(bestconf_path, ".Maternity", sep = ""),
    sep = ",",
    header = TRUE,
    fill = TRUE,
    comment.char = "?"
  )

  ##### UNKNOWN PARENTS#####
  # IDs of fathers with no samples -> vector with *XX codes (COLONY generated unknown fathers)
  ufathers <- unique(bestconfig$FatherID[!bestconfig$FatherID %in% bestconfig$OffspringID])
  # IDs of mothers with no samples -> vector with #XX codes (COLONY generated unknown mothers)
  umothers <- unique(bestconfig$MotherID[!bestconfig$MotherID %in% bestconfig$OffspringID])

  # creates DF of unknown fathers, compatible with BestConfig_ordered data
  ufathers <- data.frame(
    ufathers, as.factor(rep(0, length(ufathers))), as.factor(rep(0, length(ufathers))),
    bestconfig$ClusterIndex[match(ufathers, bestconfig$FatherID)]
  )
  # creates DF of unknown mothers, compatible with BestConfig_ordered data
  umothers <- data.frame(
    umothers, as.factor(rep(0, length(umothers))), as.factor(rep(0, length(umothers))),
    bestconfig$ClusterIndex[match(umothers, bestconfig$MotherID)]
  )

  # assign BestConfig_ordered names to DF created above
  names(ufathers) <- names(bestconfig)
  names(umothers) <- names(bestconfig)

  # adding unknown parents data to bestconfig
  bestconfig1 <- rbind(bestconfig, ufathers, umothers)

  # changes mother and father ID for unknown animals form 0 to NA
  bestconfig1$FatherID[bestconfig1$FatherID == 0] <- NA
  bestconfig1$MotherID[bestconfig1$MotherID == 0] <- NA


  ##### PROBABILITY####
  # adding probability form paternity/maternity colony outputs
  names(paternity)[3] <- "prob"
  names(maternity)[3] <- "prob"

  paternity_p <- paternity[match(bestconfig1$OffspringID, paternity$OffspringID), "prob"]
  maternity_p <- maternity[match(bestconfig1$OffspringID, maternity$OffspringID), "prob"]


  bestconfig1 <- cbind(bestconfig1, paternity_p, maternity_p)


  ##### SEX####
  # solving founders
  # adding sex to offspring (Offspring ID) that are also mothers/fathers -> sex retrieved form FatherID/MotherID columns
  # if particular offspring is included in FatherID column more or equal than 1 time it gets 1 in named vector [fvect]
  fvect <- vapply(bestconfig1$OffspringID, function(x) {
    ifelse(sum(grepl(x, bestconfig1$FatherID, fixed = TRUE)) > 0, 1, 0)
  },
  FUN.VALUE = numeric(1)
  )
  # if particular offspring is included in MotherID column more or equal than 1 time it gets 2 in named vector [mvect]
  mvect <- vapply(bestconfig1$OffspringID, function(x) {
    ifelse(sum(grepl(x, bestconfig1$MotherID, fixed = TRUE)) > 0, 2, 0)
  },
  FUN.VALUE = numeric(1)
  )

  # adding the named vector together result is named vector where all offspring have assigned sex
  # 1 = M included in FatherID
  # 2 = F included in MotherID
  # 0 = UNKNOWN sex cannot be determined based on maternity/paternity; a
  # animal is not a sampled parent -> sex data has to be added from sampler where sex was determined by genetics
  sex <- mvect + fvect

  # adding sex to bestconfig
  bestconfig1$sex <- sex

  # adding sex to animals that are not parents, sex data has to be added form sampledata table
  # vector with two levels F,M extracted from sampledata$GeneticSex for all animals where bestconfig1$sex = 0
  unknownsex <- sampledata$GeneticSex[match(bestconfig1$OffspringID[bestconfig1$sex == 0], sampledata$Sample)]
  unknownsex <- as.factor(unknownsex)
  ## matching unknowns sex levels to bestconfig levels F=2, M=1
  levels(unknownsex) <- c("2", "1")
  # assigning values of unknownsex to bestconfig$sex == 0,
  bestconfig1[bestconfig1$sex == 0, "sex"] <- as.numeric(as.character(unknownsex))
  # sex == 3 assigned to all animals that sill have unknown sex
  bestconfig1$sex[is.na(bestconfig1$sex)] <- 3

  ##### IF for OUTPUTS####
  # if you want to remove unknown parents form output of the function
  if (rm_obsolete_parents == T) {
    # removes #XX and *XX values from FatherID and MotherID column,
    # just for offspring where both parents are unknown
    bestconfig1[
      (grepl("\\*", bestconfig1$FatherID) & grepl("#", bestconfig1$MotherID)),
      c("FatherID", "MotherID")
    ] <- c(NA, NA)

    # parents that were removed from FadtherID and MotherID column are also removed form offspingID column
    bestconfig1$OffspringID[(grepl("\\*", bestconfig1$OffspringID) | grepl("#", bestconfig1$OffspringID)) &
      (!(bestconfig1$OffspringID %in% bestconfig1$FatherID) &
        !(bestconfig1$OffspringID %in% bestconfig1$MotherID))] <- NA

    # from output table removes row that have NA in OffspingID column
    bestconfig1 <- bestconfig1[!is.na(bestconfig1$OffspringID), ]
  }


  if (out == "FamAgg") {
    pedtable <- data.frame(
      family = bestconfig1$ClusterIndex,
      id = bestconfig1$OffspringID, father = bestconfig1$FatherID,
      mother = bestconfig1$MotherID, sex = bestconfig1$sex
    )

    for (i in 1:(ncol(pedtable) - 1)) pedtable[, i] <- as.character(pedtable[, i])
    return(pedtable)
  }


  if (out == "kinship2") {
    output <- data.frame(
      id = bestconfig1$OffspringID, dadid = bestconfig1$FatherID,
      momid = bestconfig1$MotherID, sex = bestconfig1$sex,
      famid = bestconfig1$ClusterIndex
    )
    return(output)
  }

  if (out == "pedtools") {
    output <- data.frame(
      id = bestconfig1$OffspringID, fid = bestconfig1$FatherID,
      mid = bestconfig1$MotherID, sex = bestconfig1$sex
    )
    output$fid[is.na(output$fid)] <- 0
    output$mid[is.na(output$mid)] <- 0
    return(output)
  }

  if (out == "table") {
    bestconfig1$sex[bestconfig1$sex == 1] <- "M"
    bestconfig1$sex[bestconfig1$sex == 2] <- "F"
    output <- bestconfig1
    return(output)
  }

  # should jump out with return statements... if not...
  stop("Unknown output selected. The function excepts
    'kinship2', 'pedtools', 'FamAgg' and 'table' as out parameter strings")
}
