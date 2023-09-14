#' Organizes Pedigree Data
#'
#' @description
#'  Offers an alternative to [`get_colony()`] function in cases where the pedigree
#'  was not reconstructed with [COLONY 2](https://www.zsl.org/about-zsl/resources/software/colony)
#'  software. It takes a pedigree dataframe and assigns sex to each individual.
#'  The function also prepares data so that the output of the function can be directly analysed with
#'  [`kinship2`](https://cran.r-project.org/web/packages/kinship2/index.html),
#'  [`pedtools`](https://cran.r-project.org/web/packages/pedtools/index.html) or
#'  [`FamAgg`](https://bioconductor.org/packages/release/bioc/html/FamAgg.html) packages.
#'
#' @details
#'  The custom pedigree specified through the `ped` parameter should mirror the
#'  structure of a COLONY2 pedigree and share the same column names.
#'  It should consist of four columns for each offspring:
#'  `OffspringID`, `FatherID`, `MotherID` and `ClusterIndex`. In the context of
#'  COLONY2 result, the `ClusterIndex` refers to a group of offspring that may
#'  share common parents or ancestors and are analyzed together. If your
#'  your pedigree does not include such information you can fill this column with
#'  the same numeric value (eg. 1, see examples) or any (numeric) information
#'  about other family structure present in your pedigree. When considering
#'  unknown parents they should be represented by `NA` values.
#'
#'
#' @param ped Data frame. Pedigree data frame with the most basic structure.
#'   Four columns corresponding to offspring, father, mother and cluster (see
#'   Details). Unknown parents should be represented by `NA` values.
#' @param sampledata Data frame. Metadata for all genetic samples that belong
#'   to the individuals included in pedigree reconstruction analysis.
#'   This data frame should adhere to the formatting and naming conventions
#'   outlined in the [`check_sampledata()`] documentation.
#' @param out Character string. For use with which package should the output be formatted?
#'   `kinship2` (out = "kinship2"), `pedtools` (out = "pedtools") or
#'   `FamAgg` (out = "FamAgg"). Defaults to "FamAgg"
#'
#' @return
#'  A data frame describing a common pedigree structure. Each individual included in
#'  pedigree represents one row. Columns describe individual identifier code, identifier code for
#'  mother and father and sex of individual. Column names and arrangement depends on selected
#'  output (`out` parameter).
#' @export
#'
#' @examples
#' #example pedigree dataframe
#' ped <- data.frame(
#'   OffspringID = c(
#'     "M273P", "M20AM", "M2757", "M2ALK", "M2ETE", "M2EUJ", "MSV00E",
#'     "MSV018", "MSV05L", "MSV0M6", "MSV0T4", "MSV0T7", "MSV0TJ", "MSV0UL"
#'   ),
#'   FatherID = c(
#'     NA, NA, "M20AM", "M20AM", "M20AM", "M20AM", "M20AM",
#'     "M20AM", "M20AM", "M20AM", "M20AM", "M20AM", "M20AM", "M20AM"
#'   ),
#'   MotherID = c(
#'     NA, NA, "M273P", "M273P", "M273P", "M273P", "M273P",
#'     "M273P", "M273P", "M273P", "M273P", "M273P", "M273P", "M273P"
#'   ),
#'   ClusterIndex = c(rep(1, 14))
#' )
#' #Get pedigree data in FamAgg format
#' get_ped(
#'     ped = ped,
#'     sampledata = wolf_samples
#'     )
#'
get_ped <- function(ped, sampledata, out = "FamAgg") {
  ##### SEX####
  # solving founders
  # adding sex to offspring (Offspring ID) that are also mothers/fathers ->
  # sex retrieved form FatherID/MotherID columns
  # if particular offspring is included in FatherID column more or equal
  # than 1 time it gets 1 in named vector [fvect]
  fvect <- vapply(ped$OffspringID, function(x) {
    ifelse(sum(grepl(x, ped$FatherID, fixed = TRUE)) > 0, 1, 0)
  },
  FUN.VALUE = numeric(1)
  )
  # if particular offsping is included in MotherID column more or equal than 1 time it gets 2 in named vector [mvect]
  mvect <- vapply(ped$OffspringID, function(x) {
    ifelse(sum(grepl(x, ped$MotherID, fixed = TRUE)) > 0, 2, 0)
  },
  FUN.VALUE = numeric(1)
  )

  # adding the named vector together result is named vector where all offspring have assigned sex
  # 1 = M included in FatherID
  # 2 = F included in MotherID
  # 0 = UNKNOWN sex cannot be determined based on maternity/paternity
  sex <- mvect + fvect

  # adding sex to ped
  ped$sex <- sex

  # adding sex to animals that are not parents, sex data has to be added form sampledata table
  # vector with two levels F,M extracted from sampledata$GeneticSex for all animals where bestconfig1$sex = 0
  unknownsex <- sampledata$GeneticSex[match(ped$OffspringID[ped$sex == 0], sampledata$Sample)]
  unknownsex <- as.factor(unknownsex)
  # matcihing unknowns sex levels to bestconfig levels F=2, M=1
  levels(unknownsex) <- c("2", "1")
  # assignig values of unknownsex to bestconfig$sex == 0
  ped[ped$sex == 0, "sex"] <- as.numeric(as.character(unknownsex))
  # sex == 3 assigned to all animals that sill have unknown sex
  ped$sex[is.na(ped$sex)] <- 3

  ##### IF for OUTPUTS####

  if (out == "FamAgg") {
    pedtable <- data.frame(
      family = ped$ClusterIndex,
      id = ped$OffspringID, father = ped$FatherID,
      mother = ped$MotherID, sex = ped$sex
    )

    for (i in 1:(ncol(pedtable) - 1)) pedtable[, i] <- as.character(pedtable[,i])
    return(pedtable)
  }


  if (out == "kinship2") {
    output <- data.frame(
      id = ped$OffspringID, dadid = ped$FatherID,
      momid = ped$MotherID, sex = ped$sex,
      famid = ped$ClusterIndex
    )
    return(output)
  }

  if (out == "pedtools") {
    output <- data.frame(
      id = ped$OffspringID, fid = ped$FatherID,
      mid = ped$MotherID, sex = ped$sex
    )
    output$fid[is.na(output$fid)] <- 0
    output$mid[is.na(output$mid)] <- 0
    return(output)
  }

  # should jump out with return statements... if not...
  stop("Unknown output selected. The function excepts
    'kinship2', 'pedtools', 'FamAgg' and 'table' as out parameter strings")
}
