#' Organizes Pedigree Data
#'
#' @description
#' `get_ped()` offers an altretnative to `get_colony()` function in cases where the pedigree
#' was not reconstruced with [COLONY 2](https://www.zsl.org/about-zsl/resources/software/colony)
#'  software. It takes a pedigree dataframe and assigns sex to each individual.
#'  The function also prepares data so that the output of the function can be directly analysed with
#'  [`kinship2`](https://cran.r-project.org/web/packages/kinship2/index.html),
#'  [`pedtools`](https://cran.r-project.org/web/packages/pedtools/index.html) or
#'  [`FamAgg`](http://bioconductor.org/packages/release/bioc/html/FamAgg.html) packages.
#'
#' @param ped data.frame. Pedigree data frame with the most basic structure.
#'   Three columns corresponding to offspirng (has to be named `OffspringID`), father
#'   (has to be named `FatherID`), mother (has to be named `MotherID`) and family
#'   (has to be named `ClusterIndex`). Unknow parents shuld be repesented by `NA` values.
#' @param sampledata Metadata for all genetic samples that belong
#'   to the individuals included in pedigree reconstruction analysis.
#'   Must have $Sample with sample names and $GeneticSex coded as M/F/NA
#' @param out tring. For use with which package should the output be formatted?
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
#'
#' ped <- data.frame(
#'        OffspringID = c("M273P", "M20AM", "M2757", "M2ALK", "M2ETE", "M2EUJ", "MSV00E",
#'                "MSV018", "MSV05L", "MSV0M6", "MSV0T4", "MSV0T7", "MSV0TJ", "MSV0UL"),
#'        FatherID = c(NA, NA, "M20AM", "M20AM", "M20AM", "M20AM", "M20AM",
#'             "M20AM", "M20AM", "M20AM", "M20AM", "M20AM", "M20AM", "M20AM"),
#'         MotherID = c(NA, NA, "M273P", "M273P", "M273P", "M273P", "M273P",
#'             "M273P", "M273P", "M273P", "M273P", "M273P", "M273P", "M273P"),
#'        ClusterIndex = c(rep(1,14))
#'             )
#'
#' get_ped(ped, wolf_samples)
#'
#'
#'
#'
get_ped <- function (ped, sampledata, out = "FamAgg") {
  #####SEX####
  #solving funders
  ## adding sex to offspring (Offspring ID) that are also mothers/fathers -> sex retrieved form FatherID/MotherID columns
  ## if particular offsping is included in FatherID column more or equal than 1 time it gets 1 in named vector [fvect]
  fvect=sapply(ped$OffspringID, function(x)
    ifelse(sum(grepl(x,ped$FatherID,fixed=T))>0,1,0))
  ## if particular offsping is included in MotherID column more or equal than 1 time it gets 2 in named vector [mvect]
  mvect=sapply(ped$OffspringID, function(x)
    ifelse(sum(grepl(x,ped$MotherID,fixed=T))>0,2,0))

  ##adding the named vectores together result is named vetor where all offspirng have assigned sex
  ## 1 = M indcluded in FatherID
  ## 2 = F included in MotherID
  ## 0 = UNKNOWN sex cannot be determined based on maternity/paternity; animal is not a sampeld parent -> sex data has to be added from sampletipe where sex was determined by genetsics
  sex=mvect+fvect
  ##hashed in the original function
  #sex[sex==0]=3

  ##adding sex to ped
  ped$sex=sex

  ##adding sex to animals that are not parents, sex data has to be added form sampledata table
  ## vector with two levels F,M extrated from sampledata$GeneticSex for all animals where bestconfig1$sex = 0
  unknownsex=sampledata$GeneticSex[match(ped$OffspringID[ped$sex==0],sampledata$Sample)]
  unknownsex <- as.factor(unknownsex)
  ##matcihing unknowns sex levels to bestconfig levels F=2, M=1
  levels(unknownsex)=c("2","1")
  ##assignig values of unknownsex to bestconfig$sex == 0, assuming that of offspringID == sampledata$Sample, I guess that match function form few lines above takes care of that
  ped[ped$sex==0,"sex"] <- as.numeric(as.character(unknownsex))
  #sex == 3 assigned to all animals that sill have unknown sex
  ped$sex[is.na(ped$sex)] = 3

  #####IF for OUTPUTS####

  ##prepares output table so that it is compatible with data structure needef for FamAgg package
  if ( out == "FamAgg") {
    pedtable = data.frame(family=ped$ClusterIndex,
                          id = ped$OffspringID, father=ped$FatherID,
                          mother=ped$MotherID, sex=ped$sex)
    ## not sure what this does, changes datatype for sume columns to character
    for (i in 1:(ncol(pedtable)-1)) pedtable[,i] = as.character(pedtable[,i])
    return(pedtable)

  }

  ##prepares output table so that it is compatible with data structure needed for kinship2 package
  if (out == "kinship2") {
    output=data.frame(id=ped$OffspringID, dadid=ped$FatherID,
                      momid=ped$MotherID,sex=ped$sex,
                      famid = ped$ClusterIndex)
    return(output)

  }

  if (out == "pedtools") {
    output = data.frame(id = ped$OffspringID, fid = ped$FatherID,
                        mid = ped$MotherID, sex = ped$sex)
    output$fid[is.na(output$fid)] = 0
    output$mid[is.na(output$mid)] = 0
    return(output)
  }

  #should jump out with return statements... if not...
  print("Unknow output selected... aborting")
}








