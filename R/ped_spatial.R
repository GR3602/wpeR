#' Get Files For Spatial Representation Of Pedigree
#'
#' @description
#' `ped_spatial` creates georeferenced data for spatial
#' pedigree representation form the output of [`fam_table`] function.
#'
#'
#' @param famtable data frame. Output of [`fam_table`] function.
#' @param na.rm logical (`TRUE`/`FALSE`). Remove samples with missing coordinates and/or dates.
#' @param output single value or vector. Type of output of function results.
#' Available outputs: list: all spatial data returned as list, gis: all spatial data
#' returned as georeferenced files.
#' @param fullsibdata data frame. COLONY2 fullsib data.
#' @param sibthreshold numeric. p value threshold for sibsihp assignment.
#' @param path string. System path where georeferenced files should be stored.
#' @param filename string. Common name for all georeferenced files.
#' @param time.limits vector of two `Date` values. Time window for movement and
#' offspring reference samples data.
#' @param time.limit.alpha logical (`TRUE`/`FALSE`). Do time limits apply to
#' reference samples of alphas?
#' @param time.limit.offspring logical (`TRUE`/`FALSE`). Do time limits also apply
#' to reference samples of offspring?
#' @param time.limit.moves logical (`TRUE`/`FALSE`). Do time limits also apply
#' to movement data?
#'
#' @return
#' Based on the `output` parameter the function can return a list of [`sf`] objects,
#' a geospatial vector data files or both.
#'
#' Most of the objects are created separately for mothers, fathers and offspring,
#' this include: reference points (`motherRpoints`, `fatherRpoints` and
#' `offspringRpoints`), movement points (`motherMovePoints`, `fatherMovePoints`
#' and `offspringMovePoints`), movement lines (`motherMoveLines`, `fatherMoveLines`
#' and `offspringMoveLines`) and movement polygons (`motherMovePoygons`,
#' `fatherMovePolygons` and `offspringMovePolygons`).
#'
#' Besides that the function also produces lines that connect mothers and
#' their offspring (`maternityLines`), fathers and their offspring
#' (`paternityLines`) and full siblings (`FullsibLines`).
#'
#'
#' @export
#' @import sf dplyr
#'
#'
#' @examples
#' animal_ts <- anim_timespan(pack21_samples$AnimalRef,
#'                                   pack21_samples$Date,
#'                                   pack21_samples$SType,
#'                                   dead = c("Tissue", "Decomposing Tissue", "Blood"))
#'
#' sampledata <- merge(pack21_samples, animal_ts, by.x = "AnimalRef", by.y = "ID", all.x = TRUE )
#'
#' path <- paste0(system.file("extdata", package = "wpeR"), "/fake_colony")
#'
#' ped_colony <- get_colony(path, sampledata, remove_obsolete_parents = TRUE, out = "FamAgg")
#'
#' org_tables <- org_fams(ped_colony, sampledata, output = "both")
#'
#' pt<-fam_table(org_tables$fams[1,],
#'               org_tables$fams,
#'               org_tables$ped,
#'               sampledata,
#'               deadSample = c("Tissue", "Decomposing Tissue", "Blood"))
#'
#'
#' ped_spatial(pt)
#'
ped_spatial <- function (famtable,
                                    na.rm = TRUE,
                                    output="list",
                                    fullsibdata = NULL,
                                    sibthreshold = 0,
                                    path="",
                                    filename="",
                                    time.limits = c(as.Date("1900-01-01"), as.Date("2100-01-01")),
                                    time.limit.alpha = FALSE,
                                    time.limit.offspring = FALSE,
                                    time.limit.moves = FALSE) {

  #TODO add selection for type of output files gpk or shp

  data = ppsList(pedplot = famtable,
                 time.limits = time.limits,
                 na.rm = na.rm,
                 time.limit.alpha = time.limit.alpha,
                 time.limit.offspring = time.limit.offspring)

  ParLines = ppsParLines(data)

  MvPoints = ppsMvPoints(ppsData = data,
                         time.limits = time.limits,
                         time.limit.moves = time.limit.moves)

  RefPoints = ppsRefPoints(ppsData = data)

  MvLines = ppsMvLines(ppsData = data,
                       pedplot = famtable,
                       time.limits = time.limits,
                       time.limit.moves = time.limit.moves)

  MvPolygons = ppsMvPolygons(ppsData = data,
                             MvPoints = MvPoints)

  if(!is.null(fullsibdata)) {
    FsLines = ppsFsLines(ppsData = data,
                         fullsibdata = fullsibdata,
                         sibthreshold = sibthreshold)
  }


  if("gis" %in% output){
    write_sf(ParLines$maternityLines, paste0(path,filename,"matLn.gpkg"))
    write_sf(ParLines$paternityLines, paste0(path,filename,"patLn.gpkg"))

    write_sf(RefPoints$motherRpoints, paste0(path, filename, "momRef.gpkg"))
    write_sf(RefPoints$fatherRpoints, paste0(path, filename, "dadRef.gpkg"))
    write_sf(RefPoints$offspringRpoints, paste0(path, filename, "offsprRef.gpkg"))

    write_sf(MvPoints$motherMovePoints, paste0(path, filename, "momMovPt.gpkg"))
    write_sf(MvPoints$fatherMovePoints, paste0(path, filename, "dadMovPt.gpkg"))
    write_sf(MvPoints$offspringMovePoints, paste0(path, filename, "offsprMovPt.gpkg"))

    write_sf(MvLines$motherMoveLines, paste0(path, filename, "momMovLn.gpkg"))
    write_sf(MvLines$fatherMoveLines, paste0(path, filename, "dadMovLn.gpkg"))
    write_sf(MvLines$offspringMoveLines, paste0(path, filename, "offsprMovLn.gpkg"))

    write_sf(MvPolygons$motherMovePolygons, paste0(path, filename, "momMovPoly.gpkg"))
    write_sf(MvPolygons$fatherMovePolygons, paste0(path, filename, "dadMovPoly.gpkg"))
    write_sf(MvPolygons$offspringMovePolygons, paste0(path, filename, "offsprMovPoly.gpkg"))

    if(!is.null(fullsibdata)) {
      write_sf(FsLines, paste0(path, filename, "FsibLn.gpkg"))
    }
  }

  if("list" %in% output) {
    list_data = list(motherRpoints = RefPoints$motherRpoints,
                     fatherRpoints = RefPoints$fatherRpoints,
                     offspringRpoints = RefPoints$offspringRpoints,
                     motherMovePoints = MvPoints$motherMovePoints,
                     fatherMovePoints = MvPoints$fatherMovePoints,
                     offspringMovePoints = MvPoints$offspringMovePoints,
                     maternityLines = ParLines$maternityLines,
                     paternityLines = ParLines$paternityLines,
                     motherMoveLines = MvLines$motherMoveLines,
                     fatherMoveLines = MvLines$fatherMoveLines,
                     offspringMoveLines = MvLines$offspringMoveLines,
                     motherMovePolygons = MvPolygons$motherMovePolygons,
                     fatherMovePolygons = MvPolygons$fatherMovePolygons,
                     offspringMovePolygons = MvPolygons$offspringMovePolygons
    )

    if(!is.null(fullsibdata)) {
      list_data$FullsibLines = FsLines
    }

    return(list_data)
  }
}
