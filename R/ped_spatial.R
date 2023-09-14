#' Get files for spatial representation of pedigree
#'
#' @description
#' Creates georeferenced data for spatial
#' pedigree representation form the output of [`plot_table()`] function.
#'
#' @details
#' The parameters `path`, `filename` and `out.format`, are used only when `output`
#' parameter is set to "gis", since they control which georeferenced files should
#' be created, where they will be saved and which common file name will they have.
#'
#'
#'
#' @param plottable Data frame. Output of [`plot_table()`] function.
#' @param na.rm Logical (`TRUE`/`FALSE`). Remove samples with missing coordinates and/or dates.
#' @param output Character vector specifying the desired output type ('list' - default or 'gis').
#' Available outputs: list: all spatial data returned as list, gis: all spatial data
#' returned as georeferenced files.
#' @param fullsibdata Data frame with COLONY2 full-sibling data.
#' @param sibthreshold Numeric. P-value threshold for sibship assignment.
#' @param path System path for storing georeferenced files.
#' @param filename Common name for all georeferenced files.
#' @param out.format Character string. Type of georeferenced files to be generated.
#' Can be ether `"geopackage"` or `"shapefile"`. Default is `"geopackage"`
#' @param time.limits Vector of two `Date` values as the time window.
#' @param time.limit.rep Logical (`TRUE`/`FALSE`). Apply time limits to
#' reference samples of reproductive animals.
#' @param time.limit.offspring Logical (`TRUE`/`FALSE`). Apply time limits to
#' reference samples of offspring.
#' @param time.limit.moves Logical (`TRUE`/`FALSE`). Apply time limits to
#' movement data.
#'
#' @return
#' Depending on the `output` parameter the function can return a list of [`sf`] objects,
#' a georeferenced vector data files or both.
#'
#' Most of the objects are created separately for mothers, fathers and offspring,
#' this include:
#'
#' - Reference Points (`motherRpoints`, `fatherRpoints`, and `offspringRpoints`).
#'    - Each point corresponds to an animal included in the 'plot_table()'
#'    function output.
#'    - For reproductive animals (mothers and fathers), a reference point is the
#'    location of their last sample within the specified time window.
#'    - For offspring, the reference point is the location of their first sample
#'    within the time window.
#'
#' - Movement Points (`motherMovePoints`, `fatherMovePoints`, and `offspringMovePoints`).
#'    - These points represent all the samples of the respective animals.
#'
#' - Movement Lines (`motherMoveLines`, `fatherMoveLines` and `offspringMoveLines`).
#'    - Movement lines connect all '...MovePoints' of a specific animal in
#'    chronological order.
#'
#' - Movement Polygons (`motherMovePolygons`, `fatherMovePolygons` and `offspringMovePolygons`):
#'    - Movement polygons represent a convex hull that encloses all the samples of an individual.
#'    -  An individual must have more than two samples for this representation.
#'
#' Besides that the function also produces lines that connect mothers and
#' their offspring (`maternityLines`), fathers and their offspring
#' (`paternityLines`), and if `fullsibdata` parameter is specified,
#'  full siblings (`FullsibLines`).
#'
#'
#' @export
#' @import sf dplyr
#'
#'
#' @examples
#' # Prepare the data for usage with ped_spatial() function.
#' # Get animal timespan data using the anim_timespan() function.
#' animal_ts <- anim_timespan(wolf_samples$AnimalRef,
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
#' # Organize families and expand pedigree data using the org_fams function.
#' org_tables <- org_fams(ped_colony, sampledata, output = "both")
#' # Prepare data for plotting.#'
#' pt <- plot_table(org_tables$fams[1, ],
#'   org_tables$fams,
#'   org_tables$ped,
#'   sampledata,
#'   deadSample = c("Tissue", "Decomposing Tissue", "Blood")
#' )
#'
#' # Run the function
#' # Get files for spatial pedigree representation in list format.
#' ped_spatial(plottable = pt)
#'
#' @aliases ped_spatial PreparePedigreeSpatial
#'
ped_spatial <- function(plottable,
                        na.rm = TRUE,
                        output = "list",
                        fullsibdata = NULL,
                        sibthreshold = 0,
                        path = "",
                        filename = "",
                        out.format = "geopackage",
                        time.limits = c(as.Date("1900-01-01"), as.Date("2100-01-01")),
                        time.limit.rep = FALSE,
                        time.limit.offspring = FALSE,
                        time.limit.moves = FALSE) {
  data <- ppsList(
    plottable = plottable,
    time.limits = time.limits,
    na.rm = na.rm,
    time.limit.rep = time.limit.rep,
    time.limit.offspring = time.limit.offspring,
    time.limit.moves = time.limit.moves
  )

  ParLines <- ppsParLines(data)

  MvPoints <- ppsMvPoints(ppsData = data)

  RefPoints <- ppsRefPoints(ppsData = data)

  MvLines <- ppsMvLines(ppsData = data)

  MvPolygons <- ppsMvPolygons(
    ppsData = data,
    MvPoints = MvPoints
  )

  if (!is.null(fullsibdata)) {
    FsLines <- ppsFsLines(
      ppsData = data,
      fullsibdata = fullsibdata,
      sibthreshold = sibthreshold
    )
  }


  if ("gis" %in% output) {
    if ("geopackage" %in% out.format) {
      ext <- ".gpkg"
    } else if ("shapefile" %in% out.format) {
      ext <- ".shp"
    } else {
      stop("Wrong out.format parameter. The out.format can either be 'geopackage' or 'shapefile'")
    }

    write_sf(ParLines$maternityLines, paste0(path, filename, "matLn", ext))
    write_sf(ParLines$paternityLines, paste0(path, filename, "patLn", ext))

    write_sf(RefPoints$motherRpoints, paste0(path, filename, "momRef", ext))
    write_sf(RefPoints$fatherRpoints, paste0(path, filename, "dadRef", ext))
    write_sf(RefPoints$offspringRpoints, paste0(path, filename, "offsprRef", ext))

    write_sf(MvPoints$motherMovePoints, paste0(path, filename, "momMovPt", ext))
    write_sf(MvPoints$fatherMovePoints, paste0(path, filename, "dadMovPt", ext))
    write_sf(MvPoints$offspringMovePoints, paste0(path, filename, "offsprMovPt", ext))

    write_sf(MvLines$motherMoveLines, paste0(path, filename, "momMovLn", ext))
    write_sf(MvLines$fatherMoveLines, paste0(path, filename, "dadMovLn", ext))
    write_sf(MvLines$offspringMoveLines, paste0(path, filename, "offsprMovLn", ext))

    write_sf(MvPolygons$motherMovePolygons, paste0(path, filename, "momMovPoly", ext))
    write_sf(MvPolygons$fatherMovePolygons, paste0(path, filename, "dadMovPoly", ext))
    write_sf(MvPolygons$offspringMovePolygons, paste0(path, filename, "offsprMovPoly", ext))

    if (!is.null(fullsibdata)) {
      write_sf(FsLines, paste0(path, filename, "FsibLn", ext))
    }
  }

  if ("list" %in% output) {
    list_data <- list(
      motherRpoints = RefPoints$motherRpoints,
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

    if (!is.null(fullsibdata)) {
      list_data$FullsibLines <- FsLines
    }

    return(list_data)
  }
}
