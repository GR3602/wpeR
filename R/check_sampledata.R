#' Check and prepare sample data
#'
#' Checks the integrity and compatibility of various columns in the sample data
#' and prepares a well-structured sample data frame. It verifies the consistency of columns
#' such as Sample, Date, AnimalRef, GeneticSex, IsAnimalReference, lat, lng, and SType.
#' The function ensures that the provided data is properly formatted and conforms to the standards
#' of functions that make up the `wpeR` package.
#'
#' @param Sample A vector of sample unique identifier codes.
#' @param Date A vector of sample collection dates in 'YYYY-MM-DD' format.
#' @param AnimalRef A vector of identifier codes of the particular individual that the sample belongs to.
#' @param GeneticSex A vector of genetic sex information ('F' for female, 'M' for male, NA for unknown).
#' @param lat A vector of latitude coordinates in the WGS84 coordinate system.
#' @param lng A vector of longitude coordinates in the WGS84 coordinate system.
#' @param SType A vector of sample types.
#'
#' @return A data frame with 8 columns and a number of rows equal to the length
#' of the input vector. Each column corresponds to one of the input parameters.
#' If the function executes without warnings or errors, the result from
#' `check_sampledata()` can be used as an input parameter for other functions
#' within this package: [`get_colony()`], [`get_ped()`], [org_fams()]
#' and [`plot_table()`].
#'
#' @examples
#' sampledata <- check_sampledata(
#'   Sample = wolf_samples$Sample,
#'   Date = wolf_samples$Date,
#'   AnimalRef = wolf_samples$AnimalRef,
#'   GeneticSex = wolf_samples$GeneticSex,
#'   lat = wolf_samples$lat,
#'   lng = wolf_samples$lng,
#'   SType = wolf_samples$SType
#' )
#'
#' @export

check_sampledata <- function(Sample,
                             Date,
                             AnimalRef,
                             GeneticSex,
                             lat,
                             lng,
                             SType) {
  n <- length(Sample)
  if (n == 0) {
    stop("sample vector has length 0")
  }

  if (length(Date) != n) {
    stop(sprintf("Incompatible input: length(Sample) = %d, but length(Date) = %d",
                 n, length(Date)))
  }

  if (length(AnimalRef) != n) {
    stop(sprintf("Incompatible input: length(Sample) = %d, but length(AnimalRef) = %d",
                 n, length(AnimalRef)))
  }

  if (length(GeneticSex) != n) {
    stop(sprintf("Incompatible input: length(Sample) = %d, but length(GeneticSex) = %d",
                 n, length(GeneticSex)))
  }


    if (length(lat) != n) {
    stop(sprintf("Incompatible input: length(Sample) = %d, but length(lat) = %d",
                 n, length(lat)))
  }

  if (length(lng) != n) {
    stop(sprintf("Incompatible input: length(Sample) = %d, but length(lng) = %d",
                 n, length(lng)))
  }

  if (length(SType) != n) {
    stop(sprintf("Incompatible input: length(Sample) = %d, but length(SType) = %d",
                 n, length(SType)))
  }


  Sample <- as.character(Sample)
  Date <- try(as.Date(Date, tryFormats = "%Y-%m-%d"), silent = TRUE)
  AnimalRef <- as.character(AnimalRef)
  GeneticSex <- as.character(GeneticSex)

  lat <- tryCatch(as.numeric(lat), warning = function(cond) {
    return(NULL)
  })
  lng <- tryCatch(as.numeric(lng), warning = function(cond) {
    return(NULL)
  })
  SType <- as.character(SType)

  if (inherits(Date, "try-error")) {
    stop("Wrong date format, date shuld be formatted as 'YYYY-MM-DD'")
  }

  if (is.null(lat)) {
    stop("Looks like there are some characters in your latitude cordinate notation.")
  }

  if (is.null(lng)) {
    stop("Looks like there are some characters in your longitude cordinate notation.")
  }

  if (anyDuplicated.default(Sample) > 0) {
    stop("Duplicated entry in `Sample` vector: ", Sample[duplicated(Sample)])
  }

  if (!inherits(Date, "Date")) {
    stop("The specified `Date` vector is not in Date format")
  }

  sex <- c("F", "M", NA)
  if (!all(GeneticSex %in% sex)) {
    stop("It seems that the GeneticSex column is not coded properly.
    You must use 'M' for males, 'F' for females and NA for unknown sex.")
  }


  if (!all(-90 <= lat & 90 >= lat)) {
    stop("It seems that the latitude column stores coordinates that are not in WGS84 coordinate system.")
  }

  if (!all(-180 <= lng & 180 >= lng)) {
    stop("It seems that the latitude column stores coordinates that are not in WGS84 coordinate system.")
  }

  sampledata <- data.frame(
    Sample = Sample,
    Date = Date,
    AnimalRef = AnimalRef,
    GeneticSex = GeneticSex,
    lat = lat,
    lng = lng,
    SType = SType
  )

  return(sampledata)
}
