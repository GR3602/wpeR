#' Check and prepare genetic sample metadata
#'
#' Verifies the consistency of columns in the genetic sample metadata and
#' prepares it for use with other functions in the `wpeR` package. The function
#' ensures that the provided data is properly formatted and conforms to the standards
#' of functions that make up the `wpeR` package.
#'
#' @details
#' By specifying the `extraCols` parameter additional information can be included
#' in the sampledata dataframe. Such additional information is not required for
#' the functioning of the `wpeR` package functions, but can be useful to the user
#' when interpreting results. When including additional columns the function
#' inputs (Sample, Date, AnimalRef...) have to be defined as a vector extracted from
#' data frame column (eg. Sample = dataframe$column) and the `extraCols` parameter is defined,
#' as a vector of column names form the same data frame (eg. extraCols = c(column1, column2,
#' column3)).
#'
#'
#' @param Sample A vector of sample unique identifier codes.
#' @param Date A vector of sample collection dates in 'YYYY-MM-DD' format.
#' @param AnimalRef A vector of identifier codes of the particular individual
#' that the sample belongs to.
#' @param GeneticSex A vector of genetic sex information
#' ('F' for female, 'M' for male, NA for unknown).
#' @param lat A vector of latitude coordinates in the WGS84 coordinate system
#' (EPSG: 4326).
#' @param lng A vector of longitude coordinates in the WGS84 coordinate system
#' (EPSG: 4326).
#' @param SType A vector of sample types eg.: scat, hair, tissue.
#' @param extraCols A vector of extra column names that the user wants to include
#' in sampledata data frame (see Details).
#'
#' @return A data frame with 7 columns and a number of rows equal to the length
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
                             SType,
                             extraCols = NULL) {
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

  #has to be sample with small s otherwise sub("\\$.*", "", deparse(substitute(sample)))
  #returns the vector below and not function parameter Sample
  sample <- as.character(Sample)
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

  if (any(is.na(Date))) {
    warning(sum(is.na(Date)),
    " sample/s are widouth collection dates. It is recomendet to exclude
            samples widouth collection dates from anaysis.")
  }

  if (is.null(lat)) {
    stop("Looks like there are some characters in your latitude cordinate notation.")
  }

  if (is.null(lng)) {
    stop("Looks like there are some characters in your longitude cordinate notation.")
  }

  if (anyDuplicated.default(sample) > 0) {
    stop("Duplicated entry in `Sample` vector: ", sample[duplicated(sample)])
  }

  if(any(is.na(AnimalRef))) {
    stop(sum(is.na(AnimalRef)), "\n samples are not assigned to any individuals.
         All samples shuld be assigened to individuals trough AnimalRef column.
         For explanation ?check_sampledata")
  }

  # if (!all(AnimalRef %in% Sample)) {
  #   stop("Individual/s: ", paste(setdiff(unique(AnimalRef), Sample),collapse = ", "), " have names
  #        that are not found as sample identifiers. Individuals have to be named
  #        by one of its samples. For explanation ?check_sampledata" )
  # }

  if (!inherits(Date, "Date")) {
    stop("The specified `Date` vector is not in Date format")
  }

  sex <- c("F", "M", NA)
  if (!all(GeneticSex %in% sex)) {
    stop("It seems that the GeneticSex column is not coded properly.
    You must use 'M' for males, 'F' for females and NA for unknown sex.")
  }


  if (!all(-90 <= lat[!is.na(lat)] & 90 >= lat[!is.na(lat)])) {
    stop("It seems that the latitude column stores coordinates that are not in WGS84 coordinate system.")
  }

  if (!all(-180 <= lng[!is.na(lng)] & 180 >= lng[!is.na(lng)])) {
    stop("It seems that the latitude column stores coordinates that are not in WGS84 coordinate system.")
  }

  if(any(is.na(lat)|any(is.na(lng)))) {
    warning("NA values in lng and/or lat column. When creating spatial
            representation of pedigree such samples will be excluded.")
  }

  sampledata <- data.frame(
    Sample = sample,
    Date = as.Date(Date),
    AnimalRef = AnimalRef,
    GeneticSex = GeneticSex,
    lat = lat,
    lng = lng,
    SType = SType
  )

  if(!is.null(extraCols)) {
    df_name <- sub("\\$.*", "", deparse(substitute(Sample)))
    addCols <- get(df_name)[extraCols]

    sampledata <- cbind(sampledata, addCols)
  }

  return(sampledata)
}
