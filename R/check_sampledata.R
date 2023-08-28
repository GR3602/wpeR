#' Check and Prepare Sample Data
#'
#' This function checks the integrity and compatibility of various columns in the sample data
#' and prepares a well-structured sample data frame. It verifies the consistency of columns
#' such as Sample, Date, AnimalRef, GeneticSex, IsAnimalReference, lat, lng, and SType.
#' The function ensures that the provided data is properly formatted and conforms to the required standards.
#'
#' @param Sample A vector of sample identifiers.
#' @param Date A vector of sample dates in 'YYYY-MM-DD' format.
#' @param AnimalRef A vector of animal reference identifiers.
#' @param GeneticSex A vector of genetic sex information ('F' for female, 'M' for male, NA for unknown).
#' @param IsAnimalReference A vector indicating whether the sample is an animal reference (0 or 1).
#' @param lat A vector of latitude coordinates in the WGS84 coordinate system.
#' @param lng A vector of longitude coordinates in the WGS84 coordinate system.
#' @param SType A vector of sample types.
#'
#' @return A well-structured sample data frame with validated and formatted columns.
#'
#' @examples
#' sampledata <- check_sampledata(Sample = c("S1", "S2"),
#'                                Date = c("2023-01-15", "2023-02-20"),
#'                                AnimalRef = c("A1", "A2"),
#'                                GeneticSex = c("F", "M"),
#'                                IsAnimalReference = c(1, 0),
#'                                lat = c(34.0522, 40.7128),
#'                                lng = c(-118.2437, -74.0060),
#'                                SType = c("Type1", "Type2"))
#'
#' @export

check_sampledata <- function(Sample, Date, AnimalRef, GeneticSex, IsAnimalReference, lat, lng, SType) {

  n = length(Sample)
  if(n == 0) {
    stop("sample vector has lenght 0")
  }

  if(length(Date) != n) {
    stop(sprintf("Incompatible input: length(Sample) = %d, but length(Date) = %d", n, length(Date)))
  }

  if(length(AnimalRef) != n) {
    stop(sprintf("Incompatible input: length(Sample) = %d, but length(AnimalRef) = %d", n, length(AnimalRef)))
  }

  if(length(GeneticSex) != n) {
    stop(sprintf("Incompatible input: length(Sample) = %d, but length(GeneticSex) = %d", n, length(GeneticSex)))
  }

  if(length(IsAnimalReference) != n) {
    stop(sprintf("Incompatible input: length(Sample) = %d, but length(IsAnimalReference) = %d", n, length(IsAnimalReference)))
  }

  if(length(lat) != n) {
    stop(sprintf("Incompatible input: length(Sample) = %d, but length(lat) = %d", n, length(lat)))
  }

  if(length(lng) != n) {
    stop(sprintf("Incompatible input: length(Sample) = %d, but length(lng) = %d", n, length(lng)))
  }

  if(length(SType) != n) {
    stop(sprintf("Incompatible input: length(Sample) = %d, but length(SType) = %d", n, length(SType)))
  }


  Sample = as.character(Sample)
  Date <- try(as.Date(Date, tryFormats ="%Y-%m-%d"), silent = TRUE)
  AnimalRef = as.character(AnimalRef)
  GeneticSex = as.character(GeneticSex)
  IsAnimalReference <- tryCatch(as.numeric(IsAnimalReference), warning=function(cond) return(NULL))
  lat <- tryCatch(as.numeric(lat), warning=function(cond) return(NULL))
  lng <- tryCatch(as.numeric(lng), warning=function(cond) return(NULL))
  SType <- as.character(SType)

  if(inherits(Date, "try-error")) {
    stop("Wrong date format, date shuld be formated as 'YYYY-MM-DD'")
  }

  if(is.null(IsAnimalReference)){
    stop("IsAnimalReference has to store numerical values, eather 0 or 1")
  }

  if(is.null(lat)){
    stop("Looks like tere are some characters in your latitude corrdinate notation.")
  }

  if(is.null(lng)){
    stop("Looks like tere are some characters in your longitude corrdinate notation.")
  }

  if(anyDuplicated.default(Sample) > 0){
    stop("Duplicated entry in `Sample` vector: ", Sample[duplicated(Sample)])
  }

  if(!inherits(Date, "Date")) {
    stop("The specified `Date` vector is not in Date format")
  }

  sex = c("F", "M", NA)
  if(!all(GeneticSex %in% sex)) {
    stop("It seems that the GeneticSex column is not coded properly.
    You must use 'M' for males, 'F' for females and NA for unknown sex.")
  }

  animrefopt =c(0,1)
  if(!all(IsAnimalReference %in% animrefopt)) {
    stop("It seems that the IsAnimRef column is not coded properly.
    You must use '0' for samples that are not reference and '1' for reference samples.")
  }

  if(!all(-90 <= lat & 90 >= lat)){
    stop("It seems that the latitude column stores coordinates that are not in WGS84 coordinate system.")
  }

  if(!all(-180 <= lng & 180 >= lng)){
    stop("It seems that the latitude column stores coordinates that are not in WGS84 coordinate system.")
  }

  sampledata = data.frame(
    Sample = Sample,
    Date = Date,
    AnimalRef = AnimalRef,
    GeneticSex = GeneticSex,
    IsAnimalReference = IsAnimalReference,
    lat = lat,
    lng = lng,
    SType = SType
  )

  return(sampledata)

}





