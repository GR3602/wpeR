#' Get dates of individuals first and last sample
#'
#' @description
#' Takes data frame of all samples and returns the dates of individuals first
#' and last sample.
#' Besides that the functions determines if animal is dead based on predefined
#' sample type eg. tissue.
#'
#'
#'
#' @param individual_id
#'   Column in the dataframe of all samples containing
#'   individual animal identifier code.
#'   Defined as `dataframe$column`.
#' @param sample_date
#'   Column in the dataframe of all samples containing
#'   the date of sample collection.
#'   Must be in `Date` format. Defined as `dataframe$column`.
#' @param sample_type
#'   Column in the dataframe of all samples containing the data
#'   on the type (eg. scat, tissue, saliva) of particular sample.
#'   Defined as `dataframe$column`.
#' @param dead
#'   Single value or vector of different lethal sample types. If no lethal
#'   samples are included in the sampledata the dead parameter can be set to `FALSE` (dead = `FALSE`).
#'   Defaults to "Tissue".
#'
#' @return
#' A data frame with four columns and one row for each `individual_id`.
#' Returned data frame columns correspond to individual identification key (`ID`),
#' date of first (`FirstSeen`) and last (`LastSeen`) sample of individual and
#' logical (`TRUE/FALSE`) value that identifies if the individual is dead (`IsDead`).
#'
#' @export
#'
#' @examples
#' anim_timespan(
#'   individual_id = wolf_samples$AnimalRef,
#'   sample_date = wolf_samples$Date,
#'   sample_type = wolf_samples$SType,
#'   dead = c("Tissue")
#' )



anim_timespan <- function(individual_id, sample_date, sample_type, dead = "Tissue") {

  if (any(dead != FALSE)){
    if (!all(dead %in% unique(sample_type), na.rm = TRUE)) {
      warning("one or more lethal sample type, defined with dead parameter,
      are not present in the sample_type vector")
    }
  }

  unique_ind <- unique(individual_id)
  if (any(is.na(unique_ind))) {
    warning("looks like you have some missing values in the individual id column.\n",
    "It is good practice that all samples have assigned individual\n. Validate your input dataframe with ?check_sampledata")
  }

  out_individual <- NULL
  out_firstseen <- NULL
  out_lastseen <- NULL
  out_isdead <- NULL

  for (individual in unique_ind) {
    ind_dates <- sample_date[individual_id == individual]
    ind_sample_types <- sample_type[individual_id == individual]

    out_individual <- c(out_individual, individual)

    first <- min(ind_dates, na.rm = TRUE) # to sort out NA-only individuals
    if (first == Inf | first == -Inf) first <- NA

    last <- max(ind_dates, na.rm = TRUE)
    if (last == Inf | last == -Inf) last <- NA

    out_firstseen <- c(out_firstseen, first)
    out_lastseen <- c(out_lastseen, last)

    out_isdead <- c(out_isdead, sum(dead %in% ind_sample_types) > 0)
  }

  out_firstseen <- as.Date(out_firstseen, origin = "1970-01-01")
  out_lastseen <- as.Date(out_lastseen, origin = "1970-01-01")

  return(data.frame(
    ID = out_individual,
    FirstSeen = out_firstseen,
    LastSeen = out_lastseen,
    IsDead = out_isdead
  ))
}
