#' Get dates of individuals first and last sample
#'
#' @description
#' Takes data frame of all samples and returns the dates of individuals first and last sample.
#' Besides that the functions determines if animal is dead based on predefined sample type eg. tissue.
#'
#'
#' @param individual_id column in the dataframe of all samples that stores individual animal identifier code.
#'   Defined as `dataframe$column`.
#' @param sample_date column in the dataframe of all samples that stores the date of sample collection.
#'   Must be in `Date` format. Defined as `dataframe$column`.
#' @param sample_type column in the dataframe of all samples that stores the data on the type (eg. scat, tissue, saliva) of particular sample.
#'   Defined as `dataframe$column`.
#' @param dead single value or vector of different lethal sample types. Defaults to "Tissue".
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
#' anim_timespan(pack21_samples$AnimalRef,
#'                      pack21_samples$Date,
#'                      pack21_samples$SType,
#'                      dead = c("Tissue", "Decomposing Tissue", "Blood"))
#'
#' @aliases make.animal.timespan
#'

anim_timespan <- function(individual_id, sample_date, sample_type, dead = "Tissue") {

  unique_ind <- unique(individual_id)
  out_individual <- NULL
  out_firstseen <- NULL
  out_lastseen <- NULL
  out_isdead <- NULL

  for (individual in unique_ind) {
    ind_dates <- sample_date[individual_id == individual]
    ind_sample_types <- sample_type[individual_id == individual]

    out_individual <- c(out_individual, individual)

    first <- min(ind_dates, na.rm = T) # to sort out NA-only individuals
    if (first == Inf | first == -Inf) first <- NA

    last <- max(ind_dates, na.rm = T)
    if (last == Inf | last == -Inf) last <- NA

    out_firstseen <- c(out_firstseen, first)
    out_lastseen <- c(out_lastseen, last)

    out_isdead <- c(out_isdead, sum(dead %in% ind_sample_types) > 0)
  }

  out_firstseen <- as.Date(out_firstseen, origin = "1970-01-01")
  out_lastseen <- as.Date(out_lastseen, origin = "1970-01-01")

  return(data.frame(ID = out_individual, FirstSeen = out_firstseen, LastSeen = out_lastseen, IsDead = out_isdead))
}
