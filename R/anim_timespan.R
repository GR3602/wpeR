#' Get dates of individuals first and last sample
#'
#' @description
#' Takes data frame of all samples and returns the dates of individuals first
#' and last sample.
#' Besides that the functions determines if animal is dead based on a column that flags
#' samples representing mortality, (eg. all tissue samples are taken from dead
#' animals)
#'
#' @param individual_id
#'   Column in the dataframe of all samples containing
#'   individual animal identifier code.
#'   Defined as `dataframe$column`.
#' @param sample_date
#'   Column in the dataframe of all samples containing
#'   the date of sample collection.
#'   Must be in `Date` format. Defined as `dataframe$column`.
#' @param mortality_sample
#'   Logical vector or column in the dataframe of all samples that
#'   identifies samples that represent a mortality event (e.g. from a dead animal).
#'   `TRUE` values mark mortality samples.
#'   Defined as `dataframe$column` or a standalone logical vector the same length
#'   as `individual_id`.
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
#'   mortality_sample = wolf_samples$IsMortality
#' )


anim_timespan <- function(individual_id, sample_date, mortality_sample) {

  if (!any(mortality_sample, na.rm = TRUE)) {
    message("No mortality samples flagged in the dataset, all individuals will be treated as still alive.")
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

    out_individual <- c(out_individual, individual)

    first <- min(ind_dates, na.rm = TRUE) # to sort out NA-only individuals
    if (first == Inf | first == -Inf) first <- NA

    last <- max(ind_dates, na.rm = TRUE)
    if (last == Inf | last == -Inf) last <- NA

    out_firstseen <- c(out_firstseen, first)
    out_lastseen <- c(out_lastseen, last)

    out_isdead <- c(out_isdead, any(mortality_sample[individual_id == individual], na.rm = TRUE))
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
