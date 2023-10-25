#' Number of detected animals between two sampling seasons
#'
#' @description Gives an numeric overview of individuals captured
#' within the second sampling season compared tho the first one.
#'
#'
#'
#' @param animal_id A column in the dataframe of all samples that stores
#' individual animal identifier code.
#' @param capture_date A column in the dataframe of all samples that stores
#' the date of sample collection. Must be in `Date` format.
#' @param season1_start String in `Date` format. Start of fist capture season.
#' Start and end date are included in the capture season.
#' @param season1_end String in `Date` format. End of fist capture season.
#' Start and end date are included in the capture season.
#' @param season2_start String in `Date` format. Start of second capture season.
#' Start and end date are included in the capture season.
#' @param season2_end String in `Date` format. End of second capture season.
#' Start and end date are included in the capture season.
#'
#' @return
#' A data frame with one row and six columns corresponding to season 1 and 2
#' start and end dates, number of detected animals in season 2 (`total_cap`),
#' number of new detentions in season 2 (`new_captures`), umber of animals from
#' season 1 detected within season 2 (`recaptured`) and number of individuals
#' skipped in season 2 but detected after the end of that season (`skipped`).
#'
#' @export
#'
#' @examples
#' # Calculate the number of animals detected between two sampling seasons.
#' nbtw_seasons(
#'  animal_id = wolf_samples$AnimalRef,
#'  capture_date = wolf_samples$Date,
#'  season1_start = as.Date("2017-01-01"),
#'  season1_end = as.Date("2017-12-31"),
#'  season2_start = as.Date("2018-01-01"),
#'  season2_end = as.Date("2018-12-31")
#' )
#'
#'
#' @aliases nbtw_seasons nBetweenSeasons
#'

nbtw_seasons <- function(animal_id, capture_date,
                         season1_start, season1_end,
                         season2_start, season2_end) {

  dates <- c(season1_start, season1_end, season2_start, season2_end)


  if (!inherits(dates, "Date")) {
    stop("Season starts and ends must be in Date format.")
  }

  if (!all(order(dates) == c(1:4))) {
    warning("Season dates are not defined in correct order.")
  }

  # not captured in the first season
  new_captures <- 0

  # captured before in the first season
  recaptures <- 0

  # captured in the first season, not captured in the second season,
  # captured later
  skipped <- 0

  capdata <- data.frame(animal_id, capture_date)

  season2data <- capdata[which(capture_date >= season2_start &
                          capture_date <= season2_end),]

  season1data <- capdata[which(capture_date >= season1_start &
                          capture_date <= season1_end),]

  post_seasondata <- capdata[which(capture_date > season2_end),]

  # total captures in season 2
  total_cap <- length(unique(season2data$animal_id))

  animals <- unique(season2data$animal_id)

  # all animals. count new captures != recaptures
  for (animal in animals) {
    if (animal %in% season1data$animal_id) {
      recaptures <- recaptures + 1
    } else {
      new_captures <- new_captures + 1
    }
  }

  # skipped season 2
  s1_animals <- unique(season1data$animal_id)

  for (animal in s1_animals) {
    if (animal %in% post_seasondata$animal_id) {
      if (!(animal %in% season2data$animal_id)) {
        skipped <- skipped + 1
      }
    }
  }

  return(data.frame(
    season1 = paste(season1_start, season1_end, sep = " - "),
    season2 = paste(season2_start, season2_end, sep = " - "),
    total_cap, new_captures, recaptures, skipped
  ))
}
