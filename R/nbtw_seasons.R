#' Number Of Detected Animals Between Two Sampling Seasons
#'
#' @description `nbtw_seasons` gives an numeric overview of individuals captured
#' within the second sampling season compared tho the first one.
#'
#'
#'
#' @param animal_id column in the dataframe of all samples that stores individual animal identifier code.
#' @param capture_date column in the dataframe of all samples that stores the date of sample collection.
#'   Must be in `Date` format.
#' @param season1_start string in `Date` format. Start of fist capture season.
#' Start and end date are included in the capture season.
#' @param season1_end string in `Date` format. End of fist capture season.
#' Start and end date are included in the capture season.
#' @param season2_start string in `Date` format. Start of second capture season.
#' Start and end date are included in the capture season.
#' @param season2_end string in `Date` format. End of second capture season.
#' Start and end date are included in the capture season.
#'
#' @return
#' A data frame with one row and six columns corresponding to season 1 and 2 start
#' and end dates, number of detected animals in season 2 (`total_cap`),
#' number of new detentions in season 2 (`new_captures`), umber of animals from
#' season 1 detected within season 2 (`recaptured`) and number of individuals
#' skipped in season 2 but detected after the end of that season (`skipped`).
#'
#' @export
#'
#' @examples
#'
#' nbtw_seasons(pack21_samples$AnimalRef, pack21_samples$Date,
#'                as.Date("2010-01-01"), as.Date("2010-12-31"),
#'                as.Date("2011-01-01"), as.Date("2011-12-31"))
#'
#'
#'
#' @aliases nbtw_seasons nBetweenSeasons
#'

nbtw_seasons = function(animal_id, capture_date,
                           season1_start, season1_end, season2_start, season2_end) {


  #no sex!
  #returns the number of animals from season 1
  #captured within season 2.
  #n_captures - number of animals captured
  #new_captures - number of animals not captured in first season
  #skipped - number of skipped
  #animal_id - vector of animal IDs
  #capture_date - vector of dates of animal captures, POSIXct
  #season_start, season_end, POSIXct
  #season start and end day is INCLUSIVE IN THE SEASON!
  #require(dplyr)

  dates = c(season1_start, season1_end, season2_start, season2_end)


  if(!inherits(dates, "Date")) {
    stop("Season starts and ends must be in Date format.")
  }

  if(!all(order(dates) == c(1:4))){
    warning("Season dates are not defined in correct order.")
  }

  #not captured in the first season
  new_captures = 0

  #captured before in the first season
  recaptures = 0

  #captured in the first season, not captured in the second season, captured later
  skipped = 0

  capdata = data.frame(animal_id, capture_date)

  season2data = capdata[capdata$capture_date >= season2_start &
                          capdata$capture_date <= season2_end, ]

  season1data = capdata[capdata$capture_date >= season1_start &
                          capdata$capture_date <= season1_end, ]

  post_seasondata = capdata[capdata$capture_date > season2_end, ]

  #total captures in season 2
  total_cap = length(unique(season2data$animal_id))

  animals = unique(season2data$animal_id)

  #all animals. count new captures != recaptures
  for (animal in animals) {
    if (animal %in% season1data$animal_id) recaptures = recaptures + 1
    else new_captures = new_captures + 1
  }

  #skipped season 2
  s1_animals = unique(season1data$animal_id)

  for (animal in s1_animals) {
    if (animal %in% post_seasondata$animal_id) {
      if(!(animal %in% season2data$animal_id))
        skipped = skipped + 1
    }
  }

  return(data.frame(season1=paste(season1_start,season1_end, sep =" - "),
                    season2=paste(season2_start,season2_end, sep =" - "),
                    total_cap, new_captures,recaptures, skipped)
  )


} #n.between.seasons END
