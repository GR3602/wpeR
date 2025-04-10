#' Get matrix ff apparent survival
#'
#' @description
#' Creates a matrix that shows number of captured animals
#' between multiple seasons.
#'
#' @param animal_id A column in the dataframe of all samples that stores
#' individual animal identifier code.
#' @param capture_date A column in the dataframe of all samples that stores
#' the date of sample collection. Must be in `Date` format.
#' @param start_dates Vector of dates in `Date` format that define the
#' start of each season.
#' @param end_dates Vector of dates in `Date` format that define the
#' end of each season.
#'
#' @return
#' A matrix with 1 + no. seasons rows and columns.
#'  * diagonal: number of new captures in each session,
#'  * above diagonal: number of recaptures from season x to season y,
#'  * below diagonal: number of animals from season y that skipped season x.
#'
#' Season x is defined in first row, season y in first column.
#' Column `Tot. Capts` gives all detected individuals in season y.
#' Row `Tot. Skipped` gives all individuals skipped in season x but detected
#' later.
#'
#' @export
#'
#' @examples
#' # Define start and end dates for sampling seasons.
#' seasons <- data.frame(
#'   start = c(
#'     as.Date("2017-01-01"),
#'     as.Date("2018-01-01"),
#'     as.Date("2019-01-01")
#'   ),
#'   end = c(
#'     as.Date("2017-12-31"),
#'     as.Date("2018-12-31"),
#'     as.Date("2019-12-31")
#'   )
#' )
#'
#' # Create a dynamics matrix for animal captures.
#' dyn_matrix(
#'   animal_id = wolf_samples$AnimalRef,
#'   capture_date = wolf_samples$Date,
#'   start_dates = seasons$start,
#'   end_dates = seasons$end
#' )
#'
#'

dyn_matrix <- function(animal_id, capture_date, start_dates, end_dates) {

  seasons <- data.frame(start_date = start_dates, end_date = end_dates)
  mtx_dim <- nrow(seasons) + 1

  seasons_names <- paste(seasons$start_date, seasons$end_date, sep = " - ")
  mtx_dimnames <- list(c(seasons_names, "Tot. Skipped"),
                       c(seasons_names, "Tot. Capts"))
  outmatrix <- matrix(NA,
                      nrow = mtx_dim,
                      ncol = mtx_dim,
                      dimnames = mtx_dimnames)


  for (i in seq_along(start_dates)) {
    # first, the total numbers for season i
    n_season <- suppressWarnings(
      nbtw_seasons(
        animal_id, capture_date,
        seasons[1, ]$start_date, seasons[i, ]$start_date - 1,
        seasons[i, ]$start_date, seasons[i, ]$end_date
      )
    )

    outmatrix[i, i] <- n_season$new_captures
    outmatrix[i, mtx_dim] <- n_season$total_cap
    outmatrix[mtx_dim, i] <- n_season$skipped

    if (i == length(start_dates)) break # get out of loop for the last element

    for (j in (i + 1):length(start_dates)) {
      n_season <- nbtw_seasons(
        animal_id, capture_date,
        seasons[i, ]$start_date, seasons[i, ]$end_date,
        seasons[j, ]$start_date, seasons[j, ]$end_date
      )

      outmatrix[i, j] <- n_season$recaptures
      outmatrix[j, i] <- n_season$skipped
    }
  }

  return(outmatrix)
}
