#' Get Matrix Of Apparent Survival
#'
#' @description
#' `MakeDynamicsMatrix()` creates a matrix that shows number of captured animals
#' between multiple seasons.
#'
#' @param animal_id column in the dataframe of all samples that stores individual animal identifier code.
#' @param capture_date column in the dataframe of all samples that stores the date of sample collection.
#'   Must be in `Date` format.
#' @param start_dates vector of dates in `Date` format that define the start of each season.
#' @param end_dates vector of dates in `Date` format that define the end of each season.
#'
#' @return
#' A matrix with 1 + no. seasons rows and columns.
#'  * diagonal: number of new captures in each session,
#'  * above diagonal: number of recaptures from season x to season y,
#'  * below diagonal: number of animals from season y that skipped season x.
#'
#' Season x is defined in first row, season y in first column. Column `Tot. Capts`
#' gives all detected individuals in season y. Row `Tot. Skipped` gives all
#' individuals skipped in season x but detected later.
#'
#' @export
#'
#' @examples
#'
#' seasons <- data.frame(start = c(as.Date("2010-01-01"),
#'                                 as.Date("2011-01-01"),
#'                                 as.Date("2012-01-01")),
#'                      end = c(as.Date("2010-12-31"),
#'                              as.Date("2011-12-31"),
#'                              as.Date("2012-12-31"))
#'                              )
#'
#'MakeDynamicsMatrix(pack21_samples$AnimalRef, pack21_samples$Date,
#'                   seasons$start, seasons$end)
#'
#'

MakeDynamicsMatrix = function(animal_id, capture_date, start_dates, end_dates){
  # Makes matrics of apparent survival
  # diagonal: number of new captures in each session
  # above diagonal: number of recaptures from season x to season y
  # below diagonal: number of animals from season y that skipped season x

  seasons = data.frame (start_date = start_dates, end_date = end_dates)
  mtx_dim = nrow(seasons)+1

  seasons_names = paste(seasons$start_date, seasons$end_date, sep = " - ")
  mtx_dimnames = list(c(seasons_names, "Tot. Skipped"), c(seasons_names, "Tot. Capts"))
  outmatrix = matrix(NA, nrow = mtx_dim, ncol = mtx_dim, dimnames=mtx_dimnames)


  for (i in 1:length(start_dates)){
    #first, the total numbers for season i
    n_season = nBetweenSeasons(animal_id, capture_date,
                               seasons[1,]$start_date, seasons[i,]$start_date - 1,
                               seasons[i,]$start_date, seasons[i,]$end_date)

    outmatrix[i,i] = n_season$new_captures
    outmatrix[i,mtx_dim] = n_season$total_cap
    outmatrix[mtx_dim,i] = n_season$skipped

    if (i == length(start_dates)) break #get out of loop for the last element

    for (j in (i+1):length(start_dates)){
      n_season = nBetweenSeasons(animal_id, capture_date,
                                 seasons[i,]$start_date, seasons[i,]$end_date,
                                 seasons[j,]$start_date, seasons[j,]$end_date)

      outmatrix[i,j] = n_season$recaptures
      outmatrix[j,i] = n_season$skipped
    }
  }

  return(outmatrix)

}

