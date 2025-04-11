#' Temporal plot of pedigree
#'
#' @description Creates "capture" history plot of individuals
#' arranged by families included in data frame created by [`plot_table()`] function.
#'
#' @param plottable Data frame. Output of [`plot_table()`] function.
#' @param famSpacing Y-axis spacing between families. Should be even number!
#' @param pClustSpacing Y-axis spacing between families. Should be even number!
#' @param xWhiteSpace Spacing on the X-axis at the beginning and end of the plot.
#' @param xlabel X-axis label.
#' @param ylabel Y-axis label.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param LegendLabel Title of the legend.
#' @param xlegend Horizontal position of the legend.
#' @param ylegend Vertical position of the legend.
#' @param text_size Plot text size.
#' @param fam_label_size Family label text size.
#'
#' @return A graphical representation of detected family members trough time.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_text geom_hline
#' @importFrom ggplot2 geom_label ggtitle labs expand_limits theme theme_bw
#' @importFrom ggplot2 scale_colour_brewer ylab xlab scale_x_date
#' @export
#'
#' @examples
#' # Prepare the data for usage with plot_table() function.
#' # Get animal timespan data using the anim_timespan() function.
#' animal_ts <- anim_timespan(wolf_samples$AnimalRef,
#'   wolf_samples$Date,
#'   wolf_samples$SType,
#'   dead = c("Tissue")
#' )
#' # Add animal timespan to the sampledata
#' sampledata <- merge(wolf_samples, animal_ts, by.x = "AnimalRef", by.y = "ID", all.x = TRUE)
#' # Define the path to the pedigree data file.
#' path <- paste0(system.file("extdata", package = "wpeR"), "/wpeR_samplePed")
#' # Retrieve the pedigree data from the get_colony function.
#' ped_colony <- get_colony(path, sampledata, rm_obsolete_parents = TRUE, out = "FamAgg")
#' # Organize families and expand pedigree data using the org_fams function.
#' org_tables <- org_fams(ped_colony, sampledata, output = "both")
#' # Prepare data for plotting.
#' pt <- plot_table(plot_fams = 1,
#'   org_tables$fams,
#'   org_tables$ped,
#'   sampledata,
#'   deadSample = c("Tissue")
#' )
#'
#' # Run the function.
#' # Get a temporal pedigree plot.
#' ped_satplot(plottable = pt)
#'
#'
#'
#'
#'
ped_satplot <- function(plottable,
                        famSpacing = 2, pClustSpacing = 2,
                        xWhiteSpace = 100,
                        xlabel = "Date", ylabel = "Animal",
                        title = "", subtitle = "",
                        LegendLabel = "Sex", xlegend = 0.2, ylegend = 0.94,
                        text_size = 2.5, fam_label_size = 2) {

  data <- data.frame(date = plottable$Date,
                     animal = plottable$AnimalRef,
                     plottingID = plottable$plottingID,
                     sex = plottable$GeneticSex,
                     fam = plottable$FamID,
                     polyCluster = plottable$polyCluster,
                     isPolygamous = plottable$isPolygamous,
                     rep = plottable$rep,
                     later_rep = plottable$later_rep,
                     dead = plottable$dead)

  # NEED THIS TO DEFINE VARIABLES NOT DEFINED BY FUNCTION
  # ELSE check() RETURNS NOTE no visible binding for global variable
  # solution found https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  first <- polyFirst <- famFirst <- Y <- yaxis <- animal <- plottingID <- sex <- NULL

  famlines <- NULL # family / polycluster separation lines



  # this gymnastics sorts and names y axis for satplot
  ## adds numerical date of fist date
  ## first sample of the animal
  data$first <- vapply(data$animal, function(x)
    as.Date(min(data$date[data$animal == x], na.rm = TRUE), origin = "1970-01-01"),
    numeric(1))
  ## first polycluster of the animal
  data$polyFirst <- vapply(data$polyCluster, function(x)
    as.Date(min(data$date[data$polyCluster == x], na.rm = TRUE), origin = "1970-01-01"),
    numeric(1))
  ## first family of the animal
  data$famFirst <- vapply(data$fam, function(x)
    as.Date(min(data$date[data$fam == x], na.rm = TRUE), origin = "1970-01-01"),
    numeric(1))

  # order data correctly, based on columns created above
  #dataOrdered <- data |>
    #dplyr::arrange(polyFirst, polyCluster, famFirst, fam, first, plottingID, date)
  dataOrdered <- dplyr::arrange(data, polyFirst, polyCluster, famFirst, fam, first, plottingID, date)

  # make Y coordinate for each plottingID
  dataOrdered$yaxis <- ordered(dataOrdered$plottingID,
                               levels = unique(dataOrdered$plottingID))
  levels(dataOrdered$yaxis) <- 1:length(levels(dataOrdered$yaxis))
  dataOrdered$yaxis <- as.numeric(dataOrdered$yaxis)

  # make line/label for the first polycluster/family
  famlines <- rbind(famlines, data.frame(Y = 1,
                                         type = "polyCluster",
                                         fam = dataOrdered$fam[1],
                                         polyCluster = dataOrdered$polyCluster[1]))
  dataOrdered$yaxis <- dataOrdered$yaxis + famSpacing

  # make Y-axis gaps between polyclusters/families, prepare data for drawing lines
  for (i in 1:(nrow(dataOrdered) - 1)) {
    if (dataOrdered$fam[i] != dataOrdered$fam[i + 1]) {
      dataOrdered$yaxis[(i + 1):nrow(dataOrdered)] <- dataOrdered$yaxis[(i + 1):nrow(dataOrdered)] + famSpacing
      plY <- dataOrdered$yaxis[i] + famSpacing / 2
      plType <- "Family"
      fam <- dataOrdered$fam[i + 1]
      polyCluster <- dataOrdered$polyCluster[i + 1]

      if (dataOrdered$polyCluster[i] != dataOrdered$polyCluster[i + 1]) {
        dataOrdered$yaxis[(i + 1):nrow(dataOrdered)] <- dataOrdered$yaxis[(i + 1):nrow(dataOrdered)] + pClustSpacing
        plY <- dataOrdered$yaxis[i] + famSpacing / 2
        plType <- "polyCluster"
        fam <- dataOrdered$fam[i + 1]
        polyCluster <- dataOrdered$polyCluster[i + 1]
      }
      famlines <- rbind(famlines, data.frame(Y = plY, type = plType, fam, polyCluster))
    }
  }

  dataOrdered$first_sample <- rep(FALSE, nrow(dataOrdered))
  # mark first samplefor labelling
  for (i in 1:(nrow(dataOrdered))) {
    minDt <- min(dataOrdered$date[dataOrdered$plottingID == dataOrdered$plottingID[i]], na.rm = TRUE)
    if (is.na(minDt)) minDt <- Inf
    if (!is.na(dataOrdered$date[i]) & minDt != Inf) {
      if (dataOrdered$date[i] == minDt) dataOrdered$first_sample[i] <- TRUE
    }
  }

  minDate <- min(dataOrdered$date) # minimum date for plotting family names

  p <- ggplot(aes(x = as.Date(date), y = yaxis), data = dataOrdered) +
    geom_line(
      aes(color = sex, group = plottingID),
      alpha = 0.5
      ) +
    geom_point(
      aes(color = sex),
      size = 1
      ) +
    geom_point(
      data = dataOrdered[dataOrdered$rep == TRUE, ],
               aes(y = yaxis, x = as.Date(date)),
               shape = 0, size = 3, color = "red"
      ) +
    geom_point(
      data = dataOrdered[dataOrdered$isPolygamous == TRUE, ],
               aes(y = yaxis, x = as.Date(date)),
               shape = 5, size = 2, color = "purple") +
    geom_point(
      data = dataOrdered[dataOrdered$later_rep == TRUE, ],
               aes(y = yaxis, x = as.Date(date)),
               shape = 1, size = 3, color = "green") +
    geom_point(
      data = dataOrdered[dataOrdered$dead == TRUE, ],
               aes(y = yaxis, x = as.Date(date)),
               shape = 4, size = 3, color = "black") +
    geom_text(
      data = dataOrdered[dataOrdered$first_sample == TRUE, ],
      aes(y = yaxis, x = as.Date(date), label = animal),
      size = text_size, hjust = 1, vjust = 0.5, nudge_x = -15
    ) + # , label.padding = unit(0.1, "lines"))+
    ggtitle(title, subtitle = subtitle) +
    labs(colour = LegendLabel) +
    expand_limits(x = c(min(dataOrdered$date) - xWhiteSpace, max(dataOrdered$date + xWhiteSpace))) +
    theme_bw() +
    scale_colour_brewer(palette = "Set1") +
    ylab(ylabel) +
    xlab(xlabel) +
    scale_x_date(date_labels = ("%m-%Y")) +
    theme(legend.position.inside = c(xlegend, ylegend), legend.direction = "horizontal")



  if (!is.null(famlines)) {
    p <- p + geom_hline(
      yintercept = famlines$Y[famlines$type == "polyCluster"],
      color = "yellow", linewidth = 1
      ) +
      geom_hline(
        yintercept = famlines$Y,
        linetype = "dashed", linewidth = 0.3
        ) +
      geom_label(
        data = famlines, aes(x = rep(minDate, nrow(famlines)),
                             y = Y,
                             label = paste("PCL:", polyCluster, "FAM:", fam, sep = "")),
        size = fam_label_size, color = "darkgreen",
        hjust = 0.5, vjust = 0.5, fontface = "bold"
      )
  }


  # print(p)
  return(p)
}
