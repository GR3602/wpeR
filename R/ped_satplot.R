#' Mark-Recapture Plot For Pedigree - Saturating
#'
#' @description Creates "capture" history plot of individuals
#' arranged by families included in data frame created by [`plot_table()`] function.
#'
#' @param date Sample collection date. Must be in `Date` format.
#' @param animal Individual animal identifier code.
#' @param plottingID ID for plotting, because of polygamy same animal can be in
#'  more families.
#' @param sex Sex, but can be other parameters - you can set up LegendLabel.
#' @param fam Family identifier number.
#' @param polyCluster Cluster for polygamous animals, includes all families in
#' which an individual is reproductive animal.
#' @param isPolygamous Logical value (`TRUE`/`FALSE`) that identifies
#' animal as polygamous (has >1 mates).
#' @param rep Logical value (`TRUE`/`FALSE`) that determines if individual
#' is reproductive or not in current family.
#' @param later_rep Logical value (`TRUE`/`FALSE`) that determines if individual
#' is reproductive or not in any other (later) family.
#' @param dead Logical value (`TRUE/FALSE`) that identifies if the individual is dead.
#' @param famSpacing Y-axis spacing between families. Should be even number!
#' @param pClustSpacing Y-axis spacing between families. Should be even number!
#' @param xWhiteSpace fill in
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
#' @return a graphical representation of detected family members trough time.
#' @import ggplot2
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
#' pt <- plot_table(org_tables$fams[1, ],
#'   org_tables$fams,
#'   org_tables$ped,
#'   sampledata,
#'   deadSample = c("Tissue", "Decomposing Tissue", "Blood")
#' )
#'
#' # Run the function.
#' # Get a temporal pedigree plot.
#' ped_satplot(
#'   pt$Date,
#'   pt$AnimalRef,
#'   pt$plottingID,
#'   pt$GeneticSex,
#'   pt$FamID,
#'   pt$polyCluster,
#'   pt$isPolygamous,
#'   pt$rep,
#'   pt$later_rep,
#'   pt$dead
#' )
#'
#' @aliases ped_satplot PedigreeCMRSatplot
#'
#'
#'
#'
ped_satplot <- function(date,
                        animal,
                        plottingID,
                        sex,
                        fam,
                        polyCluster,
                        isPolygamous,
                        rep,
                        later_rep,
                        dead,
                        famSpacing = 2, pClustSpacing = 2,
                        xWhiteSpace = 100,
                        xlabel = "Date", ylabel = "Animal",
                        title = "Pedigree CMR Graph", subtitle = "",
                        LegendLabel = "Sex", xlegend = 0.2, ylegend = 0.94,
                        text_size = 2.5, fam_label_size = 2) {

  data <- data.frame(date,
                     animal,
                     plottingID,
                     sex,
                     fam,
                     polyCluster,
                     isPolygamous,
                     rep,
                     later_rep,
                     dead)

  # NEED THIS TO DEFINE VARIABLES NOT DEFINED BY FUNCTION
  # ELSE check() RETURNS NOTE no visible binding for global variable
  # solution found https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  first <- polyFirst <- famFirst <- Y <- yaxis <- NULL

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
  dataOrdered <- data |>
    dplyr::arrange(polyFirst, polyCluster, famFirst, fam, first, plottingID, date)

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
    theme(legend.position = c(xlegend, ylegend), legend.direction = "horizontal")



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
