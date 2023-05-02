
#' Mark-Recapture Plot For Pedigree - Saturating
#'
#' @description `PedigreeCMRSatplot` creates "capture" history plot of individuals
#' arranged by families included in data frame created by [`PackTable`] function.
#'
#' @param date sample collection date. Must be in `Date` format.
#' @param animal individual animal identifier code.
#' @param plottingID ID for plotting, because of polygamy same animal can be in more families.
#' @param sex sex, but can be other parameters - you can set up LegendLabel
#' @param pack family identifier number.
#' @param polyCluster cluster for polygamous animals, includes all families in
#' which an individual is reproductive animal
#' @param isPolygamous logical value (`TRUE`/`FALSE`) that identifies
#' animal as polygamous (has >1 mates).
#' @param alpha logical value (`TRUE`/`FALSE`) that determines if individual
#' is reproductive or not in current family.
#' @param later_alpha logical value (`TRUE`/`FALSE`) that determines if individual
#' is reproductive or not in any other (later) family.
#' @param dead logical value (`TRUE/FALSE`) that identifies if the individual is dead
#' @param packSpacing Y-axis spacing between packs. Shuld be even number!
#' @param pClustSpacing Y-axis spacing between packs. Shuld be even number!
#' @param xWhiteSpace fill in
#' @param xlabel x axis label.
#' @param ylabel y axis label.
#' @param title plot title.
#' @param subtitle plot subtitle.
#' @param LegendLabel title of the legend.
#' @param xlegend horizontal position of the legend.
#' @param ylegend vertical position of the legend.
#' @param text_size plot text size.
#' @param pack_label_size family label text size.
#'
#' @return a graphical representation of detected family members trough time.
#' @import ggplot2
#' @export
#'
#' @examples
#' animal_ts <- make.animal.timespan(pack21_samples$AnimalRef,
#'                                   pack21_samples$Date,
#'                                   pack21_samples$SType,
#'                                   dead = c("Tissue", "Decomposing Tissue", "Blood"))
#'
#' sampledata <- merge(pack21_samples, animal_ts, by.x = "AnimalRef", by.y = "ID", all.x = TRUE )
#'
#' path <- paste0(system.file("extdata", package = "WildPedigreeExplorer"), "/fake_colony")
#'
#' ped_colony <- GetDigestColony(path, sampledata, remove_obsolete_parents = TRUE, out = "FamAgg")
#'
#' org_tables <- organizePacks(ped_colony, sampledata, output = "both")
#'
#' pt<-PackTable(org_tables$packs[1,],
#'               org_tables$packs,
#'               org_tables$ped,
#'               sampledata,
#'               deadSample = c("Tissue", "Decomposing Tissue", "Blood"))
#'
#'
#' PedigreeCMRSatplot(pt$Date,
#'                    pt$AnimalRef,
#'                    pt$plottingID,
#'                    pt$GeneticSex,
#'                    pt$PackID,
#'                    pt$polyCluster,
#'                    pt$isPolygamous,
#'                    pt$alpha,
#'                    pt$later_alpha,
#'                    pt$dead)
#'
#'
#'
#'
#'
#'
#'
PedigreeCMRSatplot = function(date, animal, plottingID, sex, pack, polyCluster, isPolygamous, alpha, later_alpha, dead,
                              packSpacing = 2, pClustSpacing = 2,
                              xWhiteSpace = 100,
                              xlabel="Date", ylabel="Animal",
                              title="Pedigree CMR Graph", subtitle = "",
                              LegendLabel = "Sex",xlegend=0.2, ylegend=0.94,
                              text_size = 2.5, pack_label_size = 2) {


  #Mark-Recapture plot for pedigree - saturating
  #Date - sample collection date
  #Animal- each animal (to set up a saturating graph)
  #plottingID - ID for plotting, same animal can be in more lines for polygamy plots.
  #sex - sex, but can be other parameters - you can set up LegendLabel
  #pack - pack (family = father, mother, offspring)
  #polyCluster - cluster for polygamous animals
  #isPolygamous (boolean) - is the animal polygamous (has >1 mates)
  #alpha (boolean) - is the animal Alpha animal within the pack.
  #packSpacing, pClustSpacing - Y-axis spacing between packs/polyclusters.
  #mark_pcl, mark_pack = mark polycluster and pack in the graph
  #HAS A WRAPPER THAT ORGANIZES DATA!


  #require(ggplot2, dplyr)

  data=data.frame(date, animal, plottingID, sex, pack, polyCluster, isPolygamous, alpha, later_alpha, dead)

  #NEED THIS TO DEFINE VARIABLES NOT DEFINED BY FUNCTION
  #ELSE check() RETURNS NOTE no visible binding for global variable
  #solution found https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  first <- polyFirst <- packFirst <- Y <- yaxis <- NULL

  packlines = NULL #pack / polycluster separation lines



  #this gymnastics sorts and names y axis for satplot
  ##adds numerical date of fist date
  ##first sample of the animal
  data$first = sapply(data$animal, function(x) as.Date(min(data$date[data$animal==x], na.rm=T), origin="1970-01-01"))
  ##first polycluster of the animal
  data$polyFirst = sapply(data$polyCluster, function(x) as.Date(min(data$date[data$polyCluster==x], na.rm=T), origin="1970-01-01"))
  ##first pack of the animal
  data$packFirst = sapply(data$pack, function(x) as.Date(min(data$date[data$pack == x], na.rm=T), origin="1970-01-01"))

  #order data correctly, based on clumns created above
  dataOrdered = data |>
    dplyr::arrange(polyFirst, polyCluster, packFirst, pack, first, plottingID, date)

  #make Y coordinate for each plottingID
  dataOrdered$yaxis = ordered(dataOrdered$plottingID, levels=unique(dataOrdered$plottingID))
  levels(dataOrdered$yaxis)=1:length(levels(dataOrdered$yaxis))
  dataOrdered$yaxis = as.numeric(dataOrdered$yaxis)

  #make line/label for the first polycluster/pack

  packlines = rbind(packlines, data.frame(Y = 1, type = "polyCluster", pack=dataOrdered$pack[1], polyCluster = dataOrdered$polyCluster[1]))
  dataOrdered$yaxis = dataOrdered$yaxis + packSpacing

  #make Y-axis gaps between polyclusters/packs, prepare data for drawing lines

  for (i in 1:(nrow(dataOrdered)-1)){

    if(dataOrdered$pack[i] != dataOrdered$pack[i+1]) {
      dataOrdered$yaxis[(i+1):nrow(dataOrdered)] = dataOrdered$yaxis[(i+1):nrow(dataOrdered)] + packSpacing
      plY = dataOrdered$yaxis[i]+packSpacing/2
      plType = "Pack"
      pack = dataOrdered$pack[i+1]
      polyCluster = dataOrdered$polyCluster[i+1]

      if(dataOrdered$polyCluster[i] != dataOrdered$polyCluster[i+1]) {
        dataOrdered$yaxis[(i+1):nrow(dataOrdered)] = dataOrdered$yaxis[(i+1):nrow(dataOrdered)] + pClustSpacing
        plY = dataOrdered$yaxis[i]+packSpacing/2
        plType = "polyCluster"
        pack = dataOrdered$pack[i+1]
        polyCluster = dataOrdered$polyCluster[i+1]
      }
      packlines = rbind(packlines, data.frame(Y = plY, type = plType, pack, polyCluster))
    }
  }

  dataOrdered$first_sample = rep(FALSE, nrow(dataOrdered))
  #mark first samplefor labelling
  for (i in 1:(nrow(dataOrdered))){
    minDt = min(dataOrdered$date[dataOrdered$plottingID == dataOrdered$plottingID[i]], na.rm=T)
    if(is.na(minDt)) minDt = Inf
    if(!is.na(dataOrdered$date[i]) & minDt != Inf){
      if (dataOrdered$date[i] == minDt) dataOrdered$first_sample[i] = TRUE}
  }

  minDate = min(dataOrdered$date) #minimum date for plotting pack names

  #plotting
  ##hasehd in the source code
  #increase X limit
  #mindate = as.numeric(min(dataOrdered$date))
  #maxdate = as.numeric(max(dataOrdered$date))
  #xWhiteSpace = as.numeric(xWhiteSpace)

  # xlimit1 = as.Date(mindate-xWhiteSpace*(maxdate-mindate), origin="1970-01-01")
  # xlimit2 = as.Date(maxdate+xWhiteSpace*(maxdate-mindate), origin="1970-01-01")
  ##
  p=ggplot(aes(x=as.Date(date), y=yaxis), data=dataOrdered) +
    geom_line (aes(color=sex, group=plottingID), alpha = 0.5)+
    geom_point(aes(color=sex), size=1)+
    geom_point(data = dataOrdered[dataOrdered$alpha == TRUE,], aes(y=yaxis, x=as.Date(date)), shape=0, size = 3, color = "red")+
    geom_point(data = dataOrdered[dataOrdered$isPolygamous == TRUE,], aes(y=yaxis, x=as.Date(date)), shape=5, size = 2, color = "purple")+
    geom_point(data = dataOrdered[dataOrdered$later_alpha == TRUE,], aes(y=yaxis, x=as.Date(date)), shape=1, size = 3, color = "green")+
    geom_point(data = dataOrdered[dataOrdered$dead == TRUE,], aes(y=yaxis, x=as.Date(date)), shape=4, size = 3, color = "black")+
    geom_text(data = dataOrdered[dataOrdered$first_sample == TRUE,], aes(y=yaxis, x=as.Date(date), label = animal),
              size = text_size, hjust=1, vjust=0.5, nudge_x = -15)+#, label.padding = unit(0.1, "lines"))+
    ggtitle(title, subtitle=subtitle)+
    labs(colour=LegendLabel)+
    expand_limits(x=c(min(dataOrdered$date)-xWhiteSpace, max(dataOrdered$date+xWhiteSpace)))+
    theme_bw()+scale_colour_brewer(palette = "Set1")+
    ylab(ylabel) + xlab(xlabel)+
    scale_x_date(date_labels = ("%m-%Y"))+
    theme(legend.position=c(xlegend,ylegend), legend.direction = "horizontal")



  if(!is.null(packlines)){
    p = p+geom_hline(yintercept=packlines$Y[packlines$type == "polyCluster"], color = "yellow", size = 1)+
      geom_hline(yintercept=packlines$Y, linetype = "dashed", size = 0.3)+
      geom_label(data=packlines, aes(x=rep(minDate, nrow(packlines)), y=Y, label = paste("PCL:",polyCluster, " PCK:", pack, sep="")),
                 size = pack_label_size, color="darkgreen", #label.padding = unit(1, "lines"),
                 hjust=0.5, vjust=0.5, fontface="bold")

  }


  print(p)
  return(p)

}
