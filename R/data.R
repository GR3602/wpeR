#'Wolf monitoring genetic samples metadata
#'
#'Metadata of selected genetic samples of wolves collected between 2015 and 2021,
#'in the scope of Slovenian National Wolf Monitoring
#'
#'@format
#'A data frame with 407 rows and 8 columns:
#'\describe{
#'   \item{Sample}{Sample unique identifier code}
#'   \item{Date}{Date of sample collection (format: `YYYY-MM-DD`)}
#'   \item{AnimalRef}{Identification string for particular animal}
#'   \item{GeneticSex}{Sex of animal to which the sample belong (format: `M` = male, `F` = female)}
#'   \item{lat}{latitude (N-S) of the sample (CRS: WGS84; EPSG: 4326)}
#'   \item{lng}{longitude (W-E) of the sample (CRS: WGS84; EPSG: 4326)}
#'   \item{SType}{Type of the sample. (`Direct Saliva, Scat, Urine, Saliva, Tissue, Decomposing Tissue, Blood`) }
#'   \item{IsMortality}{ A logical column (TRUE/FALSE) that marks samples collected from dead animals.
#'   TRUE = mortality event (e.g., a tissue sample from a carcass),
#'   FALSE = non-mortality sample (e.g., scat, hair, saliva from a live animal).}
#'   }
#'
#'@source Slovenian National Wolf Monitoring
"wolf_samples"


