#'Wolf Monitoring Genetic Samples Data
#'
#'Metadata of selected genetic samples of wolves collected between XXXX and YYYY,
#'in the scope of Slovenian Natioanl Wolf Monitoring
#'#'@format ##`wolf_samples`
#'A data frame with 407 rows and 8 columns:
#'\describe{
#'   \item{Sample}{Sample unique identifier code}
#'   \item{Date}{Date of sample collection (format: `YYYY-MM-DD`)}
#'   \item{AnimalRef}{Identification string for patricular animal}
#'   \item{GeneticSex}{Sex of animal to which the sample belong (format: `M` = male, `F` = female)}
#'   \item{IsAnimalReference}{Is a samle also a reference sample for particular animal (format: `1` = yes, `0` = no)}
#'   \item{lat}{latitude (N-S) of the sample (CRS: WGS84; EPSG: 4326)}
#'   \item{lng}{longitude (W-E) of the sample (CRS: WGS84; EPSG: 4326)}
#'   \item{SType}{Type of the sample. (`Direct Saliva, Scat, Urine, Saliva, Tissue, Decomposing Tissue, Blood`) }
#'}
#'
#'@source Slovenian National Wolf Monitoring
"wolf_samples"
