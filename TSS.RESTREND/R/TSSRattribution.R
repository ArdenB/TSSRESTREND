#' @title Vegetation change attribution using the Time Series Segmentation of Residual Trends (MAIN FUNCTION)
#'
#' @importFrom stats coef end frequency lm sd start time ts cor.test
#' @importFrom mblm mblm
#' @importFrom graphics abline arrows legend par plot
#' @importFrom utils tail read.csv
#' @importFrom broom glance
#'
#' @description
#' This is a wrapper function around the TSS-RESTREND main function that dows additional attribution.
#' It measues the Observed vegetation change, land use, climate change and climate varibility.
#' Unlike TSSRESTREND function, this requires both temperature and precitation data fo work.
#' @author Arden Burrell, aburrell@whrc.org
#'
#' @inheritParams TSSRESTREND
#' @param C4frac
#'        The fraction of vegetation that follows the C4 photosynthetic pathway, between 0 and 1
#' @param verbose
#'        Return all the created models as well as the original data
#' @param AnnualRes
#'        Report results in change per year. Defualt is False. Instead reports total change from the
#'        start to the end of the time series.
#'
#' @return \bold{tacp}
#'        The accumulation period for the annual max time series temperature


#' @export

TSSRattribution <- function(
	CTSR.VI, CTSR.RF, CTSR.TM, C4frac = 0,  sig = 0.05, season = "none", exclude = 0,
	allow.negative = FALSE, allowneg.retest = FALSE, h = 0.15, verbose = FALSE, AnnualRes=FALSE){

  # ========== Perfrom the data check ==========
  if (any(is.na(CTSR.VI))){
    return(NA)
  }else if (any(is.na(CTSR.RF))){
    return(NA)
  }else if (any(is.na(CTSR.TM))){
    return(NA)
  }
  # ========== Get the annual Max VI values ==========
  rawmax.df  <- AnMaxVI(CTSR.VI)
  CTSR.VIadj <- franksCO2(CTSR.VI, C4frac)
  adjmax.df  <- AnMaxVI(CTSR.VIadj)

  # +++++ the veg change due to CO2 +++++
  CO2change  <-  rawmax.df$Max - adjmax.df$Max
  ti <- time(CO2change)
  CO2cor <- cor.test(ti, CO2change, method="spearman")
  CO2theil <-  mblm(CO2change~ti, repeated=FALSE)


	# ========== Calculate the accumulated temperature and Precipitation ==========
  # Add a new temperature info so i can calculate the mean 
	browser()

}
