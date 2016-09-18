#
# load("./demo_data/stdRESTREND.Rda")
# load("./demo_data/stdRESTREND_CTSR.Rda")
# load("./demo_data/stdRESTREND_RF.Rda")
#
# load("./demo_data/segRESTREND.Rda")
# load("./demo_data/segRESTREND_CTSR.Rda")
# load("./demo_data/segRESTREND_RF.Rda")
#
# load("./demo_data/segVPRD.Rda")
# load("./demo_data/segVPRD_CTSR.Rda")
# load("./demo_data/segVPRD_RF.Rda")
#
# load("./demo_data/segVPRI.Rda")
# load("./demo_data/segVPRI_CTSR.Rda")
# load("./demo_data/segVPRI_RF.Rda")
#
#
#
# source("RF_acum.R")
# source("max_pos.R")
# source("CTSR_acp.R")
# source("Annual_precipitation.R")
# source("VPR_BFAST.R")
# source("CHOW_bptest.R")
# source("viMAX_plot.R")
# source("segVPR.R")
# source("segRESTREND.R")
# source("RESTREND.R")



#FUnctions to call functions
#Function to call the other functions
#Missing function to find optimal accumulation of the precipitation
#which will be a sperate script. If the method is segVPR, there
#is a need to recalculate precip on either side of the breakpoint
#Until it is functional  b4 and after need to be passed into this
#function.  rf.b4=FALSE, rf.af=FALSE, will be removed as soon as'
#'
#' @title A script containing all my draft paramters
#'
#' From TSS_RESTREND
#' @param CTSR.VI
#' Complete Time Series of Vegetation Index. An object of class \code{'ts'}. Monthly time series of VI values
#' @param ACP.table
#' A table of every combination of offset period and accumulation period. if ACP.table = FALSE, CTSR.RF and acu.RF must be provided
#' @param CTSR.RF
#' Complete Time Series of Rain Fall.  An object of class 'ts' and be the same length and cover the same time range as CTSR.VI.
#' If ACU.table is provided, CTSR.RF will be automitaclly calculated by \code{\link{rainfall.accumulator}}
#' @param anu.VI
#' The annual (Growing season) max VI. if anu.VI=FALSE, it will be calculated from the CTSR.VI using \code{\link{AnMax.VI}}.
#' @param acu.RF
#' The optimal accumulated rainfall for anu.VI. Mut be a object of class 'ts' and of equal length to anu.VI. if anu.RF=FALSE, it will be calculated from ACP.table. see(____)
#' @param VI.index
#' the index of the CTSR.VI ts that the anu.VI values occur at. Must be the same length as anu.VI. NOTE. R indexs from 1 rather than 0.
#' NOTE. if VI.index=FALSE, it will be calculated from the CTSR.VI using \code{\link{AnMax.VI}}.
#' @param rf.b4
#' If a breakpoint in the VPR is detected this is the optimial accumulated rainfall before the breakpoint. must be the same length as the anu.VI
#' @param rf.af
#' If a breakpoint in the VPR is detected this is the optimial accumulated rainfall after the breakpoint. must be the same length as the anu.VI
#' @param sig
#' Significance of all the functions, sig=0.05
#' @param print
#' Prints more details at every step of the procces
#' @param plot
#' creates a plots of every step
#' @param details
#' returns adational details for VPR and Residual trends
#'
#' from AnMax
#' @param CTSR.VI
#' Complete Time Series of Vegetation Index. An object of class \code{'ts'}. Monthly time series of VI values
#' @return anu.VI
#' The annual (Growing season) max VI.
#' @return Max.Month
#' a Vector that has the month number where max values wew observed
#' @return VI.index
#' the index of the CTSR.VI ts that the anu.VI values occur at.R indexs from 1 rather than 0
#'
#'
#'
#'
#' @title Antecedental Accumulation calculator for the VI Complete Time Series
#'
#' @description
#' Takes the Complete Time Series Vegetation Index and a table of every possible accumulation period and offset
#' period.  A OLS is calculated \code{\link{lm}} for every combination of VI ~ rainfall.  This Function preferences those results where
#' slope>0 (increase in rainfall causes an increase in vegetation), returning the rainfall accumulation that has the highest R-squared
#' and a positive slope. If no combinations produce a positive slope then the one with the highest Rsquared is returned.
#'
#' @param CTSR.VI
#' Complete Time Series of Vegetation Index. An object of class \code{'ts'}. Monthly time series of VI values
#' @param ACP.table
#' A table of every combination of offset period and accumulation period. if ACP.table = FALSE, CTSR.RF and acu.RF must be provided
#' @return summary
#' (To be filled in)
#' @return CTSR.precip (CTSR.RF)
#' Complete Time Series of the optimally accumulated rainfall.An object of class 'ts'.
#' @export


#' @title Antecedental Accumulation calculator for the Annual Max VI Time Series
#'
#' @description
#' Takes the Annual Max VI Time Series, the VI.index and a table of every possible accumulation period and offset
#' period.  A OLS is calculated \code{\link{lm}} for every combination of VI ~ rainfall.  This Function preferences those results where
#' slope>0 (increase in rainfall causes an increase in vegetation), returning the rainfall accumulation that has the highest R-squared
#' and a positive slope. If no combinations produce a positive slope then the one with the highest Rsquared is returned.MISSING non
#' peramtric and other variables
#'
#' @param anu.VI
#' The annual (Growing season) max VI. Which can be calculated from the CTSR.VI using \code{\link{AnMax.VI}}.
#' @param VI.index
#' the index of the CTSR.VI ts that the anu.VI values occur at. NOTE. R indexs from 1 rather than 0.
#' NOTE. if VI.index=FALSE, It can be calculated from the CTSR.VI using \code{\link{AnMax.VI}}.
#' @param ACP.table
#' A table of every combination of offset period and accumulation period. It can be calculated from ________
#'
#' @return _____(acu.RF)
#' The optimal accumulated rainfall for anu.VI. Mut be a object of class 'ts' and of equal length to anu.VI.
#' if anu.RF=FALSE, it will be calculated from ACP.table. see(____)
#' @return summary
#' (To be filled in)
#' @export
#'
#'
#'#' @title BFAST Breakpoint Detector
#'
#' @description
#' This function will perform the complete Time Series Segmented Residual Trend (TSS.RESTREND) methodology.Takes in a complete monthly time series of a VI and its corrosponding precipitation. Will caculate missing input varibles,
#' look for breakpoints using the BFAST function. The significance of the breakpoin in the residuals and the VPR is assessed using a chow test an
#' then the total time series change is calculated. Calls ____________
#' @importFrom bfast bfast


#' @param CTSR.VI
#' Complete Time Series of Vegetation Index. An object of class \code{'ts'}. Monthly time series of VI values
#' @param CTSR.RF
#' Complete Time Series of Rain Fall.  An object of class 'ts' and be the same length and cover the same time range as CTSR.VI.
#' If ACU.table is provided, CTSR.RF will be automitaclly calculated by \code{\link{rainfall.accumulator}}
#' @param sig
#' Significance of all the functions, sig=0.05
#' @param print
#' Prints more details at every step of the procces
#' @param plot
#' creates a plots of every step
#' @param details
#' returns adational details for BFAST
#'
#' @return breakpoints (____)
#' (To be filled in)
#' @export



