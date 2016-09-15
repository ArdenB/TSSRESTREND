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
  #function.  rf.b4=FALSE, rf.af=FALSE, will be removed as soon as
#' @title Time Series Segmentation of Residual Trends (MAIN FUNCTION)
#'
#' @description
#' This function will perform the complete Time Series Segmented Residual Trend (TSS.RESTREND) methodology.Takes in a complete monthly time series of a VI and its corrosponding precipitation. Will caculate missing input varibles,
#' look for breakpoints using the BFAST function. The significance of the breakpoin in the residuals and the VPR is assessed using a chow test an
#' then the total time series change is calculated. Calls ____________
#' @importFrom stats coef end frequency lm sd start time ts
#' @importFrom graphics abline arrows legend par plot
#' @importFrom utils tail
#' @importFrom broom glance
#' @author Arden Burrell, arden.burrell@unsw.edu.au
#'
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
#' @return list
#' (To be filled in)
#' @export
#' @examples
#'
#'
#' \dontrun{
#' print("Hello World")
#' }
#'
TSS.RESTREND <- function(CTSR.VI, ACP.table=FALSE, CTSR.RF=FALSE, anu.VI=FALSE, acu.RF=FALSE, VI.index=FALSE,
                         rf.b4=FALSE, rf.af=FALSE, sig=0.05, print=FALSE, plot=TRUE, details=FALSE){

  while (TRUE){ #Test the variables
    if (class(CTSR.VI) != "ts")
      stop("CTSR.VI Not a time series object. Please check the data")

    if ((!ACP.table) && (!CTSR.RF || acu.RF))
      stop("Rainfall data invalid. ACP.table or (CTSR.RF & acu.RF")

    if ((!anu.VI)||(!VI.index)){
      max.df <- AnMax.VI(CTSR.VI)
      anu.VI <- max.df$Max
      VI.index <- max.df$index
    }else{
      if (class(anu.VI) != "ts")
        stop("anu.VI Not a time series object")
    }
    if (!CTSR.RF){
      CTS.Str <- ACP.calculator(CTSR.VI, ACP.table)
      CTSR.RF <- CTS.Str$CTSR.precip
      details.CTS.VPR <- CTS.Str$summary
    }else{
      if (class(CTSR.RF) != "ts")
        stop("CTSR.RF Not a time series object")
      #get the time data out
      start.ti <- time(CTSR.VI)
      freq <- frequency(CTSR.VI)
      #check the two ts object cover the same time period
      start.ti2 <- time(CTSR.RF)
      freq2 <- frequency(CTSR.RF)
      if (!identical(start.ti, start.ti2))
        stop("ts objects do not have the same time, (CTSR.VI & CTSR.RF)")
      if (!identical(freq, freq2))
        stop("ts objects do not have the same frequency, (CTSR.VI & CTSR.RF)")
    }
    # browser()
    if (!acu.RF){
      precip.df <- AnnualRF.Cal(anu.VI, VI.index, ACP.table)
      acu.RF <- precip.df$annual.precip
      details.acu.RF <- precip.df$summary
    }else{
      if (class(acu.RF) != "ts")
        stop("acu.VI Not a time series object")
      st.ti <- time(anu.VI)
      st.f <- frequency(anu.VI)
      #check the two ts object cover the same time period
      st.ti2 <- time(acu.RF)
      st.f2 <- frequency(acu.RF)
      if (!identical(st.ti, st.ti2))
        stop("ts object do not have the same time, (acu.RF & anu.VI)")
      if (!identical(st.f, st.f2))
        stop("ts object do not have the same frequency, (acu.RF & anu.VI)")
    }
    if (class(rf.b4) != "logical"){
      if (length(rf.b4) != (length(rf.af)))
          stop("rf.b4 and rf.af are different shapes. They must be the same size and be th same lenths as acu.VI")
    }
    break
  }

  bkp = VPR.BFAST(CTSR.VI, CTSR.RF, print=print, plot=plot, details = details)
  bp <- bkp$bkps
  if (!bp){# no breakpoints detected by the BFAST
    test.Method = "RESTREND"
    if (plot){
      VImax.plot(anu.VI)
    }
  }else{
    bp<-as.numeric(bkp$bkps)
    res.chow <- CHOW(anu.VI, acu.RF, VI.index, bp, sig=sig, print=print)
    brkp <- as.integer(res.chow$bp.summary["yr.index"]) #this isn't right
    test.Method = res.chow$n.Method
    if (plot){
      VImax.plot(anu.VI, brkp=brkp)
    }
  }
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #add NDVI plot
  # browser()

  if (test.Method == "RESTREND"){
    result <- RESTREND(anu.VI, acu.RF, VI.index, sig=sig, print=print, plot=plot)
  }else if (test.Method == "seg.RESTREND"){
    breakpoint = as.integer(res.chow$bp.summary[2])
    print(brkp)
    result <- seg.RESTREND(anu.VI, acu.RF, VI.index, brkp,  sig=sig, print=print, plot=plot)
  }else if (test.Method == "seg.VPR"){
    if ((!rf.b4)||(!rf.af)){
      #Improve the figure
      VPRbp.df <-AnnualRF.Cal(anu.VI, VI.index, ACP.table, Breakpoint = brkp)
      rf.b4 <- VPRbp.df$rf.b4
      rf.af <- VPRbp.df$rf.af
      #!!!!!!!!!!!!!!!!!!!!!!!!!!#
    }
    breakpoint = as.integer(res.chow$bp.summary[2])
    print(brkp)
    result <- seg.VPR(anu.VI, acu.RF, VI.index, brkp, rf.b4, rf.af, sig=sig, print=print, plot=plot)
  }
  # print(result$summary)
  # browser()
  return(result) #add CTSR
}

