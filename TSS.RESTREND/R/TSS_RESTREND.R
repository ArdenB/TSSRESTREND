#' @title Time Series Segmentation of Residual Trends (MAIN FUNCTION)
#'
#' @importFrom stats coef end frequency lm sd start time ts
#' @importFrom graphics abline arrows legend par plot
#' @importFrom utils tail read.csv
#' @importFrom broom glance
#'
#' @description
#' Time Series Segmented Residual Trend (TSS.RESTREND) methodology.Takes in a complete monthly
#' time series of a VI and its corrosponding precipitation. Will caculate missing input varibles,
#' look for breakpoints using the BFAST function. The significance of the breakpoin in the residuals
#' and the VPR is assessed using a chow test an then the total time series change is calculated.
#'
#' @author Arden Burrell, arden.burrell@unsw.edu.au
#'
#' @param CTSR.VI
#'        Complete Monthly Time Series of Vegetation Index values.
#'        An object of class \code{'ts'} object without NA's.
#' @param ACP.table
#'        A table of every combination of offset period and accumulation period.for precipitation
#'        ACP.table can be calculated using the \code{\link{climate.accumulator}}.
#' @note  if ACP.table = FALSE, CTSR.RF and acu.RF must be provided as well as
#'        rf.b4 and rf.af for \code{'ts'} with a breakpoint in the VPR.
#' @param ACT.table
#'        A table of every combination of offset period and accumulation period.for temperature
#'        ACP.table can be calculated using the \code{\link{climate.accumulator}}.
#' @param CTSR.RF
#'        Complete Time Series of Rainfall. An object of class 'ts' object without NA's
#'        and be the same length and cover the same time range as CTSR.VI.
#'        If ACP.table is provided, CTSR.RF will be automitaclly calculated using the
#'        \code{\link{ACP.calculator}}
#' @param CTSR.TM
#'        Complete Time Series of temperature. An object of class 'ts' object without NA's
#'        and be the same length and cover the same time range as CTSR.VI.  Default (CTSR.TM=NULL).
#'        If ACT.table is provided, CTSR.RF will be automitaclly calculated using the
#'        \code{\link{ACP.calculator}}
#' @param anu.VI
#'        The annual (Growing season) max VI. Must be a object of class \code{'ts'} without NA's.
#'        if anu.VI=FALSE, it will be calculated from the CTSR.VI using \code{\link{AnMaxVI}}.
#' @param acu.RF
#'        The optimal accumulated rainfall for anu.VI. Must be a object of class \code{'ts'} without
#'        NA's and be of equal length and temporal range to anu.VI. if anu.RF=FALSE, it will be
#'        calculated from ACP.table usingthe \code{\link{AnnualRF.Cal}}
#' @param acu.TM
#'        The optimal accumulated rainfall for anu.TM. Must be a object of class \code{'ts'} without
#'        NA's and be of equal length and temporal range to anu.TM. if anu.TM=FALSE, it will be
#'        calculated from ACT.table usingthe \code{\link{AnnualRF.Cal}}
#' @param VI.index
#'        the index of the CTSR.VI ts that the anu.VI values occur at. Must be the same length
#'        as anu.VI. NOTE. R indexs from 1 rather than 0.
#'        if VI.index=FALSE, it will be calculated from the CTSR.VI using \code{\link{AnMaxVI}}.
#' @param rf.b4
#'        If a breakpoint in the VPR is detected this is the optimial accumulated rainfall before
#'        the breakpoint. must be the same length as the anu.VI. If ACP.table is provided it will
#'        be generated using \code{\link{AnnualRF.Cal}}
#' @param rf.af
#'        If a breakpoint in the VPR is detected this is the optimial accumulated rainfall after
#'        the breakpoint. must be the same length as the anu.VI. If ACP.table is provided it will
#'        be generated using \code{\link{AnnualRF.Cal}}
#' @param sig
#'        Significance of all the functions. defualt sig=0.05
#' @param season
#'        See \code{\link[bfast]{bfast}}. This season value only applies to bfast done using the CTS
#'        VPR. if a non VPR adjusted BFAST is performed.a harmonic season is used.
#' @param exclude
#'        A numberic vector containg months excluded from breakpoint detection.  This was included to
#'        allow sensor transitions to be masked.
#' @param allow.negative
#'        IF true, will not preference positive slope in either CTSR or VI calculations
#' @return
#' An object of class \code{'TSSRESTREND'} is a list with the following elements:
#' \describe{
#'   \item{\bold{\emph{summary}}}{
#'    \describe{
#'    \item{Method}{The method used to determine total change. (\emph{RESTREND} see \code{\link{RESTREND}},
#'      \emph{segmented.RESTREND} see \code{\link{seg.RESTREND}}, \emph{segmented.VPR} see
#'      \code{\link{seg.VPR}})}
#'    \item{Total.Change}{The total significant change. Residual.Change + VPR.HeightChange. }
#'    \item{Residual.Change}{The change in the VPR Residuals over the time period}
#'    \item{VPR.HeightChange}{The change in VI at mean rainfall for a "ts" with a significant
#'      breakpoint in the VPR}
#'    \item{model.p}{p value of the regression model fitted to the VPR. See \code{\link[stats]{lm}}}
#'    \item{residual.p}{p value of the regression model fitted to the VPR Residuals. See \code{\link[stats]{lm}}}
#'    \item{VPRbreak.p}{the p value associated with the break height. See \code{\link[stats]{lm}}}
#'    \item{bp.year}{The Year of the most significant breakpoint}
#'    }}
#'   \item{\bold{\emph{ts.data}}}{The Time series used in analysis. See Arguments for description
#'    \itemize{
#'      \item CTSR.VI
#'      \item CTSR.RF
#'      \item anu.VI
#'      \item VI.index
#'      \item acu.RF
#'      \item StdVar.RF see \code{\link{seg.VPR}})}
#'      }
#'   \item{\bold{\emph{ols.summary}}}{
#'      \describe{
#'        \item{chow.summary}{summary of the most significant breakpoint. }
#'        \item{chow.ind}{Summary of every detected breakpoint}
#'        \item{OLS.table}{A matrix containing the coefficents for the CTS.fit, VPR.fit, RESTREND.fit and segVPR.fit}
#'        }}
#'    \item{\bold{\emph{TSSRmodels}}}{
#'    models of class "lm" \code{\link[stats]{lm}} and class "bfast" \code{\link[bfast]{bfast}} generated.}
#'    }
#'
#' @seealso
#'  \itemize{
#'    \item \code{\link{plot.TSSRESTREND}}
#'    \item \code{\link{print.TSSRESTREND}}
#'    }
#' @export
#' @examples
#' \dontrun{
#' #To get the latest version of the package (Still in development)
#' install.packages("devtools")
#' library("devtools")
#' install_github("ArdenB/TSSRESTREND", subdir="TSS.RESTREND")
#' library(TSS.RESTREND)
#' }
#' #Find the path of the rabbitRF.csv dataset, read it in and turn it into a time series
#' rf.path<- system.file("extdata", "rabbitRF.csv", package = "TSS.RESTREND", mustWork = TRUE)
#' in.RF <- read.csv(rf.path)
#' rf.data <- ts(in.RF, end=c(2013,12), frequency = 12)
#'
#' #Find the path of the rabbitVI.csv dataset and read it in
#' vi.path <- system.file("extdata", "rabbitVI.csv", package = "TSS.RESTREND", mustWork = TRUE)
#' in.VI <- read.csv(vi.path)
#' CTSR.VI <- ts(in.VI, start=c(1982, 1), end=c(2013,12), frequency = 12)
#'
#' #Define the max accumuulation period
#' max.acp <- 12
#' #Define the max offset period
#' max.osp <- 4
#' #Create a table of every possible precipitation value given the max.acp and max.osp
#' ACP.table <- rainfall.accumulator(CTSR.VI, rf.data, max.acp, max.osp)
#' results <- TSSRESTREND(CTSR.VI, ACP.table)
#' print(results)
#'
#' \dontrun{
#' plot(results, verbose=TRUE)
#' }
#'
TSSRESTREND <- function(CTSR.VI, ACP.table=FALSE, ACT.table=NULL, CTSR.RF=FALSE, CTSR.TM=NULL,  anu.VI=FALSE, acu.RF=FALSE, acu.TM=NULL, VI.index=FALSE,
                         rf.b4=FALSE, rf.af=FALSE, sig=0.05, season="none", exclude=0, allow.negative=FALSE){

  while (TRUE){ #Test the variables for consistenty
    #Each check liiks at a different paramter. If the data fails the check will stop, else, it breaks after all the checks
    if (class(CTSR.VI) != "ts")
      stop("CTSR.VI Not a time series object. Please check the data")
    if ((class(ACP.table)=="logical") && (!CTSR.RF || acu.RF))
      stop("Insufficent Rainfall data provided. Provide either a complete ACP.table or both the CTSR.RF & acu.RF")
    if ((!anu.VI)||(!VI.index)){
      # Get the annual Max VI values
      max.df <- AnMaxVI(CTSR.VI)
      # Pull the key components from the result
      anu.VI <- max.df$Max #the VI values
      VI.index <- max.df$index #the indes values
      Max.Month <- max.df$Max.Month #month if the year the even occured
    }else{
      if (class(anu.VI) != "ts")
        stop("anu.VI Not a time series object")
    }
    # Change the allow.negative for multivariate regression
    if (!is.null(ACT.table)){allow.negative=TRUE}

    if (!CTSR.RF){
      #Calculate the Complete time seties Accumulation (rainfall)
      CTS.Str <- ACP.calculator(CTSR.VI, ACP.table, ACT.table, allow.negative=allow.negative)
      # If the allow negative is on, this is ignored, else perform a negative slope check and perform a significance check
      #This mod is will impact results comparisons before V0.1.04
      if ((!allow.negative && as.numeric(CTS.Str$summary)[1] <0)||as.numeric(CTS.Str$summary)[4] >sig){
        BFraw = TRUE
      }else{BFraw=FALSE}
      #Pull out the relevant paramters for use
      CTSR.RF <- CTS.Str$CTSR.precip #RF values
      if (is.null(CTSR.TM)){
        CTSR.TM <- CTS.Str$CTSR.tmp #Temperature values, null if temp is not considered
      }
      details.CTS.VPR <- CTS.Str$summary #Summay of the FIT between the CTS.VRP
      CTSR.osp <- CTS.Str$CTSR.osp #CTS Off set period
      CTSR.acp <- CTS.Str$CTSR.acp #CTS Sccumulation period
      CTSR.tosp <- CTS.Str$CTSR.tosp #CTS Off set period
      CTSR.tacp <- CTS.Str$CTSR.tacp #CTS Sccumulation period

    }else{
      if (class(CTSR.RF) != "ts")
        stop("CTSR.RF Not a time series object")
      #get the time data out
      start.ti <- time(CTSR.VI)
      freq <- frequency(CTSR.VI)
      #check the two ts object cover the same time period
      start.ti2 <- time(CTSR.RF)
      freq2 <- frequency(CTSR.RF)
      #Check the start dates and the frequency are correct
      if (!identical(start.ti, start.ti2))
        stop("ts objects do not have the same time, (CTSR.VI & CTSR.RF)")
      if (!identical(freq, freq2))
        stop("ts objects do not have the same frequency, (CTSR.VI & CTSR.RF)")
    }
    if (!acu.RF){
      # Calculate the Annual Accumulated Rainfall
      precip.df <- AnnualRF.Cal(anu.VI, VI.index, ACP.table, ACT.table, allow.negative=allow.negative)
      # Pull out and store key values
      osp <- precip.df$osp # offset period
      acp <- precip.df$acp # Accumulation period
      tosp <- precip.df$tosp # offset period
      tacp <- precip.df$tacp # Accumulation period
      acu.RF <- precip.df$annual.precip # precip values
      acu.TM <- precip.df$annual.temp # precip values
      details.VPR<- precip.df$summary # The summary of the lm between rainfall and Vegetation
    }else{
      # Check the passed accumulated rainfall
      if (class(acu.RF) != "ts")
        stop("acu.RF Not a time series object")
      # get the time peramaters
      st.ti <- time(anu.VI)
      st.f <- frequency(anu.VI)
      st.ti2 <- time(acu.RF)
      st.f2 <- frequency(acu.RF)
      #check the two ts object cover the same time period and frequency
      if (!identical(st.ti, st.ti2))
        stop("ts object do not have the same time, (acu.RF & anu.VI)")
      if (!identical(st.f, st.f2))
        stop("ts object do not have the same frequency, (acu.RF & anu.VI)")
    }
    # Check passed breakpoint rainfall results
    if (class(rf.b4) != "logical"){
      if (length(rf.b4) != (length(rf.af)))
          stop("rf.b4 and rf.af are different shapes. They must be the same size and be th same lenths as acu.VI")
    }
    break
  }
  # Pass the infomation about the VI and RF as well as the nfast method to the VPR.BFAST script
  bkp = VPR.BFAST(CTSR.VI, CTSR.RF, CTSR.TM, season=season, BFAST.raw = BFraw)

  bp <- bkp$bkps
  BFAST.obj <- bkp$BFAST.obj #For the models Bin
  CTS.lm <- bkp$CTS.lm #For the Models Bin
  bp <- bp[!bp %in% exclude] #remove breakpoints in the exclude list (Sensor transitions)
  BFT <-  bkp$BFAST.type #Type of BFAST used
  acum.df <- data.frame(CTSR.osp=CTSR.osp, CTSR.acp=CTSR.acp, CTSR.tosp=CTSR.tosp, CTSR.tacp=CTSR.tacp,
                        osp=osp, acp=acp, tosp=tosp, tacp=tacp,  osp.b4=NaN, acp.b4=NaN, tosp.b4=NaN, tacp.b4=NaN,
                        osp.af=NaN, acp.af = NaN, tosp.af=NaN, tacp.af = NaN)

  if (class(bp)=="logical"|length(bp)==0) { #Should catch both the false and the no breakpoints
    bp <- FALSE
    # no breakpoints detected by the BFAST
    test.Method = "RESTREND"
    chow.sum <- data.frame(abs.index=FALSE, yr.index = FALSE, reg.sig=FALSE, VPR.bpsig = FALSE)
    chow.bpi <- FALSE
  }else{
    bp<-as.numeric(bkp$bkps)
    res.chow <- CHOW(anu.VI, acu.RF, acu.TM, VI.index, bp, sig=sig)
    brkp <- as.integer(res.chow$bp.summary["yr.index"]) #this isn't right
    chow.sum <-res.chow$bp.summary
    chow.bpi <- res.chow$allbp.index
    test.Method = res.chow$n.Method
  }

  if (test.Method == "RESTREND"){
    # browser()
    result <- RESTREND(anu.VI, acu.RF, acu.TM,  VI.index, sig=sig)
    result$TSSRmodels$CTS.fit <- CTS.lm
    result$TSSRmodels$BFAST <- BFAST.obj
    result$ts.data$CTSR.VI <- CTSR.VI
    result$ts.data$CTSR.RF <- CTSR.RF
    result$ols.summary$chow.sum <- chow.sum
    result$ols.summary$chow.ind <- chow.bpi
    result$ols.summary$OLS.table["CTS.fit",] <- details.CTS.VPR

  }else if (test.Method == "seg.RESTREND"){
    breakpoint = as.integer(res.chow$bp.summary[2])
    result <- seg.RESTREND(anu.VI, acu.RF, acu.TM, VI.index, brkp,  sig=sig)
    result$TSSRmodels$CTS.fit <- CTS.lm
    result$TSSRmodels$BFAST <- BFAST.obj
    result$ts.data$CTSR.VI <- CTSR.VI
    result$ts.data$CTSR.RF <- CTSR.RF
    result$ols.summary$chow.sum <- chow.sum
    result$ols.summary$chow.ind <- chow.bpi
    result$ols.summary$OLS.table["CTS.fit",] <- details.CTS.VPR

  }else if (test.Method == "seg.VPR"){
    if ((!rf.b4)||(!rf.af)){
      VPRbp.df <-AnnualRF.Cal(anu.VI, VI.index, ACP.table, ACT.table, Breakpoint = brkp, allow.negative = allow.negative)
      rf.b4 <- VPRbp.df$rf.b4
      rf.af <- VPRbp.df$rf.af
      tm.b4 <- VPRbp.df$tm.b4
      tm.af <- VPRbp.df$tm.af
      acum.df$osp.b4 <- VPRbp.df$osp.b4
      acum.df$acp.b4 <- VPRbp.df$acp.b4
      acum.df$tosp.b4 <- VPRbp.df$tosp.b4
      acum.df$tacp.b4 <- VPRbp.df$tacp.b4
      acum.df$osp.af <- VPRbp.df$osp.af
      acum.df$acp.af <- VPRbp.df$acp.af
      acum.df$tosp.af <- VPRbp.df$tosp.af
      acum.df$tacp.af <- VPRbp.df$tacp.af
    }
    breakpoint = as.integer(res.chow$bp.summary[2])
    print(brkp)
    result <- seg.VPR(anu.VI, acu.RF, acu.TM, VI.index, brkp, rf.b4, rf.af, tm.b4, tm.af, sig=sig)
    result$TSSRmodels$CTS.fit <- CTS.lm
    result$TSSRmodels$BFAST <- BFAST.obj
    result$ts.data$CTSR.VI <- CTSR.VI
    result$ts.data$CTSR.RF <- CTSR.RF
    result$ols.summary$chow.sum <- chow.sum
    result$ols.summary$chow.ind <- chow.bpi
    result$ols.summary$OLS.table["CTS.fit",] <- details.CTS.VPR
  }
  result$acum.df <- acum.df
  #add the bfast method to the results summary
  result$summary$BFAST.Method <- BFT
  return(result)
}

