#' @title BFAST Breakpoint Detector
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

VPR.BFAST <- function(CTSR.VI, CTSR.RF, season="none") { #plot=TRUE, details=FALSE
  #functions takes the complete time series VI and rainfall (RF)

  #Check the objects are Time series
  if (class(CTSR.VI) != "ts")
    stop("CTSR.VI Not a time series object")
  if (class(CTSR.RF) != "ts")
    stop("CTSR.RF Not a time series object")
  #get the time data out
  ti <- time(CTSR.VI)
  f <- frequency(CTSR.VI)
  #check the two ts object cover the same time period
  ti2 <- time(CTSR.RF)
  f2 <- frequency(CTSR.RF)
  if (!identical(ti, ti2))
    stop("ts object do not have the same time")
  if (!identical(f, f2))
    stop("ts object do not have the same frequency")

  # Fit the two lines
  CTS.fit <- lm(CTSR.VI ~ CTSR.RF)
  #Convert to a ts object
  resid.ts<- ts(CTS.fit$residuals, start=ti[1], end=tail(ti, 1), frequency = f)

  #To Be Removed
  # #Print and plot
  # if (print) { #Need better explination and a title but it will do for the moment
  #   print(summary(CTS.fit))
  # }
  # if (plot){
  #   #This plot need serius work but will do for the moment
  #   plot(resid.ts)
  # }

  #perform the BFAST
  bf.fit <- bfast(resid.ts, h=0.15, season=season, max.iter=3, level = 0.05)

  # if (plot){
  #   #if plot is requested
  #   plot(bf.fit)
  # }

  if (bf.fit$nobp$Vt[[1]] == FALSE) {
    numiter <- length(bf.fit$output)
    tmp <- bf.fit$output[[numiter]]$bp.Vt[1]$breakpoints
    return(structure(list(bkps = tmp, BFAST.obj=bf.fit, CTS.lm = CTS.fit), class = "BFAST.Object"))
    # if (details){
    #   return(structure(list(bkps = tmp, BFAST.obj=bf.fit, CTS.lm = CTS.fit), class = "BFAST.Object"))
    # }else{return(structure(list(bkps = tmp, BFAST.obj=FALSE), class = "BFAST.Object"))
    }else {
    return(structure(list(bkps = FALSE, BFAST.obj=bf.fit, CTS.lm = CTS.fit), class = "BFAST.Object"))
    #
    # if (details){
    #   return(structure(list(bkps = FALSE, BFAST.obj=bf.fit), class = "BFAST.Object"))
    # }else{return(structure(list(bkps = FALSE, BFAST.obj=FALSE), class = "BFAST.Object"))}
  }
}

