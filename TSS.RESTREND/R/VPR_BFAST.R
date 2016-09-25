#' @title BFAST Breakpoint Detector
#'
#' @description
#'      Uses the Complete VI and Rainfall time series, calculates a \code{\link[stats]{lm}} between them
#'      And then performs a \code{\link[bfast]{bfast}}
#'
#' @author Arden Burrell, arden.burrell@unsw.edu.au
#'
#' @importFrom bfast bfast
#'
#' @inheritParams TSSRESTREND
#'
#' @return List of objects:
#' @return \bold{bkps}
#'          The index of the Breakpoints detected. If no breakpoints are detected, bkps = FASLE
#' @return \bold{BFAST.obj}
#'          See \code{\link[bfast]{bfast}}
#' @return \bold{CTS.lm}
#'          the \code{\link[stats]{lm}} of CTSR.VI and CTSR.RF
#' @export
#'
#' @examples
#' VPRBFdem <- VPR.BFAST(segVPRCTSR$cts.NDVI, segVPRCTSR$cts.precip)
#' print(VPRBFdem)

VPR.BFAST <- function(CTSR.VI, CTSR.RF, season="none") {
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
  #perform the BFAST
  bf.fit <- bfast(resid.ts, h=0.15, season=season, max.iter=3, level = 0.05)

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

