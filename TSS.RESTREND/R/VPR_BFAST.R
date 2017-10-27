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
#' @param BFAST.raw
#'        Defualt = FALSE
#'        If TRUE will perform a BFAST (season="harmonic") on the CTSR.VI
#'        If FALSE will perform BFAST on the CTSR VPR residuals
#' @return List of objects:
#' @return \bold{bkps}
#'          The index of the Breakpoints detected. If no breakpoints are detected, bkps = FASLE
#' @return \bold{BFAST.obj}
#'          See \code{\link[bfast]{bfast}}
#' @return \bold{CTS.lm}
#'          the \code{\link[stats]{lm}} of CTSR.VI and CTSR.RF
#' @return \bold{BFAST.type}
#'          the type of BFAST done (VPR residuals or on the VI timeseris itself)
#' @export
#'
#' @examples
#' VPRBFdem <- VPR.BFAST(segVPRCTSR$cts.NDVI, segVPRCTSR$cts.precip)
#' print(VPRBFdem)

VPR.BFAST <- function(CTSR.VI, CTSR.RF, season="none", BFAST.raw=FALSE) {
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
  if (BFAST.raw){
    bf.fit <- bfast(CTSR.VI, h=0.15, season="harmonic", max.iter=3, level = 0.05)
    bft <- "raw.VI"
  }else{
    bf.fit <- bfast(resid.ts, h=0.15, season=season, max.iter=3, level = 0.05)
    bft <- "CTSR.VPR"
  }


  if (bf.fit$nobp$Vt[[1]] == FALSE) {
    # Get the number of breakpoints
    numiter <- length(bf.fit$output)
    tmp <- bf.fit$output[[numiter]]$bp.Vt[1]$breakpoints
    # stack and return the data
    return(structure(list(bkps = tmp, BFAST.obj=bf.fit, CTS.lm = CTS.fit, BFAST.type=bft), class = "BFAST.Object"))
  } else {
    # Stack the data
    return(structure(list(bkps = FALSE, BFAST.obj=bf.fit, CTS.lm = CTS.fit, BFAST.type=bft), class = "BFAST.Object"))
    }
}

