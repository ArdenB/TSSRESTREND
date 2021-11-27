#' @title BFAST Breakpoint Detector
#'
#' @description
#'      takes the Complete VI and optimally accumulated Rainfall (and tmperature if included), calculates a \code{\link[stats]{lm}} between them
#'      And then performs a \code{\link[bfast]{bfast}}.in the residuals.  If BFAST.raw=TRUE, it will perform bfast on the Complete VI ts
#'
#' @author Arden Burrell, arden.burrell@unsw.edu.au
#'
#' @importfrom bfast bfast
#'
#' @inheritParams TSSRESTREND
#' @param CTSR.TM
#'        Complete Time Series of temperature. An object of class 'ts' object without NA's
#'        and be the same length and cover the same time range as CTSR.VI.  Default (CTSR.TM=NULL).
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
#'          the type of BFAST done (VPR residuals or on the VI timeseries itself)
#' @export
#'
#' @examples
#' \dontrun{
#' VPRBFdem <- VPR.BFAST(segVPRCTSR$cts.NDVI, segVPRCTSR$cts.precip)
#' print(VPRBFdem)}

VPR.BFAST <- function(CTSR.VI, CTSR.RF, CTSR.TM=NULL, season="none", BFAST.raw=FALSE, h = 0.15) {

  # ========== Complete time series breakpoint detection using BFAST ==========
  # ===========================================================================
  #functions takes the complete time series VI and rainfall (RF)
  # ===== Basic data sanity checks =====
  #Check the objects are Time series
  if (class(CTSR.VI) != "ts")
    stop("CTSR.VI Not a time series object")
  if (class(CTSR.RF) != "ts")
    stop("CTSR.RF Not a time series object")
  # get the time data out
  ti <- time(CTSR.VI)
  f <- frequency(CTSR.VI)
  # check the two ts object cover the same time period
  ti2 <- time(CTSR.RF)
  f2 <- frequency(CTSR.RF)
  if (!identical(ti, ti2))
    stop("ts object do not have the same time")
  if (!identical(f, f2))
    stop("ts object do not have the same frequency")

  # ===== Get the regression fot the VCR/VPR =====
  if (is.null(CTSR.TM)) { # no temperature data
    CTS.fit <- lm(CTSR.VI ~ CTSR.RF)
  }else{
    CTS.fit <- lm(CTSR.VI ~ CTSR.RF + CTSR.TM)
  }
  # Convert to a ts object
  resid.ts <- ts(CTS.fit$residuals, start = ti[1], end = tail(ti, 1), frequency = f)

  # ===== perform the BFAST =====
  if (BFAST.raw){
    bf.fit <- bfast(CTSR.VI, h = h, season = "harmonic", max.iter = 3, level = 0.05)
    bft <- "raw.VI"
  } else {
    bf.fit <- bfast(resid.ts, h = h, season = season, max.iter = 3, level = 0.05)
    bft <- "CTSR.VPR"
    #browser()
  }

  # ===== Create and object to return =====
  if (bf.fit$nobp$Vt[[1]] == FALSE) {
    # Get the number of breakpoints
    numiter <- length(bf.fit$output)
    tmp <- bf.fit$output[[numiter]]$bp.Vt[1]$breakpoints
    # stack and return the data
    return(structure(list(bkps = tmp, BFAST.obj = bf.fit, CTS.lm = CTS.fit, BFAST.type = bft), class = "BFAST.Object"))
  } else {
    # Stack the data
    return(structure(list(bkps = FALSE, BFAST.obj = bf.fit, CTS.lm = CTS.fit, BFAST.type = bft), class = "BFAST.Object"))
    }
}

