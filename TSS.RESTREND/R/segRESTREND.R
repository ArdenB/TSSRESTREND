#' @title Segmented RESTREND (RESidual TREND)
#'
#' @description
#' For ts with a significant breakpoints in the residuals but not in the VPR. Takes annual VI max and it associated
#' optimal accumulated precipitation and regresses it with a dummy variable that is 0 before the breakpoint and 1
#' after it
#'
#' @importFrom strucchange sctest
#' @importFrom broom glance
#' @author Arden Burrell, arden.burrell@unsw.edu.au
#'
#' @inheritParams TSSRESTREND
#' @inheritParams seg.VPR
#'
#' @return a list of class TSSRESTREND.
#'        See \code{\link{TSSRESTREND}} for details. Note. if called seperatly from TSSRESTREND,
#'        this list will be incomplete.
#' @export
#'
#' @examples
#' # brkp can be determined using VPR.BFAST and CHOW.
#' brkp <-  as.integer(11)
#' resu <- seg.RESTREND(segRESTREND$max.NDVI, segRESTREND$acc.precip, segRESTREND$index, brkp)

seg.RESTREND <- function(anu.VI, acu.RF, VI.index, breakpoint, acu.TM = NULL, sig=0.05){
  # ==============================================================================================
  # ========== Sanity check the input data ==========
  while (TRUE) {
    if (class(anu.VI) != "ts")
      stop("anu.VI Not a time series object")
    if (class(acu.RF) != "ts")
      stop("acu.VI Not a time series object")
    ti <- time(anu.VI)
    f <- frequency(anu.VI)
    #check the two ts object cover the same time period
    ti2 <- time(acu.RF)
    f2 <- frequency(acu.RF)
    if (!identical(ti, ti2))
      stop("ts object do not have the same time range")
    if (!identical(f, f2))
      stop("ts object do not have the same frequency")
    if (class(breakpoint) != "integer")
      stop("Breakpoint must be an interger")
    if (length(breakpoint) != 1)
      stop("Breakpoint must be an interger of length 1")
    break
  }
  # ==============================================================================================
  # ========== Calculate the VPR/VCR ==========
  if (is.null(acu.TM)) {# No Temp
    VPR.fit <- lm(anu.VI ~ acu.RF)
  } else {# temp
    VPR.fit <- lm(anu.VI ~ acu.RF+acu.TM)
  }
  # ===== Setup empty matrix to hold paramaters or the linear models =====
  m <- matrix(nrow = (4), ncol = 8)
  m[] <- NaN
  rownames(m) <- c("CTS.fit", "VPR.fit", "RESTREND.fit", "segVPR.fit")
  colnames(m) <- c(
    "slope", "temp.coef", "intercept", "p.value", "R^2.Value",
    "Break.Height", "Slope.Change", "Slope.ChangeTmp"
    )
  # ========== Pull out the key values ==========
  R.pval <- glance(VPR.fit)$p.value
  R.Rval <- summary(VPR.fit)$r.square
  if (is.null(acu.TM)) {
    R.tcoef <- NaN
  } else {
    R.tcoef <- as.numeric(coef(VPR.fit)[3])
  }
  R.intr <- as.numeric(coef(VPR.fit)[1])
  R.slpe <- as.numeric(coef(VPR.fit)[2])
  R.BH <- NaN
  R.SC <- NaN
  R.SCT <- NaN
  m["VPR.fit", ] <- c(R.slpe, R.tcoef, R.intr,R.pval, R.Rval, R.BH, R.SC, R.SCT)

  ####### may wat to add a nonparametric trend test here

  # ===== Critical threshold test =====
  if (R.pval > sig)
    stop("VPR significance below critical threshold. Variables should not be passed to this function")

  # ==============================================================================================
  # ========== Build results and return them significant VPR/VCR ==========
  # +++++ Pull out the key values +++++
  VPR.resid <- ts(VPR.fit$residuals, start = start(ti), end = end(ti), frequency = f)
  #Create the dummy variable
  dummy <- rep(0, length(VPR.resid))
  dummy[(breakpoint + 1):length(VPR.resid)] <- 1

  segRES.df <- data.frame(year = ti, VPR.residuals = VPR.resid, dummy.var = dummy)

  start = as.integer(start(ti)[1])
  bkp = breakpoint + start - 1

  bpanalysis <- lm(VPR.residuals ~ I(year-(bkp + 0.5)) * dummy.var, segRES.df)
  R2.pval <- glance(bpanalysis)$p.value
  R2.tcoef <- NaN
  R2.Rval <- summary(bpanalysis)$r.square
  R2.intr <- bpanalysis$coefficients[[1]]
  R2.slpe <- bpanalysis$coefficients[[2]]
  R2.BH <- bpanalysis$coefficients[[3]]
  R2.SC <- bpanalysis$coefficients[[4]]
  R2.SCT <- NaN
  m["RESTREND.fit", ] <- c(R2.slpe, R2.tcoef, R2.intr,R2.pval, R2.Rval, R2.BH, R2.SC, R2.SCT)

  # Calculate the total change
  init <- bpanalysis$fitted.values[1]
  fin <- bpanalysis$fitted.values[end(bpanalysis$fitted.values)[1]]
  change <- as.numeric(fin - init)
  if (R2.pval < 0.10) {
    tot.ch = change
  }else{
    tot.ch = 0
  }

  # ===== Build and return the results =====
  overview <- data.frame(
    Method = "segmented.RESTREND", Total.Change = tot.ch,
    Residual.Change = change, VPR.HeightChange = FALSE, model.p = glance(VPR.fit)$p.value,
    residual.p = glance(bpanalysis)$p.value, VPRbreak.p = FALSE, bp.year = bkp
    )
  models <- list(
    CTS.fit = FALSE, BFAST = FALSE, VPR.fit = VPR.fit,
    resid.fit = bpanalysis, segVPR.fit = FALSE
    )
  ts.data <- list(
    CTSR.VI = FALSE, CTSR.RF = FALSE, anu.VI = anu.VI,
    VI.index = VI.index, acu.RF = acu.RF,
    acu.TM = acu.TM, StdVar.RF = FALSE, StdVar.TM = FALSE
    )
  ols.summary <- list(chow.sum = FALSE, chow.ind = FALSE, OLS.table = m)
  acum.df <- FALSE

  return(structure(list(
    summary = overview, ts.data = ts.data, ols.summary = ols.summary,
    TSSRmodels = models, acum.df = acum.df), class = "TSSRESTREND")
    )

}
