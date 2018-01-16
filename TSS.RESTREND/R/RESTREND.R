#' @title RESTREND (RESidual TREND)
#'
#' @description
#' For ts with no significant breakpoints in the residuals or the VPR. Takes annual VI max and it associated optimal accumulated precipitation
#'
#' @importFrom stats coef end frequency lm sd start time ts
#' @importFrom graphics abline arrows legend par plot
#' @importFrom utils tail
#' @importFrom broom glance
#' @author Arden Burrell, arden.burrell@unsw.edu.au
#'
#' @inheritParams TSSRESTREND
#'
#' @return a list of class TSSRESTREND.
#'        See \code{\link{TSSRESTREND}} for details. Note. if called seperatly from TSSRESTREND,
#'        this list will be incomplete.
#' @export
#'
#' @examples
#' restrend <- RESTREND(stdRESTREND$max.NDVI, stdRESTREND$acc.precip, stdRESTREND$index)
#' print(restrend)

RESTREND <- function(anu.VI, acu.RF, acu.TM, VI.index, sig=0.05) {
  # ==============================================================================================
  # ========== Sanity check the input data ==========
  while (TRUE) {
    if (class(anu.VI) != "ts") {
      stop("anu.VI Not a time series object")}
    if (class(acu.RF) != "ts") {
      stop("acu.VI Not a time series object")}
    ti <- time(anu.VI)
    f <- frequency(anu.VI)
    #check the two ts object cover the same time period
    ti2 <- time(acu.RF)
    f2 <- frequency(acu.RF)
    if (!identical(ti, ti2)) {
      stop("ts object do not have the same time")}
    if (!identical(f, f2)) {
      stop("ts object do not have the same frequency")}
    break
  }

  # ==============================================================================================
  # ========== Calculate the VPR/VCR ==========
  # Future plans:
  #     Add the option for a non parametric version of this
  #     Add a script that just perform a RESTREND on the raw inputs
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
    "slope", "temp.coef", "intercept", "p.value", "R^2.Value", "Break.Height",
    "Slope.Change", "Slope.ChangeTmp"
    )
  # ========== Pull out the key values ==========
  R.pval <- glance(VPR.fit)$p.value #Pval
  R.Rval <- summary(VPR.fit)$r.square #r2val
  if (is.null(acu.TM)) {# No Temperature data
    R.tcoef <- NaN
  } else {# Temperature data
    R.tcoef <- as.numeric(coef(VPR.fit)[3])
  }
  R.intr <- as.numeric(coef(VPR.fit)[1])
  R.slpe <- as.numeric(coef(VPR.fit)[2])
  R.BH <- NaN
  R.SC <- NaN
  R.SCT <- NaN
  m["VPR.fit", ] <- c(R.slpe, R.tcoef, R.intr,R.pval, R.Rval, R.BH, R.SC, R.SCT)

  # ==============================================================================================
  # ========== Critical threshold and slope tests ==========
  # Catches inderterminate and agriculture results

  if (R.pval > sig) {
    # VPR/VCR failed
    print("VPR significance below critical threshold")
    tot.ch <- FALSE
    change <- FALSE
    overview <- data.frame(
      Method = "indeterminate", Total.Change = tot.ch,
      Residual.Change = change,VPR.HeightChange = FALSE,
      model.p = glance(VPR.fit)$p.value, residual.p = FALSE,
      VPRbreak.p = FALSE, bp.year = FALSE
      )
    models <- list(
      CTS.fit = FALSE, BFAST = FALSE, VPR.fit = VPR.fit,
      resid.fit = FALSE, segVPR.fit = FALSE
      )
    ts.data <- list(
      CTSR.VI = FALSE, CTSR.RF = FALSE, anu.VI = anu.VI,
      VI.index = VI.index, acu.RF = acu.RF, StdVar.RF = FALSE
      )
    ols.summary <- list(chow.sum = FALSE, chow.ind = FALSE, OLS.table = m)

    return(structure(list(
      summary = overview, ts.data = ts.data, ols.summary = ols.summary,
      TSSRmodels = models), class = "TSSRESTREND")
    )

  } else if ((R.slpe < 0) && (is.null(acu.TM))) {
    # This may need to be changes to accout for allow negatives in the future
    # VPR/VCR is significant but negative
    print("VPR slope is negative")
    tot.ch <- FALSE
    change <- FALSE

    overview <- data.frame(
      Method = "IND-agr?", Total.Change = tot.ch,
      Residual.Change = change, VPR.HeightChange = FALSE,
      model.p = glance(VPR.fit)$p.value, residual.p = FALSE,
      VPRbreak.p = FALSE, bp.year = FALSE
      )
    models <- list(
      CTS.fit = FALSE, BFAST = FALSE, VPR.fit = VPR.fit,
      resid.fit = FALSE, segVPR.fit = FALSE
      )
    ts.data <- list(
      CTSR.VI = FALSE, CTSR.RF = FALSE, anu.VI = anu.VI,
      VI.index = VI.index, acu.RF = acu.RF, StdVar.RF = FALSE
      )
    ols.summary <- list(chow.sum = FALSE, chow.ind = FALSE, OLS.table = m)

    return(structure(list(
      summary = overview, ts.data = ts.data, ols.summary = ols.summary,
      TSSRmodels = models), class = "TSSRESTREND")
    )
  }

  # ==============================================================================================
  # ========== Build results and return them significant VPR/VCR ==========
  # +++++ Pull out the key values +++++
  VPR.resid <- ts(VPR.fit$residuals, start = start(ti), end = end(ti), frequency = f)
  RES <- lm(VPR.resid ~ ti)
  R2.pval <- glance(RES)$p.value
  R2.tcoef <- NaN
  R2.Rval <- summary(RES)$r.square
  R2.intr <- RES$coefficients[[1]]
  R2.slpe <- RES$coefficients[[2]]
  R2.BH <- NaN
  R2.SC <- NaN
  R2.SCT <- NaN
  m["RESTREND.fit", ] <- c(R2.slpe,R2.tcoef, R2.intr,R2.pval, R2.Rval, R2.BH, R2.SC, R2.SCT)

  # Calculate the total change
  init <- RES$fitted.values[1]
  fin <- RES$fitted.values[end(RES$fitted.values)[1]]
  change <- fin - init
  if (R2.pval < 0.10) {
    tot.ch = change
  } else {
    tot.ch = 0
  }

  # ===== Build and return the results =====
  overview <- data.frame(
    Method = "RESTREND", Total.Change = tot.ch,
    Residual.Change = change, VPR.HeightChange = FALSE,
    model.p = glance(VPR.fit)$p.value, residual.p = glance(RES)$p.value,
    VPRbreak.p = FALSE, bp.year = FALSE
    )
  models <- list(
    CTS.fit = FALSE, BFAST = FALSE, VPR.fit = VPR.fit,
    resid.fit = RES, segVPR.fit = FALSE
    )
  ts.data <- list(
    CTSR.VI = FALSE, CTSR.RF = FALSE, anu.VI = anu.VI,
    VI.index = VI.index, acu.RF = acu.RF, acu.TM = acu.TM,
    StdVar.RF = FALSE, StdVar.TM = FALSE
    )
  ols.summary <- list(chow.sum = FALSE, chow.ind = FALSE, OLS.table = m)
  acum.df <- FALSE

  return(structure(list(
    summary = overview, ts.data = ts.data, ols.summary = ols.summary,
    TSSRmodels = models, acum.df = acum.df), class = "TSSRESTREND")
  )

}



