#' @title Segmented Vegetation Climate Relationship
#'
#' @description
#' For a ts with a significant breakpoints in the the VPR/VCR. This function takes annual VI max,
#' the optimal accumulated precipitation (& temperature) before and after the breakpoint, then
#' caculates the Standard Variance of the climate cariables.Theen an OLS is performed with a dummy
#' variable to reperesent the breakpoint (0 before the breakpoint and 1 after it)..
#'
#' @author Arden Burrell, arden.burrell@unsw.edu.au
#'
#' @importFrom strucchange sctest
#'
#' @inheritParams TSSRESTREND
#' @param tm.b4
#'        If a breakpoint in the VCR is detected this is the optimial accumulated temperature before
#'        the breakpoint. It must be the same length as the anu.VI. If ACT.table is provided it will
#'        be generated using \code{\link{AnnualClim.Cal}}
#' @param tm.af
#'        If a breakpoint in the VCR is detected this is the optimial accumulated temperature after
#'        the breakpoint. It must be the same length as the anu.VI. If ACT.table is provided it will
#'        be generated using \code{\link{AnnualClim.Cal}}
#' @param breakpoint
#'        The index of the most significant breakpoint as determined using \code{\link{CHOW}}.
#'
#' @return a list of class TSSRESTREND.
#'        See \code{\link{TSSRESTREND}} for details. Note. if called seperatly from TSSRESTREND,
#'        this list will be incomplete.
#' @export
#'
#' @examples
#' brkp <- as.integer(24) #calculated using th CHOW (DONTRUN) example
#' VPRres <- seg.VPR(segVPR$max.NDVI, segVPR$acum.RF, segVPR$index, brkp, segVPR$RFB4, segVPR$RFAF)
#' print(VPRres)

seg.VPR <- function(
  anu.VI, acu.RF, VI.index, breakpoint, rf.b4, rf.af,
  acu.TM = NULL,  tm.b4 = NULL, tm.af = NULL, sig = 0.05
  ) {
  # ==============================================================================================
  # ========== Sanity check the input data ==========
  # Description:
  #   Each check liiks at a different paramter. If the data fails
  #   the check will stop, else, it breaks after all the checks

  while (TRUE) {
    if (class(anu.VI) != "ts")
      stop("anu.VI Not a time series object")
    if (class(acu.RF) != "ts")
      stop("acu.VI Not a time series object")
    if (class(rf.b4) != "ts")
      stop("rf.b4 Not a time series object")
    if (class(rf.af) != "ts")
      stop("rf.af Not a time series object")
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
    if (class(rf.b4) != "logical") {
      if (length(rf.b4) != (length(rf.af)))
        stop("rf.b4 and rf.af are different shapes. They must be the same size and be th same lenths as acu.VI")
    }
    break
  }
  len <- length(anu.VI)
  # ==============================================================================================
  # ========== Setup for result storage and comparison ==========

  # ===== Get the OLS as if there was breapoint in the VPR/VCR =====
  if (is.null(acu.TM)) {# No Temp
    VPR.fit <- lm(anu.VI ~ acu.RF)
    R.tcoef <- NaN
  } else {# temp
    VPR.fit <- lm(anu.VI ~ acu.RF + acu.TM)
    R.tcoef <- as.numeric(coef(VPR.fit)[3])
  }
  # ===== Set up a blank table to hold the regression coefficents from the different models =====
  m <- matrix(nrow = (4), ncol = 8)
  m[] <- NaN
  rownames(m) <- c("CTS.fit", "VPR.fit", "RESTREND.fit", "segVPR.fit")
  colnames(m) <- c(
    "slope", "temp.coef", "intercept", "p.value", "R^2.Value",
    "Break.Height", "Slope.Change", "Slope.ChangeTmp"
    )
  # Pull out key values from VPR and put them in the table
  R.pval <- glance(VPR.fit)$p.value
  R.Rval <- summary(VPR.fit)$r.square
  R.intr <- as.numeric(coef(VPR.fit)[1])
  R.slpe <- as.numeric(coef(VPR.fit)[2])
  R.BH <- NaN
  R.SC <- NaN
  R.SCT <- NaN
  m["VPR.fit", ] <- c(R.slpe, R.tcoef, R.intr,R.pval, R.Rval, R.BH, R.SC, R.SCT)

  # ===== Convert accumulllated climate variables to an SV value to allow comparison across the BP =====
  #calculate the standard Variance before and after the breakpoint
  adj.rfb4 <- array((rf.b4 - mean(rf.b4)))
  sd.adjb4 <- adj.rfb4/sd(rf.b4)

  adj.rfaf <- array((rf.af - mean(rf.af)))
  sd.adjaf <- adj.rfaf/sd(rf.af)
   if (!is.null(acu.TM)) {# has temperature data
    #calculate the standard Variance before and after the breakpoint
    adj.tmb4 <- array((tm.b4 - mean(tm.b4)))
    sd.adjtmb4 <- adj.tmb4/sd(tm.b4)

    adj.tmaf <- array((tm.af - mean(tm.af)))
    sd.adjtmaf <- adj.tmaf/sd(tm.af)
   }
  # ========== Sanity check the standard variance results to avoid OLS failure ===========
  #Check and see if there is any problems with the precipitation data
  if (sd(sd.adjaf[(breakpoint + 1):len]) == 0 || sd(sd.adjb4[1:breakpoint]) == 0) {

    print("Segmented VPR failure. No variance in one segement of the precipitation data")
    tot.ch <- FALSE
    change <- FALSE
    # +++++ Create fail objects to return +++++
    overview <- data.frame(
      Method = "ind-rfsegFail", Total.Change = tot.ch, Residual.Change = change,
      VPR.HeightChange = FALSE, model.p = glance(VPR.fit)$p.value, residual.p = FALSE,
      VPRbreak.p = FALSE, bp.year = FALSE
      )
    models <- list(CTS.fit = FALSE, BFAST = FALSE, VPR.fit = VPR.fit, resid.fit = FALSE, segVPR.fit = FALSE)
    ts.data <- list(
      CTSR.VI = FALSE, CTSR.RF = FALSE, anu.VI = anu.VI, VI.index = VI.index,
      acu.RF = acu.RF, acu.TM = acu.TM, StdVar.RF = adj.RF, StdVar.TM = FALSE
      )
    ols.summary <- list(chow.sum = FALSE, chow.ind = FALSE, OLS.table = m)
    # Return the fail objects
    return(structure(list(
      summary = overview, ts.data = ts.data, ols.summary = ols.summary,TSSRmodels = models),
      class = "TSSRESTREND"))
  }
  # ==============================================================================================
  # ========== Perform a Segmented VPR/VCR calculation and store the results ==========

  # ===== Setup the regression variables =====
  # Build a single adjusted precip vector
  adj.RF <- c(sd.adjb4[1:breakpoint], sd.adjaf[(breakpoint + 1):len])
  # Create the dummy variable
  dummy <- rep(0, len)
  dummy[(breakpoint + 1):len] = 1

  # +++++ THe regression if different if temperature is included or not +++++
  if (!is.null(acu.TM)) {
    # ===== has temperature data =====
    # Build a single standard variance of temperature
    adj.tm <- c(sd.adjtmb4[1:breakpoint], sd.adjtmaf[(breakpoint + 1):len])
    # Put the regression variables in a single dataframe
    segRES.df = data.frame(
      year = ti, VI = anu.VI, sv.RF = adj.RF, sv.TM = adj.tm, breakpoint.var = dummy
      )
    # ===== perform the regression to calculate the segmented VCR =====
    segVPR.fit <-  lm(VI~(sv.RF+sv.TM)*breakpoint.var, segRES.df)

    # ===== Get the regression coefficents that are not common to the VPR =====
    R2.intr <- segVPR.fit$coefficients[[1]]
    R2.slpe <- segVPR.fit$coefficients[[2]]
    R2.tcoef <- segVPR.fit$coefficients[[3]]
    R2.BH <- segVPR.fit$coefficients[[4]]
    R2.SC <- segVPR.fit$coefficients[[5]]
    R2.SCT <- segVPR.fit$coefficients[[6]]
  }else{
    # ===== no temperature data ======
    # Put the regression variables in a single dataframe
    segRES.df = data.frame(
      year = ti, VI = anu.VI, sv.RF = adj.RF, breakpoint.var = dummy
      )
    # ===== perform the regression to calculate the segmented VPR =====
    segVPR.fit <-  lm(VI~sv.RF*breakpoint.var, segRES.df)
    # ===== Get the regression coeeficents that are not common to the VCR =====
    R2.intr <- segVPR.fit$coefficients[[1]]
    R2.slpe <- segVPR.fit$coefficients[[2]]
    R2.tcoef <- NaN
    R2.BH <- segVPR.fit$coefficients[[3]]
    R2.SC <- segVPR.fit$coefficients[[4]]
    R2.SCT <- NaN
  }

  # ===== Get the infomation that is called the same way from the VPR and the VCR =====
  start = as.integer(start(ti)[1])
  bkp = breakpoint + start - 1
  R2.pval <- glance(segVPR.fit)$p.value
  R2.Rval <- summary(segVPR.fit)$r.square
  # add to the summary
  m["segVPR.fit", ] <- c(R2.slpe, R2.tcoef, R2.intr, R2.pval, R2.Rval, R2.BH, R2.SC, R2.SCT)

  # ==============================================================================================
  # ========== Perform a Segmented VPR/VCR residual calculation and store the results ==========
  # DESCRIPTION:
  #     Using the residuals of the segmented VPR/VCR, Performs a segmenented Restrend type
  #     calculation and returns the results as an object of class TSSRESTREND

  # ===== Pull out the variables to go in the segRESTREND =====
  resid.raw <- segVPR.fit$residuals #VI Residuals
  len <- length(resid.raw) #Length of the vectors
  start = as.integer(start(ti)[1]) # start date
  end = as.integer(end(ti)[1]) #end date
  year = c(start:end) # years
  VPR.resid <- ts(resid.raw, start = start(ti), end = end(ti), frequency = f) # timseries
  # +++++ Perform a chow test on the breakpoint in the segmented residuals +++++
  RESchow <- sctest(resid.raw ~ year, type = "Chow", point = breakpoint)

  # Create the dummy variable
  dummy <- rep(0, length(VPR.resid))
  dummy[(breakpoint + 1):length(VPR.resid)] = 1

  # +++++ Put the variables in a single dataframe +++++
  segRES.df = data.frame(year = ti, VPR.residuals = VPR.resid,  dummy.var = dummy)

  start = as.integer(start(ti)[1])
  bkp = breakpoint + start - 1

  # ===== perform the regression to calculate the segmented RESTREND =====
  bpanalysis <- lm(VPR.residuals~I(year-(bkp + 0.5)) * dummy.var, segRES.df)

  # ===== Store the segRES results =====
  R3.pval <- glance(bpanalysis)$p.value
  R3.Rval <- summary(bpanalysis)$r.square
  R3.intr <- bpanalysis$coefficients[[1]]
  R3.tcoef <- NaN
  R3.slpe <- bpanalysis$coefficients[[2]]
  R3.BH <- bpanalysis$coefficients[[3]]
  R3.SC <- bpanalysis$coefficients[[4]]
  R3.SCT <- NaN
  m["RESTREND.fit", ] <- c(R3.slpe, R3.tcoef, R3.intr, R3.pval, R3.Rval, R3.BH, R3.SC, R3.SCT)

  breakheight <- m["segVPR.fit", "Break.Height"]
  bp.pval <- coef(summary(segVPR.fit))["breakpoint.var","Pr(>|t|)"]

  # work out the total  residual change (to add to get the total change)
  init <- bpanalysis$fitted.values[1]
  fin <- bpanalysis$fitted.values[end(bpanalysis$fitted.values)[1]]
  change <- as.numeric(fin - init)
  # ===== store all the relevant time series values ====
  if (!is.null(acu.TM)) {
    # has temperature data
    ts.data <- list(
      CTSR.VI = FALSE, CTSR.RF = FALSE, anu.VI = anu.VI, VI.index = VI.index,
      acu.RF = acu.RF, acu.TM = acu.TM, StdVar.RF = adj.RF, StdVar.TM = adj.tm)
  }else{
    # Has no temerature data
    ts.data <- list(
      CTSR.VI = FALSE, CTSR.RF = FALSE, anu.VI = anu.VI, VI.index = VI.index,
      acu.RF = acu.RF, acu.TM = NULL, StdVar.RF = adj.RF, StdVar.TM = NULL)
  }

  # ===== Work out the total change =====
  tc <- c(0, 0)
  # Check if residual change and the VPR/VCR breakheights meet 0.10 significance levels
  if (R3.pval < 0.10) {
    tc[1] = change
  }

  if (bp.pval < 0.10) {
    tc[2] = breakheight
  }
  # Calculate the total change
  tot.ch = sum(tc)

  # ===== put all the results in dataframes and return them =====
  overview <- data.frame(
    Method = "segmented.VPR", Total.Change = tot.ch, Residual.Change = change,
    VPR.HeightChange = breakheight, model.p = glance(segVPR.fit)$p.value,
    residual.p = R3.pval,VPRbreak.p = bp.pval, bp.year = bkp
    )
  models <- list(
    CTS.fit = FALSE, BFAST = FALSE, VPR.fit = VPR.fit,
    resid.fit = bpanalysis, segVPR.fit = segVPR.fit
    )
  ols.summary <- list(chow.sum = FALSE, chow.ind = FALSE, OLS.table = m)
  acum.df <- FALSE
  # Return a list structure of the results
  return(structure(list(
    summary = overview, ts.data = ts.data, ols.summary = ols.summary,
    TSSRmodels = models, acum.df = acum.df), class = "TSSRESTREND")
    )

}
