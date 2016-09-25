#' @title Segmented Vegetation Precipitation Relationship
#'
#' @description
#' For ts with a significant breakpoints in the the VPR. Takes annual VI max, the optimal accumulated precipitation
#' bor and after the breakpoint, then caculated the Precipitation Standard Variance, and regresses it with a dummy
#' variable that is 0 before the breakpoint and 1 after it. WARNING NEED TO REMOVE acu.RF
#'
#' @author Arden Burrell, arden.burrell@unsw.edu.au
#'
#' @importFrom strucchange sctest
#'
#' @inheritParams TSSRESTREND
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

seg.VPR <- function(anu.VI, acu.RF, VI.index, breakpoint, rf.b4, rf.af, sig=0.05){
  while (TRUE){
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
    if (class(rf.b4) != "logical"){
      if (length(rf.b4) != (length(rf.af)))
        stop("rf.b4 and rf.af are different shapes. They must be the same size and be th same lenths as acu.VI")
    }
    break
  }
  len <- length(anu.VI)

  VPR.fit <- lm(anu.VI ~ acu.RF)

  m<- matrix(nrow=(4), ncol=6)
  m[]<-NaN
  rownames(m)<- c("CTS.fit", "VPR.fit", "RESTREND.fit", "segVPR.fit")
  colnames(m)<- c("slope", "intercept", "p.value", "R^2.Value", "Break.Height", "Slope.Change")

  R.pval <- glance(VPR.fit)$p.value
  R.Rval <- summary(VPR.fit)$r.square
  R.intr <- as.numeric(coef(VPR.fit)[1])
  R.slpe <- as.numeric(coef(VPR.fit)[2])
  R.BH <- NaN
  R.SC <- NaN
  m["VPR.fit", ] <- c(R.slpe, R.intr,R.pval, R.Rval, R.BH, R.SC)

  #calculate the standard Variance before and after the breakpoint
  adj.rfb4 <- array((rf.b4-mean(rf.b4)))
  sd.adjb4 <- adj.rfb4/sd(rf.b4)

  adj.rfaf <- array((rf.af-mean(rf.af)))
  sd.adjaf <- adj.rfaf/sd(rf.af)
  adj.RF <- c(sd.adjb4[1:breakpoint], sd.adjaf[(breakpoint+1):len])

  #Create the dummy variable
  dummy <- rep(0, len)
  dummy[(breakpoint+1):len] = 1

  segRES.df = data.frame(year=ti, VI=anu.VI, sv.RF=adj.RF,   breakpoint.var=dummy)

  start = as.integer(start(ti)[1])
  bkp = breakpoint + start-1

  #perform the regression
  segVPR.fit <-  lm(VI~sv.RF*breakpoint.var, segRES.df)
  R2.pval <- glance(segVPR.fit)$p.value
  R2.Rval <- summary(segVPR.fit)$r.square
  R2.intr <- segVPR.fit$coefficients[[1]]
  R2.slpe <- segVPR.fit$coefficients[[2]]
  R2.BH <- segVPR.fit$coefficients[[3]]
  R2.SC <- segVPR.fit$coefficients[[4]]
  m["segVPR.fit", ] <- c(R2.slpe, R2.intr,R2.pval, R2.Rval, R2.BH, R2.SC)
  # browser()

  #Added new section for creating a total height
  resid.raw <- segVPR.fit$residuals
  # resid.BHadj <- c(resid.raw[1:breakpoint], resid.raw[(breakpoint+1):len] + R2.BH)
  len <- length(resid.raw)
  start = as.integer(start(ti)[1])
  end = as.integer(end(ti)[1])
  year = c(start:end)
  VPR.resid<- ts(resid.raw, start=start(ti), end=end(ti), frequency = f)
  RESchow <- sctest(resid.raw ~ year, type = "Chow", point = breakpoint)



  #Create the dummy variable
  dummy <- rep(0, length(VPR.resid))
  dummy[(breakpoint+1):length(VPR.resid)] = 1

  segRES.df = data.frame(year=ti, VPR.residuals=VPR.resid,  dummy.var=dummy)

  start = as.integer(start(ti)[1])
  bkp = breakpoint + start-1

  bpanalysis<-lm(VPR.residuals~I(year-(bkp+0.5))*dummy.var,segRES.df)
  R3.pval <- glance(bpanalysis)$p.value
  R3.Rval <- summary(bpanalysis)$r.square
  R3.intr <- bpanalysis$coefficients[[1]]
  R3.slpe <- bpanalysis$coefficients[[2]]
  R3.BH <- bpanalysis$coefficients[[3]]
  R3.SC <- bpanalysis$coefficients[[4]]
  m["RESTREND.fit", ] <- c(R3.slpe, R3.intr, R3.pval, R3.Rval, R3.BH, R3.SC)



  breakheight <- segVPR.fit$coefficients[[3]]
  bp.pval <- coef(summary(segVPR.fit))[15]

  #from the residual change
  init <- bpanalysis$fitted.values[1]
  fin <- bpanalysis$fitted.values[end(bpanalysis$fitted.values)[1]]
  change <- as.numeric(fin - init)

  ts.data <- list(CTSR.VI=FALSE, CTSR.RF=FALSE, anu.VI = anu.VI, VI.index = VI.index, acu.RF = acu.RF, StdVar.RF=adj.RF)
  #May ad a total residual change
  tc <- c(0, 0)
  if (R3.pval<0.10){
    tc[1] = change
  }

  if (bp.pval<0.10){
    tc[2] = breakheight
  }
  tot.ch = sum(tc)
  overview <- data.frame(Method = "segmented.VPR", Total.Change=tot.ch,
                         Residual.Change=change, VPR.HeightChange=breakheight, model.p = glance(segVPR.fit)$p.value,
                         residual.p = R3.pval, VPRbreak.p = bp.pval, bp.year=bkp)
  models <- list(CTS.fit=FALSE, BFAST=FALSE, VPR.fit=VPR.fit, resid.fit = bpanalysis, segVPR.fit=segVPR.fit)
  ols.summary <- list(chow.sum=FALSE, chow.ind=FALSE, OLS.table=m)

  return(structure(list(summary=overview, ts.data = ts.data, ols.summary=ols.summary,
                        TSSRmodels=models), class = "TSSRESTREND"))

}
