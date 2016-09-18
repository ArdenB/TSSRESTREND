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
#'
#' @param anu.VI
#' The annual (Growing season) max VI. if anu.VI=FALSE, it will be calculated from the CTSR.VI. See (___)
#' @param acu.RF
#' The optimal accumulated rainfall for anu.VI. Mut be a object of class 'ts' and of equal length to anu.VI. if anu.RF=FALSE, it will be calculated from ACP.table. see(____)
#' @param VI.index
#' the index of the CTSR.VI ts that the anu.VI values occur at. Must be the same length as anu.VI. NOTE. R indexs from 1 rather than 0.
#' @param sig
#' Significance of all the functions, sig=0.05
#' @param print
#' Prints more details at every step of the procces
#' @param plot
#' creates a plots of the residulas vs time
#'
#' @return summary
#' (To be filled in)
#' @return VPR
#' the lm (See lm_______)
#' @return TSS.RESTREND
#' the lm of the residuals
#' @export
#'

seg.RESTREND <- function(anu.VI, acu.RF, VI.index, breakpoint,  sig=0.05){

  while (TRUE){
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

  #Get the VPR
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


  ####### may wat to add a nonparametric trend test here

  #Critical threshold test
  if (summary(VPR.fit)$coefficients[,4][2] > sig)
    stop("VPR significance below critical threshold. Variables should not be passed to this function")

  VPR.resid<- ts(VPR.fit$residuals, start=start(ti), end=end(ti), frequency = f)
  #Create the dummy variable
  dummy <- rep(0, length(VPR.resid))
  dummy[(breakpoint+1):length(VPR.resid)] = 1

  segRES.df = data.frame(year=ti, VPR.residuals=VPR.resid,  dummy.var=dummy)

  start = as.integer(start(ti)[1])
  bkp = breakpoint + start-1

  bpanalysis<-lm(VPR.residuals~I(year-(bkp+0.5))*dummy.var,segRES.df)
  R2.pval <- glance(bpanalysis)$p.value
  R2.Rval <- summary(bpanalysis)$r.square
  R2.intr <- bpanalysis$coefficients[[1]]
  R2.slpe <- bpanalysis$coefficients[[2]]
  R2.BH <- bpanalysis$coefficients[[3]]
  R2.SC <- bpanalysis$coefficients[[4]]
  m["RESTREND.fit", ] <- c(R2.slpe, R2.intr,R2.pval, R2.Rval, R2.BH, R2.SC)


  init <- bpanalysis$fitted.values[1]
  fin <- bpanalysis$fitted.values[end(bpanalysis$fitted.values)[1]]
  change <- as.numeric(fin - init)
  overview <- data.frame(Method = "segmented.RESTREND",
                         Residual.Change=change, VPR.HeightChange =FALSE, model.p = glance(VPR.fit)$p.value,
                         residual.p = glance(bpanalysis)$p.value, VPRbreak.p = FALSE, bp.year=bkp)
  models <- list(CTS.fit=FALSE, BFAST=FALSE, VPR.fit=VPR.fit, resid.fit = bpanalysis, segVPR.fit=FALSE)
  ts.data <- list(CTSR.VI=FALSE, CTSR.RF=FALSE, anu.VI = anu.VI, VI.index = VI.index, acu.RF = acu.RF, StdVar.RF=FALSE)
  ols.summary <- list(chow.sum=FALSE, OLS.table=m)
  return(structure(list(summary=overview, ts.data = ts.data, ols.summary=ols.summary,
                        TSSRmodels=models), class = "TSSRESTREND"))

}
