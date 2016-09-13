#' @title Segmented Vegetation Precipitation Relationship
#'
#' @description
#' For ts with a significant breakpoints in the the VPR. Takes annual VI max, the optimal accumulated precipitation
#' bor and after the breakpoint, then caculated the Precipitation Standard Variance, and regresses it with a dummy
#' variable that is 0 before the breakpoint and 1 after it. WARNING NEED TO REMOVE acu.RF
#'
#' @importFrom strucchange sctest
#' @author Arden Burrell, arden.burrell@unsw.edu.au
#'
#'
#' @param anu.VI
#' The annual (Growing season) max VI. if anu.VI=FALSE, it will be calculated from the CTSR.VI. See (___)
#' @param acu.RF
#' The optimal accumulated rainfall for anu.VI. Mut be a object of class 'ts' and of equal length to anu.VI. if anu.RF=FALSE, it will be calculated from ACP.table. see(____)
#' @param VI.index
#' the index of the CTSR.VI ts that the anu.VI values occur at. Must be the same length as anu.VI. NOTE. R indexs from 1 rather than 0.
#' @param breakpoint
#' TO BE ADDED
#' @param sig
#' Significance of all the functions, sig=0.05
#' @param print
#' Prints more details at every step of the procces
#' @param plot
#' creates a plots of the VPR vs rainfall standard variance
#'
#' @return summary
#' (To be filled in)
#' @return VPR
#' the lm (See lm_______)
#' @return TSS.RESTREND
#' FALSE
#' @export
#'
seg.VPR <- function(anu.VI, acu.RF, VI.index, breakpoint, rf.b4, rf.af, sig=0.05, print=FALSE, plot=FALSE){
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
  if (plot){
    # df = data.frame(year=(t-1982), RF=adj.RF, NDVI,  breakpoint=z)
    fit0 <- lm(anu.VI[1:breakpoint] ~ sd.adjb4[1:breakpoint])
    fit1 <- lm(anu.VI[(breakpoint+1):len] ~ sd.adjaf[(breakpoint+1):len])
    fitRES <- lm(segRES.df$VI ~ segRES.df$sv.RF)
    # chow1 <- sctest(bpanalysis$residuals ~ t, type = "Chow", point = bkp)
    plt.ymin <- min(segRES.df$VI)
    plt.ymax <- max(segRES.df$VI)
    plt.xmin <- min(segRES.df$sv.RF)
    plt.xmax <- max(segRES.df$sv.RF)
    plot(segRES.df$sv.RF[1:breakpoint], segRES.df$VI[1:breakpoint], pch=16,
         xlab="Rainfall Standard Variance", ylab="Annual max VI", col="orange",
         xlim=c(plt.xmin, plt.xmax), ylim=c(plt.ymin, plt.ymax))
    par(new=T)
    plot(segRES.df$sv.RF[(breakpoint+1):len], segRES.df$VI[(breakpoint+1):len], pch=16,
         xlab="", ylab="", col="purple", main="",
         xlim=c(plt.xmin, plt.xmax), ylim=c(plt.ymin, plt.ymax))
    par(new=T)
    # abline(fit, col = "red",lwd = 2, lty = "dashed")
    abline(fit0, col = "orange", lwd = 2)
    abline(fit1, col = "purple", lwd = 2)
    abline(fitRES, col="red", lwd=2, lty="dashed")
    top <- segVPR.fit$coefficients[[1]]
    bh <-  segVPR.fit$coefficients[[3]]
    bot <- top+bh
    # lines(x=c(0, 0), y=c(top, bot), col="red", lwd=2, pch=0)
    arrows(0, bot, x1=0, y1=top, length =.075,  angle = 90, code=3, col="red", lwd=2)
    R.Fval = summary(segVPR.fit)$f[[1]]
    R.pval = glance(segVPR.fit)$p.value
    R.Rval = summary(segVPR.fit)$r.squared
    rp = vector('expression',3)
    rp[1] = substitute(expression(italic(R^2) == R.Rval),
                       list(R.Rval = format(R.Rval,dig=3)))[2]
    rp[2] = substitute(expression(italic(F) == R.Fval),
                       list(R.Fval  = format(R.Fval,dig=3)))[2]

    rp[3] = substitute(expression(italic(p) == R.pval),
                       list(R.pval  = format(R.pval, digits = 3)))[2]
    legend('topleft', legend = rp, bty = 'n')

    # browser()
  }
  breakheight <- segVPR.fit$coefficients[[3]]
  bp.pval <- coef(summary(segVPR.fit))[15]

  overview <- data.frame(Method = "segmented.VPR",
                         Total.Change=breakheight, model.p = glance(segVPR.fit)$p.value,
                         residual.p = FALSE, VPRbreak.p = bp.pval)
  return(structure(list(summary=overview,
                        VPR = segVPR.fit, TSS.RESTREND = FALSE), class = "RESTREND.Object"))

  browser()
}
