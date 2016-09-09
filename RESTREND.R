

RESTREND <- function(anu.VI, acu.RF, VI.index, sig=0.05, print=FALSE, plot=FALSE) {
  #check the data
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
      stop("ts object do not have the same time")
    if (!identical(f, f2))
      stop("ts object do not have the same frequency")
    break
  }
  
  #Get the VPR
  VPR.fit <- lm(anu.VI ~ acu.RF)
  
  #may wat to add a nonparametric trend test here
  
  #Critical threshold test 
  if (summary(VPR.fit)$coefficients[,4][2] > sig){
    if (print){ print("VPR significance below critical threshold")}
    overview <- data.frame(Method = FALSE, Total.Change=FALSE, model.p = glance(VPR.fit)$p.value,
                           residual.p = FALSR, VPRbreak.p = FALSE)
    return(structure(list(summary=overview, 
                          VPR = VPR.fit, TSS.RESTREND = FALSE), class = "RESTREND.Object"))
  } 
  VPR.resid<- ts(VPR.fit$residuals, start=start(ti), end=end(ti), frequency = f)
  RES <- lm(VPR.resid ~ ti)
  
  #!!!!!!!!!!!!!!!!!!!!!!!! Plots need to be tested and improved !!!!!!!!!!!!!!!!!!!!!!!
  if (plot){
    m.range = 2*max(abs(RES$fitted.values))
    plot(c(start(ti)[1]:end(ti)[1]), VPR.fit$residuals, pch=16,xlab="Accumulated Rainfall", 
         ylab="Annual VI", col="orange",main="VI vs Precip", ylim = c(-m.range, m.range))
    par(new=T)   
    plot(c(start(ti)[1]:end(ti)[1]), RES$fitted.values, type = "l", lwd = 2, pch=16,
         xlab="", ylab="", col="red", main="", ylim = c(-m.range, m.range))
    
    R.Fval = summary(RES)$f[[1]]
    R.Rval = summary(RES)$r.squared
    R.pval = glance(RES)$p.value
    
    rp = vector('expression',3)
    rp[1] = substitute(expression(italic(R^2) == R.Rval), 
                       list(R.Rval = format(R.Rval,dig=3)))[2]
    rp[2] = substitute(expression(italic(F) == R.Fval), 
                       list(R.Fval = format(R.Fval,dig=3)))[2]
    
    rp[3] = substitute(expression(italic(p) == R.pval), 
                       list(R.pval = format(R.pval, digits = 3)))[2]
    legend('topleft', legend = rp, bty = 'n')
    
    # abline(RES, col = "red",lwd = 2)#, lty = "dashed")# need to add a trend line 
    
  }
  init <- RES$fitted.values[1]
  fin <- RES$fitted.values[end(RES$fitted.values)[1]]
  change <- fin - init
  overview <- data.frame(Method = "RESTREND", Total.Change=change, model.p = glance(VPR.fit)$p.value,
                         residual.p = glance(RES)$p.value, VPRbreak.p = FALSE)
  return(structure(list(summary=overview, 
                        VPR = VPR.fit, TSS.RESTREND = RES), class = "RESTREND.Object"))
}

