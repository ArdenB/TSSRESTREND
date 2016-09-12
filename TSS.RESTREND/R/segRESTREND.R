
seg.RESTREND <- function(anu.VI, acu.RF, VI.index, breakpoint,  sig=0.05, print=FALSE, plot=FALSE){
  
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
  #Convert to a ts object
  #may wat to add a nonparametric trend test here
  
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
  if (print){
    print(summary(bpanalysis))
  }
  
  if (plot){
    len <- length(VPR.resid)
    R.fit <- lm(VPR.residuals ~ year, segRES.df)
    R.fit0 <- lm(VPR.residuals[1:breakpoint] ~ year[1:breakpoint], segRES.df)
    R.fit1 <- lm(VPR.residuals[breakpoint+1:len] ~ year[breakpoint+1:len], segRES.df)
    
    start = as.integer(start(ti)[1])
    end = as.integer(end(ti)[1])
    xlim = c(start, end)
    m.range = 2*max(abs(bpanalysis$fitted.values))
    
    RESchow <- sctest(segRES.df$VPR.residuals ~ segRES.df$year, type = "Chow", point = breakpoint)
    plot(segRES.df$year[1:breakpoint], segRES.df$VPR.residuals[1:breakpoint], pch=16,
         xlab="time", ylab="Residuals", col="orange", xlim=xlim, ylim = c(-m.range, m.range))
    par(new=T)   
    plot(segRES.df$year[breakpoint+1:len], segRES.df$VPR.residuals[breakpoint+1:len], pch=16,
         xlab="", ylab="", col="purple", main="", 
         xlim=c(start, end), ylim = c(-m.range, m.range))
    abline(R.fit, col = "red",lwd = 2, lty = "dashed")
    
    
    par(new=T)   
    bpa.fitts <- ts(bpanalysis$fitted.values, start=ti[1], end=tail(ti, 1), frequency = f) 
    b4.bp = bpanalysis$coefficients[[1]]
    af.bp = bpanalysis$coefficients[[1]] + bpanalysis$coefficients[[3]]
    bpats2 <- append(bpa.fitts, c(b4.bp, af.bp), after = breakpoint)
    t2 <- append(ti, c(start+breakpoint-0.50001, start+breakpoint-0.49999), after=breakpoint)
    plot(t2, bpats2, pch=16, type = "l", lwd = 2,
         xlab="", ylab="", col="red", main="", 
         xlim=c(start, end), ylim = c(-m.range, m.range))
    #add a breakpoint band
    abline(v=(breakpoint-0.5+start), col="white", lwd = 3, lty = "dotted")
    #Need to change the stastics that is shows here
    
    R.Fval = summary(bpanalysis)$f[[1]]
    R.pval = glance(bpanalysis)$p.value
    R.Rval = summary(bpanalysis)$r.squared
    rp = vector('expression',3)
    rp[1] = substitute(expression(italic(R^2) == R.Rval), 
                       list(R.Rval = format(R.Rval,dig=3)))[2]
    rp[2] = substitute(expression(italic(F) == R.Fval), 
                       list(R.Fval  = format(R.Fval,dig=3)))[2]
    
    rp[3] = substitute(expression(italic(p) == R.pval), 
                       list(R.pval  = format(R.pval, digits = 3)))[2]
    legend('topleft', legend = rp, bty = 'n')
  }
  
  init <- bpanalysis$fitted.values[1]
  fin <- bpanalysis$fitted.values[end(bpanalysis$fitted.values)[1]]
  change <- as.numeric(fin - init)
  overview <- data.frame(Method = "segmented.RESTREND", 
                         Total.Change=change, model.p = glance(VPR.fit)$p.value, 
                         residual.p = glance(bpanalysis)$p.value, VPRbreak.p = FALSE)
  return(structure(list(summary=overview, 
                        VPR = VPR.fit, TSS.RESTREND = bpanalysis), class = "RESTREND.Object"))
  
}
