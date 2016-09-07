

VPR.BFAST <- function(CTSR.VI, CTSR.RF, print=FALSE, plot=FALSE, details=FALSE) {
  #functions takes the complete time series VI and rainfall (RF)
  
  #Check the objects are Time series
  if (class(CTSR.VI) != "ts") 
    stop("CTSR.VI Not a time series object")
  if (class(CTSR.RF) != "ts") 
    stop("CTSR.VI Not a time series object")
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
  
  #Print and plot
  if (print) { #Need better explination and a title but it will do for the moment
    print(summary(CTS.fit))
  } 
  if (plot){ 
    #This plot need serius work but will do for the moment
    plot(resid.ts)
  }
  
  #perform the BFAST
  bf.fit <- bfast(resid.ts, h=0.15, season="none", max.iter=3, level = 0.05)
  
  if (plot){
    #if plot is requested
    plot(bf.fit)
  }
  
  if (bf.fit$nobp$Vt[[1]] == FALSE) {
    numiter <- length(bf.fit$output)
    tmp <- bf.fit$output[[numiter]]$bp.Vt[1]$breakpoints
    if (details){
      return(structure(list(bkps = tmp, BFAST.obj=bf.fit), class = "BFAST.Object"))
    }else{return(structure(list(bkps = tmp, BFAST.obj=FALSE), class = "BFAST.Object"))
    }
  }else {
    if (details){
      return(structure(list(bkps = FALSE, BFAST.obj=bf.fit), class = "BFAST.Object"))
    }else{return(structure(list(bkps = FALSE, BFAST.obj=FALSE), class = "BFAST.Object"))}
  }
}
