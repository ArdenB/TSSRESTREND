library("bfast")
#library("forecast")
library("RcppCNPy")
#library("gap", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
library("strucchange")

setwd("/mnt/FCBE3028BE2FD9C2/Users/user/Documents/segres_demo")

#load the data

load("./demo_data/stdRESTREND.Rda")
load("./demo_data/stdRESTREND_CTSR.Rda")



BFAST.RESID <- function(CTSR.VI, CTSR.RF, print=FALSE, plot=FALSE) {
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
    return(tmp)
  }else {return(FALSE)}
  
}
# Missing the chow test

RESTREND <- function(anu.VI, acu.RF, VI.index, sig=0.05, breakpoint=FALSE, print=FALSE, plot=FALSE) {
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
  if (!breakpoint){
    VPR.fit <- lm(anu.VI ~ acu.RF)
    #Convert to a ts object
    #may wat to add a nonparametric trend test here
    
    
    #Critical threshold test 
    if (summary(VPR.fit)$coefficients[,4][2] > sig){
      print("VPR significance below critical threshold")
      return(structure(list(Method = FALSE, 
                            VPR = VPR.fit, TSS.RESTREND = FALSE, total.change = FALSE), class = "RESTREND Object"))
    } 
    VPR.resid<- ts(VPR.fit$residuals, start=start(ti), end=end(ti), frequency = f)
    RES <- lm(VPR.resid ~ ti)
    if (plot){
      plot(c(start(ti)[1]:end(ti)[1]), VPR.fit$residuals, pch=16,xlab="Accumulated Rainfall", 
           ylab="Annual VI", col="orange",main="VI vs Precip")
      abline(RES, col = "red",lwd = 2, lty = "dashed")# need to add a trend line 
    }
    init <- RES$fitted.values[1]
    fin <- RES$fitted.values[end(RES$fitted.values)[1]]
    change <- fin - init
    return(structure(list(Method = "RESTREND", 
                          VPR = VPR.fit, TSS.RESTREND = RES, total.change = change), class = "RESTREND Object"))
    
  } else{
    print("segres")
  }
}

#CTSR.VI <- stdRESTREND.CTSR$cts.NDVI
#CTSR.RF <- stdRESTREND.CTSR$cts.precip
a<- BFAST.RESID(stdRESTREND.CTSR$cts.NDVI, stdRESTREND.CTSR$cts.precip, print=TRUE, plot=TRUE)
print(a)

#anu.VI <- stdRESTREND$max.NDVI
#acu.RF <- stdRESTREND$acc.precip
#VI.index <- stdRESTREND$index

res <- RESTREND(stdRESTREND$max.NDVI, stdRESTREND$acc.precip, stdRESTREND$index, sig=0.05, breakpoint=FALSE, print=TRUE, plot=TRUE) 

