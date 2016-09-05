library("bfast")
#library("forecast")
library("RcppCNPy")
library("strucchange")
library("broom")

setwd("/mnt/FCBE3028BE2FD9C2/Users/user/Documents/segres_demo") #needs to be replaced witha variable function

#load the data

load("./demo_data/stdRESTREND.Rda")
load("./demo_data/stdRESTREND_CTSR.Rda")
load("./demo_data/segRESTREND.Rda")
load("./demo_data/segRESTREND_CTSR.Rda")
load("./demo_data/segVPRD.Rda")
load("./demo_data/segVPRD_CTSR.Rda")
load("./demo_data/segVPRI.Rda")
load("./demo_data/segVPRI_CTSR.Rda")



BFAST.RESID <- function(CTSR.VI, CTSR.RF, print=FALSE, plot=FALSE, details=FALSE) {
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

CHOW <- function(anu.VI, acu.RF, VI.index, breakpoints, sig=0.05, print=FALSE){
  #test the data to make sure its valid
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
  # need to test the breakpoints to make sure they are numeric
  if (class(breakpoints) != "numeric") 
    stop("Breakpoints are not class numeric")
  
  #count of all the breakpoints
  len <- length(breakpoints)
  # convert the breakpoints into year posistion
    #the breakpoint loc will be the last anaual max before the breakpoint
  
  #empty variables to add to a datafram
  empty.1 <- NaN
  empty.2 <- NaN 
  empty.3 <- NaN
  ind.df <- data.frame(abs.index=breakpoints, yr.index = empty.1, reg.sig=empty.2, VPR.bpsig = empty.3)
  for (bp in 1:length(breakpoints)){
    bpv = ind.df$abs.index[bp]
    for (n in 1:length(VI.index)){
      if (bpv>=VI.index[n] & bpv<=VI.index[n+1]){
        #print(n)}
        ind.df$yr.index[bp] = n
        }
    }
  }
  #create the lm for the VPR and test is significance, 
  #split VPR sig from VPR insignificant 
  VPR.fit <- lm(anu.VI ~ acu.RF)
  if (summary(VPR.fit)$coefficients[,4][2] > sig){
    if (print){print("VPR significance below critical threshold, Testing breakpoints in the VPR")}
    ind <- acu.RF #independent variable
    dep <- anu.VI #dependent variable
    Method = "seg.VPR"
  } else {
    ind <- ti
    dep <- VPR.fit$residuals
    Method = "seg.RESTREND"
  } 
  #Iterate over each of the breakpoints
  while (TRUE){
    for (bp.num in nrow(ind.df)){ #the breakpoints number, first bp is 1,  etc
      bp = ind.df$yr.index[bp.num]
      #start and ends
      if (identical(ind.df$yr.index[bp.num-1], numeric(0))){
        bp.start = 1
        #print("here0")
      }else{
        bp.start = ind.df$yr.index[bp.num-1]
        #print("here1")
        }
      if (is.na(ind.df$yr.index[bp.num+1])){
        bp.end = length(dep)
      }else{
        bp.end = ind.df$yr.index[bp.num+1]}
      #perform the chow test
      bkp = bp - (bp.start-1)
      print(bkp)
      chow <- sctest(dep[bp.start:bp.end] ~ ind[bp.start:bp.end], type = "Chow", point = bkp)

      
      ind.df$reg.sig[bp.num] = chow$p.value
      
    }
    if (nrow(ind.df)>1){ #*****This needs to be tested with mutiple breakpoints*****
      #delete breakpoint with the largest p values (lowest significance)
      ind.df <- ind.df[!(1:nrow(ind.df) %in% (which.max(ind.df$reg.sig))),]
    }else if (nrow(ind.df)==1){
      if (print){
        print(chow)
      }
      VPR.chow<- sctest(anu.VI ~ acu.RF, type = "Chow", point = ind.df$yr.index[1])
      ind.df$VPR.bpsig[1] = VPR.chow$p.value
      
      if (Method == "seg.VPR" & ind.df$reg.sig[1] > sig){ # cant chow non-sig residulas (bpRESID.chow = FALSE)
        ind.df$reg.sig[1] = NaN
        return(structure(list(n.Method = FALSE, bp.summary = ind.df,
                              bpRESID.chow = FALSE, bpVPR.chow=VPR.chow), class = "CHOW.Object"))
      }else if (Method == "seg.VPR" & ind.df$reg.sig[1] <= sig){
        ind.df$reg.sig[1] = NaN
        return(structure(list(n.Method = "seg.VPR", bp.summary = ind.df,
                              bpRESID.chow = FALSE, bpVPR.chow=VPR.chow), class = "CHOW.Object"))
      }else if (Method == "seg.RESTREND" & ind.df$reg.sig[1] > sig){
        return(structure(list(n.Method = "RESTREND", bp.summary = ind.df,
                              bpRESID.chow = chow, bpVPR.chow=VPR.chow), class = "CHOW.Object"))
      }else if (Method == "seg.RESTREND" & ind.df$reg.sig[1] <= sig){
        if (VPR.chow$p.value>sig){
          return(structure(list(n.Method = "seg.RESTREND", bp.summary = ind.df,
                              bpRESID.chow = chow, bpVPR.chow=VPR.chow), class = "CHOW.Object"))
        }else if(VPR.chow$p.value<=sig){
          return(structure(list(n.Method = "seg.VPR", bp.summary = ind.df,
                                bpRESID.chow = chow, bpVPR.chow=VPR.chow), class = "CHOW.Object"))
        }
      }else{
        print("Error in Method and df, exit point 1")
        return(FALSE)
      }
    }else{
      print("ind.df shape is wrong, exited to avoid infinit loop, exit point 2")
      return(FALSE)
    }
  }
}



#final methods Section
seg.VPR <- function(anu.VI, acu.RF, VI.index, breakpoint, rf.b4, rf.af, sig=0.05, print=FALSE, plot=FALSE){
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
    fit0 <- lm(NDVI[1:breakpoint] ~ sd.adjb4[1:breakpoint])
    fit1 <- lm(NDVI[(breakpoint+1):len] ~ sd.adjaf[(breakpoint+1):len])
    fitRES <- lm(segRES.df$VI ~ segRES.df$sv.RF)
    # chow1 <- sctest(bpanalysis$residuals ~ t, type = "Chow", point = bkp)
    plt.ymin <- min(segRES.df$VI)
    plt.ymax <- max(segRES.df$VI)
    plt.xmin <- min(segRES.df$sv.RF)
    plt.xmax <- max(segRES.df$sv.RF)
    plot(segRES.df$sv.RF[1:breakpoint], segRES.df$VI[1:breakpoint], pch=16,
         xlab="Rainfall Standard Variance", ylab="Annual max NDVI", col="orange", 
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
    abline(v=0,lty="dashed") #to be improved asap
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
    abline(v=(breakpoint-0.5+start), col="black", lwd = 2, lty = "dashed")
    par(new=T)   
    bpa.fitts <- ts(bpanalysis$fitted.values, start=ti[1], end=tail(ti, 1), frequency = f) 
    b4.bp = bpanalysis$coefficients[[1]]
    af.bp = bpanalysis$coefficients[[1]] + bpanalysis$coefficients[[3]]
    bpats2 <- append(bpa.fitts, c(b4.bp, af.bp), after = breakpoint)
    t2 <- append(ti, c(start+breakpoint-0.50001, start+breakpoint-0.49999), after=breakpoint)
    plot(t2, bpats2, pch=16, type = "l", lwd = 2,
         xlab="", ylab="", col="red", main="", 
         xlim=c(start, end), ylim = c(-m.range, m.range))
    
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


#FUnctions to call functions
TSS.RESTREND <- function(CTSR.VI, CTSR.RF, anu.VI, acu.RF, VI.index, rf.b4=FALSE, rf.af=FALSE, 
                         sig=0.05, print=FALSE, plot=FALSE, details=FALSE){
  #Function to call the other functions
  #Missing function to find optimal accumulation of the precipitation
  #which will be a sperate script. If the method is segVPR, there 
  #is a need to recalculate precip on either side of the breakpoint
  #Until it is functional  b4 and after need to be passed into this
  #function.  rf.b4=FALSE, rf.af=FALSE, will be removed as soon as 
  while (TRUE){ #Test the variables 
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
    if (class(rf.b4) != "logical"){
      if (length(rf.b4) != (length(rf.af)))
          stop("rf.b4 and rf.af are different shapes. They must be the same size and be th same lenths as acu.VI")
    }
    break
  }
  bkp = BFAST.RESID(CTSR.VI, CTSR.RF, print=print, plot=plot)
  bp<-as.numeric(bkp$bkps)
  if (!bp){# no breakpoints detected by the BFAST
    test.Method = "RESTREND"
  }else{
    res.chow <- CHOW(anu.VI, acu.RF, VI.index, bp, sig=sig, print=print)
    test.Method = res.chow$n.Method
  }
  #add NDVI plot
  
  if (test.Method == "RESTREND"){
    result <- RESTREND(anu.VI, acu.RF, VI.index, sig=sig, print=print, plot=plot) 
  }else if (test.Method == "seg.RESTREND"){
    breakpoint = as.integer(res.chow$bp.summary[2])
    result <- seg.RESTREND(anu.VI, acu.RF, VI.index, breakpoint,  sig=sig, print=print, plot=plot)
  }else if (test.Method == "seg.VPR"){
    breakpoint = as.integer(res.chow$bp.summary[2])
    result <- seg.VPR(anu.VI, acu.RF, VI.index, breakpoint, rf.b4, rf.af, sig=sig, print=print, plot=plot)
  }
  browser()
  print(result$summary)
}


demo.stdRESTEND <- function(sig=0.05, print=TRUE, plot=TRUE, details=FALSE, mode="TSS.RESTREND"){
  #set the environment variables 
  CTSR.VI <- stdRESTREND.CTSR$cts.NDVI
  CTSR.RF <- stdRESTREND.CTSR$cts.precip
  anu.VI <- stdRESTREND$max.NDVI
  acu.RF <- stdRESTREND$acc.precip
  VI.index <- stdRESTREND$index
  if (mode == "TSS.RESTREND"){
    TSSR.result <-TSS.RESTREND(CTSR.VI, CTSR.RF, anu.VI, acu.RF, VI.index, sig=sig, print=print, plot=plot, details = details)
  }else{
    #drops into a browser so the user can call the functions individually 
    print("loading standard RESTREND environment variables")
    browser()
  }
}

demo.segRESTEND <- function(sig=0.05, print=TRUE, plot=TRUE, details=FALSE, mode="TSS.RESTREND"){
  #set the environment variables 
  CTSR.VI <- segRESTREND.CTSR$cts.NDVI
  CTSR.RF <- segRESTREND.CTSR$cts.precip
  anu.VI <- segRESTREND$max.NDVI
  acu.RF <- segRESTREND$acc.precip
  VI.index <- segRESTREND$index
  if (mode == "TSS.RESTREND"){
    TSSR.result <-TSS.RESTREND(CTSR.VI, CTSR.RF, anu.VI, acu.RF, VI.index, sig=sig, print=print, plot=plot, details = details)
  }else{
    #drops into a browser so the user can call the functions individually 
    print("loading segmented RESTREND environment variables")
    browser()
    #need to add the other modes 
  }
}

demo.segVPRD <- function(sig=0.05, print=TRUE, plot=TRUE, details=FALSE, mode="TSS.RESTREND"){
  #set the environment variables 
  CTSR.VI <- segVPRD.CTSR$cts.NDVI
  CTSR.RF <- segVPRD.CTSR$cts.precip
  anu.VI <- segVPRD$max.NDVI
  acu.RF <- segVPRD$acc.precip
  VI.index <- segVPRD$index
  rf.b4 <- segVPRD$acp.b4
  rf.af <- segVPRD$acp.af
  
  if (mode == "TSS.RESTREND"){
    TSSR.result <-TSS.RESTREND(CTSR.VI, CTSR.RF, anu.VI, acu.RF, VI.index, rf.b4, rf.af, sig=sig, print=print, plot=plot, details = details)
  }else{
    #drops into a browser so the user can call the functions individually 
    print("loading segmented RESTREND environment variables")
    browser()
    #need to add the other modes 
  }
}


#need to figure add a return 
# res <- demo.stdRESTEND()
# res <- demo.segRESTEND()
res <- demo.segVPRD()

print("hello World")



# a<- BFAST.RESID(stdRESTREND.CTSR$cts.NDVI, stdRESTREND.CTSR$cts.precip, print=TRUE, plot=TRUE)
# print(a)
# 
# 
# #fin functions use class(a) to determine if its numeric or logical
# 
# se<- BFAST.RESID(segRESTREND.CTSR$cts.NDVI, segRESTREND.CTSR$cts.precip, print=TRUE, plot=TRUE)
# print(se)
# #CHOW <- function(anu.VI, acu.RF, VI.index, se, sig=0.05, print=FALSE)
# 
# #for testing the chow test
# anu.VI <- segRESTREND$max.NDVI
# acu.RF <- segRESTREND$acc.precip
# VI.index <- segRESTREND$index
# 
# se.chow <- CHOW(anu.VI, acu.RF, VI.index, se, sig=0.05, print=TRUE)
# breakpoint = as.integer(se.chow$bp.summary[2])
# se.RES <- seg.RESTREND(anu.VI, acu.RF, VI.index, breakpoint,  sig=0.05, print=TRUE, plot=TRUE)
# 
# anu.VI <- stdRESTREND$max.NDVI
# acu.RF <- stdRESTREND$acc.precip
# VI.index <- stdRESTREND$index
# res <- RESTREND(anu.VI, acu.RF, VI.index, sig=0.05, print=TRUE, plot=TRUE) 
# 
# res <- RESTREND(stdRESTREND$max.NDVI, stdRESTREND$acc.precip, stdRESTREND$index, sig=0.05, print=TRUE, plot=TRUE) 
# 
