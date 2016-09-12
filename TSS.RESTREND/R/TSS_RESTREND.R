library("bfast")
#library("forecast")
#library("RcppCNPy")
library("strucchange")
library("broom")
library("RcppRoll")

this.dir <- dirname(parent.frame(2)$ofile)
former.dir <- getwd()
setwd(this.dir)

#load the test data

load("./demo_data/stdRESTREND.Rda")
load("./demo_data/stdRESTREND_CTSR.Rda")
load("./demo_data/stdRESTREND_RF.Rda")

load("./demo_data/segRESTREND.Rda")
load("./demo_data/segRESTREND_CTSR.Rda")
load("./demo_data/segRESTREND_RF.Rda")

load("./demo_data/segVPRD.Rda")
load("./demo_data/segVPRD_CTSR.Rda")
load("./demo_data/segVPRD_RF.Rda")

load("./demo_data/segVPRI.Rda")
load("./demo_data/segVPRI_CTSR.Rda")
load("./demo_data/segVPRI_RF.Rda")



source("RF_acum.R")
source("max_pos.R")
source("CTSR_acp.R")
source("Annual_precipitation.R")
source("VPR_BFAST.R")
source("CHOW_bptest.R")
source("viMAX_plot.R")
source("segVPR.R")
source("segRESTREND.R")
source("RESTREND.R")

#return to the past directory 
setwd(former.dir)


#FUnctions to call functions
  #Function to call the other functions
  #Missing function to find optimal accumulation of the precipitation
  #which will be a sperate script. If the method is segVPR, there 
  #is a need to recalculate precip on either side of the breakpoint
  #Until it is functional  b4 and after need to be passed into this
  #function.  rf.b4=FALSE, rf.af=FALSE, will be removed as soon as 
TSS.RESTREND <- function(CTSR.VI, ACP.table=FALSE, CTSR.RF=FALSE, anu.VI=FALSE, acu.RF=FALSE, VI.index=FALSE, rf.b4=FALSE, rf.af=FALSE, 
                         sig=0.05, print=FALSE, plot=TRUE, details=FALSE){

  while (TRUE){ #Test the variables 
    if (class(CTSR.VI) != "ts") 
      stop("CTSR.VI Not a time series object. Please check the data")
    
    if ((!ACP.table) && (!CTSR.RF || acu.RF))
      stop("Rainfall data invalid. ACP.table or (CTSR.RF & acu.RF")
    
    if ((!anu.VI)||(!VI.index)){
      max.df <- AnMax.VI(CTSR.VI)
      anu.VI <- max.df$Max
      VI.index <- max.df$index
    }else{
      if (class(anu.VI) != "ts") 
        stop("anu.VI Not a time series object")
    }
    if (!CTSR.RF){
      CTS.Str <- ACP_calculator(CTSR.VI, ACP.table)
      CTSR.RF <- CTS.Str$CTSR.precip
      details.CTS.VPR <- CTS.Str$summary
    }else{
      if (class(CTSR.RF) != "ts") 
        stop("CTSR.RF Not a time series object")
      #get the time data out
      start.ti <- time(CTSR.VI)
      freq <- frequency(CTSR.VI)
      #check the two ts object cover the same time period
      start.ti2 <- time(CTSR.RF)
      freq2 <- frequency(CTSR.RF)
      if (!identical(start.ti, startti2))
        stop("ts objects do not have the same time, (CTSR.VI & CTSR.RF)")
      if (!identical(f, f2))
        stop("ts objects do not have the same frequency, (CTSR.VI & CTSR.RF)")
    }
    # browser()
    if (!acu.RF){
      precip.df <- AnnualRF.Cal(anu.VI, VI.index, ACP.table)
      acu.RF <- precip.df$annual.precip
      details.acu.RF <- precip.df$summary
    }else{
      if (class(acu.RF) != "ts") 
        stop("acu.VI Not a time series object")
      st.ti <- time(anu.VI)
      st.f <- frequency(anu.VI)
      #check the two ts object cover the same time period
      st.ti2 <- time(acu.RF)
      st.f2 <- frequency(acu.RF)
      if (!identical(st.ti, st.ti2))
        stop("ts object do not have the same time, (acu.RF & anu.VI)")
      if (!identical(st.f, st.f2))
        stop("ts object do not have the same frequency, (acu.RF & anu.VI)")
    }
    if (class(rf.b4) != "logical"){
      if (length(rf.b4) != (length(rf.af)))
          stop("rf.b4 and rf.af are different shapes. They must be the same size and be th same lenths as acu.VI")
    }
    break
  }
  
  bkp = VPR.BFAST(CTSR.VI, CTSR.RF, print=print, plot=plot, details = details)
  bp <- bkp$bkps
  if (!bp){# no breakpoints detected by the BFAST
    test.Method = "RESTREND"
    if (plot){
      VImax.plot(anu.VI)
    }
  }else{
    bp<-as.numeric(bkp$bkps)
    res.chow <- CHOW(anu.VI, acu.RF, VI.index, bp, sig=sig, print=print)
    brkp <- as.integer(res.chow$bp.summary["yr.index"])
    test.Method = res.chow$n.Method
    if (plot){
      VImax.plot(anu.VI, brkp=brkp)
    }
  }
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #add NDVI plot
  # browser()
  
  if (test.Method == "RESTREND"){
    result <- RESTREND(anu.VI, acu.RF, VI.index, sig=sig, print=print, plot=plot) 
  }else if (test.Method == "seg.RESTREND"){
    breakpoint = as.integer(res.chow$bp.summary[2])
    result <- seg.RESTREND(anu.VI, acu.RF, VI.index, brkp,  sig=sig, print=print, plot=plot)
  }else if (test.Method == "seg.VPR"){
    if ((!rf.b4)||(!rf.af)){
      #Improve the figure
      VPRbp.df <-AnnualRF.Cal(anu.VI, VI.index, ACP.table, Breakpoint = brkp)
      rf.b4 <- VPRbp.df$rf.b4
      rf.af <- VPRbp.df$rf.af
      #!!!!!!!!!!!!!!!!!!!!!!!!!!#
    }
    breakpoint = as.integer(res.chow$bp.summary[2])
    result <- seg.VPR(anu.VI, acu.RF, VI.index, brkp, rf.b4, rf.af, sig=sig, print=print, plot=plot)
  }
  # print(result$summary)
  # browser()
  return(result) #add CTSR 
}

