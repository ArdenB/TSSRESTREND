# Title: TSS-RESTREND change detection and attribution demonstration
# Author: Arden Burrell

library(foreach)
library(xts)
library(zoo)
library(lubridate)
library(TSS.RESTREND)

# ========== read in the vegetation and climate data ==========
fnVI <- "./demo_dataframe_ndvi.csv"
VIdf <- read.csv(fnVI, row.names=1, check.names=FALSE) # NDVI
fnPP <- "./demo_dataframe_ppt.csv"
PPdf <- read.csv(fnPP, row.names=1, check.names=FALSE) # precip
fnTM <- "./demo_dataframe_tmean.csv"
TMdf <- read.csv(fnTM, row.names=1, check.names=FALSE) # temperature


tssr.attr <- function(line, VI, PP, TM){
  # =========== Function is applied to one pixel at a time ===========
  # ========== Perfrom the data check ==========
  if (any(is.na(VI))){
    return(NA)
  }else if (any(is.na(PP))){
    return(NA)
  }else if (any(is.na(TM))){
    return(NA)
  }
  print(line)
  # ========== Deal with Dates ==========
  # +++++ COnvert the VI to a TS object +++++
  VIdates <- as.POSIXlt(colnames(VI))
  VIys    <- VIdates[1]$year + 1900
  VIms    <- month(VIdates[1])
  VIyf    <- tail(VIdates, n=1)$year+ 1900
  VImf    <- month(tail(VIdates, n=1))
  CTSR.VI <- ts(as.numeric(VI), start=c(VIys, VIms), end=c(VIyf,VImf), frequency = 12)

  # +++++ COnvert the PPT to a TS object +++++
  PPdates <- as.POSIXlt(colnames(PP))
  PPys    <- PPdates[1]$year + 1900
  PPms    <- month(PPdates[1])
  PPyf    <- tail(PPdates, n=1)$year+ 1900
  PPmf    <- month(tail(PPdates, n=1))
  CTSR.RF <- ts(as.numeric(PP), start=c(PPys, PPms), end=c(PPyf,PPmf), frequency = 12)

  # +++++ COnvert the tMANS to a TS object +++++
  TMdates <- as.POSIXlt(colnames(TM))
  TMys    <- TMdates[1]$year + 1900
  TMms    <- month(TMdates[1])
  TMyf    <- tail(TMdates, n=1)$year+ 1900
  TMmf    <- month(tail(TMdates, n=1))
  CTSR.TM <- ts(as.numeric(TM), start=c(TMys, TMms), end=c(TMyf,TMmf), frequency = 12)


  TSSRattribution(CTSR.VI, CTSR.RF, CTSR.TM)
  # browser()
}

# ========== Calculate the number of rows i need to loop over ==========
# nline  <- dim(VIdf)[1]

# line = rownames(VIdf)[1]

ptime  <- system.time(
  tss.atdf <- foreach(
    line=rownames(VIdf), .combine = rbind) %do% {
      tssr.attr(line, VIdf[line, ], PPdf[line,], TMdf[line,])
    })

