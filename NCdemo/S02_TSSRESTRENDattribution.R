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


tssr.attr <- function(line, VI, PP, TM, max.acp, max.osp){
  # =========== Function is applied to one pixel at a time ===========
  # ========== Perfrom the data checks and if anything fails skipp processing ==========
  # There is a data check for NANs in the TSSRattribution function, If SkipError is True
  # It then returns an opject of the same structure as actual results but filled with NaN
  # Usefull stacking using the foreac::do command.
  if (any(is.na(VI))){
    results = TSSRattribution(c(NA, NA), c(NA, NA), c(NA, NA), max.acp, max.osp, SkipError=TRUE)
  }else if (any(is.na(PP))){
    results = TSSRattribution(c(1,1), c(NA, NA), c(NA, NA), max.acp, max.osp, SkipError=TRUE)
  }else if (any(is.na(TM))){
    results = TSSRattribution(c(1,1), c(1,1), c(NA, NA), max.acp, max.osp, SkipError=TRUE)
  }else{

  	# +++++ Nothing has failed +++++
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

  	# ========== get the results ==========
  	results = TSSRattribution(CTSR.VI, CTSR.RF, CTSR.TM, max.acp, max.osp)
  }
  # ========== return the results ==========
  ret <- results$summary
  rownames(ret) <- line          # add the row name back
  ret["errors"] = results$errors # add the reason fro failure
  return(ret)
}

# ========== Calculate the number of rows i need to loop over ==========
max.acp <- 12
max.osp <- 4

ptime  <- system.time(
  tss.atdf <- foreach(
    line=rownames(VIdf), .combine = rbind) %do% {
      tssr.attr(line, VIdf[line, ], PPdf[line,], TMdf[line,], max.acp, max.osp)
    })

# ========== name to save the file ==========
fnout <- "./AttributionResults.csv"
write.csv(tss.atdf, fnout)

# browser()

