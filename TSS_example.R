this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
#
load("./demo_data/stdRESTREND.Rda")
load("./demo_data/stdRESTREND_CTSR.Rda")
load("./demo_data/stdRESTREND_RF.Rda")
CTSR.VI <- stdRESTREND.CTSR$cts.NDVI
rf.data <- stdRES.CTSRrf.TS$precip
#
# load("./demo_data/segRESTREND.Rda")
# load("./demo_data/segRESTREND_CTSR.Rda")
# load("./demo_data/segRESTREND_RF.Rda")
# CTSR.VI <- segRESTREND.CTSR$cts.NDVI
# rf.data <- segRES.CTSRrf.TS$precip
# # #
# load("./demo_data/segVPRD.Rda")
# load("./demo_data/segVPRD_CTSR.Rda")
# load("./demo_data/segVPRD_RF.Rda")
# CTSR.VI <- segVPRD.CTSR$cts.NDVI
# rf.data <- segVPRD.CTSRrf.TS$precip



# load("./demo_data/segVPRI.Rda")
# load("./demo_data/segVPRI_CTSR.Rda")
# load("./demo_data/segVPRI_RF.Rda")
# CTSR.VI <- segVPRI.CTSR$cts.NDVI
# rf.data <- segVPRI.CTSRrf.TS$precip

# install.packages("devtools")

# library("devtools")

# install_github("ArdenB/SegmentedRESTREND_pub", subdir="TSS.RESTREND")
# library(TSS.RESTREND)


# start.time <- Sys.time()
#
#
# source("TSS_RESTREND.R")

# Set the Complete VI time series



#Set the matching rainfall time series


# Define the max accumuulation period
max.acp <- 12

#Define the max offset period
max.osp <- 4

#Create a table of every possible precipitation value given the max.acp and max.osp
ACP.table <- rainfall.accumulator(CTSR.VI, rf.data, max.acp, max.osp)

# check the ACP.table is the right shape
print(dim(ACP.table))

#Pass the ACP.table and the CTSR.VI to the TSS.RESTREND
# exc = c(125:127)
results <- TSSRESTREND(CTSR.VI, ACP.table)
plot(results)
# stdRESTRENDctRF <- stdRES.CTSRrf.TS
# stdRESTRENDCTSR <- stdRESTREND.CTSR
# save(stdRESTREND, file = "./TSS.RESTREND/data/stdRESTREND.rda")
# save(stdRESTRENDCTSR, file = "./TSS.RESTREND/data/stdRESTREND_CTSR.rda")
# save(stdRESTRENDctRF, file = "./TSS.RESTREND/data/stdRESTREND_RF.rda")
# stdRESTRENDrfTab <- ACP.table
# save(stdRESTRENDrfTab, file = "./TSS.RESTREND/data/stRESTREND_ACPtable.rda")

# # segVPR <- data.frame(max.NDVI = results$ts.data$anu.VI, acum.RF=results$ts.data$acu.RF,
# #                      index=results$ts.data$VI.index, RFB4=ARCseg$rf.b4, RFAF= ARCseg$rf.af)
# # save(segVPR, file = "./TSS.RESTREND/data/segVPR.rda")
# segVPRctRF <- segVPRD.CTSRrf.TS
# save(segVPRctRF, file = "./TSS.RESTREND/data/segVPRctRF.rda")
# segVPRrfTab <- ACP.table
# save(segVPRrfTab, file = "./TSS.RESTREND/data/segVPR_ACPtable.rda")
# segVPRCTSR <- segVPRD.CTSR
# save(segVPRCTSR, file = "./TSS.RESTREND/data/segVPR_CTSR.rda")


#
# segRESTRENDctRF <- segRES.CTSRrf.TS
# segRESTRENDCTSR <- segRESTREND.CTSR
# save(segRESTREND, file = "./TSS.RESTREND/data/segRESTREND.rda")
# save(segRESTRENDCTSR, file = "./TSS.RESTREND/data/segRESTREND_CTSR.rda")
# save(segRESTRENDctRF, file = "./TSS.RESTREND/data/segRESTREND_RF.rda")
# segRESTRENDrfTab <- ACP.table
# save(segRESTRENDrfTab, file = "./TSS.RESTREND/data/segRESTREND_ACPtable.rda")



# # end.time <- Sys.time()
# time.taken <- end.time - start.time
# print(time.taken)
