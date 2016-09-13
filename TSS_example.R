this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
#
# # load("./demo_data/stdRESTREND.Rda")
# # load("./demo_data/stdRESTREND_CTSR.Rda")
# # load("./demo_data/stdRESTREND_RF.Rda")
#
# # load("./demo_data/segRESTREND.Rda")
# # load("./demo_data/segRESTREND_CTSR.Rda")
# # load("./demo_data/segRESTREND_RF.Rda")
# # #
load("./demo_data/segVPRD.Rda")
load("./demo_data/segVPRD_CTSR.Rda")
load("./demo_data/segVPRD_RF.Rda")



# load("./demo_data/segVPRI.Rda")
# load("./demo_data/segVPRI_CTSR.Rda")
# load("./demo_data/segVPRI_RF.Rda")

# install.packages("devtools")

# library("devtools")

# install_github("ArdenB/SegmentedRESTREND_pub", subdir="TSS.RESTREND")
# library(TSS.RESTREND)


# start.time <- Sys.time()
#
#
# source("TSS_RESTREND.R")

# Set the Complete VI time series


# CTSR.VI <- stdRESTREND.CTSR$cts.NDVI
# CTSR.VI <- segRESTREND.CTSR$cts.NDVI
CTSR.VI <- segVPRD.CTSR$cts.NDVI
# CTSR.VI <- segVPRI.CTSR$cts.NDVI

#Set the matching rainfall time series

# rf.data <- stdRES.CTSRrf.TS$precip
# rf.data <- segRES.CTSRrf.TS$precip
rf.data <- segVPRD.CTSRrf.TS$precip
# rf.data <- segVPRI.CTSRrf.TS$precip

# Define the max accumuulation period
max.acp <- 12

#Define the max offset period
max.osp <- 4

#Create a table of every possible precipitation value given the max.acp and max.osp
ACP.table <- rainfall.accumulator(CTSR.VI, rf.data, max.acp, max.osp)

# check the ACP.table is the right shape
print(dim(ACP.table))

#Pass the ACP.table and the CTSR.VI to the TSS.RESTREND
results <- TSS.RESTREND(CTSR.VI, ACP.table, print=TRUE, plot=TRUE)



# end.time <- Sys.time()
# time.taken <- end.time - start.time
# print(time.taken)
