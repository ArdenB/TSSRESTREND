


load("./demo_data/segVPRD_RF.Rda")
load("./demo_data/segVPRD_CTSR.Rda")
source("TSS_RESTREND.R")

# Set the Complete VI time series 
CTSR.VI <- segVPRD.CTSR$cts.NDVI

#Set the matching rainfall time series
rf.data <- CTSRrf.TS$precip

# Define the max accumuulation period
max.acp <- 12

#Define the max offset period
max.osp <- 4

#Create a table of every possible precipitation value given the max.acp and max.osp
ACP.table <- rainfall.accumulator(CTSR.VI, rf.data, max.acp, max.osp)

# check the ACP.table is the right shape
print(dim(ACP.table))

#Pass the ACP.table and the CTSR.VI to the TSS.RESTREND
TSS.RESTREND(CTSR.VI, ACP.table, print=TRUE, plot=TRUE)
