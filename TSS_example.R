


load("./demo_data/segVPRD_RF.Rda")
load("./demo_data/segVPRD_CTSR.Rda")
CTSR.VI <- segVPRD.CTSR$cts.NDVI
rf.data <- CTSRrf.TS$precip
max.acp <- 12
max.osp <- 4
acum.table <-RF.accumulator(CTSR.VI, rf.data, max.acp, max.osp )
dim(acum.table)