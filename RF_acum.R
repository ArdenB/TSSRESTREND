


RF.accumulator <- function(CTSR.VI, rf.data, max.acp, max.osp){
  
  if (class(CTSR.VI) != "ts") 
    stop("CTSR.VI Not a time series object")
  if (class(rf.data) != "ts") 
    stop("rf.data Not a time series object")
  yst <- start(CTSR.VI)[1]
  mst <-  start(CTSR.VI)[2] 
  y.en <- end(CTSR.VI)[1]
  m.en <- end(CTSR.VI)[2]
  rf.yend <- end(rf.data)[1]
  rf.mend <-  end(rf.data)[2]
  if ((y.en != rf.yend) || rf.mend != m.en)
    stop("rf.data does not end at the same time as CTSR.VI")
  if (length(rf.data)<(length(CTSR.VI)+max.acp+max.osp))
    stop("rf.data is not long enough for the set max.acp and max.ops")
  

  row.nm <- rep(0, max.acp)
  if (max.osp>1){
    for (n in 1:(max.osp-1)){
      row.nm <- c(row.nm, rep(n, max.acp))
    }
  }
  
  len <- length(CTSR.VI)
  
  #Set up a blank matrix to write into 
  m<- matrix(nrow=(max.acp*max.osp), ncol=len)

  rownames(m)<- paste(row.nm, rep(1:max.acp, max.osp), sep = "-")
  colnames(m)<- c(1:len)
  # index <- 1
  
  m2 <- matrix(nrow=(max.acp), ncol=length(rf.data))
  rev.rf = rev(rf.data)
  for (n in 1:max.acp){
    roll= (roll_sum(rev.rf, n))
    # print(length(roll))
    if (n>1){
      roll = c(roll, rep(NaN, (n-1)))
    }
    m2[n, ] = roll
  }
  

  #turn the suplied ts in a table
  for (osp in 0:(max.osp-1)){
    
    m3 <- matrix(m2[,(1+osp):(len+osp)])
    m4 <- matrix(m3[, ncol(m3):1])
    # m4 <- apply(m3, 3, rev)
    ind <- 1 + (osp*max.acp)
    # print(ind)
    # print(dim(m4))
    m[ind:(ind+max.acp-1),] <- m4

  }
  return(m)
}
# 
# setwd("/mnt/FCBE3028BE2FD9C2/Users/user/Documents/segres_demo")
# 
load("./demo_data/segVPRD_RF.Rda")
load("./demo_data/segVPRD_CTSR.Rda")
CTSR.VI <- segVPRD.CTSR$cts.NDVI
rf.data <- CTSRrf.TS$precip
max.acp <- 12
max.osp <- 4
acum.table <-RF.accumulator(CTSR.VI, rf.data, max.acp, max.osp )
dim(res)
# 
