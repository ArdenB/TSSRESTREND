# library("bfast")
# #library("forecast")
# library("RcppCNPy")
# library("strucchange")
# library("broom")
# 
# setwd("/mnt/FCBE3028BE2FD9C2/Users/user/Documents/segres_demo") #needs to be replaced witha variable function
# 
# #load the data
# 
# load("./demo_data/stdRESTREND.Rda")
# load("./demo_data/stdRESTREND_CTSR.Rda")
# load("./demo_data/segRESTREND.Rda")
# load("./demo_data/segRESTREND_CTSR.Rda")
# load("./demo_data/segVPRD.Rda")
# load("./demo_data/segVPRD_CTSR.Rda")
# load("./demo_data/segVPRI.Rda")
# load("./demo_data/segVPRI_CTSR.Rda")

AnMax.VI <- function(CTSR.VI){
  
  sty <- start(CTSR.VI)[1]
  stm <-  start(CTSR.VI)[2] 
  eny <- end(CTSR.VI)[1]
  enm <- end(CTSR.VI)[2]
  
  m<- matrix(nrow=(eny-sty+1), ncol=12)
  rownames(m)<- c(sty:eny)
  colnames(m)<- c(month.abb[1:12])
  index <- 1
  for (yr in sty:eny){
    # m[toString(yr), 13:14] = 0
    # print(yr)  
    for (mon in 1:12){
      if (yr==sty & mon<stm){
        m[toString(yr), month.abb[mon]] = NaN
      }else if (yr==eny & mon>enm){
        m[toString(yr), month.abb[mon]] = NaN
      }else{
        m[toString(yr), month.abb[mon]] = CTSR.VI[index]
      }
      index <- (index+1)
    }
  }
  anmax.ts <- ts(apply(m,1,max, na.rm=TRUE), start=sty, frequency = 1)
  whmax.ts <- ts(apply(m,1,which.max), start=sty, frequency = 1)
  index.ts <- ts(c(sty:eny), start=sty, frequency = 1)
  df <- data.frame( Max=anmax.ts, Max.month=whmax.ts, index = index.ts)
  for (row in 1:(dim(df)[1]-1)){
    if ((df$Max.month[row] == 11 || df$Max.month[row] == 12) && (df$Max.month[row+1]==1 || df$Max.month[row+1]==2)){
      print(df[row,])
      print(row)
      test <- c(m[row, 11:12], m[row+1, 1:2])
      # print(test)
      df$Max.month[row] = which.max(test)+10
      df$Max[row] = max(test, na.rm=TRUE)
      nxtyr <- m[row+2, 2:12]
      df$Max.month[row+1] = which.max(nxtyr)+2
      df$Max[row+1] = max(nxtyr)
    }
  }
  df$index <- ((df$index-sty)*12+df$Max.month)
  #need to test this
  if (stm != 1){
    df$index <- (df$index-(stm-1))
    print("this feature requires further testing")
  }
  return(df)
}
CTSR.VI <- segVPRD.CTSR$cts.NDVI
t1 <-AnMax.VI(CTSR.VI)
