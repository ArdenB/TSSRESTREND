#' @title Annual max VI Calculator
#'
#' @description
#' Takes the montly time series of the VI and calculates the growing season max VI. In series where the peak occurs in November or December,
#' an interannual growing season is assessed.
#'@param CTSR.VI
#' Complete Time Series of Vegetation Index. An object of class \code{'ts'}. Monthly time series of VI values
#' @return Max(anu.VI)
#' The annual (Growing season) max VI.
#' @return Max.Month
#' a Vector that has the month number where max values wew observed
#' @return index(VI.index)
#' the index of the CTSR.VI ts that the anu.VI values occur at.R indexs from 1 rather than 0.
#' @export
AnMax.VI <- function(CTSR.VI){
  if (class(CTSR.VI) != "ts")
      stop("CTSR.VI Not a time series object")

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
      # print(df[row,])
      # print(row)
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
# CTSR.VI <- segVPRD.CTSR$cts.NDVI
# t1 <-AnMax.VI(CTSR.VI)
