#' @title Rainfall Accumulator
#'
#' @description
#' Takes the time series of rainfall and returns a rainfall accumulation table of every possible combination
#' of the max accumulation period and the max offset period.
#' @importFrom RcppRoll roll_sum
#' @param CTSR.VI
#' Complete Time Series of Vegetation Index. An object of class \code{'ts'}. Monthly time series of VI values
#' @param rf.data
#' Complete Time Series of monthly rainfall. An object of class \code{'ts'}. Must have the same end as CTSR.VI
#' and be longer that the max acsumuation period and the max offset period plus the start CTSR.VI
#' @param max.acp
#' The max accumuation period. Must be an integer > 1.
#' @param max.osp
#' The max offset period. Must be an integer >1
#'
#' @return ACP.table
#' A matrix with ever possible accumuated rainfall combination
#'
#' @export
#'
#' @examples
#' # Define the max accumuulation period
#' acp <- 12
#' #Define the max offset period
#' osp <- 4
#' rftable <- rainfall.accumulator(segRESTRENDCTSR$cts.NDVI, segRESTRENDctRF$precip, acp, osp)

rainfall.accumulator <- function(CTSR.VI, rf.data, max.acp, max.osp){

  if (class(CTSR.VI) != "ts")
    stop("CTSR.VI Not a time series object")
  if (class(rf.data) != "ts")
    stop("rf.data Not a time series object")
  if (sd(rf.data)==0)
    stop("The precipitation data has identical values (SD=0)")
  # Get the start and end dates of the precip and the VI
  yst <- start(CTSR.VI)[1]
  mst <-  start(CTSR.VI)[2]
  y.en <- end(CTSR.VI)[1]
  m.en <- end(CTSR.VI)[2]
  rf.yend <- end(rf.data)[1]
  rf.mend <-  end(rf.data)[2]
  #Check to make sure they have no issues
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
    #defence against index 1 i missing
    m3 <- (m2[,(1+osp):(len+osp)])
    m4 <- (m3[, ncol(m3):1])
    ind <- 1 + (osp*max.acp)
    m[ind:(ind+max.acp-1),] <- m4

  }
  # Return the table
  return(m)
}
