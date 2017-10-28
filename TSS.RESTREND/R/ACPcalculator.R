#' @title Antecedental Rainfall Accumulation calculator for the VI Complete Time Series
#'
#' @description
#' Takes the Complete Time Series Vegetation Index and a table of every possible accumulation period and offset
#' period.  A OLS is calculated \code{\link{lm}} for every combination of VI ~ rainfall.  This Function preferences those results where
#' slope>0 (increase in rainfall causes an increase in vegetation), returning the rainfall accumulation that has the highest R-squared
#' and a positive slope. If no combinations produce a positive slope then the one with the highest Rsquared is returned.
#'
#' @inheritParams TSSRESTREND
#' @return A list containing:
#' @return \bold{summary}
#'        a Matrix containing "slope", "intercept", "p.value", "R^2.Value", "Break.Height", "Slope.Change"
#'        of the \code{\link[stats]{lm}} between Antecedental Rainfall Accumulation (CTSR.RF) and the CTSR.VI
#' @return \bold{CTSR.precip}
#'        see CTSR.RF in \code{\link{TSSRESTREND}} for description
#' @return \bold{CTSR.osp}
#'        The offest period for the complete time series rainfall
#' @return \bold{CTSR.acp}
#'        The accumulation period for the complete time series rainfall
#' @export
#' @examples
#' #Find the data
#' vi.path <- system.file("extdata", "rabbitVI.csv", package = "TSS.RESTREND", mustWork = TRUE)
#' in.VI <- read.csv(vi.path)
#' CTSR.VI <- ts(in.VI, start=c(1982, 1), end=c(2013,12), frequency = 12)
#' data(rabbitACPtable)
#' ACPres <- ACP.calculator(CTSR.VI, rabbitACPtable)
#' print(ACPres$summary)

ACP.calculator <- function(CTSR.VI, ACP.table, allow.negative=FALSE){
  #Check classes of input data
  if (class(CTSR.VI) != "ts")
    stop("CTSR.VI Not a time series object")
  if (length(CTSR.VI) != dim(ACP.table)[2])
    stop("ACP.table size does not match CTSR.VI")

  #Get the start year and start month of the data
  yst <- start(CTSR.VI)[1]
  mst <-  start(CTSR.VI)[2]

  #Get the dimensions of the data for indexing
  len <- dim(ACP.table)[2]
  lines <- dim(ACP.table)[1]
  #empy matrix to store lm results
  m<- matrix(nrow=(lines), ncol=6)

  #Pull out the row names
  rownames(m)<- rownames(ACP.table)
  colnames(m)<- c("slope", "intercept", "p.value", "R^2.Value", "Break.Height", "Slope.Change")

  # fit a LM to every combination of rainfall and vegetation
  for (n in 1:lines){
    #get the lm
    fit <- lm(CTSR.VI ~ ACP.table[n, ])
    #Get the key values [pval, Rsquared, intercept, slope]
    R.pval <- glance(fit)$p.value
    R.Rval <- summary(fit)$r.square
    R.intr <- as.numeric(coef(fit)[1])
    R.slpe <- as.numeric(coef(fit)[2])
    R.BH <- NaN
    R.SC <- NaN
    # Stack the results in the empyt matryx m
    m[n, ] <- c(R.slpe, R.intr,R.pval, R.Rval, R.BH, R.SC)


  }
  # Look for the number of combinations the had positive slopes
  mx <- matrix(m[m[, "slope"] > 0,], ncol=6)
  colnames(mx)<- c("slope", "intercept", "p.value", "R^2.Value", "Break.Height", "Slope.Change")

  if (dim(mx)[1] <= 1||allow.negative){
    # warning("<2 CTS positve slopes exist. Returing most significant negative slope")
    # Get the max line and values
    max.line <- which.max(m[, "R^2.Value"])
    suma <- m[max.line,]
    # Get the values to return
    CTSR.ARF <- ts(ACP.table[max.line, ], start=c(yst, mst), frequency = 12)
    namestr <- rownames(ACP.table)[max.line]
    nmsplit <- strsplit(namestr, "\\-")[[1]]
    osp <- as.numeric(nmsplit[1])
    acp <- as.numeric(nmsplit[2])
    return(structure(list(summary=suma, CTSR.precip = CTSR.ARF, CTSR.osp = osp, CTSR.acp = acp)))
  }else{
    #Filter the negative results out
    rfx <- ACP.table[m[, "slope"] > 0,]
    # Get the max line and values
    max.line <- which.max(mx[, "R^2.Value"])
    suma <- mx[max.line,]
    # get the values to return
    CTSR.ARF <- ts(rfx[max.line, ], start=c(yst, mst), frequency = 12)
    namestr <- rownames(rfx)[max.line]
    nmsplit <- strsplit(namestr, "\\-")[[1]]
    osp <- as.numeric(nmsplit[1])
    acp <- as.numeric(nmsplit[2])
    return(structure(list(summary=suma, CTSR.precip = CTSR.ARF, CTSR.osp = osp, CTSR.acp = acp)))
  }
}
