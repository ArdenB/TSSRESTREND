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
#' @export
#' @examples
#' #Find the data
#' vi.path <- system.file("extdata", "rabbitVI.csv", package = "TSS.RESTREND", mustWork = TRUE)
#' in.VI <- read.csv(vi.path)
#' CTSR.VI <- ts(in.VI, start=c(1982, 1), end=c(2013,12), frequency = 12)
#' data(rabbitACPtable)
#' ACPres <- ACP.calculator(CTSR.VI, rabbitACPtable)
#' print(ACPres$summary)

ACP.calculator <- function(CTSR.VI, ACP.table){

  if (class(CTSR.VI) != "ts")
    stop("CTSR.VI Not a time series object")
  if (length(CTSR.VI) != dim(ACP.table)[2])
    stop("ACP.table size does not match CTSR.VI")


  yst <- start(CTSR.VI)[1]
  mst <-  start(CTSR.VI)[2]

  len <- dim(ACP.table)[2]
  lines <- dim(ACP.table)[1]
  m<- matrix(nrow=(lines), ncol=6)

  rownames(m)<- rownames(ACP.table)
  colnames(m)<- c("slope", "intercept", "p.value", "R^2.Value", "Break.Height", "Slope.Change")

  for (n in 1:lines){
    fit <- lm(CTSR.VI ~ ACP.table[n, ])
    R.pval <- glance(fit)$p.value
    R.Rval <- summary(fit)$r.square
    R.intr <- as.numeric(coef(fit)[1])
    R.slpe <- as.numeric(coef(fit)[2])
    R.BH <- NaN
    R.SC <- NaN
    m[n, ] <- c(R.slpe, R.intr,R.pval, R.Rval, R.BH, R.SC)


  }
  mx <- matrix(m[m[, "slope"] > 0,], ncol=6)
  colnames(mx)<- c("slope", "intercept", "p.value", "R^2.Value", "Break.Height", "Slope.Change")
  #
  if (dim(mx)[1] == 0){
    # warning("No positve slopes exist. Returing most significant negative slope")
    max.line <- which.max(m[, "R^2.Value"])
    suma <- m[max.line,]
    CTSR.ARF <- ts(ACP.table[max.line, ], start=c(yst, mst), frequency = 12)
    # browser()
    return(structure(list(summary=suma, CTSR.precip = CTSR.ARF)))
  }else{
    rfx <- ACP.table[m[, "slope"] > 0,]
    max.line <- which.max(mx[, "R^2.Value"])
    suma <- mx[max.line,]
    CTSR.ARF <- ts(rfx[max.line, ], start=c(yst, mst), frequency = 12)
    # browser()
    return(structure(list(summary=suma, CTSR.precip = CTSR.ARF)))
  }
}
