#' @title Antecedental Accumulation calculator for the VI Complete Time Series
#'
#' @description
#' Takes the Complete Time Series Vegetation Index and a table of every possible accumulation period and offset
#' period.  A OLS is calculated \code{\link{lm}} for every combination of VI ~ rainfall.  This Function preferences those results where
#' slope>0 (increase in rainfall causes an increase in vegetation), returning the rainfall accumulation that has the highest R-squared
#' and a positive slope. If no combinations produce a positive slope then the one with the highest Rsquared is returned.
#'
#' @param CTSR.VI
#' Complete Time Series of Vegetation Index. An object of class \code{'ts'}. Monthly time series of VI values
#' @param ACP.table
#' A table of every combination of offset period and accumulation period. if ACP.table = FALSE, CTSR.RF and acu.RF must be provided
#' @return summary
#' (To be filled in)
#' @return CTSR.precip (CTSR.RF)
#' Complete Time Series of the optimally accumulated rainfall.An object of class 'ts'.
#' @export
#'
ACP.calculator <- function(CTSR.VI, ACP.table){

  if (class(CTSR.VI) != "ts")
    stop("CTSR.VI Not a time series object")
  if (length(CTSR.VI) != dim(ACP.table)[2])
    stop("ACP.table size does not match CTSR.VI")


  yst <- start(CTSR.VI)[1]
  mst <-  start(CTSR.VI)[2]

  len <- dim(ACP.table)[2]
  lines <- dim(ACP.table)[1]
  m<- matrix(nrow=(lines), ncol=4)

  rownames(m)<- rownames(ACP.table)
  colnames(m)<- c("slope", "intercept", "p-value", "R^2.Value")

  for (n in 1:lines){
    fit <- lm(CTSR.VI ~ ACP.table[n, ])
    R.pval <- glance(fit)$p.value
    R.Rval <- summary(fit)$r.square
    R.intr <- as.numeric(coef(fit)[1])
    R.slpe <- as.numeric(coef(fit)[2])
    m[n, ] <- c(R.slpe, R.intr,R.pval, R.Rval)

  }
  mx <- m[m[, "slope"] > 0,]
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
