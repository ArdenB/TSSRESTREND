# antecedental Accumulation calculator

ACP_calculator <- function(CTSR.VI, ACP.table){
  
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
  
#   max.line <- which.max(m[, "R^2.Value"])
#   suma <- m[max.line,]
#   CTSR.ARF <- ts(ACP.table[max.line, ], start=c(yst, mst), frequency = 12)
#   
#   return(structure(list(summary=suma, CTSR.precip = CTSR.ARF)))
}
