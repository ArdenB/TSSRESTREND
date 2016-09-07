# antecedental Accumulation calculator

ACP_calculator <- function(CTSR.VI, acum.table){
  
  if (class(CTSR.VI) != "ts") 
    stop("CTSR.VI Not a time series object")
  if (length(CTSR.VI) != dim(acum.table)[2])
    stop("acum.table size does not match CTSR.VI")
  
  
  yst <- start(CTSR.VI)[1]
  mst <-  start(CTSR.VI)[2] 
  
  len <- dim(acum.table)[2]
  lines <- dim(acum.table)[1]
  m<- matrix(nrow=(lines), ncol=4)
  
  rownames(m)<- rownames(acum.table)
  colnames(m)<- c("slope", "intercept", "p-value", "R^2.Value")
  # n <- 1
  for (n in 1:lines){
    # print(n)
    fit <- lm(CTSR.VI ~ acum.table[n, ])
    R.pval <- glance(fit)$p.value
    R.Rval <- summary(fit)$r.square
    R.intr <- as.numeric(coef(fit)[1])
    R.slpe <- as.numeric(coef(fit)[2])
    m[n, ] <- c(R.slpe, R.intr,R.pval, R.Rval)
    
  }
  max.line <- which.max(m[, "R^2.Value"])
  suma <- m[max.line,]
  CTSR.ARF <- ts(acum.table[max.line, ], start=c(yst, mst), frequency = 12)
  
  return(structure(list(summary=suma, CTSR.precip = CTSR.ARF)))
}
