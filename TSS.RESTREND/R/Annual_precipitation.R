#' @title Antecedental Accumulation calculator for the Annual Max VI Time Series
#'
#' @description
#' Takes the Annual Max VI Time Series, the VI.index and a table of every possible accumulation period and offset
#' period.  A OLS is calculated \code{\link{lm}} for every combination of VI ~ rainfall.  This Function preferences those results where
#' slope>0 (increase in rainfall causes an increase in vegetation), returning the rainfall accumulation that has the highest R-squared
#' and a positive slope. If no combinations produce a positive slope then the one with the highest Rsquared is returned.MISSING non
#' peramtric and other variables
#'
#' @param anu.VI
#' The annual (Growing season) max VI. Which can be calculated from the CTSR.VI using \code{\link{AnMax.VI}}.
#' @param VI.index
#' the index of the CTSR.VI ts that the anu.VI values occur at. NOTE. R indexs from 1 rather than 0.
#' NOTE. if VI.index=FALSE, It can be calculated from the CTSR.VI using \code{\link{AnMax.VI}}.
#' @param ACP.table
#' A table of every combination of offset period and accumulation period. It can be calculated from ________
#'
#' @return _____(acu.RF)
#' The optimal accumulated rainfall for anu.VI. Mut be a object of class 'ts' and of equal length to anu.VI.
#' if anu.RF=FALSE, it will be calculated from ACP.table. see(____)
#' @return summary
#' (To be filled in)
#' @export
#'
AnnualRF.Cal <- function(anu.VI, VI.index, ACP.table, Breakpoint = FALSE, allow.negative=FALSE){
  if (class(anu.VI) != "ts")
    stop("anu.VI Not a time series object")
  if (length(VI.index) != length(anu.VI))
    stop("acu.VI and VI.index are not of the same length")


  yst <- start(anu.VI)[1]
  mst <- start(anu.VI)[2]

  anu.ACUP <- ACP.table[, VI.index]
  lines <- dim(anu.ACUP)[1]
  len <- dim(anu.ACUP)[2]

  m<- matrix(nrow=(lines), ncol=4)
  rownames(m)<- rownames(ACP.table)
  colnames(m)<- c("slope", "intercept", "p-value", "R^2.Value")

  p <- matrix(nrow=(lines), ncol=4)
  rownames(p)<- rownames(ACP.table)
  colnames(p)<- c("slope", "intercept", "p-value", "R^2.Value")

  # n <- 1
  for (n in 1:lines){
    # print(n)
    if (!Breakpoint){
      fit <- lm(anu.VI ~ anu.ACUP[n, ])
      # fit <- lm(anu.VI ~ segVPRD$acc.precip)
      R.pval <- glance(fit)$p.value
      R.Rval <- summary(fit)$r.square
      R.intr <- as.numeric(coef(fit)[1])
      R.slpe <- as.numeric(coef(fit)[2])
      m[n, ] <- c(R.slpe, R.intr,R.pval, R.Rval)
    }else{
      # print("Breakpoint")
      fit <- lm(anu.VI[1:Breakpoint] ~ anu.ACUP[n, 1:Breakpoint])
      R.pval <- glance(fit)$p.value
      R.Rval <- summary(fit)$r.square
      R.intr <- as.numeric(coef(fit)[1])
      R.slpe <- as.numeric(coef(fit)[2])
      m[n, ] <- c(R.slpe, R.intr,R.pval, R.Rval)

      fit2 <- lm(anu.VI[Breakpoint:len] ~ anu.ACUP[n, Breakpoint:len])
      R.pval2 <- glance(fit2)$p.value
      R.Rval2 <- summary(fit2)$r.square
      R.intr2 <- as.numeric(coef(fit2)[1])
      R.slpe2 <- as.numeric(coef(fit2)[2])
      p[n, ] <- c(R.slpe2, R.intr2,R.pval2, R.Rval2)
    }
  }
  if (!Breakpoint){
    if (allow.negative){

      max.line <- which.max(m[, "R^2.Value"])
      suma <- m[max.line,]
      anu.ARF <- ts(anu.ACUP[max.line, ], start=c(yst, mst), frequency = 1)
      return(structure(list(summary=suma, annual.precip = anu.ARF)))
    }else{

      mx <- m[m[, "slope"] > 0,]
      if (dim(mx)[1] == 0){
        warning("No positve slopes exist. Returing most significant negative slope")
        max.line <- which.max(m[, "R^2.Value"])
        suma <- m[max.line,]
        anu.ARF <- ts(anu.ACUP[max.line, ], start=c(yst, mst), frequency = 1)
        return(structure(list(summary=suma, annual.precip = anu.ARF)))
      }else{
        rfx <- anu.ACUP[m[, "slope"] > 0,]
        max.line <- which.max(mx[, "R^2.Value"])
        suma <- mx[max.line,]
        anu.ARF <- ts(rfx[max.line, ], start=c(yst, mst), frequency = 1)
        return(structure(list(summary=suma, annual.precip = anu.ARF)))
      }
    }
  }else{
    if (allow.negative){
      max.line <- which.max(m[, "R^2.Value"])
      suma <- m[max.line,]
      anu.ARF <- ts(anu.ACUP[max.line, ], start=yst, frequency = 1)

      pmax.line <- which.max(p[, "R^2.Value"])
      p.suma <- p[pmax.line,]
      panu.ARF <- ts(anu.ACUP[pmax.line, ], start=yst, frequency = 1)
      summ <- c(suma, p.suma)
      return(structure(list(summary=suma, rf.b4 = anu.ARF, rf.af= panu.ARF)))
    }else{
      mx <- m[m[, "slope"] >= 0,]
      if (dim(mx)[1] == 0){
        warning("No positve slopes exist before the bp. Returing most significant negative slope")
        max.line <- which.max(m[, "R^2.Value"])
        suma <- m[max.line,]
        anu.ARF <- ts(anu.ACUP[max.line, ], start=c(yst, mst), frequency = 1)
      }else{
        rfx <- anu.ACUP[m[, "slope"] > 0,]
        max.line <- which.max(mx[, "R^2.Value"])
        suma <- mx[max.line,]
        anu.ARF <- ts(rfx[max.line, ], start=c(yst, mst), frequency = 1)
      }
      px <- p[p[, "slope"] >= 0,]
      if (dim(px)[1] == 0){
        warning("No positve slopes exist after the bp. Returing most significant negative slope")
        pmax.line <- which.max(p[, "R^2.Value"])
        p.suma <- p[pmax.line,]
        panu.ARF <- ts(anu.ACUP[pmax.line, ], start=c(yst, mst), frequency = 1)
      }else{
        p.rfx <- anu.ACUP[p[, "slope"] > 0,]
        pmax.line <- which.max(px[, "R^2.Value"])
        p.suma <- px[pmax.line,]
        panu.ARF <- ts(p.rfx[pmax.line, ], start=c(yst, mst), frequency = 1)
      }
      summ <- c(suma, p.suma)
      return(structure(list(summary=suma, rf.b4 = anu.ARF, rf.af= panu.ARF)))




    }

  }

}

