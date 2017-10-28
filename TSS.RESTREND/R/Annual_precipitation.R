#' @title Antecedental Accumulation calculator for the Annual Max VI Time Series
#'
#' @description
#' Takes the Annual Max VI Time Series, the VI.index and a table of every possible accumulation period and offset
#' period.  A OLS is calculated \code{\link[stats]{lm}} for every combination of VI ~ rainfall.  This Function preferences those results where
#' slope>0 (increase in rainfall causes an increase in vegetation), returning the rainfall accumulation that has the highest R-squared
#' and a positive slope. If no combinations produce a positive slope then the one with the highest Rsquared is returned.MISSING non
#' peramtric and other variables
#' @author Arden Burrell, arden.burrell@unsw.edu.au
#' @inheritParams TSSRESTREND
#' @param Breakpoint
#'        Used when calcualting rf.bf and rf.af for ts with breakpoints in the VPR.  See \code{\link{CHOW}}
#'
#' @return \bold{summary}
#'        a Matrix containing "slope", "intercept", "p.value", "R^2.Value", "Break.Height", "Slope.Change"
#'        of the \code{\link[stats]{lm}} of VI ~ rainfall. If Breakpoint, summary covers both rf.b4 and rf.af.
#' @return \bold{acu.RF}
#'        The optimal accumulated rainfall for anu.VI. Mut be a object of class 'ts' and of equal length to anu.VI. IT is caculated
#'        from the ACP.table by finding the acp and osp that has the largest R^2 value. \code{\link[stats]{lm}}(anu.VI ~ rainfall)
#' @return \bold{rf.b4}
#'        The optimal acumulated rainfall before the Breakpoint
#' @return \bold{rf.af}
#'        The Optimally accumulated rainfall after the Breakpoint
#' @return \bold{osp}
#'        The offest period for the annual max time series rainfall
#' @return \bold{acp}
#'        The accumulation period for the annual max time series rainfall
#' @export
#'
#' @examples
#' ARC <- AnnualRF.Cal(stdRESTREND$max.NDVI, stdRESTREND$index, stdRESTRENDrfTab)
#' print(ARC)
#' \dontrun{
#'
#' #Test the complete time series for breakpoints
#' VPRBFdem <- VPR.BFAST(segVPRCTSR$cts.NDVI, segVPRCTSR$cts.precip)
#' bp<-as.numeric(VPRBFdem$bkps)
#'
#' #test the significance of the breakpoints
#' reschow <- CHOW(segVPR$max.NDVI, segVPR$acum.RF, segVPR$index, bp)
#' brkp <- as.integer(reschow$bp.summary["yr.index"])
#' ARCseg <-AnnualRF.Cal(segVPR$max.NDVI, segVPR$index, segVPRrfTab, Breakpoint = brkp)
#' print(ARCseg)
#' }
AnnualRF.Cal <- function(anu.VI, VI.index, ACP.table, Breakpoint = FALSE, allow.negative=FALSE){
  #Perform the basic sanity checks
  if (class(anu.VI) != "ts")
    stop("anu.VI Not a time series object")
  if (length(VI.index) != length(anu.VI))
    stop("acu.VI and VI.index are not of the same length")

  #Get the start year and month
  yst <- start(anu.VI)[1]
  mst <- start(anu.VI)[2]

  # Subset the complete ACP.table using the indexs of the an mav values and get dims
  anu.ACUP <- ACP.table[, VI.index]
  lines <- dim(anu.ACUP)[1]
  len <- dim(anu.ACUP)[2]

  #create a blank matrix and assign it colum headings
  m<- matrix(nrow=(lines), ncol=6)
  rownames(m)<- rownames(ACP.table)
  colnames(m)<- c("slope", "intercept", "p.value", "R^2.Value", "Break.Height", "Slope.Change")
  # Second matrix in case of breakpoint in the VPR
  p <- matrix(nrow=(lines), ncol=6)
  rownames(p)<- rownames(ACP.table)
  colnames(p)<- c("slope", "intercept", "p.value", "R^2.Value", "Break.Height", "Slope.Change")

  # For each line of the accumulation table calculate the coefficents from a lm
  for (n in 1:lines){
    if (!Breakpoint){ #No Breakpoint
      if (sd(anu.ACUP[n, ])== 0 ){
        #if a combination of acp and osp leads to SD=0 rainfall, this will catch it
        # All values are bad
        m[n, ] <- c(-1, -1, 1, 0, NaN, NaN)
      }else{
        #Get the lm between VI and the precip on line n of the table
        fit <- lm(anu.VI ~ anu.ACUP[n, ])
        #Pull out key infomation
        R.pval <- glance(fit)$p.value
        R.Rval <- summary(fit)$r.square
        R.intr <- as.numeric(coef(fit)[1])
        R.slpe <- as.numeric(coef(fit)[2])
        R.BH <- NaN
        R.SC <- NaN
        #add the infomation the the m data table
        m[n, ] <- c(R.slpe, R.intr,R.pval, R.Rval, R.BH, R.SC)
      }

    }else{
      # Before the breakpoint
      if (sd(anu.ACUP[n, 1:Breakpoint]) ==0){
        #if a combination of acp and osp leads to SD=0 rainfall, this will catch it
        # All values are bad
        m[n, ] <- c(-1, -1, 1, 0, NaN, NaN)
      }else{
        #get the fit before the breakpoint
        fit <- lm(anu.VI[1:Breakpoint] ~ anu.ACUP[n, 1:Breakpoint])
        #Pull out key infomation before the bp
        R.pval <- glance(fit)$p.value
        R.Rval <- summary(fit)$r.square
        R.intr <- as.numeric(coef(fit)[1])
        R.slpe <- as.numeric(coef(fit)[2])
        R.BH <- NaN
        R.SC <- NaN
        m[n, ] <- c(R.slpe, R.intr,R.pval, R.Rval, R.BH, R.SC)
      }
      if (sd(anu.ACUP[n, (Breakpoint+1):len])==0) {
        #if a combination of acp and osp leads to SD=0 rainfall, this will catch it
        # All values are bad
        p[n, ] <- c(-1, -1, 1, 0, NaN, NaN)
      }else{
        #get the fit after the breakpoint
        fit2 <- lm(anu.VI[(Breakpoint+1):len] ~ anu.ACUP[n, (Breakpoint+1):len])
        #Pull out key infomation after the bp
        R.pval2 <- glance(fit2)$p.value
        R.Rval2 <- summary(fit2)$r.square
        R.intr2 <- as.numeric(coef(fit2)[1])
        R.slpe2 <- as.numeric(coef(fit2)[2])
        R.BH <- NaN
        R.SC <- NaN
        p[n, ] <- c(R.slpe2, R.intr2,R.pval2, R.Rval2, R.BH, R.SC)
      }
    }
  }

  # Figure out which acp&osp to use
  if (!Breakpoint){
    if (allow.negative){

      max.line <- which.max(m[, "R^2.Value"])
      #Catch failures
      if (max(m[, "R^2.Value"])==0){
        #Means that all the values after the breakpoint are the same
        warning("There is no variance (SD=0). Precipitation data failure")
      }
      suma <- m[max.line,]
      anu.ARF <- ts(anu.ACUP[max.line, ], start=c(yst, mst), frequency = 1)

      namestr <- rownames(m)[max.line]
      nmsplit <- strsplit(namestr, "\\-")[[1]]
      osp <- as.numeric(nmsplit[1])
      acp <- as.numeric(nmsplit[2])

      return(structure(list(summary=suma, annual.precip = anu.ARF, osp=osp, acp=acp)))

    }else{

      mx <- matrix(m[m[, "slope"] >= 0,],  ncol=6)
      colnames(mx) <- c("slope", "intercept", "p.value", "R^2.Value", "Break.Height", "Slope.Change")
      if (dim(mx)[1] <= 1){
        warning("Insufficent positve slopes exist. Returing most significant negative slope")
        max.line <- which.max(m[, "R^2.Value"])
        suma <- m[max.line,]
        anu.ARF <- ts(anu.ACUP[max.line, ], start=c(yst, mst), frequency = 1)

        namestr <- rownames(m)[max.line]
        nmsplit <- strsplit(namestr, "\\-")[[1]]
        osp <- as.numeric(nmsplit[1])
        acp <- as.numeric(nmsplit[2])

        return(structure(list(summary=suma, annual.precip = anu.ARF, osp=osp, acp=acp)))
      }else{
        rfx <- anu.ACUP[m[, "slope"] >= 0,]
        max.line <- which.max(mx[, "R^2.Value"])
        suma <- mx[max.line,]
        anu.ARF <- ts(rfx[max.line, ], start=c(yst, mst), frequency = 1)

        namestr <- rownames(rfx)[max.line]
        nmsplit <- strsplit(namestr, "\\-")[[1]]
        osp <- as.numeric(nmsplit[1])
        acp <- as.numeric(nmsplit[2])

        return(structure(list(summary=suma, annual.precip = anu.ARF, osp=osp, acp=acp)))
      }
    }
  }else{

    if (allow.negative){
      max.line <- which.max(m[, "R^2.Value"])
      suma <- m[max.line,]
      anu.ARF <- ts(anu.ACUP[max.line, ], start=yst, frequency = 1)

      namestr <- rownames(m)[max.line]
      # browser()
      nmsplit <- strsplit(namestr, "\\-")[[1]]
      osp.b4 <- as.numeric(nmsplit[1])
      acp.b4 <- as.numeric(nmsplit[2])

      pmax.line <- which.max(p[, "R^2.Value"])
      p.suma <- p[pmax.line,]
      panu.ARF <- ts(anu.ACUP[pmax.line, ], start=yst, frequency = 1)
      summ <- c(suma, p.suma)

      namestr.af <- rownames(p)[pmax.line]
      nmsplit.af <- strsplit(namestr.af, "\\-")[[1]]
      osp.af <- as.numeric(nmsplit[1])
      acp.af <- as.numeric(nmsplit[2])

      return(structure(list(summary=suma, rf.b4 = anu.ARF, rf.af= panu.ARF, osp.b4 = osp.b4, acp.b4 = acp.b4, osp.af = osp.af, acp.af = acp.af)))
    }else{
      mx <- matrix(m[m[, "slope"] >= 0,],  ncol=6)
      colnames(mx)<- c("slope", "intercept", "p.value", "R^2.Value", "Break.Height", "Slope.Change")
      if (dim(mx)[1] <= 1){
        warning("<2 positve slopes exist before the bp. Returing most significant negative slope")
        max.line <- which.max(m[, "R^2.Value"])
        suma <- m[max.line,]
        anu.ARF <- ts(anu.ACUP[max.line, ], start=c(yst, mst), frequency = 1)
        namestr <- rownames(m)[max.line]
        nmsplit <- strsplit(namestr, "\\-")[[1]]
        osp.b4 <- as.numeric(nmsplit[1])
        acp.b4 <- as.numeric(nmsplit[2])
      }else{
        rfx <- anu.ACUP[m[, "slope"] > 0,]
        max.line <- which.max(mx[, "R^2.Value"])
        #Catch failures
        if (max(m[, "R^2.Value"])==0){
          #Means that all the values after the breakpoint are the same
          warning("There is no variance (SD=0) before the breakpoint. Precipitation data failure")
        }

        suma <- mx[max.line,]
        anu.ARF <- ts(rfx[max.line, ], start=c(yst, mst), frequency = 1)
        namestr <- rownames(rfx)[max.line]
        nmsplit <- strsplit(namestr, "\\-")[[1]]
        osp.b4 <- as.numeric(nmsplit[1])
        acp.b4 <- as.numeric(nmsplit[2])
      }

      px <- matrix(p[p[, "slope"] > 0,],  ncol=6)
      colnames(px)<- c("slope", "intercept", "p.value", "R^2.Value", "Break.Height", "Slope.Change")
      if (dim(px)[1] <= 1){
        warning("<2 positve slopes exist after the bp. Returing most significant negative slope")
        pmax.line <- which.max(p[, "R^2.Value"])
        # catch failures past the breakpoint
        if (max(p[, "R^2.Value"])==0){
          #Means that all the values after the breakpoint are the same
          warning("There is no variance (SD=0) after the breakpoint. Precipitation data failure")
        }
        p.suma <- p[pmax.line,]
        panu.ARF <- ts(anu.ACUP[pmax.line, ], start=c(yst, mst), frequency = 1)
        namestr.af <- rownames(p)[pmax.line]
        nmsplit.af <- strsplit(namestr.af, "\\-")[[1]]
        osp.af <- as.numeric(nmsplit[1])
        acp.af <- as.numeric(nmsplit[2])
      }else{
        p.rfx <- anu.ACUP[p[, "slope"] > 0,]
        pmax.line <- which.max(px[, "R^2.Value"])
        p.suma <- px[pmax.line,]
        panu.ARF <- ts(p.rfx[pmax.line, ], start=c(yst, mst), frequency = 1)
        namestr.af <- rownames(p.rfx)[pmax.line]
        nmsplit.af <- strsplit(namestr.af, "\\-")[[1]]
        osp.af <- as.numeric(nmsplit[1])
        acp.af <- as.numeric(nmsplit[2])
      }
      summ <- c(suma, p.suma)
      return(structure(list(summary=suma, rf.b4 = anu.ARF, rf.af= panu.ARF, osp.b4 = osp.b4,
                            acp.b4 = acp.b4, osp.af = osp.af, acp.af = acp.af)))
      # return(structure(list(summary=suma, rf.b4 = anu.ARF, rf.af= panu.ARF)))

    }

  }

}

