#' @title Antecedental Rainfall (and temperature) Accumulation calculator for the VI Complete Time Series
#'
#' @description
#' Takes the Complete Time Series Vegetation Index and a table of every possible accumulation period and offset
#' period for precipitation and temperature(optional).  A OLS is calculated \code{\link{lm}} for every combination
#' of VI ~ rainfall (and temperature if that is included).
#' if only the VPR is being calculated, this Function preferences those results where slope>0
#' (increase in rainfall causes an increase in vegetation),returning the rainfall accumulation
#' that has the highest R-squared and a positive slope. If no combinationsproduce a
#' positive slope then the one with the highest Rsquared is returned.
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
#' @return \bold{CTSR.tmp}
#'        The optimally accumulated CTS temperature
#' @return \bold{CTSR.tosp}
#'        The offest period for the complete time series temperature
#' @return \bold{CTSR.tacp}
#'        The accumulation period for the complete time series temperature
#' @export
#' @examples
#' #Find the data
#' vi.path <- system.file("extdata", "rabbitVI.csv", package = "TSS.RESTREND", mustWork = TRUE)
#' in.VI <- read.csv(vi.path)
#' CTSR.VI <- ts(in.VI, start=c(1982, 1), end=c(2013,12), frequency = 12)
#' data(rabbitACPtable)
#' ACPres <- ACP.calculator(CTSR.VI, rabbitACPtable)
#' print(ACPres$summary)

ACP.calculator <- function(CTSR.VI, ACP.table, ACT.table=NULL, allow.negative=FALSE, allowneg.retest=FALSE){
  #Check classes of input data
  if (class(CTSR.VI) != "ts")
    stop("CTSR.VI Not a time series object")
  if (length(CTSR.VI) != dim(ACP.table)[2])
    stop("ACP.table size does not match CTSR.VI")
  if ((!is.null(ACT.table)) && (dim(ACT.table)[2] != dim(ACP.table)[2]))
    stop("ACT.table size does not match CTSR.VI") #Check the temp table (ACT.table)

  #Get the start year and start month of the data
  yst <- start(CTSR.VI)[1]
  mst <-  start(CTSR.VI)[2]

    # Define a Function to perform linear regression and get out values
  linreg <- function(CTSR.VI, ACP.N, ACT.N = NULL, simple=TRUE){
    #get the lm
    # Simple tells me which values to return
    if (sd(ACP.table[n, ])== 0){
      #if a combination of acp and osp leads to SD=0 rainfall, this will catch it
      # All values are bad
      if (is.null(ACT.N)){return(c(0, -1))}else{return(0)}
      # return(c(-1, -1, 1, 0, NaN, NaN))
    }else{
      if (is.null(ACT.N)) { # No temperature data
        fit <- lm(CTSR.VI ~ ACP.N)
        R.Rval <- summary(fit)$r.square
        R.slpe <- as.numeric(coef(fit)[2])
        if (simple) {return(c(R.Rval, R.slpe))} #To speed up looping over all the pixels (simple = TRUE)
        R.pval <- glance(fit)$p.value
        R.intr <- as.numeric(coef(fit)[1])
        R.Tcoef <- NaN
        # pass a true value for the tempurature sig test
        t.sig <- TRUE

      }else{ # temperature data
      # do the multivariate regression with precip and temperature
        fit <- lm(CTSR.VI ~ ACP.N + ACT.N)
        R.Rval <- summary(fit)$r.square
        if (simple) {return(R.Rval)}#To speed up looping over all the pixels (simple = TRUE)
        R.pval <- glance(fit)$p.value
        R.intr <- as.numeric(coef(fit)[1])
        R.slpe <- as.numeric(coef(fit)[2]) #precip Coefficent
        R.Tcoef <- as.numeric(coef(fit)[3]) #Temperatur coefficent
        # test to see if temperature should be left in tegression
        t.sig <- (coef(summary(fit))["ACT.N","Pr(>|t|)"] < 0.05)

      }

      R.BH <- NaN
      R.SC <- NaN
      R.SCT <- NaN
      #Return the results when simple = FALSE
      return(structure(list(lm.sum = c(R.slpe, R.Tcoef, R.intr, R.pval, R.Rval, R.BH, R.SC, R.SCT), temp.sig = t.sig)))
    }}

  while (TRUE) {
    # While loop to handles what happens if temp is not a significant variable
      #On the second loop temperature will be excluded and the analysis run again
    #Get the dimensions of the data for indexing
    len <- dim(ACP.table)[2]
    if (is.null(ACT.table)){
      lines <- dim(ACP.table)[1]
    }else{
      #Add the climate table
      lines <- dim(ACP.table)[1]*dim(ACT.table)[1] #Should allow for variable size tables
    }
    #empy matrix to store lm results
    if (is.null(ACT.table)){ # no Temperature data
      # Build an empy array for the numbers
      m<- matrix(nrow=(lines), ncol=2)
      rownames(m)<- rownames(ACP.table)
      colnames(m)<- c("R^2.Value", "slope")
    }else{
      # Build an empy array for the numbers
      m<- matrix(nrow=(lines), ncol=1)
      rn.names <- NULL
      for (rnm in rownames(ACP.table)){
        rnms <- cbind(rnm, rownames(ACT.table))
        rmx <- paste(rnms[,1] , rnms[,2], sep = ":")
        rn.names <- c(rn.names, rmx)}
      rownames(m)<- rn.names
      colnames(m)<- c("R^2.Value")
    }
    # For loops to call the linreg function and fit a LM to every combination of rainfall and vegetation
    for (n in 1:dim(ACP.table)[1]){
      # browser()
      if (is.null(ACT.table)){
        # Stack the results in the empyt matryx m
        m[n, ] <- linreg(CTSR.VI, ACP.table[n, ])
      }else{
        for (nx in 1:dim(ACT.table)[1]){
          #Get the row name of the correct line
          rn.loop <- paste(rownames(ACP.table)[n], rownames(ACT.table)[nx], sep = ":")
          # Stack the results in the empyt matryx m
          m[rn.loop, ] <- linreg(CTSR.VI, ACP.table[n, ], ACT.table[nx,])
        }}
    }
    # figure out if i need to be checking positive slopes
    if (allow.negative) {
      #To avoild memory duplication
      mx <- m
    }else{ # all the things that enter this else should be single variate regressions
      mx <- matrix(m[m[, "slope"] > 0,], ncol=2)
      colnames(mx)<- c("R^2.Value", "slope")
      rownames(mx)<- rownames(m[m[, "slope"] > 0,])
    }

    #Perform all the regressions
    if (dim(mx)[1] <= 2||allow.negative){ # No positve slpoes or allow negative = TRUE
      # Get the max line and values
      max.line <- which.max(m[, "R^2.Value"])
      #calulate a new linear regression and pull out the summay params
      fulname <- rownames(m)[max.line]
      if (is.null(ACT.table)){ # no temperature data
        sum.res <- linreg(CTSR.VI, ACP.table[fulname, ], simple = FALSE)
        suma <- sum.res$lm.sum
        # Setup for below the else
        precip.nm <- fulname
        CTSR.ATM <- NULL
        t.osp <- NaN
        t.acp <- NaN
      }else{ # include Temperature data
        # namesplit the precip part
        part.nm <- strsplit(fulname, "\\:")[[1]]
        precip.nm <- part.nm[1]
        # browser()
        sum.res <- linreg(CTSR.VI, ACP.table[precip.nm, ], ACT.table[part.nm[2],], simple = FALSE)
        suma <- sum.res$lm.sum
        #pull out tose parms that are only temp related
        CTSR.ATM <- ts(ACT.table[part.nm[2], ], start=c(yst, mst), frequency = 12)
        Tnmsplit <- strsplit(part.nm[2], "\\-")[[1]]
        t.osp <- as.numeric(Tnmsplit[1])
        t.acp <- as.numeric(Tnmsplit[2])

      }
      if (!sum.res$temp.sig) {
        ACT.table=NULL
        allow.negative=allowneg.retest
        # browser()
      }else{
        CTSR.ARF <- ts(ACP.table[precip.nm, ], start=c(yst, mst), frequency = 12)
        # Get the values to return
        nmsplit <- strsplit(precip.nm, "\\-")[[1]]
        osp <- as.numeric(nmsplit[1])
        acp <- as.numeric(nmsplit[2])
        browser()
        return(structure(list(summary=suma, CTSR.precip = CTSR.ARF, CTSR.tmp=CTSR.ATM,
                              CTSR.osp = osp, CTSR.acp = acp, CTSR.tosp = t.osp, CTSR.tacp = t.acp)))
      }

    }else{
      #Filter the negative results out
      rfx <- ACP.table[m[, "slope"] > 0,]
      # Get the max line and values
      max.line <- which.max(mx[, "R^2.Value"])
      fulname <- rownames(mx)[max.line]

      #get all the details needed from the lm()
      sum.res <- linreg(CTSR.VI, ACP.table[fulname, ], simple = FALSE)
      suma <- sum.res$lm.sum

      # get the values to return
      CTSR.ARF <- ts(rfx[max.line, ], start=c(yst, mst), frequency = 12)
      namestr <- rownames(rfx)[max.line]
      nmsplit <- strsplit(namestr, "\\-")[[1]]
      osp <- as.numeric(nmsplit[1])
      acp <- as.numeric(nmsplit[2])
      return(structure(list(summary=suma, CTSR.precip = CTSR.ARF, CTSR.tmp=NULL,
                            CTSR.osp = osp, CTSR.acp = acp, CTSR.tosp = NaN, CTSR.tacp = NaN)))
    }
  }
  browser()
}
