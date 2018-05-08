#' @title Antecedental Rainfall (and temperature) Accumulation calculator for the VI Complete Time Series
#'
#' @description
#' Takes the Complete Time Series Vegetation Index and a table of every possible accumulation period and offset
#' period for precipitation and temperature(optional).  A OLS is calculated \code{\link{lm}} for every combination
#' of VI ~ rainfall (and temperature if that is included).
#' if only the VPR is being calculated, this Function preferences those results where slope>0
#' (increase in rainfall causes an increase in vegetation),returning the rainfall accumulation
#' that has the highest R-squared and a positive slope. If no combination produces a
#' positive slope then the one with the highest R-squared is returned.
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

ACP.calculator <- function(CTSR.VI, ACP.table, ACT.table = NULL, allow.negative = FALSE, allowneg.retest = FALSE) {
  # ==============================================================================================
  # ========== Sanity check the input data ==========
  if (class(CTSR.VI) != "ts") {
    stop("CTSR.VI Not a time series object")}
  if (length(CTSR.VI) != dim(ACP.table)[2]) {
    stop("ACP.table size does not match CTSR.VI")}
  if ((!is.null(ACT.table)) && (dim(ACT.table)[2] != dim(ACP.table)[2])) {
    # Check the temp table (ACT.table)
    stop("ACT.table size does not match CTSR.VI")
    }

  # +++++ Get the start year and start month of the data +++++
  yst <- start(CTSR.VI)[1]
  mst <-  start(CTSR.VI)[2]

  # ==============================================================================================
  # ========== Linear regression function ==========
  # Define a Function to perform linear regression and get out values
  # DESCRIPTION:
  #   This script uses OLS many time. This function allows the usee of
  #   apply to spped up the process

  linreg <- function(CTSR.VI, ACP.N, ACT.N = NULL, simple = TRUE) {
    #   +++++ Simple tells function which values to return, if Simple = TRUE,
    #   only the slope and the R^2 values are returned, if FALSE, all the details

    # +++++ Sanity check the data +++++
    if (sd(ACP.table[n, ]) == 0) {
      #if a combination of acp and osp leads to SD=0 rainfall, this will catch it
      # All values are bad
      if (is.null(ACT.N)) {
        return(c(0, -1))}else{return(0)
        }
    }else{
      # +++++ Check if temperature needs to be considered +++++
      if (is.null(ACT.N)) {# No temperature data
        # perform the regression between VI and precipitation
        fit <- lm(CTSR.VI ~ ACP.N)
        R.Rval <- summary(fit)$r.square
        R.slpe <- as.numeric(coef(fit)[2])

        if (simple) {
          #To speed up looping over all the pixels (simple = TRUE)
          return(c(R.Rval, R.slpe))
          }
        #  +++++ Pull of the rest of the key values +++++
        R.pval <- glance(fit)$p.value
        R.intr <- as.numeric(coef(fit)[1])
        R.Tcoef <- NaN
        # pass a true value for the tempurature sig test
        t.sig <- TRUE

      }else{# temperature data
        # +++++ do the multivariate regression with precip and temperature +++++
        fit <- lm(CTSR.VI ~ ACP.N + ACT.N)
        R.Rval <- summary(fit)$r.square
        R.slpe <- as.numeric(coef(fit)[2])

        if (simple) {
          #To speed up looping over all the pixels (simple = TRUE)
          return(c(R.Rval, R.slpe))
        }
        #  +++++ Pull of the rest of the key values +++++
        R.pval <- glance(fit)$p.value
        R.intr <- as.numeric(coef(fit)[1])
        R.slpe <- as.numeric(coef(fit)[2]) #precip Coefficent
        R.Tcoef <- as.numeric(coef(fit)[3]) #Temperatur coefficent
        # test to see if temperature should be left in regression
        t.sig <- (coef(summary(fit))["ACT.N","Pr(>|t|)"] < 0.05)

      }

      R.BH <- NaN
      R.SC <- NaN
      R.SCT <- NaN
      # +++++ Return the results when simple = FALSE +++++
      return(structure(list(lm.sum = c(R.slpe, R.Tcoef, R.intr, R.pval, R.Rval, R.BH, R.SC, R.SCT), temp.sig = t.sig)))
    }}

  # ==============================================================================================
  # ========== Calculate the opimial climate accumulation  ==========

  while (TRUE) {
    # While loop to handles what happens if temp is not a significant variable
      #On the second loop temperature will be excluded and the analysis run again
    
    # ========== get the size of the tables to be created ==========
    len <- dim(ACP.table)[2]
    if (is.null(ACT.table)) {
      lines <- dim(ACP.table)[1]
    }else{
      #Add the climate table
      lines <- dim(ACP.table)[1]*dim(ACT.table)[1] #Should allow for variable size tables
    }
    # ========== Create a matrix to hold the regression results for all acp/osp ==========
    m <- matrix(nrow = (lines), ncol = 2)
    colnames(m) <- c("R^2.Value", "slope")

    # +++++ Set the row names (vary for temperature and no temp) +++++
    if (is.null(ACT.table)) { # no Temperature data
      rownames(m) <- rownames(ACP.table)
    }else{
      # clear the existing row names
      rn.names <- NULL
      # create new row names that combine precip and temp
      for (rnm in rownames(ACP.table)) {
        rnms <- cbind(rnm, rownames(ACT.table))
        rmx <- paste(rnms[,1] , rnms[,2], sep = ":")
        rn.names <- c(rn.names, rmx)}

      rownames(m) <- rn.names
    }

    # ========== Loop over the acculuation tables ==========
    # For loops to call the linreg function and fit a LM to every combination of rainfall and vegetation
    for (n in 1:dim(ACP.table)[1]) {
      # No temperature, loop over precip
      if (is.null(ACT.table)) {
        # Stack the results in the empyt matryx m
        m[n, ] <- linreg(CTSR.VI, ACP.table[n, ])
      }else{
        # Nested for lopp to allow every combination of temp and precipitation
        for (nx in 1:dim(ACT.table)[1]) {
          #Get the row name of the correct line
          rn.loop <- paste(rownames(ACP.table)[n], rownames(ACT.table)[nx], sep = ":")
          # Stack the results in the empyt matryx m
          m[rn.loop, ] <- linreg(CTSR.VI, ACP.table[n, ], ACT.table[nx,])
        }}
    }

    # ===== Exclude negative precip if allow.negative=False =====
    # figure out if i need to be checking positive slopes
    if (allow.negative) {
      #To avoild memory duplication
      mx <- m
    }else{
      mx <- matrix(m[m[, "slope"] > 0,], ncol = 2)
      colnames(mx) <- c("R^2.Value", "slope")
      rownames(mx) <- rownames(m[m[, "slope"] > 0,])
    }

    # ========== Find the Max R^2 =========
    # +++++ Check and see if negative should be considered +++++
    if (dim(mx)[1] <= 2 || allow.negative) {
      # ++++++++ No positve slpoes or allow negative = TRUE ++++++
      # ++++++ Includes negative precip relationships +++++
      # Get the max line and values
      max.line <- which.max(m[, "R^2.Value"])
      #calulate a new linear regression and pull out the summay params
      fulname <- rownames(m)[max.line]
      if (is.null(ACT.table)) { # no temperature data
        # Pull all infomation from the linreg function
        sum.res <- linreg(CTSR.VI, ACP.table[fulname, ], simple = FALSE)
        suma <- sum.res$lm.sum
        # Setup for below the else
        precip.nm <- fulname
        CTSR.ATM <- NULL
        t.osp <- NaN
        t.acp <- NaN
      } else {# include Temperature data
        # namesplit the precip part
        part.nm <- strsplit(fulname, "\\:")[[1]]
        precip.nm <- part.nm[1]
        # get the regression results
        sum.res <- linreg(CTSR.VI, ACP.table[precip.nm, ], ACT.table[part.nm[2],], simple = FALSE)
        suma <- sum.res$lm.sum
        # pull out tose parms that are only temp related
        CTSR.ATM <- ts(ACT.table[part.nm[2], ], start = c(yst, mst), frequency = 12)
        Tnmsplit <- strsplit(part.nm[2], "\\-")[[1]]
        t.osp <- as.numeric(Tnmsplit[1])
        t.acp <- as.numeric(Tnmsplit[2])

      }
      # ===== If temperature is provided but not significant, setup to do the WHile loop again =====
      if (!sum.res$temp.sig) {
        # +++++ goes round again +++++
        ACT.table = NULL
        allow.negative = allowneg.retest
      }else{
        # +++++ temperature is significant +++++
        CTSR.ARF <- ts(ACP.table[precip.nm, ], start = c(yst, mst), frequency = 12)
        # Get the values to return
        nmsplit <- strsplit(precip.nm, "\\-")[[1]]
        osp <- as.numeric(nmsplit[1])
        acp <- as.numeric(nmsplit[2])
        return(structure(list(
          summary = suma, CTSR.precip = CTSR.ARF, CTSR.tmp = CTSR.ATM,CTSR.osp = osp,
          CTSR.acp = acp, CTSR.tosp = t.osp, CTSR.tacp = t.acp))
          )
      }

    }else{
      # +++++ excludes negative precipitation relationships +++++
      # Get the max line and values
      max.line <- which.max(mx[, "R^2.Value"])
      fulname <- rownames(mx)[max.line]

      if (is.null(ACT.table)) {
        # Filter the negative results out
        rfx <- ACP.table[m[, "slope"] > 0,]
        #get all the details needed from the lm()
        sum.res <- linreg(CTSR.VI, ACP.table[fulname, ], simple = FALSE)
        suma <- sum.res$lm.sum

        # get the values to return
        CTSR.ARF <- ts(rfx[max.line, ], start = c(yst, mst), frequency = 12)
        namestr <- rownames(rfx)[max.line]
        nmsplit <- strsplit(namestr, "\\-")[[1]]
        osp <- as.numeric(nmsplit[1])
        acp <- as.numeric(nmsplit[2])

        precip.nm <- fulname
        CTSR.ATM <- NULL
        t.osp <- NaN
        t.acp <- NaN
      } else {
        # +++++ Temperature data included +++++
        # namesplit the precip part
        part.nm <- strsplit(fulname, "\\:")[[1]]
        precip.nm <- part.nm[1]
        sum.res <- linreg(CTSR.VI, ACP.table[precip.nm, ], ACT.table[part.nm[2],], simple = FALSE)
        suma <- sum.res$lm.sum
        # Get the values to return
        nmsplit <- strsplit(precip.nm, "\\-")[[1]]
        osp <- as.numeric(nmsplit[1])
        acp <- as.numeric(nmsplit[2])
        #pull out tose parms that are only temp related
        CTSR.ATM <- ts(ACT.table[part.nm[2], ], start = c(yst, mst), frequency = 12)
        CTSR.ARF <- ts(ACP.table[precip.nm, ], start = c(yst, mst), frequency = 12)
        Tnmsplit <- strsplit(part.nm[2], "\\-")[[1]]
        t.osp <- as.numeric(Tnmsplit[1])
        t.acp <- as.numeric(Tnmsplit[2])
      }

      if (!sum.res$temp.sig) {
        ACT.table = NULL
        allow.negative = allowneg.retest

      }else {
        return(structure(list(
          summary = suma, CTSR.precip = CTSR.ARF, CTSR.tmp = CTSR.ATM,CTSR.osp = osp,
          CTSR.acp = acp, CTSR.tosp = t.osp, CTSR.tacp = t.acp))
        )
      }

    }
  }
  browser()
}
