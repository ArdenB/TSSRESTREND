#' @title Chow test on detected breakpoints
#'
#' @description
#' Takes the breakpoints detected by \code{\link{VPR.BFAST}}, finds the most significant one then tests it in both
#' the residuals and the VPR.
#'
#' @author Arden Burrell, arden.burrell@unsw.edu.au
#'
#' @importFrom strucchange sctest
#' @importFrom broom glance
#'
#' @inheritParams TSSRESTREND
#'
#' @param breakpoints
#' vector containing the breakpoints detected by \code{\link{VPR.BFAST}} (bkps)
#'
#' @return \bold{n.Method}
#'          The method that the ts should be tested with.  TSSRESTREND internal communication.
#' @return \bold{bp.summary}
#'          Summary of the most signifcant breakpoint in the residuals and VPR.
#'          see \code{\link[strucchange]{sctest}}
#' @return \bold{allbp.index}
#'          the Annual index of every breakpoint. Used by \code{\link{plot.TSSRESTREND}}
#' @return bpRESID.chow
#'          Chow test in the VPR residuals. See \code{\link[strucchange]{sctest}}
#' @return bpVPR.chow
#'          Chow test in the VPR. See \code{\link[strucchange]{sctest}}
#' @export
#'
#' @examples
#' #Test the complete time series for breakpoints
#' VPRBFdem <- VPR.BFAST(segVPRCTSR$cts.NDVI, segVPRCTSR$cts.precip)
#' bp <- as.numeric(VPRBFdem$bkps)
#' #test the significance of the breakpoints
#' reschow <- CHOW(segVPR$max.NDVI, segVPR$acum.RF, segVPR$index, bp)
#' print(reschow)
#'

CHOW <- function(anu.VI, acu.RF, VI.index, breakpoints, acu.TM = NULL, sig = 0.05) {
  # ==============================================================================================
  # ========== Sanity check the input data ==========
  if (class(anu.VI) != "ts")
    stop("anu.VI Not a time series object")
  if (class(acu.RF) != "ts")
    stop("acu.RF Not a time series object")
  ti <- time(anu.VI)
  f <- frequency(anu.VI)
  #check the two ts object cover the same time period
  ti2 <- time(acu.RF)
  f2 <- frequency(acu.RF)
  if (!identical(ti, ti2))
    stop("ts object do not have the same time")
  if (!identical(f, f2))
    stop("ts object do not have the same frequency")
  # need to test the breakpoints to make sure they are numeric
  if (class(breakpoints) != "numeric")
    stop("Breakpoints are not class numeric")

  # ==============================================================================================
  # ========== Look at the breakpoints and set up variables ==========

  # +++++ count of all the breakpoints +++++
  len <- length(breakpoints)

  # convert the breakpoints into year posistion
  # the breakpoint loc will be the last anaual max before the breakpoint

  #empty variables to add to a datafram
  empty.1 <- NaN
  empty.2 <- NaN
  empty.3 <- NaN
  # Setup empty data frames to stor infomation
  ind.df <- data.frame(
    abs.index = breakpoints, yr.index = empty.1, reg.sig = empty.2, VPR.bpsig = empty.3
    )
  bp.ind <- data.frame(abs.index = breakpoints, yr.index = NaN)
  # ===== Get the year indexs of the breakpoints =====
  for (bp in 1:length(breakpoints)) {
    bpv = ind.df$abs.index[bp]
    for (n in 1:length(VI.index)) {
      if (bpv >= VI.index[n] & bpv <= VI.index[n + 1]) {
        ind.df$yr.index[bp] = n
        bp.ind$yr.index[bp] = n
      }
    }
  }
  # ===== create the lm for the VPR and test is significance, =====
  #   split VPR sig from VPR insignificant
  if (is.null(acu.TM)) {# no temp data
    VPR.fit <- lm(anu.VI ~ acu.RF)
  }else{# temp data
    VPR.fit <- lm(anu.VI ~ acu.RF+acu.TM)
  }
  if (glance(VPR.fit)$p.value > sig) {
    # VPR significance below critical threshold
    ind <- acu.RF #independent variable
    dep <- anu.VI #dependent variable
    Method = "seg.VPR"
    if (!is.null(acu.TM)) {# temperature data
      ind2 <- acu.TM
    }else {# No temperature data
        ind2 = NULL
    }
  } else {
    # VPR is significant
    ind <- ti
    dep <- VPR.fit$residuals
    Method = "seg.RESTREND"
    ind2 = NULL
  }
  # ========== Iterate over each of the breakpoints ==========
  # Description:
  #   Will loop over each breakpoint, removing the least significant each time
  #   until there is only one. Then testing that in the residuals and the VPR

  while (TRUE) {
    # ===== loop over the breakpoints. The first bp is 1,etc =====
    for (bp.num in 1:nrow(ind.df)) {
      bp = ind.df$yr.index[bp.num]
      # +++++ set the start and ends indexes  +++++
      # Start
      if (identical(ind.df$yr.index[bp.num - 1], numeric(0))) {
        bp.start = 1
      }else{
        bp.start = ind.df$yr.index[bp.num - 1]
      }
      # End
      if (is.na(ind.df$yr.index[bp.num + 1])) {
        bp.end = length(dep)
      }else{
        bp.end = ind.df$yr.index[bp.num + 1]}
      # Locate the tested breakpoint
      bkp = bp - (bp.start - 1)

      #  =========== Perform the CHOW test on the Residuals  ==========
      if (is.null(ind2)) { # no temp
        chow <- sctest(dep[bp.start:bp.end] ~ ind[bp.start:bp.end], type = "Chow", point = bkp)
        }else{
          # This is here to catch VPR fails that have a breakpoint
          #  Currently this shouldn't be used because when the VCR fails it reverts to a VPR
          #  This may need to be changes in the future
          browser("This may need further development")
          chow <- sctest(dep[bp.start:bp.end] ~ ind[bp.start:bp.end]+ind2[bp.start:bp.end], type = "Chow", point = bkp)
        }
      # Pull out the Chow results
      ind.df$reg.sig[bp.num] = chow$p.value

    }
    # ========== Test pixels with more than one breakpoints ==========
    if (nrow(ind.df) > 1) {
      # ===== More than one BP =====
      # delete breakpoint with the largest p values (lowest significance)
      ind.df <- ind.df[!(1:nrow(ind.df) %in% (which.max(ind.df$reg.sig))),]
    }else if (nrow(ind.df) == 1) {
      # ===== last breakpoint standing =====
      # Check if temperature needs to be testeds
      if (is.null(acu.TM)) {
        VPR.chow <- sctest(anu.VI ~ acu.RF, type = "Chow", point = ind.df$yr.index[1])
      }else{
        browser()
        VPR.chow <- sctest(anu.VI ~ acu.RF + acu.TM, type = "Chow", point = bkp)
      }
      # Pull out the p value
      ind.df$VPR.bpsig[1] = VPR.chow$p.value

      # ========== Set the test method and return chow test results ==========

      if (Method == "seg.VPR" & ind.df$reg.sig[1] > sig) { # cant chow non-sig residulas (bpRESID.chow = FALSE)
        ind.df$reg.sig[1] = NaN
        return(structure(list(
          # +++++ Passing to the RESTREND function to catch broken VPR +++++
          # VPR failed and significance of breakpoint failed
          n.Method = "RESTREND", bp.summary = ind.df, allbp.index = bp.ind,
          bpRESID.chow = FALSE, bpVPR.chow = VPR.chow), class = "CHOW.Object")
          )
      } else if (Method == "seg.VPR" & ind.df$reg.sig[1] <= sig) {
        # VPR failed and breakpoint is significant
        ind.df$reg.sig[1] = NaN
        return(structure(list(
          n.Method = "seg.VPR", bp.summary = ind.df, allbp.index = bp.ind,
          bpRESID.chow = FALSE, bpVPR.chow = VPR.chow), class = "CHOW.Object")
          )
      } else if (Method == "seg.RESTREND" & ind.df$reg.sig[1] > sig) {
        # Breakpoint not significant in the residuals
        return(structure(list(
          n.Method = "RESTREND", bp.summary = ind.df, allbp.index = bp.ind,
          bpRESID.chow = chow, bpVPR.chow = VPR.chow), class = "CHOW.Object")
          )
      } else if (Method == "seg.RESTREND" & ind.df$reg.sig[1] <= sig) {
        # Breakpoint is significant in the residuals
        if (VPR.chow$p.value > sig) {
          # Breakpoint is significant in the residuals not in the VPR/VCR
          return(structure(list(
            n.Method = "seg.RESTREND", bp.summary = ind.df, allbp.index = bp.ind,
            bpRESID.chow = chow, bpVPR.chow = VPR.chow), class = "CHOW.Object")
            )
        }else if (VPR.chow$p.value <= sig) {
          # Breakpoint is significant in the VPR
          return(structure(list(
            n.Method = "seg.VPR", bp.summary = ind.df, allbp.index = bp.ind,
            bpRESID.chow = chow, bpVPR.chow = VPR.chow), class = "CHOW.Object"))
        }
      }else{
        warning("Error in Method and df, exit point 1")
        return(FALSE)
      }
    }else{
      warning("ind.df shape is wrong, exited to avoid infinite loop, exit point 2")
      return(FALSE)
    }
  }
}
