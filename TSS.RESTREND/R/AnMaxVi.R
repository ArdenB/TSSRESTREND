#' @title Annual max VI Calculator
#'
#' @description
#'        Takes the montly time series of the VI and calculates the growing season max VI.
#'        In series where the peak occurs in November or December,
#'        an interannual growing season is assessed.
#' @author Arden Burrell, arden.burrell@unsw.edu.au
#' @inheritParams TSSRESTREND
#' @return Max(anu.VI)
#'          The annual (Growing season) max VI. See \code{\link{TSSRESTREND}}
#' @return  Max.Month
#'          The month number where max values were observed (1 for January).
#'          if month > 12, the peak was detected in Nov, Dec, Jan.
#'          In this case the peak seasonal value and position is used.
#' @return index(VI.index)
#'          the index of the CTSR.VI ts that the anu.VI values occur at.
#'          See \code{\link{TSSRESTREND}}. Note.R indexs from 1 rather than 0.
#' @export
#' @examples
#' anmax <- AnMaxVI(stdRESTRENDCTSR$cts.NDVI)
#' print(anmax)


AnMaxVI <- function(CTSR.VI) {
  # ==============================================================================================
  # ========== Sanity check the input data ==========
  if (class(CTSR.VI) != "ts") {
    stop("CTSR.VI Not a time series object")
  }
  # ==============================================================================================
  # =========== Organise the data and key variables ==========
  # +++++ Work out the start and end dates +++++
  sty <- start(CTSR.VI)[1]
  stm <- start(CTSR.VI)[2]
  eny <- end(CTSR.VI)[1]
  enm <- end(CTSR.VI)[2]
  # ++++++ Create a blank matrix+++++
  # matrix is sorted by year so its easier to find the annual max values
  m <- matrix(nrow = (eny - sty + 1), ncol = 12)
  rownames(m) <- c(sty:eny)
  colnames(m) <- c(month.abb[1:12])
  index <- 1
  # +++++  loop over each year To set the names +++++
  for (yr in sty:eny) {
    for (mon in 1:12) {
      if (yr == sty & mon < stm) {# Ignore months before the start
        m[toString(yr), month.abb[mon]] = NaN
      }else if (yr == eny & mon > enm) { # Ignore months After the end of the data
        m[toString(yr), month.abb[mon]] = NaN
      }else{
        m[toString(yr), month.abb[mon]] = CTSR.VI[index]
      }
      index <- (index + 1)
    }
  }
  # ==============================================================================================
  # =========== Work out the seasonal max NDVI values ==========

  # =========== Get the annual max NDVI values ==========
  # +++++ build time series of max values, max months and max index +++++
  anmax.ts <- ts(apply(m, 1, max, na.rm = TRUE), start = sty, frequency = 1)
  whmax.ts <- ts(apply(m, 1, which.max), start = sty, frequency = 1)
  index.ts <- ts(c(sty:eny), start = sty, frequency = 1)
  # +++++ add these ts' to a dataframe
  df <- data.frame( Max = anmax.ts, Max.month = whmax.ts, index = index.ts)

  # ========== Test seasonality ============
  # Loop over each year
  for (row in 1:(dim(df)[1] - 1)) {
    # Check and see of the max values are occuring in nov and december
    if (
      (df$Max.month[row] == 11 || df$Max.month[row] == 12) &&
      (df$Max.month[row + 1] == 1 || df$Max.month[row + 1] == 2)
      ) {# Value occurs in the NDJF period
      # Add Jan anf Feb data and re check the max
      test <- c(m[row, 11:12], m[row + 1, 1:2])
      df$Max.month[row] = which.max(test) + 10
      df$Max[row] = max(test, na.rm = TRUE)
      # exclude jan and feb for the subsuquent years calculation
      try({# Should only fail on the last year, this is to catch that
        nxtyr <- m[row + 2, 2:12]
        df$Max.month[row + 1] = which.max(nxtyr) + 2
        df$Max[row + 1] = max(nxtyr)
      }, silent = TRUE)

    }
  }
  df$index <- ((df$index - sty) * 12 + df$Max.month)
  if (stm != 1) {# the first month of data is not january
    df$index <- (df$index - (stm - 1))
    warning(
      "TSS-RESTREND should work with any start month.
      However this has not been tested. Additional
      caution should be excercised when considering
      the results ")
  }
  return(df)
}
