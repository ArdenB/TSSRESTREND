#' @title  Print function for class TSSRESTREND
#'
#' @description print function for objects of class TSSRESTREND
#' @inheritParams plot.TSSRESTREND
#' @export

print.TSSRESTREND <- function(x, ...) {
  print(x$summary) #SUmmary table
  print(x$ols.summary) # all of the regression and stat test summary
}
