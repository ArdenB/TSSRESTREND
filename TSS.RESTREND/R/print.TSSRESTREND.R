#' @title  Print function for class TSSRESTREND class
#'
#' @description print function for class TSSRESTREND
#' @inheritParams plot.TSSRESTREND
#' @export

print.TSSRESTREND <- function(x){
  print(x$summary)
  print(x$ols.summary)
}
