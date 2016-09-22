#' @title  Print function for class TSSRESTREND class
#'
#' @description print function
#'
#' @export

print.TSSRESTREND <- function(x , ...){
  print(x$summary)
  print(x$ols.summary)
}
