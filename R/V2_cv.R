#' V2_cv
#'
#' This function is used internally to calculate the coefficient of variation.
#'
#' @param x The vector for which to calculate the coefficient of variation.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{FALSE}.
#'
#' @return The calculated coefficient of variation.
#'
#' @export V2_cv
#'
#' @author Bhargab Chattopadhyay \email{bhargab@iimv.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#' @examples
#' x <- c(1:5)
#' V2_cv(x)
#'
V2_cv <- function(x, na.rm=FALSE){
  if(na.rm){
    x <-  x[!is.na( x)]
  }
  t1 <- mu4(x)/(4*mean(x)^2*stats::var(x))
  t2 <- stats::var(x)/(4*mean(x)^2)
  t3 <- mu3(x)/mean(x)^3
  t4 <- stats::var(x)^2/(mean(x)^4)
  n <- length(x)
  vn2 <- max(t1 - t2 - t3 + t4, n^(-3))
  return(vn2);
}
