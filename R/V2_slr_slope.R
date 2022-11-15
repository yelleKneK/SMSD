#' V2_slr_slope
#'
#' This function is used internally to calculate the estimate of the asymptotic variance of the slope of two vectors.
#'
#' @param x The independent variable data vector.
#' @param y The dependent variable data vector.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{FALSE}.
#'
#' @return The estimated asymptotic variance of the calculated slope.
#'
#' @author Bhargab Chattopadhyay \email{bhargab@iimv.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#' @examples
#' x <- c(1:4)
#' y <- c(2:5)
#' V2_slr_slope(x, y)
#'
#' @export V2_slr_slope
#'
V2_slr_slope <- function(x, y, na.rm=FALSE){
  if(na.rm){
    x <-  x[!is.na( x)]
    y <-  y[!is.na( y)]
  }
  n <- length(x)
  t1 <- mu22(x,y)/stats::var(x)^2
  t2 <- 2*stats::cov(x,y)*mu31(x,y)/stats::var(x)^3
  t3 <- stats::cov(x,y)^2*mu4(x)/stats::var(x)^4
  vn2 <- max(t1-t2+t3, n^(-3));
  return(vn2)
}
