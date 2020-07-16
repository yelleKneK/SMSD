#' V2_gini_index
#'
#' This function is used internally to calculate the Gini index.
#'
#' @param x The vector for which to calculate the Gini index.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{FALSE}.
#'
#' @return The estimated estimated asymptotic variance of the Gini index. 
#'
#' @author Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#' @examples
#' x <- c(1:4)
#' V2_gini_index(x)
#'
#' @export V2_gini_index
#'
V2_gini_index <- function(x, na.rm=FALSE){
  if(na.rm){
    x <-  x[!is.na( x)]
  }
  n <- length(x)
  xbar <- mean(x)
  GMD <- gini_mean_diff(x)
  tn <- 0.5*gini_mean_diff(x^2)
  u <- sapply(1:n, function(i) gini_mean_diff(x[-i]))
  z <- n*GMD - (n - 2)*u
  sw2 <- stats::var(z)

  t1 <- GMD^2*stats::var(x)/(4*xbar^4)
  t2 <- GMD*tn/xbar^3
  t3 <- GMD^2/xbar^2
  t4 <- sw2/(4*xbar^2)
  vn2 <- t1 - t2 + t3 + t4;
  return(vn2)
}
