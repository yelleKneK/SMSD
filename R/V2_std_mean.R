#' V2_std_mean
#' 
#' This function is used internally to calculate the estimate of the asymptotic variance of the standardized mean.
#'
#' @param x The vector for which to calculate the standardized mean.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{FALSE}.
#'
#' @return The estimated asymptotic variance of the standardized mean.
#'
#' @references
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in parameter estimation for a general class of effect sizes: A sequential approach. \emph{Psychological Methods}, \emph{23}, 226–243.
#' 
#'@author Bhargab Chattopadhyay \email{bhargab@iimv.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#' @examples
#' x <- c(1:4)
#' V2_std_mean(x)
#'
#' @export V2_std_mean
#'
V2_std_mean <- function(x, na.rm=FALSE){
  if(na.rm){
    x <-  x[!is.na( x)]
  }
  xbar <- mean(x)
  sn2 <- stats::var(x)

  vn2 <- 1 - xbar^2/(4*sn2) - xbar*mu3(x)/sn2^2 + xbar^2*mu4(x)/(4*sn2^3)
  return(vn2)
}
