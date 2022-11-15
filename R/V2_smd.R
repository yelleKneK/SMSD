#' V2_smd
#'
#' This function is used internally to calculate the estimate of the asymptotic variance of the standardized mean difference.
#'
#' @param x The first vector.
#' @param y The second vector.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{FALSE}.
#'
#' @return The estimated asymptotic variance of the standardized mean.
#'
#'@author Bhargab Chattopadhyay \email{bhargab@iimv.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#' @examples
#' x <- c(1:4)
#' y <- c(2:5)
#' V2_smd(x, y)
#'
#' @import stats
#'
#' @export V2_smd
#'
V2_smd <- function(x, y, na.rm=FALSE)
{
  if(na.rm){
    x <-  x[!is.na( x)]
    y <-  y[!is.na( y)]
  }
  Spn <- sqrt(0.5*(stats::var(x) + stats::var(y)))
  
  #t1 is correct here. In the article, the 0.5 (i.e., 1/2) was missing.
  t1 <- 0.5*(mean(x) - mean(y))*(mu3(x) - mu3(y))/Spn^4
  t2 <- ((mean(x) - mean(y))^2)/(4*Spn^6)
  t3 <- ((mu4(x) + mu4(y))/4 - (Spn^4)/2) # correct formula;
  n <- length(x)
  vn2 = max(2 - t1 + t2*t3, n^(-3))
  return(vn2)
}
