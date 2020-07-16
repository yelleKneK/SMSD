#' Gini Mean Difference
#'
#' This function calculates the Gini Mean Difference value. The input is
#' a vector of values for which the Gini Mean Difference is calculated.
#'
#' @param x The vector of values for which the Gini Mean Difference will
#' be calculated. This vector should be a numeric vector.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}.
#' @return The calculated Gini Mean Difference.
#'
#' @examples
#' # Household money from income brackets (2010) from Table C of
#' # https://en.wikipedia.org/wiki/Gini_coefficient. 
#' x <- c(13.7, 12.0, 10.9, 13.9, 17.7, 11.4, 12.1, 4.5, 3.9)
#' gini_mean_diff(x)
#'
#' @export gini_mean_diff
#' 
gini_mean_diff <- function(x, na.rm=TRUE)
{
  if(na.rm){
    x <- x[!is.na(x)]
  }
  y <- sort(x)
  n <- length(x)
  i <- seq(1:n)
  s11 <- 2*sum(as.numeric(i*y)) - (n + 1)*sum(y)
  Gn <- s11/choose(n, 2)
  return(Gn)
}
