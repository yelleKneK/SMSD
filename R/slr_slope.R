#' slr_slope
#'
#' This function calcualtes the slope between two vectors.
#'
#' @param x The first vector
#' @param y The second vector
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}.
#'
#'
#' @return The calculated slope.
#'
#' @examples
#' x <- c(1:4)
#' y <- c(2:5)
#' slr_slope(x, y)
#'
#' @export slr_slope
#'
slr_slope <- function(x, y, na.rm=TRUE){
  if(na.rm){
    x <-  x[!is.na(x)]
    y <-  y[!is.na(y)]
  }
  return(stats::cov(x, y)/stats::var(x))
}
