#' mu31
#'
#' This an internal function called by other functions.
#'
#' @param x The first vector
#' @param y The second vector
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{FALSE}.
#'
#' @return This function returns the calculated value for use in
#' other functions.
#'
#' @examples
#' x <- c(1:4)
#' y <- c(2:5)
#' mu31(x,y)
#'
#' @export mu31
#'
mu31 <- function(x, y,na.rm=FALSE)
{
  if(na.rm){
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
  }
  return(mu13(y, x))
}
