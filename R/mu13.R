#' mu13
#'
#' This an internal function called by other functions.
#'
#' @param x The first vector
#' @param y The second vector
#' @param na.rm This parameter controls whether NA values are removed from
#' the vector. Default is \code{FALSE}.
#'
#' @return The calculated value for use in other functions.
#'
#'
#' @examples
#' x <- c(1:4)
#' y <- c(2:5)
#' mu13(x,y)
#'
#' @export mu13
#'
mu13 <- function(x,y, na.rm = FALSE) {
  if(na.rm){
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
  }
  dat <- cbind(x,y)
  n <- length(x)

  k02 <- stats::var(y)
  k11 <- stats::cov(x, y)
  k13 <- n/((n-1)*(n-2)*(n-3)) * (
    (n+1)*s(dat, 1,3) - (n+1)/n*s(dat, 0,3)*s(dat, 1,0) -
      3*(n-1)/n*s(dat, 1,1)*s(dat, 0,2) - 3*(n+1)/n*s(dat,1,2)*s(dat,0,1) +
      6/n*s(dat, 1,1)*s(dat, 0,1)^2 + 6/n*s(dat, 0,2)*s(dat, 0,1)*s(dat, 1,0)
    - 6/n^2*s(dat,1,0)*s(dat,0,1)^3
  )

  return(k13 + 3*k02*k11)
}
