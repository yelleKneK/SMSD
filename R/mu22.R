#' mu22
#'
#' This an internal function called by other functions.
#'
#' @param x The first vector
#' @param y The second vector
#' @param na.rm This parameter controls whether NA values are removed from
#' the vector. Default is \code{FALSE}.
#'
#' @return The calcualted value for use in other functions.
#'
#' @author Ken Kelley \email{kkelley@nd.edu},
#' Francis Bilson Darku \email{fbilsond@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@iimv.ac.in}
#'
#' @examples
#' x <- c(1:4)
#' y <- c(2:5)
#' mu22(x, y)
#'
#' @export mu22
#'
mu22 <- function(x,y, na.rm = FALSE){
  if(na.rm){
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
  }
  dat <- cbind(x,y)
  n <- length(x)

  k20 <- stats::var(dat[,1])
  k02 <- stats::var(dat[,2])
  k11 <- 1/(n-1) * (s(dat,1,1)-1/n*s(dat,1,0)*s(dat,0,1))
  k22 <- n/((n-1)*(n-2)*(n-3)) * ((n+1)*s(dat,2,2) -
                                    2*(n+1)/n*s(dat,2,1)*s(dat,0,1)
                                  - 2*(n+1)/n*s(dat,1,2)*s(dat,1,0)
                                  - (n-1)/n*s(dat,2,0)*s(dat,0,2)
                                  - 2*(n-1)/n*s(dat,1,1)^2 +
                                    8/n*s(dat,1,1)*s(dat,0,1)*s(dat,1,0)
                                  + 2/n*s(dat,0,2)*s(dat,1,0)^2 +
                                    2/n*s(dat,2,0)*s(dat,0,1)^2
                                  - 6/n^2*s(dat,1,0)^2*s(dat,0,1)^2)

  return(k22 + k20*k02 + 2*k11^2)
}
