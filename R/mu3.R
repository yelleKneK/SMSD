#' mu3
#' 
#' The function computes the unbiased estimate of the third central moment
#'
#' The function \code{mu3} computes the unbiased estimate of the third central moment
#' of a given vector. The input must numeric (as determined by
#' \code{is.numeric}) and of length at least 3
#' (excluding \code{NA}'s. If \code{na.rm} is \code{TRUE}, the \code{NA}
#' entries in the given vector will be ignored.
#'
#' @param x numeric vector of at least length 3
#' @param na.rm logical value indicating whether NA values should be stripped
#' before the computation proceeds.
#'
#' @return The unbiased sample central moment of \code{X}
#'
#' @author Ken Kelley \email{kkelley@nd.edu},
#' Francis Bilson Darku \email{fbilsond@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@utdallas.edu}
#'
#' @examples
#' x <- rnorm(10)
#' mu3(x)
#' y <- c(1, 3, 1, 5, 6, 3, NA, 3 ,4, NA, 2)
#' mu3(y, na.rm = TRUE)
#'
#' @export mu3
#'
mu3 <- function(x, na.rm = FALSE) {
  if (!is.numeric(x)){
    stop("The argument 'x' needs to be numeric")
  }
  if (!is.logical(na.rm)){
    stop("The argument 'na.rm' needs to be logical")
  }
  if (is.data.frame(x)){
    x <- as.vector(x)
  }
  if (na.rm){
    x <- x[!is.na(x)]
  }
  if (length(x) < 3) {
    stop("Not enough numeric values; needs at least 3")
  }

  n <- length(x)
  w <- sum((x - mean(x))^3)
  mu3n <- n * w/((n - 1) * (n - 2))
  return(mu3n)
}
