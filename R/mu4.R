#' mu4
#'
#' A function to compute the unbiased sample central fourth moment.
#'
#' The function \code{mu4} computes the unbiased sample central moment
#' of a given vector. The input must numeric (as determined by
#' \code{is.numeric}) and of length at least 4
#' (excluding \code{NA}'s. If \code{na.rm} is \code{TRUE}, the \code{NA}
#' entries in the given vector will be ignored.
#'
#' @param x A numeric vector of at least length 4.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{FALSE}.
#'
#' @return The calculated moment.
#'
#' @examples
#' x <- c(1:5)
#' mu4(x)
#'
#' @export mu4
#'
mu4 <- function(x, na.rm = FALSE) {
  if (!is.numeric(x)){
    stop("The argument 'x' needs to be numeric")
  }

  if (!is.logical(na.rm)){
    stop("The argument 'na.rm' needs to be logical")
  }

  if (is.data.frame(x)){
    x <- as.vector(x)
  }
  if (na.rm == TRUE){
    x <- x[!is.na(x)]
  }
  if (length(x) < 4){
    stop("Not enough numeric values; needs at least 4")
  }


  n <- length(x)
  w1 <- n^2 * sum((x - mean(x))^4)
  w2 <- (-2*n + 3) * sum(x^4)
  w3 <- (8*n - 12) * mean(x) * sum(x^3)
  w4 <- (-6 + 9/n) * (sum(x^2))^2

  mu4n <- (w1 + w2 + w3 + w4)/((n - 1) * (n - 2) * (n - 3))
  return(mu4n)
}
