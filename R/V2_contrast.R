#' V2_contrast
#'
#' This is an internal function called to calculate the estimate of the asymptotic variance of a contrast.
#'
#' @param data The data set for which to calculate the contrast.
#' @param coef The set of coefficients.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{FALSE}.
#'
#' @return The coefficients multiplied by the variance of each column
#'
#' @import matrixStats
#'
#' @export V2_contrast
#'
#' @author Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#' @examples
#' dat_1 <- matrix(1:10, ncol=2)
#' V2_contrast(dat_1, coef=1)
#'
V2_contrast <- function(data, coef, na.rm=FALSE)
{
  S2 <- matrixStats::colVars(data, na.rm = na.rm)
  vn2 <- sum(coef^2*S2)
}
