#' Simple Linear Regression Slope
#'
#' Calculates the ordinary least squares slope of the simple linear
#' regression of \code{y} on \code{x}: the expected change in \code{y} for a
#' one unit increase in \code{x}. This is the point estimate used by the
#' sequential procedure \code{\link{seq_aipe_reg_coef}}.
#'
#' The slope is \eqn{b_1 = s_{xy} / s_x^2}, the sample covariance of
#' \code{x} and \code{y} divided by the sample variance of \code{x}
#' (Kelley, Darku, and Chattopadhyay, 2018, Equation 35). It is the same
#' number as \code{coef(lm(y ~ x))[2]}.
#'
#' @param x A numeric vector, the predictor.
#' @param y A numeric vector, the outcome. It must have the same length as
#'   \code{x}; the two are paired by position.
#' @param na.rm Should a pair be removed when either of its values is
#'   missing? Default \code{TRUE}. With \code{na.rm = FALSE} a missing value
#'   is an error.
#'
#' @return A single number, the estimated slope, in units of \code{y} per
#'   unit of \code{x}. At least two complete pairs are needed, and \code{x}
#'   must vary: a predictor with no variation is an error, because the slope
#'   is then undefined. Values of \code{x} that agree to ten significant
#'   digits count as no variation (such differences are usually rounding
#'   error); if they are real, subtract a constant from \code{x} first. A
#'   slope too large for a double (which takes units of \code{y} and
#'   \code{x} more than about 300 orders of magnitude apart) is an error
#'   that says so.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#' Francis Bilson Darku \email{fbilsond@@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @references
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#' parameter estimation for a general class of effect sizes: A sequential
#' approach. \emph{Psychological Methods}, \emph{23}(2), 226--243.
#' \doi{10.1037/met0000127}
#'
#' @examples
#' # Hours of training and a performance rating for eight employees.
#' hours  <- c(2, 4, 5, 7, 8, 10, 12, 15)
#' rating <- c(48, 55, 53, 61, 66, 64, 73, 80)
#' slr_slope(hours, rating)
#'
#' @seealso \code{\link{seq_aipe_reg_coef}}; in \pkg{DMAR},
#'   \code{\link[DMAR]{ci_reg_coef}} for a fixed sample size confidence
#'   interval.
#' @export
slr_slope <- function(x, y, na.rm = TRUE) {
  na.rm <- .check_flag(na.rm, "na.rm")
  pair <- .aipe_check_pair(x, y, "x", "y", min_n = 2L, na.rm = na.rm,
                           plain_estimator = TRUE)
  if (.is_constant(pair$x)) {
    .smsd_stop("'x' has no variation (its values agree to ten significant ",
               "digits; the first is ", format(pair$x[1L]), "), so the ",
               "slope is undefined. If the variation is real but small ",
               "beside the size of the values, subtract a constant from ",
               "'x' first.")
  }
  slope <- .slr_slope_value(pair$x, pair$y)
  if (!is.finite(slope)) {
    .smsd_stop("The slope is too large to represent. It is in units of 'y' ",
               "per unit of 'x', and the values of 'y' (the largest is ",
               format(max(abs(pair$y))), " in absolute value) are too ",
               "large beside those of 'x' (the largest is ",
               format(max(abs(pair$x))), "). Express 'x' or 'y' in other ",
               "units.")
  }
  slope
}

# The slope of y on x for clean paired data in which x varies: the value
# slr_slope() returns, or Inf or -Inf when the slope is too large for a
# double (slr_slope() stops then; seq_aipe_reg_coef() reports an estimate
# that is not finite).
#
# The covariance and variance are formed on x and y each divided by a power
# of two of its magnitude, so that the sums of squares cannot overflow or
# underflow (see .aipe_scale()), and the slope on that scale, a number of
# ordinary size, is returned to the units of the data. A slope of exactly
# zero is zero in any units; multiplying it by the ratio of the two scales
# would give NaN when that ratio is infinite (x subnormal, say). The ratio
# is applied as a power of two in steps, because the ratio itself can
# overflow or underflow when the slope in the units of the data does not.
#' @noRd
.slr_slope_value <- function(x, y) {
  scale_x <- .aipe_scale(x)
  scale_y <- .aipe_scale(y)
  slope <- stats::cov(x / scale_x, y / scale_y) / stats::var(x / scale_x)
  if (isTRUE(slope == 0)) {
    return(0)
  }
  # Both scales are powers of two, so their base two logarithms are whole
  # numbers (rounded in case log2() is off in the last place).
  e <- round(log2(scale_y)) - round(log2(scale_x))
  step <- 2^(sign(e) * 1000)
  while (abs(e) > 1000) {
    slope <- slope * step
    e <- e - sign(e) * 1000
  }
  slope * 2^e
}
