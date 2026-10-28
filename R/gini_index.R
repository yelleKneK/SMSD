#' Gini Index
#'
#' Estimates the Gini index of inequality from a sample of nonnegative
#' values such as incomes, wealth, sales, or citation counts. The Gini index
#' is 0 when every unit has the same amount and approaches 1 when a single
#' unit holds everything.
#'
#' The population Gini index is \eqn{G = \Delta / (2\mu)}, where
#' \eqn{\Delta = E|X_1 - X_2|} is Gini's mean difference (see
#' \code{\link{gini_mean_diff}}) and \eqn{\mu} is the mean. It is defined for
#' nonnegative data with a positive mean. The estimate returned is
#' \deqn{G_n = \frac{\hat{\Delta}_n}{2 \bar{X}_n},}
#' the sample Gini mean difference (the average of \eqn{|x_i - x_j|} over all
#' pairs) divided by twice the sample mean, which is the estimator used in
#' the sequential procedure of Chattopadhyay and De (2016) that
#' \code{\link{seq_aipe_gini}} implements. Because the average is over pairs
#' of distinct observations, \eqn{G_n} equals \eqn{n / (n - 1)} times the
#' "plug in" Gini index that some software reports; the two agree as the
#' sample grows. The estimate does not change when every value is multiplied
#' by the same positive constant.
#'
#' @param x A numeric vector of nonnegative values with at least two
#'   observations, not all zero.
#' @param na.rm Should missing values be removed before the calculation?
#'   Default \code{TRUE}. With \code{na.rm = FALSE} a missing value is an
#'   error.
#'
#' @return A single number between 0 and 1, the estimated Gini index.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#' Francis Bilson Darku \email{fbilsond@@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @references
#' Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
#' pre-specified error bound. \emph{Econometrics}, \emph{4}(3), 30.
#' \doi{10.3390/econometrics4030030}
#'
#' @examples
#' # Annual salaries (in thousands of dollars) of ten employees of a firm.
#' salary <- c(28, 31, 35, 38, 42, 47, 55, 64, 90, 170)
#' gini_index(salary)
#'
#' # Equal salaries give 0; the index does not depend on the currency unit.
#' gini_index(rep(50, 10))
#' gini_index(1000 * salary)
#'
#' @seealso \code{\link{seq_aipe_gini}} for the sequential confidence
#'   interval; \code{\link{gini_mean_diff}}.
#' @export
gini_index <- function(x, na.rm = TRUE) {
  if (missing(x)) {
    .smsd_stop("'x' must be supplied: a numeric vector of nonnegative ",
               "values.")
  }
  na.rm <- .check_flag(na.rm, "na.rm")
  x <- .check_data(x, "x", min_n = 2L, na.rm = na.rm, nonnegative = TRUE,
                   advice = NULL)
  y <- sort(x)
  n <- as.numeric(length(y))
  if (y[length(y)] <= 0) {
    .smsd_stop("Every value in 'x' is zero, so the Gini index is ",
               "undefined: it needs a positive mean.")
  }
  # The index does not depend on the scale of the data, which are divided
  # by a power of two of their magnitude so that the total cannot overflow;
  # see .aipe_scale().
  y <- y / .aipe_scale(y)
  rank_weight <- 2 * seq_along(y) - n - 1
  # Sum over pairs of |x_i - x_j| written with the order statistics; the
  # weights sum to zero, so the centered values give the same sum without
  # cancellation.
  index <- sum(rank_weight * .aipe_center(y)) / ((n - 1) * sum(y))
  # The index is at most 1 for nonnegative data, with equality when one unit
  # holds everything. Rounding in the last place can carry the computed
  # value past that (to 1.0000000000000002 for c(0, 0, 0, 0, 0, 0.1)), so it
  # is held to 1; this is a bound, not rounding.
  min(1, index)
}
