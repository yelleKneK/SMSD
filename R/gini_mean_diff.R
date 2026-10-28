#' Gini's Mean Difference
#'
#' Estimates Gini's mean difference, the average absolute difference between
#' two randomly chosen members of a population. It is a measure of spread on
#' the scale of the data (like the standard deviation, but based on absolute
#' rather than squared differences) and is the numerator of the Gini index.
#'
#' The population quantity is \eqn{\Delta = E|X_1 - X_2|} for two independent
#' draws. The estimate returned is the U-statistic
#' \deqn{\hat{\Delta}_n = \binom{n}{2}^{-1} \sum_{i < j} |x_i - x_j|,}
#' the average of \eqn{|x_i - x_j|} over all pairs of distinct observations,
#' which is unbiased for \eqn{\Delta}. It does not change when a constant is
#' added to every value. Kelley, Darku, and Chattopadhyay (2018, Example 7)
#' name Gini's mean difference among the parameters covered by their general
#' sequential procedure; Chattopadhyay and De (2016) use it to estimate the
#' Gini index, \eqn{\hat{\Delta}_n / (2 \bar{X}_n)} (see
#' \code{\link{gini_index}}).
#'
#' @param x A numeric vector with at least two observations. Any real values
#'   are allowed.
#' @param na.rm Should missing values be removed before the calculation?
#'   Default \code{TRUE}. With \code{na.rm = FALSE} a missing value is an
#'   error.
#'
#' @return A single nonnegative number, the estimated Gini mean difference,
#'   in the units of \code{x}.
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
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#' parameter estimation for a general class of effect sizes: A sequential
#' approach. \emph{Psychological Methods}, \emph{23}(2), 226--243.
#' \doi{10.1037/met0000127}
#'
#' @examples
#' # Annual salaries (in thousands of dollars) of ten employees of a firm:
#' # two employees chosen at random differ by this much on average.
#' salary <- c(28, 31, 35, 38, 42, 47, 55, 64, 90, 170)
#' gini_mean_diff(salary)
#'
#' # The Gini index is the mean difference over twice the mean.
#' gini_mean_diff(salary) / (2 * mean(salary))
#' gini_index(salary)
#'
#' @seealso \code{\link{gini_index}}, \code{\link{seq_aipe_gini}}
#' @export
gini_mean_diff <- function(x, na.rm = TRUE) {
  if (missing(x)) {
    .smsd_stop("'x' must be supplied: a numeric vector.")
  }
  na.rm <- .check_flag(na.rm, "na.rm")
  x <- .check_data(x, "x", min_n = 2L, na.rm = na.rm, advice = NULL)
  y <- sort(x)
  n <- as.numeric(length(y))
  # The sum is formed on data divided by a power of two of their magnitude
  # and returned to the units of the data, so that it cannot overflow; see
  # .aipe_scale().
  scale <- .aipe_scale(y)
  y <- y / scale
  rank_weight <- 2 * seq_along(y) - n - 1
  # Sum over pairs of |x_i - x_j| written with the order statistics; the
  # weights sum to zero, so the values are centered first (the statistic is
  # location invariant and the centered sum does not cancel).
  scale * (2 * sum(rank_weight * .aipe_center(y)) / (n * (n - 1)))
}
