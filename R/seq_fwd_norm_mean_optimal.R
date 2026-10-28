#' Known Variance Sample Size for a Fixed Width Interval for a Normal Mean
#'
#' Returns the sample size at which a confidence interval for the mean of
#' a normal distribution has half width \code{d} when the population
#' variance is known. This is the benchmark that the purely sequential
#' (\code{\link{seq_fwd_norm_mean}}) and two stage
#' (\code{\link{ts_fwd_norm_mean}}) fixed width procedures aim for without
#' knowing the variance. It is useful for planning with an assumed
#' variance and for judging how much a sequential procedure costs relative
#' to knowing the variance in advance.
#'
#' With known variance \eqn{\sigma^2} the interval
#' \eqn{\bar{X}_n \pm d} covers the mean with probability
#' \eqn{2\Phi(d \sqrt{n} / \sigma) - 1}, which reaches the confidence
#' level \eqn{1 - \alpha} when \eqn{n} is at least
#' \deqn{C = \frac{z_{\alpha/2}^2 \, \sigma^2}{d^2},}
#' where \eqn{z_{\alpha/2}} is the upper \eqn{\alpha/2} point of the
#' standard normal distribution. \eqn{C} is rarely a whole number; the
#' function returns the smallest whole sample size at which the coverage
#' is at least the confidence level, which is \eqn{C} rounded up and never
#' less than 1. (When \eqn{C} is within rounding error of a whole number,
#' the coverage formula decides between that number and the next.) The
#' unrounded value is
#' \code{qnorm((1 - conf_level) / 2, lower.tail = FALSE)^2 * variance / d^2}.
#'
#' Because the result is proportional to the variance, a planning value
#' that is too small gives a sample that is too small by the same factor.
#' When the variance is not known, use one of the procedures that
#' estimate it as the data arrive.
#'
#' @param variance The population variance \eqn{\sigma^2} (the square of
#'   the standard deviation, not the standard deviation itself). A single
#'   positive number, in the squared units of the data.
#' @param d The half width of the confidence interval (the margin of
#'   error), in the units of the data. A single positive number. The full
#'   width is \eqn{2d}.
#' @param conf_level The confidence level, strictly between 0 and 1.
#'   Default 0.95.
#'
#' @return A length one integer of at least 1: the smallest sample size
#'   \eqn{n} for which the known variance interval \eqn{\bar{X}_n \pm d}
#'   has confidence level at least \code{conf_level}, that is,
#'   \eqn{\lceil C \rceil}.
#'
#' @references
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#' Francis Bilson Darku \email{fbilsond@@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # Mean score on a 100 point skills assessment to within 2 points, with
#' # 95% confidence, if the standard deviation were known to be 8.
#' seq_fwd_norm_mean_optimal(variance = 8^2, d = 2)
#'
#' # The same target at 99% confidence.
#' seq_fwd_norm_mean_optimal(variance = 8^2, d = 2, conf_level = 0.99)
#'
#' @seealso \code{\link{seq_fwd_norm_mean}} and
#'   \code{\link{ts_fwd_norm_mean}}, which reach the same target without
#'   knowing the variance; \code{\link{seq_budget_optimal}} for the
#'   corresponding known parameter allocation in the budget constrained
#'   procedure.
#' @family seq_fwd_norm_mean
#' @export
seq_fwd_norm_mean_optimal <- function(variance, d, conf_level = 0.95) {
  if (missing(variance)) {
    .smsd_stop("'variance' must be supplied: the population variance (the ",
               "square of the standard deviation).")
  }
  variance <- .check_positive(variance, "variance")
  .ts_check_width(d)
  conf_level <- .check_conf_level(conf_level)
  sigma <- sqrt(variance)

  # Does the known variance interval reach the level at sample size n? The
  # coverage formula of the help page, as it stands.
  at_level <- function(n) {
    if (1 - conf_level < 16 * .Machine$double.eps) {
      2 * stats::pnorm(d * sqrt(n) / sigma, lower.tail = FALSE) <=
        1 - conf_level
    } else {
      2 * stats::pnorm(d * sqrt(n) / sigma) - 1 >= conf_level
    }
  }

  # C = z^2 sigma^2 / d^2, computed as the square of z sigma / d so that
  # neither a very large nor a very small scale of measurement overflows or
  # underflows on the way.
  z <- .aipe_z(conf_level)
  advice <- paste0("Check that 'd' is in the units of the data and ",
                   "'variance' in their square; the target as given calls ",
                   "for more than two billion observations.")
  n <- .as_sample_size((z * sigma / d)^2, advice = advice)
  # C and the coverage it inverts can disagree by one when C is within
  # rounding error of a whole number. The coverage decides, so the result is
  # the smallest n at which the coverage formula is at least conf_level. The
  # second branch is also what makes the result at least 1: C rounds to 0
  # when it underflows (a target far wider than the spread of the data), and
  # no sample size below 1 is at the level.
  if (n > 1L && at_level(n - 1L)) {
    n <- n - 1L
  } else if (!at_level(n)) {
    n <- .as_sample_size(n + 1, advice = advice)
  }
  n
}
