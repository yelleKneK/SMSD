# Unbiased estimator of the fourth central moment (the h-statistic of order
# four). Kelley, Darku, and Chattopadhyay (2018, Equations 76 and 79; 2019,
# Equations 71 and 72). Internal; not exported.
#
# x: numeric vector with at least 4 (non missing) values.
# na.rm: remove missing values first? Default FALSE.
#
# The articles print the estimator in terms of raw power sums,
#
#   { n^2 sum((x - xbar)^4) - (2n - 3) sum(x^4) + (8n - 12) xbar sum(x^3)
#     - (6 - 9/n) (sum(x^2))^2 } / ((n - 1)(n - 2)(n - 3)).
#
# That expression is location invariant, so it may be evaluated on the
# centered data d = x - xbar, where the term in sum(d^3) drops out (it is
# multiplied by the mean of d, which is zero) and the rest collapses to
#
#   { (n^2 - 2n + 3) sum(d^4) - 3 (2n - 3) / n (sum(d^2))^2 }
#     / ((n - 1)(n - 2)(n - 3)).
#
# The two forms are algebraically identical; the centered one is what is
# computed, because the raw power sums lose all precision when the mean is
# large relative to the standard deviation (for standard normal data shifted
# by 1e6 the raw form returned 0 where the centered form returns 2.17). The
# estimator is exactly unbiased for any distribution with a finite fourth
# moment.
#
# Authors: Ken Kelley, Francis Bilson Darku, Bhargab Chattopadhyay.
#' @noRd
.mu4 <- function(x, na.rm = FALSE) {
  if (!is.numeric(x)) {
    .smsd_stop("'x' must be a numeric vector.")
  }
  .check_flag(na.rm, "na.rm")
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  if (length(x) < 4) {
    .smsd_stop("'x' must contain at least 4 observations.")
  }
  n <- as.numeric(length(x))
  d <- .aipe_center(x)
  ((n^2 - 2 * n + 3) * sum(d^4) - 3 * (2 * n - 3) / n * sum(d^2)^2) /
    ((n - 1) * (n - 2) * (n - 3))
}
