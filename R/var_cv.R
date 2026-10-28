# Estimator of the asymptotic variance xi^2 of the sample coefficient of
# variation k_n = s_n / xbar_n: sqrt(n) (k_n - kappa) is asymptotically
# N(0, xi^2), with
#
#   xi^2 = mu_4 / (4 mu^2 sigma^2) - sigma^2 / (4 mu^2) - mu_3 / mu^3
#          + sigma^4 / mu^4
#
# (Kelley, Darku, and Chattopadhyay, 2018, Equation 30). The sample version
# (their Equation 32) substitutes the sample mean and variance and the
# unbiased central moment estimators .mu3() and .mu4() (their Equations 78
# and 79). Used by seq_aipe_cv() and by seq_mr_cv(). Internal; not exported.
#
# x: numeric vector (at least 4 values) with a positive mean.
# na.rm: remove missing values first? Default FALSE.
# floor: apply the lower bound n^-3? Default TRUE.
# Returns a single number: max{V_n^2, n^-3} (their Equation 31), or, with
# floor = FALSE, V_n^2 itself, which can be zero or negative, and is NaN for
# a sample with no variation (see .is_constant()). A value whose terms
# cancel to ten significant digits is rounding noise and is returned as
# exactly zero (see .aipe_zero_if_cancelled()); every other value is the
# number the formula gives.
#
# This is the distribution free variance. It differs, for data that are not
# normal, from Equations 9 and 22 of Chattopadhyay and Kelley (2016), which
# agree with it only under normality.
#
# Authors: Bhargab Chattopadhyay, Ken Kelley.
#' @noRd
.var_cv <- function(x, na.rm = FALSE, floor = TRUE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  if (.is_constant(x)) {
    # No variation: the ratio mu4_hat / s^2 is 0 / 0.
    return(.var_floor(NaN, n, floor))
  }
  # The estimate does not depend on the scale of the data; see
  # .aipe_scale().
  x <- x / .aipe_scale(x)
  mx <- mean(x)
  vx <- stats::var(x)
  t1 <- .mu4(x) / (4 * mx^2 * vx)
  t2 <- vx / (4 * mx^2)
  t3 <- .mu3(x) / mx^3
  t4 <- vx^2 / (mx^4)
  vn2 <- .aipe_zero_if_cancelled(t1 - t2 - t3 + t4, c(t1, t2, t3, t4))
  .var_floor(vn2, n, floor)
}
