# Estimator of the asymptotic variance xi^2 of the sample standardized mean
# xbar_n / s_n: sqrt(n) (xbar_n / s_n - mu / sigma) is asymptotically
# N(0, xi^2), with
#
#   xi^2 = 1 - mu mu_3 / sigma^4 + mu^2 (mu_4 - sigma^4) / (4 sigma^6).
#
# The standardized mean is Example 3 of Kelley, Darku, and Chattopadhyay
# (2018, p. 229). The article does not print this variance; it follows from
# the article's Theorem 1 (Equation 55), xi^2 = D' Sigma D, with
# theta_1 = mu, theta_2 = sigma^2, g_1(u) = u, g_2(u) = sqrt(u), so that
# D = (1 / sigma, -mu / (2 sigma^3))' and Sigma has entries sigma^2, mu_3,
# and mu_4 - sigma^4. Under normality it is 1 + (mu / sigma)^2 / 2. The
# sample version substitutes the sample mean and variance and the unbiased
# central moment estimators .mu3() and .mu4(). Used by seq_aipe_sm().
# Internal; not exported.
#
# x: numeric vector (at least 4 values).
# na.rm: remove missing values first? Default FALSE.
# floor: apply the lower bound n^-3? Default TRUE.
# Returns a single number: max{V_n^2, n^-3}, or, with floor = FALSE, V_n^2
# itself, which can be zero or negative, and is NaN for a sample with no
# variation (see .is_constant()). A value whose terms cancel to ten
# significant digits is rounding noise and is returned as exactly zero (see
# .aipe_zero_if_cancelled()).
#
# Authors: Bhargab Chattopadhyay, Ken Kelley.
#' @noRd
.var_sm <- function(x, na.rm = FALSE, floor = TRUE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  if (.is_constant(x)) {
    # No variation: the standardized mean and its variance are undefined.
    return(.var_floor(NaN, n, floor))
  }
  # The estimate does not depend on the scale of the data; see
  # .aipe_scale().
  x <- x / .aipe_scale(x)
  xbar <- mean(x)
  sn2 <- stats::var(x)
  t1 <- xbar^2 / (4 * sn2)
  t2 <- xbar * .mu3(x) / sn2^2
  t3 <- xbar^2 * .mu4(x) / (4 * sn2^3)
  vn2 <- .aipe_zero_if_cancelled(1 - t1 - t2 + t3, c(1, t1, t2, t3))
  .var_floor(vn2, n, floor)
}
