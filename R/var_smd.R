# Estimator of the asymptotic variance xi^2 of the standardized mean
# difference for two independent groups of equal size: sqrt(n) (d_n - delta)
# is asymptotically N(0, xi^2), with n the per group sample size. Used by
# seq_aipe_smd() in the sequential accuracy in parameter estimation stopping
# rule of Kelley, Darku, and Chattopadhyay (2018). Internal; not exported.
#
# x, y: numeric vectors, the two groups (at least 4 values each).
# na.rm: remove missing values from each group first? Default FALSE.
# floor: apply the lower bound n^-3? Default TRUE.
# Returns a single number: max{V_n^2, n^-3} (their Equation 27), with n the
# length of x, or, with floor = FALSE, V_n^2 itself, which can be zero or
# negative, and is NaN when either group has no variation (see
# .is_constant()). A value whose terms cancel to ten significant digits is
# rounding noise and is returned as exactly zero (see
# .aipe_zero_if_cancelled()); every other value is the number the formula
# gives.
#
# One group with no variation. The formula has a value then (the pooled
# variance is half the variance of the other group), but the variance of
# the constant group has not been estimated, and the value understates the
# sampling variance: the pooled standard deviation is too small and the
# standardized mean difference too large. NaN is returned, which
# seq_aipe_smd() reports as a variance estimate that is not yet usable;
# seq_aipe_c() and seq_budget() read the same data the same way. (The test
# "one group with no variation makes the variance not yet usable" in
# tests/testthat/test-seq_aipe_smd.R has an example worked by hand; the
# help of seq_aipe_smd() gives the coverage and the sample size at stake.)
#
# Derivation. With theta_1 = mu_1 - mu_2, theta_2 = sigma^2 (the common
# variance), g_1(u) = u and g_2(u) = sqrt(u), delta = g_1(theta_1) /
# g_2(theta_2), and the delta method (their Theorem 1, Equation 55, and
# Theorem 3) gives xi^2 = D' Sigma D with
#
#   D     = ( 1 / sigma,  -(mu_1 - mu_2) / (2 sigma^3) )'
#   Sigma = | 2 sigma^2                (mu_13 - mu_23) / 2            |
#           | (mu_13 - mu_23) / 2      (mu_14 + mu_24 - 2 sigma^4) / 4 |
#
# where mu_kj is the j-th central moment of group k. Expanding,
#
#   xi^2 = 2 - t_1 + t_2 t_3,
#   t_1  = (1/2) (mu_1 - mu_2) (mu_13 - mu_23) / sigma^4,
#   t_2  = (mu_1 - mu_2)^2 / (4 sigma^6),
#   t_3  = (mu_14 + mu_24) / 4 - sigma^4 / 2.
#
# Correction to Equation 26. Equation 26 of the article (p. 232), and the
# sample version printed under its Equation 27, omit the factor 1/2 in t_1.
# The factor follows from the matrix above (the cross term of D' Sigma D is
# 2 * (1 / sigma) * (-(mu_1 - mu_2) / (2 sigma^3)) * (mu_13 - mu_23) / 2).
# The corrected form is implemented here. The authors' correction is at
# https://kenkelley.org/wp-content/uploads/articles/Kelley_Darku_Chattopadhyay_Equation_26.pdf
# and ships with the package as inst/eq26_correction.pdf, with the algebra
# written out in inst/eq26_correction.md. tests/testthat/test-var_smd.R pins
# the corrected form against the matrix D' Sigma D and against the sampling
# variance of d in a case (groups skewed in opposite directions) where the
# two forms differ by a factor of about 2.5.
#
# The sample version substitutes the sample means, the pooled variance
# s_p^2 = (s_1^2 + s_2^2) / 2, and the unbiased central moment estimators
# .mu3() and .mu4() (their Equations 75 and 76).
#
# Authors: Bhargab Chattopadhyay, Ken Kelley.
#' @noRd
.var_smd <- function(x, y, na.rm = FALSE, floor = TRUE) {
  if (na.rm) {
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
  }
  n <- length(x)
  if (.is_constant(x) || .is_constant(y)) {
    # No variation in a group: its variance has not been estimated (and
    # with none in either group the pooled variance is zero and the
    # estimate is undefined).
    return(.var_floor(NaN, n, floor))
  }
  # The estimate does not depend on the common scale of the two groups;
  # see .aipe_scale().
  scale <- .aipe_scale(x, y)
  x <- x / scale
  y <- y / scale

  s_p <- sqrt(0.5 * (stats::var(x) + stats::var(y)))
  mean_diff <- mean(x) - mean(y)

  t_1 <- 0.5 * mean_diff * (.mu3(x) - .mu3(y)) / s_p^4
  t_2 <- mean_diff^2 / (4 * s_p^6)
  mu4_sum <- .mu4(x) + .mu4(y)
  t_3 <- mu4_sum / 4 - s_p^4 / 2

  vn2 <- .aipe_zero_if_cancelled(
    2 - t_1 + t_2 * t_3,
    c(2, t_1, t_2 * mu4_sum / 4, t_2 * s_p^4 / 2)
  )
  .var_floor(vn2, n, floor)
}
