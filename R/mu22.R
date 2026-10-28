# Estimator of the joint central moment mu_22 = E[(X - mu_X)^2 (Y - mu_Y)^2]
# used by the Pearson correlation and regression slope variance estimators:
#
#   mu22_hat = k22 + k20 k02 + 2 k11^2,
#
# Kelley, Darku, and Chattopadhyay (2019, Equation 73), where k20 and k02 are
# the sample variances, k11 is the sample covariance, and k22 is the
# bivariate k-statistic (the unbiased estimator of the joint cumulant
# kappa_22; their Equation 74, after Cook, 1951, and Fisher, 1930).
# Internal; not exported.
#
# x, y: paired numeric vectors with at least 4 pairs.
# na.rm: drop a pair when either member is missing? Default FALSE.
#
# Equation 74 is printed in raw power sums W_pq = sum(x^p y^q). The
# k-statistics are location invariant, so they are evaluated here on the
# centered data, where every term carrying W_10 or W_01 vanishes and
#
#   k22 = n { (n + 1) S22 - (n - 1)/n S20 S02 - 2 (n - 1)/n S11^2 }
#           / ((n - 1)(n - 2)(n - 3)),
#
# with S_pq the sums of products of the centered values. The raw form loses
# all precision when a mean is large relative to its standard deviation.
# (As printed, Equation 74 gives the coefficient of W_12 W_10 in k22 as
# (n + 1)/n; symmetry in x and y, and location invariance, require
# 2 (n + 1)/n, the coefficient of W_21 W_01. The term vanishes for centered
# data, so the misprint has no bearing on the value computed here.)
#
# k22 is exactly unbiased for kappa_22; mu22_hat is consistent but not
# unbiased, because the products k20 k02 and k11^2 are biased for the
# products of the moments they estimate.
#
# Authors: Ken Kelley, Francis Bilson Darku, Bhargab Chattopadhyay.
#' @noRd
.mu22 <- function(x, y, na.rm = FALSE) {
  if (na.rm) {
    ok <- !is.na(x) & !is.na(y)
    x <- x[ok]
    y <- y[ok]
  }
  n <- as.numeric(length(x))
  dx <- .aipe_center(x)
  dy <- .aipe_center(y)
  s20 <- sum(dx^2)
  s02 <- sum(dy^2)
  s11 <- sum(dx * dy)
  s22 <- sum(dx^2 * dy^2)

  k20 <- s20 / (n - 1)
  k02 <- s02 / (n - 1)
  k11 <- s11 / (n - 1)
  k22 <- n / ((n - 1) * (n - 2) * (n - 3)) *
    ((n + 1) * s22 - (n - 1) / n * s20 * s02 - 2 * (n - 1) / n * s11^2)

  k22 + k20 * k02 + 2 * k11^2
}
