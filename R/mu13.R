# Estimator of the joint central moment mu_13 = E[(X - mu_X) (Y - mu_Y)^3]
# used by the Pearson correlation variance estimator:
#
#   mu13_hat = k13 + 3 k02 k11,
#
# Kelley, Darku, and Chattopadhyay (2019, Equation 73), where k02 is the
# sample variance of y, k11 the sample covariance, and k13 the bivariate
# k-statistic (the unbiased estimator of the joint cumulant kappa_13; their
# Equation 74). Internal; not exported.
#
# x, y: paired numeric vectors with at least 4 pairs.
# na.rm: drop a pair when either member is missing? Default FALSE.
#
# As in .mu22(), the k-statistic is evaluated on the centered data, where
# Equation 74 reduces to
#
#   k13 = n { (n + 1) S13 - 3 (n - 1)/n S11 S02 } / ((n - 1)(n - 2)(n - 3)),
#
# with S_pq the sums of products of the centered values; the raw power sum
# form loses all precision when a mean is large relative to its standard
# deviation. k13 is exactly unbiased for kappa_13; mu13_hat is consistent
# but not unbiased.
#
# Authors: Ken Kelley, Francis Bilson Darku, Bhargab Chattopadhyay.
#' @noRd
.mu13 <- function(x, y, na.rm = FALSE) {
  if (na.rm) {
    ok <- !is.na(x) & !is.na(y)
    x <- x[ok]
    y <- y[ok]
  }
  n <- as.numeric(length(x))
  dx <- .aipe_center(x)
  dy <- .aipe_center(y)
  s02 <- sum(dy^2)
  s11 <- sum(dx * dy)
  s13 <- sum(dx * dy^3)

  k02 <- s02 / (n - 1)
  k11 <- s11 / (n - 1)
  k13 <- n / ((n - 1) * (n - 2) * (n - 3)) *
    ((n + 1) * s13 - 3 * (n - 1) / n * s11 * s02)

  k13 + 3 * k02 * k11
}
