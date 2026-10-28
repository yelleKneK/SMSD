# Estimator of the asymptotic variance xi^2 of the least squares slope
# b_1n = s_XY / s_X^2 of the simple linear regression of y on x:
# sqrt(n) (b_1n - beta_1) is asymptotically N(0, xi^2), with
#
#   xi^2 = mu_22 / sigma_X^4 - 2 sigma_XY mu_31 / sigma_X^6
#          + sigma_XY^2 mu_40 / sigma_X^8
#
# (Kelley, Darku, and Chattopadhyay, 2018, Equation 37 and Theorem 6,
# Equation 81), where mu_ij = E[(X - mu_X)^i (Y - mu_Y)^j]. It makes no
# assumption of normal or homoscedastic errors. The sample version (their
# Equation 85) substitutes the sample variance and covariance, .mu22(),
# .mu31(), and .mu4(). Used by seq_aipe_reg_coef(). Internal; not exported.
#
# x: the predictor; y: the response; paired, at least 4 pairs.
# na.rm: drop a pair when either member is missing? Default FALSE.
# floor: apply the lower bound n^-3? Default TRUE.
# Returns a single number: max{V_n^2, n^-3} (their Equation 84), or, with
# floor = FALSE, V_n^2 itself (their Equation 85).
#
# Being a difference of moment estimates, the raw value can be negative,
# and at small n it often is; see the help of seq_aipe_reg_coef() for how
# often and what the check function does about it. Two cases are settled
# without computing it (see .is_constant() for what counts as no
# variation). A predictor with no variation has no slope, and the raw value
# is NaN. A response with no variation has a slope of zero and a variance
# estimate of zero (every term carries a moment of y), and zero is
# returned: arithmetic on values that are equal on paper leaves noise of
# the order of 1e-30, which would otherwise pass as a positive estimate.
# For the same reason a value whose three terms cancel to ten significant
# digits is returned as zero (.aipe_zero_if_cancelled()): the terms cancel
# exactly on paper for collinear pairs whose sample kurtosis is
# (n^2 - 3) / (n - 1)^2, which discrete data can produce, and the sign of
# what is computed then is the sign of rounding noise.
#
# The slope has units (those of y per unit of x) and so does xi^2 (their
# square). The estimate is computed from x and y each divided by a power of
# two of its own magnitude and then returned to the units of the data, so
# that the powers of the predictor's variance stay inside the range of a
# double whatever those units are; see .aipe_scale(). A raw value that is
# positive but below n^-3 is replaced by the bound, as for the other effect
# sizes, and seq_aipe_reg_coef() adds a note (see its help).
#
# Authors: Bhargab Chattopadhyay, Ken Kelley.
#' @noRd
.var_reg_coef <- function(x, y, na.rm = FALSE, floor = TRUE) {
  if (na.rm) {
    ok <- !is.na(x) & !is.na(y)
    x <- x[ok]
    y <- y[ok]
  }
  n <- length(x)
  if (.is_constant(x)) {
    return(.var_floor(NaN, n, floor))
  }
  if (.is_constant(y)) {
    return(.var_floor(0, n, floor))
  }
  scale_x <- .aipe_scale(x)
  scale_y <- .aipe_scale(y)
  x <- x / scale_x
  y <- y / scale_y
  vx <- stats::var(x)
  cxy <- stats::cov(x, y)
  t1 <- .mu22(x, y) / vx^2
  t2 <- 2 * cxy * .mu31(x, y) / vx^3
  t3 <- cxy^2 * .mu4(x) / vx^4
  vn2 <- .aipe_zero_if_cancelled(t1 - t2 + t3, c(t1, t2, t3))
  .var_floor(vn2 * (scale_y / scale_x)^2, n, floor)
}
