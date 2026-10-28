# Estimator of the asymptotic variance xi^2 of the sample Gini index
# G_n = Delta_n / (2 xbar_n): sqrt(n) (G_n - G_F) is asymptotically
# N(0, xi^2), with
#
#   xi^2 = Delta^2 sigma^2 / (4 mu^4) - Delta tau / mu^3 + Delta^2 / mu^2
#          + sigma_1^2 / mu^2
#
# (Chattopadhyay and De, 2016, Equation 6), where Delta = E|X_1 - X_2| is
# Gini's mean difference, tau = E(X_1 |X_1 - X_2|), and sigma_1^2 is the
# variance of E(|X_1 - X_2| | X_1). Used by seq_aipe_gini(). Internal; not
# exported.
#
# x: numeric vector of nonnegative values with a positive mean (at least 4).
# na.rm: remove missing values first? Default FALSE.
# floor: apply the lower bound n^-3? Default TRUE.
# Returns a single number: max{V_n^2, n^-3}, or, with floor = FALSE, V_n^2
# itself, which is zero for a sample with no variation (see .is_constant()).
#
# The estimator is their Equation 11,
#
#   V_n^2 = Delta_n^2 s_n^2 / (4 xbar^4) - Delta_n tau_n / xbar^3
#           + Delta_n^2 / xbar^2 + s_wn^2 / (4 xbar^2),
#
# where Delta_n is the sample Gini mean difference, tau_n is the U-statistic
# with kernel (1/2)(x_1 + x_2)|x_1 - x_2| (their Equation 10), and s_wn^2 is
# the sample variance of W_jn = n Delta_n - (n - 2) Delta_n^(j), with
# Delta_n^(j) the Gini mean difference of the sample without observation j
# (their Equation 9; s_wn^2 estimates 4 sigma_1^2). Their Equation 11
# truncates V_n^2 below at zero; here, as in the rest of the family, a
# positive value is bounded below by n^-3, and a value that is zero or
# negative is treated by seq_aipe_gini() as not yet usable. For a sample
# with no variation every term is zero, and zero is returned without
# computing them: arithmetic on values that are equal on paper leaves
# noise of the order of 1e-30, which would otherwise pass as a positive
# estimate.
#
# Computation. Everything is obtained from the order statistics y_(1) <=
# ... <= y_(n) in O(n log n):
#
#   Delta_n = 2 sum_j (2j - n - 1) y_(j) / (n (n - 1)),
#   tau_n   = sum_j (2j - n - 1) y_(j)^2 / (n (n - 1)),
#   W_(j)n  = 2 { (2j - n - 1) y_(j) + T - C_j - C_(j-1) } / (n - 1),
#
# with C_j the cumulative sums and T the total. The last line uses the
# identity n Delta_n - (n - 2) Delta_n^(j) = (2 / (n - 1)) sum_i |x_j - x_i|,
# so no leave one out loop is needed. Delta_n and W are location invariant
# and are computed from the centered order statistics.
#
# Authors: Bhargab Chattopadhyay, Ken Kelley.
#' @noRd
.var_gini <- function(x, na.rm = FALSE, floor = TRUE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  n <- as.numeric(length(x))
  if (.is_constant(x)) {
    return(.var_floor(0, n, floor))
  }
  # The estimate does not depend on the scale of the data; see
  # .aipe_scale().
  y <- sort(x)
  y <- y / .aipe_scale(y)
  xbar <- mean(y)
  d <- .aipe_center(y)
  rank_weight <- 2 * seq_along(y) - n - 1

  delta <- 2 * sum(rank_weight * d) / (n * (n - 1))
  tau <- sum(rank_weight * y^2) / (n * (n - 1))
  cum <- cumsum(d)
  w <- 2 * (rank_weight * d + cum[length(cum)] - cum -
              c(0, cum[-length(cum)])) / (n - 1)
  s_w2 <- stats::var(w)

  vn2 <- delta^2 * stats::var(y) / (4 * xbar^4) - delta * tau / xbar^3 +
    delta^2 / xbar^2 + s_w2 / (4 * xbar^2)
  .var_floor(vn2, n, floor)
}
