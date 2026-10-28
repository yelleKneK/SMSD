# Value tests for the internal asymptotic variance estimators (.var_r,
# .var_cv, .var_sm, .var_reg_coef, .var_c, .var_gini; .var_smd has its own
# file). Each estimator is checked against (a) the formula as printed in its
# source, transcribed literally here (loops and raw power sums, on data
# where that route is exact), and (b) a population value known in closed
# form or reported in the source, by large sample evaluation. Tolerances for
# the large sample checks are at least four standard deviations of the
# estimate, measured over 40 to 100 replications and stated with each test.

# Kelley, Darku, and Chattopadhyay (2019), Equation 70: the Pearson
# estimator, with the moment estimators of Equations 71 to 74 in raw power
# sums.
pearson_as_printed <- function(x, y) {
  kx <- k_as_printed(x, y)
  ky <- k_as_printed(y, x)
  k20 <- var(x); k02 <- var(y); k11 <- cov(x, y)
  m22 <- kx[["k22"]] + k20 * k02 + 2 * k11^2
  m13 <- kx[["k13"]] + 3 * k02 * k11
  m31 <- ky[["k13"]] + 3 * k20 * k11
  cor(x, y)^2 / 4 *
    (mu4_as_printed(x) / k20^2 + mu4_as_printed(y) / k02^2 +
       2 * m22 / (k20 * k02) + 4 * m22 / k11^2 -
       4 * m31 / (k11 * k20) - 4 * m13 / (k11 * k02))
}

# Kelley, Darku, and Chattopadhyay (2018), Equation 85: the slope estimator.
slope_as_printed <- function(x, y) {
  kx <- k_as_printed(x, y)
  ky <- k_as_printed(y, x)
  k20 <- var(x); k02 <- var(y); k11 <- cov(x, y)
  m22 <- kx[["k22"]] + k20 * k02 + 2 * k11^2
  m31 <- ky[["k13"]] + 3 * k20 * k11
  m22 / k20^2 - 2 * k11 * m31 / k20^3 + k11^2 * mu4_as_printed(x) / k20^4
}

# Kelley, Darku, and Chattopadhyay (2019), Equations 27 to 29 (Kendall) and
# 51 to 53 (Spearman), as loops over the ranks.
kendall_as_printed <- function(x, y) {
  n <- length(x); rx <- rank(x); ry <- rank(y); w <- numeric(n)
  for (i in seq_len(n)) {
    w[i] <- 2 / n * sum(rx <= rx[i] & ry <= ry[i]) -
      rx[i] / (n + 1) - ry[i] / (n + 1)
  }
  16 / (n - 1) * sum((w - mean(w))^2)
}
spearman_as_printed <- function(x, y) {
  n <- length(x); rx <- rank(x); ry <- rank(y); zz <- numeric(n)
  for (i in seq_len(n)) {
    zz[i] <- rx[i] / (n + 1) * ry[i] / (n + 1) +
      mean((rx[i] <= rx) * ry / (n + 1)) +
      mean((ry[i] <= ry) * rx / (n + 1))
  }
  144 / (n - 1) * sum((zz - mean(zz))^2)
}

# Chattopadhyay and De (2016), Equations 9 to 11, by brute force over pairs
# and with an explicit leave one out loop.
gini_var_as_printed <- function(x) {
  n <- length(x)
  pairwise <- function(v, kernel) {
    p <- utils::combn(length(v), 2)
    mean(kernel(v[p[1, ]], v[p[2, ]]))
  }
  delta <- pairwise(x, function(a, b) abs(a - b))
  tau <- pairwise(x, function(a, b) 0.5 * (a + b) * abs(a - b))
  delta_without <- vapply(seq_len(n), function(j) {
    pairwise(x[-j], function(a, b) abs(a - b))
  }, numeric(1))
  w <- n * delta - (n - 2) * delta_without
  xbar <- mean(x)
  delta^2 * var(x) / (4 * xbar^4) - delta * tau / xbar^3 +
    delta^2 / xbar^2 + var(w) / (4 * xbar^2)
}

# Kelley, Darku, and Chattopadhyay (2018), Equation 32: the coefficient of
# variation estimator, with the third central moment estimator of their
# Equation 78 and the fourth of Equation 79.
mu3_as_printed <- function(x) {
  n <- length(x)
  n / ((n - 1) * (n - 2)) * sum((x - mean(x))^3)
}
cv_var_as_printed <- function(x) {
  xbar <- mean(x); s2 <- var(x)
  mu4_as_printed(x) / (4 * xbar^2 * s2) - s2 / (4 * xbar^2) -
    mu3_as_printed(x) / xbar^3 + s2^2 / xbar^4
}
# The standardized mean: 1 - mu mu3 / sigma^4 + mu^2 (mu4 - sigma^4) /
# (4 sigma^6) (the delta method variance; see R/var_sm.R), with the same
# moment estimators.
sm_var_as_printed <- function(x) {
  xbar <- mean(x); s2 <- var(x)
  1 - xbar * mu3_as_printed(x) / s2^2 +
    xbar^2 * (mu4_as_printed(x) - s2^2) / (4 * s2^3)
}

# Bisection for the point at which an estimate changes sign along a path of
# samples: f(0) < 0 < f(1). It stops at the first t where f is exactly zero
# and returns c(t, t); if f is never exactly zero it ends with two adjacent
# doubles, f negative at the first and positive at the second. With real
# arithmetic the estimate is zero at a single t; in floating point the
# values computed next to that t are rounding noise of either sign, and an
# estimator that hands the noise back never returns an exact zero.
sign_change <- function(f) {
  lo <- 0
  hi <- 1
  repeat {
    mid <- (lo + hi) / 2
    if (mid <= lo || mid >= hi) return(c(lo, hi))
    value <- f(mid)
    if (value == 0) return(c(mid, mid))
    if (value < 0) lo <- mid else hi <- mid
  }
}

test_that(".var_r (Pearson) equals Equation 70 as printed, bounded below as in Equation 8", {
  # Integer data, on which the raw power sums of the printed formula are
  # exact. In samples this small the printed value is often zero or
  # negative, and the estimator then returns the bound n^-3 of the
  # article's Equation 8, so the comparison is with max{printed, n^-3}. The
  # count at the end makes sure that the comparison was not only between
  # bounds (over 300 other seeds at least 8 of the 18 samples were above the
  # bound). A sample covariance of exactly zero, where the printed formula
  # is 0 / 0, has its own test below.
  withr::local_seed(113)
  above_bound <- 0
  for (n in c(8, 15, 40)) {
    for (i in 1:6) {
      x <- sample(-6:9, n, replace = TRUE)
      y <- x + sample(-5:5, n, replace = TRUE)
      if (cov(x, y) == 0) next
      printed <- pearson_as_printed(x, y)
      expect_equal(.var_r(x, y, method = "pearson"), max(printed, n^(-3)),
                   tolerance = 1e-9)
      above_bound <- above_bound + (printed > n^(-3))
    }
  }
  expect_gte(above_bound, 4)
})

test_that(".var_r (Pearson) has its limiting value when the sample covariance is zero", {
  # x = 1, ..., 5 twice and y below: the deviations of x are -2, -1, 0, 1,
  # 2 (twice), those of y are 0, -2, 2, -2, 0, 0, 2, -2, 2, 0, and the sum
  # of their products is 0, so r = 0 and Equation 70 as printed is 0 / 0.
  # With r^2 multiplied through, every term carrying r vanishes and
  # V^2 = mu22_hat / (s_x^2 s_y^2), where mu22_hat = k22 + k20 k02 (the
  # covariance k11 being zero). By hand, with n = 10: S20 = 20, S02 = 24,
  # S22 = 16, so k20 k02 = (20 / 9)(24 / 9) = 480 / 81 and
  #   k22 = 10 (11 * 16 - (9 / 10) 20 * 24) / (9 * 8 * 7) = -2560 / 504,
  # giving V^2 = 1 - (2560 / 504) / (480 / 81) = 1 / 7.
  x <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
  y <- c(3, 1, 5, 1, 3, 3, 5, 1, 5, 3)
  expect_identical(cov(x, y), 0)
  expect_true(is.nan(pearson_as_printed(x, y)))
  expect_equal(.var_r(x, y, method = "pearson"), 1 / 7)
  # The limiting value is approached continuously: nudging one value makes
  # the covariance nonzero and the printed formula usable, and the two
  # agree.
  y2 <- replace(y, 1, 3.001)
  expect_equal(.var_r(x, y2, method = "pearson"), pearson_as_printed(x, y2),
               tolerance = 1e-6)
  expect_equal(.var_r(x, y2, method = "pearson"), 1 / 7, tolerance = 0.01)
})

test_that(".var_r (Pearson) approaches (1 - rho^2)^2 for bivariate normal data", {
  # Their Equation 3 reduces to (1 - rho^2)^2 under bivariate normality.
  # Standard deviation of the estimate at n = 6,000: 0.016.
  withr::local_seed(113)
  d <- aipe_pairs(6000, rho = 0.6)
  expect_lt(abs(.var_r(d$x, d$y, method = "pearson") - (1 - 0.6^2)^2), 0.07)
})

test_that(".var_r (Kendall) equals Equations 27 to 29 as printed", {
  withr::local_seed(113)
  for (n in c(6, 25, 80)) {
    d <- aipe_pairs(n, rho = 0.5)
    expect_equal(.var_r(d$x, d$y, method = "kendall"),
                 kendall_as_printed(d$x, d$y), tolerance = 1e-12)
  }
})

test_that(".var_r (Kendall) approaches the bivariate normal value", {
  # For bivariate normal data the asymptotic variance of Kendall's tau is
  # 4 (1/9 - ((2/pi) asin(rho/2))^2). That closed form reproduces the
  # theoretical sample sizes n_KT printed in Table 8 of Kelley, Darku, and
  # Chattopadhyay (2019) for a 95% interval: 524 and 131 at rho = 0.5, and
  # 677 and 170 at rho = 0.1, for widths 0.1 and 0.2.
  xi2 <- function(rho) 4 * (1 / 9 - (2 / pi * asin(rho / 2))^2)
  z <- qnorm(0.975)
  expect_equal(ceiling((2 * z / c(0.1, 0.2))^2 * xi2(0.5)), c(524, 131))
  expect_equal(ceiling((2 * z / c(0.1, 0.2))^2 * xi2(0.1)), c(677, 170))
  # Standard deviation of the estimate at n = 1,500: 0.0137 at independence
  # (4,000 samples, mean 0.0007 above 4/9) and 0.0160 at rho = 0.5 (2,000
  # samples, mean 0.0008 above the value), so each tolerance, less that
  # offset, is 4.3 of them. (The samples were drawn as here, the 4,000
  # first, after set.seed(20261012).)
  withr::local_seed(113)
  d <- aipe_pairs(1500, rho = 0.5)
  expect_lt(abs(.var_r(d$x, d$y, method = "kendall") - xi2(0.5)), 0.07)
  # At independence the value is 4/9.
  expect_lt(abs(.var_r(rnorm(1500), rnorm(1500), method = "kendall") - 4 / 9),
            0.06)
})

test_that(".var_r (Spearman) equals Equations 51 to 53 as printed", {
  withr::local_seed(113)
  for (n in c(6, 25, 80)) {
    d <- aipe_pairs(n, rho = 0.5)
    expect_equal(.var_r(d$x, d$y, method = "spearman"),
                 spearman_as_printed(d$x, d$y), tolerance = 1e-12)
  }
})

test_that(".var_r (Spearman) approaches 1 at independence", {
  # sqrt(n) r_s is asymptotically standard normal under independence.
  # Standard deviation of the estimate at n = 1,500: 0.0218 (200 samples),
  # so 0.09 is a little over four of them.
  withr::local_seed(113)
  expect_lt(abs(.var_r(rnorm(1500), rnorm(1500), method = "spearman") - 1),
            0.09)
})

test_that(".var_cv approaches Equation 30 of Kelley, Darku, and Chattopadhyay (2018)", {
  # xi^2 = mu4/(4 mu^2 sigma^2) - sigma^2/(4 mu^2) - mu3/mu^3 + sigma^4/mu^4.
  eq30 <- function(mu, s2, mu3, mu4) {
    mu4 / (4 * mu^2 * s2) - s2 / (4 * mu^2) - mu3 / mu^3 + s2^2 / mu^4
  }
  # Gamma with shape 4, scale 1: mean 4, variance 4, mu3 = 8, mu4 = 72.
  expect_equal(eq30(4, 4, 8, 72), 0.15625)
  withr::local_seed(113)
  # Standard deviation of the estimate at n = 20,000: 0.0073.
  expect_lt(abs(.var_cv(rgamma(20000, shape = 4)) - 0.15625), 0.032)
  # Normal with mean 4, sd 1 (kappa = 0.25): kappa^2 (kappa^2 + 1/2).
  expect_equal(eq30(4, 1, 0, 3), 0.25^2 * (0.25^2 + 0.5))
  # Standard deviation of the estimate at n = 20,000: 0.0007.
  expect_lt(abs(.var_cv(rnorm(20000, 4, 1)) - 0.25^2 * (0.25^2 + 0.5)), 0.003)
})

test_that(".var_sm approaches the delta method variance of the standardized mean", {
  # xi^2 = 1 - mu mu3 / sigma^4 + mu^2 (mu4 - sigma^4) / (4 sigma^6).
  # Normal, mu / sigma = 0.5: 1 + 0.25 / 2. Standard deviation of the
  # estimate at n = 20,000: 0.010.
  withr::local_seed(113)
  expect_lt(abs(.var_sm(rnorm(20000, mean = 0.5)) - 1.125), 0.045)
  # Standardized chi square (3 df) shifted to mean 0.4: mu3 = sqrt(8/3),
  # mu4 = 7, so xi^2 = 1 - 0.4 sqrt(8/3) + 0.16 * 6 / 4 = 0.5868. Standard
  # deviation of the estimate at n = 20,000: 0.0065 (200 samples), so 0.027
  # is a little over four of them.
  expect_lt(abs(.var_sm(aipe_skewed(20000) + 0.4) -
                  (1 - 0.4 * sqrt(8 / 3) + 0.16 * 6 / 4)), 0.027)
})

test_that(".var_reg_coef equals Equation 85 as printed, bounded below as in Equation 84", {
  # As for the Pearson estimator: the printed value is often zero or
  # negative in samples this small, so the comparison is with
  # max{printed, n^-3}, and the unbounded value is compared as well (over
  # 300 other seeds at least 10 of the 18 samples were above the bound).
  withr::local_seed(113)
  above_bound <- 0
  for (n in c(8, 15, 40)) {
    for (i in 1:6) {
      x <- sample(-6:9, n, replace = TRUE)
      y <- x + sample(-5:5, n, replace = TRUE)
      printed <- slope_as_printed(x, y)
      expect_equal(.var_reg_coef(x, y), max(printed, n^(-3)),
                   tolerance = 1e-9)
      expect_equal(.var_reg_coef(x, y, floor = FALSE), printed,
                   tolerance = 1e-9)
      above_bound <- above_bound + (printed > n^(-3))
    }
  }
  expect_gte(above_bound, 4)
})

test_that(".var_reg_coef approaches the distribution free variance of the slope", {
  # Equation 37 is E[(X - mu_X)^2 e^2] / sigma_X^4 for the model
  # y = b0 + b1 x + e with E(e | x) = 0. With x and e independent standard
  # normal it is 1. Standard deviation of the estimate at n = 20,000: 0.021.
  withr::local_seed(113)
  x <- rnorm(20000)
  expect_lt(abs(.var_reg_coef(x, 0.5 * x + rnorm(20000)) - 1), 0.09)
  # With heteroscedastic errors e = x * u (u standard normal, independent
  # of x) it is E[x^4] = 3, three times the value 1 that the usual constant
  # variance formula sigma_e^2 / sigma_X^2 would give. Standard deviation of
  # the estimate at n = 20,000: 0.11.
  expect_lt(abs(.var_reg_coef(x, 0.5 * x + x * rnorm(20000)) - 3), 0.5)
})

test_that(".var_c is the weighted sum of the group variances", {
  # Groups 1, 2, 3, 4 and 2, 4, 6, 8 have variances 5/3 and 20/3; with
  # weights 1 and -1 the contrast variance is 25/3.
  g <- list(c(1, 2, 3, 4), c(2, 4, 6, 8))
  expect_equal(.var_c(g, c(1, -1)), 25 / 3)
  expect_equal(.var_c(g, c(0.5, 2)), 0.25 * 5 / 3 + 4 * 20 / 3)
  # A matrix or data frame with one column per group gives the same value.
  expect_equal(.var_c(cbind(g[[1]], g[[2]]), c(1, -1)), 25 / 3)
  expect_equal(.var_c(data.frame(a = g[[1]], b = g[[2]]), c(1, -1)), 25 / 3)
  # Missing values are removed group by group on request.
  expect_equal(.var_c(list(c(1, 2, 3, 4, NA), c(2, 4, NA, 6, 8)), c(1, -1),
                      na.rm = TRUE), 25 / 3)
})

test_that(".var_c is bounded below by n^-3 with n the smallest group size", {
  expect_identical(.var_c(list(rep(1, 10), rep(2, 10)), c(1, -1)), 10^(-3))
  expect_identical(.var_c(list(rep(1, 10), rep(2, 6)), c(1, -1)), 6^(-3))
  # A group with no variation and a weight that is not zero: its variance
  # has not been estimated, so neither has the variance of the contrast.
  # There is no raw value (NaN), whether the other group counts or not, and
  # the bound is returned. (Summing what there is would give 0 with weights
  # 1 and 0, and 6, the variance of 1, ..., 8, with weights 1 and -1: the
  # contrast's variance with the constant group left out.)
  flat_weighted <- list(rep(1, 8), c(1, 2, 3, 4, 5, 6, 7, 8))
  expect_identical(.var_c(flat_weighted, c(1, 0)), 8^(-3))
  expect_true(is.nan(.var_c(flat_weighted, c(1, 0), floor = FALSE)))
  expect_true(is.nan(.var_c(flat_weighted, c(1, -1), floor = FALSE)))
  expect_identical(.var_c(flat_weighted, c(1, -1)), 8^(-3))
  # A group with no variation whose weight is zero takes no part in the
  # contrast. floor = FALSE then returns the weighted sum itself, whatever
  # its size: the variance of 1, ..., 8 is 6, so with weights 0 and 0.01 it
  # is 0.0006, below 8^-3 = 0.00195.
  expect_equal(.var_c(flat_weighted, c(0, 0.01), floor = FALSE), 0.0006)
  expect_identical(.var_c(flat_weighted, c(0, 0.01)), 8^(-3))
})

test_that(".var_gini equals Equations 9 to 11 of Chattopadhyay and De (2016)", {
  # Their Equation 11 bounds the value below at zero; the package bounds it
  # at n^-3, so the comparison is with max{printed, n^-3}. (The printed
  # value was above the bound in 96% of 3,000 samples of 5 and in every one
  # of 3,000 samples of 12 and of 30.)
  withr::local_seed(113)
  above_bound <- 0
  for (n in c(5, 12, 30)) {
    x <- rlnorm(n, meanlog = 0, sdlog = 0.7)
    printed <- gini_var_as_printed(x)
    expect_equal(.var_gini(x), max(printed, n^(-3)), tolerance = 1e-10)
    above_bound <- above_bound + (printed > n^(-3))
  }
  expect_gte(above_bound, 1)
  # With ties.
  x <- c(1, 1, 2, 2, 2, 5, 9, 9, 14)
  expect_equal(.var_gini(x), gini_var_as_printed(x), tolerance = 1e-10)
})

test_that(".var_gini approaches the population values", {
  withr::local_seed(113)
  # Exponential: G = 1/2 and, from their Equation 6 with Delta = 1,
  # tau = 3/2, and sigma_1^2 = 1/3, xi^2 = 1/3 + 1/4 - 3/2 + 1 = 1/12.
  # Standard deviation of the estimate at n = 20,000: 0.0015.
  expect_lt(abs(.var_gini(rexp(20000)) - 1 / 12), 0.007)
  # The two simulation populations of Chattopadhyay and De (2016), with the
  # values they report: lognormal(2.185, 0.562), xi^2 = 0.0532; and
  # gamma(shape 2.649, rate 0.84), xi^2 = 0.0468. Standard deviations of
  # the estimate at n = 50,000: 0.0012 and 0.0005.
  expect_lt(abs(.var_gini(rlnorm(50000, 2.185, 0.562)) - 0.0532), 0.005)
  expect_lt(abs(.var_gini(rgamma(50000, shape = 2.649, rate = 0.84)) -
                  0.0468), 0.0025)
})

test_that("every estimator returns its lower bound, silently, when the data have no variation", {
  k <- rep(3, 12)
  withr::local_seed(113)
  v <- rnorm(12)
  bound_12 <- 12^(-3)
  expect_silent(expect_identical(.var_cv(k), bound_12))
  expect_silent(expect_identical(.var_sm(k), bound_12))
  expect_silent(expect_identical(.var_gini(k), bound_12))
  expect_silent(expect_identical(.var_reg_coef(k, v), bound_12))
  expect_silent(expect_identical(.var_reg_coef(v, k), bound_12))
  for (m in c("pearson", "kendall", "spearman")) {
    expect_silent(expect_identical(.var_r(k, v, method = m), bound_12))
    expect_silent(expect_identical(.var_r(v, k, method = m), bound_12))
  }
})

test_that("with floor = FALSE the estimators say why there is no estimate", {
  # NaN where the effect size itself is undefined (a ratio 0 / 0) or the
  # variance of a group has not been estimated (a contrast with a constant
  # group), zero where the estimate is a true zero (no inequality; a slope
  # of zero on a constant response). Either way the check functions report
  # the variance as not yet usable.
  k <- rep(3, 12)
  withr::local_seed(113)
  v <- rnorm(12)
  expect_true(is.nan(.var_cv(k, floor = FALSE)))
  expect_true(is.nan(.var_sm(k, floor = FALSE)))
  expect_true(is.nan(.var_smd(k, k + 1, floor = FALSE)))
  expect_true(is.nan(.var_reg_coef(k, v, floor = FALSE)))
  for (m in c("pearson", "kendall", "spearman")) {
    expect_true(is.nan(.var_r(k, v, method = m, floor = FALSE)))
    expect_true(is.nan(.var_r(v, k, method = m, floor = FALSE)))
  }
  expect_identical(.var_gini(k, floor = FALSE), 0)
  expect_identical(.var_reg_coef(v, k, floor = FALSE), 0)
  expect_true(is.nan(.var_c(list(k, k + 2), c(1, -1), floor = FALSE)))
  expect_true(is.nan(.var_smd(k, v, floor = FALSE)))
  expect_true(is.nan(.var_smd(v, k, floor = FALSE)))
  # For data that vary, floor = FALSE and the default differ only by the
  # bound: here every estimate is above it, so the two agree exactly.
  # (Fixed data, so that this does not depend on a random draw.)
  x <- c(2.1, 3.4, 5.0, 1.7, 4.2, 6.3, 2.8, 3.9, 7.5, 3.1, 4.6, 2.4)
  v <- c(-0.3, 1.2, 0.4, -1.1, 0.8, 2.0, -0.6, 0.1, 1.5, -0.9, 0.3, 0.7)
  expect_identical(.var_cv(x, floor = FALSE), .var_cv(x))
  expect_identical(.var_sm(x, floor = FALSE), .var_sm(x))
  expect_identical(.var_gini(x, floor = FALSE), .var_gini(x))
  expect_identical(.var_smd(x, v, floor = FALSE), .var_smd(x, v))
  for (m in c("kendall", "spearman")) {
    expect_identical(.var_r(x, v, method = m, floor = FALSE),
                     .var_r(x, v, method = m))
  }
  expect_gt(min(.var_cv(x), .var_sm(x), .var_gini(x), .var_smd(x, v),
                .var_r(x, v, method = "kendall")), 12^(-3))
  # The internal estimator accepts a method by partial match and names the
  # argument when it cannot.
  expect_identical(.var_r(x, v, method = "ken"), .var_r(x, v, method = "kendall"))
  expect_error(.var_r(x, v, method = "tau"), "'method' must be one of")
})

test_that("the rank based estimators are exactly zero where they are zero on paper", {
  # Perfectly discordant pairs: every W_i of the Kendall estimator
  # (Equation 28) is 2 / n - 1, because only the pair itself is at or below
  # it on both margins and its two ranks sum to n + 1. The variance of
  # equal numbers is zero. Computed from ranks divided by n + 1 it came out
  # as about 1e-32, positive at 30 of the 37 sample sizes from 4 to 40 and
  # exactly zero at the other 7 (4, 7, 8, 15, 16, 31, 32).
  for (n in c(4, 5, 9, 10, 12, 14, 15, 17, 20, 33)) {
    x <- seq_len(n) / 7
    expect_identical(.var_r(x, -x, method = "kendall", floor = FALSE), 0)
    expect_identical(.var_r(rev(x), exp(x), method = "kendall", floor = FALSE),
                     0)
    expect_identical(.var_r(x, -x, method = "kendall"), n^(-3))
  }
  # Two dichotomous variables in perfect opposition, in equal halves: the
  # two halves have the same W_i and the same Z_i, so both estimates are
  # zero on paper (kendall_as_printed() and spearman_as_printed() above,
  # evaluated in exact rational arithmetic, give 0 for n = 4, 6, 8, 10, 12,
  # and 20).
  for (n in c(4, 6, 10, 12, 20)) {
    b <- rep(c(0, 1), each = n / 2)
    for (m in c("kendall", "spearman")) {
      expect_identical(.var_r(b, 1 - b, method = m, floor = FALSE), 0)
    }
  }
  # The Z_i of the Spearman estimator (Equation 52) are all 0.77 for
  # 1, 2, 3, 4 against 3, 4, 1, 2. For the first pair: 1 * 3 / 25, plus
  # (3 + 4 + 1 + 2) / (4 * 5), plus (1 + 2) / (4 * 5), which is
  # 0.12 + 0.5 + 0.15. For the second: 2 * 4 / 25 + (4 + 1 + 2) / 20 +
  # 2 / 20 = 0.32 + 0.35 + 0.1. The third and fourth mirror them.
  expect_identical(.var_r(c(1, 2, 3, 4), c(3, 4, 1, 2), method = "spearman",
                          floor = FALSE), 0)
  expect_identical(.var_r(c(1, 2, 3, 4), c(2, 1, 4, 3), method = "spearman",
                          floor = FALSE), 0)
  # Among all 46,224 pairings of 4 to 8 untied values (1, ..., n against
  # each permutation, the same printed forms in exact rational arithmetic)
  # these are the only zeros: the reversal for Kendall, the two above for
  # Spearman.
})

test_that("the rank based estimators equal their exact values, with and without ties", {
  # The values of kendall_as_printed() and spearman_as_printed() above
  # (Equations 27 to 29 and 51 to 53) evaluated in exact rational
  # arithmetic, with midranks for ties. Several are positive but below
  # n^-3, so they also show that a small estimate is not rounded to zero.
  cases <- list(
    list(x = c(1, 2, 3, 4, 5), y = c(2, 1, 5, 4, 3),
         kendall = 2 / 375, spearman = 107 / 450),
    list(x = c(1, 2, 3, 4, 5, 6, 7), y = c(3, 2, 1, 7, 6, 5, 4),
         kendall = 1 / 686, spearman = 243 / 1792),
    list(x = c(1, 2, 3, 4, 5, 6), y = c(3, 2, 1, 6, 5, 4),
         kendall = 24 / 245, spearman = 16 / 36015),
    list(x = c(1, 2, 3, 4, 5, 6, 7, 8), y = c(4, 3, 2, 1, 8, 7, 6, 5),
         kendall = 32 / 567, spearman = 2 / 5103),
    list(x = c(1, 2, 3, 4, 5, 6), y = c(2, 4, 1, 6, 3, 5),
         kendall = 1408 / 6615, spearman = 1040 / 7203),
    # With ties (average ranks).
    list(x = c(1, 1, 1, 1, 2), y = c(1, 1, 1, 2, 1),
         kendall = 1 / 750, spearman = 245 / 96),
    list(x = c(1, 1, 2, 3), y = c(2, 3, 1, 1),
         kendall = 12 / 25, spearman = 27 / 2500),
    list(x = c(1, 1, 2, 2, 2), y = c(2, 2, 1, 1, 2),
         kendall = 2209 / 1125, spearman = 4 / 1125),
    list(x = c(1, 1, 2, 2, 3, 3, 3), y = c(2, 1, 1, 3, 3, 2, 3),
         kendall = 19529 / 16464, spearman = 290805 / 702464)
  )
  for (case in cases) {
    for (m in c("kendall", "spearman")) {
      expect_equal(.var_r(case$x, case$y, method = m, floor = FALSE),
                   case[[m]], tolerance = 1e-13)
      # The value does not depend on the units or the order of the pairs.
      expect_equal(.var_r(rev(case$x) * 0.3, rev(case$y) + 11, method = m,
                          floor = FALSE), case[[m]], tolerance = 1e-13)
    }
  }
  # Perfectly concordant pairs. W_i = 2 i / n - 2 i / (n + 1) =
  # 2 i / (n (n + 1)), and the variance of 1, ..., n is n (n + 1) / 12, so
  # the Kendall estimate is 16 / (3 n (n + 1)): small, but not zero.
  for (n in c(5, 12, 40)) {
    x <- seq_len(n) / 7
    expect_equal(.var_r(x, 3 * x + 1, method = "kendall", floor = FALSE),
                 16 / (3 * n * (n + 1)), tolerance = 1e-13)
  }
  # The printed forms, as loops, on data with many ties.
  withr::local_seed(113)
  for (n in c(7, 30, 90)) {
    x <- sample(1:5, n, replace = TRUE)
    y <- pmin(5, pmax(1, x + sample(-2:2, n, replace = TRUE)))
    expect_equal(.var_r(x, y, method = "kendall", floor = FALSE),
                 kendall_as_printed(x, y), tolerance = 1e-12)
    expect_equal(.var_r(x, y, method = "spearman", floor = FALSE),
                 spearman_as_printed(x, y), tolerance = 1e-12)
  }
})

test_that("a Pearson or slope estimate that is zero on paper is exactly zero", {
  # For collinear pairs, y = a + b x, both estimators reduce to a multiple of
  # (n^2 - 3) S2^2 / (n (n - 1)) - (n - 1) S4, with S2 and S4 the sums of
  # squared and fourth power deviations of x: zero when the sample kurtosis
  # n S4 / S2^2 equals (n^2 - 3) / (n - 1)^2. For x below, the deviations
  # are 1.5 (three), 0.5, -1.5 (three), -0.5, so S2 = 14 and S4 = 30.5, and
  # 8 * 30.5 / 196 = 61 / 49 = (64 - 3) / 49. pearson_as_printed() and
  # slope_as_printed() above, evaluated in exact rational arithmetic, give
  # 0 for each of the five lines below. In floating point the printed forms
  # give rounding noise of either sign: 4.4e-16 for
  # y = x and -2.2e-16 for y = 0.3 x (Pearson), so that the first was used
  # through the bound, and the check met, where the second was not usable.
  x <- c(3, 3, 3, 2, 0, 0, 0, 1)
  expect_equal(8 * sum((x - 1.5)^4) / sum((x - 1.5)^2)^2, 61 / 49)
  lines <- list(function(x) x, function(x) 0.3 * x, function(x) 10 - x,
                function(x) 2 * x + 1, function(x) 1.8 * x + 32)
  for (f in lines) {
    y <- f(x)
    # Zero up to the noise of the raw power sums of the printed forms
    # (collinear pairs with another kurtosis give 0.02 and more; see below).
    expect_lt(abs(pearson_as_printed(x, y)), 1e-9)
    expect_lt(abs(slope_as_printed(x, y)), 1e-9 * max(1, var(y)))
    expect_identical(.var_r(x, y, method = "pearson", floor = FALSE), 0)
    expect_identical(.var_reg_coef(x, y, floor = FALSE), 0)
    expect_identical(.var_r(x, y, method = "pearson"), 8^(-3))
    expect_identical(.var_reg_coef(x, y), 8^(-3))
  }
  # Collinear pairs with any other kurtosis have a clearly negative
  # estimate, which is returned as it is: -901 / 41070 (Pearson) and
  # -901 / 20535 (slope) for y = x after the last value is changed to 2,
  # and -19 / 120 and -19 / 60 for 1, ..., 8 (the same printed forms in
  # exact rational arithmetic; the slope of 0.3 x is 0.09 times -19 / 60).
  x2 <- replace(x, 8, 2)
  expect_equal(.var_r(x2, x2, method = "pearson", floor = FALSE),
               -901 / 41070, tolerance = 1e-12)
  expect_equal(.var_reg_coef(x2, x2, floor = FALSE), -901 / 20535,
               tolerance = 1e-12)
  expect_equal(.var_r(1:8, 10 - 1:8, method = "pearson", floor = FALSE),
               -19 / 120, tolerance = 1e-12)
  expect_equal(.var_reg_coef(1:8, 0.3 * (1:8), floor = FALSE), -57 / 2000,
               tolerance = 1e-12)
})

test_that(".aipe_zero_if_cancelled returns zero for cancellation noise only", {
  # A sum that is no more than 1e-10 of the size of its terms is noise.
  expect_identical(.aipe_zero_if_cancelled(4.4e-16, c(1.5, -2, 0.5)), 0)
  expect_identical(.aipe_zero_if_cancelled(-2.2e-16, c(1.5, -2, 0.5)), 0)
  expect_identical(.aipe_zero_if_cancelled(0, c(1.5, -2, 0.5)), 0)
  expect_identical(.aipe_zero_if_cancelled(3.9e-10, c(1.5, -2, 0.5)), 0)
  # Anything larger is returned untouched, whatever its sign or size: a
  # small estimate is not a cancelled one unless its terms are large.
  expect_identical(.aipe_zero_if_cancelled(4.1e-10, c(1.5, -2, 0.5)), 4.1e-10)
  expect_identical(.aipe_zero_if_cancelled(-0.3, c(1.5, -2, 0.2)), -0.3)
  expect_identical(.aipe_zero_if_cancelled(2.5e-30, c(1e-30, 1.5e-30)), 2.5e-30)
  expect_identical(.aipe_zero_if_cancelled(1e-300, 1e-300), 1e-300)
  # Values that are not finite pass through for the caller to report.
  expect_true(is.nan(.aipe_zero_if_cancelled(NaN, c(NaN, 1))))
  expect_identical(.aipe_zero_if_cancelled(Inf, c(Inf, 1)), Inf)
  expect_identical(.aipe_zero_if_cancelled(-Inf, c(-Inf, 1)), -Inf)
})

test_that("the sign of rounding noise never decides whether an estimate is usable", {
  # The moment based estimators change sign along a path from a sample
  # with a negative estimate to one with a positive estimate. On paper the
  # estimate is zero at one point of the path; the values computed next to
  # it are noise of either sign, and a positive one would be used through
  # the bound n^-3 while a negative one is not usable. Each estimator must
  # return exactly zero there (sign_change() then stops at that point;
  # without the cancellation step it ends with two adjacent doubles and
  # noise at both). Fixed data, so nothing depends on a random draw.
  crossing <- function(estimate, as_printed, sample_at) {
    f <- function(t) estimate(sample_at(t))
    expect_lt(f(0), 0)
    expect_gt(f(1), 0)
    at <- sign_change(f)
    expect_identical(at[1], at[2])               # an exact zero was found
    expect_identical(f(at[1]), 0)
    # It is the crossing of the printed formula: there the printed value
    # is within 1e-10 of the size of its terms (which are of the order of
    # 1 to 20 here) of zero, and a little to either side the package
    # returns the printed values, of opposite signs: nothing is set to
    # zero except at the crossing itself.
    expect_lt(abs(as_printed(sample_at(at[1]))), 1e-8)
    for (step in c(1e-6, 1e-3)) {
      below <- sample_at(at[1] - step)
      above <- sample_at(at[1] + step)
      expect_lt(f(at[1] - step), 0)
      expect_gt(f(at[1] + step), 0)
      expect_equal(f(at[1] - step), as_printed(below), tolerance = 1e-6)
      expect_equal(f(at[1] + step), as_printed(above), tolerance = 1e-6)
    }
    at[1]
  }
  # One sample: five scores of 1 and one of 10 (negative estimates for the
  # coefficient of variation and the standardized mean, -0.41 and -0.088)
  # moved toward 1, ..., 6 (0.108 and 1.325).
  one <- function(t) (1 - t) * c(1, 1, 1, 1, 1, 10) + t * c(1, 2, 3, 4, 5, 6)
  t_cv <- crossing(function(x) .var_cv(x, floor = FALSE), cv_var_as_printed,
                   one)
  t_sm <- crossing(function(x) .var_sm(x, floor = FALSE), sm_var_as_printed,
                   one)
  # The two estimates differ by the positive factor (s / mean)^4, so they
  # cross at the same sample.
  expect_equal(t_cv, t_sm, tolerance = 1e-8)
  # And the checks say "not yet usable" there, not "met".
  expect_true(attr(seq_aipe_cv(one(t_cv), width = 2), "variance_unusable"))
  expect_false(attr(seq_aipe_cv(one(t_cv), width = 2), "is_satisfied"))
  expect_true(attr(seq_aipe_sm(one(t_sm), width = 2), "variance_unusable"))
  expect_false(attr(seq_aipe_sm(one(t_sm), width = 2), "is_satisfied"))

  # Two groups: 20, 20, 20, 30 against 10, 10, 10, 0 (skewed in opposite
  # directions, with means far apart) moved toward 21, 24, 26, 29 against
  # 1, 4, 6, 9 (5.14). The first pair by hand: both variances are 25, the
  # mean difference is 15, the third moment estimates are 250 and -250
  # (4 * 375 / 6), and the fourth moment estimates are both
  # (11 * 3281.25 - 3.75 * 75^2) / 6 = 2500, so that
  #   t_1 = (1/2) 15 * 500 / 625 = 6,
  #   t_2 t_3 = (225 / (4 * 15625)) (5000 / 4 - 625 / 2) = 3.375,
  # and the estimate is 2 - 6 + 3.375 = -0.625.
  expect_equal(.var_smd(c(20, 20, 20, 30), c(10, 10, 10, 0), floor = FALSE),
               -0.625)
  two <- function(t) {
    list((1 - t) * c(20, 20, 20, 30) + t * c(21, 24, 26, 29),
         (1 - t) * c(10, 10, 10, 0) + t * c(1, 4, 6, 9))
  }
  smd_as_printed <- function(g) {
    x <- g[[1]]; y <- g[[2]]
    sp2 <- (var(x) + var(y)) / 2
    d <- mean(x) - mean(y)
    2 - 0.5 * d * (mu3_as_printed(x) - mu3_as_printed(y)) / sp2^2 +
      d^2 / (4 * sp2^3) *
      ((mu4_as_printed(x) + mu4_as_printed(y)) / 4 - sp2^2 / 2)
  }
  t_smd <- crossing(function(g) .var_smd(g[[1]], g[[2]], floor = FALSE),
                    smd_as_printed, two)
  g <- two(t_smd)
  out <- seq_aipe_smd(g[[1]], g[[2]], width = 3)
  expect_true(attr(out, "variance_unusable"))
  expect_false(attr(out, "is_satisfied"))
})

test_that("the variance estimators are location invariant where the effect size is", {
  # Dyadic data, so that the shift itself is exact. With raw power sums the
  # Pearson estimate went from 0.58 to 52 at a shift of 1e4.
  withr::local_seed(113)
  x <- round(rnorm(200) * 1024) / 1024
  y <- round((0.5 * x + rnorm(200)) * 1024) / 1024
  for (shift in c(1e3, 1e6, 1.76e9)) {
    for (m in c("pearson", "kendall", "spearman")) {
      expect_equal(.var_r(x + shift, y - shift, method = m),
                   .var_r(x, y, method = m), tolerance = 1e-8)
    }
    expect_equal(.var_reg_coef(x + shift, y + shift), .var_reg_coef(x, y),
                 tolerance = 1e-8)
  }
  # The contrast variance does not depend on the group means either.
  g <- list(x, y, x + y)
  expect_equal(.var_c(lapply(g, function(v) v + 1e6), c(1, -0.5, -0.5)),
               .var_c(g, c(1, -0.5, -0.5)), tolerance = 1e-9)
})

test_that(".var_cv, .var_sm, and .var_gini remove missing values when na.rm = TRUE", {
  # No check function passes na.rm (they refuse missing values first), so
  # the argument is tested here: the estimate is that of the complete
  # values.
  withr::local_seed(113)
  x <- stats::rgamma(40, shape = 4)
  with_na <- c(x[1:20], NA, x[21:40], NA)
  for (f in list(.var_cv, .var_sm, .var_gini)) {
    expect_identical(f(with_na, na.rm = TRUE), f(x))
    expect_identical(f(with_na, na.rm = TRUE, floor = FALSE),
                     f(x, floor = FALSE))
  }
})
