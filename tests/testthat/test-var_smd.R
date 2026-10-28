# Tests for .var_smd(), the estimator of the asymptotic variance of the
# standardized mean difference, which implements the corrected Equation 26
# of Kelley, Darku, and Chattopadhyay (2018). The published Equation 26
# omits a factor of 1/2 in the term carrying the third central moments; the
# correction is at
# https://kenkelley.org/wp-content/uploads/articles/Kelley_Darku_Chattopadhyay_Equation_26.pdf
# and ships as inst/eq26_correction.pdf.
#
# Every value here is independent of the three term expression the function
# computes: the delta method quadratic form D' Sigma D from the article's
# Theorem 3 matrix, the population value in a case where the corrected and
# published forms differ by a factor of about 2.5, and the sampling variance
# of d itself.

# D' Sigma D with the sample moments, from the matrix of Theorem 3 (2018,
# p. 242): Sigma = [2 s^2, (m13 - m23)/2; (m13 - m23)/2, (m14 + m24 - 2 s^4)/4]
# and D = (1/s, -(xbar1 - xbar2)/(2 s^3)).
quadratic_form <- function(x, y) {
  s2 <- (var(x) + var(y)) / 2
  s <- sqrt(s2)
  D <- c(1 / s, -(mean(x) - mean(y)) / (2 * s^3))
  off <- (.mu3(x) - .mu3(y)) / 2
  Sigma <- matrix(c(2 * s2, off, off,
                    (.mu4(x) + .mu4(y) - 2 * s2^2) / 4), 2, 2)
  drop(t(D) %*% Sigma %*% D)
}

# Equation 26 exactly as published (no factor 1/2 in the second term), with
# the same sample moments.
published_form <- function(x, y) {
  s2 <- (var(x) + var(y)) / 2
  diff <- mean(x) - mean(y)
  2 - diff * (.mu3(x) - .mu3(y)) / s2^2 +
    diff^2 / (4 * s2^3) * ((.mu4(x) + .mu4(y)) / 4 - s2^2 / 2)
}

test_that(".var_smd equals the delta method quadratic form D' Sigma D", {
  withr::local_seed(113)
  # Groups skewed in opposite directions, so the third moment term is large.
  x <- aipe_skewed(200) + 0.5
  y <- -aipe_skewed(200)
  expect_equal(.var_smd(x, y), quadratic_form(x, y), tolerance = 1e-12)
  # The published form is a different number. The population values differ
  # by 0.5 sqrt(8/3) = 0.816; with 2,000 per group the difference of the
  # two sample forms has mean 0.814 and standard deviation 0.060 (1,000
  # samples), so 0.4 is almost seven standard deviations below it.
  x <- aipe_skewed(2000) + 0.5
  y <- -aipe_skewed(2000)
  expect_equal(.var_smd(x, y), quadratic_form(x, y), tolerance = 1e-12)
  expect_gt(.var_smd(x, y) - published_form(x, y), 0.4)

  # Symmetric groups with equal third moments: the two forms then agree, so
  # this case alone could not detect the missing factor.
  x <- rnorm(60, mean = 0.5)
  y <- rnorm(60)
  expect_equal(.var_smd(x, y), quadratic_form(x, y), tolerance = 1e-12)
})

test_that(".var_smd converges to the corrected population variance, not the published one", {
  # X = S + 0.5 and Y = -S with S a standardized chi square (3 df): delta =
  # 0.5, sigma = 1, third central moments +/- sqrt(8/3), fourth central
  # moments 7. Corrected Equation 26:
  #   2 - (1/2)(0.5)(2 sqrt(8/3)) + (0.25/4)(14/4 - 1/2) = 1.3710.
  # As published (no 1/2): 2 - 0.5 * 2 sqrt(8/3) + 0.1875 = 0.5545.
  corrected <- 2 - 0.5 * 0.5 * 2 * sqrt(8 / 3) + 0.25 / 4 * (14 / 4 - 0.5)
  published <- 2 - 0.5 * 2 * sqrt(8 / 3) + 0.25 / 4 * (14 / 4 - 0.5)
  expect_equal(corrected, 1.3710, tolerance = 1e-4)
  expect_equal(published, 0.5545, tolerance = 1e-4)
  withr::local_seed(113)
  x <- aipe_skewed(20000) + 0.5
  y <- -aipe_skewed(20000)
  # The standard deviation of the estimate at n = 20,000 is about 0.010
  # (200 replications), so 0.05 is more than four standard errors.
  expect_lt(abs(.var_smd(x, y) - corrected), 0.05)
})

test_that("the corrected variance matches the sampling variance of d", {
  # The same populations: n Var(d_n) over 4,000 samples of n = 200 per
  # group. The Monte Carlo standard error of n Var(d) is about
  # 1.37 * sqrt(2 / 3999) = 0.031, so 0.13 is more than four standard
  # errors; the published form (0.55) is 26 standard errors away.
  withr::local_seed(113)
  n <- 200
  d <- replicate(4000, {
    x <- aipe_skewed(n) + 0.5
    y <- -aipe_skewed(n)
    (mean(x) - mean(y)) / sqrt((var(x) + var(y)) / 2)
  })
  expect_lt(abs(n * var(d) - 1.3710), 0.13)
})

test_that(".var_smd under normality is near 2 + delta^2 / 4", {
  # For normal data the third moments vanish and the fourth are 3 sigma^4:
  # xi^2 = 2 + delta^2 / 4 (the large sample variance of Cohen's d with
  # equal groups, in per group n).
  withr::local_seed(113)
  x <- rnorm(20000, mean = 0.8)
  y <- rnorm(20000)
  # The standard deviation of the estimate at n = 20,000 is about 0.011
  # (100 replications), so 0.05 is more than four standard errors.
  expect_lt(abs(.var_smd(x, y) - (2 + 0.8^2 / 4)), 0.05)
})

test_that(".var_smd is location and scale invariant", {
  withr::local_seed(113)
  x <- round(rnorm(50, 0.5) * 1024) / 1024
  y <- round(rnorm(50) * 1024) / 1024
  v <- .var_smd(x, y)
  expect_equal(.var_smd(x + 1e6, y + 1e6), v, tolerance = 1e-9)
  expect_equal(.var_smd(8 * x, 8 * y), v, tolerance = 1e-12)
  # The raw power sum form returned 71.4 instead of about 2 for groups near
  # 100,000 with unit variance.
  withr::local_seed(113)
  a <- rnorm(200, 1e5, 1)
  b <- rnorm(200, 1e5 - 0.5, 1)
  expect_equal(.var_smd(a, b), .var_smd(a - 1e5, b - 1e5), tolerance = 1e-6)
})

test_that(".var_smd has no raw value, and returns its lower bound, when there is no variation", {
  expect_identical(.var_smd(rep(1, 10), rep(2, 10)), 10^(-3))
  expect_true(is.nan(.var_smd(rep(1, 10), rep(2, 10), floor = FALSE)))
  # One constant group: the formula has a value (the pooled variance is
  # half the variance of the other group), but the variance of the constant
  # group has not been estimated, so there is no raw value either, whichever
  # group it is, and the bound is returned.
  withr::local_seed(113)
  v <- rnorm(30)
  expect_gt(stats::var(v), 0)
  expect_true(is.nan(.var_smd(rep(1, 30), v, floor = FALSE)))
  expect_true(is.nan(.var_smd(v, rep(1, 30), floor = FALSE)))
  expect_identical(.var_smd(rep(1, 30), v), 30^(-3))
  # A group that is constant only up to rounding counts as constant.
  noise <- 1 + c(0, 1e-15, -1e-15)[rep(1:3, 10)]
  expect_gt(stats::var(noise), 0)
  expect_true(is.nan(.var_smd(noise, v, floor = FALSE)))
})

test_that(".var_smd removes missing values when na.rm = TRUE", {
  withr::local_seed(113)
  x <- c(rnorm(40), NA, NA)
  y <- c(rnorm(40), NA, NA)
  expect_equal(.var_smd(x, y, na.rm = TRUE),
               .var_smd(x[!is.na(x)], y[!is.na(y)]))
})
