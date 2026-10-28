# Tests for the central moment estimators (.mu3, .mu4, .mu13, .mu22, .mu31).
# Values are pinned against hand computation, against the formulas as
# printed in the source articles (in raw power sums, on small integer data
# where those sums are exact), against the exact unbiasedness the articles
# claim (by full enumeration of a small discrete population), and for
# location invariance to working precision.

test_that(".mu3 and .mu4 equal hand computed values", {
  # x = 1, 2, 3, 4, 10: mean 4, deviations -3, -2, -1, 0, 6, so
  # sum d^2 = 50, sum d^3 = 180, sum d^4 = 1394, and with n = 5
  #   mu3 = 5 * 180 / (4 * 3) = 75,
  #   mu4 = ((25 - 10 + 3) * 1394 - 3 * 7 / 5 * 50^2) / (4 * 3 * 2) = 608.
  x <- c(1, 2, 3, 4, 10)
  expect_equal(.mu3(x), 75)
  expect_equal(.mu4(x), 608)
})

test_that(".mu4 equals the raw power sum form printed in the articles", {
  withr::local_seed(113)
  for (n in c(4, 7, 25)) {
    x <- sample(-5:9, n, replace = TRUE)       # integers: raw sums are exact
    expect_equal(.mu4(x), mu4_as_printed(x), tolerance = 1e-12)
  }
})

test_that(".mu3 and .mu4 are exactly unbiased (full enumeration, n = 4)", {
  # Every sample of size 4 from a five point population, weighted by its
  # probability: the expectation of each estimator must equal the population
  # central moment exactly.
  support <- c(0, 1, 3, 7, 8)
  prob <- c(0.30, 0.25, 0.20, 0.15, 0.10)
  mu <- sum(prob * support)
  central <- function(r) sum(prob * (support - mu)^r)
  idx <- as.matrix(expand.grid(rep(list(1:5), 4)))
  w <- apply(idx, 1, function(i) prod(prob[i]))
  e3 <- sum(w * apply(idx, 1, function(i) .mu3(support[i])))
  e4 <- sum(w * apply(idx, 1, function(i) .mu4(support[i])))
  expect_equal(e3, central(3), tolerance = 1e-10)   # 19.686
  expect_equal(e4, central(4), tolerance = 1e-10)   # 148.2197
})

test_that(".mu22, .mu13, and .mu31 equal hand computed values", {
  # x = 1, 2, 3, 4, 10 and y = 2, 1, 4, 3, 5: deviations dx = -3, -2, -1, 0,
  # 6 and dy = -1, -2, 1, 0, 2, so S20 = 50, S02 = 10, S11 = 18, S22 = 170,
  # S13 = 66, S31 = 474; k20 = 12.5, k02 = 2.5, k11 = 4.5; and with n = 5
  #   k22 = 5/24 (6 * 170 - 0.8 * 500 - 1.6 * 324) = 21.1666...
  #   k13 = 5/24 (6 * 66 - 2.4 * 18 * 10)          = -7.5
  #   k31 = 5/24 (6 * 474 - 2.4 * 18 * 50)         = 142.5
  # giving mu22 = k22 + k20 k02 + 2 k11^2, mu13 = k13 + 3 k02 k11,
  # mu31 = k31 + 3 k20 k11 (2019, Equation 73).
  x <- c(1, 2, 3, 4, 10)
  y <- c(2, 1, 4, 3, 5)
  expect_equal(.mu22(x, y), 5 / 24 * 101.6 + 12.5 * 2.5 + 2 * 4.5^2)
  expect_equal(.mu13(x, y), -7.5 + 3 * 2.5 * 4.5)
  expect_equal(.mu31(x, y), 142.5 + 3 * 12.5 * 4.5)
})

test_that(".mu22, .mu13, and .mu31 equal the raw power sum forms printed in the article", {
  withr::local_seed(113)
  for (n in c(4, 6, 20)) {
    x <- sample(-4:8, n, replace = TRUE)
    y <- sample(-6:6, n, replace = TRUE)
    kx <- k_as_printed(x, y)
    ky <- k_as_printed(y, x)                   # k31(x, y) = k13(y, x)
    expect_equal(.mu22(x, y),
                 kx[["k22"]] + var(x) * var(y) + 2 * cov(x, y)^2,
                 tolerance = 1e-10)
    expect_equal(.mu13(x, y), kx[["k13"]] + 3 * var(y) * cov(x, y),
                 tolerance = 1e-10)
    expect_equal(.mu31(x, y), ky[["k13"]] + 3 * var(x) * cov(x, y),
                 tolerance = 1e-10)
  }
})

test_that("the k-statistics inside .mu22, .mu13, and .mu31 are exactly unbiased", {
  # Full enumeration of samples of size 4 from a four point bivariate
  # population. k22 = mu22_hat - k20 k02 - 2 k11^2 must have expectation
  # kappa_22 = mu_22 - mu_20 mu_02 - 2 mu_11^2, and likewise for k13 and k31
  # (kappa_13 = mu_13 - 3 mu_02 mu_11).
  bx <- c(0, 1, 2, 5)
  by <- c(1, 0, 4, 3)
  prob <- c(0.4, 0.3, 0.2, 0.1)
  ex <- sum(prob * bx)
  ey <- sum(prob * by)
  m <- function(i, j) sum(prob * (bx - ex)^i * (by - ey)^j)
  idx <- as.matrix(expand.grid(rep(list(1:4), 4)))
  w <- apply(idx, 1, function(i) prod(prob[i]))
  expectation <- function(f) {
    sum(w * apply(idx, 1, function(i) f(bx[i], by[i])))
  }
  expect_equal(
    expectation(function(x, y) .mu22(x, y) - var(x) * var(y) - 2 * cov(x, y)^2),
    m(2, 2) - m(2, 0) * m(0, 2) - 2 * m(1, 1)^2, tolerance = 1e-10)   # -4.02
  expect_equal(
    expectation(function(x, y) .mu13(x, y) - 3 * var(y) * cov(x, y)),
    m(1, 3) - 3 * m(0, 2) * m(1, 1), tolerance = 1e-10)               # -4.73
  expect_equal(
    expectation(function(x, y) .mu31(x, y) - 3 * var(x) * cov(x, y)),
    m(3, 1) - 3 * m(2, 0) * m(1, 1), tolerance = 1e-10)               # 0.412
})

test_that("the bivariate estimators converge to the joint central moments", {
  # Direct definition: mu_ij = E[(X - mu_X)^i (Y - mu_Y)^j], here estimated
  # by the plain average of products of deviations. At n = 50,000 the two
  # differ only by terms of order 1/n.
  withr::local_seed(113)
  n <- 50000
  x <- rexp(n)
  y <- 0.5 * x + rchisq(n, 2)
  plain <- function(i, j) mean((x - mean(x))^i * (y - mean(y))^j)
  expect_equal(.mu22(x, y), plain(2, 2), tolerance = 1e-3)
  expect_equal(.mu13(x, y), plain(1, 3), tolerance = 1e-3)
  expect_equal(.mu31(x, y), plain(3, 1), tolerance = 1e-3)
})

test_that(".mu22(x, x) and .mu13(x, x) reduce to k4 + 3 k2^2", {
  # With y = x both estimate the fourth central moment as k4 + 3 k2^2, where
  # k4 = n^2 ((n + 1) m4 - 3 (n - 1) m2^2) / ((n - 1)(n - 2)(n - 3)) is the
  # univariate fourth k-statistic (m_r the average r-th power of the
  # deviations).
  withr::local_seed(113)
  x <- rgamma(15, shape = 2)
  n <- length(x)
  d <- x - mean(x)
  k4 <- n^2 * ((n + 1) * mean(d^4) - 3 * (n - 1) * mean(d^2)^2) /
    ((n - 1) * (n - 2) * (n - 3))
  expect_equal(.mu22(x, x), k4 + 3 * var(x)^2)
  expect_equal(.mu13(x, x), k4 + 3 * var(x)^2)
  expect_equal(.mu31(x, x), k4 + 3 * var(x)^2)
})

test_that("the moment estimators are location invariant to working precision", {
  # Dyadic data (multiples of 1/1024), so that adding the shift is exact in
  # double precision and any change in the result is the estimator's own
  # rounding. The raw power sum forms failed this at a shift of 1e4.
  withr::local_seed(113)
  x <- round(rnorm(200) * 1024) / 1024
  y <- round((0.5 * x + rnorm(200)) * 1024) / 1024
  # 1.76e9 is about the number of seconds from 1970 to 2025 (a POSIXct
  # timestamp used as a variable).
  for (shift in c(1e3, 1e6, 1.76e9)) {
    expect_equal(.mu3(x + shift), .mu3(x), tolerance = 1e-9)
    expect_equal(.mu4(x + shift), .mu4(x), tolerance = 1e-9)
    expect_equal(.mu22(x + shift, y + shift), .mu22(x, y), tolerance = 1e-9)
    expect_equal(.mu13(x + shift, y - shift), .mu13(x, y), tolerance = 1e-9)
    expect_equal(.mu31(x - shift, y + shift), .mu31(x, y), tolerance = 1e-9)
  }
})

test_that(".mu3 and .mu4 check their input", {
  expect_error(.mu3(c(1, 2)), "at least 3")
  expect_error(.mu3("a"), "must be a numeric vector")
  expect_error(.mu4(c(1, 2, 3)), "at least 4")
  expect_error(.mu4("a"), "must be a numeric vector")
  expect_error(.mu3(1:5, na.rm = NA), "'na.rm' must be TRUE or FALSE")
})

test_that("missing values are removed, pairwise for the bivariate estimators, when na.rm = TRUE", {
  x <- c(1, 2, 3, NA, 4, 6)
  expect_equal(.mu3(x, na.rm = TRUE), .mu3(c(1, 2, 3, 4, 6)))
  expect_equal(.mu4(x, na.rm = TRUE), .mu4(c(1, 2, 3, 4, 6)))
  x <- c(1, 2, NA, 4, 5, 6, 7)
  y <- c(2, NA, 3, 5, 6, 7, 9)
  ok <- !is.na(x) & !is.na(y)
  expect_equal(.mu22(x, y, na.rm = TRUE), .mu22(x[ok], y[ok]))
  expect_equal(.mu13(x, y, na.rm = TRUE), .mu13(x[ok], y[ok]))
  expect_equal(.mu31(x, y, na.rm = TRUE), .mu31(x[ok], y[ok]))
  expect_equal(.var_r(x, y, na.rm = TRUE), .var_r(x[ok], y[ok]))
  expect_equal(.var_reg_coef(x, y, na.rm = TRUE), .var_reg_coef(x[ok], y[ok]))
})
