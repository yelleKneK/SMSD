# Tests for seq_aipe_cv() and seq_aipe_cv_pilot(): the sequential AIPE
# procedure for the coefficient of variation of Kelley, Darku, and
# Chattopadhyay (2018).

row_value <- function(out, term) out$value[out$term == term]
cv_stat <- function(d) c(sd(d) / mean(d), .var_cv(d, floor = FALSE))
take_first <- function(d, n) d[1:n]

test_that("seq_aipe_cv_pilot is max{m_0, ceiling(2 z / width)} (2018, Equation 21)", {
  # Width 0.04, the target of the article's Tables 3 and 4:
  # 2 * 1.959964 / 0.04 = 98.0 (just under 98), and 2 * 1.644854 / 0.04 =
  # 82.2, so 83.
  expect_identical(seq_aipe_cv_pilot(width = 0.04), 98L)
  expect_identical(seq_aipe_cv_pilot(width = 0.04, conf_level = 0.90), 83L)
  expect_identical(seq_aipe_cv_pilot(width = 0.30), 14L)
  expect_identical(seq_aipe_cv_pilot(width = 0.30, m_0 = 25), 25L)
})

test_that("seq_aipe_cv_pilot checks its arguments", {
  expect_pilot_checks(seq_aipe_cv_pilot)
})

test_that("seq_aipe_cv returns the documented rows and attributes", {
  withr::local_seed(113)
  x <- rgamma(60, shape = 4, scale = 2)
  out <- seq_aipe_cv(x, width = 0.30)
  expect_s3_class(out, c("smsd_check", "data.frame"), exact = TRUE)
  expect_identical(out$term, c("lower_limit", "cv", "upper_limit",
                               "xi_squared", "criterion_n", "current_n"))
  expect_true(all(is.finite(out$value)))
  expect_identical(attr(out, "method"), "coefficient of variation")
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_identical(attr(out, "notes"), character(0))
  expect_identical(attr(out, "conf_level"), 0.95)
  expect_identical(attr(out, "width"), 0.30)
  expect_equal(attr(out, "current_n"), 60)
  expect_silent(seq_aipe_cv(x, width = 0.30))
})

test_that("seq_aipe_cv computes the estimate, the interval (2018, Equation 33), and the criterion (Equation 22)", {
  withr::local_seed(113)
  x <- rgamma(60, shape = 4, scale = 2)
  for (cl in c(0.90, 0.95)) {
    out <- seq_aipe_cv(x, width = 0.30, conf_level = cl)
    k <- sd(x) / mean(x)
    xi2 <- .var_cv(x)
    z <- qnorm(1 - (1 - cl) / 2)
    expect_equal(row_value(out, "cv"), k)
    expect_equal(row_value(out, "xi_squared"), xi2)
    expect_equal(row_value(out, "lower_limit"), k - z * sqrt(xi2 / 60))
    expect_equal(row_value(out, "upper_limit"), k + z * sqrt(xi2 / 60))
    expect_equal(row_value(out, "criterion_n"),
                 ceiling((2 * z / 0.30)^2 * (xi2 + 1 / 60)))
    expect_identical(attr(out, "is_satisfied"),
                     60 >= row_value(out, "criterion_n"))
  }
})

test_that("seq_aipe_cv meets the rule exactly when n reaches the criterion", {
  withr::local_seed(113)
  x <- rgamma(2000, shape = 4, scale = 2)
  wide <- seq_aipe_cv(x, width = 1.00)
  expect_true(attr(wide, "is_satisfied"))
  expect_false(attr(seq_aipe_cv(x[1:30], width = 0.02), "is_satisfied"))
})

test_that("the stopping rule cannot be met below the pilot size, because of its 1 / n term", {
  # n >= (2z/w)^2 (xi2 + 1/n) >= (2z/w)^2 / n forces n >= 2z/w, which is why
  # the pilot is at least ceiling(2z/w). The coefficient of variation shows
  # it, because its variance is small: for normal data with kappa = 0.1,
  # xi^2 = kappa^2 (kappa^2 + 1/2) = 0.0051. One observation short of the
  # pilot of 40 for width 0.1, the rule without the 1 / n term would be met
  # many times over ((2z/w)^2 xi^2 is about 8), and the rule itself is not.
  w <- 0.1
  z <- qnorm(0.975)
  m <- seq_aipe_cv_pilot(w)
  expect_identical(m, 40L)
  withr::local_seed(113)
  for (i in 1:20) {
    out <- seq_aipe_cv(rnorm(m - 1, mean = 10, sd = 1), width = w)
    xi2 <- row_value(out, "xi_squared")
    expect_false(attr(out, "variance_unusable"))
    expect_lte(ceiling((2 * z / w)^2 * xi2), m - 1)  # met without 1 / n
    expect_false(attr(out, "is_satisfied"))          # not met with it
    expect_gte(row_value(out, "criterion_n"), m)
  }
})

test_that("seq_aipe_cv does not depend on the unit of measurement", {
  withr::local_seed(113)
  x <- rgamma(80, shape = 4, scale = 2)
  expect_equal(seq_aipe_cv(1000 * x, width = 0.3)$value,
               seq_aipe_cv(x, width = 0.3)$value, tolerance = 1e-12)
  # Integer data (counts) give the same answer as the same values as doubles.
  counts <- as.integer(round(10 * x))
  expect_equal(seq_aipe_cv(counts, width = 0.3)$value,
               seq_aipe_cv(as.numeric(counts), width = 0.3)$value)
})

test_that("seq_aipe_cv allows negative observations when the mean is positive", {
  # The article's Table 4 uses normal populations, which have no lower
  # bound; only the mean must be positive.
  x <- c(-1, 4, 6, 3, 8, 5, 7, 2, 9, 4)
  out <- seq_aipe_cv(x, width = 0.5)
  expect_equal(row_value(out, "cv"), sd(x) / mean(x))
  expect_true(is.finite(row_value(out, "lower_limit")))
})

test_that("seq_aipe_cv checks its arguments", {
  withr::local_seed(113)
  good <- rgamma(20, shape = 4)
  expect_error(seq_aipe_cv(width = 0.3), "'x' must be supplied")
  expect_data_checks(function(v) seq_aipe_cv(v, width = 0.3), "x", good)
  expect_scalar_checks(function(...) seq_aipe_cv(good, ...))
})

test_that("the exported function agrees with the test search on whole paths", {
  withr::local_seed(113)
  check <- function(d) seq_aipe_cv(d, width = 0.15)
  pilot <- seq_aipe_cv_pilot(0.15)
  for (i in 1:3) {
    d <- rlnorm(1500, 0, 0.5)
    run <- aipe_run(d, take_first, cv_stat, check, width = 0.15, pilot = pilot)
    expect_equal(aipe_walk(d, take_first, check, pilot = pilot), run[["n"]])
  }
})

test_that("operating characteristics: Kelley, Darku, and Chattopadhyay (2018), Table 4", {
  # Row N(10, 4) (mean 10, variance 4, kappa = 0.2), width 0.04, 95%: mean
  # final sample size 241.2 (theoretical 208), coverage 0.9422 (se 0.0033),
  # mean width 0.0363, from 5,000 replications. In 5,000 replications of
  # this package: 242.9 (standard deviation 32.1), 0.936 (se 0.004), 0.0363;
  # in 20,000 more (200 other seeds of this test): 241.8, 0.935, 0.0363.
  # With the 100 replications here the standard error of the mean is
  # 32.1 / sqrt(100) = 3.21, so its window is 4 * 3.21 + 1.7 = 14.5. The
  # standard error of the coverage is sqrt(0.935 * 0.065 / 100) = 0.0247 (and
  # 0.0249 over the 200 seeds), so its window is 4 * 0.0249 + 0.007 = 0.107.
  # The mean width had a standard deviation of 0.0001 over the 200 seeds,
  # so 0.0005 is five of them.
  withr::local_seed(113)
  w <- 0.04
  check <- function(d) seq_aipe_cv(d, width = w)
  runs <- replicate(100, {
    aipe_run(rnorm(700, 10, 2), take_first, cv_stat, check, width = w,
             pilot = seq_aipe_cv_pilot(w))
  })
  s <- aipe_summary(runs, truth = 0.2)
  expect_lt(abs(s[["mean_n"]] - 241.2), 14.5)
  expect_lt(abs(s[["coverage"]] - 0.9422), 0.107)
  expect_lt(abs(s[["mean_width"]] - 0.0363), 0.0005)
  expect_lte(s[["max_width"]], w)              # the width guarantee is exact
})

test_that("operating characteristics at 90%: Kelley, Darku, and Chattopadhyay (2018), Table 3", {
  # Row N(10, 4), width 0.04, 90%: mean final sample size 178.8
  # (theoretical 147), coverage 0.8844 (se 0.0049), mean width 0.0351. In
  # 5,000 replications of this package: 179.1 (standard deviation 25.0),
  # 0.882 (se 0.005), 0.0351; in 12,000 more (200 other seeds of this
  # test): 179.1, 0.879, 0.0351. The mean of the 60 replications here had a
  # standard deviation of 3.28 over the 200 seeds, so its window is
  # 4 * 3.28 + 0.3 = 13.4. The standard error of the coverage is
  # sqrt(0.879 * 0.121 / 60) = 0.0421, so its window is
  # 4 * 0.0421 + 0.006 = 0.175. The mean width had a standard deviation of
  # 0.00019 over the 200 seeds, so its window is 4 * 0.00019 = 0.0008.
  withr::local_seed(113)
  w <- 0.04
  check <- function(d) seq_aipe_cv(d, width = w, conf_level = 0.90)
  runs <- replicate(60, {
    aipe_run(rnorm(600, 10, 2), take_first, cv_stat, check, width = w,
             pilot = seq_aipe_cv_pilot(w, conf_level = 0.90),
             conf_level = 0.90)
  })
  s <- aipe_summary(runs, truth = 0.2)
  expect_lt(abs(s[["mean_n"]] - 178.8), 13.4)
  expect_lt(abs(s[["coverage"]] - 0.8844), 0.175)
  expect_lt(abs(s[["mean_width"]] - 0.0351), 0.0008)
  expect_lte(s[["max_width"]], w)
})
