# Smoke tests for low level helpers: gini_index, gini_mean_diff, and
# slr_slope.

test_that("gini_index equals hand computed values", {
  # 1, 2, 4: the three absolute differences are 1, 3, 2 (mean 2) and the
  # mean is 7/3, so G = 2 / (2 * 7/3) = 3/7.
  expect_equal(gini_index(c(1, 2, 4)), 3 / 7)
  # 1, 2, 3, 4: mean absolute difference 10/6, mean 5/2, G = 1/3.
  expect_equal(gini_index(1:4), 1 / 3)
  # Complete equality and complete concentration.
  expect_identical(gini_index(rep(50, 10)), 0)
  expect_equal(gini_index(c(0, 0, 0, 5)), 1)
  # When one unit holds everything the index is exactly 1 and must not come
  # out above it: rounding used to give 1.0000000000000002 for the first
  # case here, and a value above 1 for about one in seven such samples.
  expect_lte(gini_index(c(0, 0, 0, 0, 0, 0.1)), 1)
  cases <- expand.grid(n = 2:30,
                       v = c(0.1, 0.3, 1 / 3, 0.7, 1, pi, 7, 13, 99, 1e-7, 1e9))
  g <- mapply(function(n, v) gini_index(c(rep(0, n - 1), v)), cases$n, cases$v)
  expect_true(all(g <= 1))
  expect_equal(g, rep(1, nrow(cases)))
  # Order does not matter.
  expect_equal(gini_index(c(4, 1, 2)), 3 / 7)
})

test_that("gini_index is the mean difference over twice the mean, n / (n - 1) times the plug in index", {
  withr::local_seed(113)
  x <- rexp(50)
  n <- length(x)
  mean_abs_diff <- sum(abs(outer(x, x, "-"))) / (n * (n - 1))
  expect_equal(gini_index(x), mean_abs_diff / (2 * mean(x)))
  plug_in <- sum(abs(outer(x, x, "-"))) / (2 * n^2 * mean(x))
  expect_equal(gini_index(x), n / (n - 1) * plug_in)
  expect_true(gini_index(x) > 0 && gini_index(x) < 1)
})

test_that("gini_index does not depend on the unit and does not overflow on integers", {
  withr::local_seed(113)
  x <- rlnorm(40, 3, 0.8)
  expect_equal(gini_index(1000 * x), gini_index(x), tolerance = 1e-13)
  # Very small positive values used to be rejected as "essentially zero".
  expect_equal(gini_index(1e-20 * c(1, 2, 4)), 3 / 7)
  # 1, ..., n has G = 1/3 for every n. As integers the sums overflow the
  # integer range from n of about 65,000; the old code returned NA.
  expect_silent(g <- gini_index(1:100000))
  expect_equal(g, 1 / 3)
  expect_silent(d <- gini_mean_diff(1:100000))
  expect_equal(d, 100001 / 3)
})

test_that("gini_index approaches the population values of Chattopadhyay and De (2016)", {
  # Their simulation populations, with the Gini indexes they report:
  # lognormal(2.185, 0.562), G = 0.3089 (in closed form
  # 2 pnorm(0.562 / sqrt(2)) - 1), and gamma(shape 2.649, rate 0.84),
  # G = 0.3308. The standard deviation of the estimate at n = 50,000 is
  # 0.0010 in both cases, so 0.0045 is more than four.
  expect_equal(2 * pnorm(0.562 / sqrt(2)) - 1, 0.3089, tolerance = 1e-3)
  withr::local_seed(113)
  expect_lt(abs(gini_index(rlnorm(50000, 2.185, 0.562)) - 0.3089), 0.0045)
  expect_lt(abs(gini_index(rgamma(50000, shape = 2.649, rate = 0.84)) -
                  0.3308), 0.0045)
})

test_that("gini_index checks its input", {
  expect_error(gini_index(), "'x' must be supplied")
  expect_error(gini_index(letters[1:5]), "'x' must be a numeric vector")
  expect_error(gini_index(matrix(1:6, 2)), "'x'.*not a matrix")
  # A plain estimator: the message says what is wrong and stops there (the
  # shared data check otherwise closes with advice for a sequential check).
  expect_error(gini_index(3),
               "^'x' needs at least 2 observations; it has 1\\.$")
  # A negative value used to give an "index" of 8.33.
  expect_error(gini_index(c(-5, 1, 2, 3)), "no negative values.*1 of 4")
  expect_error(gini_index(rep(0, 20)), "Every value in 'x' is zero")
  expect_error(gini_index(c(1, 2, Inf)), "infinite")
  expect_error(gini_index(1:4, na.rm = NA), "'na.rm' must be TRUE or FALSE")
})

test_that("gini_index and gini_mean_diff handle missing values as documented", {
  # Removed by default.
  expect_equal(gini_index(c(1, 2, NA, 4)), 3 / 7)
  expect_equal(gini_mean_diff(c(1, 2, NA, 4)), 2)
  # With na.rm = FALSE a missing value is an error. It used to return a
  # wrong number (0.3333 and 1.1667) with only a recycling warning.
  expect_error(gini_index(c(1, 2, NA, 4), na.rm = FALSE), "missing values")
  expect_error(gini_mean_diff(c(1, 2, NA, 4), na.rm = FALSE),
               "missing values")
  expect_equal(gini_index(c(1, 2, 4), na.rm = FALSE), 3 / 7)
  expect_error(gini_index(c(NA_real_, NA_real_, 3)), "at least 2")
})

test_that("gini_mean_diff equals hand computed values and the pairwise definition", {
  # 1, 2, 4: differences 1, 3, 2, mean 2. 1, 2, 3, 4: 10/6.
  expect_equal(gini_mean_diff(c(1, 2, 4)), 2)
  expect_equal(gini_mean_diff(1:4), 5 / 3)
  expect_identical(gini_mean_diff(rep(7, 5)), 0)
  # Never negative, and 0 for a sample of identical values whatever they
  # are (values whose sum is not exact in binary included).
  for (v in c(0.1, 1 / 3, 1e-300, 1e300, -2.7)) {
    expect_identical(gini_mean_diff(rep(v, 7)), 0)
    expect_identical(gini_index(rep(abs(v), 7)), 0)
  }
  withr::local_seed(113)
  x <- rnorm(60)                               # negative values are fine
  n <- length(x)
  expect_equal(gini_mean_diff(x), sum(abs(outer(x, x, "-"))) / (n * (n - 1)))
  # For the standard normal the population value is 2 / sqrt(pi) = 1.128;
  # the standard deviation of the estimate at n = 20,000 is about 0.005.
  expect_lt(abs(gini_mean_diff(rnorm(20000)) - 2 / sqrt(pi)), 0.025)
})

test_that("gini_mean_diff is location invariant and scale equivariant", {
  withr::local_seed(113)
  x <- round(rexp(80) * 1024) / 1024           # dyadic: the shift is exact
  expect_equal(gini_mean_diff(x + 1.76e9), gini_mean_diff(x),
               tolerance = 1e-9)
  expect_equal(gini_mean_diff(-3 * x), 3 * gini_mean_diff(x))
})

test_that("gini_mean_diff checks its input", {
  expect_error(gini_mean_diff(), "'x' must be supplied")
  expect_error(gini_mean_diff(letters[1:5]), "'x' must be a numeric vector")
  expect_error(gini_mean_diff(data.frame(a = 1:5)), "'x'.*not a data frame")
  expect_error(gini_mean_diff(3),
               "^'x' needs at least 2 observations; it has 1\\.$")
  expect_error(gini_mean_diff(c(1, -Inf, 2)), "infinite")
})

test_that("slr_slope equals a hand computed value and lm()", {
  # x = 1, 2, 3, 4 and y = 2, 4, 5, 9: sum of cross products 11, sum of
  # squares of x 5, slope 2.2.
  expect_equal(slr_slope(1:4, c(2, 4, 5, 9)), 2.2)
  withr::local_seed(113)
  x <- rnorm(40)
  y <- 1.5 + 0.7 * x + rnorm(40)
  expect_equal(slr_slope(x, y), unname(stats::coef(stats::lm(y ~ x))[2]),
               tolerance = 1e-10)
  # Two points determine a line.
  expect_equal(slr_slope(c(1, 3), c(10, 4)), -3)
})

test_that("slr_slope is right, or says why not, when the units are far apart", {
  # A slope of exactly zero is zero in any units. With x in units that make
  # it subnormal the ratio of the two scales is infinite, and zero times
  # that ratio used to be returned: NaN.
  expect_identical(slr_slope(c(1, 2, 3) * 1e-321, c(1, 1, 1)), 0)
  expect_identical(slr_slope(c(1, 2, 3, 4), c(2, 5, 5, 2)), 0)
  # A slope too large for a double is an error that names the units; it
  # used to be Inf or -Inf with no message.
  expect_error(slr_slope(c(1, 2, 4) * 1e-200, c(1, 2, 3) * 1e200),
               paste0("The slope is too large to represent\\. It is in units ",
                      "of 'y' per unit of 'x'.*largest is 3e\\+200.*largest ",
                      "is 4e-200.*Express 'x' or 'y' in other units"))
  expect_error(slr_slope(c(1, 2, 4) * 1e-200, -c(1, 2, 3) * 1e200),
               "The slope is too large to represent")
  # Units far apart are not in themselves a problem: the slope of
  # (1, 2, 4) on (1, 2, 3) is 4.5 / 7 * 3 = 9 / 14 (cross products 4.5,
  # squares 14 / 3), here times 1e300, and with the units exchanged times
  # 1e-320, a number that a double holds only as a subnormal.
  expect_equal(slr_slope(c(1, 2, 4) * 1e-150, c(1, 2, 3) * 1e150),
               9 / 14 * 1e300)
  expect_equal(slr_slope(c(1, 2, 4) * 1e300, c(1, 2, 3) * 1e-20) / 1e-320,
               9 / 14, tolerance = 1e-3)
  # The scales themselves can be further apart than a double can hold
  # (2^-1074 and 2^1000) while the slope is an ordinary number.
  expect_equal(.slr_slope_value(c(1, 2, 4) * 2^-1074, c(1, 2, 3) * 2^-1074),
               9 / 14)
  # The check function reports a slope it cannot represent as an estimate
  # that is not finite; it does not stop.
  out <- seq_aipe_reg_coef(c(1, 2, 4, 7) * 1e-200, c(1, 2, 3, 5) * 1e200,
                           width = 1)
  expect_true(is.na(out$value[out$term == "reg_coef"]))
  expect_false(attr(out, "is_satisfied"))
})

test_that("slr_slope removes incomplete pairs, or refuses them", {
  x <- c(1, 2, NA, 3, 4, 7)
  y <- c(2, 4, 8, NA, 5, 9)
  # Complete pairs: (1, 2), (2, 4), (4, 5), (7, 9).
  expect_equal(slr_slope(x, y), slr_slope(c(1, 2, 4, 7), c(2, 4, 5, 9)))
  expect_error(slr_slope(x, y, na.rm = FALSE), "'x' contains missing values")
  expect_error(slr_slope(c(1, 2, 3), c(1, NA, 2), na.rm = FALSE),
               "'y' contains missing values")
  # Too few complete pairs. (No advice about a stopping rule: this is a
  # plain estimator.)
  expect_error(slr_slope(c(1, NA, 3), c(NA, 2, 3)),
               "^'x' needs at least 2 observations; it has 1\\.$")
})

test_that("slr_slope checks its input", {
  expect_error(slr_slope(y = 1:4), "'x' must be supplied")
  expect_error(slr_slope(1:4), "'y' must be supplied")
  expect_error(slr_slope(1:10, 1:9), "same length.*10 and 9")
  expect_error(slr_slope(letters[1:4], 1:4), "'x' must be a numeric vector")
  # A matrix x used to be recycled against y and return a wrong slope (1)
  # with only a warning.
  expect_error(slr_slope(matrix(1:8, ncol = 2) * 1.0,
                         c(1, 3, 2, 5, 4, 7, 6, 9)), "'x'.*not a matrix")
  expect_error(slr_slope(1:4, data.frame(y = 1:4)), "'y'.*not a data frame")
  expect_error(slr_slope(c(1, 2, Inf), c(1, 2, 3)), "infinite")
  # A predictor with no variation used to return NaN silently.
  expect_error(slr_slope(rep(1, 5), 1:5), "'x' has no variation")
  # So is a predictor whose values are equal on paper but not in floating
  # point (one dose of 0.3 obtained as a difference): its "slope" would be
  # rounding noise, 2.9e14 here. The package wide test for no variation
  # (agreement to ten significant digits) catches it, and the message says
  # what to do when the variation is real.
  dose <- c(0.3, 0.3, 0.3, (63.9 + 0.3) - 63.9, 0.3)
  expect_gt(max(dose) - min(dose), 0)
  expect_equal(cov(dose, 1:5) / var(dose), 2.924e14, tolerance = 1e-3)
  expect_error(slr_slope(dose, 1:5),
               "'x' has no variation.*ten significant digits.*subtract a constant")
  # Real variation that is small beside the size of the values is a slope
  # once the predictor is centered: 1e9 + 0, 0.01, ..., 0.04 against 1:5 is
  # 100 per unit.
  expect_error(slr_slope(1e9 + (0:4) / 100, 1:5), "'x' has no variation")
  expect_equal(slr_slope((0:4) / 100, 1:5), 100)
  expect_error(slr_slope(1, 2), "at least 2")
  expect_error(slr_slope(1:4, 1:4, na.rm = "yes"), "'na.rm' must be TRUE")
})
