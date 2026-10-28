# Tests for what the seq_aipe_* family shares through R/aipe_internals.R:
# the normal quantile, the refusal of a criterion that cannot be
# represented, and the exact power of two scaling that keeps every
# estimator right whatever the units of the data.

row_value <- function(out, term) out$value[out$term == term]

test_that(".aipe_z is the upper alpha / 2 normal quantile", {
  # Textbook values.
  expect_equal(.aipe_z(0.90), 1.644854, tolerance = 1e-6)
  expect_equal(.aipe_z(0.95), 1.959964, tolerance = 1e-6)
  expect_equal(.aipe_z(0.99), 2.575829, tolerance = 1e-6)
  # The round trip through the normal distribution function holds to
  # working precision however close conf_level is to 1. The form
  # qnorm(1 - alpha / 2) loses accuracy there (8.0414 instead of 8.0269 at
  # 1 - 1e-15) and is infinite at the largest number below 1.
  for (cl in c(0.5, 0.999, 1 - 1e-9, 1 - 1e-12, 1 - 1e-15, 1 - 2^(-53))) {
    z <- .aipe_z(cl)
    expect_true(is.finite(z))
    expect_equal(2 * pnorm(z, lower.tail = FALSE), 1 - cl, tolerance = 1e-10)
  }
  expect_equal(.aipe_z(1 - 1e-15), 8.0269, tolerance = 1e-4)
})

test_that("a confidence level just below 1 gives finite pilots and finite checks", {
  cl <- 1 - 2^(-53)                            # the largest double below 1
  z <- .aipe_z(cl)
  for (pilot in list(seq_aipe_smd_pilot, seq_aipe_r_pilot, seq_aipe_cv_pilot,
                     seq_aipe_reg_coef_pilot, seq_aipe_sm_pilot,
                     seq_aipe_c_pilot, seq_aipe_gini_pilot)) {
    expect_identical(pilot(width = 0.3, conf_level = cl),
                     as.integer(ceiling(2 * z / 0.3)))
  }
  withr::local_seed(113)
  x <- rnorm(30, mean = 1)
  y <- 0.5 * x + rnorm(30)
  checks <- list(
    seq_aipe_sm(x, 0.3, conf_level = cl),
    seq_aipe_smd(x, y, 0.3, conf_level = cl),
    seq_aipe_r(x, y, 0.3, conf_level = cl),
    seq_aipe_reg_coef(x, y, 0.3, conf_level = cl),
    seq_aipe_cv(x + 5, 0.3, conf_level = cl),
    seq_aipe_gini(x + 5, 0.3, conf_level = cl),
    seq_aipe_c(list(x, y), c(1, -1), 0.3, conf_level = cl)
  )
  for (out in checks) {
    # Previously the limits were -Inf and Inf and the criterion Inf.
    expect_true(all(is.finite(out$value)))
    half <- (row_value(out, "upper_limit") - row_value(out, "lower_limit")) / 2
    expect_equal(half, z * sqrt(row_value(out, "xi_squared") / 30))
  }
})

test_that("a width too small for the criterion to be represented is an error", {
  withr::local_seed(113)
  x <- rnorm(30, mean = 1)
  y <- 0.5 * x + rnorm(30)
  msg <- paste0("sample size required by the stopping rule \\(Inf\\) is too ",
                "large to represent\\. Check 'width'")
  # (2 z / width)^2 overflows. Previously criterion_n was Inf.
  expect_error(seq_aipe_smd(x, y, width = 1e-200), msg)
  expect_error(seq_aipe_r(x, y, width = 1e-200), msg)
  expect_error(seq_aipe_r(x, y, width = 1e-200, method = "kendall"), msg)
  expect_error(seq_aipe_reg_coef(x, y, width = 1e-200), msg)
  expect_error(seq_aipe_sm(x, width = 1e-200), msg)
  expect_error(seq_aipe_cv(x + 5, width = 1e-200), msg)
  expect_error(seq_aipe_gini(x + 5, width = 1e-200), msg)
  expect_error(seq_aipe_c(list(x, y), c(1, -1), width = 1e-200), msg)
  # A criterion that is very large but finite is still reported: with a
  # width of 1e-5 it is (2 z / width)^2 (xi2 + 1 / n), beyond the range of an
  # integer but an ordinary double.
  out <- seq_aipe_sm(x, width = 1e-5)
  crit <- ceiling((2 * qnorm(0.975) / 1e-5)^2 *
                    (row_value(out, "xi_squared") + 1 / 30))
  expect_gt(crit, .Machine$integer.max)
  expect_equal(row_value(out, "criterion_n"), crit)
  expect_false(attr(out, "is_satisfied"))
})

test_that(".aipe_scale is a power of two of the order of the largest absolute value", {
  expect_identical(.aipe_scale(c(3, -5)), 4)
  expect_identical(.aipe_scale(0.3), 0.25)
  expect_identical(.aipe_scale(c(1, 2), 9), 8)       # several vectors at once
  expect_identical(.aipe_scale(c(0, 0, 0)), 1)       # nothing to scale by
  below <- 1 - .Machine$double.eps / 2               # the double before 1
  for (v in c(4.9e-324, 1e-300, 1e-17, 1, 7, 8 * below, 8, 1e17, 1e300,
              2^1023 * below, 2^1023, .Machine$double.xmax)) {
    s <- .aipe_scale(c(-v, v / 3))
    expect_true(is.finite(s) && s > 0)
    expect_identical(s, 2^round(log2(s)))            # a power of two
    expect_gte(v / s, 1)
    expect_lt(v / s, 2)
  }
  # log2() of a value within rounding of a power of two is that power's
  # exponent, which would give a scale above the value (and an infinite one
  # for the largest double, whose base 2 logarithm is returned as 1024).
  expect_identical(.aipe_scale(8 * below), 4)
  expect_identical(.aipe_scale(.Machine$double.xmax), 2^1023)
  # A first guess of the exponent that is one too high or one too low is
  # put right; a correct one is left alone.
  expect_identical(.aipe_power_of_two(8 * below, e = 3), 4)
  expect_identical(.aipe_power_of_two(8, e = 2), 8)
  expect_identical(.aipe_power_of_two(8, e = 3), 8)
  expect_identical(.aipe_power_of_two(4.9e-324, e = -1075), 2^(-1074))
  expect_equal(gini_index(c(1, 2, 3, 4) * (.Machine$double.xmax / 4)), 1 / 3)
})

test_that("the scale free check functions give the same answer at any magnitude", {
  # Before the data were scaled, the sixth power of the pooled standard
  # deviation overflowed for data of magnitude 1e52 (in one sample the
  # variance of the standardized mean difference came out as 2.11 where the
  # value is 2.17, with no flag), and ordinary data were reported as
  # "variance not yet usable" by three of the procedures at a magnitude of
  # 1e-60 and by all of them at 1e80.
  withr::local_seed(113)
  a <- rnorm(40, mean = 1)
  b <- rnorm(40, mean = 0.5)
  u <- rnorm(40)
  v <- 0.5 * u + rnorm(40)
  amount <- rgamma(40, shape = 4)
  base <- list(
    smd = seq_aipe_smd(a, b, 0.3),
    pearson = seq_aipe_r(u, v, 0.3),
    kendall = seq_aipe_r(u, v, 0.3, method = "kendall"),
    spearman = seq_aipe_r(u, v, 0.3, method = "spearman"),
    cv = seq_aipe_cv(amount, 0.3),
    sm = seq_aipe_sm(a, 0.3),
    gini = seq_aipe_gini(amount, 0.3),
    slope = seq_aipe_reg_coef(u, v, 0.3)
  )
  for (out in base) expect_false(attr(out, "variance_unusable"))
  for (s in c(1e-300, 1e-150, 1e-60, 1e-41, 1e39, 1e52, 1e150, 1e300)) {
    scaled <- list(
      smd = seq_aipe_smd(a * s, b * s, 0.3),
      # The two variables of a correlation need not share a scale.
      pearson = seq_aipe_r(u * s, v / s, 0.3),
      kendall = seq_aipe_r(u * s, v / s, 0.3, method = "kendall"),
      spearman = seq_aipe_r(u * s, v / s, 0.3, method = "spearman"),
      cv = seq_aipe_cv(amount * s, 0.3),
      sm = seq_aipe_sm(a * s, 0.3),
      gini = seq_aipe_gini(amount * s, 0.3),
      # A slope is unchanged when both variables are rescaled alike.
      slope = seq_aipe_reg_coef(u * s, v * s, 0.3)
    )
    for (nm in names(base)) {
      expect_equal(scaled[[nm]]$value, base[[nm]]$value, tolerance = 1e-9,
                   info = paste(nm, "at scale", s))
      expect_identical(attr(scaled[[nm]], "is_satisfied"),
                       attr(base[[nm]], "is_satisfied"))
    }
  }
  # Degenerate data produce no warning at these magnitudes either (cor()
  # used to warn that the standard deviation is zero at 1e-300).
  expect_silent(seq_aipe_r(u * 1e-300, v * 1e-300, 0.3))
})

test_that("each scale free variance estimator is right at any magnitude on its own", {
  # The estimators are also called outside the check functions
  # (seq_mr_cv() uses .var_cv()), so each scales its own data.
  withr::local_seed(113)
  a <- rnorm(40, mean = 1)
  b <- rnorm(40, mean = 0.5)
  amount <- rgamma(40, shape = 4)
  for (s in c(1e-150, 1e-60, 1e52, 1e150)) {
    expect_equal(.var_smd(a * s, b * s), .var_smd(a, b), tolerance = 1e-9)
    expect_equal(.var_cv(amount * s), .var_cv(amount), tolerance = 1e-9)
    expect_equal(.var_sm(a * s), .var_sm(a), tolerance = 1e-9)
    expect_equal(.var_gini(amount * s), .var_gini(amount), tolerance = 1e-9)
    expect_equal(.var_r(a * s, b / s), .var_r(a, b), tolerance = 1e-9)
    expect_equal(.var_reg_coef(a * s, b * s), .var_reg_coef(a, b),
                 tolerance = 1e-9)
  }
  # None of these is a value replaced by the bound.
  for (v in c(.var_smd(a, b), .var_cv(amount), .var_sm(a), .var_gini(amount),
              .var_r(a, b), .var_reg_coef(a, b))) {
    expect_gt(v, 40^(-3))
  }
})

test_that("the estimators that have units carry them exactly", {
  withr::local_seed(113)
  x <- rnorm(40)
  y <- 0.5 * x + rnorm(40)
  # Powers of two, so that the rescaling of the data is itself exact. The
  # slope is in units of y per unit of x and its variance in their square.
  fx <- 2^200
  fy <- 2^(-150)
  expect_identical(slr_slope(x * fx, y * fy), slr_slope(x, y) * fy / fx)
  expect_identical(.var_reg_coef(x * fx, y * fy, floor = FALSE),
                   .var_reg_coef(x, y, floor = FALSE) * (fy / fx)^2)
  # Hand value at a magnitude where the sum of squares of x overflows: x =
  # 1, 2, 3, 4 and y = 2, 4, 5, 9 have slope 2.2.
  expect_equal(slr_slope(c(1, 2, 3, 4) * 1e200, c(2, 4, 5, 9) * 1e-100),
               2.2e-300)
  # The contrast and its variance: units of the outcome and their square.
  g <- list(c(1, 2, 3, 4), c(2, 4, 6, 8), c(3, 3, 5, 9))
  cw <- c(1, -0.5, -0.5)
  big <- lapply(g, function(v) v * 2^400)
  raw_c <- function(groups, weights) .var_c(groups, weights, floor = FALSE)
  expect_equal(raw_c(g, cw), 16 / 3)                 # 5/3 + 20/12 + 8/4
  expect_identical(raw_c(big, cw), raw_c(g, cw) * 2^800)
  out <- seq_aipe_c(big, cw, width = 2^402)
  expect_identical(row_value(out, "contrast"), -2.5 * 2^400)
  expect_identical(row_value(out, "xi_squared"), raw_c(g, cw) * 2^800)
  # Near the top of the range of a double the sum of squares inside a
  # variance can overflow although the variance itself does not: the
  # variance of 1, 2, 3, 4 times 2^511 is (5 / 3) 2^1022, while its sum of
  # squared deviations, 5 * 2^1022, is beyond the largest double (on a
  # platform without extended precision var() returns Inf).
  high <- list(c(1, 2, 3, 4) * 2^511, c(1, 2, 3, 4) * 2^511)
  expect_equal(raw_c(high, c(1, -1)) / 2^1022, 10 / 3)
  # Gini's mean difference of 1, 1.2, 1.5, 1.7 is 2.4 / 6 = 0.4, and their
  # Gini index 0.4 / (2 * 1.35); the sums overflowed at this magnitude
  # (the index came back NaN and the mean difference Inf).
  top <- c(1, 1.2, 1.5, 1.7) * 1e308
  expect_equal(gini_mean_diff(top), 0.4e308)
  expect_equal(gini_index(top), 0.4 / 2.7)
  expect_equal(gini_index(c(1, 2, 3, 4) * 1e307), 1 / 3)
  expect_equal(gini_index(c(1, 2, 4) * 1e-310), 3 / 7, tolerance = 1e-9)
})
