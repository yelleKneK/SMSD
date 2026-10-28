# Tests for seq_aipe_reg_coef() and seq_aipe_reg_coef_pilot(): the
# sequential AIPE procedure for the slope of a simple linear regression of
# Kelley, Darku, and Chattopadhyay (2018).

row_value <- function(out, term) out$value[out$term == term]

test_that("seq_aipe_reg_coef_pilot is max{m_0, ceiling(2 z / width)} (2018, Equation 21)", {
  expect_identical(seq_aipe_reg_coef_pilot(width = 0.20), 20L)
  expect_identical(seq_aipe_reg_coef_pilot(width = 0.40), 10L)
  expect_identical(seq_aipe_reg_coef_pilot(width = 0.20, conf_level = 0.99),
                   as.integer(ceiling(2 * qnorm(0.995) / 0.20)))   # 26
  expect_identical(seq_aipe_reg_coef_pilot(width = 0.20, m_0 = 30), 30L)
})

test_that("seq_aipe_reg_coef_pilot checks its arguments", {
  expect_pilot_checks(seq_aipe_reg_coef_pilot)
})

test_that("seq_aipe_reg_coef returns the documented rows and attributes", {
  withr::local_seed(113)
  x <- rnorm(40)
  y <- 1.5 + 0.7 * x + rnorm(40)
  out <- seq_aipe_reg_coef(x, y, width = 0.50)
  expect_s3_class(out, c("smsd_check", "data.frame"), exact = TRUE)
  expect_identical(out$term, c("lower_limit", "reg_coef", "upper_limit",
                               "xi_squared", "criterion_n", "current_n"))
  expect_true(all(is.finite(out$value)))
  expect_identical(attr(out, "method"), "simple linear regression slope")
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_identical(attr(out, "notes"), character(0))
  expect_identical(attr(out, "conf_level"), 0.95)
  expect_identical(attr(out, "width"), 0.50)
  expect_equal(attr(out, "current_n"), 40)
  expect_silent(seq_aipe_reg_coef(x, y, width = 0.50))
})

test_that("seq_aipe_reg_coef computes the slope, the interval (2018, Equation 38), and the criterion (Equation 22)", {
  withr::local_seed(113)
  x <- rnorm(40)
  y <- 1.5 + 0.7 * x + rnorm(40)
  for (cl in c(0.90, 0.95)) {
    out <- seq_aipe_reg_coef(x, y, width = 0.50, conf_level = cl)
    b <- unname(coef(lm(y ~ x))[2])
    xi2 <- .var_reg_coef(x, y)
    z <- qnorm(1 - (1 - cl) / 2)
    expect_equal(row_value(out, "reg_coef"), b)
    expect_equal(row_value(out, "xi_squared"), xi2)
    expect_equal(row_value(out, "lower_limit"), b - z * sqrt(xi2 / 40))
    expect_equal(row_value(out, "upper_limit"), b + z * sqrt(xi2 / 40))
    expect_equal(row_value(out, "criterion_n"),
                 ceiling((2 * z / 0.50)^2 * (xi2 + 1 / 40)))
    expect_identical(attr(out, "is_satisfied"),
                     40 >= row_value(out, "criterion_n"))
  }
})

test_that("seq_aipe_reg_coef meets the rule exactly when n reaches the criterion", {
  withr::local_seed(113)
  x <- rnorm(500)
  y <- 1 + 0.5 * x + rnorm(500)
  wide <- seq_aipe_reg_coef(x, y, width = 1.50)
  expect_true(attr(wide, "is_satisfied"))
  expect_lte(row_value(wide, "upper_limit") - row_value(wide, "lower_limit"),
             1.50)
  expect_false(attr(seq_aipe_reg_coef(x[1:20], y[1:20], width = 0.05),
                    "is_satisfied"))
})

test_that("seq_aipe_reg_coef is unchanged by shifting either variable", {
  # With raw power sums the variance estimate went from 0.87 to 62 when both
  # variables were shifted by 1e4.
  withr::local_seed(113)
  x <- round(rnorm(80) * 1024) / 1024
  y <- round((0.5 * x + rnorm(80)) * 1024) / 1024
  base <- seq_aipe_reg_coef(x, y, width = 0.4)
  moved <- seq_aipe_reg_coef(x + 1e6, y - 1e4, width = 0.4)
  expect_equal(moved$value, base$value, tolerance = 1e-8)
  expect_identical(attr(moved, "is_satisfied"), attr(base, "is_satisfied"))
  # Rescaling x by 10 divides the slope by 10 and its variance by 100.
  scaled <- seq_aipe_reg_coef(10 * x, y, width = 0.4)
  expect_equal(row_value(scaled, "reg_coef"), row_value(base, "reg_coef") / 10)
  expect_equal(row_value(scaled, "xi_squared"),
               row_value(base, "xi_squared") / 100)
})

test_that("a positive variance estimate below n^-3 is used through the bound, with a note about units", {
  # The slope has units; the bound n^-3 has none. Here is one study with
  # the predictor recorded in hours and in seconds: the slope and the width
  # divide by 3600 and the variance estimate by 3600^2, which takes it below
  # n^-3 (with xi^2 = 1 in hours the estimate in seconds stays below n^-3
  # until n is about 235). The bound then stands in for the estimate, as
  # the article prescribes, the rule is evaluated, and the interval is
  # built; a note says that the bound was used and that the units are the
  # usual reason. (An earlier version reported such a check as "not yet
  # usable", with NA limits.)
  withr::local_seed(113)
  hours <- rnorm(200)
  y <- 0.5 * hours + rnorm(200)
  seconds <- 3600 * hours
  z <- qnorm(0.975)
  for (n in c(60, 200)) {
    i <- seq_len(n)
    in_hours <- seq_aipe_reg_coef(hours[i], y[i], width = 0.36)
    in_seconds <- seq_aipe_reg_coef(seconds[i], y[i], width = 0.36 / 3600)
    expect_false(attr(in_hours, "variance_unusable"))
    expect_identical(attr(in_hours, "notes"), character(0))
    raw <- row_value(in_hours, "xi_squared") / 3600^2
    expect_gt(raw, 0)
    expect_lt(raw, n^(-3))
    expect_equal(.var_reg_coef(seconds[i], y[i], floor = FALSE), raw,
                 tolerance = 1e-10)

    expect_false(attr(in_seconds, "variance_unusable"))
    expect_identical(row_value(in_seconds, "xi_squared"), n^(-3))
    b <- row_value(in_hours, "reg_coef") / 3600
    expect_equal(row_value(in_seconds, "reg_coef"), b)
    expect_equal(row_value(in_seconds, "lower_limit"),
                 b - z * sqrt(n^(-3) / n))
    expect_equal(row_value(in_seconds, "upper_limit"),
                 b + z * sqrt(n^(-3) / n))
    # The criterion is the article's with the bound: far above n here,
    # because 1 / n is huge beside a variance of this size.
    crit <- ceiling((2 * z * 3600 / 0.36)^2 * (n^(-3) + 1 / n))
    expect_equal(row_value(in_seconds, "criterion_n"), crit)
    expect_gt(crit, 1e6)
    expect_false(attr(in_seconds, "is_satisfied"))
    note <- attr(in_seconds, "notes")
    expect_length(note, 1L)
    expect_match(note, "positive but below its lower bound")
    expect_match(note, "the bound was used in its place")
    expect_match(note, "Rescale x or y, and width with it")
    expect_match(note, "seq_aipe_reg_coef")
    expect_false(grepl("not yet usable", note))
    expect_match(paste(format(in_seconds), collapse = "\n"), "Rescale x or y")
  }

  # The note informs; it does not change the decision. With a target wide
  # enough for the rule with the bound to be met, it is met, and the note
  # is still there: (2 z / 0.8)^2 (60^-3 + 1/60) = 0.4, so the criterion
  # is 1.
  i <- seq_len(60)
  wide <- seq_aipe_reg_coef(seconds[i], y[i], width = 0.8)
  expect_identical(row_value(wide, "xi_squared"), 60^(-3))
  expect_equal(row_value(wide, "criterion_n"), 1)
  expect_true(attr(wide, "is_satisfied"))
  expect_length(attr(wide, "notes"), 1L)
  expect_match(attr(wide, "notes"), "the bound was used in its place")

  # A raw estimate that is zero or negative, the usual reason for the bound
  # in a small sample, is not an estimate at all: the check is "not yet
  # usable" with the standard note. Pilot samples of 10 pairs are drawn
  # until one has a negative raw estimate (about one in six does).
  found <- FALSE
  for (i in 1:200) {
    x <- rnorm(10)
    yy <- 0.5 * x + rnorm(10)
    if (.var_reg_coef(x, yy, floor = FALSE) < 0) {
      found <- TRUE
      break
    }
  }
  expect_true(found)
  out <- seq_aipe_reg_coef(x, yy, width = 0.4)
  expect_true(attr(out, "variance_unusable"))
  expect_false(attr(out, "is_satisfied"))
  expect_true(is.na(row_value(out, "lower_limit")))
  expect_identical(attr(out, "notes"), .note_variance_unusable())
  # So is a constant response (raw estimate exactly zero).
  expect_identical(.var_reg_coef(x, rep(2, 10), floor = FALSE), 0)
  flat <- seq_aipe_reg_coef(x, rep(2, 10), width = 0.4)
  expect_true(attr(flat, "variance_unusable"))
  expect_identical(attr(flat, "notes"), .note_variance_unusable())
  # A predictor with no variation has no raw estimate at all.
  expect_true(is.nan(.var_reg_coef(rep(2, 10), yy, floor = FALSE)))
})

test_that("seq_aipe_reg_coef checks its arguments", {
  withr::local_seed(113)
  good <- rnorm(20)
  other <- rnorm(20)
  expect_error(seq_aipe_reg_coef(y = good, width = 0.3),
               "'x' must be supplied")
  expect_error(seq_aipe_reg_coef(good, width = 0.3), "'y' must be supplied")
  expect_data_checks(function(v) seq_aipe_reg_coef(v, other, width = 0.3),
                     "x", good, short = FALSE)
  expect_data_checks(function(v) seq_aipe_reg_coef(other, v, width = 0.3),
                     "y", good, short = FALSE)
  expect_error(seq_aipe_reg_coef(good, other[1:15], width = 0.50),
               "same length.*20 and 15")
  # The whole message: a sequential check closes with what to do next (the
  # plain estimator slr_slope(), which shares the pair check, does not).
  expect_error(seq_aipe_reg_coef(good[1:3], other[1:3], width = 0.5),
               paste0("^'x' needs at least 4 observations; it has 3\\. ",
                      "Collect more data, then check the stopping rule ",
                      "again\\.$"))
  expect_scalar_checks(function(...) seq_aipe_reg_coef(good, other, ...))
})

test_that("operating characteristics: the slope procedure reaches its target sample size", {
  # No published table exists for the slope, so the target is the
  # theoretical sample size. Normal x and errors with unit variances and a
  # slope of 0.5: xi^2 = 1, so at width 0.4 the theoretical sample size is
  # ceiling((2 * 1.96 / 0.4)^2) = 97. The minimum pilot is set to 50
  # because with the default (a pilot of 10 here) the variance estimate is
  # often unusable at the pilot (see test-var_floor.R and the help page). In
  # 4,000 replications of this condition the mean final sample size was
  # 90.3 (standard deviation 27.2; the estimator runs low at small n, so
  # the procedure stops short of the theoretical size on average) and
  # coverage 0.919 (se 0.004), below the nominal 0.95, as the help page
  # reports. The mean of the 40 replications here had a standard deviation
  # of at most 4.48 over 200 and 300 other seeds, so the window around 97
  # is 4 * 4.48 plus the 6.7 by which the procedure's own mean falls short
  # of it: 24.6.
  z <- qnorm(0.975)
  expect_equal(ceiling((2 * z / 0.4)^2), 97)
  withr::local_seed(113)
  w <- 0.4
  take <- function(d, n) list(x = d$x[1:n], y = d$y[1:n])
  stat <- function(d) {
    c(cov(d$x, d$y) / var(d$x), .var_reg_coef(d$x, d$y, floor = FALSE))
  }
  check <- function(d) seq_aipe_reg_coef(d$x, d$y, width = w)
  runs <- replicate(40, {
    x <- rnorm(600)
    d <- list(x = x, y = 0.5 * x + rnorm(600))
    aipe_run(d, take, stat, check, width = w,
             pilot = seq_aipe_reg_coef_pilot(w, m_0 = 50))
  })
  s <- aipe_summary(runs, truth = 0.5)
  expect_lt(abs(s[["mean_n"]] - 97), 24.6)
  expect_lte(s[["max_width"]], w)
  expect_gte(min(runs["n", ]), 50)             # never below the pilot
})
