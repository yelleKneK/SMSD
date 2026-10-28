# Tests for seq_aipe_smd() and seq_aipe_smd_pilot(): the sequential AIPE
# procedure for the standardized mean difference of Kelley, Darku, and
# Chattopadhyay (2018), with the corrected Equation 26.

row_value <- function(out, term) out$value[out$term == term]

test_that("seq_aipe_smd_pilot is max{m_0, ceiling(2 z / width)} (2018, Equation 21)", {
  # 2 * 1.959964 / 0.2 = 19.6, so 20; / 0.3 = 13.07, so 14.
  expect_identical(seq_aipe_smd_pilot(width = 0.2), 20L)
  expect_identical(seq_aipe_smd_pilot(width = 0.3), 14L)
  # 90%: 2 * 1.644854 / 0.2 = 16.4, so 17.
  expect_identical(seq_aipe_smd_pilot(width = 0.2, conf_level = 0.90), 17L)
  # The minimum applies when it is larger, and only then.
  expect_identical(seq_aipe_smd_pilot(width = 0.2, m_0 = 50), 50L)
  expect_identical(seq_aipe_smd_pilot(width = 0.2, m_0 = 10), 20L)
  expect_identical(seq_aipe_smd_pilot(width = 5), 4L)
  # Integer valued m_0 stored as integer or double.
  expect_identical(seq_aipe_smd_pilot(width = 5, m_0 = 6L), 6L)
})

test_that("seq_aipe_smd_pilot checks its arguments", {
  expect_pilot_checks(seq_aipe_smd_pilot)
})

test_that("seq_aipe_smd returns the documented rows and attributes", {
  withr::local_seed(113)
  group_1 <- rnorm(20, mean = 0.5)
  group_2 <- rnorm(20)
  out <- seq_aipe_smd(group_1, group_2, width = 0.30)

  expect_s3_class(out, c("smsd_check", "data.frame"), exact = TRUE)
  expect_identical(out$term, c("lower_limit", "smd", "upper_limit",
                               "xi_squared", "criterion_n", "current_n"))
  expect_type(out$value, "double")
  expect_true(all(is.finite(out$value)))
  expect_identical(attr(out, "method"), "standardized mean difference")
  expect_identical(attr(out, "is_satisfied"), FALSE)
  expect_equal(attr(out, "current_n"), 20)
  expect_equal(attr(out, "criterion_n"), row_value(out, "criterion_n"))
  expect_identical(attr(out, "conf_level"), 0.95)
  expect_identical(attr(out, "width"), 0.30)
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_identical(attr(out, "notes"), character(0))
  expect_silent(seq_aipe_smd(group_1, group_2, width = 0.30))
  capture.output(expect_invisible(print(out)))
})

test_that("seq_aipe_smd computes d, the interval, and the criterion as published", {
  withr::local_seed(113)
  group_1 <- rnorm(40, mean = 0.5)
  group_2 <- rnorm(40)
  for (cl in c(0.90, 0.95)) {
    out <- seq_aipe_smd(group_1, group_2, width = 0.50, conf_level = cl)
    # d with the pooled variance of two equal groups (2018, Equation 24).
    d <- (mean(group_1) - mean(group_2)) /
      sqrt((var(group_1) + var(group_2)) / 2)
    expect_equal(row_value(out, "smd"), d)
    xi2 <- row_value(out, "xi_squared")
    expect_equal(xi2, .var_smd(group_1, group_2))
    z <- qnorm(1 - (1 - cl) / 2)
    # Interval (Equation 28) and criterion (Equation 22).
    expect_equal(row_value(out, "lower_limit"), d - z * sqrt(xi2 / 40))
    expect_equal(row_value(out, "upper_limit"), d + z * sqrt(xi2 / 40))
    expect_equal(row_value(out, "criterion_n"),
                 ceiling((2 * z / 0.50)^2 * (xi2 + 1 / 40)))
    expect_identical(attr(out, "is_satisfied"),
                     40 >= row_value(out, "criterion_n"))
  }
  # Swapping the groups changes the sign of d and nothing else.
  a <- seq_aipe_smd(group_1, group_2, width = 0.5)
  b <- seq_aipe_smd(group_2, group_1, width = 0.5)
  expect_equal(row_value(b, "smd"), -row_value(a, "smd"))
  expect_equal(row_value(b, "xi_squared"), row_value(a, "xi_squared"))
})

test_that("seq_aipe_smd meets the rule exactly when n reaches the criterion", {
  withr::local_seed(113)
  group_1 <- rnorm(500)
  group_2 <- rnorm(500)
  wide <- seq_aipe_smd(group_1, group_2, width = 1.50)
  expect_true(attr(wide, "is_satisfied"))
  expect_lte(row_value(wide, "upper_limit") - row_value(wide, "lower_limit"),
             1.50)
  narrow <- seq_aipe_smd(group_1[1:10], group_2[1:10], width = 0.10)
  expect_false(attr(narrow, "is_satisfied"))
})

test_that("seq_aipe_smd is unchanged by a common shift and a common rescaling", {
  withr::local_seed(113)
  g1 <- round(rnorm(60, 0.5) * 1024) / 1024
  g2 <- round(rnorm(60) * 1024) / 1024
  base <- seq_aipe_smd(g1, g2, width = 0.4)
  shifted <- seq_aipe_smd(g1 + 1e6, g2 + 1e6, width = 0.4)
  scaled <- seq_aipe_smd(16 * g1, 16 * g2, width = 0.4)
  expect_equal(shifted$value, base$value, tolerance = 1e-9)
  expect_equal(scaled$value, base$value, tolerance = 1e-12)
  expect_identical(attr(shifted, "is_satisfied"), attr(base, "is_satisfied"))
})

test_that("seq_aipe_smd accepts integer data", {
  g1 <- c(4L, 6L, 5L, 7L, 9L, 3L, 6L, 8L)
  g2 <- c(2L, 5L, 3L, 4L, 6L, 2L, 5L, 3L)
  expect_equal(seq_aipe_smd(g1, g2, width = 1)$value,
               seq_aipe_smd(as.numeric(g1), as.numeric(g2), width = 1)$value)
})

test_that("one group with no variation makes the variance not yet usable", {
  # Engagement ratings of 10 trained and 10 control employees; every
  # control rating is 3. The formula has a value (the pooled variance is
  # half the variance of the trained group), and at a width of 2.6 the
  # rule as printed is met with it; but the variance of the control group
  # has not been estimated, the pooled standard deviation is too small, and
  # the standardized mean difference too large.
  trained <- c(4, 5, 3, 4, 4, 5, 3, 4, 5, 4)
  control <- rep(3, 10)
  # By hand: means 4.1 and 3; variances 4.9 / 9 and 0, pooled 4.9 / 18.
  d <- (4.1 - 3) / sqrt(4.9 / 18)
  expect_equal(d, 2.108293, tolerance = 1e-6)
  out <- seq_aipe_smd(trained, control, width = 2.6)
  expect_equal(row_value(out, "smd"), d)
  expect_false(attr(out, "is_satisfied"))
  expect_true(attr(out, "variance_unusable"))
  expect_true(is.na(row_value(out, "lower_limit")))
  expect_true(is.na(row_value(out, "upper_limit")))
  expect_identical(row_value(out, "xi_squared"), 10^(-3))
  expect_identical(attr(out, "criterion_n"), 11)
  expect_identical(attr(out, "notes")[1L], .note_variance_unusable())
  expect_identical(
    attr(out, "notes")[2L],
    paste("The observations in 'group_2' are all equal, or differ too",
          "little relative to their size for a variance to be computed, so",
          "the variance of that group has not been estimated. If a group",
          "shows no variation well past the pilot sample, reconsider the",
          "measure (a ceiling or a floor).")
  )
  # The other way around the note names the first group, and the sign of
  # the estimate changes.
  swapped <- seq_aipe_smd(control, trained, width = 2.6)
  expect_equal(row_value(swapped, "smd"), -d)
  expect_true(attr(swapped, "variance_unusable"))
  expect_match(attr(swapped, "notes")[2L], "^The observations in 'group_1'")
  # The budget procedure reads the same two groups the same way.
  expect_true(attr(seq_budget(list(trained, control), c(1, -1), c(1, 1),
                              budget = 100), "variance_unusable"))
  # Neither group varies: no estimate at all, and the standard notes alone.
  neither <- seq_aipe_smd(control + 1, control, width = 2.6)
  expect_true(is.na(row_value(neither, "smd")))
  expect_identical(attr(neither, "notes"),
                   c(.note_not_estimable(), .note_variance_unusable()))
  # One control rating that differs and the check is an ordinary one.
  control[10] <- 4
  usual <- seq_aipe_smd(trained, control, width = 2.6)
  expect_false(attr(usual, "variance_unusable"))
  expect_identical(attr(usual, "notes"), character(0))
  expect_true(is.finite(row_value(usual, "lower_limit")))
})

test_that("seq_aipe_smd checks its arguments", {
  withr::local_seed(113)
  good <- rnorm(20)
  expect_error(seq_aipe_smd(group_2 = good, width = 0.3),
               "'group_1' must be supplied")
  expect_error(seq_aipe_smd(good, width = 0.3), "'group_2' must be supplied")
  expect_data_checks(function(v) seq_aipe_smd(v, good, width = 0.3),
                     "group_1", good)
  expect_data_checks(function(v) seq_aipe_smd(good, v, width = 0.3),
                     "group_2", good)
  expect_error(seq_aipe_smd(good, good[1:15], width = 0.30),
               "same length.*20 and 15")
  expect_scalar_checks(function(...) seq_aipe_smd(good, good + 1, ...))
})

test_that("the exported function agrees with the test search on whole paths", {
  withr::local_seed(113)
  take <- function(d, n) list(x = d$x[1:n], y = d$y[1:n])
  stat <- function(d) {
    c((mean(d$x) - mean(d$y)) / sqrt((var(d$x) + var(d$y)) / 2),
      .var_smd(d$x, d$y, floor = FALSE))
  }
  check <- function(d) seq_aipe_smd(d$x, d$y, width = 0.7)
  m <- seq_aipe_smd_pilot(0.7)
  for (i in 1:3) {
    d <- list(x = aipe_skewed(400) + 0.5, y = -aipe_skewed(400))
    run <- aipe_run(d, take, stat, check, width = 0.7, pilot = m)
    expect_equal(aipe_walk(d, take, check, pilot = m), run[["n"]])
  }
})

test_that("operating characteristics: Kelley, Darku, and Chattopadhyay (2018), Table 2", {
  # Row N(10, 1) versus N(9.5, 1), delta = 0.5, width 0.2, 95%: the article
  # reports a mean final sample size of 793.6 per group (theoretical 793)
  # and a mean width of 0.1999, from 5,000 replications. In 1,500
  # replications of this package the mean was 793.3 with a standard
  # deviation of 12.8 (and coverage 0.951, se 0.006, against the article's
  # 0.946). The mean of the 8 replications here had a standard deviation of
  # 4.85 over 200 other seeds (mean 793.4), so the window is 4 * 4.85 plus
  # the 0.3 between the article's mean and this package's: 19.7.
  withr::local_seed(113)
  w <- 0.2
  take <- function(d, n) list(x = d$x[1:n], y = d$y[1:n])
  stat <- function(d) {
    c((mean(d$x) - mean(d$y)) / sqrt((var(d$x) + var(d$y)) / 2),
      .var_smd(d$x, d$y, floor = FALSE))
  }
  check <- function(d) seq_aipe_smd(d$x, d$y, width = w)
  runs <- replicate(8, {
    d <- list(x = rnorm(1100, 10, 1), y = rnorm(1100, 9.5, 1))
    aipe_run(d, take, stat, check, width = w, pilot = seq_aipe_smd_pilot(w))
  })
  s <- aipe_summary(runs, truth = 0.5)
  expect_lt(abs(s[["mean_n"]] - 793.6), 19.7)
  # The width guarantee is exact, and at this sample size it is tight.
  expect_lte(s[["max_width"]], w)
  expect_gt(s[["mean_width"]], 0.199)
  # The variance estimate is never zero or negative here, so the procedure
  # as shipped and the rule as printed stop at the same point.
  expect_identical(s[["mean_n_printed"]], s[["mean_n"]])
})
