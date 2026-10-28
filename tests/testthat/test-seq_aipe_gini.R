# Tests for seq_aipe_gini() and seq_aipe_gini_pilot(): the sequential fixed
# width interval for the Gini index of Chattopadhyay and De (2016), with the
# full width (twice their half width d) as the target.

row_value <- function(out, term) out$value[out$term == term]

test_that("seq_aipe_gini_pilot is the pilot of Chattopadhyay and De (2016) in terms of the full width", {
  # Their pilot is max{4, ceiling(z / d)} with d the half width. For the
  # settings of their simulation study (d = 0.01 and 0.02; alpha = 0.10 and
  # 0.05): z / d = 164.5, 196.0 (just under), 82.2, 98.0 (just under).
  expect_identical(seq_aipe_gini_pilot(width = 0.02, conf_level = 0.90), 165L)
  expect_identical(seq_aipe_gini_pilot(width = 0.02, conf_level = 0.95), 196L)
  expect_identical(seq_aipe_gini_pilot(width = 0.04, conf_level = 0.90), 83L)
  expect_identical(seq_aipe_gini_pilot(width = 0.04, conf_level = 0.95), 98L)
  # A target close to the whole range: 2 * 1.96 / 0.99 = 3.96, so the
  # minimum of 4 is the pilot.
  expect_identical(seq_aipe_gini_pilot(width = 0.99), 4L)
  expect_identical(seq_aipe_gini_pilot(width = 0.2, m_0 = 40), 40L)
})

test_that("seq_aipe_gini_pilot checks its arguments", {
  expect_pilot_checks(seq_aipe_gini_pilot)
})

test_that("seq_aipe_gini returns the documented rows and attributes", {
  withr::local_seed(113)
  x <- rexp(60, rate = 0.5)
  out <- seq_aipe_gini(x, width = 0.30)
  expect_s3_class(out, c("smsd_check", "data.frame"), exact = TRUE)
  expect_identical(out$term, c("lower_limit", "gini", "upper_limit",
                               "xi_squared", "criterion_n", "current_n"))
  expect_true(all(is.finite(out$value)))
  expect_identical(attr(out, "method"), "Gini index")
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_identical(attr(out, "notes"), character(0))
  expect_identical(attr(out, "conf_level"), 0.95)
  expect_identical(attr(out, "width"), 0.30)
  expect_equal(attr(out, "current_n"), 60)
  expect_silent(seq_aipe_gini(x, width = 0.30))
})

test_that("seq_aipe_gini computes the estimate, the interval, and the stopping rule of Chattopadhyay and De (2016, Equation 12)", {
  withr::local_seed(113)
  x <- rlnorm(60, meanlog = 4, sdlog = 0.6)
  n <- 60
  # The Gini index as the mean absolute difference over twice the mean.
  g <- sum(abs(outer(x, x, "-"))) / (n * (n - 1)) / (2 * mean(x))
  for (cl in c(0.90, 0.95)) {
    d <- 0.06                                  # their half width
    out <- seq_aipe_gini(x, width = 2 * d, conf_level = cl)
    z <- qnorm(1 - (1 - cl) / 2)
    v2 <- row_value(out, "xi_squared")
    expect_equal(row_value(out, "gini"), g)
    expect_equal(v2, .var_gini(x))
    expect_equal(row_value(out, "lower_limit"), g - z * sqrt(v2 / n))
    expect_equal(row_value(out, "upper_limit"), g + z * sqrt(v2 / n))
    # Equation 12 in their parameterization: n >= (z / d)^2 (V_n^2 + 1/n).
    expect_equal(row_value(out, "criterion_n"),
                 ceiling((z / d)^2 * (v2 + 1 / n)))
    expect_identical(attr(out, "is_satisfied"),
                     n >= row_value(out, "criterion_n"))
  }
})

test_that("seq_aipe_gini handles integer data without overflow", {
  # 6,000 incomes stored as integers: the sums involved exceed the integer
  # range. Previously gini_index() returned NA with a warning and the check
  # reported the rule as met with a criterion of 2.
  withr::local_seed(113)
  income <- as.integer(round(rlnorm(6000, 11, 0.8)))
  expect_silent(out <- seq_aipe_gini(income, width = 0.05))
  expect_equal(out$value, seq_aipe_gini(as.numeric(income), width = 0.05)$value)
  expect_true(is.finite(row_value(out, "gini")))
  expect_gt(row_value(out, "criterion_n"), 500)
})

test_that("seq_aipe_gini does not depend on the unit of measurement", {
  withr::local_seed(113)
  x <- rlnorm(80, 4, 0.6)
  expect_equal(seq_aipe_gini(1000 * x, width = 0.2)$value,
               seq_aipe_gini(x, width = 0.2)$value, tolerance = 1e-12)
  expect_equal(seq_aipe_gini(1e-9 * x, width = 0.2)$value,
               seq_aipe_gini(x, width = 0.2)$value, tolerance = 1e-12)
})

test_that("seq_aipe_gini checks its arguments", {
  withr::local_seed(113)
  good <- rexp(20)
  expect_error(seq_aipe_gini(width = 0.3), "'x' must be supplied")
  expect_data_checks(function(v) seq_aipe_gini(v, width = 0.3), "x", good)
  # The Gini index is defined for nonnegative data: previously a sample
  # with negative values gave an "index" of 1.31 without complaint.
  expect_error(seq_aipe_gini(c(good, -0.5), width = 0.3),
               "'x' must contain no negative values.*1 of 21")
  # Every value negative: refused by the check on the data, before the
  # index is looked at (gini_index() refuses a mix of signs in the same
  # words, so only a sample with no positive value shows which check
  # stopped it).
  expect_error(seq_aipe_gini(c(-4, -3, -2, -1), width = 0.2),
               "'x' must contain no negative values.*4 of 4 are negative")
  # Zeros are allowed (households with no income).
  expect_silent(seq_aipe_gini(c(good, 0, 0), width = 0.3))
  expect_scalar_checks(function(...) seq_aipe_gini(good, ...))
})

test_that("a target width as long as the range of the Gini index is refused", {
  # The Gini index lies between 0 and 1: an interval 1 or more wide says
  # nothing. A width of 1.5 used to give a pilot of 4 and, for four
  # observations, a rule that was met with an upper limit above 1.
  withr::local_seed(113)
  x <- rexp(20)
  for (w in c(1, 1.5)) {
    msg <- paste0("'width' must be a single finite number greater than 0 ",
                  "and less than 1; got ", w, ".")
    expect_error(seq_aipe_gini_pilot(width = w), msg, fixed = TRUE)
    expect_error(seq_aipe_gini(x, width = w), msg, fixed = TRUE)
    expect_error(seq_aipe_gini_sensitivity(0.3, width = w, G = 2, seed = 113),
                 msg, fixed = TRUE)
  }
  expect_s3_class(seq_aipe_gini(x, width = 0.99), "smsd_check")
})

test_that("operating characteristics: the Gini index procedure reaches its target", {
  # Chattopadhyay and De (2016) report simulations, but their table is not
  # on hand, so the targets are theoretical. Exponential data: G = 1/2 and
  # xi^2 = 1/12, so at full width 0.1 (their d = 0.05) the stopping rule
  # n >= (z / d)^2 (xi^2 + 1 / n) is first met, with xi^2 known, at n = 140
  # (the theoretical size without the 1 / n term is 129). In 3,000
  # replications of this condition the mean final sample size was 141.9
  # (standard deviation 23.1) and coverage 0.951 (se 0.004); in 8,000 more
  # (200 other seeds of this test) they were 141.0 and 0.944. The standard
  # error of the mean of the 40 replications here is 23.1 / sqrt(40) = 3.65,
  # so its window is 4 * 3.65 + 1.9 = 16.5. The standard error of the
  # coverage is sqrt(0.944 * 0.056 / 40) = 0.0364, so its window is
  # 4 * 0.0364 + 0.006 = 0.152.
  z <- qnorm(0.975)
  n <- 100:200
  expect_equal(min(n[n >= (z / 0.05)^2 * (1 / 12 + 1 / n)]), 140)
  expect_equal(ceiling((z / 0.05)^2 / 12), 129)
  withr::local_seed(113)
  w <- 0.1
  take <- function(d, n) d[1:n]
  stat <- function(d) c(gini_index(d), .var_gini(d, floor = FALSE))
  check <- function(d) seq_aipe_gini(d, width = w)
  runs <- replicate(40, {
    aipe_run(rexp(700), take, stat, check, width = w,
             pilot = seq_aipe_gini_pilot(w))
  })
  s <- aipe_summary(runs, truth = 0.5)
  expect_lt(abs(s[["mean_n"]] - 140), 16.5)
  expect_lt(abs(s[["coverage"]] - 0.95), 0.152)
  expect_lte(s[["max_width"]], w)
})
