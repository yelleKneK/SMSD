# Tests for seq_aipe_sm() and seq_aipe_sm_pilot(): the sequential AIPE
# procedure of Kelley, Darku, and Chattopadhyay (2018) applied to the
# standardized mean (their Example 3).

row_value <- function(out, term) out$value[out$term == term]

test_that("seq_aipe_sm_pilot is max{m_0, ceiling(2 z / width)} (2018, Equation 21)", {
  expect_identical(seq_aipe_sm_pilot(width = 0.20), 20L)
  expect_identical(seq_aipe_sm_pilot(width = 0.50), 8L)
  expect_identical(seq_aipe_sm_pilot(width = 0.50, conf_level = 0.90), 7L)
  expect_identical(seq_aipe_sm_pilot(width = 0.50, m_0 = 15), 15L)
})

test_that("seq_aipe_sm_pilot checks its arguments", {
  expect_pilot_checks(seq_aipe_sm_pilot)
})

test_that("seq_aipe_sm returns the documented rows and attributes", {
  withr::local_seed(113)
  x <- rnorm(40, mean = 1, sd = 2)
  out <- seq_aipe_sm(x, width = 0.50)
  expect_s3_class(out, c("smsd_check", "data.frame"), exact = TRUE)
  expect_identical(out$term, c("lower_limit", "sm", "upper_limit",
                               "xi_squared", "criterion_n", "current_n"))
  expect_true(all(is.finite(out$value)))
  expect_identical(attr(out, "method"), "standardized mean")
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_identical(attr(out, "notes"), character(0))
  expect_identical(attr(out, "conf_level"), 0.95)
  expect_identical(attr(out, "width"), 0.50)
  expect_equal(attr(out, "current_n"), 40)
  expect_silent(seq_aipe_sm(x, width = 0.50))
})

test_that("seq_aipe_sm computes the estimate, the interval (2018, Equation 23), and the criterion (Equation 22)", {
  withr::local_seed(113)
  x <- rnorm(40, mean = 1, sd = 2)
  for (cl in c(0.90, 0.95)) {
    out <- seq_aipe_sm(x, width = 0.50, conf_level = cl)
    est <- mean(x) / sd(x)
    xi2 <- .var_sm(x)
    z <- qnorm(1 - (1 - cl) / 2)
    expect_equal(row_value(out, "sm"), est)
    expect_equal(row_value(out, "xi_squared"), xi2)
    expect_equal(row_value(out, "lower_limit"), est - z * sqrt(xi2 / 40))
    expect_equal(row_value(out, "upper_limit"), est + z * sqrt(xi2 / 40))
    expect_equal(row_value(out, "criterion_n"),
                 ceiling((2 * z / 0.50)^2 * (xi2 + 1 / 40)))
    expect_identical(attr(out, "is_satisfied"),
                     40 >= row_value(out, "criterion_n"))
  }
})

test_that("seq_aipe_sm handles a negative mean and is scale invariant", {
  withr::local_seed(113)
  x <- rnorm(60, mean = 1, sd = 2)
  pos <- seq_aipe_sm(x, width = 0.5)
  neg <- seq_aipe_sm(-x, width = 0.5)
  # Reversing the sign of the data reverses the estimate and the interval;
  # the variance is unchanged (the third moment changes sign with the mean).
  expect_equal(row_value(neg, "sm"), -row_value(pos, "sm"))
  expect_equal(row_value(neg, "xi_squared"), row_value(pos, "xi_squared"))
  expect_equal(row_value(neg, "lower_limit"), -row_value(pos, "upper_limit"))
  expect_equal(seq_aipe_sm(250 * x, width = 0.5)$value, pos$value,
               tolerance = 1e-12)
})

test_that("seq_aipe_sm meets the rule exactly when n reaches the criterion", {
  withr::local_seed(113)
  x <- rnorm(500, mean = 1, sd = 2)
  wide <- seq_aipe_sm(x, width = 1.50)
  expect_true(attr(wide, "is_satisfied"))
  expect_lte(row_value(wide, "upper_limit") - row_value(wide, "lower_limit"),
             1.50)
  expect_false(attr(seq_aipe_sm(x[1:20], width = 0.05), "is_satisfied"))
})

test_that("seq_aipe_sm checks its arguments", {
  withr::local_seed(113)
  good <- rnorm(20, 1)
  expect_error(seq_aipe_sm(width = 0.3), "'x' must be supplied")
  expect_data_checks(function(v) seq_aipe_sm(v, width = 0.3), "x", good)
  expect_scalar_checks(function(...) seq_aipe_sm(good, ...))
})

test_that("operating characteristics: the standardized mean procedure reaches its target", {
  # No published table exists for the standardized mean, so the targets are
  # the theoretical sample size and the nominal coverage. Normal data with
  # mu / sigma = 0.4: xi^2 = 1 + 0.4^2 / 2 = 1.08, so at width 0.5 the
  # theoretical sample size is ceiling((2 * 1.96 / 0.5)^2 * 1.08) = 67. In
  # 3,000 replications of this condition the mean final sample size was
  # 67.6 (standard deviation 8.0) and coverage 0.946 (se 0.004). The mean
  # of the 80 replications here had a standard deviation of 0.99 over 200
  # other seeds, so its window is 4 * 0.99 + 0.6 = 4.6. The standard error
  # of the coverage is sqrt(0.946 * 0.054 / 80) = 0.0253, so its window is
  # 4 * 0.0253 + 0.004 = 0.106.
  z <- qnorm(0.975)
  expect_equal(ceiling((2 * z / 0.5)^2 * 1.08), 67)
  withr::local_seed(113)
  w <- 0.5
  take <- function(d, n) d[1:n]
  stat <- function(d) c(mean(d) / sd(d), .var_sm(d, floor = FALSE))
  check <- function(d) seq_aipe_sm(d, width = w)
  runs <- replicate(80, {
    aipe_run(rnorm(400, mean = 0.4), take, stat, check, width = w,
             pilot = seq_aipe_sm_pilot(w))
  })
  s <- aipe_summary(runs, truth = 0.4)
  expect_lt(abs(s[["mean_n"]] - 67), 4.6)
  expect_lt(abs(s[["coverage"]] - 0.95), 0.106)
  expect_lte(s[["max_width"]], w)
})
