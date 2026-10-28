# Tests for seq_aipe_r() and seq_aipe_r_pilot(): the sequential AIPE
# procedure for the Pearson, Kendall, and Spearman correlations of Kelley,
# Darku, and Chattopadhyay (2019).

row_value <- function(out, term) out$value[out$term == term]
take_pairs <- function(d, n) list(x = d$x[1:n], y = d$y[1:n])

test_that("seq_aipe_r_pilot reproduces the pilots printed in the article", {
  # Equation 9 and the worked example on p. 495: a 95% interval of width
  # 0.10 gives m = max{4, ceiling(2 * 1.96 / 0.10)} = 40.
  expect_identical(seq_aipe_r_pilot(width = 0.10), 40L)
  # Simulation section, p. 496: width 0.1 at 90% gives m = 33.
  expect_identical(seq_aipe_r_pilot(width = 0.10, conf_level = 0.90), 33L)
  expect_identical(seq_aipe_r_pilot(width = 0.40), 10L)
  expect_identical(seq_aipe_r_pilot(width = 0.40, m_0 = 30), 30L)
  # A target close to the whole range: 2 * 1.96 / 1.9 = 2.06, so the
  # minimum of 4 is the pilot.
  expect_identical(seq_aipe_r_pilot(width = 1.9), 4L)
})

test_that("a target width as long as the range of a correlation is refused", {
  # A correlation lies between -1 and 1: an interval 2 or more wide says
  # nothing. Such a width (a percentage typed for a proportion, say) used
  # to give a pilot of 4 and a rule met almost at once.
  withr::local_seed(113)
  d <- aipe_pairs(20, rho = 0.5)
  for (w in c(2, 20)) {
    msg <- paste0("'width' must be a single finite number greater than 0 ",
                  "and less than 2; got ", w, ".")
    expect_error(seq_aipe_r_pilot(width = w), msg, fixed = TRUE)
    expect_error(seq_aipe_r(d$x, d$y, width = w), msg, fixed = TRUE)
    expect_error(seq_aipe_r(d$x, d$y, width = w, method = "kendall"), msg,
                 fixed = TRUE)
    expect_error(seq_aipe_r_sensitivity(0.4, width = w, G = 2, seed = 113),
                 msg, fixed = TRUE)
  }
  expect_s3_class(seq_aipe_r(d$x, d$y, width = 1.99), "smsd_check")
})

test_that("seq_aipe_r_pilot checks its arguments", {
  expect_pilot_checks(seq_aipe_r_pilot)
})

test_that("seq_aipe_r returns the documented rows and attributes", {
  withr::local_seed(113)
  d <- aipe_pairs(40, rho = 0.5)
  for (m in c("pearson", "kendall", "spearman")) {
    out <- seq_aipe_r(d$x, d$y, width = 0.50, method = m)
    expect_s3_class(out, c("smsd_check", "data.frame"), exact = TRUE)
    expect_identical(out$term, c("lower_limit", "r", "upper_limit",
                                 "xi_squared", "criterion_n", "current_n"))
    expect_true(all(is.finite(out$value)))
    expect_identical(attr(out, "cor_method"), m)
    # The label carries the proper name; cor_method the argument as given.
    expect_identical(attr(out, "method"),
                     c(pearson = "Pearson correlation",
                       kendall = "Kendall correlation",
                       spearman = "Spearman correlation")[[m]])
    expect_identical(attr(out, "variance_unusable"), FALSE)
    expect_identical(attr(out, "notes"), character(0))
    expect_identical(attr(out, "conf_level"), 0.95)
    expect_identical(attr(out, "width"), 0.50)
    expect_equal(attr(out, "current_n"), 40)
    expect_silent(seq_aipe_r(d$x, d$y, width = 0.50, method = m))
  }
  # The default method is Pearson, and partial matching is allowed.
  expect_identical(attr(seq_aipe_r(d$x, d$y, width = 0.5), "cor_method"),
                   "pearson")
  expect_identical(
    attr(seq_aipe_r(d$x, d$y, width = 0.5, method = "ken"), "cor_method"),
    "kendall")
})

test_that("seq_aipe_r computes the estimate, the interval (2019, Equation 11), and the criterion (Equation 10)", {
  withr::local_seed(113)
  d <- aipe_pairs(60, rho = 0.5)
  for (m in c("pearson", "kendall", "spearman")) {
    for (cl in c(0.90, 0.95)) {
      out <- seq_aipe_r(d$x, d$y, width = 0.50, conf_level = cl, method = m)
      r <- cor(d$x, d$y, method = m)
      xi2 <- .var_r(d$x, d$y, method = m)
      z <- qnorm(1 - (1 - cl) / 2)
      expect_equal(row_value(out, "r"), r)
      expect_equal(row_value(out, "xi_squared"), xi2)
      # Symmetric on the correlation scale: no Fisher z transformation.
      expect_equal(row_value(out, "lower_limit"), r - z * sqrt(xi2 / 60))
      expect_equal(row_value(out, "upper_limit"), r + z * sqrt(xi2 / 60))
      expect_equal(row_value(out, "criterion_n"),
                   ceiling((2 * z / 0.50)^2 * (xi2 + 1 / 60)))
      expect_identical(attr(out, "is_satisfied"),
                       60 >= row_value(out, "criterion_n"))
    }
  }
})

test_that("seq_aipe_r meets the rule exactly when n reaches the criterion", {
  withr::local_seed(113)
  d <- aipe_pairs(500, rho = 0.5)
  wide <- seq_aipe_r(d$x, d$y, width = 1.50)
  expect_true(attr(wide, "is_satisfied"))
  expect_lte(row_value(wide, "upper_limit") - row_value(wide, "lower_limit"),
             1.50)
  expect_false(attr(seq_aipe_r(d$x[1:20], d$y[1:20], width = 0.10),
                    "is_satisfied"))
})

test_that("seq_aipe_r is unchanged by shifting and rescaling either variable", {
  # A predictor recorded as a timestamp in seconds made the Pearson variance
  # estimate 52 instead of 0.33 with the raw power sum moment estimators.
  # Dyadic data (multiples of 1/1024) for Pearson, so that the shift itself
  # is exact; continuous data for the rank based coefficients, which depend
  # on the data only through the ranks.
  withr::local_seed(113)
  x <- rnorm(80)
  y <- 0.5 * x + rnorm(80)
  xd <- round(x * 1024) / 1024
  yd <- round(y * 1024) / 1024
  base <- seq_aipe_r(xd, yd, width = 0.3)
  moved <- seq_aipe_r(1.76e9 + 3600 * xd, 4 * yd - 1e6, width = 0.3)
  expect_equal(moved$value, base$value, tolerance = 1e-7)
  expect_identical(attr(moved, "is_satisfied"), attr(base, "is_satisfied"))
  for (m in c("kendall", "spearman")) {
    base <- seq_aipe_r(x, y, width = 0.3, method = m)
    moved <- seq_aipe_r(1.76e9 + 3600 * x, exp(y), width = 0.3, method = m)
    expect_equal(moved$value, base$value, tolerance = 1e-12)
  }
})

test_that("the Pearson check is defined when the sample covariance is exactly zero", {
  # Two 5 point items on ten respondents with a sample covariance of
  # exactly zero. Equation 70 as printed divides by the covariance (0 / 0);
  # with the factor r^2 multiplied through its value is
  # mu22_hat / (s_x^2 s_y^2) = 1/7 for these data (worked by hand in
  # test-var_estimators.R). Previously the check reported "variance not yet
  # usable" here.
  x <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
  y <- c(3, 1, 5, 1, 3, 3, 5, 1, 5, 3)
  expect_identical(sum((x - mean(x)) * (y - mean(y))), 0)
  out <- seq_aipe_r(x, y, width = 1.5)
  expect_false(attr(out, "variance_unusable"))
  expect_identical(attr(out, "notes"), character(0))
  expect_identical(row_value(out, "r"), 0)
  expect_equal(row_value(out, "xi_squared"), 1 / 7)
  z <- qnorm(0.975)
  expect_equal(row_value(out, "upper_limit"), z * sqrt(1 / 7 / 10))
  expect_equal(row_value(out, "lower_limit"), -z * sqrt(1 / 7 / 10))
  # (2 z / 1.5)^2 (1/7 + 1/10) = 1.66, so the criterion is 2.
  expect_equal(row_value(out, "criterion_n"), 2)
  expect_true(attr(out, "is_satisfied"))
})

test_that("seq_aipe_r warns once per call about ties for Kendall and Spearman only", {
  # Ratings on a 5 point scale.
  a <- c(1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 1, 2, 3, 4, 5, 3, 2, 4, 3, 3)
  b <- c(1, 1, 2, 3, 2, 4, 4, 5, 5, 4, 2, 2, 3, 5, 4, 3, 3, 4, 2, 3)
  for (m in c("kendall", "spearman")) {
    res <- aipe_warnings(seq_aipe_r(a, b, width = 0.5, method = m))
    expect_length(res$warnings, 1L)
    expect_match(res$warnings, "tied values")
    expect_match(res$warnings, m)
    # The estimate is still the usual tie corrected coefficient.
    expect_equal(row_value(res$value, "r"), cor(a, b, method = m))
    # For Kendall's tau that is tau-b, and the warning says so; the
    # Spearman warning has no such sentence.
    expect_identical(grepl("Kendall's tau-b", res$warnings, fixed = TRUE),
                     m == "kendall")
    # Ties in only one of the two variables are enough.
    cont <- a + seq_along(a) / 1000
    expect_warning(seq_aipe_r(cont, b, width = 0.5, method = m),
                   "tied values")
    expect_silent(seq_aipe_r(cont, b + seq_along(b) / 977, width = 0.5,
                             method = m))
  }
  # The Pearson estimator makes no continuity assumption.
  expect_silent(seq_aipe_r(a, b, width = 0.5, method = "pearson"))

  # With ties, tau-b is not the statistic of the article's Equation 23,
  # 2 / (n (n - 1)) times the sum over pairs of the products of signs
  # (tau-a), to which the variance estimator belongs. Coded directly:
  n <- length(a)
  pairs <- utils::combn(n, 2)
  tau_a <- 2 * sum(sign(a[pairs[1, ]] - a[pairs[2, ]]) *
                     sign(b[pairs[1, ]] - b[pairs[2, ]])) / (n * (n - 1))
  tau_b <- suppressWarnings(
    row_value(seq_aipe_r(a, b, width = 0.5, method = "kendall"), "r")
  )
  expect_gt(tau_b - tau_a, 0.05)
  # Without ties the two are the same number.
  cont_a <- a + seq_along(a) / 1000
  cont_b <- b + seq_along(b) / 977
  tau_a <- 2 * sum(sign(cont_a[pairs[1, ]] - cont_a[pairs[2, ]]) *
                     sign(cont_b[pairs[1, ]] - cont_b[pairs[2, ]])) /
    (n * (n - 1))
  expect_equal(row_value(seq_aipe_r(cont_a, cont_b, width = 0.5,
                                    method = "kendall"), "r"), tau_a)
})

test_that("a perfect sample correlation is not yet usable, for all three methods", {
  # A sample correlation of 1 or -1 says nothing about its own sampling
  # variance. The Pearson estimator is zero for collinear pairs and the
  # Kendall estimator for perfectly discordant ones, but for perfectly
  # concordant pairs the Kendall estimator is positive, and so is the
  # Spearman estimator for perfectly monotone pairs in either direction:
  # the rule as printed would be met at once, with a limit outside
  # [-1, 1].
  up <- (1:12)^2
  down <- 12:1
  z <- stats::qnorm(0.975)
  raw_spearman <- .var_r(1:12, down, method = "spearman", floor = FALSE)
  raw_kendall <- .var_r(1:12, up, method = "kendall", floor = FALSE)
  expect_gt(raw_spearman, 12^(-3))
  expect_gt(raw_kendall, 12^(-3))
  expect_true(aipe_rule_met(12, raw_spearman, 0.4))
  expect_true(aipe_rule_met(12, raw_kendall, 0.4))
  # The limits the printed rule would report: below -1, and above 1.
  expect_lt(-1 - z * sqrt(raw_spearman / 12), -1.03)
  expect_gt(1 + z * sqrt(raw_kendall / 12), 1.10)

  cases <- list(
    list(y = 3 * (1:12) + 1, method = "pearson", r = 1),
    list(y = 1 - 0.3 * (1:12), method = "pearson", r = -1),
    list(y = up, method = "kendall", r = 1),
    list(y = down, method = "kendall", r = -1),
    list(y = up, method = "spearman", r = 1),
    list(y = down, method = "spearman", r = -1)
  )
  for (case in cases) {
    out <- seq_aipe_r(1:12, case$y, width = 0.4, method = case$method)
    label <- paste(case$method, case$r)
    expect_equal(row_value(out, "r"), case$r, label = label)
    expect_false(attr(out, "is_satisfied"), label = label)
    expect_true(attr(out, "variance_unusable"), label = label)
    expect_true(is.na(row_value(out, "lower_limit")), label = label)
    expect_true(is.na(row_value(out, "upper_limit")), label = label)
    expect_identical(row_value(out, "xi_squared"), 12^(-3), label = label)
    expect_identical(attr(out, "criterion_n"), 13, label = label)
    expect_identical(attr(out, "notes"), .note_variance_unusable(),
                     label = label)
  }
  # One pair out of order and the correlation is no longer perfect: the
  # check is an ordinary one again, for a sample as large as 2,000 pairs,
  # where Spearman's rho is within 2e-9 of 1.
  for (n in c(12, 2000)) {
    x <- seq_len(n)
    y <- x
    y[1:2] <- y[2:1]
    for (m in c("kendall", "spearman")) {
      out <- seq_aipe_r(x, y, width = 0.4, method = m)
      expect_lt(row_value(out, "r"), 1)
      expect_false(attr(out, "variance_unusable"), label = paste(m, n))
      expect_true(is.finite(row_value(out, "lower_limit")))
    }
  }
})

test_that("seq_aipe_r checks its arguments", {
  withr::local_seed(113)
  good <- rnorm(20)
  other <- rnorm(20)
  # The data arguments are x and y (paired scores), as in
  # seq_aipe_reg_coef(); they were group_1 and group_2.
  expect_identical(names(formals(seq_aipe_r)),
                   c("x", "y", "width", "conf_level", "method"))
  expect_identical(seq_aipe_r(x = good, y = other, width = 0.3)$value,
                   seq_aipe_r(good, other, width = 0.3)$value)
  expect_error(seq_aipe_r(group_1 = good, group_2 = other, width = 0.3),
               "unused argument")
  expect_error(seq_aipe_r(y = good, width = 0.3), "'x' must be supplied")
  expect_error(seq_aipe_r(good, width = 0.3), "'y' must be supplied")
  expect_data_checks(function(v) seq_aipe_r(v, other, width = 0.3),
                     "x", good, short = FALSE)
  expect_data_checks(function(v) seq_aipe_r(other, v, width = 0.3),
                     "y", good, short = FALSE)
  expect_error(seq_aipe_r(good, other[1:15], width = 0.30),
               "'x' and 'y' must have the same length.*20 and 15")
  # The whole message: a sequential check closes with what to do next (a
  # plain estimator such as slr_slope(), which shares the pair check, does
  # not).
  expect_error(seq_aipe_r(good[1:3], other[1:3], width = 0.5),
               paste0("^'x' needs at least 4 observations; it has 3\\. ",
                      "Collect more data, then check the stopping rule ",
                      "again\\.$"))
  expect_error(seq_aipe_r(good, replace(other, 4, NA), width = 0.3),
               "'y'.*missing values")
  expect_error(seq_aipe_r(good, matrix(other, ncol = 2), width = 0.3,
                          method = "kendall"),
               "'y'.*not a matrix")
  expect_scalar_checks(function(...) seq_aipe_r(good, other, ...))
  # The shared choice check names the argument, lists the choices, and
  # shows what was received.
  expect_error(seq_aipe_r(good, other, width = 0.3, method = "foo"),
               paste0("'method' must be one of \"pearson\", \"kendall\", ",
                      "\"spearman\"; got \"foo\""), fixed = TRUE)
  expect_error(seq_aipe_r(good, other, width = 0.3, method = 2),
               "'method' must be one of.*got 2")
  expect_error(seq_aipe_r(good, other, width = 0.3,
                          method = c("kendall", "spearman")),
               "'method' must be one of.*vector of length 2")
  expect_error(seq_aipe_r(good, other, width = 0.3, method = NA_character_),
               "'method' must be one of")
})

test_that("the exported function agrees with the test search on whole paths", {
  withr::local_seed(113)
  for (m in c("pearson", "kendall", "spearman")) {
    stat <- function(d) {
      c(cor(d$x, d$y, method = m), .var_r(d$x, d$y, m, floor = FALSE))
    }
    check <- function(d) seq_aipe_r(d$x, d$y, width = 0.5, method = m)
    pilot <- seq_aipe_r_pilot(0.5)
    for (i in 1:2) {
      d <- aipe_pairs(400, rho = 0.5)
      run <- aipe_run(d, take_pairs, stat, check, width = 0.5, pilot = pilot)
      expect_equal(aipe_walk(d, take_pairs, check, pilot = pilot), run[["n"]])
    }
  }
})

test_that("operating characteristics, Pearson: Kelley, Darku, and Chattopadhyay (2019), Table 2", {
  # Row N2(0, 0, 1, 1, .5), width 0.2, 95%: mean final sample size 204.0
  # (se 0.7435 from 5,000 replications, so a standard deviation of 52.6),
  # theoretical 217, coverage 0.9012, mean width 0.1938. The article's rule
  # stops even when the raw variance estimate is zero or negative (it is
  # replaced by the bound n^-3), so the published mean is compared with the
  # stopping time of the rule as printed: in 5,000 replications with this
  # package's estimator that rule gave 205.4 (standard deviation 52.7),
  # coverage 0.905, width 0.1942. The procedure as shipped does not stop on
  # such an estimate (1.5% of pilots here) and gave 207.1, coverage 0.913.
  # The mean of the 20 replications here had a standard deviation of 11.9
  # over 200 other seeds (52.7 / sqrt(20) is 11.8), and this package's mean
  # under the printed rule was within 1.4 of the article's in three runs
  # (205.4, 203.4, 202.9), so the window is 4 * 11.9 + 1.4 = 49.0.
  withr::local_seed(113)
  w <- 0.2
  stat <- function(d) {
    c(cor(d$x, d$y), .var_r(d$x, d$y, "pearson", floor = FALSE))
  }
  check <- function(d) seq_aipe_r(d$x, d$y, width = w)
  runs <- replicate(20, {
    aipe_run(aipe_pairs(900, rho = 0.5), take_pairs, stat, check, width = w,
             pilot = seq_aipe_r_pilot(w))
  })
  s <- aipe_summary(runs, truth = 0.5)
  expect_lt(abs(s[["mean_n_printed"]] - 204.0), 49.0)
  # The shipped procedure never stops before the printed rule would.
  expect_true(all(runs["n", ] >= runs["n_printed", ]))
  expect_lte(s[["max_width"]], w)
})

test_that("operating characteristics, Spearman: Kelley, Darku, and Chattopadhyay (2019), Table 10", {
  # Row N2(0, 0, 1, 1, .5), width 0.2, 95%: mean final sample size 238.8
  # (se 0.5289 from 5,000 replications, so a standard deviation of 37.4),
  # theoretical 243, coverage 0.9300, mean width 0.1986. In 900
  # replications of this package: 238.3 (standard deviation 34.2), coverage
  # 0.941 (se 0.008), width 0.1987; in 1,000 more: 237.5 (standard deviation
  # 33.4). The final sample size has a long right tail, so its standard
  # deviation is itself estimated loosely: the mean of the 8 replications
  # here had a standard deviation of 14.1 over 200 other seeds (the
  # article's 37.4 gives 37.4 / sqrt(8) = 13.2), and the window uses the
  # larger: 4 * 14.1 + 1.3 = 57.7.
  withr::local_seed(113)
  w <- 0.2
  stat <- function(d) c(cor(d$x, d$y, method = "spearman"),
                        .var_r(d$x, d$y, "spearman", floor = FALSE))
  check <- function(d) seq_aipe_r(d$x, d$y, width = w, method = "spearman")
  runs <- replicate(8, {
    aipe_run(aipe_pairs(700, rho = 0.5), take_pairs, stat, check, width = w,
             pilot = seq_aipe_r_pilot(w))
  })
  s <- aipe_summary(runs, truth = 6 / pi * asin(0.25))
  expect_lt(abs(s[["mean_n"]] - 238.8), 57.7)
  expect_lte(s[["max_width"]], w)
})

test_that("operating characteristics, Kendall: Kelley, Darku, and Chattopadhyay (2019), Table 8", {
  # Row N2(0, 0, 1, 1, .5), width 0.2, 95%: mean final sample size 133.2
  # (se 0.1514 from 5,000 replications, so a standard deviation of 10.7),
  # theoretical 131, coverage 0.9318, mean width 0.1972. In 10,000
  # replications of the package's validation study
  # (data-raw/validation_study.R): 135.9 (se 0.20, standard deviation
  # 20.4), coverage 0.937 (se 0.002), width 0.1968. The package's mean is
  # 2% above the article's and its standard deviation nearly twice as
  # large, although the estimator equals the article's Equations 27 to 29
  # exactly (see test-var_estimators.R); the difference is unresolved (the
  # article's simulation code is not available) and is recorded here rather
  # than hidden. With the 10 replications here the standard error of the
  # mean is 20.4 / sqrt(10) = 6.5, so the window of 31.3 around the
  # article's mean is the 2.7 by which this package's mean exceeds it plus
  # 4.4 standard errors.
  withr::local_seed(113)
  w <- 0.2
  stat <- function(d) c(cor(d$x, d$y, method = "kendall"),
                        .var_r(d$x, d$y, "kendall", floor = FALSE))
  check <- function(d) seq_aipe_r(d$x, d$y, width = w, method = "kendall")
  runs <- replicate(10, {
    aipe_run(aipe_pairs(500, rho = 0.5), take_pairs, stat, check, width = w,
             pilot = seq_aipe_r_pilot(w))
  })
  s <- aipe_summary(runs, truth = 2 / pi * asin(0.5))
  expect_lt(abs(s[["mean_n"]] - 133.2), 31.3)
  expect_lte(s[["max_width"]], w)
})
