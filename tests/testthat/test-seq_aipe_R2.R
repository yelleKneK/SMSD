# Tests for seq_aipe_R2_pilot(), seq_aipe_R2(), and seq_aipe_R2_sensitivity():
# sequential AIPE for the squared multiple correlation coefficient, Kelley,
# Darku, and Chattopadhyay (2019, pp. 503--507 and Appendix C). The article
# writes k for the number of predictors; the package writes p.
#
# Hard coded expected values were computed by hand from the article's
# equations (base R arithmetic, without the package), not taken from the
# functions under test.

z90 <- 1.644853627   # qnorm(0.95)
z95 <- 1.959963985   # qnorm(0.975)

# Data whose sample squared multiple correlation equals `R2` exactly: the
# response is a combination of one direction inside the span of the centered
# predictors and one direction orthogonal to them and to the constant.
exact_R2_data <- function(n, p, R2, shift = 0, scale = 1) {
  withr::local_seed(113)
  x <- matrix(stats::rnorm(n * p), n, p)
  q <- qr.Q(qr(cbind(1, x, stats::rnorm(n))))
  y <- sqrt(R2) * q[, 2] + sqrt(1 - R2) * q[, p + 2]
  list(x = x, y = shift + scale * y)
}

row_value <- function(out, term) out$value[out$term == term]

survey_data <- function(n, p = 3, seed = 113) {
  withr::local_seed(seed)
  x <- matrix(stats::rnorm(n * p), n, p)
  y <- drop(x %*% rep(0.4, p)) + stats::rnorm(n)
  list(x = x, y = y)
}

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

test_that("the arguments are in the order the package uses everywhere", {
  # Pilot: the method's own arguments, then m_0. Check: data, the precision
  # target, conf_level, then what is specific to the method. Sensitivity:
  # the design, then G, seed, m_0, max_n, batch.
  expect_identical(names(formals(seq_aipe_R2_pilot)),
                   c("width", "p", "conf_level", "ci_method", "m_0"))
  expect_identical(names(formals(seq_aipe_R2)),
                   c("x", "y", "width", "conf_level", "ci_method"))
  expect_identical(names(formals(seq_aipe_R2_sensitivity)),
                   c("true_value", "width", "p", "conf_level", "ci_method",
                     "G", "seed", "m_0", "max_n", "batch"))
  for (f in list(seq_aipe_R2_pilot, seq_aipe_R2, seq_aipe_R2_sensitivity)) {
    expect_identical(formals(f)$conf_level, 0.95)
    expect_identical(eval(formals(f)$ci_method),
                     c("wald", "log_transformed"))
  }
  expect_identical(formals(seq_aipe_R2_sensitivity)$G, 1000)
  expect_null(formals(seq_aipe_R2_sensitivity)$seed)
})

# ---------------------------------------------------------------------------
# Pilot sample size
# ---------------------------------------------------------------------------

test_that("the Wald pilot is max{p, 4 z / omega} rounded up (2019, p. 505 and Equation 97)", {
  # The article prints max{k, 4 z / omega} with no rounding; a sample size is
  # a whole number, and rounding up (not down) is what reproduces the P2 = .90
  # rows of Tables 13 and 14, whose mean final sample sizes are just above
  # the pilot (68.66 and 84.17 at omega = .10).
  # 4 z / omega by hand: 65.79 (90%, .10), 78.40 (95%, .10), 131.59 (90%,
  # .05), 156.80 (95%, .05).
  expect_identical(seq_aipe_R2_pilot(0.10, p = 2, conf_level = 0.90), 66L)
  expect_identical(seq_aipe_R2_pilot(0.10, p = 2, conf_level = 0.95), 79L)
  expect_identical(seq_aipe_R2_pilot(0.05, p = 10, conf_level = 0.90), 132L)
  expect_identical(seq_aipe_R2_pilot(0.05, p = 10, conf_level = 0.95), 157L)
  # p governs when it exceeds 4 z / omega (8.71 at omega = .9), but never
  # below p + 2, where the sample R2 is identically one.
  expect_identical(seq_aipe_R2_pilot(0.90, p = 30), 32L)
  expect_identical(seq_aipe_R2_pilot(0.90, p = 8), 10L)
})

test_that("the log transformed pilot is p + 2 (2019, Equation 95)", {
  for (p in c(1L, 2L, 5L, 10L)) {
    expect_identical(
      seq_aipe_R2_pilot(0.10, p = p, ci_method = "log_transformed"), p + 2L)
  }
  # It does not depend on the width or the confidence level.
  expect_identical(
    seq_aipe_R2_pilot(0.02, p = 4, conf_level = 0.99,
                      ci_method = "log_transformed"), 6L)
})

test_that("m_0 raises the pilot and cannot be below p + 2", {
  expect_identical(seq_aipe_R2_pilot(0.25, p = 3), 32L)   # 4 z / omega = 31.36
  expect_identical(seq_aipe_R2_pilot(0.25, p = 3, m_0 = 50), 50L)
  expect_identical(seq_aipe_R2_pilot(0.25, p = 3, m_0 = 20), 32L)
  expect_identical(
    seq_aipe_R2_pilot(0.25, p = 3, ci_method = "log_transformed", m_0 = 30),
    30L)
  expect_error(seq_aipe_R2_pilot(0.25, p = 3, m_0 = 4), "'m_0'.*at least 5")
  expect_error(seq_aipe_R2_pilot(0.25, p = 3, m_0 = 10.5), "'m_0'")
  expect_error(seq_aipe_R2_pilot(0.25, p = 3, m_0 = NA), "'m_0'")
})

test_that("seq_aipe_R2_pilot checks its arguments", {
  expect_error(seq_aipe_R2_pilot(p = 3), "'width' must be supplied")
  expect_error(seq_aipe_R2_pilot(0, p = 3), "'width'")
  expect_error(seq_aipe_R2_pilot(-0.1, p = 3), "'width'")
  expect_error(seq_aipe_R2_pilot(1, p = 3), "'width'.*less than 1")
  expect_error(seq_aipe_R2_pilot(NA_real_, p = 3), "'width'")
  expect_error(seq_aipe_R2_pilot(c(0.1, 0.2), p = 3), "'width'")
  expect_error(seq_aipe_R2_pilot("0.1", p = 3), "'width'")
  expect_error(seq_aipe_R2_pilot(0.2), "'p' must be supplied")
  expect_error(seq_aipe_R2_pilot(0.2, p = 0), "'p'")
  expect_error(seq_aipe_R2_pilot(0.2, p = 2.5), "'p'")
  expect_error(seq_aipe_R2_pilot(0.2, p = NA), "'p'")
  # A confidence level passed by position where p belongs is caught.
  expect_error(seq_aipe_R2_pilot(0.2, 0.95), "'p'")
  expect_error(seq_aipe_R2_pilot(0.2, p = 3, conf_level = 1), "'conf_level'")
  expect_error(seq_aipe_R2_pilot(0.2, p = 3, conf_level = 0), "'conf_level'")
  expect_error(seq_aipe_R2_pilot(0.2, p = 3, conf_level = 95), "'conf_level'")
  expect_error(seq_aipe_R2_pilot(0.2, p = 3, ci_method = "exact"),
               paste0("'ci_method' must be one of \"wald\", ",
                      "\"log_transformed\"; got \"exact\""))
  expect_error(seq_aipe_R2_pilot(0.2, p = 3, ci_method = "Wald"),
               "'ci_method' must be one of")
  expect_error(seq_aipe_R2_pilot(0.2, p = 3, ci_method = 1),
               "'ci_method' must be one of.*got 1")
  expect_error(seq_aipe_R2_pilot(0.2, p = 3,
                                 ci_method = c("log_transformed", "wald")),
               "'ci_method' must be one of.*character vector of length 2")
  expect_error(seq_aipe_R2_pilot(0.2, p = 3, ci_method = NA_character_),
               "'ci_method' must be one of")
  expect_error(seq_aipe_R2_pilot(0.2, p = 3, ci_method = NULL),
               "'ci_method' must be one of.*NULL")
})

test_that("ci_method may be abbreviated, and the default is the Wald interval", {
  expect_identical(seq_aipe_R2_pilot(0.2, p = 3, ci_method = "log"),
                   seq_aipe_R2_pilot(0.2, p = 3,
                                     ci_method = "log_transformed"))
  expect_identical(seq_aipe_R2_pilot(0.2, p = 3, ci_method = "w"),
                   seq_aipe_R2_pilot(0.2, p = 3, ci_method = "wald"))
  expect_identical(seq_aipe_R2_pilot(0.2, p = 3),
                   seq_aipe_R2_pilot(0.2, p = 3, ci_method = "wald"))
  # The two pilots differ here (40 and 5), so the comparisons above can fail.
  expect_identical(seq_aipe_R2_pilot(0.2, p = 3, ci_method = "wald"), 40L)
  expect_identical(seq_aipe_R2_pilot(0.2, p = 3, ci_method = "log"), 5L)
})

test_that("a number of predictors too large for p + 2 to be an integer is refused by name", {
  # No warning about integer overflow may escape, and the message is about p.
  expect_no_warning(
    expect_error(seq_aipe_R2_pilot(0.2, p = .Machine$integer.max),
                 "'p' must be a single whole number of at most 2147483645")
  )
  expect_no_warning(
    expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = .Machine$integer.max - 1,
                                         G = 2),
                 "'p' must be a single whole number of at most")
  )
  # The largest p accepted gives a pilot of p + 2 without overflow.
  expect_identical(seq_aipe_R2_pilot(0.2, p = .Machine$integer.max - 2L),
                   .Machine$integer.max)
})

# ---------------------------------------------------------------------------
# The estimate
# ---------------------------------------------------------------------------

test_that("the estimate is the R squared of lm(), for a matrix, a data frame, and a vector", {
  d <- survey_data(45)
  expected <- summary(stats::lm(d$y ~ d$x))$r.squared
  out <- seq_aipe_R2(d$x, d$y, width = 0.20)
  expect_equal(row_value(out, "R2"), expected, tolerance = 1e-12)
  expect_identical(attr(out, "p"), 3L)
  expect_identical(attr(out, "current_n"), 45L)

  as_frame <- data.frame(autonomy = d$x[, 1], support = d$x[, 2],
                         workload = d$x[, 3])
  expect_equal(seq_aipe_R2(as_frame, d$y, width = 0.20)$value, out$value,
               tolerance = 1e-12)

  # One predictor: R2 is the squared correlation.
  one <- seq_aipe_R2(d$x[, 1], d$y, width = 0.20)
  expect_equal(row_value(one, "R2"), stats::cor(d$x[, 1], d$y)^2,
               tolerance = 1e-12)
  expect_identical(attr(one, "p"), 1L)
  expect_equal(seq_aipe_R2(d$x[, 1, drop = FALSE], d$y, width = 0.20)$value,
               one$value)
})

test_that("the estimate does not change with the location or scale of any variable", {
  d <- survey_data(40)
  base <- seq_aipe_R2(d$x, d$y, width = 0.20)
  moved <- seq_aipe_R2(sweep(d$x %*% diag(c(1e-4, 250, 3)), 2,
                             c(1000, -5000, 70), `+`),
                       1e4 + 1e3 * d$y, width = 0.20)
  # The first predictor, scaled to a spread of 1e-4 and shifted to a level
  # of 1000, keeps about nine significant digits of its spread when it is
  # stored, so the agreement is to about that precision. (The largest
  # relative difference in a row was 5.4e-9 over this data set and 60
  # others drawn by survey_data() with seeds 2001 to 2060; at a tolerance
  # of 1e-10 the test passed for this data set and failed for 6 of the 61.)
  expect_equal(moved$value, base$value, tolerance = 1e-8)
  # With a transformation that is exact in binary floating point there is
  # nothing to lose: data rounded to multiples of 1/1024, scales and shifts
  # that are powers of two. The agreement is then to rounding error (at
  # most 5.8e-14 over the same 61 data sets).
  xd <- round(d$x * 1024) / 1024
  yd <- round(d$y * 1024) / 1024
  exact <- sweep(xd %*% diag(c(2^-13, 256, 4)), 2, c(1024, -4096, 64), `+`)
  expect_identical((exact[, 1] - 1024) * 2^13, xd[, 1])
  expect_equal(seq_aipe_R2(exact, 8192 + 1024 * yd, width = 0.20)$value,
               seq_aipe_R2(xd, yd, width = 0.20)$value, tolerance = 1e-12)
  # A shift of 1.7e9 (seconds since 1970) leaves about seven significant
  # digits for the variation in the data themselves, so the agreement is to
  # that precision and no worse.
  far <- seq_aipe_R2(d$x + 1.7e9, d$y - 1.7e9, width = 0.20)
  expect_equal(far$value, base$value, tolerance = 1e-5)
  expect_identical(attr(far, "variance_unusable"), FALSE)
  expect_identical(attr(far, "is_satisfied"), attr(base, "is_satisfied"))
  # The order of the predictors does not matter either.
  expect_equal(seq_aipe_R2(d$x[, 3:1], d$y, width = 0.20)$value, base$value,
               tolerance = 1e-12)
})

test_that("the result does not depend on the units of a variable, however extreme", {
  # Multiplying a variable by a power of two is exact in binary floating
  # point, so the whole result must be the same to the last bit. Computed
  # from the values as supplied, the sums of squares of a response of the
  # order of 1e-160 underflow (R2 wrong from the seventh digit at 2^-530,
  # not a number from 2^-600) and those of a response of the order of 1e155
  # overflow (not a number from 2^510), and the mean of a predictor of the
  # order of 1e307 overflows.
  d <- survey_data(50)
  base <- seq_aipe_R2(d$x, d$y, width = 0.20)
  expect_equal(row_value(base, "R2"), summary(stats::lm(d$y ~ d$x))$r.squared,
               tolerance = 1e-12)
  expect_identical(attr(base, "notes"), .note_R2_small_sample(50L, 3L))
  # On an ordinary scale the common binary scale changes nothing at all: the
  # estimate is, to the last bit, the one computed here from the values as
  # supplied (QR decomposition of the centered predictors; model sum of
  # squares over model plus residual sums of squares).
  x_qr <- qr(sweep(d$x, 2, colMeans(d$x)))
  y_centered <- d$y - mean(d$y)
  ss_model <- sum(qr.fitted(x_qr, y_centered)^2)
  ss_resid <- sum(qr.resid(x_qr, y_centered)^2)
  expect_identical(row_value(base, "R2"), ss_model / (ss_model + ss_resid))
  for (power in c(-1000, -600, -530, 510, 600, 1000)) {
    for (ci_method in c("wald", "log_transformed")) {
      reference <- seq_aipe_R2(d$x, d$y, width = 0.20, ci_method = ci_method)
      out <- seq_aipe_R2(d$x, d$y * 2^power, width = 0.20,
                         ci_method = ci_method)
      expect_identical(out$value, reference$value)
      expect_identical(attr(out, "notes"), attr(reference, "notes"))
      out <- seq_aipe_R2(d$x * 2^power, d$y, width = 0.20,
                         ci_method = ci_method)
      expect_identical(out$value, reference$value)
    }
  }
  # Each variable in units of its own, at both ends of the range at once.
  scaled <- sweep(d$x, 2, 2^c(1020, -1000, 300), `*`)
  expect_identical(seq_aipe_R2(scaled, d$y * 2^-900, width = 0.20)$value,
                   base$value)
  # A unit that is not a power of two rounds every value in its last place,
  # so the agreement is to rounding error, and the decision and the notes
  # are the same.
  for (unit in c(1e-170, 1e160)) {
    out <- seq_aipe_R2(d$x, d$y * unit, width = 0.20)
    expect_equal(out$value, base$value, tolerance = 1e-12)
    expect_identical(attr(out, "is_satisfied"), attr(base, "is_satisfied"))
    expect_identical(attr(out, "variance_unusable"), FALSE)
    expect_identical(attr(out, "notes"), attr(base, "notes"))
  }
})

test_that("only a response that does not vary can be said to differ beyond the tenth digit", {
  # .R2_estimate() reports response_differs for the note that explains a
  # response that "does not vary" although its values are unequal. It is
  # FALSE for a response that varies, whatever its units, and for one whose
  # values are all equal.
  d <- survey_data(40)
  expect_false(.R2_estimate(d$x, d$y)$response_differs)
  expect_false(.R2_estimate(d$x, d$y * 1e-170)$response_differs)
  expect_false(.R2_estimate(d$x, d$y * 1e160)$response_differs)
  expect_false(.R2_estimate(d$x, rep(2.5, 40))$response_differs)
  rounded <- rep(0.3, 40)
  rounded[1] <- 0.1 * 3
  expect_true(.R2_estimate(d$x, rounded)$response_differs)
  # So the sentence about the tenth significant digit is never attached to
  # a response that varies.
  for (unit in c(1e-170, 1, 1e160)) {
    notes <- attr(seq_aipe_R2(d$x, d$y * unit, width = 0.2), "notes")
    expect_false(any(grepl("tenth significant digit", notes)))
  }
})

test_that("integer input gives the lm() value", {
  withr::local_seed(113)
  x <- matrix(sample.int(7L, 90L, replace = TRUE), 30L, 3L)
  y <- sample.int(5L, 30L, replace = TRUE) + x[, 1L]
  out <- seq_aipe_R2(x, y, width = 0.30)
  expect_equal(row_value(out, "R2"), summary(stats::lm(y ~ x))$r.squared,
               tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Intervals and criteria against the article's formulas
# ---------------------------------------------------------------------------

test_that("the Wald interval is Equation 59 with R in place of the population value", {
  # n = 100, R2 = .5 exactly, 95%: half width = z * 2 * sqrt(.5) * .5 / 10
  # = 0.1385903824.
  d <- exact_R2_data(100, 2, 0.5)
  out <- seq_aipe_R2(d$x, d$y, width = 0.20)
  expect_equal(row_value(out, "R2"), 0.5, tolerance = 1e-12)
  expect_equal(row_value(out, "lower_limit"), 0.3614096176, tolerance = 1e-9)
  expect_equal(row_value(out, "upper_limit"), 0.6385903824, tolerance = 1e-9)
  # Olkin and Finn's asymptotic variance, 4 P2 (1 - P2)^2 = 4 * .5 * .25.
  expect_equal(row_value(out, "xi_squared"), 0.5, tolerance = 1e-12)
})

test_that("the Wald stopping rule adds 1/n to both factors and so reproduces Tables 13 and 14", {
  # n = 100, R2 = .5, omega = .2, 95%:
  # 16 z^2 / omega^2 * (.51)(.51)^2 = 203.83, so the criterion is 204.
  d <- exact_R2_data(100, 2, 0.5)
  out <- seq_aipe_R2(d$x, d$y, width = 0.20)
  expect_identical(attr(out, "criterion_n"), 204L)
  expect_equal(row_value(out, "criterion_n"), 204)
  expect_false(attr(out, "is_satisfied"))

  # n = 50, R2 = .7, omega = .1, 90%: (.72)(.32)^2 gives 319.16, so 320.
  # Equation 63 as printed, (.72)(1 - .72)^2, would give 244.36, so 245; the
  # printed form does not reproduce the article's tables (see the operating
  # characteristic tests below and the help page).
  d <- exact_R2_data(50, 3, 0.7)
  out <- seq_aipe_R2(d$x, d$y, width = 0.10, conf_level = 0.90)
  expect_identical(attr(out, "criterion_n"), 320L)

  # The corrected criterion across the range of R2 at n = 40, omega = .10,
  # 90%: 16 z^2 / omega^2 = 4328.87, times (R2 + 1/40)(1 - R2 + 1/40)^2.
  # By hand: 196.75, 662.93, 741.39, 488.69, 62.57, 5.38. Without the
  # correction these would be 83.15, 554.10, 641.12, 415.57, 38.96, 0.43, and
  # with Equation 63 as printed 177.66, 585.01, 637.11, 380.47, 22.52, 0.99,
  # so each expected value below separates the three forms.
  R2_values <- c(0.02, 0.2, 0.34, 0.6, 0.9, 0.99)
  by_hand <- c(197L, 663L, 742L, 489L, 63L, 6L)
  for (i in seq_along(R2_values)) {
    d <- exact_R2_data(40, 2, R2_values[i])
    out <- seq_aipe_R2(d$x, d$y, width = 0.10, conf_level = 0.90)
    expect_identical(attr(out, "criterion_n"), by_hand[i])
  }
})

test_that("the Wald rule is met at the first n that reaches the criterion, not one observation earlier or later", {
  # R2 = .5 exactly, omega = .2, 95%: 16 z^2 / omega^2 = 1536.584.
  # n = 197: 1536.584 * (.5 + 1/197)^3 = 197.98 > 197, so not met.
  # n = 198: 1536.584 * (.5 + 1/198)^3 = 197.95 <= 198, so met.
  # (Equation 63 as printed would be met from n = 191: 190.04 <= 191.)
  d <- exact_R2_data(197, 2, 0.5)
  out <- seq_aipe_R2(d$x, d$y, width = 0.20)
  expect_identical(attr(out, "criterion_n"), 198L)
  expect_false(attr(out, "is_satisfied"))
  d <- exact_R2_data(198, 2, 0.5)
  out <- seq_aipe_R2(d$x, d$y, width = 0.20)
  expect_identical(attr(out, "criterion_n"), 198L)
  expect_true(attr(out, "is_satisfied"))
  # At the stopping point the interval is within the target: half width
  # z * 2 * sqrt(.5) * .5 / sqrt(198) = 0.09849, so a width of 0.19698.
  expect_equal(row_value(out, "upper_limit") - row_value(out, "lower_limit"),
               0.1969838, tolerance = 1e-6)
  d <- exact_R2_data(191, 2, 0.5)
  expect_false(attr(seq_aipe_R2(d$x, d$y, width = 0.20), "is_satisfied"))
})

test_that("the log transformed rule is met at n_BW and not one observation earlier", {
  # Equation 62 at R2 = .9, k = 2, omega = .10, 90%:
  # 2 + 2 + 4 (.9) z^2 / ln[(1/2)(1 + sqrt(5))]^2 = 46.06, so n_BW = 47, the
  # value printed in Table 11.
  d <- exact_R2_data(46, 2, 0.9)
  out <- seq_aipe_R2(d$x, d$y, 0.10, 0.90, "log_transformed")
  expect_identical(attr(out, "criterion_n"), 47L)
  expect_false(attr(out, "is_satisfied"))
  d <- exact_R2_data(47, 2, 0.9)
  out <- seq_aipe_R2(d$x, d$y, 0.10, 0.90, "log_transformed")
  expect_identical(attr(out, "criterion_n"), 47L)
  expect_true(attr(out, "is_satisfied"))
})

test_that("with a very large n the Wald criterion is the article's optimal sample size (Equation 61, Tables 13 and 14)", {
  # n_OF as printed in the tables: 1,403 (P2 = .1, omega = .05, 90%),
  # 56 (P2 = .9, omega = .10, 95%), 3,074 (P2 = .5, omega = .05, 95%),
  # 637 (P2 = .3, omega = .10, 90%).
  n_of <- function(P2, width, z) {
    ceiling(.R2_rule(P2, 1 - P2, n = 1e12, p = 2, width = width, z = z,
                     ci_method = "wald")$criterion)
  }
  expect_equal(n_of(0.1, 0.05, z90), 1403)
  expect_equal(n_of(0.9, 0.10, z95), 56)
  expect_equal(n_of(0.5, 0.05, z95), 3074)
  expect_equal(n_of(0.3, 0.10, z90), 637)
})

test_that("the log transformed interval is Equation 58 (Bonett and Wright, 2011)", {
  # n = 100, p = 2, R2 = .5, 95%:
  # 1 - exp(ln(.5) +/- z * 2 * sqrt(.5) / sqrt(96)) = 0.3365161373 and
  # 0.6232010844.
  d <- exact_R2_data(100, 2, 0.5)
  out <- seq_aipe_R2(d$x, d$y, width = 0.20, ci_method = "log_transformed")
  expect_equal(row_value(out, "lower_limit"), 0.3365161373, tolerance = 1e-9)
  expect_equal(row_value(out, "R2"), 0.5, tolerance = 1e-12)
  expect_equal(row_value(out, "upper_limit"), 0.6232010844, tolerance = 1e-9)
  expect_equal(row_value(out, "xi_squared"), 0.5, tolerance = 1e-12)
  # Equation 62 at R2 = .5, p = 2, omega = .2: 2 + 2 + 4 * .5 * z^2 / L^2
  # with L = ln[(1/2)(.4 + sqrt(.16 + 4))] is 198.61, so 199.
  expect_identical(attr(out, "criterion_n"), 199L)
})

test_that("the log transformed criterion is the article's optimal sample size at R2 = P2 (Equation 60, Tables 11 and 12)", {
  # n_BW as printed: 1,407 (k = 2, omega = .05, P2 = .1, 90%); 72 (k = 10,
  # omega = .10, P2 = .9, 95%); 3,083 (k = 5, omega = .05, P2 = .5, 95%);
  # 47 (k = 2, omega = .10, P2 = .9, 90%). The rule has no correction term,
  # so the criterion at a sample R2 equal to P2 is n_BW at any n.
  d <- exact_R2_data(40, 2, 0.1)
  expect_identical(
    attr(seq_aipe_R2(d$x, d$y, 0.05, 0.90, "log_transformed"),
         "criterion_n"), 1407L)
  d <- exact_R2_data(40, 10, 0.9)
  expect_identical(
    attr(seq_aipe_R2(d$x, d$y, 0.10, 0.95, "log_transformed"),
         "criterion_n"), 72L)
  d <- exact_R2_data(40, 5, 0.5)
  expect_identical(
    attr(seq_aipe_R2(d$x, d$y, 0.05, 0.95, "log_transformed"),
         "criterion_n"), 3083L)
  d <- exact_R2_data(60, 2, 0.9)
  out <- seq_aipe_R2(d$x, d$y, 0.10, 0.90, "log_transformed")
  expect_identical(attr(out, "criterion_n"), 47L)
  expect_true(attr(out, "is_satisfied"))   # 60 >= 47
})

test_that("the limits are not truncated to [0, 1]", {
  # Wald, n = 20, R2 = .05, 95%: .05 -/+ z * 2 * sqrt(.05) * .95 / sqrt(20).
  d <- exact_R2_data(20, 2, 0.05)
  out <- seq_aipe_R2(d$x, d$y, width = 0.20)
  expect_equal(row_value(out, "lower_limit"), -0.1361965785, tolerance = 1e-9)
  expect_equal(row_value(out, "upper_limit"), 0.2361965785, tolerance = 1e-9)
  # Wald upper limit above one, which needs n < 4 z^2 R2 (14.6 here):
  # n = 8, R2 = .95, 95%: .95 + z * 2 * sqrt(.95) * .05 / sqrt(8).
  d <- exact_R2_data(8, 2, 0.95)
  out <- seq_aipe_R2(d$x, d$y, width = 0.20)
  expect_equal(row_value(out, "upper_limit"), 1.017540598, tolerance = 1e-9)
  # Log transformed, n = 8, p = 2, R2 = .3, 95%:
  # 1 - .7 exp(+/- z * 2 * sqrt(.3) / 2).
  d <- exact_R2_data(8, 2, 0.3)
  out <- seq_aipe_R2(d$x, d$y, width = 0.20, ci_method = "log_transformed")
  expect_equal(row_value(out, "lower_limit"), -1.047954607, tolerance = 1e-9)
  expect_equal(row_value(out, "upper_limit"), 0.7607368843, tolerance = 1e-9)
})

test_that("the rule is met exactly when n reaches the criterion, and the interval is then no wider than the target", {
  withr::local_seed(113)
  met <- 0L
  for (i in 1:200) {
    p <- sample(1:4, 1)
    n <- sample(seq.int(p + 3L, 400L), 1)
    width <- stats::runif(1, 0.1, 0.6)
    ci_method <- sample(c("wald", "log_transformed"), 1)
    x <- matrix(stats::rnorm(n * p), n, p)
    y <- drop(x %*% rep(stats::runif(1, 0, 1), p)) + stats::rnorm(n)
    out <- seq_aipe_R2(x, y, width = width, ci_method = ci_method)
    expect_identical(attr(out, "is_satisfied"),
                     n >= attr(out, "criterion_n"))
    if (attr(out, "is_satisfied")) {
      met <- met + 1L
      expect_lte(row_value(out, "upper_limit") - row_value(out, "lower_limit"),
                 width + 1e-12)
    }
  }
  expect_gt(met, 30L)   # the property was exercised
})

# ---------------------------------------------------------------------------
# The returned object
# ---------------------------------------------------------------------------

test_that("seq_aipe_R2 returns the standard smsd_check", {
  d <- survey_data(60)
  out <- seq_aipe_R2(d$x, d$y, width = 0.25)
  expect_s3_class(out, c("smsd_check", "data.frame"), exact = TRUE)
  expect_identical(out$term, c("lower_limit", "R2", "upper_limit",
                               "xi_squared", "criterion_n", "current_n"))
  expect_type(out$value, "double")
  expect_identical(attr(out, "method"),
                   "squared multiple correlation, Wald interval")
  expect_identical(attr(out, "ci_method"), "wald")
  expect_identical(attr(out, "p"), 3L)
  expect_identical(attr(out, "conf_level"), 0.95)
  expect_identical(attr(out, "width"), 0.25)
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_identical(attr(out, "notes"), character(0))
  expect_identical(row_value(out, "current_n"), 60)

  log_out <- seq_aipe_R2(d$x, d$y, width = 0.25,
                         ci_method = "log_transformed")
  expect_identical(attr(log_out, "method"),
                   "squared multiple correlation, log transformed interval")
  expect_identical(attr(log_out, "ci_method"), "log_transformed")
  # The two intervals share the estimate and differ in their limits.
  expect_identical(row_value(log_out, "R2"), row_value(out, "R2"))
  expect_false(isTRUE(all.equal(row_value(log_out, "lower_limit"),
                                row_value(out, "lower_limit"))))

  tidied <- generics::tidy(out)
  expect_identical(tidied$term, "R2")
  expect_identical(tidied$estimate, row_value(out, "R2"))
  expect_identical(tidied$ci_lower, row_value(out, "lower_limit"))
  expect_identical(tidied$ci_upper, row_value(out, "upper_limit"))
  expect_output(print(out), "squared multiple correlation, Wald interval")
})

test_that("seq_path() needs the pilot sample size for this procedure, as the help page says", {
  d <- survey_data(200)
  # The pilot sample size depends on the number of predictors, which the
  # check function reads from the data and seq_path() cannot hand to the
  # pilot function: the default start is an error that names the argument.
  expect_error(seq_path(seq_aipe_R2, list(x = d$x, y = d$y), width = 0.25),
               "seq_aipe_R2_pilot.*'p'")
  # The call the help page gives.
  path <- seq_path(seq_aipe_R2, list(x = d$x, y = d$y), width = 0.25,
                   start = seq_aipe_R2_pilot(0.25, p = ncol(d$x)))
  expect_identical(attr(path, "start"), 32L)   # 4 z / omega = 31.36
  # It stops where a loop of seq_aipe_R2() calls from the pilot sample size
  # first meets the rule.
  n <- 32L
  while (n < 200L &&
         !attr(seq_aipe_R2(d$x[seq_len(n), ], d$y[seq_len(n)], width = 0.25),
               "is_satisfied")) {
    n <- n + 1L
  }
  expect_lt(n, 200L)
  expect_identical(attr(path, "stop_n"), n)
  # It stops at n = 151. That many observations for 3 predictors are not
  # few (the threshold is 60), so the check at the stopping point has no
  # note.
  expect_gte(n, 60L)
  expect_identical(attr(attr(path, "final_check"), "notes"), character(0))

  # A path that stops while the sample is small relative to the number of
  # predictors: the help page sends the reader to the check at the stopping
  # point, attr(path, "final_check"), for the note. Ten predictors at a
  # width of .30 start at 27 (4 z / omega = 26.13) and stop here at n = 65,
  # far below 10 * 20 = 200 observations.
  d10 <- survey_data(200, p = 10)
  path <- seq_path(seq_aipe_R2, list(x = d10$x, y = d10$y), width = 0.30,
                   start = seq_aipe_R2_pilot(0.30, p = 10))
  stop_n <- attr(path, "stop_n")
  expect_identical(attr(path, "start"), 27L)
  expect_lt(stop_n, 200L)
  final <- attr(path, "final_check")
  expect_s3_class(final, "smsd_check")
  expect_identical(attr(final, "current_n"), stop_n)
  expect_true(attr(final, "is_satisfied"))
  expect_identical(attr(final, "notes"), .note_R2_small_sample(stop_n, 10L))
  expect_identical(final$value,
                   seq_aipe_R2(d10$x[seq_len(stop_n), ], d10$y[seq_len(stop_n)],
                               width = 0.30)$value)
})

# ---------------------------------------------------------------------------
# Guards
# ---------------------------------------------------------------------------

test_that("too few observations for the number of predictors is an error", {
  d <- survey_data(20)
  # p + 1 observations make the sample R2 identically one.
  expect_error(seq_aipe_R2(d$x[1:4, ], d$y[1:4], width = 0.2),
               "needs at least 5 observations")
  expect_error(seq_aipe_R2(d$x[1:2, ], d$y[1:2], width = 0.2),
               "needs at least 5 observations")
  expect_error(seq_aipe_R2(d$x[1:2, 1], d$y[1:2], width = 0.2),
               "needs at least 3 observations")
  # p + 2 is accepted.
  expect_s3_class(seq_aipe_R2(d$x[1:5, ], d$y[1:5], width = 0.2),
                  "smsd_check")
  # No data at all is the same condition, reported in the same words.
  expect_error(seq_aipe_R2(matrix(numeric(0), 0, 3), numeric(0), width = 0.2),
               "needs at least 5 observations.*there are 0")
  expect_error(seq_aipe_R2(numeric(0), numeric(0), width = 0.2),
               "needs at least 3 observations.*there are 0")
})

test_that("collinear or constant predictors are an error", {
  d <- survey_data(30)
  expect_error(seq_aipe_R2(cbind(d$x, d$x[, 1] + d$x[, 2]), d$y, width = 0.2),
               "collinear: 4 columns were supplied but only 3")
  expect_error(seq_aipe_R2(cbind(d$x, d$x[, 1]), d$y, width = 0.2),
               "collinear")
  expect_error(seq_aipe_R2(cbind(a = d$x[, 1], b = 2), d$y, width = 0.2),
               "does not vary in the current data \\(b of 'x'\\)")
  expect_error(seq_aipe_R2(cbind(d$x[, 1], 0.1), d$y, width = 0.2),
               "does not vary in the current data \\(column 2 of 'x'\\)")
  # A matrix with names for some of its columns only, as cbind(1, x) gives
  # for a named x. The column without a name is given by its position (the
  # message used to name nothing: "( of 'x')"), and a column of ones gets a
  # word about the intercept, which the function fits itself.
  named <- d$x[, 1:2]
  colnames(named) <- c("a", "b")
  expect_identical(colnames(cbind(1, named)), c("", "a", "b"))
  err <- tryCatch(seq_aipe_R2(cbind(1, named), d$y, width = 0.2),
                  error = function(e) conditionMessage(e))
  expect_match(err, "does not vary in the current data (column 1 of 'x')",
               fixed = TRUE)
  expect_match(err, paste("(The function fits the intercept itself, so a",
                          "column of ones is not needed.)"), fixed = TRUE)
  expect_error(seq_aipe_R2(cbind(named, 1, 2), d$y, width = 0.2),
               "\\(column 3, column 4 of 'x'\\)")
  expect_error(seq_aipe_R2(cbind(named, 3, late = 2), d$y, width = 0.2),
               "\\(column 3, late of 'x'\\)")
  # A constant that is not one says nothing about an intercept.
  err <- tryCatch(seq_aipe_R2(cbind(a = d$x[, 1], b = 2), d$y, width = 0.2),
                  error = function(e) conditionMessage(e))
  expect_false(grepl("intercept", err, fixed = TRUE))
  expect_identical(.R2_column_labels(matrix(0, 2, 3), c(1, 3)),
                   "column 1, 3")
  # A predictor that is constant except for the rounding of its values (0.3
  # and 0.1 * 3 differ by one unit in the last place) does not vary either.
  rounded <- rep(0.3, 30)
  rounded[1] <- 0.1 * 3
  expect_false(identical(rounded[1], rounded[2]))
  expect_error(seq_aipe_R2(cbind(d$x[, 1], rounded), d$y, width = 0.2),
               "does not vary in the current data \\(rounded of 'x'\\)")
  # An indicator for every category of a factor is collinear with the
  # intercept.
  group <- rep(1:3, 10)
  expect_error(seq_aipe_R2(cbind(group == 1, group == 2, group == 3) + 0,
                           d$y, width = 0.2), "collinear")
})

test_that("missing, infinite, and mismatched data are errors", {
  d <- survey_data(30)
  x_na <- d$x
  x_na[3, 2] <- NA
  expect_error(seq_aipe_R2(x_na, d$y, width = 0.2),
               "'x' contains missing values")
  y_na <- d$y
  y_na[7] <- NA
  expect_error(seq_aipe_R2(d$x, y_na, width = 0.2),
               "'y' contains missing values")
  x_inf <- d$x
  x_inf[1, 1] <- Inf
  expect_error(seq_aipe_R2(x_inf, d$y, width = 0.2),
               "'x' contains infinite values")
  y_inf <- d$y
  y_inf[1] <- -Inf
  expect_error(seq_aipe_R2(d$x, y_inf, width = 0.2),
               "'y' contains infinite values")
  expect_error(seq_aipe_R2(d$x, d$y[-1], width = 0.2),
               "'x' has 30 rows and 'y' has 29 values")
  expect_error(seq_aipe_R2(d$x[-1, ], d$y, width = 0.2),
               "'x' has 29 rows and 'y' has 30 values")
  expect_error(seq_aipe_R2(d$x, numeric(0), width = 0.2), "'y'")
})

test_that("inputs of the wrong kind are errors that name the argument", {
  d <- survey_data(30)
  expect_error(seq_aipe_R2(y = d$y, width = 0.2), "must be supplied")
  expect_error(seq_aipe_R2(d$x, width = 0.2), "must be supplied")
  expect_error(seq_aipe_R2(as.character(d$x[, 1]), d$y, width = 0.2), "'x'")
  expect_error(seq_aipe_R2(matrix(letters[1:30], 30, 1), d$y, width = 0.2),
               "'x' must be a numeric matrix")
  expect_error(seq_aipe_R2(array(1, c(30, 2, 2)), d$y, width = 0.2), "'x'")
  expect_error(seq_aipe_R2(list(d$x[, 1]), d$y, width = 0.2), "'x'")
  frame <- data.frame(score = d$x[, 1], group = rep(c("a", "b"), 15))
  expect_error(seq_aipe_R2(frame, d$y, width = 0.2),
               "Every column of 'x' must be numeric; not numeric: group")
  expect_error(seq_aipe_R2(data.frame(), d$y, width = 0.2),
               "at least one predictor column")
  expect_error(seq_aipe_R2(d$x, matrix(d$y, ncol = 1), width = 0.2),
               "'y' must be a numeric vector")
  expect_error(seq_aipe_R2(d$x, as.character(d$y), width = 0.2), "'y'")
  expect_error(seq_aipe_R2(d$x, d$y), "'width' must be supplied")
  expect_error(seq_aipe_R2(d$x, d$y, width = 0), "'width'")
  expect_error(seq_aipe_R2(d$x, d$y, width = 1.2), "'width'.*less than 1")
  expect_error(seq_aipe_R2(d$x, d$y, width = NA), "'width'")
  expect_error(seq_aipe_R2(d$x, d$y, width = 0.2, conf_level = 1.5),
               "'conf_level'")
  expect_error(seq_aipe_R2(d$x, d$y, width = 0.2, conf_level = NA),
               "'conf_level'")
  expect_error(seq_aipe_R2(d$x, d$y, width = 0.2, ci_method = "fisher"),
               paste0("'ci_method' must be one of \"wald\", ",
                      "\"log_transformed\"; got \"fisher\""))
})

# ---------------------------------------------------------------------------
# When no decision can be made the rule is "not met", never NA and never
# "met": the limits are NA, the criterion is at least n + 1, and a note says
# why (as in the other seq_aipe_* check functions)
# ---------------------------------------------------------------------------

test_that("a response that does not vary is not estimable and not met", {
  d <- survey_data(40)
  for (ci_method in c("wald", "log_transformed")) {
    out <- seq_aipe_R2(d$x, rep(2.5, 40), width = 0.9, ci_method = ci_method)
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_identical(attr(out, "variance_unusable"), TRUE)
    # The two notes of the whole package: the estimate cannot be computed,
    # and so neither can its variance.
    expect_identical(attr(out, "notes"),
                     c(.note_not_estimable(), .note_variance_unusable()))
    expect_true(all(is.na(
      out$value[!out$term %in% c("criterion_n", "current_n")]
    )))
    expect_identical(row_value(out, "current_n"), 40)
    # At least one more observation is needed: the criterion is n + 1, an
    # integer, never NA.
    expect_identical(attr(out, "criterion_n"), 41L)
    expect_identical(row_value(out, "criterion_n"), 41)
    expect_identical(attr(out, "current_n"), 40L)
    expect_output(print(out), "Stopping rule met: FALSE")
  }
})

test_that("a response that is constant except for rounding is not estimable either", {
  # 0.3 and 0.1 * 3 differ by one unit in the last place (5.6e-17). Taken at
  # face value the response would "vary", and the sample R2 of that rounding
  # noise on two predictors would be reported as an ordinary estimate.
  d <- survey_data(200, p = 2)
  y <- rep(0.3, 200)
  y[1] <- 0.1 * 3
  expect_gt(diff(range(y)), 0)
  for (ci_method in c("wald", "log_transformed")) {
    out <- seq_aipe_R2(d$x, y, width = 0.6, ci_method = ci_method)
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_identical(attr(out, "variance_unusable"), TRUE)
    expect_true(is.na(row_value(out, "R2")))
    # Because the values are not all equal, a note between the two standard
    # ones says what "does not vary" means.
    notes <- attr(out, "notes")
    expect_length(notes, 3L)
    expect_identical(notes[c(1L, 3L)],
                     c(.note_not_estimable(), .note_variance_unusable()))
    expect_match(notes[2L],
                 paste0("^The values of the response differ only beyond the ",
                        "tenth significant digit, which counts as not ",
                        "varying"))
  }
})

test_that("a perfect fit has a zero variance estimate and is not met", {
  # As printed, the rules would be met at once here, with an interval of
  # zero width at 1: with R2 = 1, n = 40, omega = .5, and 95% the Wald
  # criterion is (4 z / .5)^2 (1 + 1/40) (1/40)^2 = 0.16, and the log
  # transformed criterion is p + 2 = 5.
  expect_lt((4 * z95 / 0.5)^2 * (1 + 1 / 40) * (1 / 40)^2, 1)
  d <- survey_data(40)
  y <- 3 + drop(d$x %*% c(1, -2, 0.5))
  for (ci_method in c("wald", "log_transformed")) {
    out <- seq_aipe_R2(d$x, y, width = 0.5, ci_method = ci_method)
    expect_equal(row_value(out, "R2"), 1)
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_identical(attr(out, "variance_unusable"), TRUE)
    expect_identical(attr(out, "notes"), .note_variance_unusable())
    # No interval is built on an unusable variance estimate, and at least
    # one more observation is needed.
    expect_true(is.na(row_value(out, "lower_limit")))
    expect_true(is.na(row_value(out, "upper_limit")))
    expect_identical(attr(out, "criterion_n"), 41L)
    expect_identical(row_value(out, "criterion_n"), 41)
    # The variance estimate is reported as computed; nothing bounds it below.
    expect_lt(row_value(out, "xi_squared"), 40^(-3))
  }
  # The criterion is raised to n + 1 only when the rule's own criterion is
  # below it. At a width of .01 the Wald criterion for R2 = 1 at n = 40 is
  # (4 z / .01)^2 (1 + 1/40) (1/40)^2 = 393.75, so 394 is reported, not 41.
  out <- seq_aipe_R2(d$x, y, width = 0.01)
  expect_identical(attr(out, "is_satisfied"), FALSE)
  expect_identical(attr(out, "variance_unusable"), TRUE)
  expect_true(is.na(row_value(out, "lower_limit")))
  expect_identical(attr(out, "criterion_n"), 394L)
  expect_identical(row_value(out, "criterion_n"), 394)
})

test_that("a perfect fit is not met when the variables carry a large offset", {
  # With an offset of 1e10 the stored values are rounded to about 2e-6, so
  # the sample R2 of an exact linear response is 1 only to within rounding:
  # 1 - R2 is anywhere from exactly 0 to about 1e-11, depending on the data
  # and on whether R accumulates sums in extended precision. As printed,
  # both rules would be met at once (criterion 1 and 5) with an interval of
  # essentially zero width.
  withr::local_seed(113)
  spreads <- replicate(5, stats::rnorm(40), simplify = FALSE)
  for (spread in spreads) {
    for (offset in c(1e8, 1e9, 1e10)) {
      x <- offset + spread
      expect_false(.is_constant(x))
      # A slope that is a power of two: the stored response is an exact
      # linear function of the stored predictor (a true residual share of
      # 0), and what is left is the last-bit error of colMeans() where R
      # has no long double, and nothing where it has. And a slope that is
      # not (0.3): the stored response is itself rounded, on every platform.
      for (slope in c(2, 0.3)) {
        y <- slope * x + 1
        for (ci_method in c("wald", "log_transformed")) {
          out <- seq_aipe_R2(x, y, width = 0.2, ci_method = ci_method)
          expect_identical(attr(out, "is_satisfied"), FALSE)
          expect_identical(attr(out, "variance_unusable"), TRUE)
          expect_identical(attr(out, "notes"), .note_variance_unusable())
          expect_identical(attr(out, "criterion_n"), 41L)
          expect_true(is.na(row_value(out, "lower_limit")))
        }
      }
    }
  }
  # One of these as numbers: for the first draw at an offset of 1e10 with a
  # slope of 2, 1 - R2 is 3.5e-12 where capabilities("long.double") is
  # FALSE (macOS arm64, a build without long double) and exactly 0 where
  # it is TRUE (x86_64). Only the upper bound is asserted.
  x <- 1e10 + spreads[[1L]]
  out <- seq_aipe_R2(x, 2 * x + 1, width = 0.2)
  expect_lt(1 - row_value(out, "R2"), 1e-10)
  # Several predictors, each with its own offset and scale.
  x <- cbind(1e9 + stats::rnorm(60), -3e10 + 50 * stats::rnorm(60),
             stats::rnorm(60))
  y <- 7e9 + drop(x %*% c(1.5, -0.25, 3))
  out <- seq_aipe_R2(x, y, width = 0.2)
  expect_identical(attr(out, "is_satisfied"), FALSE)
  expect_identical(attr(out, "variance_unusable"), TRUE)
})

test_that("a variable whose values agree to ten significant digits does not vary (the package's one definition)", {
  d <- survey_data(40, p = 2)
  # An offset of 1e11 on a spread of about 4: the range is 4e-11 of the
  # largest value, inside the package's tolerance of 1e-10. The spread is
  # 40 further draws of the stream of seed 113 (survey_data() used its first
  # 120 for the two predictors and the error), so that it is not the first
  # predictor over again.
  far <- 1e11 + withr::with_seed(113, stats::rnorm(160))[121:160]
  expect_gt(diff(range(far)), 1)
  expect_true(.is_constant(far))

  # As a predictor it is an error that says why and what to do.
  expect_error(
    seq_aipe_R2(cbind(d$x[, 1], far), d$y, width = 0.2),
    paste0("does not vary in the current data \\(far of 'x'\\).*The values ",
           "of far differ only beyond the tenth significant digit, which ",
           "counts as not varying\\. If it is a large constant plus a small ",
           "spread, subtract the constant")
  )
  # Subtracting the constant, as the message says, gives lm()'s value.
  near <- far - 1e11
  out <- seq_aipe_R2(cbind(d$x[, 1], near), d$y, width = 0.2)
  expect_equal(row_value(out, "R2"),
               summary(stats::lm(d$y ~ d$x[, 1] + near))$r.squared,
               tolerance = 1e-12)
  # A predictor whose values are all equal gets the same error without the
  # explanation, which would not apply.
  message_equal <- tryCatch(
    seq_aipe_R2(cbind(d$x[, 1], 0.1), d$y, width = 0.2),
    error = function(e) conditionMessage(e)
  )
  expect_match(message_equal, "does not vary in the current data")
  expect_false(grepl("tenth significant digit", message_equal))

  # As the response it is "not estimable", with a note that says why.
  out <- seq_aipe_R2(d$x, far, width = 0.2)
  expect_identical(attr(out, "is_satisfied"), FALSE)
  expect_true(is.na(row_value(out, "R2")))
  expect_length(attr(out, "notes"), 3L)
  expect_match(attr(out, "notes")[2L],
               "^The values of the response differ only beyond the tenth")
  # One tenth of that offset is ordinary data.
  usual <- seq_aipe_R2(d$x, far - 9e10, width = 0.2)
  expect_identical(attr(usual, "variance_unusable"), FALSE)
  expect_equal(row_value(usual, "R2"),
               summary(stats::lm(I(far - 1e11) ~ d$x))$r.squared,
               tolerance = 1e-6)

  # The boundary is that of .is_constant(): a range of 0.9e-10 of the
  # largest value does not vary, a range of 1.5e-10 does.
  just_under <- 1 + 0.9e-10 * (0:39) / 39
  just_over <- 1 + 1.5e-10 * (0:39) / 39
  expect_error(seq_aipe_R2(just_under, d$y, width = 0.2), "does not vary")
  expect_s3_class(seq_aipe_R2(just_over, d$y, width = 0.2), "smsd_check")
})

test_that("a response with real error is an ordinary estimate at a large offset", {
  # The same offsets with a residual standard deviation of .01 (a residual
  # share of about 2.5e-5): the estimate is usable and must equal the value
  # computed after removing the offset, which is exact in floating point.
  withr::local_seed(113)
  shift <- stats::rnorm(40)
  error <- 0.01 * stats::rnorm(40)
  x <- 1e10 + shift
  y <- 2 * x + 1 + error
  expected <- summary(stats::lm(I(y - (2e10 + 1)) ~ I(x - 1e10)))$r.squared
  expect_gt(1 - expected, 1e-6)
  for (ci_method in c("wald", "log_transformed")) {
    out <- seq_aipe_R2(x, y, width = 0.2, ci_method = ci_method)
    expect_identical(attr(out, "variance_unusable"), FALSE)
    expect_identical(attr(out, "notes"), character(0))
    expect_equal(1 - row_value(out, "R2"), 1 - expected, tolerance = 1e-4)
    expect_true(attr(out, "is_satisfied"))
  }
})

test_that("the precision of the data is that of the arithmetic unless a mean is far from zero", {
  # R2 or 1 - R2 at or below `resolution` is treated as zero. It is
  # .Machine$double.eps for ordinary data, including a predictor on the scale
  # of a calendar year, and grows as (100 eps kappa)^2 once the largest
  # kappa = max|v| / sd(v) passes sqrt(eps) / (100 eps) = 671,089.
  eps <- .Machine$double.eps
  d <- survey_data(60)
  expect_identical(.R2_estimate(d$x, d$y)$resolution, eps)
  year <- d$x
  year[, 1] <- 2020 + year[, 1]
  expect_identical(.R2_estimate(year, d$y)$resolution, eps)
  expect_identical(.R2_estimate(d$x + 1e5, d$y)$resolution, eps)

  far <- d$x
  far[, 2] <- 1e9 + far[, 2]
  kappa <- max(apply(cbind(far, d$y), 2,
                     function(v) max(abs(v)) / stats::sd(v)))
  expect_gt(kappa, 1e8)
  expect_equal(.R2_estimate(far, d$y)$resolution, (100 * eps * kappa)^2,
               tolerance = 1e-6)
  # It does not depend on which variable carries the offset, or on the units
  # of the variables.
  expect_equal(.R2_estimate(d$x, 1e9 * stats::sd(d$y) / stats::sd(far[, 2]) +
                              d$y)$resolution,
               .R2_estimate(far, d$y)$resolution, tolerance = 0.05)
  expect_equal(.R2_estimate(1000 * far, d$y / 1000)$resolution,
               .R2_estimate(far, d$y)$resolution, tolerance = 1e-6)
})

test_that("a sample R2 of exactly zero has a zero variance estimate and is not met", {
  # As printed, the rules would be met: with R2 = 0, n = 60, omega = .5, and
  # 95% the Wald criterion is (4 z / .5)^2 (1/60) (1 + 1/60)^2 = 4.24, and
  # the log transformed criterion is p + 2 = 4.
  expect_lt((4 * z95 / 0.5)^2 * (1 / 60) * (1 + 1 / 60)^2, 5)
  d <- exact_R2_data(60, 2, 0)
  for (ci_method in c("wald", "log_transformed")) {
    out <- seq_aipe_R2(d$x, d$y, width = 0.5, ci_method = ci_method)
    expect_lt(row_value(out, "R2"), 1e-20)
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_identical(attr(out, "variance_unusable"), TRUE)
    expect_identical(attr(out, "notes"), .note_variance_unusable())
    expect_true(is.na(row_value(out, "lower_limit")))
    expect_true(is.na(row_value(out, "upper_limit")))
    expect_identical(attr(out, "criterion_n"), 61L)
  }
  # A narrower target, whose criterion is above n + 1, keeps the criterion
  # of the rule: (4 z / .05)^2 (1/60) (1 + 1/60)^2 = 423.53, so 424, not 61.
  out <- seq_aipe_R2(d$x, d$y, width = 0.05)
  expect_identical(attr(out, "is_satisfied"), FALSE)
  expect_identical(attr(out, "variance_unusable"), TRUE)
  expect_true(is.na(row_value(out, "lower_limit")))
  expect_identical(attr(out, "criterion_n"), 424L)
})

test_that("a small positive variance estimate is used as it stands: no n^(-3) bound", {
  # R2 = .999 at n = 12 has 4 R2 (1 - R2)^2 = 4.0e-6, far below
  # 12^(-3) = 5.8e-4, and is an ordinary estimate here: the interval is
  # built and the rule is evaluated. (With 12 observations and 2 predictors
  # the only note is the one about a small sample.)
  d <- exact_R2_data(12, 2, 0.999)
  out <- seq_aipe_R2(d$x, d$y, width = 0.2, ci_method = "log_transformed")
  expect_equal(row_value(out, "xi_squared"), 4 * 0.999 * 0.001^2,
               tolerance = 1e-9)
  expect_lt(row_value(out, "xi_squared"), 12^(-3))
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_true(all(is.finite(out$value)))
  expect_false(.note_variance_unusable() %in% attr(out, "notes"))
  expect_identical(attr(out, "notes"), .note_R2_small_sample(12L, 2L))
  # Equation 62 at R2 = .999: 4 + 4 (.999) z^2 / asinh(.2 / .002)^2 = 4.55,
  # so the criterion is 5 and the rule is met at n = 12.
  expect_identical(attr(out, "criterion_n"), 5L)
  expect_true(attr(out, "is_satisfied"))
})

test_that("the log transformed interval is not defined at n = p + 2 and is not met there", {
  d <- survey_data(5)
  out <- seq_aipe_R2(d$x, d$y, width = 0.9, ci_method = "log_transformed")
  expect_true(is.na(row_value(out, "lower_limit")))
  expect_true(is.na(row_value(out, "upper_limit")))
  expect_true(is.finite(row_value(out, "R2")))
  expect_identical(attr(out, "is_satisfied"), FALSE)
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_length(attr(out, "notes"), 1L)
  expect_match(attr(out, "notes"), "not defined with p \\+ 2 observations")
  # The criterion is that of the rule, Equation 62, computed here from the
  # R squared of lm(): p + 2 + 4 R2 z^2 / L^2 = 10.35 at R2 = .790, so 11.
  # It is above n + 1 = 6, so it is reported as it is.
  R2 <- summary(stats::lm(d$y ~ d$x))$r.squared
  L <- log(0.5 * (0.9 / (1 - R2) + sqrt(0.9^2 / (1 - R2)^2 + 4)))
  expect_equal(5 + 4 * R2 * z95^2 / L^2, 10.3531, tolerance = 1e-5)
  expect_identical(attr(out, "criterion_n"), 11L)
  # The Wald interval is defined there (and carries only the note about a
  # small sample: 5 observations for 3 predictors).
  wald <- seq_aipe_R2(d$x, d$y, width = 0.9)
  expect_true(all(is.finite(wald$value)))
  expect_identical(attr(wald, "notes"), .note_R2_small_sample(5L, 3L))
  # One more observation and the log transformed interval is defined.
  d6 <- survey_data(6)
  expect_true(all(is.finite(
    seq_aipe_R2(d6$x, d6$y, width = 0.9, ci_method = "log_transformed")$value
  )))
})

test_that("is_satisfied is TRUE exactly when current_n reaches criterion_n", {
  # The criterion sample size is rounded up by the shared .as_sample_size(),
  # which takes a criterion a few parts in 1e11 above a whole number for
  # that whole number. The decision must use the same rounding.
  # R2 = .5 exactly at n = 198: the Wald criterion is
  # (4 z / omega)^2 (.5 + 1/198)^3, so the width that gives a chosen
  # criterion is 4 z sqrt((.5 + 1/198)^3 / criterion).
  d <- exact_R2_data(198, 2, 0.5)
  z <- stats::qnorm(0.975)
  width_for <- function(criterion) {
    4 * z * sqrt((0.5 + 1 / 198)^3 / criterion)
  }
  # 198 (1 + 5e-11) = 198.0000000099: inside the allowance, so 198, met.
  out <- seq_aipe_R2(d$x, d$y, width = width_for(198 * (1 + 5e-11)))
  expect_identical(attr(out, "criterion_n"), 198L)
  expect_true(attr(out, "is_satisfied"))
  # 198 (1 + 1e-8) = 198.00000198: outside it, so 199, not met.
  out <- seq_aipe_R2(d$x, d$y, width = width_for(198 * (1 + 1e-8)))
  expect_identical(attr(out, "criterion_n"), 199L)
  expect_false(attr(out, "is_satisfied"))
  # 198 (1 - 1e-8): 198, met.
  out <- seq_aipe_R2(d$x, d$y, width = width_for(198 * (1 - 1e-8)))
  expect_identical(attr(out, "criterion_n"), 198L)
  expect_true(attr(out, "is_satisfied"))

  # The comparison the engine and the check function share, against the
  # shared rounding, around whole numbers of every size.
  compared <- 0L
  for (whole in c(5, 56, 198, 3074, 1e5, 2e9)) {
    for (delta in c(-1e-6, -5e-9, -5e-10, -5e-11, -5e-12, 0, 5e-12, 5e-11,
                    5e-10, 5e-9, 1e-6)) {
      criterion <- whole * (1 + delta)
      for (n in whole + (-1:1)) {
        expect_identical(.R2_reached(n, criterion),
                         n >= .as_sample_size(criterion))
        compared <- compared + 1L
      }
    }
  }
  expect_identical(compared, 198L)
  # Both answers occur at n = whole (so the comparison above can fail).
  expect_true(.R2_reached(56, 56 * (1 + 5e-11)))
  expect_false(.R2_reached(56, 56 * (1 + 5e-10)))

  # And in every state in which no decision can be made (see the tests
  # above): the criterion is above n, so the identity holds there as well.
  d <- survey_data(40)
  states <- list(
    seq_aipe_R2(d$x, rep(1, 40), width = 0.5),
    seq_aipe_R2(d$x, drop(d$x %*% c(1, 2, 3)), width = 0.5),
    seq_aipe_R2(d$x, drop(d$x %*% c(1, 2, 3)), width = 0.5,
                ci_method = "log_transformed"),
    seq_aipe_R2(d$x[1:5, ], d$y[1:5], width = 0.9,
                ci_method = "log_transformed")
  )
  for (out in states) {
    expect_false(attr(out, "is_satisfied"))
    expect_gt(attr(out, "criterion_n"), attr(out, "current_n"))
  }
})

test_that("a width too narrow to give a sample size is refused with advice about the width", {
  # 4 z / 1e-9 = 7839855938.16 and, at R2 = .5 and n = 60,
  # (4 z / 1e-5)^2 (.5 + 1/60)^3 = 8.5e10: neither is a sample size that can
  # be represented (or collected).
  expect_error(
    seq_aipe_R2_pilot(1e-9, p = 3),
    paste0("The implied pilot sample size \\(7839855938\\) is too large to ",
           "represent\\. Check 'width': a target this narrow calls for ",
           "more observations than a study can collect\\.")
  )
  d <- exact_R2_data(60, 2, 0.5)
  expect_error(
    seq_aipe_R2(d$x, d$y, width = 1e-5),
    paste0("The implied criterion sample size .* is too large to ",
           "represent\\. Check 'width': a target this narrow")
  )
  expect_error(seq_aipe_R2_sensitivity(0.5, 1e-9, p = 3, G = 2),
               "pilot sample size .* too large to represent\\. Check 'width'")
  # A width ten times as large is inside the range: the criterion is
  # (4 z / 1e-4)^2 (.5 + 1/60)^3 = 8.48e8, an integer.
  out <- seq_aipe_R2(d$x, d$y, width = 1e-4)
  expect_identical(
    attr(out, "criterion_n"),
    as.integer(ceiling((4 * stats::qnorm(0.975) / 1e-4)^2 *
                         (0.5 + 1 / 60)^3))
  )
})

# ---------------------------------------------------------------------------
# The note on a sample that is small relative to the number of predictors
# ---------------------------------------------------------------------------

test_that("a sample is small relative to the predictors below p * max(20, 2p) observations", {
  # By hand: 20 per predictor up to 10 predictors (20, 40, 100, 200), then
  # 2 p per predictor (2 * 11^2 = 242, 2 * 20^2 = 800, 1,800, 3,200).
  expect_identical(
    vapply(c(1L, 2L, 5L, 10L, 11L, 20L, 30L, 40L), .R2_small_sample_below,
           numeric(1)),
    c(20, 40, 100, 200, 242, 800, 1800, 3200)
  )
  # No integer overflow for the largest number of predictors accepted.
  largest <- .Machine$integer.max - 2L
  expect_no_warning(
    expect_equal(.R2_small_sample_below(largest), 2 * as.numeric(largest)^2)
  )
})

test_that("the small sample note appears below the threshold and not from the threshold on", {
  for (p in c(1L, 2L, 5L, 10L, 11L, 20L)) {
    threshold <- as.integer(p * max(20L, 2L * p))
    d <- survey_data(threshold, p = p)
    under <- seq_len(threshold - 1L)
    for (ci_method in c("wald", "log_transformed")) {
      below <- seq_aipe_R2(d$x[under, , drop = FALSE], d$y[under],
                           width = 0.2, ci_method = ci_method)
      at <- seq_aipe_R2(d$x, d$y, width = 0.2, ci_method = ci_method)
      expect_length(attr(below, "notes"), 1L)
      expect_match(
        attr(below, "notes"),
        paste0("^With ", threshold - 1L, " observations and ", p,
               " predictors? the sample is small relative to the number of ",
               "predictors \\(fewer than ", threshold, " observations, ",
               max(20L, 2L * p), " per predictor\\)")
      )
      expect_identical(attr(at, "notes"), character(0))
      # The note is the only thing that differs in kind between the two.
      expect_identical(attr(below, "variance_unusable"), FALSE)
      expect_identical(attr(at, "variance_unusable"), FALSE)
    }
  }
})

test_that("the small sample note is worded in full and does not change the decision", {
  # R2 = .9 exactly with 5 predictors at n = 27, omega = .30, 95%: the Wald
  # criterion is (4 z / .3)^2 (.9 + 1/27) (.1 + 1/27)^2 = 12.02, so 13, and
  # the rule is met; 27 is below 100 = 5 * 20, so the note is attached.
  d <- exact_R2_data(27, 5, 0.9)
  out <- seq_aipe_R2(d$x, d$y, width = 0.30)
  expect_identical(attr(out, "criterion_n"), 13L)
  expect_true(attr(out, "is_satisfied"))
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_identical(
    attr(out, "notes"),
    paste0("With 27 observations and 5 predictors the sample is small ",
           "relative to the number of predictors (fewer than 100 ",
           "observations, 20 per predictor). The sample squared multiple ",
           "correlation then overestimates the population value, and in ",
           "simulations the interval of this procedure contained the ",
           "population value less often than the confidence level, far ",
           "less often the smaller the sample. This note does not change ",
           "the stopping rule. See 'Small Samples Relative to the Number ",
           "of Predictors' in the help of seq_aipe_R2().")
  )
  expect_output(print(out), "Stopping rule met: TRUE")
  expect_output(print(out), "Note: With 27 observations and 5 predictors")

  # R2 = .5 with the same n and p: (4 z / .3)^2 (.5 + 1/27)^3 = 105.78, so
  # 106: not met, and the note is attached all the same (it depends on n
  # and p alone).
  d <- exact_R2_data(27, 5, 0.5)
  out <- seq_aipe_R2(d$x, d$y, width = 0.30)
  expect_identical(attr(out, "criterion_n"), 106L)
  expect_false(attr(out, "is_satisfied"))
  expect_match(attr(out, "notes"), "^With 27 observations and 5 predictors")

  # One predictor is singular, and a large threshold is written with a
  # comma.
  expect_match(.note_R2_small_sample(12L, 1L),
               "^With 12 observations and 1 predictor the sample")
  expect_match(.note_R2_small_sample(100L, 30L),
               "fewer than 1,800 observations, 60 per predictor")
})

test_that("the numbers in the small sample note are written in full, whatever the options of the session", {
  # By hand: 500 predictors give 500 * max(20, 1000) = 500,000 observations,
  # 1,000 per predictor; 1,000 predictors give 1000 * 2000 = 2,000,000 and
  # 2,000. Left to the options of the session, a round number of 100,000 or
  # more is written in scientific notation (5e+05).
  expect_match(
    .note_R2_small_sample(600L, 500L),
    paste0("^With 600 observations and 500 predictors the sample is small ",
           "relative to the number of predictors \\(fewer than 500,000 ",
           "observations, 1,000 per predictor\\)\\. ")
  )
  expect_match(
    .note_R2_small_sample(250000L, 1000L),
    paste0("^With 250,000 observations and 1,000 predictors the sample is ",
           "small relative to the number of predictors \\(fewer than ",
           "2,000,000 observations, 2,000 per predictor\\)\\. ")
  )
  expect_false(grepl("e\\+", .note_R2_small_sample(100L, 50000L)))

  # The note is the same under options that change how numbers are written,
  # and no warning escapes: with OutDec = "," format() warns when it is left
  # to take the decimal mark from the option, and with a negative scipen it
  # writes 8e+02 for 800.
  d <- survey_data(30)
  usual_note <- .note_R2_small_sample(27L, 20L)
  expect_match(usual_note,
               "\\(fewer than 800 observations, 40 per predictor\\)")
  usual <- seq_aipe_R2(d$x, d$y, width = 0.3)
  expect_match(attr(usual, "notes"),
               paste0("^With 30 observations and 3 predictors .*\\(fewer ",
                      "than 60 observations, 20 per predictor\\)"))
  withr::local_options(OutDec = ",", scipen = -5)
  expect_no_warning(changed_note <- .note_R2_small_sample(27L, 20L))
  expect_identical(changed_note, usual_note)
  expect_no_warning(large <- .note_R2_small_sample(250000L, 1000L))
  expect_match(large, "^With 250,000 observations and 1,000 predictors ")
  expect_match(large,
               "fewer than 2,000,000 observations, 2,000 per predictor")
  expect_no_warning(changed <- seq_aipe_R2(d$x, d$y, width = 0.3))
  expect_identical(attr(changed, "notes"), attr(usual, "notes"))
  expect_identical(changed$value, usual$value)
})

test_that("twenty predictors unrelated to the response at the pilot sample size: met in many data sets, never silently", {
  # The case that prompted the note. The default pilot for a width of .30 is
  # 27 (4 z / omega = 26.13); with 20 pure noise predictors the sample R2 is
  # about 20/26 = .77 and the rule is met at once in more than a third of
  # data sets, with an interval nowhere near the population value of 0:
  # 3,703 and 3,623 of 10,000 data sets drawn as in the loop below with
  # seeds 1 and 2, a proportion of .366.
  expect_identical(seq_aipe_R2_pilot(0.30, p = 20), 27L)
  withr::local_seed(113)
  met <- 0L
  for (i in 1:90) {
    x <- matrix(stats::rnorm(27 * 20), 27, 20)
    y <- stats::rnorm(27)
    out <- seq_aipe_R2(x, y, width = 0.30)
    expect_match(attr(out, "notes"),
                 "^With 27 observations and 20 predictors .*fewer than 800")
    if (attr(out, "is_satisfied")) {
      met <- met + 1L
      expect_gt(row_value(out, "lower_limit"), 0)
    }
  }
  # .366 of 90 is 33.0 on average and four binomial standard errors are
  # 4 * sqrt(90 * .366 * .634) = 18.3, so at least 14.
  expect_gte(met, 14L)
})

test_that("at the threshold the fixed sample coverage of the Wald interval is about .90 (the calibration of the note)", {
  # The sample R2 is drawn from its exact distribution under multivariate
  # normality with random predictors, as Gatsonis and Sampson (1989, p. 517,
  # Equations 2 and 3) print it, which shares nothing with the package:
  #   R2 / (1 - R2) is distributed as chisq_{p + 2K} / chisq_{n - 1 - p},
  # the two chi square variables independent given K, and K negative
  # binomial with
  #   P(K = k) = Gamma(k + (n - 1)/2) / [k! Gamma((n - 1)/2)]
  #              P2^k (1 - P2)^((n - 1)/2),
  # which is stats::rnbinom() with size (n - 1)/2 and prob 1 - P2.
  # The interval is the package's (.R2_rule()). The expected values are the
  # coverage in 400,000 draws per condition of the same distribution written
  # in another form, conditionally on the predictors, with the two intervals
  # coded from their formulas and not from the package:
  #   R2 / (1 - R2) = ((Z + sqrt(phi W))^2 + U) / V,  phi = P2 / (1 - P2),
  # with Z standard normal and W, U, and V chi square on n - 1, p - 1, and
  # n - p - 1 degrees of freedom, all independent. Another 400,000 draws of
  # that form per condition reproduce the nine values within .0014, and
  # 400,000 draws of draw_R2() below within .0021 (at most 2.5 and 2.8
  # standard errors of the difference). Tolerance: four standard errors of
  # the difference. This run has 20,000 draws, so its standard error is at
  # most sqrt(.21 * .79 / 20000) = .0029, and that of the expected value at
  # most .0007; four times the root sum of squares is .012.
  draw_R2 <- function(G, n, p, P2) {
    K <- stats::rnbinom(G, size = (n - 1) / 2, prob = 1 - P2)
    ratio <- stats::rchisq(G, p + 2 * K) / stats::rchisq(G, n - 1 - p)
    ratio / (1 + ratio)
  }
  coverage <- function(n, p, P2, ci_method = "wald") {
    R2 <- draw_R2(20000, n, p, P2)
    rule <- .R2_rule(R2, 1 - R2, n, p, width = 0.2, z = z95,
                     ci_method = ci_method)
    mean(rule$lower <= P2 & P2 <= rule$upper)
  }
  withr::local_seed(113)
  # At the threshold, n = p * max(20, 2p):
  expect_lt(abs(coverage(40, 2, 0.7) - 0.9084), 0.012)
  expect_lt(abs(coverage(100, 5, 0.5) - 0.9145), 0.012)
  expect_lt(abs(coverage(800, 20, 0.3) - 0.9049), 0.012)
  expect_lt(abs(coverage(3200, 40, 0.3) - 0.9040), 0.012)
  # Half and a quarter of it: lower. Twice: higher.
  expect_lt(abs(coverage(100, 10, 0.5) - 0.8530), 0.012)
  expect_lt(abs(coverage(200, 20, 0.5) - 0.7922), 0.012)
  expect_lt(abs(coverage(400, 10, 0.5) - 0.9288), 0.012)
  # The log transformed interval at the threshold and at a quarter of it.
  expect_lt(abs(coverage(800, 20, 0.3, "log_transformed") - 0.9208), 0.012)
  expect_lt(abs(coverage(200, 20, 0.5, "log_transformed") - 0.8701), 0.012)
})

# ---------------------------------------------------------------------------
# Sensitivity
# ---------------------------------------------------------------------------

sensitivity_value <- function(out, term) out$value[out$term == term]

test_that("seq_aipe_R2_sensitivity returns the summary rows of the AIPE sensitivity family", {
  out <- seq_aipe_R2_sensitivity(true_value = 0.4, width = 0.3, p = 3, G = 40,
                                 seed = 113)
  # The class of every simulation summary in the package, printed by its
  # shared method; underneath, a term and value data frame.
  expect_s3_class(out, c("smsd_sensitivity", "data.frame"), exact = TRUE)
  expect_identical(capture.output(print(out)), format(out))
  expect_length(format(out), 14L)
  expect_identical(names(out), c("term", "value"))
  expect_type(out$value, "double")
  expect_identical(out$term,
                   c("true_value", "target_width", "conf_level",
                     "replications", "mean_n", "median_n", "sd_n",
                     "mean_width", "median_width", "pct_width_within",
                     "coverage", "pct_hit_max_n", "pct_variance_unusable"))
  expect_identical(sensitivity_value(out, "true_value"), 0.4)
  expect_identical(sensitivity_value(out, "target_width"), 0.3)
  expect_identical(sensitivity_value(out, "conf_level"), 0.95)
  expect_identical(sensitivity_value(out, "replications"), 40)
  expect_true(all(is.finite(out$value)))
  # The width guarantee holds in every replication; no cap; continuous data.
  expect_identical(sensitivity_value(out, "pct_width_within"), 100)
  expect_identical(sensitivity_value(out, "pct_hit_max_n"), 0)
  expect_identical(sensitivity_value(out, "pct_variance_unusable"), 0)
  # No replication stops before the pilot sample size.
  expect_gte(sensitivity_value(out, "mean_n"), seq_aipe_R2_pilot(0.3, 3))
})

test_that("a seeded call reproduces and leaves the caller's random number stream alone", {
  a <- seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 25, seed = 113)
  b <- seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 25, seed = 113)
  expect_identical(a, b)
  other <- seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 25, seed = 114)
  expect_false(identical(a, other))

  # The caller's stream is seeded with a value other than 113 (42) so that
  # it cannot coincide with the seed passed to the function.
  untouched <- withr::with_seed(42, stats::runif(3))
  after_call <- withr::with_seed(42, {
    seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 5, seed = 113)
    stats::runif(3)
  })
  expect_identical(after_call, untouched)

  # seed = NULL draws from the caller's stream (again seeded with a value
  # other than 113, so the result must differ from the seeded call above).
  first <- withr::with_seed(42,
                            seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 25))
  second <- withr::with_seed(42,
                             seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 25))
  expect_identical(first, second)
  expect_false(identical(first, a))
})

test_that("seq_aipe_R2_sensitivity checks its arguments", {
  expect_error(seq_aipe_R2_sensitivity(width = 0.3, p = 2),
               "'true_value' must be supplied")
  expect_error(seq_aipe_R2_sensitivity(1, 0.3, p = 2), "'true_value'")
  expect_error(seq_aipe_R2_sensitivity(-0.1, 0.3, p = 2), "'true_value'")
  expect_error(seq_aipe_R2_sensitivity(NA, 0.3, p = 2), "'true_value'")
  expect_error(seq_aipe_R2_sensitivity(0.4, p = 2), "'width' must be supplied")
  expect_error(seq_aipe_R2_sensitivity(0.4, 1.5, p = 2), "'width'")
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3), "'p' must be supplied")
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 0), "'p'")
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, conf_level = 2),
               "'conf_level'")
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, ci_method = "x"),
               paste0("'ci_method' must be one of \"wald\", ",
                      "\"log_transformed\"; got \"x\""))
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 0), "'G'")
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 10.5), "'G'")
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 5, seed = 1.5),
               "'seed'")
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 5, seed = "a"),
               "'seed'")
  # A seed too large for set.seed() is refused in the package's own words,
  # with no warning about coercion.
  expect_no_warning(
    expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 5, seed = 1e10),
                 "'seed' must be NULL or a single whole number.*got 1e\\+10")
  )
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 5, seed = NA),
               "'seed' must be NULL or a single whole number")
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 5, seed = 1:2),
               "'seed' must be NULL or a single whole number")
  # The seed is checked with the other arguments, one at a time, before the
  # checks that combine arguments (here max_n against the pilot sample size
  # of 27), so a bad seed is not masked by them.
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 5, seed = 1.5,
                                       max_n = 20),
               "'seed' must be NULL or a single whole number")
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, m_0 = 3), "'m_0'")
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, batch = 0), "'batch'")
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, batch = 2.5),
               "'batch'")
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, max_n = 0), "'max_n'")
  # The pilot is 27 here (4 z / omega = 26.13).
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, max_n = 20),
               "'max_n' \\(20\\) is below the pilot sample size \\(27\\)")
  expect_error(seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, max_n = 4,
                                       ci_method = "log_transformed"),
               "'max_n' must be at least p \\+ 3")
})

test_that("max_n is an exact cap and batch sets the spacing of the checks", {
  capped <- seq_aipe_R2_sensitivity(0.5, 0.05, p = 2, G = 20, seed = 113,
                                    max_n = 300)
  expect_identical(sensitivity_value(capped, "mean_n"), 300)
  expect_identical(sensitivity_value(capped, "sd_n"), 0)
  expect_identical(sensitivity_value(capped, "pct_hit_max_n"), 100)
  expect_identical(sensitivity_value(capped, "pct_width_within"), 0)

  # The cap is exact when it falls inside a batch: the pilot is 53
  # (4 z / omega = 52.27), the next check would be at 103, and the rule
  # needs about 340 observations here, so every replication stops at 80.
  inside <- seq_aipe_R2_sensitivity(0.5, 0.15, p = 2, G = 20, seed = 113,
                                    max_n = 80, batch = 50)
  expect_identical(sensitivity_value(inside, "mean_n"), 80)
  expect_identical(sensitivity_value(inside, "pct_hit_max_n"), 100)

  # A cap equal to the pilot sample size stops every replication there.
  at_pilot <- seq_aipe_R2_sensitivity(0.5, 0.3, p = 2, G = 20, seed = 113,
                                      max_n = 27)
  expect_identical(sensitivity_value(at_pilot, "mean_n"), 27)

  # With batch = 5 every realized sample size is the pilot plus a multiple
  # of 5.
  withr::local_seed(113)
  draw <- function(m) {
    new <- matrix(stats::rnorm(m * 3), m, 3)
    new[, 3] <- 0.5 * (new[, 1] + new[, 2]) + sqrt(0.5) * new[, 3]
    new
  }
  realized <- vapply(1:60, function(i) {
    .R2_sensitivity_one(draw, pilot_n = 27L, p = 2L, width = 0.3, z = z95,
                        ci_method = "wald", max_n = 1e5L, batch = 5L)[["n"]]
  }, numeric(1))
  expect_true(all((realized - 27) %% 5 == 0))
  expect_gt(length(unique(realized)), 3L)

  # Checking less often cannot stop earlier, so the mean sample size grows.
  one <- seq_aipe_R2_sensitivity(0.5, 0.3, p = 2, G = 300, seed = 113)
  twenty <- seq_aipe_R2_sensitivity(0.5, 0.3, p = 2, G = 300, seed = 113,
                                    batch = 20)
  expect_gt(sensitivity_value(twenty, "mean_n"),
            sensitivity_value(one, "mean_n"))
})

test_that("a negative seed is accepted, as set.seed() accepts it", {
  a <- seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 10, seed = -7)
  b <- seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 10, seed = -7)
  expect_identical(a, b)
  expect_false(identical(a, seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 10,
                                                    seed = 7)))
})

test_that("m_0 sets the first check in seq_aipe_R2_sensitivity", {
  # Wald, P2 = .5, omega = .3, 95%: the rule needs about 85 observations
  # (16 z^2 (.5)(.25) / .09 = 85.4) and the default pilot is 27
  # (4 z / omega = 26.13), so with the default the replications stop near 85,
  # and with a pilot of 200 or 500 every one of them stops at the pilot.
  default <- seq_aipe_R2_sensitivity(0.5, 0.3, p = 2, G = 30, seed = 113)
  expect_gt(sensitivity_value(default, "mean_n"), 60)
  expect_lt(sensitivity_value(default, "mean_n"), 110)
  for (m_0 in c(200, 500)) {
    out <- seq_aipe_R2_sensitivity(0.5, 0.3, p = 2, G = 30, seed = 113,
                                   m_0 = m_0)
    expect_identical(sensitivity_value(out, "mean_n"), m_0)
    expect_identical(sensitivity_value(out, "median_n"), m_0)
    expect_identical(sensitivity_value(out, "sd_n"), 0)
  }
  # Log transformed, P2 = .9, omega = .3: Equation 62 at R2 = .9 is 13.7 and
  # the default pilot is 4, so a pilot of 60 is where every replication
  # stops, and without it most stop far earlier.
  out <- seq_aipe_R2_sensitivity(0.9, 0.3, p = 2, G = 30, seed = 113,
                                 m_0 = 60, ci_method = "log_transformed")
  expect_identical(sensitivity_value(out, "mean_n"), 60)
  default <- seq_aipe_R2_sensitivity(0.9, 0.3, p = 2, G = 30, seed = 113,
                                     ci_method = "log_transformed")
  expect_lt(sensitivity_value(default, "mean_n"), 30)
})

test_that("the summary rows are the mean, median, and standard deviation of the replications", {
  # The same 60 replications, produced one at a time with the generator the
  # help page describes and summarized here in base R.
  p <- 2L
  draw <- function(m) {
    obs <- matrix(stats::rnorm(m * (p + 1L)), nrow = m, ncol = p + 1L,
                  byrow = TRUE)
    obs[, p + 1L] <- sqrt(0.5 / p) * rowSums(obs[, seq_len(p), drop = FALSE]) +
      sqrt(0.5) * obs[, p + 1L]
    obs
  }
  runs <- withr::with_seed(113, vapply(1:60, function(g) {
    .R2_sensitivity_one(draw, 27L, p, 0.3, z95, "wald", 1e5L, 1L)
  }, numeric(5)))
  out <- seq_aipe_R2_sensitivity(0.5, 0.3, p = 2, G = 60, seed = 113)
  widths <- runs["upper", ] - runs["lower", ]
  expect_gt(stats::sd(runs["n", ]), 1)
  expect_equal(sensitivity_value(out, "mean_n"), mean(runs["n", ]))
  expect_equal(sensitivity_value(out, "median_n"), stats::median(runs["n", ]))
  expect_equal(sensitivity_value(out, "sd_n"), stats::sd(runs["n", ]))
  expect_false(isTRUE(all.equal(stats::sd(runs["n", ]),
                                stats::var(runs["n", ]))))
  expect_equal(sensitivity_value(out, "mean_width"), mean(widths))
  expect_equal(sensitivity_value(out, "median_width"), stats::median(widths))
  expect_false(isTRUE(all.equal(stats::median(widths), mean(widths))))
  covered <- runs["lower", ] <= 0.5 & 0.5 <= runs["upper", ]
  expect_gt(sum(!covered), 0)
  expect_identical(sensitivity_value(out, "coverage"), sum(covered) / 60)
  expect_identical(sensitivity_value(out, "pct_width_within"),
                   100 * sum(widths <= 0.3) / 60)
})

test_that("the summary is formed from counts and is not broken by a replication without an interval", {
  # Five replications written out by hand, for a population value of .50 and
  # a target width of .20. The first three stop by the rule; the fourth is
  # stopped by the cap with an interval wider than the target; the fifth is
  # stopped by the cap with no interval at all (limits of NA).
  runs <- rbind(n       = c(100, 120, 90, 300, 300),
                lower   = c(0.42, 0.51, 0.40, 0.30, NA),
                upper   = c(0.60, 0.69, 0.59, 0.62, NA),
                capped  = c(0, 0, 0, 1, 1),
                unusable = c(0, 1, 0, 0, 1))
  out <- .R2_sensitivity_summary(runs, true_value = 0.5, width = 0.2,
                                 conf_level = 0.9)
  expect_identical(out$term,
                   c("true_value", "target_width", "conf_level",
                     "replications", "mean_n", "median_n", "sd_n",
                     "mean_width", "median_width", "pct_width_within",
                     "coverage", "pct_hit_max_n", "pct_variance_unusable"))
  expect_identical(sensitivity_value(out, "true_value"), 0.5)
  expect_identical(sensitivity_value(out, "target_width"), 0.2)
  expect_identical(sensitivity_value(out, "conf_level"), 0.9)
  expect_identical(sensitivity_value(out, "replications"), 5)
  expect_identical(sensitivity_value(out, "mean_n"), 182)
  expect_identical(sensitivity_value(out, "median_n"), 120)
  # Deviations from 182: -82, -62, -92, 118, 118; squares sum to 46880, so
  # the standard deviation is sqrt(46880 / 4) = 108.2589488.
  expect_equal(sensitivity_value(out, "sd_n"), 108.2589488, tolerance = 1e-9)
  # Widths .18, .18, .19, .32 and none: mean .2175, median .185.
  expect_equal(sensitivity_value(out, "mean_width"), 0.2175, tolerance = 1e-12)
  expect_equal(sensitivity_value(out, "median_width"), 0.185,
               tolerance = 1e-12)
  # Three of the five are within the target; the one without an interval is
  # not.
  expect_identical(sensitivity_value(out, "pct_width_within"), 60)
  # A width equal to the target is within it. (Limits that are exact in
  # binary, so that the widths are exactly .25 and .375.)
  at_target <- rbind(n = c(50, 60), lower = c(0.25, 0.375),
                     upper = c(0.5, 0.75), capped = c(0, 0),
                     unusable = c(0, 0))
  expect_identical(
    sensitivity_value(.R2_sensitivity_summary(at_target, true_value = 0.4,
                                              width = 0.25,
                                              conf_level = 0.9),
                      "pct_width_within"),
    50
  )
  # .50 lies in the first, third, and fourth intervals; not in the second;
  # the fifth has no interval and does not cover.
  expect_identical(sensitivity_value(out, "coverage"), 0.6)
  expect_identical(sensitivity_value(out, "pct_hit_max_n"), 40)
  expect_identical(sensitivity_value(out, "pct_variance_unusable"), 40)
  expect_true(all(is.finite(out$value)))

  # No replication has an interval: the width summaries are NA, the rest are
  # not.
  none <- runs[, c(5, 5), drop = FALSE]
  out <- .R2_sensitivity_summary(none, true_value = 0.5, width = 0.2,
                                 conf_level = 0.9)
  expect_true(is.na(sensitivity_value(out, "mean_width")))
  expect_true(is.na(sensitivity_value(out, "median_width")))
  expect_identical(sensitivity_value(out, "coverage"), 0)
  expect_identical(sensitivity_value(out, "pct_width_within"), 0)
  expect_identical(sensitivity_value(out, "pct_hit_max_n"), 100)
})

test_that("one replication is allowed and has no standard deviation", {
  out <- seq_aipe_R2_sensitivity(0.4, 0.3, p = 2, G = 1, seed = 113)
  expect_identical(sensitivity_value(out, "replications"), 1)
  expect_true(is.na(sensitivity_value(out, "sd_n")))
  expect_true(all(is.finite(out$value[out$term != "sd_n"])))
  expect_identical(sensitivity_value(out, "mean_n"),
                   sensitivity_value(out, "median_n"))
})

# Serve the rows of a fixed data matrix as the engine asks for them.
serve_rows <- function(data) {
  served <- new.env()
  served$rows <- 0L
  function(m) {
    rows <- served$rows + seq_len(m)
    served$rows <- served$rows + m
    data[rows, , drop = FALSE]
  }
}

test_that("the engine draws a pilot sample larger than one chunk in pieces", {
  # The matrix of running sums is kept near 2e6 cells: with 100 predictors
  # (5,151 sums of products) a chunk is floor(2e6 / 5151) = 388 rows, so a
  # pilot sample of 400 is drawn in two pieces before the first check. The
  # result is the one a loop of seq_aipe_R2() calls gives on the same rows.
  p <- 100L
  expect_identical(floor(2e6 / ((p + 1) * (p + 2) / 2)), 388)
  withr::local_seed(113)
  x <- matrix(stats::rnorm(450 * p), 450, p)
  y <- 0.5 * rowSums(x[, 1:3]) + stats::rnorm(450)
  one <- .R2_sensitivity_one(serve_rows(cbind(x, y)), pilot_n = 400L, p = p,
                             width = 0.3, z = z95, ci_method = "wald",
                             max_n = 450L, batch = 1L)
  n <- 400L
  repeat {
    direct <- seq_aipe_R2(x[seq_len(n), ], y[seq_len(n)], width = 0.3)
    if (attr(direct, "is_satisfied") || n == 450L) break
    n <- n + 1L
  }
  expect_identical(one[["n"]], as.numeric(n))
  expect_identical(one[["capped"]], 0)
  expect_equal(one[["lower"]], row_value(direct, "lower_limit"),
               tolerance = 1e-10)
  expect_equal(one[["upper"]], row_value(direct, "upper_limit"),
               tolerance = 1e-10)
})

test_that("collinear predictors are a check that is not met in the engine, and an error in the check", {
  # In a simulated study a sample whose predictors are collinear is a check
  # at which the rule is not met (R2 is NaN), and sampling goes on; the
  # check function stops with an error that says what to do.
  withr::local_seed(113)
  x <- matrix(stats::rnorm(80), 40, 2)
  x[1:6, 2] <- 2 * x[1:6, 1]
  y <- drop(x %*% c(0.5, 0.5)) + stats::rnorm(40)
  expect_error(seq_aipe_R2(x[1:6, ], y[1:6], width = 0.9),
               "The predictors in 'x' are collinear")
  at_six <- .R2_estimate(x[1:6, ], y[1:6], stop_on_degenerate = FALSE)
  expect_true(is.nan(at_six$R2))
  expect_true(is.nan(at_six$one_minus_R2))
  one <- .R2_sensitivity_one(serve_rows(cbind(x, y)), pilot_n = 6L, p = 2L,
                             width = 0.9, z = z95, ci_method = "wald",
                             max_n = 40L, batch = 1L)
  expect_gt(one[["n"]], 6)
  n <- one[["n"]]
  expect_true(attr(seq_aipe_R2(x[seq_len(n), ], y[seq_len(n)], width = 0.9),
                   "is_satisfied"))
  # A matrix with no predictor column at all.
  expect_error(seq_aipe_R2(matrix(numeric(0), 10, 0), stats::rnorm(10),
                           width = 0.3),
               "'x' must have at least one predictor column")
})

test_that("a replication that meets the rule exactly at max_n is not counted as stopped by the cap", {
  # The served data have R2 = .5 exactly at n = 198, where the Wald rule is
  # first met (see the boundary test above); the only checks are at the
  # pilot (40) and at max_n.
  d <- exact_R2_data(198, 2, 0.5)
  met <- .R2_sensitivity_one(serve_rows(cbind(d$x, d$y)), pilot_n = 40L,
                             p = 2L, width = 0.20, z = z95,
                             ci_method = "wald", max_n = 198L, batch = 1000L)
  expect_identical(met[["n"]], 198)
  expect_identical(met[["capped"]], 0)
  # One observation fewer and the rule is not met there (197.98 > 197): the
  # cap stopped the replication.
  d <- exact_R2_data(197, 2, 0.5)
  capped <- .R2_sensitivity_one(serve_rows(cbind(d$x, d$y)), pilot_n = 40L,
                                p = 2L, width = 0.20, z = z95,
                                ci_method = "wald", max_n = 197L,
                                batch = 1000L)
  expect_identical(capped[["n"]], 197)
  expect_identical(capped[["capped"]], 1)
})

test_that("the engine treats an unusable variance estimate as seq_aipe_R2() does", {
  # The first 30 responses are an exact linear function of the predictors
  # (sample R2 = 1 at the pilot check); error enters afterwards. As printed
  # the rule is met at the pilot with an interval of zero width at 1; with
  # the package's rule the pilot check is "not met", and the rule is first met at
  # n = 31, where seq_aipe_R2() on the same rows reports it met.
  withr::local_seed(113)
  x <- matrix(stats::rnorm(600), 300, 2)
  y <- drop(x %*% c(1, -1)) + c(rep(0, 30), stats::rnorm(270))
  data <- cbind(x, y)
  for (ci_method in c("wald", "log_transformed")) {
    one <- .R2_sensitivity_one(serve_rows(data), pilot_n = 30L, p = 2L,
                               width = 0.3, z = z95, ci_method = ci_method,
                               max_n = 300L, batch = 1L)
    at_pilot <- seq_aipe_R2(x[1:30, ], y[1:30], 0.3, ci_method = ci_method)
    expect_identical(attr(at_pilot, "variance_unusable"), TRUE)
    expect_identical(attr(at_pilot, "is_satisfied"), FALSE)
    at_stop <- seq_aipe_R2(x[1:31, ], y[1:31], 0.3, ci_method = ci_method)
    expect_identical(attr(at_stop, "is_satisfied"), TRUE)
    expect_identical(one[["n"]], 31)
    expect_identical(one[["unusable"]], 1)
    expect_identical(one[["capped"]], 0)
    expect_equal(one[["lower"]], row_value(at_stop, "lower_limit"),
                 tolerance = 1e-10)
    expect_equal(one[["upper"]], row_value(at_stop, "upper_limit"),
                 tolerance = 1e-10)
  }

  # The same with an offset of 1e10 on every variable, where the residual
  # share of the perfect fit is rounding error (about 1e-12) rather than
  # zero: the pilot check is still "not met", as it is in seq_aipe_R2().
  # (The slope is not a power of two, so the products are rounded.)
  shifted <- cbind(1e10 + x[, 1], 0.3 * (1e10 + x[, 1]) + 1)
  shifted[31:300, 2] <- shifted[31:300, 2] + stats::rnorm(270)
  at_pilot <- seq_aipe_R2(shifted[1:30, 1], shifted[1:30, 2], 0.3)
  expect_identical(attr(at_pilot, "variance_unusable"), TRUE)
  expect_gt(1 - row_value(at_pilot, "R2"), .Machine$double.eps)
  one <- .R2_sensitivity_one(serve_rows(shifted), pilot_n = 30L, p = 1L,
                             width = 0.3, z = z95, ci_method = "wald",
                             max_n = 300L, batch = 1L)
  expect_gt(one[["n"]], 30)
  expect_identical(one[["unusable"]], 1)

  # A perfect fit throughout: never met, stopped by the cap with no interval,
  # counted as unusable.
  perfect <- cbind(x[1:60, ], 3 + x[1:60, 1] - 2 * x[1:60, 2])
  one <- .R2_sensitivity_one(serve_rows(perfect), pilot_n = 30L, p = 2L,
                             width = 0.3, z = z95, ci_method = "wald",
                             max_n = 60L, batch = 1L)
  expect_identical(one[["n"]], 60)
  expect_identical(one[["capped"]], 1)
  expect_identical(one[["unusable"]], 1)
  expect_true(is.na(one[["lower"]]))
  expect_true(is.na(one[["upper"]]))
  at_cap <- seq_aipe_R2(perfect[, 1:2], perfect[, 3], 0.3)
  expect_true(is.na(row_value(at_cap, "lower_limit")))

  # A response that does not vary: no estimate, so no limits, stopped by the
  # cap.
  flat <- cbind(x[1:60, ], 0.1)
  one <- .R2_sensitivity_one(serve_rows(flat), pilot_n = 30L, p = 2L,
                             width = 0.3, z = z95, ci_method = "wald",
                             max_n = 60L, batch = 1L)
  expect_identical(one[["n"]], 60)
  expect_true(is.na(one[["lower"]]))
  expect_true(is.na(one[["upper"]]))
  expect_identical(one[["capped"]], 1)
  expect_identical(one[["unusable"]], 1)

  # A sample on which seq_aipe_R2() stops with an error (here the second
  # predictor does not vary in the first 35 rows) is a check that is not
  # met, and sampling goes on.
  stuck <- cbind(x[, 1], c(rep(1, 35), x[36:300, 2]),
                 drop(x %*% c(1, -1)) + stats::rnorm(300))
  expect_error(seq_aipe_R2(stuck[1:35, 1:2], stuck[1:35, 3], 0.3),
               "does not vary")
  expect_s3_class(seq_aipe_R2(stuck[1:36, 1:2], stuck[1:36, 3], 0.3),
                  "smsd_check")
  one <- .R2_sensitivity_one(serve_rows(stuck), pilot_n = 30L, p = 2L,
                             width = 0.3, z = z95, ci_method = "wald",
                             max_n = 300L, batch = 1L)
  expect_gt(one[["n"]], 35)
  expect_identical(one[["unusable"]], 1)
  reference <- 36L
  repeat {
    out <- seq_aipe_R2(stuck[seq_len(reference), 1:2],
                       stuck[seq_len(reference), 3], 0.3)
    if (attr(out, "is_satisfied") || reference >= 300L) break
    reference <- reference + 1L
  }
  expect_identical(one[["n"]], as.numeric(reference))
  expect_equal(one[["lower"]], row_value(out, "lower_limit"),
               tolerance = 1e-8)

  # Ordinary continuous data are never counted as unusable.
  ordinary <- cbind(x, drop(x %*% c(1, -1)) + stats::rnorm(300))
  one <- .R2_sensitivity_one(serve_rows(ordinary), pilot_n = 30L, p = 2L,
                             width = 0.3, z = z95, ci_method = "wald",
                             max_n = 300L, batch = 1L)
  expect_identical(one[["unusable"]], 0)
  expect_identical(one[["capped"]], 0)
})

test_that("the sensitivity engine gives what a loop of seq_aipe_R2() calls gives on the same observations", {
  withr::local_seed(113)
  loop <- function(data, pilot_n, p, width, conf_level, ci_method, max_n,
                   batch) {
    n <- pilot_n
    repeat {
      out <- seq_aipe_R2(data[seq_len(n), seq_len(p), drop = FALSE],
                         data[seq_len(n), p + 1L], width, conf_level,
                         ci_method)
      if (attr(out, "is_satisfied") || n >= max_n) break
      n <- min(n + batch, max_n)
    }
    c(n = n, lower = row_value(out, "lower_limit"),
      upper = row_value(out, "upper_limit"),
      capped = as.numeric(!attr(out, "is_satisfied")))
  }
  compared <- 0L
  for (i in 1:24) {
    p <- sample(c(1L, 2L, 3L, 6L), 1)
    P2 <- stats::runif(1, 0, 0.9)
    width <- stats::runif(1, 0.2, 0.5)
    conf_level <- sample(c(0.90, 0.95), 1)
    ci_method <- sample(c("wald", "log_transformed"), 1)
    batch <- sample(c(1L, 3L, 10L), 1)
    max_n <- sample(c(60L, 400L), 1)
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    pilot_n <- seq_aipe_R2_pilot(width, p, conf_level, ci_method)
    if (max_n < pilot_n) next
    data <- matrix(stats::rnorm(400 * (p + 1L)), 400, p + 1L)
    data[, p + 1L] <- sqrt(P2 / p) * rowSums(data[, seq_len(p), drop = FALSE]) +
      sqrt(1 - P2) * data[, p + 1L]
    data <- sweep(data, 2, stats::runif(p + 1L, -20, 20), `+`)
    engine <- .R2_sensitivity_one(serve_rows(data), pilot_n, p, width, z,
                                  ci_method, max_n, batch)
    reference <- loop(data, pilot_n, p, width, conf_level, ci_method, max_n,
                      batch)
    expect_identical(engine[["n"]], reference[["n"]])
    expect_identical(engine[["capped"]], reference[["capped"]])
    expect_equal(engine[["lower"]], reference[["lower"]], tolerance = 1e-8)
    expect_equal(engine[["upper"]], reference[["upper"]], tolerance = 1e-8)
    compared <- compared + 1L
  }
  expect_gt(compared, 18L)

  # The same holds when the variables carry an offset large enough to defeat
  # running sums (the engine then computes from the observations themselves).
  for (offset in c(1e3, 1e6, 1e9)) {
    data <- matrix(stats::rnorm(300 * 3L), 300, 3L)
    data[, 3L] <- 0.5 * (data[, 1L] + data[, 2L]) + sqrt(0.5) * data[, 3L]
    data <- data + offset
    engine <- .R2_sensitivity_one(serve_rows(data), 27L, 2L, 0.3, z95,
                                  "wald", 300L, 1L)
    reference <- loop(data, 27L, 2L, 0.3, 0.95, "wald", 300L, 1L)
    expect_identical(engine[["n"]], reference[["n"]])
    expect_identical(engine[["capped"]], reference[["capped"]])
    expect_equal(engine[["lower"]], reference[["lower"]], tolerance = 1e-8)
    expect_equal(engine[["upper"]], reference[["upper"]], tolerance = 1e-8)
  }
})

# ---------------------------------------------------------------------------
# Operating characteristics against the article's tables
# ---------------------------------------------------------------------------

# Each article value is the mean of 5,000 replications. For the mean final
# sample size and the coverage the article prints the standard error, and
# the tolerance is four standard errors of the difference between its value
# and this simulation (both Monte Carlo errors combined). The article prints
# no standard error for the mean interval width; its tolerance is four
# standard errors of the difference computed from the standard deviation of
# the realized width (the width of the final interval) in 20,000
# replications of the same design, stated in each test, plus .00005 for the
# rounding of the printed value. So the tests cannot flake.

test_that("the Wald procedure reproduces Table 13 (90%, k = 2, omega = .10, P2 = .90)", {
  # Article: mean N 68.66 (se 0.0963), n_OF 39, ratio 1.7605, coverage .8678
  # (se .0048), mean width .0733. Equation 63 as printed gives a mean N of
  # about 66.4 here, more than 20 standard errors away.
  G <- 4000
  out <- seq_aipe_R2_sensitivity(0.90, width = 0.10, p = 2, conf_level = 0.90,
                                 G = G, seed = 113)
  mean_n <- sensitivity_value(out, "mean_n")
  se_n <- sqrt(0.0963^2 + sensitivity_value(out, "sd_n")^2 / G)   # about 0.14
  expect_lt(abs(mean_n - 68.66), 4 * se_n)
  expect_lt(abs(mean_n / 39 - 1.7605), 4 * se_n / 39)
  coverage <- sensitivity_value(out, "coverage")
  se_p <- sqrt(0.0048^2 + 0.8678 * (1 - 0.8678) / G)              # about .0072
  expect_lt(abs(coverage - 0.8678), 4 * se_p)
  # Standard deviation of the realized width .01274, so four standard errors
  # of the difference are 4 * .01274 * sqrt(1/4000 + 1/5000) = .00108.
  expect_lt(abs(sensitivity_value(out, "mean_width") - 0.0733), 0.0012)
  expect_identical(sensitivity_value(out, "median_n"), 66)   # the pilot
  expect_identical(sensitivity_value(out, "pct_width_within"), 100)
})

test_that("the log transformed procedure reproduces Table 11 (90%, k = 2, omega = .10, P2 = .90)", {
  # Article: mean N 34.84 (se 0.3222), n_BW 47, ratio .7412, coverage .7688
  # (se .0060), mean width .0908.
  G <- 4000
  out <- seq_aipe_R2_sensitivity(0.90, width = 0.10, p = 2, conf_level = 0.90,
                                 ci_method = "log_transformed", G = G,
                                 seed = 113)
  mean_n <- sensitivity_value(out, "mean_n")
  se_n <- sqrt(0.3222^2 + sensitivity_value(out, "sd_n")^2 / G)   # about 0.49
  expect_lt(abs(mean_n - 34.84), 4 * se_n)
  expect_lt(abs(mean_n / 47 - 0.7412), 4 * se_n / 47)
  coverage <- sensitivity_value(out, "coverage")
  se_p <- sqrt(0.0060^2 + 0.7688 * (1 - 0.7688) / G)              # about .0090
  expect_lt(abs(coverage - 0.7688), 4 * se_p)
  # Standard deviation of the realized width .01664, so four standard errors
  # of the difference are 4 * .01664 * sqrt(1/4000 + 1/5000) = .00141.
  expect_lt(abs(sensitivity_value(out, "mean_width") - 0.0908), 0.0015)
})

test_that("the Wald procedure reproduces Table 14 (95%, k = 5, omega = .10, P2 = .50)", {
  skip_on_cran()
  # Article: mean N 766.17 (se 0.5611), n_OF 769, ratio .9963, coverage .9450
  # (se .0032), mean width .0995.
  G <- 600
  out <- seq_aipe_R2_sensitivity(0.50, width = 0.10, p = 5, conf_level = 0.95,
                                 G = G, seed = 113)
  mean_n <- sensitivity_value(out, "mean_n")
  se_n <- sqrt(0.5611^2 + sensitivity_value(out, "sd_n")^2 / G)   # about 1.7
  expect_lt(abs(mean_n - 766.17), 4 * se_n)
  expect_lt(abs(mean_n / 769 - 0.9963), 4 * se_n / 769)
  coverage <- sensitivity_value(out, "coverage")
  se_p <- sqrt(0.0032^2 + 0.9450 * (1 - 0.9450) / G)              # about .0098
  expect_lt(abs(coverage - 0.9450), 4 * se_p)
  # Standard deviation of the realized width .00009 (every replication stops
  # just inside the target), so four standard errors of the difference are
  # .00002; with the rounding of the printed value, .0001.
  expect_lt(abs(sensitivity_value(out, "mean_width") - 0.0995), 0.0001)
})

test_that("where the note is attached at every stop, the procedure's coverage is far below the confidence level", {
  # Ten predictors, a population value of .50, a target width of .30, 95%:
  # coverage .6658 in seq_aipe_R2_sensitivity(0.50, width = 0.30, p = 10,
  # G = 10000, seed = 1), with a mean final sample size of 71. The rule
  # calls for about 85 observations at the population value, far below
  # 10 * 20 = 200, so every replication stops with the note attached.
  # Tolerance: four standard errors of the difference between 400 and 10,000
  # replications, 4 * sqrt(.6658 * .3342 * (1/400 + 1/10000)) = .096.
  out <- seq_aipe_R2_sensitivity(0.50, width = 0.30, p = 10, G = 400,
                                 seed = 113)
  expect_lt(abs(sensitivity_value(out, "coverage") - 0.6658), 0.096)
  expect_lt(sensitivity_value(out, "mean_n") +
              4 * sensitivity_value(out, "sd_n"), 200)
  # The same design with a pilot sample at the threshold, where the note
  # goes away: the rule is already met there, so every replication stops at
  # 200, and coverage is that of the interval at a fixed sample size of 200:
  # .9046 in 400,000 draws of the second representation of the sample R2
  # given in the test of the threshold above (.9039 in another 400,000).
  # Tolerance: four standard errors at 400 replications,
  # 4 * sqrt(.905 * .095 / 400) = .059.
  raised <- seq_aipe_R2_sensitivity(0.50, width = 0.30, p = 10, G = 400,
                                    seed = 113, m_0 = 200)
  expect_identical(sensitivity_value(raised, "mean_n"), 200)
  expect_lt(abs(sensitivity_value(raised, "coverage") - 0.9046), 0.059)
})

test_that("a one dimensional array of predictor values is a single predictor", {
  # tapply() returns one; .drop_one_dim() (R/checks.R) treats it as the
  # vector it is, as every other data argument of the package does.
  withr::local_seed(113)
  x <- stats::rnorm(30)
  y <- x + stats::rnorm(30)
  expect_identical(seq_aipe_R2(array(x, 30), y, width = 0.3),
                   seq_aipe_R2(x, y, width = 0.3))
  expect_error(seq_aipe_R2(table(rep(1:3, 10)), y[1:3], width = 0.3),
               "'x' must be a numeric matrix or data frame of predictors")
})
