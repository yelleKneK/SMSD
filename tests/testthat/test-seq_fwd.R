# Tests for the purely sequential fixed width confidence interval procedures:
# seq_fwd_norm_mean(), seq_fwd_dfr_mean(), seq_fwd_nexp_mean(), and their
# pilots. Expected values are hand computed, or derived in the test from the
# coverage probability that defines each criterion; none is the package's own
# formula typed again.

# ---- pilots -----------------------------------------------------------------

test_that("fixed width pilots equal max(m_0, ceiling(K^(1 / (1 + gamma))))", {
  # z_.025 = 1.959964. With d = 0.1, z / d = 19.60, so the pilot is 20: the
  # pilot of Kelley, Darku, and Chattopadhyay (2018, Equation 21),
  # ceiling(2 z / omega), at omega = 2 d = 0.2. (Not 21: an earlier version
  # added 1 after the ceiling.)
  expect_identical(seq_fwd_norm_mean_pilot(d = 0.1), 20L)
  expect_identical(seq_fwd_dfr_mean_pilot(d = 0.1), 20L)
  # d = 0.5: z / d = 3.92 -> 4, below the default minimum m_0 = 10.
  expect_identical(seq_fwd_norm_mean_pilot(d = 0.5), 10L)
  expect_identical(seq_fwd_dfr_mean_pilot(d = 0.5), 10L)
  expect_identical(seq_fwd_norm_mean_pilot(d = 0.5, m_0 = 2), 4L)
  # gamma = 0.49: (19.59964^2)^(1 / 1.49) = 54.27 -> 55.
  expect_identical(seq_fwd_norm_mean_pilot(d = 0.1, gamma = 0.49), 55L)
  # 90%: z_.05 = 1.644854; d = 0.02 gives 82.24 -> 83.
  expect_identical(seq_fwd_dfr_mean_pilot(d = 0.02, conf_level = 0.90), 83L)
  # Negative exponential: K = log(1 / alpha) / d. At 95% and d = 0.02,
  # K = 2.995732 / 0.02 = 149.79 and sqrt(K) = 12.24 -> 13. At 99% and
  # d = 0.01, K = 4.605170 / 0.01 = 460.52 and sqrt(K) = 21.46 -> 22.
  expect_identical(seq_fwd_nexp_mean_pilot(d = 0.02), 13L)
  expect_identical(seq_fwd_nexp_mean_pilot(d = 0.01, conf_level = 0.99), 22L)
  expect_identical(seq_fwd_nexp_mean_pilot(d = 2), 10L)   # sqrt(1.50) -> 2
})

test_that("fixed width pilots are of lower order than the required sample size", {
  # The required sample size is n* = (z / d)^2 sigma^2. A pilot equal to n*
  # at sigma = 1 (the defect this replaces) is 1537 at d = 0.05, sixteen
  # times the 96 needed when sigma = 0.25. The pilot must stay below that 96.
  z <- stats::qnorm(0.975)
  expect_lt(seq_fwd_norm_mean_pilot(d = 0.05), (z / 0.05)^2 * 0.25^2)
  expect_lt(seq_fwd_dfr_mean_pilot(d = 0.05), (z / 0.05)^2 * 0.25^2)
  # Same for the negative exponential: n* = log(1 / alpha) sigma / d is 37
  # at d = 0.02 and sigma = 0.25.
  expect_lt(seq_fwd_nexp_mean_pilot(d = 0.02), -log(0.05) / 0.02 * 0.25)
  # And the pilot becomes a vanishing fraction of n* (sigma = 1) as d shrinks.
  d <- c(0.2, 0.05, 0.0125)
  frac <- vapply(d, function(di) seq_fwd_norm_mean_pilot(di) / (z / di)^2,
                 numeric(1))
  expect_true(all(diff(frac) < 0))
  expect_lt(frac[3], 0.01)
})

test_that("fixed width pilots check every argument", {
  for (pilot in list(seq_fwd_norm_mean_pilot, seq_fwd_dfr_mean_pilot,
                     seq_fwd_nexp_mean_pilot)) {
    expect_error(pilot(), "'d' must be supplied")
    expect_error(pilot(d = NA_real_), "'d' must be a single finite number")
    expect_error(pilot(d = -1), "'d' must be a single finite number")
    expect_error(pilot(d = Inf), "'d' must be a single finite number")
    expect_error(pilot(d = c(0.1, 0.2)), "'d' must be a single finite number")
    expect_error(pilot(d = 0.1, conf_level = 1), "'conf_level'")
    expect_error(pilot(d = 0.1, conf_level = NA_real_), "'conf_level'")
    expect_error(pilot(d = 0.1, conf_level = 95), "'conf_level'")
    expect_error(pilot(d = 0.1, m_0 = 1), "'m_0' must be a single whole number")
    expect_error(pilot(d = 0.1, m_0 = 4.5), "'m_0'")
    expect_error(pilot(d = 0.1, m_0 = NA), "'m_0'")
    expect_error(pilot(d = 0.1, gamma = 0), "'gamma'")
    expect_error(pilot(d = 0.1, gamma = NA_real_), "'gamma'")
    expect_error(pilot(d = 0.1, gamma = "a"), "'gamma'")
    # A target so narrow that the pilot cannot be represented is an error,
    # not an NA.
    expect_error(pilot(d = 1e-20),
                 "implied pilot sample size .* is too large to represent")
    expect_type(pilot(d = 0.1), "integer")
  }
  # The message says what d is in each procedure.
  expect_error(seq_fwd_norm_mean_pilot(d = 0), "half width")
  expect_error(seq_fwd_dfr_mean_pilot(d = 0), "half width")
  expect_error(seq_fwd_nexp_mean_pilot(d = 0), "full width")
  # So does the advice that closes the too large message: it names the
  # argument to check, not "the precision target and the other arguments".
  expect_error(seq_fwd_norm_mean_pilot(d = 1e-20),
               "represent\\. Check 'd' \\(the half width of the confidence interval\\)\\.$")
  expect_error(seq_fwd_dfr_mean_pilot(d = 1e-20),
               "represent\\. Check 'd' \\(the half width of the confidence interval\\)\\.$")
  expect_error(seq_fwd_nexp_mean_pilot(d = 1e-20),
               "represent\\. Check 'd' \\(the full width of the confidence interval\\)\\.$")
})

# ---- criteria against the coverage probability ------------------------------

test_that("seq_fwd_norm_mean's criterion is the smallest n with the stated coverage", {
  # With the variance known to equal v, xbar +/- d covers the mean with
  # probability 2 * pnorm(d * sqrt(n / v)) - 1. The criterion is the smallest
  # n for which that is at least conf_level when v is the sample variance.
  set.seed(113)
  for (cfg in list(c(d = 0.5, conf = 0.95, v = 4), c(d = 0.2, conf = 0.90, v = 0.36),
                   c(d = 3, conf = 0.99, v = 150))) {
    x <- seqc_scale_var(rnorm(25, 50, 1), cfg[["v"]])
    out <- seq_fwd_norm_mean(x, d = cfg[["d"]], conf_level = cfg[["conf"]])
    n_opt <- out$value[out$term == "optimal_n"]
    cover <- function(n) 2 * stats::pnorm(cfg[["d"]] * sqrt(n / cfg[["v"]])) - 1
    expect_gte(cover(n_opt), cfg[["conf"]])
    expect_lt(cover(n_opt - 1), cfg[["conf"]])
    expect_identical(attr(out, "criterion_n"), as.integer(n_opt))
    expect_identical(attr(out, "is_satisfied"), 25 >= n_opt)
  }
})

test_that("seq_fwd_dfr_mean adds 1 / n to the variance in its criterion", {
  # Same coverage identity, with the variance v + 1 / n at the current n.
  set.seed(113)
  for (n in c(12L, 40L)) {
    v <- 2.25
    x <- seqc_scale_var(rnorm(n, 5, 1), v)
    out <- seq_fwd_dfr_mean(x, d = 0.4, conf_level = 0.95)
    n_opt <- out$value[out$term == "optimal_n"]
    cover <- function(k) 2 * stats::pnorm(0.4 * sqrt(k / (v + 1 / n))) - 1
    expect_gte(cover(n_opt), 0.95)
    expect_lt(cover(n_opt - 1), 0.95)
    # It can never ask for fewer observations than the normal theory rule.
    expect_gte(n_opt, seq_fwd_norm_mean(x, d = 0.4)$value[4])
  }
  # 12 observations with variance 2.25 and d = 0.4: z^2 (2.25 + 1 / 12) / 0.16
  # = 56.02, so 57; the normal theory rule gives 54.02, so 55.
  x <- seqc_scale_var(rnorm(12, 5, 1), 2.25)
  expect_equal(seq_fwd_dfr_mean(x, d = 0.4)$value[4], 57)
  expect_equal(seq_fwd_norm_mean(x, d = 0.4)$value[4], 55)
})

test_that("seq_fwd_nexp_mean's criterion is the smallest n with the stated coverage", {
  # X_(1) - theta is exponential with scale sigma / n, so [X_(1) - d, X_(1)]
  # covers theta with probability pexp(d, rate = n / sigma). The criterion is
  # the smallest n for which that is at least conf_level when sigma is
  # U = sum(x - min(x)) / (n - 1).
  set.seed(113)
  x <- 3 + rexp(30, rate = 1 / 2)
  U <- sum(x - min(x)) / 29
  for (cfg in list(c(d = 0.4, conf = 0.95), c(d = 0.1, conf = 0.99))) {
    out <- seq_fwd_nexp_mean(x, d = cfg[["d"]], conf_level = cfg[["conf"]])
    n_opt <- out$value[out$term == "optimal_n"]
    expect_gte(stats::pexp(cfg[["d"]], rate = n_opt / U), cfg[["conf"]])
    expect_lt(stats::pexp(cfg[["d"]], rate = (n_opt - 1) / U), cfg[["conf"]])
    expect_identical(attr(out, "is_satisfied"), 30 >= n_opt)
  }
})

test_that("fixed width criteria round up, so a fractional criterion is not met early", {
  # A sample variance chosen so that z^2 v / d^2 = 10.7: eleven observations
  # are needed, and ten must not satisfy the rule.
  set.seed(113)
  d <- 0.5
  z2 <- stats::qnorm(0.975)^2
  x <- seqc_scale_var(rnorm(10), 10.7 * d^2 / z2)
  out <- seq_fwd_norm_mean(x, d = d)
  expect_equal(out$value[out$term == "optimal_n"], 11)
  expect_false(attr(out, "is_satisfied"))
  x11 <- seqc_scale_var(rnorm(11), 10.7 * d^2 / z2)
  expect_true(attr(seq_fwd_norm_mean(x11, d = d), "is_satisfied"))
})

# ---- what is returned -------------------------------------------------------

test_that("fixed width intervals for a mean are the sample mean plus and minus d", {
  set.seed(113)
  x <- rnorm(200, 10, 2)
  for (f in list(seq_fwd_norm_mean, seq_fwd_dfr_mean)) {
    out <- f(x, d = 0.5, conf_level = 0.90)
    expect_s3_class(out, c("smsd_check", "data.frame"), exact = TRUE)
    expect_identical(out$term, c("lower_limit", "mean", "upper_limit",
                                 "optimal_n", "current_n"))
    expect_equal(out$value[1:3], mean(x) + c(-0.5, 0, 0.5))
    expect_equal(out$value[5], 200)
    expect_equal(attr(out, "width"), 1)             # full width 2 d
    expect_equal(attr(out, "conf_level"), 0.90)
    expect_identical(attr(out, "current_n"), 200L)
    expect_identical(attr(out, "variance_unusable"), FALSE)
    expect_identical(attr(out, "notes"), character(0))
  }
  expect_identical(attr(seq_fwd_norm_mean(x, d = 0.5), "method"),
                   "fixed width interval (normal mean)")
  expect_identical(attr(seq_fwd_dfr_mean(x, d = 0.5), "method"),
                   "fixed width interval (distribution free mean)")
})

test_that("seq_fwd_nexp_mean reports [min - d, min] with d the full width", {
  set.seed(113)
  x <- 3 + rexp(80, rate = 2)
  out <- seq_fwd_nexp_mean(x, d = 0.05, conf_level = 0.99)
  expect_identical(out$term, c("lower_limit", "location", "upper_limit",
                               "optimal_n", "current_n"))
  expect_equal(out$value[1:3], min(x) + c(-0.05, 0, 0))
  expect_equal(attr(out, "width"), 0.05)            # d itself, not 2 d
  expect_equal(attr(out, "conf_level"), 0.99)
  expect_identical(attr(out, "method"),
                   "fixed width interval (negative exponential location)")
  # A negative location is legitimate: shifting the data shifts the interval
  # and changes nothing else.
  shifted <- seq_fwd_nexp_mean(x - 10, d = 0.05, conf_level = 0.99)
  expect_equal(shifted$value[1:3], out$value[1:3] - 10)
  expect_equal(shifted$value[4:5], out$value[4:5])
})

test_that("the data argument is x, as for every one sample procedure", {
  set.seed(113)
  obs <- rnorm(30, 10, 2)
  for (f in list(seq_fwd_norm_mean, seq_fwd_dfr_mean, seq_fwd_nexp_mean)) {
    expect_identical(names(formals(f)), c("x", "d", "conf_level", "na.rm"))
    # By name and by position give the same object.
    expect_identical(f(x = obs, d = 0.5), f(obs, 0.5))
    # The name used before 1.0.0 is not accepted. (The error is R's own
    # "unused argument"; its wording is not matched because R translates it.)
    expect_error(f(data = obs, d = 0.5))
  }
  for (pilot in list(seq_fwd_norm_mean_pilot, seq_fwd_dfr_mean_pilot,
                     seq_fwd_nexp_mean_pilot)) {
    expect_identical(names(formals(pilot)), c("d", "conf_level", "m_0", "gamma"))
  }
})

test_that("missing values are removed or refused as na.rm says", {
  set.seed(113)
  x <- rnorm(30, 10, 2)
  for (f in list(seq_fwd_norm_mean, seq_fwd_dfr_mean, seq_fwd_nexp_mean)) {
    expect_equal(f(c(x, NA), d = 0.5)$value, f(x, d = 0.5)$value)
    expect_error(f(c(x, NA), d = 0.5, na.rm = FALSE), "missing values")
    expect_error(f(x, d = 0.5, na.rm = NA), "'na.rm' must be TRUE or FALSE")
  }
})

# ---- degenerate and invalid input -------------------------------------------

test_that("constant data are 'not met' with a note, never an optimal n of 0", {
  for (f in list(seq_fwd_norm_mean, seq_fwd_dfr_mean, seq_fwd_nexp_mean)) {
    for (x in list(rep(5, 2), rep(5, 40), rep(0.1, 7))) {
      out <- f(x, d = 0.5)
      expect_identical(attr(out, "is_satisfied"), FALSE)
      expect_identical(attr(out, "variance_unusable"), TRUE)
      expect_true(is.na(out$value[out$term == "optimal_n"]))
      expect_true(is.na(attr(out, "criterion_n")))
      # Exactly equal observations: the package's standard note, alone.
      expect_length(attr(out, "notes"), 1L)
      expect_match(attr(out, "notes"),
                   "not yet usable .*zero, negative, or not finite")
      expect_output(print(out), "Stopping rule met: FALSE")
    }
  }
  # One differing observation and the variance is usable again.
  expect_identical(attr(seq_fwd_norm_mean(c(rep(5, 40), 5.001), d = 0.5),
                        "variance_unusable"), FALSE)
})

test_that("data that are constant up to rounding error are 'not met' as well", {
  # Ten gain scores of 0.3, each a posttest minus a pretest. On paper they
  # are equal; in floating point they differ by about 4e-15 and have a
  # "variance" near 2e-30, which satisfied every rule with a criterion of 1.
  pre <- c(101, 57, 12, 80, 33, 64, 90, 45, 71, 28)
  gain <- (pre + 0.3) - pre
  expect_gt(stats::var(gain), 0)
  expect_lt(max(gain) - min(gain), 1e-14)
  # The same thing with no large intermediate at all.
  small <- c(1.3 - 1, 0.5 - 0.2, 0.9 - 0.6, 0.3)
  expect_gt(stats::var(small), 0)
  for (f in list(seq_fwd_norm_mean, seq_fwd_dfr_mean, seq_fwd_nexp_mean)) {
    for (x in list(gain, small)) {
      out <- f(x, d = 0.05)
      expect_identical(attr(out, "is_satisfied"), FALSE)
      expect_identical(attr(out, "variance_unusable"), TRUE)
      expect_true(is.na(out$value[out$term == "optimal_n"]))
      # Two notes: the standard one, and one that says why observations that
      # are not all equal were treated as constant and what to do about it.
      notes <- attr(out, "notes")
      expect_length(notes, 2L)
      expect_match(notes[1], "not yet usable .*zero, negative, or not finite")
      expect_match(notes[2], "agree to ten significant digits")
      expect_match(notes[2], "no more than 1e-10 times their largest absolute value")
      expect_match(notes[2], "subtract a common value from every observation")
      expect_output(print(out), "Note: The observations are not all equal")
    }
  }
})

test_that("a small real spread is not mistaken for rounding error", {
  # Four values near 5 that differ in the sixth decimal place. In units of
  # 1e-6 the deviations from their mean (-0.025) are -0.975, 0.325, 1.025,
  # and -0.375, with squares summing to 2.2475, so the variance is
  # 0.749167e-12. With d = 1e-7 the criterion is
  # 3.841459 * 0.749167e-12 / 1e-14 = 287.8, so 288.
  x <- 5 + c(-1, 0.3, 1, -0.4) * 1e-6
  out <- seq_fwd_norm_mean(x, d = 1e-7)
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_equal(out$value[out$term == "optimal_n"], 288)
  # Times of day recorded in seconds since 1970, a few seconds apart: a
  # range of 12 on values near 1.7e9 is a relative spread of 7e-9.
  stamps <- 1.7e9 + c(0, 3, 7, 12, 4)
  expect_identical(attr(seq_fwd_dfr_mean(stamps, d = 1), "variance_unusable"),
                   FALSE)
  # The documented threshold: a range of at most 1e-10 times the largest
  # absolute value counts as no spread; anything above it is used.
  expect_identical(attr(seq_fwd_norm_mean(c(1, 1 + 5e-11, 1 + 2e-11), d = 1),
                        "variance_unusable"), TRUE)
  expect_identical(attr(seq_fwd_norm_mean(c(1, 1 + 2e-10, 1 + 1e-10), d = 1),
                        "variance_unusable"), FALSE)
  # The threshold is relative: the same pattern at a size of 1e6.
  expect_identical(attr(seq_fwd_norm_mean(1e6 * c(1, 1 + 5e-11, 1 + 2e-11), d = 1),
                        "variance_unusable"), TRUE)
  expect_identical(attr(seq_fwd_norm_mean(1e6 * c(1, 1 + 2e-10, 1 + 1e-10), d = 1),
                        "variance_unusable"), FALSE)
  # One tolerance for the whole family (and the package): the same two
  # samples on either side of it in the other two procedures.
  for (f in list(seq_fwd_dfr_mean, seq_fwd_nexp_mean)) {
    expect_identical(attr(f(c(1, 1 + 5e-11, 1 + 2e-11), d = 1),
                          "variance_unusable"), TRUE)
    expect_identical(attr(f(c(1, 1 + 2e-10, 1 + 1e-10), d = 1),
                          "variance_unusable"), FALSE)
  }
  # Real variation that small beside the size of the data is reported the
  # same way, and the second note says what to do. After a common value is
  # subtracted the same data are an ordinary case: 1e15 + (0, 1, 2, 1, 0)
  # has variance 0.7, and 3.841459 * 0.7 / 0.25 = 10.76, so 11.
  big <- 1e15 + c(0, 1, 2, 1, 0)
  out <- seq_fwd_norm_mean(big, d = 0.5)
  expect_identical(attr(out, "variance_unusable"), TRUE)
  expect_match(attr(out, "notes")[2], "subtract a common value")
  centered <- seq_fwd_norm_mean(big - 1e15, d = 0.5)
  expect_identical(attr(centered, "variance_unusable"), FALSE)
  expect_equal(centered$value[centered$term == "optimal_n"], 11)
})

test_that("a spread that cannot be represented is an error, not a 'not yet' result", {
  x <- c(1.3, 0.2, 2.9, 0.7, 1.1)
  # Squares of data near 1e160 overflow, so the sample variance is Inf. More
  # data cannot change that; the message says to change the units.
  expect_identical(stats::var(x * 1e160), Inf)
  expect_error(seq_fwd_norm_mean(x * 1e160, d = 1e158),
               "sample variance of 'x' is too large to represent")
  expect_error(seq_fwd_dfr_mean(x * 1e160, d = 1e158),
               "sample variance of 'x' is too large to represent")
  # Squares of data near 1e-200 underflow to zero although the data differ.
  expect_identical(stats::var(x * 1e-200), 0)
  expect_error(seq_fwd_norm_mean(x * 1e-200, d = 1),
               "sample variance of 'x' is too small to represent")
  expect_error(seq_fwd_dfr_mean(x * 1e-200, d = 1),
               "sample variance of 'x' is too small to represent")
  # The negative exponential scale estimate is a sum of differences; it
  # overflows only at the very top of the range.
  expect_error(seq_fwd_nexp_mean(x * 6e307, d = 1e307),
               "scale estimate of 'x' is too large to represent")
})

test_that("a limit or a width that cannot be represented is an error, never Inf", {
  # The largest double is about 1.797e308. Two observations of -1.7e308 and
  # d = 5e307 put the lower limit at -2.2e308, which is not a double; it was
  # returned as -Inf. (The full width 2 d = 1e308 is representable, so this
  # is the limit alone.) The data are constant, so no variance error comes
  # first.
  low <- rep(-1.7e308, 2)
  high <- rep(1.7e308, 2)
  expect_error(seq_fwd_norm_mean(low, d = 5e307),
               "'lower_limit' row of the result is too large to represent \\(-Inf\\)\\. Express the data in larger units")
  expect_error(seq_fwd_dfr_mean(low, d = 5e307),
               "'lower_limit' row of the result is too large to represent")
  expect_error(seq_fwd_nexp_mean(low, d = 5e307),
               "'lower_limit' row of the result is too large to represent")
  expect_error(seq_fwd_norm_mean(high, d = 5e307),
               "'upper_limit' row of the result is too large to represent \\(Inf\\)")
  expect_error(seq_fwd_dfr_mean(high, d = 5e307),
               "'upper_limit' row of the result is too large to represent")
  # The negative exponential interval ends at the sample minimum, so it has
  # no limit above the data: 1.7e308 - 5e307 = 1.2e308 is an ordinary row.
  out <- seq_fwd_nexp_mean(high, d = 5e307)
  expect_identical(out$value[1:3], c(1.7e308 - 5e307, 1.7e308, 1.7e308))
  # With d = 5e306 every limit is representable (1.7e308 + 5e306 = 1.75e308)
  # and the result is the ordinary one for constant data: not met, a note,
  # finite rows.
  for (f in list(seq_fwd_norm_mean, seq_fwd_dfr_mean)) {
    out <- f(high, d = 5e306)
    expect_identical(out$value[1:3], c(1.7e308 - 5e306, 1.7e308, 1.7e308 + 5e306))
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_identical(attr(out, "variance_unusable"), TRUE)
  }
  # A half width beyond half of the largest double has a full width 2 d that
  # is not a double; the width attribute was Inf. Here the limits themselves
  # (2.5 - 1e308 and 2.5 + 1e308) are representable.
  x <- c(1, 2, 3, 4)
  expect_error(seq_fwd_norm_mean(x, d = 1e308),
               "full width of the interval, twice 'd', is too large to represent \\(Inf\\)\\. Express the data in larger units \\(divide them by a constant\\) and state 'd' in the same units\\.$")
  expect_error(seq_fwd_dfr_mean(x, d = 1e308),
               "full width of the interval, twice 'd', is too large to represent")
  # 2 * 8.9e307 = 1.78e308 is still a double.
  expect_identical(attr(seq_fwd_norm_mean(x, d = 8.9e307), "width"), 1.78e308)
  # In the negative exponential interval d is the full width, so 1e308 is
  # an (absurd but) representable target.
  out <- seq_fwd_nexp_mean(x, d = 1e308)
  expect_identical(attr(out, "width"), 1e308)
  expect_identical(out$value[1:3], c(1 - 1e308, 1, 1))
})

test_that("fixed width checks refuse invalid arguments and data", {
  set.seed(113)
  x <- rnorm(30, 10, 2)
  for (f in list(seq_fwd_norm_mean, seq_fwd_dfr_mean, seq_fwd_nexp_mean)) {
    expect_error(f(x), "'d' must be supplied")
    expect_error(f(x, d = NA_real_), "'d' must be a single finite number")
    expect_error(f(x, d = Inf), "'d' must be a single finite number")
    expect_error(f(x, d = -0.5), "'d' must be a single finite number")
    expect_error(f(x, d = "a"), "'d' must be a single finite number")
    expect_error(f(x, d = 0.5, conf_level = 0), "'conf_level'")
    expect_error(f(x, d = 0.5, conf_level = NA_real_), "'conf_level'")
    expect_error(f(d = 0.5), "'x' must be supplied")
    # Too few observations: the message names the argument, counts in the
    # singular, and says what to do in a sequential procedure.
    expect_error(f(1, d = 0.5),
                 "'x' needs at least 2 observations; it has 1\\. Collect more data, then check the stopping rule again\\.")
    expect_error(f(letters, d = 0.5), "'x' must be a numeric vector")
    expect_error(f(matrix(x, 15), d = 0.5), "'x' must be a numeric vector, not a matrix")
    expect_error(f(data.frame(x), d = 0.5), "not a data frame")
    expect_error(f(c(x, Inf), d = 0.5), "'x' contains infinite values")
    # A target too narrow for the sample size to be represented is an error
    # rather than an NA in the stopping flag.
    expect_error(f(x, d = 1e-12), "implied sample size .* is too large to represent")
  }
  expect_error(seq_fwd_norm_mean(x, d = 0), "half width")
  expect_error(seq_fwd_dfr_mean(x, d = 0), "half width")
  expect_error(seq_fwd_nexp_mean(x, d = 0), "full width")
  # The advice that closes the too large message names what to check.
  expect_error(seq_fwd_norm_mean(x, d = 1e-12),
               "Check 'd' \\(the half width of the confidence interval\\) and that it is in the units of 'x'\\.$")
  expect_error(seq_fwd_dfr_mean(x, d = 1e-12),
               "Check 'd' \\(the half width of the confidence interval\\) and that it is in the units of 'x'\\.$")
  expect_error(seq_fwd_nexp_mean(x, d = 1e-12),
               "Check 'd' \\(the full width of the confidence interval\\) and that it is in the units of 'x'\\.$")
})

test_that("a criterion is never below one observation", {
  # Data with a spread near 1e-160 and a target of d = 1000: the variance is
  # positive, but z^2 s^2 / d^2 underflows to exactly zero. A target that
  # loose is met at once, with a criterion of 1, not 0.
  x <- c(1.3, 0.2, 2.9, 0.7, 1.1) * 1e-160
  expect_gt(stats::var(x), 0)
  expect_identical(stats::qnorm(0.975)^2 * stats::var(x) / 1000^2, 0)
  out <- seq_fwd_norm_mean(x, d = 1000)
  expect_equal(out$value[out$term == "optimal_n"], 1)
  expect_true(attr(out, "is_satisfied"))
  expect_identical(attr(out, "variance_unusable"), FALSE)
  # A merely loose target gives a criterion of 1 by ordinary rounding up.
  expect_equal(seq_fwd_nexp_mean(c(1, 1.2, 1.5, 1.1), d = 1e6)$value[4], 1)
})

test_that("integer data cannot overflow in the range or the scale estimate", {
  # Balances stored as integers: -2e9, 2e9, 0. Their range, 4e9, is beyond
  # the largest integer, so in integer arithmetic max - min (and x - min(x))
  # is NA with a warning. The data must be handled as doubles.
  x <- c(-2000000000L, 2000000000L, 0L)
  expect_warning(max(x) - min(x), "integer overflow")
  # Negative exponential: U = (0 + 4e9 + 2e9) / 2 = 3e9, and with d = 1e9
  # the criterion is log(20) * 3 = 8.99, so 9.
  out <- expect_no_warning(seq_fwd_nexp_mean(x, d = 1e9))
  expect_equal(out$value[out$term == "optimal_n"], 9)
  expect_identical(attr(out, "is_satisfied"), FALSE)
  expect_identical(attr(out, "variance_unusable"), FALSE)
  # Mean: the variance is (4e18 + 4e18 + 0) / 2 = 4e18, and with d = 1e9 the
  # criterion is 3.841459 * 4 = 15.37, so 16.
  out <- expect_no_warning(seq_fwd_norm_mean(x, d = 1e9))
  expect_equal(out$value[out$term == "optimal_n"], 16)
  expect_identical(attr(out, "variance_unusable"), FALSE)
})

# ---- operating characteristics ----------------------------------------------
# Each test runs the complete procedure 2,000 times from the package's own
# pilot and compares coverage and the mean final sample size with their
# targets. The coverage of each procedure at its test target is a little
# under the nominal 0.95 (the guarantee is asymptotic), so the comparison is
# with the coverage found by an independent simulation of 400,000 runs
# (Monte Carlo SE 0.0004; an engine that never calls the package),
# with a tolerance of four standard errors of the 2,000 run estimate on
# either side: 4 * sqrt(0.95 * 0.05 / 2000) = 0.0195. The seed is fixed, so
# the tests cannot flake.

test_that("seq_fwd_norm_mean attains its confidence level at a moderate target", {
  # N(10, 2^2), d = 0.25, 95%: n* = z^2 sigma^2 / d^2 = 245.9.
  set.seed(113)
  d <- 0.25
  sigma <- 2
  z <- stats::qnorm(0.975)
  n_star <- z^2 * sigma^2 / d^2
  m <- seq_fwd_norm_mean_pilot(d)
  sim <- seqc_simulate(
    2000, m, rgen = function(k) rnorm(k, 10, sigma),
    ok = function(x, n) n >= (z / d)^2 * seqc_run_mv(x - 10, n)$var,
    n_guess = n_star,
    measure = function(x, N) abs(mean(x[seq_len(N)]) - 10) <= d,
    check = function(x) seq_fwd_norm_mean(x, d = d))
  expect_true(sim$agree)
  # Coverage at this target is 0.949 (independent simulation); SE of a
  # coverage near 0.95 from 2,000 runs is 0.0049, so 4 SE = 0.0195.
  expect_lt(abs(sim$mean_y - 0.949), 0.0195)
  # E[N] is about one observation below n* and the SE of the mean of N is
  # about 0.5, so the ratio is within 0.02 of 1 by a wide margin.
  expect_lt(abs(sim$mean_n / n_star - 1), 0.02)
})

test_that("seq_fwd_norm_mean does not presume a scale: sd 0.25, 1, and 3", {
  # d = 0.1, 95%: n* = 384.1 sigma^2, that is 24, 384, and 3457. The pilot is
  # 20 whatever the scale. 1,000 runs at each scale; an independent 100,000
  # run simulation gave E[N] / n* = 1.026, 0.997, 1.000 (the first is above 1
  # because the pilot is most of the 24 needed) and coverage 0.947, 0.949,
  # and 0.951.
  d <- 0.1
  z <- stats::qnorm(0.975)
  m <- seq_fwd_norm_mean_pilot(d)
  expected_ratio <- c(1.026, 0.997, 1.000)
  expected_cover <- c(0.947, 0.949, 0.951)
  # sd(N) is about 5, 28, and 83, so 4 SE of the mean of N over n* is 0.026,
  # 0.009, and 0.003.
  tol_ratio <- c(0.03, 0.015, 0.01)
  set.seed(113)
  for (i in 1:3) {
    sigma <- c(0.25, 1, 3)[i]
    n_star <- (z / d)^2 * sigma^2
    sim <- seqc_simulate(
      1000, m, rgen = function(k) rnorm(k, 0, sigma),
      ok = function(x, n) n >= (z / d)^2 * seqc_run_mv(x, n)$var,
      n_guess = n_star,
      measure = function(x, N) abs(mean(x[seq_len(N)])) <= d,
      check = function(x) seq_fwd_norm_mean(x, d = d),
      n_confirm = 20L, n_full = 1L)
    expect_true(sim$agree)
    expect_lt(abs(sim$mean_n / n_star - expected_ratio[i]), tol_ratio[i])
    # SE of a coverage near 0.95 from 1,000 runs is 0.0069; 4 SE = 0.028.
    expect_lt(abs(sim$mean_y - expected_cover[i]), 0.028)
  }
})

test_that("seq_fwd_dfr_mean attains its confidence level for skewed data at a narrow target", {
  # Exponential data with mean and sd 2, d = 0.125, 95%: n* = 983.4. (At
  # wider targets coverage for skewed data is well below nominal; the help
  # page gives the figures.)
  set.seed(113)
  d <- 0.125
  sigma <- 2
  z <- stats::qnorm(0.975)
  n_star <- z^2 * sigma^2 / d^2
  m <- seq_fwd_dfr_mean_pilot(d)
  sim <- seqc_simulate(
    2000, m, rgen = function(k) rexp(k, rate = 1 / sigma),
    ok = function(x, n) n >= (z / d)^2 * (seqc_run_mv(x, n)$var + 1 / n),
    n_guess = n_star,
    measure = function(x, N) abs(mean(x[seq_len(N)]) - sigma) <= d,
    check = function(x) seq_fwd_dfr_mean(x, d = d),
    n_confirm = 40L, n_full = 1L)
  expect_true(sim$agree)
  # Coverage at this target is 0.946 (independent simulation), still a
  # little short of 0.95 for skewed data; 4 SE = 0.0195 as above.
  expect_lt(abs(sim$mean_y - 0.946), 0.0195)
  # The SE of the mean of N is about 2 (sd of N near 90), 0.002 of n*; the
  # mean falls about 0.6% short of n*. Within 0.03 of 1.
  expect_lt(abs(sim$mean_n / n_star - 1), 0.03)
})

test_that("seq_fwd_nexp_mean attains its confidence level at a moderate target", {
  # Location 3, scale 2, d = 0.025, 95%: n* = log(20) * 2 / 0.025 = 239.7.
  set.seed(113)
  d <- 0.025
  sigma <- 2
  n_star <- log(20) * sigma / d
  m <- seq_fwd_nexp_mean_pilot(d)
  sim <- seqc_simulate(
    2000, m, rgen = function(k) 3 + rexp(k, rate = 1 / sigma),
    ok = function(x, n) n >= (log(20) / d) * seqc_run_u(x, n),
    n_guess = n_star,
    measure = function(x, N) {
      x1 <- min(x[seq_len(N)])
      x1 - d <= 3 && 3 <= x1
    },
    check = function(x) seq_fwd_nexp_mean(x, d = d))
  expect_true(sim$agree)
  # Coverage at this target is 0.948 (independent simulation); 4 SE = 0.0195.
  expect_lt(abs(sim$mean_y - 0.948), 0.0195)
  # sd of N is about sqrt(n*) = 15.5, so the SE of its mean is 0.35.
  expect_lt(abs(sim$mean_n / n_star - 1), 0.02)
})
