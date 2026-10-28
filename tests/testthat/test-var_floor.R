# Tests for the lower bound n^-3 on the variance estimators and for what the
# seq_aipe_* check functions do with an estimate at or below it.
#
# The rule (R/var_floor.R). A raw estimate that is positive is used as the
# source articles prescribe, max{raw, n^-3}: the stopping rule is evaluated
# and the interval is built, whether or not the bound replaced the estimate.
# A raw estimate that is zero, negative, or not finite is not an estimate of
# a variance: the stopping rule is reported as not met, the confidence
# limits are NA, criterion_n is at least current_n + 1, the
# variance_unusable attribute is TRUE, and a note says what happened. The
# same holds when the point estimate cannot be computed. This is a
# deliberate departure from the stopping rule as printed in the articles,
# which stops at once in that situation with a near zero width interval.

row_value <- function(out, term) out$value[out$term == term]

# The contract shared by every "not yet estimable" result, as one logical
# (for use inside loops) and as a set of expectations (for single cases).
# `unusable`: is the variance estimate reported as unusable? (It is in every
# such result of this family: when the point estimate cannot be computed,
# neither can its variance.)
is_not_estimable <- function(out, unusable = TRUE) {
  ok <- identical(attr(out, "is_satisfied"), FALSE) &&
    identical(attr(out, "variance_unusable"), unusable) &&
    is.na(row_value(out, "lower_limit")) &&
    is.na(row_value(out, "upper_limit")) &&
    attr(out, "criterion_n") > attr(out, "current_n") &&
    identical(row_value(out, "criterion_n"),
              as.numeric(attr(out, "criterion_n"))) &&
    length(attr(out, "notes")) >= 1L
  if (unusable) {
    # xi_squared holds the bound, and the standard note is attached.
    ok <- ok &&
      identical(row_value(out, "xi_squared"), attr(out, "current_n")^(-3)) &&
      .note_variance_unusable() %in% attr(out, "notes")
  }
  ok
}
expect_not_estimable <- function(out, unusable = TRUE) {
  expect_false(attr(out, "is_satisfied"))
  expect_identical(attr(out, "variance_unusable"), unusable)
  expect_true(is.na(row_value(out, "lower_limit")))
  expect_true(is.na(row_value(out, "upper_limit")))
  expect_gt(attr(out, "criterion_n"), attr(out, "current_n"))
  expect_true(is_not_estimable(out, unusable))
}

test_that(".var_floor bounds below, replaces non-finite values, and can be switched off", {
  expect_equal(.var_floor(0.5, 10), 0.5)          # finite value passes through
  expect_equal(.var_floor(1e-9, 10), 10^(-3))     # bounded from below
  expect_equal(.var_floor(-0.3, 10), 10^(-3))     # negative raw estimate
  expect_equal(.var_floor(NaN, 10), 10^(-3))      # non-finite -> bound
  expect_equal(.var_floor(Inf, 10), 10^(-3))
  # floor = FALSE hands back the raw value, whatever it is.
  expect_identical(.var_floor(1e-9, 10, floor = FALSE), 1e-9)
  expect_identical(.var_floor(-0.3, 10, floor = FALSE), -0.3)
  expect_true(is.nan(.var_floor(NaN, 10, floor = FALSE)))
})

test_that("only a raw estimate that is zero, negative, or not finite is unusable", {
  for (raw in c(0, -1e-300, -0.3, NaN, NA_real_, Inf, -Inf)) {
    expect_true(.var_unusable(raw))
  }
  # Any positive number is an estimate, however small; whether it is below
  # the bound is a separate question.
  for (raw in c(1e-300, 1e-9, 10^(-3), 0.5)) {
    expect_false(.var_unusable(raw))
  }
  expect_true(.aipe_below_bound(1e-9, 10))
  expect_false(.aipe_below_bound(10^(-3), 10))     # at the bound: unchanged
  expect_false(.aipe_below_bound(0.5, 10))
  expect_false(.aipe_below_bound(0, 10))
  expect_false(.aipe_below_bound(-0.3, 10))
  expect_false(.aipe_below_bound(NaN, 10))
})

test_that("a positive variance estimate below n^-3 is used through the bound", {
  # Ten scores with a coefficient of variation of about 0.015: 98, 99, 100,
  # 101, 102, twice. By hand: mean 100, sum of squared deviations 20 (so
  # s^2 = 20/9), third central moment 0 by symmetry, sum of fourth powers
  # 68, and the unbiased fourth moment estimate (Kelley, Darku, and
  # Chattopadhyay, 2018, Equation 79)
  #   {(n^2 - 2n + 3) 68 - 3 (2n - 3)/n 20^2} / ((n - 1)(n - 2)(n - 3))
  #     = (83 * 68 - 2040) / 504 = 3604 / 504.
  # Their Equation 32 then gives
  #   V^2 = (3604/504) / (4 * 100^2 * 20/9) - (20/9) / (4 * 100^2)
  #         + (20/9)^2 / 100^4 = 0.0000249,
  # a genuine estimate (the population value for normal data with this
  # coefficient of variation is about 0.0001) that is below 10^-3.
  x <- rep(c(98, 99, 100, 101, 102), 2)
  by_hand <- (3604 / 504) / (4 * 100^2 * 20 / 9) - (20 / 9) / (4 * 100^2) +
    (20 / 9)^2 / 100^4
  expect_equal(by_hand, 2.494e-05, tolerance = 1e-3)
  expect_equal(.var_cv(x, floor = FALSE), by_hand, tolerance = 1e-9)
  expect_identical(.var_cv(x), 10^(-3))
  z <- qnorm(0.975)
  cv <- sd(x) / 100

  # The bound stands in for the estimate, the rule is evaluated, and the
  # interval is built: (2 z / 0.5)^2 (0.001 + 1/10) = 6.2, so the criterion
  # is 7 and ten observations meet it. (When every estimate at or below the
  # bound was treated as unusable this check was "not met" with NA limits,
  # and an estimate of this size would have stayed so until n^-3 fell below
  # it, at 35 observations.)
  out <- seq_aipe_cv(x, width = 0.5)
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_identical(attr(out, "notes"), character(0))
  expect_identical(row_value(out, "xi_squared"), 10^(-3))
  expect_equal(row_value(out, "cv"), cv)
  expect_equal(row_value(out, "lower_limit"), cv - z * sqrt(10^(-3) / 10))
  expect_equal(row_value(out, "upper_limit"), cv + z * sqrt(10^(-3) / 10))
  expect_equal(row_value(out, "criterion_n"), 7)
  expect_true(attr(out, "is_satisfied"))

  # A narrower target is not met, and the criterion is the article's, not
  # one forced above the sample size: (2 z / 0.2)^2 (0.001 + 0.1) = 38.8.
  narrow <- seq_aipe_cv(x, width = 0.2)
  expect_identical(attr(narrow, "variance_unusable"), FALSE)
  expect_false(attr(narrow, "is_satisfied"))
  expect_equal(row_value(narrow, "criterion_n"), 39)
  expect_true(is.finite(row_value(narrow, "lower_limit")))
  expect_equal(row_value(narrow, "upper_limit") -
                 row_value(narrow, "lower_limit"),
               2 * z * sqrt(10^(-3) / 10))

  # The result does not depend on the units of the data.
  expect_equal(seq_aipe_cv(x / 1000, width = 0.5)$value, out$value)
})

test_that("a positive rank based estimate below n^-3 is used through the bound too", {
  # 1, ..., 7 against 3, 2, 1, 7, 6, 5, 4: the 12 pairs between the two
  # blocks are concordant and the 9 within them discordant, so Kendall's
  # tau is 3/21. Its variance estimate is 1/686 in exact rational
  # arithmetic (see test-var_estimators.R): small because the pattern is
  # so regular, but an estimate, and below the bound 7^-3 = 1/343. It is
  # not rounded to zero; the bound stands in for it and the rule is
  # evaluated: (2 z / 0.6)^2 (1/343 + 1/7) = 6.2, so seven pairs meet it.
  x <- c(1, 2, 3, 4, 5, 6, 7)
  y <- c(3, 2, 1, 7, 6, 5, 4)
  out <- seq_aipe_r(x, y, width = 0.6, method = "kendall")
  expect_equal(row_value(out, "r"), 3 / 21)
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_identical(row_value(out, "xi_squared"), 7^(-3))
  expect_equal(row_value(out, "criterion_n"), 7)
  expect_true(attr(out, "is_satisfied"))
  expect_identical(attr(out, "notes"), character(0))
  expect_equal(row_value(out, "upper_limit") - row_value(out, "lower_limit"),
               2 * qnorm(0.975) * sqrt(7^(-3) / 7))
})

test_that("a sample with no variation never reports the stopping rule as met", {
  # The old behavior, pinned by an earlier version of this test, was
  # is_satisfied = TRUE with an estimate of Inf, NaN, or NA.
  # First one case of each kind, in full.
  k <- rep(3, 10)
  v <- seq_len(10) / 7
  no_estimate <- list(
    sm = seq_aipe_sm(k, width = 0.4),
    smd = seq_aipe_smd(k, k + 1, width = 0.4),
    reg_coef = seq_aipe_reg_coef(k, v, width = 0.4),
    pearson = seq_aipe_r(k, v, width = 0.4),
    kendall = seq_aipe_r(v, k, width = 0.4, method = "kendall"),
    spearman = seq_aipe_r(k, v, width = 0.4, method = "spearman")
  )
  for (out in no_estimate) {
    expect_not_estimable(out)
    expect_true(is.na(out$value[2]))
    expect_true(.note_not_estimable() %in% attr(out, "notes"))
  }
  # An estimate exists, but no usable variance: the coefficient of
  # variation and the Gini index of a constant sample are 0, a contrast of
  # constant groups is the contrast of the constants, and the slope on a
  # constant response is 0.
  cv <- seq_aipe_cv(k, width = 0.4)
  expect_not_estimable(cv)
  expect_identical(row_value(cv, "cv"), 0)
  gi <- seq_aipe_gini(k, width = 0.4)
  expect_not_estimable(gi)
  expect_identical(row_value(gi, "gini"), 0)
  ct <- seq_aipe_c(list(k, k + 2, k - 1), c(1, -0.5, -0.5), width = 0.4)
  expect_not_estimable(ct)
  expect_equal(row_value(ct, "contrast"), 3 - 0.5 * 5 - 0.5 * 2)
  rc <- seq_aipe_reg_coef(v, k, width = 0.4)
  expect_not_estimable(rc)
  expect_equal(row_value(rc, "reg_coef"), 0)

  # Then every procedure over a grid of sample sizes and widths, including
  # widths so wide that the published criterion is below the sample size.
  # (The widest target a correlation accepts is below 2 and the widest for
  # the Gini index below 1, the lengths of their ranges; with the bound
  # n^-3 in place of the estimate the published criterion at those widths
  # is 2 for the correlation at every n here, and 2 or less for the Gini
  # index from n = 10 on.)
  ok <- TRUE
  for (n in c(4, 10, 40, 200)) {
    k <- rep(0.1, n)
    v <- seq_len(n) / 7
    for (w in c(0.05, 0.4, 5)) {
      w_r <- min(w, 1.9)
      w_gini <- min(w, 0.99)
      ok <- ok &&
        is_not_estimable(seq_aipe_sm(k, width = w)) &&
        is_not_estimable(seq_aipe_smd(k, k + 1, width = w)) &&
        is_not_estimable(seq_aipe_reg_coef(k, v, width = w)) &&
        is_not_estimable(seq_aipe_reg_coef(v, k, width = w)) &&
        is_not_estimable(seq_aipe_r(k, v, width = w_r)) &&
        is_not_estimable(seq_aipe_r(v, k, width = w_r, method = "kendall")) &&
        is_not_estimable(seq_aipe_r(k, v, width = w_r,
                                    method = "spearman")) &&
        is_not_estimable(seq_aipe_cv(k, width = w)) &&
        is_not_estimable(seq_aipe_gini(k, width = w_gini)) &&
        is_not_estimable(seq_aipe_c(list(k, k + 2), c(1, -1), width = w))
    }
  }
  expect_true(ok)
  expect_true(aipe_rule_met(4, 4^(-3), 1.9))
  expect_true(aipe_rule_met(10, 10^(-3), 0.99))
})

test_that("values that are equal on paper but not in floating point count as no variation", {
  # Twelve gains of 0.3, each computed as a posttest minus a pretest
  # between 3 and 411. They take five distinct values that differ by up to
  # 1.4e-14, and their "variance" is about 3e-29. Computed from that noise,
  # the raw variance estimates of the coefficient of variation, the Gini
  # index, and the contrast come out as positive numbers of that order; a
  # positive estimate is used through the bound n^-3, so with an exact test
  # for equal values these checks would report the rule as met (for the
  # coefficient of variation at width 0.4: criterion 9 with 12
  # observations). The package wide .is_constant() (agreement to ten
  # significant digits) recognizes the sample for what it is.
  pre <- c(3.1, 7.4, 12.6, 25.3, 31.7, 55.2, 63.9, 82.1, 96.6, 120.4, 200.7,
           410.9)
  gain <- (pre + 0.3) - pre
  expect_gt(length(unique(gain)), 1L)             # not exactly equal
  expect_lt(max(gain) - min(gain), 1e-13)
  expect_true(.is_constant(gain))
  v <- seq_along(gain) / 7
  # The published rule with the bound in place of the estimate would be met.
  expect_true(aipe_rule_met(12, 12^(-3), 0.4))

  # No variance estimate is manufactured from the noise.
  expect_true(is.nan(.var_cv(gain, floor = FALSE)))
  expect_true(is.nan(.var_sm(gain, floor = FALSE)))
  expect_identical(.var_gini(gain, floor = FALSE), 0)
  expect_true(is.nan(.var_c(list(gain, gain + 2), c(1, -1), floor = FALSE)))
  expect_identical(.var_reg_coef(v, gain, floor = FALSE), 0)
  expect_true(is.nan(.var_reg_coef(gain, v, floor = FALSE)))
  expect_true(is.nan(.var_smd(gain, gain + 1, floor = FALSE)))
  for (m in c("pearson", "kendall", "spearman")) {
    expect_true(is.nan(.var_r(gain, v, method = m, floor = FALSE)))
  }

  # And no check reports the rule as met.
  expect_not_estimable(seq_aipe_cv(gain, width = 0.4))
  # The coefficient of variation is reported as it was computed, s / xbar
  # (here the ratio of the noise to 0.3, about 1.7e-14), which is what
  # seq_mr_cv() reports for the same sample: the two procedures for the
  # coefficient of variation give the same estimate for every sample. It
  # is exactly zero for a sample of equal values (tested above).
  cv_gain <- row_value(seq_aipe_cv(gain, width = 0.4), "cv")
  by_hand <- sqrt(sum((gain - mean(gain))^2) / 11) / mean(gain)
  expect_gt(by_hand, 0)
  expect_lt(by_hand, 1e-13)
  expect_equal(cv_gain / by_hand, 1, tolerance = 1e-12)
  mr_gain <- seq_mr_cv(gain, A = 1000, c = 1)
  expect_identical(row_value(mr_gain, "cv"), cv_gain)
  # A sample that varies, though too little beside its size to count: both
  # report its coefficient of variation, s / xbar, not zero.
  withr::local_seed(113)
  far <- 1e11 + stats::rnorm(12)
  expect_true(.is_constant(far))
  cv_far <- row_value(seq_aipe_cv(far, width = 0.3), "cv")
  expect_equal(cv_far / (stats::sd(far) / mean(far)), 1, tolerance = 1e-12)
  expect_gt(cv_far, 1e-12)
  expect_identical(row_value(seq_mr_cv(far, A = 1000, c = 1), "cv"), cv_far)
  expect_not_estimable(seq_aipe_cv(far, width = 0.3))
  expect_not_estimable(seq_aipe_gini(gain, width = 0.4))
  expect_not_estimable(seq_aipe_c(list(gain, gain + 2), c(1, -1), width = 0.4))
  expect_not_estimable(seq_aipe_reg_coef(v, gain, width = 0.4))
  # For the Gini index and the slope, exactly zero is reported for a sample
  # with no variation, not the noise (the plain estimators return 9.0e-15
  # for the Gini index and 4.5e-15 for the slope here).
  expect_gt(gini_index(gain), 0)
  expect_false(slr_slope(v, gain) == 0)
  expect_identical(row_value(seq_aipe_gini(gain, width = 0.4), "gini"), 0)
  expect_identical(row_value(seq_aipe_reg_coef(v, gain, width = 0.4),
                             "reg_coef"), 0)
  no_estimate <- list(
    seq_aipe_sm(gain, width = 0.4),
    seq_aipe_smd(gain, gain + 1, width = 0.4),
    seq_aipe_reg_coef(gain, v, width = 0.4),
    seq_aipe_r(gain, v, width = 0.4),
    seq_aipe_r(v, gain, width = 0.4, method = "kendall"),
    seq_aipe_r(gain, v, width = 0.4, method = "spearman")
  )
  for (out in no_estimate) {
    expect_not_estimable(out)
    expect_true(is.na(out$value[2]))
  }
  # One group with no variation and one that varies: the standardized mean
  # difference has an estimate (the pooled standard deviation is positive),
  # but the variance of the constant group has not been estimated, so the
  # variance estimate is not usable and the note names the group.
  withr::local_seed(113)
  other <- rnorm(12)
  one_flat <- seq_aipe_smd(gain, other, width = 0.4)
  expect_not_estimable(one_flat)
  expect_equal(row_value(one_flat, "smd"),
               (mean(gain) - mean(other)) / sqrt(stats::var(other) / 2))
  expect_match(attr(one_flat, "notes")[2L],
               "^The observations in 'group_1' are all equal")
})

test_that("degenerate data produce no warnings", {
  k <- rep(3, 12)
  v <- seq_len(12) / 7
  expect_silent(seq_aipe_r(k, v, width = 0.4))
  expect_silent(seq_aipe_r(k, v, width = 0.4, method = "kendall"))
  expect_silent(seq_aipe_r(v, k, width = 0.4, method = "spearman"))
  expect_silent(seq_aipe_sm(k, width = 0.4))
  expect_silent(seq_aipe_smd(k, k, width = 0.4))
  expect_silent(seq_aipe_reg_coef(k, v, width = 0.4))
})

test_that("a sample mean that is not positive is 'not estimable' for the CV, not an error", {
  # Negative mean: sd / mean would be a negative number, which is not a
  # coefficient of variation.
  x <- c(-4, -6, -5, -7, -3, -5, -6, -4, -5, -5)
  out <- seq_aipe_cv(x, width = 0.3)
  expect_not_estimable(out)
  expect_true(is.na(row_value(out, "cv")))
  expect_match(attr(out, "notes")[1], "sample mean is not positive")
  # A mean of exactly zero (previously an error).
  out0 <- seq_aipe_cv(c(-2, -1, 1, 2, -3, 3), width = 0.3)
  expect_not_estimable(out0)
  expect_match(attr(out0, "notes")[1], "sample mean is not positive")
  # The standardized mean, in contrast, is defined for any mean.
  sm <- seq_aipe_sm(x, width = 0.3)
  expect_equal(row_value(sm, "sm"), mean(x) / sd(x))
  expect_false(attr(sm, "variance_unusable"))
  expect_true(is.finite(row_value(sm, "lower_limit")))
})

test_that("a sample mean that is positive only by rounding is not a positive mean for the CV", {
  # Six scores that sum to 14 units in the last place of 3: the mean is
  # 1.04e-15, positive, and below what the arithmetic can tell from zero
  # for values of this size (ten units in the last place of the largest
  # value, the tolerance seq_mr_cv() applies to the same question). Taken
  # as a real mean it gives a coefficient of variation of 2.3e15, with
  # limits of -/+ 4e30, a criterion of 2.6e63, and no note.
  x <- c(-1, 1, -2, 2, -3, 3 + 14 * 2^-51)
  expect_gt(mean(x), 0)
  expect_lt(mean(x), 2e-15)
  out <- seq_aipe_cv(x, width = 0.4)
  expect_not_estimable(out)
  expect_true(is.na(row_value(out, "cv")))
  expect_match(attr(out, "notes")[1],
               "not distinguishable from zero at machine precision")
  # Centered scores: their mean is rounding noise of either sign, and both
  # signs get the same answer.
  scores <- c(2.1, 3.4, 5.0, 1.7, 4.2, 6.3, 2.8, 3.9, 7.5, 3.1, 4.6, 2.4)
  centered <- scores - mean(scores)
  expect_lt(abs(mean(centered)), 1e-15)
  for (sign in c(1, -1)) {
    out <- seq_aipe_cv(sign * centered, width = 0.4)
    expect_not_estimable(out)
    expect_true(is.na(row_value(out, "cv")))
  }
  # A mean that is small but real is a mean: the same scores with a last
  # value of 3.000006 have a mean of 1e-6 and a coefficient of variation in
  # the millions, which is reported.
  small <- c(-1, 1, -2, 2, -3, 3.000006)
  out <- seq_aipe_cv(small, width = 0.4)
  expect_equal(row_value(out, "cv"), sd(small) / mean(small))
  expect_gt(row_value(out, "cv"), 1e6)
  expect_false(attr(out, "variance_unusable"))
  # .aipe_mean_positive() is the test, in units of the largest value.
  expect_true(.aipe_mean_positive(c(1, 2, 3)))
  expect_false(.aipe_mean_positive(c(-1, 0, 1)))
  expect_false(.aipe_mean_positive(c(-3, -2, -1)))
  expect_false(.aipe_mean_positive(c(-1, 1, 9 * .Machine$double.eps)))
  expect_true(.aipe_mean_positive(c(-1, 1, 40 * .Machine$double.eps)))
  expect_false(.aipe_mean_positive(1e6 * c(-1, 1, 9 * .Machine$double.eps)))
})

test_that("an estimate that is zero on paper is not yet usable, whatever the rounding", {
  # Perfectly discordant pairs: Kendall's tau is -1 and its variance
  # estimate is zero on paper. Computed as rounding noise it was positive
  # at most sample sizes, the bound n^-3 stood in for it, and the check
  # reported the rule as met with limits of -1.020 and -0.980 for a
  # requested width of 0.4 (and as not yet usable at the sample sizes where
  # the noise happened to be exactly zero: 15 and 16, but not 14 or 17).
  x <- c(2.1, 3.4, 1.2, 5.6, 4.4, 3.3, 2.8, 6.1, 0.4, 3.9)
  expect_true(aipe_rule_met(10, 10^(-3), 0.4))   # what the bound would do
  out <- seq_aipe_r(x, exp(-x), width = 0.4, method = "kendall")
  expect_not_estimable(out)
  expect_equal(row_value(out, "r"), -1)
  expect_identical(attr(out, "criterion_n"), 11)
  expect_identical(attr(out, "notes"), .note_variance_unusable())
  ok <- TRUE
  for (n in 4:40) {
    v <- seq_len(n) / 7
    out <- seq_aipe_r(v, -v, width = 0.4, method = "kendall")
    ok <- ok && is_not_estimable(out) &&
      isTRUE(all.equal(row_value(out, "r"), -1))
  }
  expect_true(ok)
  # Two dichotomous variables in perfect opposition: zero on paper for
  # Kendall and for Spearman (met at n = 10 and 20 and unusable at 8, 12,
  # and 30 when the Spearman estimate was computed as noise).
  ok <- TRUE
  for (n in c(8, 10, 12, 20, 30)) {
    b <- rep(c(0, 1), each = n / 2)
    for (m in c("kendall", "spearman")) {
      out <- suppressWarnings(seq_aipe_r(b, 1 - b, width = 0.4, method = m))
      ok <- ok && is_not_estimable(out) &&
        isTRUE(all.equal(row_value(out, "r"), -1))
    }
  }
  expect_true(ok)
  out <- seq_aipe_r(c(1, 2, 3, 4), c(3, 4, 1, 2), width = 1.2,
                    method = "spearman")
  expect_not_estimable(out)
  expect_equal(row_value(out, "r"), -0.6)

  # Collinear pairs whose Pearson and slope variance estimates are zero on
  # paper (see test-var_estimators.R): the same answer for every line
  # through the points. The noise had made the line with slope 1 "met",
  # with limits of 0.969 and 1.031 for the correlation, and the line with
  # slope 0.3 "not yet usable".
  x <- c(3, 3, 3, 2, 0, 0, 0, 1)
  expect_true(aipe_rule_met(8, 8^(-3), 0.6))
  for (b in c(1, 0.3, -1, 2, 1.8)) {
    y <- 10 + b * x
    out <- seq_aipe_r(x, y, width = 0.6)
    expect_not_estimable(out)
    expect_equal(row_value(out, "r"), sign(b))
    out <- seq_aipe_reg_coef(x, y, width = 0.6)
    expect_not_estimable(out)
    expect_equal(row_value(out, "reg_coef"), b)
    # Only the standard note: the bound did not replace a positive estimate.
    expect_identical(attr(out, "notes"), .note_variance_unusable())
  }
})

test_that("a sample of zeros is 'not estimable' for the Gini index, not an error", {
  out <- seq_aipe_gini(rep(0, 10), width = 0.2)
  expect_not_estimable(out)
  expect_true(is.na(row_value(out, "gini")))
  expect_true(.note_not_estimable() %in% attr(out, "notes"))
})

test_that("the unusable path at the pilot: the slope", {
  # At width 0.4 the pilot is 10 pairs, and the raw moment based variance
  # estimate of the slope is zero or negative in about 16% of pilot samples
  # of bivariate normal data (0.1610, se 0.0008, in 200,000 pilots; it is
  # positive but below the bound in another 0.05%). With 300 pilots the
  # Monte Carlo standard error of that proportion is 0.0212, so the window
  # below is four standard errors on each side (0.085).
  withr::local_seed(113)
  w <- 0.4
  z <- qnorm(0.975)
  m <- seq_aipe_reg_coef_pilot(width = w)
  expect_identical(m, 10L)
  # The rule as printed stops at the pilot whenever the estimate is replaced
  # by the bound, with an interval of width 2 z sqrt(n^-3 / n) = 0.039, a
  # tenth of the target.
  expect_true(aipe_rule_met(m, m^(-3), w))
  expect_equal(2 * z * sqrt(m^(-3) / m), 0.0392, tolerance = 1e-3)

  unusable <- logical(300)
  unusable_ok <- usable_ok <- TRUE
  for (i in seq_along(unusable)) {
    x <- rnorm(m)
    y <- 0.5 * x + rnorm(m)
    out <- seq_aipe_reg_coef(x, y, width = w)
    raw <- .var_reg_coef(x, y, floor = FALSE)
    unusable[i] <- attr(out, "variance_unusable")
    if (raw <= 0) {
      # Not met, limits NA, criterion above n, note attached; the slope
      # itself is still reported.
      unusable_ok <- unusable_ok && is_not_estimable(out) &&
        isTRUE(all.equal(row_value(out, "reg_coef"), cov(x, y) / var(x)))
    } else {
      # A positive estimate: the rule as printed, with max{raw, n^-3}.
      xi2 <- row_value(out, "xi_squared")
      usable_ok <- usable_ok && !unusable[i] &&
        isTRUE(all.equal(xi2, max(raw, m^(-3)))) &&
        identical(attr(out, "is_satisfied"), aipe_rule_met(m, xi2, w)) &&
        isTRUE(all.equal(row_value(out, "upper_limit") -
                           row_value(out, "lower_limit"),
                         2 * z * sqrt(xi2 / m))) &&
        length(attr(out, "notes")) == as.integer(raw < m^(-3))
    }
  }
  expect_true(unusable_ok)
  expect_true(usable_ok)
  expect_lt(abs(mean(unusable) - 0.161), 0.085)
})

test_that("the unusable path at the pilot: the Pearson correlation", {
  # Bivariate normal, rho = 0.4, width 0.4 (pilot 10): the raw variance
  # estimate is zero or negative in about 14% of pilots (0.1383, se 0.0008,
  # in 200,000 pilots; positive but below the bound in another 0.08%). The
  # standard error at 300 pilots is 0.0199; four of them: 0.080.
  withr::local_seed(113)
  w <- 0.4
  m <- seq_aipe_r_pilot(width = w)
  unusable <- logical(300)
  unusable_ok <- TRUE
  for (i in seq_along(unusable)) {
    d <- aipe_pairs(m, rho = 0.4)
    out <- seq_aipe_r(d$x, d$y, width = w)
    unusable[i] <- attr(out, "variance_unusable")
    unusable_ok <- unusable_ok &&
      identical(unusable[i], .var_r(d$x, d$y, floor = FALSE) <= 0)
    if (unusable[i]) {
      unusable_ok <- unusable_ok && is_not_estimable(out) &&
        isTRUE(all.equal(row_value(out, "r"), cor(d$x, d$y)))
    }
  }
  expect_true(unusable_ok)
  expect_lt(abs(mean(unusable) - 0.138), 0.080)
})

test_that("is_satisfied is identical to current_n >= criterion_n in every state", {
  withr::local_seed(113)
  invariant <- function(out) {
    identical(attr(out, "is_satisfied"),
              attr(out, "current_n") >= attr(out, "criterion_n"))
  }
  ok <- TRUE
  n_unusable <- 0
  for (i in 1:60) {
    n <- sample(c(4, 6, 10, 25, 80), 1)
    w <- sample(c(0.2, 0.4, 0.8, 2), 1)
    d <- aipe_pairs(n, rho = 0.5)
    # The correlation takes a width below 2 and the Gini index one below 1.
    outs <- list(
      seq_aipe_r(d$x, d$y, width = min(w, 1.9)),
      seq_aipe_reg_coef(d$x, d$y, width = w),
      seq_aipe_smd(d$x, d$y, width = w),
      seq_aipe_sm(d$x, width = w),
      seq_aipe_cv(abs(d$x) + 0.1, width = w),
      seq_aipe_gini(abs(d$x), width = min(w, 0.9)),
      seq_aipe_c(list(d$x, d$y), c(1, -1), width = w)
    )
    ok <- ok && all(vapply(outs, invariant, logical(1)))
    n_unusable <- n_unusable +
      sum(vapply(outs, attr, logical(1), "variance_unusable"))
  }
  expect_true(ok)
  expect_gt(n_unusable, 0)      # the unusable state was among those checked
  # The state in which the bound replaced a positive estimate (the scores
  # of the test above), met and not met.
  x <- rep(c(98, 99, 100, 101, 102), 2)
  expect_true(invariant(seq_aipe_cv(x, width = 0.5)))
  expect_true(invariant(seq_aipe_cv(x, width = 0.2)))
})

test_that("a loop that collects the reported shortfall always moves forward", {
  # criterion_n - current_n is at least 1 whenever the rule is not met, in
  # particular when the variance estimate is unusable, so a driver that
  # adds that many observations cannot stall, and it ends with an interval
  # no wider than the target.
  withr::local_seed(113)
  w <- 0.4
  min_shortfall <- Inf
  max_width <- 0
  n_unusable <- 0
  for (i in 1:40) {
    d <- aipe_pairs(seq_aipe_r_pilot(w), rho = 0.6)
    steps <- 0
    repeat {
      out <- seq_aipe_r(d$x, d$y, width = w)
      if (attr(out, "is_satisfied")) break
      n_unusable <- n_unusable + attr(out, "variance_unusable")
      shortfall <- attr(out, "criterion_n") - attr(out, "current_n")
      min_shortfall <- min(min_shortfall, shortfall)
      more <- aipe_pairs(max(shortfall, 1), rho = 0.6)
      d <- list(x = c(d$x, more$x), y = c(d$y, more$y))
      steps <- steps + 1
      if (steps > 500) stop("the driver did not terminate")
    }
    max_width <- max(max_width, row_value(out, "upper_limit") -
                       row_value(out, "lower_limit"))
  }
  expect_gte(min_shortfall, 1)
  expect_lte(max_width, w)
  expect_gt(n_unusable, 0)      # the driver passed through unusable checks
})

test_that("print() shows the notes of a check that is not yet estimable", {
  out <- seq_aipe_sm(rep(3, 10), width = 0.4)
  printed <- paste(format(out), collapse = "\n")
  expect_match(printed, "Stopping rule met: FALSE", fixed = TRUE)
  expect_match(printed, "not yet usable")
  expect_match(printed, "not finite")
})
