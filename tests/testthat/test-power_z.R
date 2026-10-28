# Tests for power_z(): the one sample, upper tailed z test with a known
# standard deviation. Expected values are closed forms evaluated from
# printed normal quantiles (z_{0.95} = 1.644854, z_{0.80} = 0.841621), the
# defining inequality checked directly with pnorm(), or a simulation of the
# test itself.

test_that("power_z returns the Type II error rate at a given sample size", {
  # Phi(z_{0.95} - sqrt(16) * (1 - 0) / 2) = Phi(1.644854 - 2) = Phi(-0.355146).
  out <- power_z(mu_0 = 0, mu_1 = 1, sigma = 2, n = 16, method = "type_2")
  expect_equal(out, 0.3612400, tolerance = 1e-6)
  expect_equal(out, stats::pnorm(1.644854 - 2), tolerance = 1e-6)
  # Only the standardized effect matters: (104 - 100) / 8 = 0.5 as well.
  expect_equal(power_z(100, 104, 8, n = 16, method = "type_2"), out)
  # conf_level sets the one-sided level: 0.99 -> z = 2.326348,
  # Phi(2.326348 - 2) = Phi(0.326348) = 0.627920.
  expect_equal(power_z(0, 1, 2, conf_level = 0.99, n = 16, method = "type_2"),
               0.627920, tolerance = 1e-5)
  # More observations, smaller Type II error; n = 1 is allowed.
  expect_lt(power_z(0, 1, 2, n = 64, method = "type_2"), out)
  expect_equal(power_z(0, 1, 2, n = 1, method = "type_2"),
               stats::pnorm(1.644854 - 0.5), tolerance = 1e-6)
  # 1.644854 is the printed z_{0.95} that the lines above rely on
  # (documented number, not a check of the code).
  expect_equal(stats::qnorm(0.95), 1.644854, tolerance = 1e-6)
})

test_that("the Type II error rate agrees with a simulation of the test", {
  withr::local_seed(113)
  reps <- 20000; n <- 16
  # Sample means under the alternative: N(1, 2^2 / 16). The test accepts H0
  # when sqrt(n) (xbar - 0) / 2 <= z_{0.95}.
  xbar <- rnorm(reps, mean = 1, sd = 2 / sqrt(n))
  accepted <- mean(sqrt(n) * xbar / 2 <= stats::qnorm(0.95))
  expected <- power_z(0, 1, 2, n = n, method = "type_2")
  # Tolerance: four Monte Carlo standard errors (about 0.014).
  expect_lt(abs(accepted - expected),
            4 * sqrt(expected * (1 - expected) / reps))
})

test_that("power_z returns the smallest whole sample size for a target Type II error", {
  # Unrounded: ((z_{0.95} + z_{0.80}) / 0.5)^2 = (2.486475 / 0.5)^2 = 24.73,
  # so 25 (the textbook value for a one-sided 0.05 test, power 0.80, and a
  # standardized effect of one half).
  out <- power_z(mu_0 = 0, mu_1 = 1, sigma = 2, beta = 0.20,
                 method = "sample_size")
  expect_identical(out, 25L)
  # beta = 0.361: unrounded ((1.644854 + 0.355735) * 2)^2 = 16.01, so 17.
  # (Before 1.0.0 the unrounded 16.01 was returned.)
  expect_identical(power_z(0, 1, 2, beta = 0.361, method = "sample_size"),
                   17L)
  # The defining property, checked with pnorm() alone: the returned n meets
  # the target and n - 1 does not.
  type_2 <- function(n, effect, alpha) {
    stats::pnorm(stats::qnorm(1 - alpha) - sqrt(n) * effect)
  }
  grid <- expand.grid(effect = c(0.1, 0.27, 0.5, 1.3),
                      beta = c(0.05, 0.2, 0.5, 0.8),
                      conf_level = c(0.9, 0.95, 0.999))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    n <- power_z(10, 10 + g$effect * 3, 3, conf_level = g$conf_level,
                 beta = g$beta, method = "sample_size")
    expect_type(n, "integer")
    expect_lte(type_2(n, g$effect, 1 - g$conf_level), g$beta)
    if (n > 1) {
      expect_gt(type_2(n - 1, g$effect, 1 - g$conf_level), g$beta)
    }
  }
})

test_that("the sample size inverts the Type II error rate that power_z reports", {
  # If n observations give a Type II error rate of b, then the smallest
  # sample size whose rate is at most b is n: the rate falls as n grows, so
  # n - 1 misses the target, and n meets it by definition. The closed form
  # for the sample size is then a whole number up to rounding error, and
  # rounding it up used to return n + 1 in about 40% of cases (61 for the
  # help page's scenario at n = 60).
  b_60 <- power_z(mu_0 = 100, mu_1 = 104, sigma = 15, n = 60)
  expect_identical(power_z(100, 104, 15, beta = b_60), 60L)
  for (effect in c(0.10, 0.25, 4 / 15)) {
    back <- vapply(1:400, function(n) {
      power_z(0, effect, 1, beta = power_z(0, effect, 1, n = n))
    }, integer(1))
    expect_identical(back, 1:400)
  }
  # A target the smallest step below the rate at n is not met by n, so the
  # answer is n + 1 (the rate at n + 1 is lower by far more than that step).
  up <- vapply(1:400, function(n) {
    b <- power_z(0, 0.25, 1, n = n)
    power_z(0, 0.25, 1, beta = b * (1 - .Machine$double.eps))
  }, integer(1))
  expect_identical(up, 2:401)
  # At other levels too.
  for (conf_level in c(0.90, 0.99)) {
    back <- vapply(1:200, function(n) {
      b <- power_z(0, 0.2, 1, conf_level = conf_level, n = n)
      power_z(0, 0.2, 1, conf_level = conf_level, beta = b)
    }, integer(1))
    expect_identical(back, 1:200)
  }
})

test_that("power_z keeps the quantile of a very small Type II error rate", {
  # z for an upper tail of 1e-17 is 8.493793 and for 1e-30 is 11.464025
  # (computed in the upper tail; 1 - 1e-17 is 1 in floating point, whose
  # quantile is infinite). With z_{0.95} = 1.644854 and an effect of 0.5:
  # ((1.644854 + 8.493793) / 0.5)^2 = 411.17, so 412;
  # ((1.644854 + 11.464025) / 0.5)^2 = 687.37, so 688.
  # (Before 1.0.0 both stopped with "(Inf) is too large to represent".)
  expect_identical(power_z(0, 0.5, 1, beta = 1e-17), 412L)
  expect_identical(power_z(0, 0.5, 1, beta = 1e-30), 688L)
  # The defining property, with pnorm() alone.
  type_2 <- function(n) stats::pnorm(stats::qnorm(0.95) - sqrt(n) * 0.5)
  expect_lte(type_2(412), 1e-17)
  expect_gt(type_2(411), 1e-17)
  expect_lte(type_2(688), 1e-30)
  expect_gt(type_2(687), 1e-30)
  # Just above the point where 1 - beta rounds to 1 the old code was one
  # short: 3e-16 gave 379 where 380 is needed (z for an upper tail of 3e-16
  # is 8.089317, and ((1.644854 + 8.089317) / 0.5)^2 = 379.02).
  n_small <- power_z(0, 0.5, 1, beta = 3e-16)
  expect_lte(type_2(n_small), 3e-16)
  expect_gt(type_2(n_small - 1L), 3e-16)
  expect_identical(n_small, 380L)
})

test_that("power_z returns 1 when every sample size meets the target", {
  # beta >= conf_level asks for a power no greater than the Type I error
  # rate. The Type II error at n = 1 is Phi(1.644854 - 0.5) = 0.874, already
  # below these targets. (Before 1.0.0, squaring the negative sum of
  # quantiles returned 17.2 for beta = 0.9999 and 0.045 for beta = 0.96.)
  expect_identical(power_z(0, 1, 2, beta = 0.9999), 1L)
  expect_identical(power_z(0, 1, 2, beta = 0.96), 1L)
  expect_identical(power_z(0, 1, 2, beta = 0.95), 1L)
  # Documented number (not a check of the code): the rate at n = 1 that the
  # comment above quotes is below these targets.
  expect_lt(stats::pnorm(1.644854 - 0.5), 0.95)
  # Just inside: beta = 0.94 still needs more than nothing but only n = 1.
  expect_identical(power_z(0, 1, 2, beta = 0.94), 1L)
})

test_that("power_z tabulates Type II error rates over whole sample sizes", {
  out <- power_z(mu_0 = 0, mu_1 = 1, sigma = 2, n_min = 5, n_max = 50,
                 n_values = 10, method = "type_2_seq")
  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("n", "type_2_error"))
  expect_identical(out$n, seq(5L, 50L, by = 5L))
  expect_equal(out$type_2_error,
               stats::pnorm(1.644854 - sqrt(seq(5, 50, by = 5)) * 0.5),
               tolerance = 1e-6)
  expect_equal(out$type_2_error[out$n == 15],
               power_z(0, 1, 2, n = 15, method = "type_2"))
  # A grid that is not whole is taken to the nearest whole number, a half
  # going up, and duplicates are dropped:
  # 1, 1.5, 2, 2.5, 3 -> 1, 2, 2, 3, 3 -> 1, 2, 3.
  small <- power_z(0, 1, 2, n_min = 1, n_max = 3, n_values = 5)
  expect_identical(small$n, 1:3)
  # Nearest, not truncated: 1, 3.33, 5.67, 8 -> 1, 3, 6, 8.
  expect_identical(power_z(0, 1, 2, n_min = 1, n_max = 8, n_values = 4)$n,
                   c(1L, 3L, 6L, 8L))
  # A half goes up every time (round() would send 2.5 to 2 but 3.5 to 4):
  # 1, 2.5, 4 -> 1, 3, 4;  2, 3.5, 5 -> 2, 4, 5;  4, 5.5, 7 -> 4, 6, 7.
  expect_identical(power_z(0, 1, 2, n_min = 1, n_max = 4, n_values = 3)$n,
                   c(1L, 3L, 4L))
  expect_identical(power_z(0, 1, 2, n_min = 2, n_max = 5, n_values = 3)$n,
                   c(2L, 4L, 5L))
  expect_identical(power_z(0, 1, 2, n_min = 4, n_max = 7, n_values = 3)$n,
                   c(4L, 6L, 7L))
  # Tenths from 1 to 3: 21 points, the halves 1.5 and 2.5 go up, so one 1,
  # ten 2s, and ten 3s before duplicates are dropped; and a long grid of
  # whole numbers comes back as itself.
  expect_identical(power_z(0, 1, 2, n_min = 1, n_max = 3, n_values = 21)$n,
                   1:3)
  expect_identical(power_z(0, 1, 2, n_min = 10, n_max = 5000,
                           n_values = 4991)$n, 10:5000)
  # The table's rates are the rates at those whole sample sizes.
  expect_equal(power_z(0, 1, 2, n_min = 1, n_max = 8,
                       n_values = 4)$type_2_error,
               stats::pnorm(1.644854 - sqrt(c(1, 3, 6, 8)) * 0.5),
               tolerance = 1e-6)
  # A single value, and n_min = n_max.
  expect_identical(power_z(0, 1, 2, n_min = 7, n_max = 9, n_values = 1)$n, 7L)
  expect_identical(power_z(0, 1, 2, n_min = 7, n_max = 7, n_values = 4)$n, 7L)
})

test_that("the table costs what its rows cost, not what n_values asks for", {
  # No table has more rows than there are whole sample sizes from n_min to
  # n_max, and from that many points on, the points are no more than 1
  # apart, so the table is every sample size in the range. A grid built with
  # n_values points and reduced afterward would, at the largest n_values the
  # argument check admits, take 2^31 - 1 doubles (16 GiB) for a table of 10
  # rows; the time limit makes that a failure here and not a session that
  # swaps for minutes.
  setTimeLimit(elapsed = 30, transient = TRUE)
  withr::defer(setTimeLimit(elapsed = Inf))
  expect_identical(.power_z_grid(1L, 10L, .Machine$integer.max), 1:10)
  huge <- power_z(100, 104, 15, n_min = 1, n_max = 10,
                  n_values = .Machine$integer.max)
  expect_identical(huge, power_z(100, 104, 15, n_min = 1, n_max = 10,
                                 n_values = 10))
  expect_identical(huge$n, 1:10)
  # At and above the number of whole sample sizes the grid is all of them;
  # one point fewer and it is not (9 points from 1 to 10 are 1.125 apart:
  # 1, 2.125, 3.25, 4.375, 5.5, 6.625, 7.75, 8.875, 10, so 5 is missing).
  for (extra in c(0L, 1L, 7L, 1000L)) {
    expect_identical(.power_z_grid(1L, 10L, 10L + extra), 1:10)
    expect_identical(.power_z_grid(37L, 61L, 25L + extra), 37:61)
  }
  expect_identical(.power_z_grid(1L, 10L, 9L),
                   c(1L, 2L, 3L, 4L, 6L, 7L, 8L, 9L, 10L))
  # The widest range the argument checks admit does not overflow the count.
  expect_identical(.power_z_grid(1L, .Machine$integer.max, 2L),
                   c(1L, .Machine$integer.max))
})

test_that("power_z chooses the method from the arguments supplied", {
  expect_identical(power_z(0, 1, 2, beta = 0.2),
                   power_z(0, 1, 2, beta = 0.2, method = "sample_size"))
  expect_identical(power_z(0, 1, 2, n = 16),
                   power_z(0, 1, 2, n = 16, method = "type_2"))
  expect_identical(power_z(0, 1, 2, n_min = 5, n_max = 50, n_values = 10),
                   power_z(0, 1, 2, n_min = 5, n_max = 50, n_values = 10,
                           method = "type_2_seq"))
  expect_error(power_z(0, 1, 2), "Nothing to compute")
  expect_error(power_z(0, 1, 2, n = 16, beta = 0.2), "More than one result")
  # With method given, arguments it does not use are ignored.
  expect_identical(power_z(0, 1, 2, n = 16, beta = 0.2, method = "type_2"),
                   power_z(0, 1, 2, n = 16))
  # A partial grid selects the table and then says what is missing.
  expect_error(power_z(0, 1, 2, n_min = 5), "'n_min', 'n_max', and 'n_values'")
})

test_that("power_z checks 'method' itself and names it", {
  # Before 1.0.0 these were R's own match.arg() errors about 'arg'. The
  # message is the package's shared one for a choice among strings.
  expect_error(power_z(0, 1, 2, n = 16, method = "nonsense"),
               "'method' must be one of \"sample_size\", \"type_2\", \"type_2_seq\"; got \"nonsense\"\\.")
  expect_error(power_z(0, 1, 2, n = 16, method = c("type_2", "sample_size")),
               "'method' must be one of .* got a character vector of length 2")
  expect_error(power_z(0, 1, 2, n = 16, method = 1), "'method' must be one of")
  expect_error(power_z(0, 1, 2, n = 16, method = NA), "'method' must be one of")
  expect_error(power_z(0, 1, 2, n = 16, method = NA_character_),
               "'method' must be one of")
  expect_error(power_z(0, 1, 2, n = 16, method = TRUE), "'method' must be one of")
  # "type" abbreviates two of the choices, so it names neither.
  expect_error(power_z(0, 1, 2, n = 16, method = "type"),
               "'method' must be one of")
  expect_error(power_z(0, 1, 2, n = 16, method = ""), "'method' must be one of")
  # A unique abbreviation is accepted; "type_2" is itself, not the table.
  expect_identical(power_z(0, 1, 2, beta = 0.2, n = 16, method = "sample"), 25L)
  expect_identical(power_z(0, 1, 2, beta = 0.2, n = 16, method = "type_2"),
                   power_z(0, 1, 2, n = 16))
  expect_s3_class(power_z(0, 1, 2, n_min = 5, n_max = 50, n_values = 10,
                          method = "type_2_s"), "data.frame")
  # NULL, and the default vector of choices passed along by a wrapper, mean
  # that no method was named: it is chosen from the arguments.
  expect_identical(power_z(0, 1, 2, n = 16, method = NULL),
                   power_z(0, 1, 2, n = 16))
  wrapper <- function(..., method = c("sample_size", "type_2", "type_2_seq")) {
    power_z(..., method = method)
  }
  expect_identical(wrapper(0, 1, 2, n = 16), power_z(0, 1, 2, n = 16))
  expect_identical(wrapper(0, 1, 2, beta = 0.2), 25L)
  expect_identical(wrapper(0, 1, 2, n = 16, beta = 0.2, method = "type_2"),
                   power_z(0, 1, 2, n = 16))
  expect_error(power_z(0, 1, 2, method = NULL), "Nothing to compute")
  # The choices shown in the usage are the ones the function accepts.
  expect_identical(eval(formals(power_z)$method),
                   c("sample_size", "type_2", "type_2_seq"))
})

test_that("power_z checks every argument and says what to do", {
  expect_error(power_z(0, 1), "must all be supplied")
  # Before 1.0.0: mu_1 < mu_0 returned 0.9999 where the lower tailed answer
  # is 0.361, and mu_1 = mu_0 returned Inf.
  expect_error(power_z(0, -1, 2, n = 16), "'mu_1' must be greater than 'mu_0'")
  expect_error(power_z(0, -1, 2, n = 16), "pass -mu_0 and -mu_1")
  expect_error(power_z(0, 0, 2, beta = 0.2), "'mu_1' must be greater than 'mu_0'")
  # The reflection the message recommends gives the lower tailed answer.
  expect_equal(power_z(0, 1, 2, n = 16), 0.3612400, tolerance = 1e-6)
  expect_error(power_z(NA_real_, 1, 2, n = 16), "'mu_0'")
  expect_error(power_z(0, Inf, 2, n = 16), "'mu_1'")
  expect_error(power_z(0, "a", 2, n = 16), "'mu_1'")
  expect_error(power_z(0, 1, -2, n = 16), "'sigma'")
  expect_error(power_z(0, 1, 0, n = 16), "'sigma'")
  expect_error(power_z(0, 1, NA_real_, n = 16), "'sigma'")
  expect_error(power_z(0, 1, 2, conf_level = 1.2, n = 16), "'conf_level'")
  expect_error(power_z(0, 1, 2, conf_level = NA_real_, n = 16), "'conf_level'")
  expect_error(power_z(0, 1, 2, beta = 1.5), "'beta'")
  expect_error(power_z(0, 1, 2, beta = 0), "'beta'")
  expect_error(power_z(0, 1, 2, beta = 1), "'beta'")
  expect_error(power_z(0, 1, 2, beta = NA_real_), "'beta'")
  expect_error(power_z(0, 1, 2, n = 0), "'n'")
  expect_error(power_z(0, 1, 2, n = 16.5), "'n' must be a single whole number")
  expect_error(power_z(0, 1, 2, n = NA), "'n'")
  expect_error(power_z(0, 1, 2, n = c(16, 25)), "'n'")
  expect_error(power_z(0, 1, 2, method = "sample_size"), "needs 'beta'")
  expect_error(power_z(0, 1, 2, method = "type_2"), "needs 'n'")
  expect_error(power_z(0, 1, 2, method = "type_2_seq"), "needs 'n_min'")
  expect_error(power_z(0, 1, 2, n_min = 0, n_max = 5, n_values = 3), "'n_min'")
  expect_error(power_z(0, 1, 2, n_min = 10, n_max = 5, n_values = 3),
               "'n_max' must be a single whole number of at least 10")
  expect_error(power_z(0, 1, 2, n_min = 5, n_max = 50, n_values = 0), "'n_values'")
  expect_error(power_z(0, 1, 2, n_min = 5, n_max = 50, n_values = 2.5), "'n_values'")
  # A sample size too large to represent is an error, not NA, and the
  # message says what to look at.
  expect_error(power_z(0, 1e-9, 1, beta = 0.2),
               "implied sample size .* is too large to represent")
  expect_error(power_z(0, 1e-9, 1, beta = 0.2),
               "Check 'mu_0', 'mu_1', and 'sigma'; a difference this small")
})

test_that("power_z at the largest representable sample size is that size or a clear error", {
  # The closed form can round up to 2147483647, the largest integer, while
  # the Type II error rate at that size is still above the target by
  # rounding error; the answer is then 2147483648, which cannot be
  # represented, and that second step must stop with the same advice as the
  # first. The inputs are the effects within a few doubles of the boundary
  # (the closed form is used here only to find them, never as an expected
  # value). Which of them need the second step depends on how the platform
  # rounds pnorm() and qnorm(), so every outcome is checked and no
  # particular count is required: the result is 2147483647, or the error
  # with the advice. (On the platform this was written on, 4 inputs in 100
  # stop.)
  withr::local_seed(113)
  largest <- .Machine$integer.max
  sizes <- integer(0)
  messages <- character(0)
  for (i in 1:150) {
    conf_level <- stats::runif(1, 0.6, 0.999)
    beta <- stats::runif(1, 0.01, 0.5)
    z_sum <- stats::qnorm(conf_level) + stats::qnorm(beta, lower.tail = FALSE)
    effect <- z_sum / sqrt(largest)
    for (step in 0:12) {
      closed_form <- (z_sum / effect)^2
      if (closed_form <= largest && closed_form > largest - 0.5) {
        res <- tryCatch(power_z(0, effect, 1, conf_level = conf_level,
                                beta = beta), error = conditionMessage)
        if (is.character(res)) {
          messages <- c(messages, res)
        } else {
          sizes <- c(sizes, res)
        }
      }
      effect <- effect * (1 + .Machine$double.eps)
    }
  }
  expect_gt(length(sizes) + length(messages), 1000L)
  expect_true(all(sizes == largest))
  expect_true(all(messages == paste0(
    "The implied sample size (2147483648) is too large to represent. ",
    "Check 'mu_0', 'mu_1', and 'sigma'; a difference this small next to ",
    "the standard deviation calls for more than two billion observations.")))
})
