# Tests for seq_fwd_norm_mean_optimal(), the known variance sample size for
# a fixed width interval for a normal mean (calculate_c() before 1.0.0), and
# for the functions held out of 1.0.0.

test_that("seq_fwd_norm_mean_optimal returns the smallest sufficient whole sample size", {
  # C = z^2 sigma^2 / d^2 with z_{0.025} = 1.959964.
  # variance 1, d = 0.75: (1.959964 / 0.75)^2 = 6.829, so 7.
  expect_identical(seq_fwd_norm_mean_optimal(variance = 1, d = 0.75), 7L)
  # variance 64, d = 2: 3.841459 * 64 / 4 = 61.46, so 62.
  expect_identical(seq_fwd_norm_mean_optimal(variance = 64, d = 2), 62L)
  # 99%: z_{0.005} = 2.575829; 6.634897 * 64 / 4 = 106.16, so 107.
  expect_identical(
    seq_fwd_norm_mean_optimal(variance = 64, d = 2, conf_level = 0.99), 107L)
  # It is the variance that is supplied, not the standard deviation: four
  # times the variance needs four times the (unrounded) sample size.
  expect_identical(seq_fwd_norm_mean_optimal(variance = 256, d = 2), 246L)
  expect_identical(names(formals(seq_fwd_norm_mean_optimal)),
                   c("variance", "d", "conf_level"))
})

test_that("the known variance interval reaches the level exactly at the returned n", {
  # Coverage of xbar +/- d with known sigma is 2 pnorm(d sqrt(n) / sigma) - 1.
  coverage <- function(n, variance, d) {
    2 * stats::pnorm(d * sqrt(n) / sqrt(variance)) - 1
  }
  grid <- expand.grid(variance = c(0.3, 1, 64, 1500), d = c(0.05, 0.4, 2),
                      conf_level = c(0.8, 0.95, 0.999))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    n <- seq_fwd_norm_mean_optimal(g$variance, g$d, g$conf_level)
    expect_type(n, "integer")
    expect_gte(coverage(n, g$variance, g$d), g$conf_level)
    if (n > 1) {
      expect_lt(coverage(n - 1, g$variance, g$d), g$conf_level)
    }
  }
})

test_that("the returned n has the defining property where C is next to a whole number", {
  # Choosing the variance as n d^2 / z^2 puts C = z^2 variance / d^2 within
  # rounding error of the whole number n, on one side or the other. The
  # coverage formula of the help page, evaluated with pnorm() alone, then
  # decides: the returned sample size reaches the level and the one below
  # it does not. (Rounding C up alone returns n + 1 in many of these cases
  # although n already reaches the level.)
  coverage <- function(n, variance, d) {
    2 * stats::pnorm(d * sqrt(n) / sqrt(variance)) - 1
  }
  for (case in list(c(0.95, 0.5), c(0.99, 2), c(0.80, 0.1))) {
    conf_level <- case[1]; d <- case[2]
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    for (n in 1:500) {
      variance <- n * d^2 / z^2
      got <- seq_fwd_norm_mean_optimal(variance, d, conf_level)
      ok <- coverage(got, variance, d) >= conf_level &&
        (got == 1L || coverage(got - 1L, variance, d) < conf_level) &&
        abs(got - n) <= 1L
      if (!ok) {
        fail(sprintf("n = %d, conf_level = %g: returned %d", n, conf_level,
                     got))
      }
    }
  }
  succeed()
})

test_that("the returned n changes exactly where the coverage formula crosses the level", {
  # For a whole number n, bisection over doubles finds two adjacent values
  # of the variance, lo < hi, such that n observations reach the level at lo
  # and miss it at hi (pnorm() alone). By the definition in the help page the
  # answer is n at lo and n + 1 at hi; the closed form C cannot tell the two
  # apart.
  coverage <- function(n, variance, d) {
    2 * stats::pnorm(d * sqrt(n) / sqrt(variance)) - 1
  }
  for (conf_level in c(0.80, 0.95, 0.99)) for (d in c(0.1, 2)) {
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    for (n in c(1, 2, 7, 62, 500, 12345)) {
      lo <- n * d^2 / z^2 * (1 - 1e-9)
      hi <- n * d^2 / z^2 * (1 + 1e-9)
      expect_gte(coverage(n, lo, d), conf_level)
      expect_lt(coverage(n, hi, d), conf_level)
      repeat {
        mid <- (lo + hi) / 2
        if (mid == lo || mid == hi) break
        if (coverage(n, mid, d) >= conf_level) lo <- mid else hi <- mid
      }
      expect_identical(seq_fwd_norm_mean_optimal(lo, d, conf_level),
                       as.integer(n))
      expect_identical(seq_fwd_norm_mean_optimal(hi, d, conf_level),
                       as.integer(n) + 1L)
    }
  }
})

test_that("seq_fwd_norm_mean_optimal never returns a sample size below 1", {
  # A very wide target needs one observation, not none. C here is
  # 3.84 * 1e200 / 1e320 = 3.8e-120 and 3.84 * 1e300 / 1e310 = 3.8e-10;
  # computing d^2 on its own overflowed and the function returned 0.
  expect_identical(seq_fwd_norm_mean_optimal(variance = 1e200, d = 1e160), 1L)
  expect_identical(seq_fwd_norm_mean_optimal(variance = 1e300, d = 1e155), 1L)
  expect_identical(seq_fwd_norm_mean_optimal(variance = 1, d = 100), 1L)
  # A confidence level so low that z is 0 in floating point.
  expect_identical(seq_fwd_norm_mean_optimal(4, 1, conf_level = 1e-17), 1L)
  # The scale of measurement does not matter: variance 64 s^2 with d = 2 s
  # is the 62 of the first test at any s.
  for (s in c(1e-150, 1e-3, 1e3, 1e150)) {
    expect_identical(seq_fwd_norm_mean_optimal(variance = 64 * s^2, d = 2 * s),
                     62L)
  }
  # Nor do the largest representable magnitudes: a standard deviation of
  # 1e154 with d = 1e153 is sigma / d = 10, so C = 3.841459 * 100 = 384.15
  # and the answer is 385 (z^2 * variance overflows if formed on its own).
  expect_identical(seq_fwd_norm_mean_optimal(variance = 1e308, d = 1e153),
                   385L)
})

test_that("seq_fwd_norm_mean_optimal checks every argument", {
  expect_error(seq_fwd_norm_mean_optimal(d = 0.5), "'variance' must be supplied")
  expect_error(seq_fwd_norm_mean_optimal(variance = -1, d = 0.5), "'variance'")
  expect_error(seq_fwd_norm_mean_optimal(variance = 0, d = 0.5), "'variance'")
  # Before 1.0.0 a missing value surfaced as R's own
  # "missing value where TRUE/FALSE needed".
  expect_error(seq_fwd_norm_mean_optimal(variance = NA_real_, d = 0.5),
               "'variance' must be a single finite number")
  expect_error(seq_fwd_norm_mean_optimal(variance = c(1, 2), d = 0.5), "'variance'")
  expect_error(seq_fwd_norm_mean_optimal(variance = "1", d = 0.5), "'variance'")
  expect_error(seq_fwd_norm_mean_optimal(variance = 1), "'d' must be supplied")
  expect_error(seq_fwd_norm_mean_optimal(variance = 1, d = 0), "'d'")
  expect_error(seq_fwd_norm_mean_optimal(variance = 1, d = NA_real_), "'d'")
  expect_error(seq_fwd_norm_mean_optimal(variance = 1, d = Inf), "'d'")
  expect_error(seq_fwd_norm_mean_optimal(variance = 1, d = 0.5, conf_level = 1),
               "'conf_level'")
  expect_error(seq_fwd_norm_mean_optimal(variance = 1, d = 0.5, conf_level = NA_real_),
               "'conf_level'")
  expect_error(seq_fwd_norm_mean_optimal(variance = 1, d = 1e-6),
               "implied sample size .* is too large to represent")
  expect_error(seq_fwd_norm_mean_optimal(variance = 1, d = 1e-6),
               "Check that 'd' is in the units of the data and 'variance' in their square")
})

test_that("seq_fwd_norm_mean_optimal at the largest representable sample size is that size or a clear error", {
  # C can round up to 2147483647, the largest integer, while the coverage
  # at that size is still below the level by rounding error; the answer is
  # then 2147483648, which cannot be represented, and that second step must
  # stop with the same advice as the first. The inputs are the values of d
  # within a few doubles of the boundary (C is computed here only to find
  # them, never as an expected value). Which of them need the second step
  # depends on how the platform rounds pnorm() and qnorm(), so every
  # outcome is checked and no particular count is required: the result is
  # 2147483647, or the error with the advice. (On the platform this was
  # written on, 5 inputs in 100 stop.)
  withr::local_seed(113)
  largest <- .Machine$integer.max
  sizes <- integer(0)
  messages <- character(0)
  for (i in 1:150) {
    conf_level <- stats::runif(1, 0.6, 0.999)
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    d <- z / sqrt(largest)
    for (step in 0:12) {
      C <- (z / d)^2
      if (C <= largest && C > largest - 0.5) {
        res <- tryCatch(seq_fwd_norm_mean_optimal(variance = 1, d = d,
                                                  conf_level = conf_level),
                        error = conditionMessage)
        if (is.character(res)) {
          messages <- c(messages, res)
        } else {
          sizes <- c(sizes, res)
        }
      }
      d <- d * (1 + .Machine$double.eps)
    }
  }
  expect_gt(length(sizes) + length(messages), 1000L)
  expect_true(all(sizes == largest))
  # Upper-tail evaluation of the normal quantile can put the intermediate
  # a fraction above the largest integer; its displayed value then rounds
  # to 2147483647 although the unrounded number is not representable as an
  # integer. Both paths must retain the same advice.
  expected_messages <- paste0(
    "The implied sample size (", c("2147483647", "2147483648"),
    ") is too large to represent. ",
    "Check that 'd' is in the units of the data and 'variance' in their ",
    "square; the target as given calls for more than two billion ",
    "observations.")
  expect_true(all(messages %in% expected_messages))
})

test_that("the functions held out of 1.0.0 are gone, with no alias for the renamed one", {
  exported <- getNamespaceExports("SMSD")
  gone <- c("calculate_c", "prob_gn", "prob_hn", "prob_kn", "seven_1",
            "ts_mr_dfr_mean", "ts_mr_dfr_mean_pilot")
  expect_false(any(gone %in% exported))
  expect_true("seq_fwd_norm_mean_optimal" %in% exported)
  ns <- asNamespace("SMSD")
  for (helper in c(".chi_prob", ".chi_prob_2", ".chi_prob_3")) {
    expect_false(exists(helper, envir = ns, inherits = FALSE))
  }
})
