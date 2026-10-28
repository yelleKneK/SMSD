# Tests for the purely sequential bounded risk point estimator of a normal
# mean, seq_br_norm_mean(), and its pilot. The loss is A |xbar - mu|^k. The
# expected values come from numerical integration of that loss against the
# normal distribution of the sample mean, or from hand arithmetic; none is
# the package's closed form typed again.

# Risk of the mean of n observations with variance v, by quadrature.
br_risk <- function(A, k, v, n) {
  stats::integrate(function(e) A * abs(e)^k * stats::dnorm(e, 0, sqrt(v / n)),
                   -Inf, Inf, rel.tol = 1e-10)$value
}

# ---- criterion and risk -----------------------------------------------------

test_that("optimal_n is the smallest n whose risk is at most w, and risk is the risk", {
  set.seed(113)
  for (cfg in list(c(A = 2, k = 1, w = 0.5, v = 9), c(A = 2, k = 2, w = 0.5, v = 8.9),
                   c(A = 1, k = 1, w = 0.25, v = 1.44), c(A = 30, k = 2, w = 4, v = 2.7))) {
    A <- cfg[["A"]]; k <- cfg[["k"]]; w <- cfg[["w"]]; v <- cfg[["v"]]
    x <- seqc_scale_var(rnorm(20, 10, 1), v)
    out <- seq_br_norm_mean(x, A = A, k = k, w = w)
    n_opt <- out$value[out$term == "optimal_n"]
    expect_lte(br_risk(A, k, v, n_opt), w)
    expect_gt(br_risk(A, k, v, n_opt - 1), w)
    # The risk row is the risk at the current n = 20 with sigma^2 = v.
    expect_equal(out$value[out$term == "risk"], br_risk(A, k, v, 20),
                 tolerance = 1e-7)
    expect_identical(attr(out, "is_satisfied"), 20 >= n_opt)
    # Meeting the rule is the same thing as an estimated risk of at most w.
    expect_identical(attr(out, "is_satisfied"),
                     out$value[out$term == "risk"] <= w)
  }
})

test_that("hand computed criteria", {
  set.seed(113)
  # Squared error, A = 2, w = 0.5, variance 8.9: n* = A v / w = 35.6 -> 36.
  x <- seqc_scale_var(rnorm(30, 10, 1), 8.9)
  out <- seq_br_norm_mean(x, A = 2, k = 2, w = 0.5)
  expect_equal(out$value[out$term == "optimal_n"], 36)
  expect_false(attr(out, "is_satisfied"))                 # 30 < 36
  expect_equal(out$value[out$term == "risk"], 2 * 8.9 / 30)
  # Absolute error, A = 1, w = 0.25, variance 1.44 (sd 1.2):
  # E|xbar - mu| = sd sqrt(2 / (pi n)), so n* = (2 / pi) (1.2 / 0.25)^2
  # = 14.67 -> 15.
  y <- seqc_scale_var(rnorm(15, 4.8, 1), 1.44)
  out <- seq_br_norm_mean(y, A = 1, k = 1, w = 0.25)
  expect_equal(out$value[out$term == "optimal_n"], 15)
  expect_true(attr(out, "is_satisfied"))                  # 15 >= 15
  expect_equal(out$value[out$term == "risk"], 1.2 * sqrt(2 / (pi * 15)))
})

test_that("the risk constant is in closed form: exactly A for squared error", {
  # B = A E|Z|^k for a standard normal Z. E|Z|^k by numerical integration of
  # its definition (no gamma function): 1 for k = 2, sqrt(2 / pi) for k = 1.
  abs_moment <- function(k) {
    2 * stats::integrate(function(z) z^k * stats::dnorm(z), 0, Inf,
                         rel.tol = 1e-12)$value
  }
  expect_equal(abs_moment(2), 1, tolerance = 1e-9)
  expect_equal(.br_constant(3, 1), 3 * abs_moment(1), tolerance = 1e-9)
  expect_equal(.br_constant(3, 1), 3 * 0.7978846, tolerance = 1e-7)
  # k = 2: B is A itself, with no rounding error at all. The general
  # expression A 2^(k / 2) Gamma((k + 1) / 2) / sqrt(pi), used before 1.0.0,
  # evaluated to A (1 + 2.2e-16).
  for (A in c(1, 2, 3, 5, 7.25, 10, 100, 0.1, 1e-3)) {
    expect_identical(.br_constant(A, 2), A)
    expect_identical(.seq_br_scale(A, 2, 4), A / 4)
  }
  expect_identical(.seq_br_scale(3, 1, 0.5), (3 * sqrt(2 / pi) / 0.5)^2)
})

test_that("whole number criteria and the risk at the bound are exact (k = 2)", {
  # Five observations with mean 5 and variance exactly 7 (deviations -3, -1,
  # -1, 1, 4; squares sum to 28). With A = 5 and w = 7 the criterion
  # A S^2 / w is exactly 5 and the estimated risk at n = 5, A S^2 / n, is
  # exactly 7: the rule is met at n = 5 and the risk row equals the bound.
  # With the general gamma function expression for B the risk row was
  # 7.0000000000000009, above the bound beside a met rule, and (before the
  # package's rounding up allowed for representation error) optimal_n was 6:
  # one observation more than the exact answer.
  x <- c(2, 4, 4, 6, 9)
  expect_identical(stats::var(x), 7)
  out <- seq_br_norm_mean(x, A = 5, k = 2, w = 7)
  expect_identical(attr(out, "criterion_n"), 5L)
  expect_true(attr(out, "is_satisfied"))
  expect_identical(out$value[out$term == "risk"], 7)
  expect_lte(out$value[out$term == "risk"], 7)
  # The same case with no help from the package's rounding: the unrounded
  # criterion (B / w) S^2 is exactly 5, so a plain ceiling gives 5 as well.
  # (From the general expression it was 5.000000000000001, whose ceiling is
  # 6: the case in which the two forms of B differ by one observation.)
  expect_identical(ceiling(.seq_br_scale(5, 2, 7) * stats::var(x)), 5)
  # Other whole number cases, every step exact in floating point:
  # (A / w) S^2 = (2 / 0.5) 7 = 28, (1 / 1) 7 = 7, (3 / 2) 7 = 10.5 -> 11,
  # (10 / 2) 7 = 35; risk A S^2 / n = 2 * 7 / 5 = 2.8.
  expect_identical(attr(seq_br_norm_mean(x, A = 2, k = 2, w = 0.5),
                        "criterion_n"), 28L)
  expect_identical(attr(seq_br_norm_mean(x, A = 1, k = 2, w = 1),
                        "criterion_n"), 7L)
  expect_identical(attr(seq_br_norm_mean(x, A = 3, k = 2, w = 2),
                        "criterion_n"), 11L)
  expect_identical(attr(seq_br_norm_mean(x, A = 10, k = 2, w = 2),
                        "criterion_n"), 35L)
  expect_identical(seq_br_norm_mean(x, A = 2, k = 2, w = 0.5)$value[3], 2.8)
  # A / w a perfect square r^2: the pilot ceiling(sqrt(A / w)) is r itself,
  # not r + 1 (which the general expression gave for 138 of these 397 roots
  # before the package's rounding up allowed for representation error).
  roots <- 4:400
  got <- vapply(roots, function(r) {
    seq_br_norm_mean_pilot(A = r^2, k = 2, w = 1, m_0 = 2)
  }, integer(1))
  expect_identical(got, roots)
  expect_identical(seq_br_norm_mean_pilot(A = 100, k = 2, w = 0.01), 100L)
})

test_that("a loss constant near the largest double does not overflow on the way (k = 2)", {
  # A = w = 1e308: the criterion (A / w) S^2 is the sample variance itself.
  # The deviations of x from its mean 1.24 are 0.06, -1.04, 1.66, -0.54, and
  # -0.14, with squares summing to 4.152, so S^2 = 1.038 and the criterion
  # is 2; the estimated risk is A S^2 / n = 1e308 * 1.038 / 5 = 2.076e307.
  # (The general expression for B formed 2 A first, which is Inf, and the
  # function stopped with "implied sample size (Inf)".)
  x <- c(1.3, 0.2, 2.9, 0.7, 1.1)
  out <- seq_br_norm_mean(x, A = 1e308, k = 2, w = 1e308)
  expect_identical(attr(out, "criterion_n"), 2L)
  expect_true(attr(out, "is_satisfied"))
  expect_equal(out$value[out$term == "risk"], 2.076e307)
  # Pilot: sqrt(A / w) = sqrt(10) = 3.16, so 4.
  expect_identical(seq_br_norm_mean_pilot(A = 1e308, k = 2, w = 1e307,
                                          m_0 = 2), 4L)
})

test_that("seq_br_norm_mean returns the documented rows and attributes", {
  set.seed(113)
  x <- rnorm(50, 3, 2)
  out <- seq_br_norm_mean(x, A = 2, k = 1, w = 0.5)
  expect_s3_class(out, c("smsd_check", "data.frame"), exact = TRUE)
  expect_identical(out$term, c("mean", "sd", "risk", "optimal_n", "current_n"))
  expect_equal(out$value[1:2], c(mean(x), stats::sd(x)))
  expect_equal(out$value[5], 50)
  expect_null(attr(out, "conf_level"))            # risk based
  expect_null(attr(out, "width"))
  expect_identical(attr(out, "loss_A"), 2)
  expect_identical(attr(out, "k"), 1)
  expect_identical(attr(out, "risk_bound"), 0.5)
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_identical(attr(out, "notes"), character(0))
  expect_identical(attr(out, "method"), "bounded risk estimation (normal mean)")
  expect_equal(seq_br_norm_mean(c(x, NA), A = 2, k = 1, w = 0.5)$value, out$value)
  # The data argument is x, as for every one sample procedure.
  expect_identical(names(formals(seq_br_norm_mean)),
                   c("x", "A", "k", "w", "na.rm"))
  expect_identical(seq_br_norm_mean(x = x, A = 2, k = 1, w = 0.5), out)
  # The name used before 1.0.0 is R's "unused argument" error (its wording
  # is not matched because R translates it).
  expect_error(seq_br_norm_mean(data = x, A = 2, k = 1, w = 0.5))
  expect_identical(names(formals(seq_br_norm_mean_pilot)),
                   c("A", "k", "w", "m_0", "gamma"))
})

# ---- pilot ------------------------------------------------------------------

test_that("the pilot is max(m_0, ceiling(K^(1 / (1 + gamma)))) with K = (B / w)^(2 / k)", {
  # k = 1, A = 2, w = 0.05: B = 2 sqrt(2 / pi) = 1.595769, K = 31.91538^2
  # = 1018.59, sqrt(K) = 31.92 -> 32.
  expect_identical(seq_br_norm_mean_pilot(A = 2, k = 1, w = 0.05), 32L)
  # k = 2, A = 2, w = 0.0049: B = 2, K = 408.16, sqrt(K) = 20.20 -> 21.
  expect_identical(seq_br_norm_mean_pilot(A = 2, k = 2, w = 0.0049), 21L)
  # gamma = 0.49: 408.16^(1 / 1.49) = 56.52 -> 57.
  expect_identical(seq_br_norm_mean_pilot(A = 2, k = 2, w = 0.0049,
                                          gamma = 0.49), 57L)
  # A loose bound: K = 4 (A = 2, k = 2, w = 0.5), so the minimum m_0
  # governs. It was 4 under the old rule, with risk 1.66 times the bound.
  expect_identical(seq_br_norm_mean_pilot(A = 2, k = 2, w = 0.5), 10L)
  expect_identical(seq_br_norm_mean_pilot(A = 2, k = 2, w = 0.5, m_0 = 3), 3L)
  # The old rule returned 1 here, fewer than the check accepts.
  expect_identical(seq_br_norm_mean_pilot(A = 0.1, k = 1, w = 5), 10L)
})

test_that("the pilot is of lower order than the required sample size", {
  # n* = K sigma^2. With K = 408 the old pilot (K itself) was sixteen times
  # the 25.5 needed at sigma = 0.25; the pilot must stay below that.
  expect_lt(seq_br_norm_mean_pilot(A = 2, k = 2, w = 0.0049), 408.16 * 0.25^2)
  w <- c(0.5, 0.05, 0.005)
  frac <- vapply(w, function(wi) seq_br_norm_mean_pilot(A = 2, k = 2, w = wi) /
                   (2 / wi), numeric(1))
  expect_true(all(diff(frac) < 0))
})

# ---- degenerate and invalid input -------------------------------------------

test_that("constant data are 'not met' with a note, never an optimal n of 0", {
  for (k in 1:2) {
    out <- seq_br_norm_mean(rep(1, 2), A = 2, k = k, w = 0.5)
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_identical(attr(out, "variance_unusable"), TRUE)
    expect_true(all(is.na(out$value[out$term %in% c("risk", "optimal_n")])))
    expect_equal(out$value[out$term == "sd"], 0)
    # Exactly equal observations: the package's standard note, alone.
    expect_length(attr(out, "notes"), 1L)
    expect_match(attr(out, "notes"),
                 "not yet usable .*zero, negative, or not finite")
  }
})

test_that("data that are constant up to rounding error are 'not met' as well", {
  # Ten gain scores of 0.3, each a posttest minus a pretest. On paper they
  # are equal; in floating point they differ by about 4e-15 and have a
  # "variance" near 2e-30, which satisfied the rule with a criterion of 1.
  pre <- c(101, 57, 12, 80, 33, 64, 90, 45, 71, 28)
  gain <- (pre + 0.3) - pre
  expect_gt(stats::var(gain), 0)
  expect_lt(max(gain) - min(gain), 1e-14)
  out <- seq_br_norm_mean(gain, A = 1, k = 2, w = 0.001)
  expect_identical(attr(out, "is_satisfied"), FALSE)
  expect_identical(attr(out, "variance_unusable"), TRUE)
  expect_true(all(is.na(out$value[out$term %in% c("risk", "optimal_n")])))
  # Two notes: the standard one, and one that says why observations that are
  # not all equal were treated as constant.
  expect_length(attr(out, "notes"), 2L)
  expect_match(attr(out, "notes")[1], "not yet usable")
  expect_match(attr(out, "notes")[2], "agree to ten significant digits")
})

test_that("the sd row is a number when the variance of constant data overflows", {
  # Every row holds a finite number or a documented NA, never Inf or NaN.
  rows_ok <- function(out) {
    all(is.finite(out$value) | (is.na(out$value) & !is.nan(out$value)))
  }
  # Four observations near 2^1000 (about 1e301) that differ in the twelfth
  # significant digit: 2^1000 (1 + j 2^-40) for j = 0, 1, 3, 1. They are
  # constant up to rounding (range 3 * 2^-40 = 2.7e-12 of their size), and
  # the squares in their variance overflow, so stats::var() is Inf and its
  # square root was reported in the sd row. The standard deviation itself is
  # an ordinary number: in units of 2^960 the observations are j, with mean
  # 1.25 and squared deviations 1.5625 + 0.0625 + 3.0625 + 0.0625 = 4.75, so
  # it is 2^960 sqrt(4.75 / 3) = 2^960 sqrt(19 / 12), about 1.23e289. Every
  # step but the last division and the square root is exact in binary, so
  # the comparison is to twelve digits.
  j <- c(0, 1, 3, 1)
  x <- 2^1000 * (1 + j * 2^-40)
  expect_identical(stats::var(x), Inf)
  for (k in 1:2) {
    out <- seq_br_norm_mean(x, A = 1, k = k, w = 1)
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_identical(attr(out, "variance_unusable"), TRUE)
    expect_length(attr(out, "notes"), 2L)
    expect_equal(out$value[out$term == "sd"], 2^960 * sqrt(19 / 12),
                 tolerance = 1e-12)
    expect_equal(out$value[out$term == "mean"], 2^1000 * (1 + 1.25 * 2^-40),
                 tolerance = 1e-14)
    expect_true(all(is.na(out$value[out$term %in% c("risk", "optimal_n")])))
    expect_true(rows_ok(out))
    expect_output(print(out), "Stopping rule met: FALSE")
  }
  # The same pattern in decimal near 1e300 and 1e170: 1e287 sqrt(19 / 12)
  # and 1e157 sqrt(19 / 12). The stored doubles carry a relative difference
  # of 1e-13 only to about three digits, hence the tolerance.
  for (size in c(1e300, 1e170)) {
    out <- seq_br_norm_mean(size * (1 + j * 1e-13), A = 1, k = 2, w = 1)
    expect_identical(attr(out, "variance_unusable"), TRUE)
    expect_equal(out$value[out$term == "sd"], size * 1e-13 * sqrt(19 / 12),
                 tolerance = 1e-3)
    expect_true(rows_ok(out))
  }
  # At the other end the squares underflow: the variance of
  # 2^-1000 (1 + j 2^-40) is zero in floating point, and the standard
  # deviation is 2^-1040 sqrt(19 / 12), about 1.07e-313, not zero. (The
  # comparison is of the ratio with 1, because expect_equal() compares
  # numbers this small by their absolute difference, which a zero would
  # pass. A number of this size is stored to about ten digits.)
  y <- 2^-1000 * (1 + j * 2^-40)
  expect_identical(stats::var(y), 0)
  out <- seq_br_norm_mean(y, A = 1, k = 2, w = 1)
  expect_identical(attr(out, "variance_unusable"), TRUE)
  expect_gt(out$value[out$term == "sd"], 0)
  expect_equal(out$value[out$term == "sd"] / (2^-1040 * sqrt(19 / 12)), 1,
               tolerance = 1e-9)
  expect_true(rows_ok(out))
  # Exactly equal observations have a standard deviation of exactly zero at
  # any magnitude: near the largest double (where the sum of the data
  # overflows, so that stats::var() need not be zero) and at the smallest.
  for (v in c(1.7e308, -1.7e308, 2^1023, 5, 0, 1e-310, 2^-1074)) {
    for (n in 2:4) {
      out <- seq_br_norm_mean(rep(v, n), A = 1, k = 2, w = 1)
      expect_identical(out$value[out$term == "sd"], 0)
      expect_length(attr(out, "notes"), 1L)
      expect_true(rows_ok(out))
    }
  }
  # Two copies of the largest double: their mean is that double (half of
  # each is exact), and their standard deviation is zero, not NaN (log2()
  # of the largest double rounds up to 1024, and 2^1024 is not a double).
  out <- seq_br_norm_mean(rep(.Machine$double.xmax, 2), A = 1, k = 2, w = 1)
  expect_identical(out$value[1:2], c(.Machine$double.xmax, 0))
  expect_true(rows_ok(out))
  # Three copies of the largest double itself: where mean() cannot form
  # their mean (it returns Inf on a platform without extended precision),
  # the function stops and says to change the units; it never returns Inf
  # in a row.
  top <- tryCatch(seq_br_norm_mean(rep(.Machine$double.xmax, 3), A = 1, k = 2,
                                   w = 1), error = conditionMessage)
  if (is.character(top)) {
    expect_match(top, "'mean' row of the result is too large to represent")
  } else {
    expect_true(rows_ok(top))
  }
  # Ordinary data are untouched: the row is the square root of the sample
  # variance. For c(2, 4, 4, 6, 9) the variance is exactly 7.
  expect_identical(seq_br_norm_mean(c(2, 4, 4, 6, 9), A = 5, k = 2,
                                    w = 7)$value[2], sqrt(7))
})

test_that("the family's check object refuses a row that is Inf or NaN", {
  # .seq_check() is new_smsd_check(), which guards the value column, behind
  # a guard on the full width. NA is allowed (rows documented as NA); Inf,
  # -Inf, and NaN are not, and the message names the first such row.
  make <- function(value, width = NULL) {
    .seq_check(term = c("mean", "sd", "risk"), value = value, method = "m",
               is_satisfied = FALSE, current_n = 3, criterion_n = NA_real_,
               width = width)
  }
  ok <- make(c(1, 0, NA))
  expect_s3_class(ok, "smsd_check")
  expect_identical(ok$value, c(1, 0, NA_real_))
  expect_null(attr(ok, "width"))
  expect_identical(attr(make(c(1, 0, NA), width = 4), "width"), 4)
  expect_error(make(c(1, Inf, NA)),
               "^The 'sd' row of the result is too large to represent \\(Inf\\)\\.")
  expect_error(make(c(-Inf, 2, Inf)),
               "^The 'mean' row of the result is too large to represent \\(-Inf\\)\\.")
  expect_error(make(c(1, 2, NaN)),
               "^The 'risk' row of the result is too large to represent \\(NaN\\)\\.")
  expect_error(make(c(1, 2, 3), width = Inf),
               "^The full width of the interval, twice 'd', is too large to represent")
})

test_that("is_satisfied is the authority when the risk and w agree to ten digits", {
  # c(2, 4, 4, 6, 9) has variance exactly 7, so with A = 5, k = 2, and w = 7
  # the criterion A S^2 / w is exactly 5 and the estimated risk A S^2 / n at
  # n = 5 is exactly 7. Multiplying the data by sqrt(1 + e) multiplies the
  # variance, the criterion, and the risk by 1 + e. The package's rounding
  # up ignores an excess of one part in 1e10 (so that an exact 5 computed as
  # 5.000000000000001 is not 6). Hence:
  #   e = 5e-11: criterion 5.00000000025, within the allowance, so optimal_n
  #     is 5 and the rule is met, although the risk row, 7.00000000035, is
  #     above w in the eleventh significant digit. This is the case the help
  #     page describes.
  #   e = 1e-9: criterion 5.000000005, beyond the allowance: optimal_n 6,
  #     not met, risk above w.
  #   e = -1e-9: criterion 4.999999995: optimal_n 5, met, risk below w.
  base <- c(2, 4, 4, 6, 9)
  at <- function(e) seq_br_norm_mean(base * sqrt(1 + e), A = 5, k = 2, w = 7)
  near <- at(5e-11)
  expect_identical(attr(near, "criterion_n"), 5L)
  expect_true(attr(near, "is_satisfied"))
  expect_gt(near$value[near$term == "risk"], 7)
  expect_equal(near$value[near$term == "risk"], 7 * (1 + 5e-11),
               tolerance = 1e-13)
  above <- at(1e-9)
  expect_identical(attr(above, "criterion_n"), 6L)
  expect_false(attr(above, "is_satisfied"))
  expect_gt(above$value[above$term == "risk"], 7)
  below <- at(-1e-9)
  expect_identical(attr(below, "criterion_n"), 5L)
  expect_true(attr(below, "is_satisfied"))
  expect_lt(below$value[below$term == "risk"], 7)
})

test_that("a variance or a risk that cannot be represented is an error", {
  x <- c(1.3, 0.2, 2.9, 0.7, 1.1)
  # The squares of data near 1e160 overflow; those of data near 1e-200
  # underflow to zero although the observations differ.
  expect_error(seq_br_norm_mean(x * 1e160, A = 1, k = 1, w = 1e158),
               "sample variance of 'x' is too large to represent")
  expect_error(seq_br_norm_mean(x * 1e-200, A = 1, k = 1, w = 1),
               "sample variance of 'x' is too small to represent")
  # Squared error with A = w = 8e307: the criterion is the sample variance
  # of x * 10, about 104, a perfectly good sample size, but the risk
  # A s^2 / n = 8e307 * 104 / 5 overflows. It used to be reported as Inf
  # beside a definite stopping decision. The message shows the arguments.
  expect_lt(stats::var(x * 10), 200)
  expect_error(seq_br_norm_mean(x * 10, A = 8e307, k = 2, w = 8e307),
               "estimated risk is too large to represent \\(Inf\\) with A = 8e\\+307, k = 2")
})

test_that("seq_br_norm_mean and its pilot check every argument", {
  set.seed(113)
  x <- rnorm(30)
  expect_error(seq_br_norm_mean(x, k = 1, w = 0.5), "'A' must be supplied")
  expect_error(seq_br_norm_mean(x, A = 2, w = 0.5), "'k' must be supplied")
  expect_error(seq_br_norm_mean(x, A = 2, k = 1), "'w' must be supplied")
  expect_error(seq_br_norm_mean(x, A = NA_real_, k = 1, w = 0.5), "'A' must be")
  expect_error(seq_br_norm_mean(x, A = -2, k = 1, w = 0.5), "'A' must be")
  expect_error(seq_br_norm_mean(x, A = 2, k = 3, w = 0.5), "'k' must be 1")
  expect_error(seq_br_norm_mean(x, A = 2, k = NA, w = 0.5), "'k' must be 1")
  expect_error(seq_br_norm_mean(x, A = 2, k = 1.5, w = 0.5), "'k' must be 1")
  expect_error(seq_br_norm_mean(x, A = 2, k = 1, w = 0), "'w' must be")
  expect_error(seq_br_norm_mean(x, A = 2, k = 1, w = Inf), "'w' must be")
  expect_error(seq_br_norm_mean(A = 2, k = 1, w = 0.5), "'x' must be supplied")
  expect_error(seq_br_norm_mean(1, A = 2, k = 1, w = 0.5),
               "'x' needs at least 2 observations; it has 1\\. Collect more data")
  expect_error(seq_br_norm_mean("a", A = 2, k = 1, w = 0.5),
               "'x' must be a numeric vector")
  expect_error(seq_br_norm_mean(c(x, NA), A = 2, k = 1, w = 0.5, na.rm = FALSE),
               "missing values")
  # The advice that closes the too large message names what to check.
  expect_error(seq_br_norm_mean(x, A = 2, k = 2, w = 1e-12),
               "implied sample size .* is too large to represent\\. Check 'A', the risk bound 'w', and the units of 'x'\\.$")

  expect_error(seq_br_norm_mean_pilot(k = 1, w = 0.5), "'A' must be supplied")
  expect_error(seq_br_norm_mean_pilot(A = 2, k = 0, w = 0.5), "'k' must be 1")
  expect_error(seq_br_norm_mean_pilot(A = 2, k = 1, w = NA_real_), "'w' must be")
  expect_error(seq_br_norm_mean_pilot(A = 2, k = 1, w = 0.5, m_0 = 1), "'m_0'")
  expect_error(seq_br_norm_mean_pilot(A = 2, k = 1, w = 0.5, m_0 = 2.5), "'m_0'")
  expect_error(seq_br_norm_mean_pilot(A = 2, k = 1, w = 0.5, gamma = 0), "'gamma'")
  expect_error(seq_br_norm_mean_pilot(A = 2, k = 2, w = 1e-25),
               "implied pilot sample size .* is too large to represent\\. Check 'A' and the risk bound 'w'\\.$")
  expect_type(seq_br_norm_mean_pilot(A = 2, k = 1, w = 0.5), "integer")
})

# ---- operating characteristics ----------------------------------------------

test_that("the risk is at the bound, within Monte Carlo error, at a moderate target", {
  # N(10, 3^2), A = 2, bound chosen so that n* = 400. The complete procedure
  # is run 2,000 times from the package's pilot.
  for (k in 1:2) {
    set.seed(113)
    A <- 2
    sigma <- 3
    n_star <- 400
    # w solving A E|Z|^k (sigma^2 / n*)^(k / 2) = w, with E|Z| = sqrt(2 / pi)
    # and E Z^2 = 1.
    w <- if (k == 1) A * sqrt(2 / pi) * sigma / sqrt(n_star) else A * sigma^2 / n_star
    expect_equal(br_risk(A, k, sigma^2, n_star), w, tolerance = 1e-7)
    K <- n_star / sigma^2
    m <- seq_br_norm_mean_pilot(A, k, w)
    sim <- seqc_simulate(
      2000, m, rgen = function(n) rnorm(n, 10, sigma),
      ok = function(x, n) n >= K * seqc_run_mv(x - 10, n)$var,
      n_guess = n_star,
      measure = function(x, N) A * abs(mean(x[seq_len(N)]) - 10)^k,
      check = function(x) seq_br_norm_mean(x, A = A, k = k, w = w))
    expect_true(sim$agree)
    # The loss has coefficient of variation sqrt(2) for k = 2 and 0.756 for
    # k = 1, so the SE of risk / w from 2,000 runs is 0.032 and 0.017. The
    # risk must be within 4 SE of the bound on either side.
    se <- if (k == 2) sqrt(2 / 2000) else sqrt((pi / 2 - 1) / 2000)
    expect_lt(sim$mean_y / w, 1 + 4 * se)
    expect_gt(sim$mean_y / w, 1 - 4 * se)
    # sd(N) is about sqrt(2 n*) = 28, so the SE of the mean of N is 0.6.
    expect_lt(abs(sim$mean_n / n_star - 1), 0.02)
    if (k == 2) {
      # For squared error and normal data the final sample size depends on
      # the data only through the sample variances, which are independent of
      # the sample means, so the risk is exactly A sigma^2 E(1 / N) and
      # risk / w = n* E(1 / N). That estimate has an SE near 0.002. It must
      # sit just above 1 (the bound is approached from above), not below it
      # and not far above.
      expect_gt(n_star * mean(1 / sim$N), 1)
      expect_lt(n_star * mean(1 / sim$N), 1.03)
    }
  }
})

test_that("seq_br_norm_mean does not presume a scale: sd 0.25, 1, and 3", {
  # Squared error, A = 2, w = 0.005: n* = 400 sigma^2, that is 25, 400, and
  # 3600. The pilot is 20 whatever the scale (the old pilot was 400: sixteen
  # times the 25 needed at sd 0.25). 1,000 runs at each scale. For squared
  # error and normal data risk / w = n* E(1 / N) exactly (see the test above),
  # which is estimated far more precisely than the mean loss. An independent
  # 40,000 run simulation gave E[N] / n* = 1.014, 0.998, 1.000 and
  # risk / w = 1.03, 1.00, 1.00.
  A <- 2
  w <- 0.005
  m <- seq_br_norm_mean_pilot(A, 2, w)
  expect_identical(m, 20L)
  expected_ratio <- c(1.014, 0.998, 1.000)
  # sd(N) is about 5, 28, and 85, so 4 SE of the mean of N over n* is 0.027,
  # 0.009, and 0.003.
  tol_ratio <- c(0.03, 0.015, 0.01)
  set.seed(113)
  for (i in 1:3) {
    sigma <- c(0.25, 1, 3)[i]
    n_star <- A * sigma^2 / w
    sim <- seqc_simulate(
      1000, m, rgen = function(n) rnorm(n, 0, sigma),
      ok = function(x, n) n >= (A / w) * seqc_run_mv(x, n)$var,
      n_guess = n_star,
      measure = function(x, N) A * mean(x[seq_len(N)])^2,
      check = function(x) seq_br_norm_mean(x, A = A, k = 2, w = w),
      n_confirm = 20L, n_full = 1L)
    expect_true(sim$agree)
    expect_lt(abs(sim$mean_n / n_star - expected_ratio[i]), tol_ratio[i])
    # At or just above the bound at every scale, never far above it: the SE
    # of n* mean(1 / N) is at most 0.007.
    expect_gt(n_star * mean(1 / sim$N), 0.99)
    expect_lt(n_star * mean(1 / sim$N), 1.08)
  }
})
