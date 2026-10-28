# Tests for the two stage bounded risk procedures: ts_br_norm_mean(),
# ts_br_exp_mean(), and their pilot functions.
#
# Every expected value is derived here, independently of the package code: a
# closed form, arithmetic shown in a comment, a numerical integration of the
# defining expectation, or an exact operating characteristic computed from
# the distribution of the first stage statistic that fixes N (helpers
# ts_n_dist() and ts_br_exp_exact() in helper-ts.R).

# The stopping constant by numerical integration of its definition,
# E[(S^2 / sigma^2)^(-k/2)] with (m - 1) S^2 / sigma^2 chi square on m - 1
# degrees of freedom (no gamma functions).
b_by_integration <- function(m, k) {
  stats::integrate(
    function(v) (v / (m - 1))^(-k / 2) * stats::dchisq(v, m - 1),
    lower = stats::qchisq(1e-13, m - 1),
    upper = stats::qchisq(1e-13, m - 1, lower.tail = FALSE),
    rel.tol = 1e-10, subdivisions = 2000L)$value
}

# ---- ts_br_norm_mean: the stopping constant ------------------------------

test_that("the stopping constant b matches its closed forms and its definition", {
  # k = 2: b = (m - 1) / (m - 3), to the last bit where that is a whole
  # number (a value one unit in the last place high would push a whole
  # number criterion to the next integer).
  expect_identical(.ts_br_norm_b(4, 2), 3)
  expect_identical(.ts_br_norm_b(5, 2), 2)
  expect_identical(.ts_br_norm_b(7, 2), 1.5)
  expect_equal(.ts_br_norm_b(10, 2), 9 / 7)
  expect_equal(.ts_br_norm_b(500, 2), 499 / 497)
  expect_equal(.ts_br_norm_b(1019, 2), 1018 / 1016)
  expect_equal(.ts_br_norm_b(1e6, 2), (1e6 - 1) / (1e6 - 3))
  # k = 1: b(3) = sqrt(pi), b(5) = sqrt(pi / 2).
  expect_equal(.ts_br_norm_b(3, 1), sqrt(pi))
  expect_equal(.ts_br_norm_b(5, 1), sqrt(pi / 2))
  # k = 1, any m: b(m) b(m + 1) = sqrt(m (m - 1)) / (m - 2), an identity
  # free of gamma functions. Checked where gamma() itself overflows.
  for (m in c(3, 10, 344, 345, 500, 1019)) {
    expect_equal(.ts_br_norm_b(m, 1) * .ts_br_norm_b(m + 1, 1),
                 sqrt(m * (m - 1)) / (m - 2))
  }
  # Against the defining expectation, by numerical integration.
  for (k in 1:2) for (m in c(k + 2, 6, 30, 500)) {
    expect_equal(.ts_br_norm_b(m, k), b_by_integration(m, k),
                 tolerance = 1e-6)
  }
})

test_that("the risk constant B is A E|Z|^k, and exactly A for squared error", {
  # E|Z|^k for a standard normal Z by numerical integration of its
  # definition (no gamma function): sqrt(2 / pi) for k = 1, 1 for k = 2.
  abs_moment <- function(k) {
    2 * stats::integrate(function(z) z^k * stats::dnorm(z), 0, Inf,
                         rel.tol = 1e-12)$value
  }
  expect_equal(abs_moment(1), sqrt(2 / pi), tolerance = 1e-9)
  expect_equal(abs_moment(2), 1, tolerance = 1e-9)
  expect_equal(.br_constant(3, 1), 3 * abs_moment(1), tolerance = 1e-9)
  expect_equal(.br_constant(3, 1), 3 * 0.7978846, tolerance = 1e-7)
  # k = 2: B is A itself, with no rounding error at all. (Before 1.0.0 the
  # general gamma function expression gave A * (1 + 2.2e-16).)
  for (A in c(1, 2, 3, 5, 10, 100, 0.1, 1e-3, 7.25)) {
    expect_identical(.br_constant(A, 2), A)
  }
})

# ---- ts_br_norm_mean: criterion and guards -------------------------------

test_that("ts_br_norm_mean criterion matches the closed form (small first stage)", {
  x <- c(2, 4, 4, 6, 9)           # m = 5, mean 5, variance 28 / 4 = 7
  # k = 2: B = A, b = (m - 1) / (m - 3) = 2, so N = ceiling((2 A / w) * 7).
  # A = 2, w = 0.5: (4 / 0.5) * 7 = 56 exactly, so N = 56, not 57.
  # A = 1, w = 1:   (2 / 1) * 7   = 14 exactly, so N = 14, not 15.
  # Every step of these is exact in floating point, so the whole number
  # must come back as itself.
  expect_identical(attr(ts_br_norm_mean(x, A = 2, k = 2, w = 0.5),
                        "criterion_n"), 56L)
  expect_identical(attr(ts_br_norm_mean(x, A = 1, k = 2, w = 1),
                        "criterion_n"), 14L)
  expect_identical(attr(ts_br_norm_mean(x, A = 3, k = 2, w = 2),
                        "criterion_n"), 21L)       # (6 / 2) * 7
  expect_identical(attr(ts_br_norm_mean(x, A = 5, k = 2, w = 1),
                        "criterion_n"), 70L)       # (10 / 1) * 7
  expect_identical(attr(ts_br_norm_mean(x, A = 10, k = 2, w = 2),
                        "criterion_n"), 70L)       # (20 / 2) * 7
  # A = 2, w = 0.6:
  # N = ceiling((b B / w) S^2) = ceiling((4 / 0.6) * 7) = ceiling(46.67) = 47.
  # (Without b the size would be ceiling(23.33) = 24.)
  out <- ts_br_norm_mean(x, A = 2, k = 2, w = 0.6)
  v <- stats::setNames(out$value, out$term)
  expect_equal(v[["optimal_n"]], 47)
  expect_equal(v[["mean"]], 5)
  expect_equal(v[["current_n"]], 5)
  # Estimated risk at the current n: B S^2 / n = 2 * 7 / 5.
  expect_equal(v[["risk"]], 2.8)
  expect_false(attr(out, "is_satisfied"))
  expect_null(attr(out, "conf_level"))
  expect_null(attr(out, "width"))
  expect_equal(attr(out, "loss_A"), 2)
  expect_equal(attr(out, "k"), 2)
  expect_equal(attr(out, "risk_bound"), 0.6)
  expect_identical(out$term, c("mean", "risk", "optimal_n", "current_n"))
  expect_identical(attr(out, "method"),
                   "two stage bounded risk estimation (normal mean)")
  # k = 1, A = 2, w = 0.6: B = 2 sqrt(2 / pi), b(5) = sqrt(pi / 2), so
  # b B = 2 and N = ceiling((2 / 0.6)^2 * 7) = ceiling(77.78) = 78.
  out1 <- ts_br_norm_mean(x, A = 2, k = 1, w = 0.6)
  expect_equal(attr(out1, "criterion_n"), 78)
  # Risk row for k = 1: B sqrt(S^2 / n) = 2 sqrt(2 / pi) sqrt(7 / 5).
  expect_equal(out1$value[out1$term == "risk"],
               2 * sqrt(2 / pi) * sqrt(7 / 5))
  # k = 1, A = 2, w = 0.5: b B = 2 as above, so the criterion is
  # (2 / 0.5)^2 * 7 = 112 exactly. In floating point b B is 2 only up to
  # rounding error (sqrt(pi / 2) * 2 * sqrt(2 / pi)), the product lands a
  # few units in the last place above 112, and a plain ceiling gave 113. A
  # whole number criterion comes back as itself.
  expect_identical(attr(ts_br_norm_mean(x, A = 2, k = 1, w = 0.5),
                        "criterion_n"), 112L)
  # k = 1, A = 1, w = 0.5: (1 / 0.5)^2 * 7 = 28 exactly.
  expect_identical(attr(ts_br_norm_mean(x, A = 1, k = 1, w = 0.5),
                        "criterion_n"), 28L)
  # N is fixed by the first stage; the rule is met exactly at N = 47.
  short <- ts_br_norm_mean(x, rep(c(4, 6), 21)[-1], A = 2, k = 2, w = 0.6)
  met <- ts_br_norm_mean(x, rep(c(4, 6), 21), A = 2, k = 2, w = 0.6)
  expect_equal(attr(short, "current_n"), 46)
  expect_equal(attr(short, "criterion_n"), 47)
  expect_false(attr(short, "is_satisfied"))
  expect_equal(attr(met, "current_n"), 47)
  expect_equal(attr(met, "criterion_n"), 47)
  expect_true(attr(met, "is_satisfied"))
  # The estimate is the mean of both stages together.
  expect_equal(met$value[met$term == "mean"], mean(c(x, rep(c(4, 6), 21))))
})

test_that("ts_br_norm_mean reports the rows for both stages together", {
  x <- c(2, 4, 4, 6, 9)
  second <- c(1, 8, 5)
  # All eight observations: sum 39, sum of squares 243, so the mean is
  # 39 / 8 = 4.875 and the variance (243 - 39^2 / 8) / 7 = 52.875 / 7.
  out <- ts_br_norm_mean(x, second, A = 2, k = 2, w = 0.6)
  v <- stats::setNames(out$value, out$term)
  expect_equal(v[["mean"]], 4.875)
  # Risk at the current n from all observations: B S^2 / n with B = A = 2,
  # (2 * 52.875 / 7) / 8 = 1.888393. (From the first stage alone it would be
  # 2 * 7 / 8 = 1.75, and at the first stage size 2 * 7 / 5 = 2.8.)
  expect_equal(v[["risk"]], 2 * (52.875 / 7) / 8)
  expect_equal(v[["risk"]], 1.888393, tolerance = 1e-6)
  expect_equal(v[["optimal_n"]], 47)               # the first stage alone
  expect_equal(v[["current_n"]], 8)
  expect_identical(attr(out, "current_n"), 8L)
  expect_identical(attr(out, "criterion_n"), 47L)
  # k = 1: B sqrt(S^2 / n) = 2 sqrt(2 / pi) sqrt(52.875 / 56) = 1.550605.
  out1 <- ts_br_norm_mean(x, second, A = 2, k = 1, w = 0.6)
  expect_equal(out1$value[out1$term == "risk"],
               2 * sqrt(2 / pi) * sqrt(52.875 / 56))
  expect_equal(out1$value[out1$term == "current_n"], 8)
})

test_that("ts_br_norm_mean stops when the first stage is below k + 2", {
  withr::local_seed(113)
  # Before 1.0.0: m = 2, k = 2 gave a negative constant and 'satisfied' at
  # n = 2; m = k + 1 gave NaN.
  expect_error(ts_br_norm_mean(rnorm(2), A = 1, k = 2, w = 1),
               "at least k \\+ 2 = 4")
  expect_error(ts_br_norm_mean(rnorm(3), A = 1, k = 2, w = 1),
               "at least k \\+ 2 = 4")                      # m = k + 1
  expect_error(ts_br_norm_mean(rnorm(2), A = 1, k = 1, w = 1),
               "at least k \\+ 2 = 3")                      # m = k + 1
  # The count is of usable observations (missing values removed), in the
  # wording of the shared data check, followed by the reason and the remedy.
  expect_error(ts_br_norm_mean(c(1, NA, 2, NA, 3), A = 1, k = 2, w = 1),
               "'pilot_data' needs at least 4 observations; it has 3\\.")
  expect_error(ts_br_norm_mean(numeric(0), A = 1, k = 1, w = 1),
               "'pilot_data' needs at least 3 observations; it has 0\\.")
  expect_error(ts_br_norm_mean(7, A = 1, k = 1, w = 1),
               "'pilot_data' needs at least 3 observations; it has 1\\.")
  expect_error(ts_br_norm_mean(7, A = 1, k = 1, w = 1),
               "standard deviation to the power -k\\) exists only for a first stage of at least k \\+ 2 = 3 observations")
  expect_error(ts_br_norm_mean(7, A = 1, k = 1, w = 1),
               "Enlarge the first stage; ts_br_norm_mean_pilot\\(\\) returns a suitable size")
  expect_false(grepl("check the stopping rule again",
                     tryCatch(ts_br_norm_mean(7, A = 1, k = 1, w = 1),
                              error = conditionMessage)))
  # The smallest permitted first stage works for each k.
  expect_s3_class(ts_br_norm_mean(rnorm(4), A = 1, k = 2, w = 1), "smsd_check")
  expect_s3_class(ts_br_norm_mean(rnorm(3), A = 1, k = 1, w = 1), "smsd_check")
})

test_that("ts_br_norm_mean is finite and correct for a large first stage", {
  withr::local_seed(113)
  x <- rnorm(500, mean = 20, sd = 1.5)
  # k = 2, A = 2, w = 0.004 (the configuration that returned NaN before
  # 1.0.0, because gamma() overflows beyond m = 344): B = 2,
  # b = 499 / 497, N = ceiling((499 / 497) * 500 * S^2).
  out <- ts_br_norm_mean(x, A = 2, k = 2, w = 0.004)
  expect_true(is.finite(attr(out, "criterion_n")))
  expect_equal(attr(out, "criterion_n"),
               max(500, ceiling((499 / 497) * (2 / 0.004) * stats::var(x))))
  expect_false(is.na(attr(out, "is_satisfied")))
  # k = 1, A = 2, w = 0.05, with b from numerical integration.
  out1 <- ts_br_norm_mean(x, A = 2, k = 1, w = 0.05)
  b <- b_by_integration(500, 1)
  expect_equal(attr(out1, "criterion_n"),
               max(500, ceiling((b * 2 * sqrt(2 / pi) / 0.05)^2 *
                                  stats::var(x))))
  # And well beyond: m = 5000.
  big <- ts_br_norm_mean(rnorm(5000), A = 2, k = 2, w = 0.004)
  expect_true(is.finite(attr(big, "criterion_n")))
})

test_that("ts_br_norm_mean checks every argument and degenerate data", {
  x <- c(2, 4, 4, 6, 9, 5, 7)
  expect_error(ts_br_norm_mean(A = 2, k = 1, w = 0.5), "'pilot_data' must be supplied")
  expect_error(ts_br_norm_mean(x, k = 1, w = 0.5), "'A' must be supplied")
  expect_error(ts_br_norm_mean(x, A = 2, w = 0.5), "'k' must be supplied")
  expect_error(ts_br_norm_mean(x, A = 2, k = 1), "'w' must be supplied")
  expect_error(ts_br_norm_mean(x, A = -2, k = 1, w = 0.5), "'A'")
  expect_error(ts_br_norm_mean(x, A = NA_real_, k = 1, w = 0.5), "'A'")
  expect_error(ts_br_norm_mean(x, A = 2, k = 3, w = 0.5), "'k' must be 1")
  expect_error(ts_br_norm_mean(x, A = 2, k = NA_real_, w = 0.5), "'k'")
  expect_error(ts_br_norm_mean(x, A = 2, k = 1, w = 0), "'w'")
  expect_error(ts_br_norm_mean(x, A = 2, k = 1, w = Inf), "'w'")
  expect_error(ts_br_norm_mean(letters, A = 2, k = 1, w = 0.5), "'pilot_data' must be a numeric vector")
  expect_error(ts_br_norm_mean(matrix(1:8, 4), A = 2, k = 1, w = 0.5), "'pilot_data'")
  expect_error(ts_br_norm_mean(x, "b", A = 2, k = 1, w = 0.5), "'additional_data'")
  expect_error(ts_br_norm_mean(c(x, NA), A = 2, k = 1, w = 0.5, na.rm = FALSE), "missing values")
  expect_error(ts_br_norm_mean(x, A = 2, k = 1, w = 0.5, na.rm = "yes"), "'na.rm'")
  expect_equal(ts_br_norm_mean(c(x, NA), A = 2, k = 1, w = 0.5)$value,
               ts_br_norm_mean(x, A = 2, k = 1, w = 0.5)$value)
  # Constant first stage: no second stage size can be determined.
  expect_error(ts_br_norm_mean(rep(5, 6), A = 2, k = 1, w = 0.5),
               "first stage variance is zero")
  # Constant apart from rounding error is the same case (see test-ts_fwd.R).
  expect_error(ts_br_norm_mean(c(0.1 + 0.2, 0.3, 0.3, 0.3, 0.3), A = 1, k = 2,
                               w = 1e-9),
               "are equal up to rounding error")
  # A criterion too large to represent is an error, never NA, and the
  # message says what to look at.
  expect_error(ts_br_norm_mean(x, A = 2, k = 2, w = 1e-12),
               "implied total sample size .* is too large to represent")
  expect_error(ts_br_norm_mean(x, A = 2, k = 2, w = 1e-12),
               "Check that the loss constant 'A' and the risk bound 'w'")
})

# ---- data on a scale whose squares cannot be represented -----------------

test_that("observations whose squares cannot be represented are an error that says to change the units", {
  msg <- function(expr) tryCatch(expr, error = conditionMessage)
  larger <- paste0(" is too large to represent. Express the data in larger ",
                   "units (divide them by a constant) and restate the loss ",
                   "constant 'A' for those units.")
  smaller <- paste0("). Express the data in smaller units (multiply them by ",
                    "a constant) and restate the loss constant 'A' for ",
                    "those units.")
  # x = (2, 4, 4, 6, 9, 5, 7, 3, 8, 6): m = 10, variance 44.4 / 9 = 4.9333.
  x <- c(2, 4, 4, 6, 9, 5, 7, 3, 8, 6)
  # Absolute error (k = 1), A = 1, w = 0.25: b = sqrt(4.5) Gamma(4) /
  # Gamma(4.5) = 2.12132 * 6 / 11.6317 = 1.09424, B = sqrt(2 / pi) = 0.79788,
  # so N = (1.09424 * 0.79788 / 0.25)^2 * 4.9333 = 60.2, that is 61, and the
  # risk row at the first stage is 0.79788 * sqrt(4.9333 / 10) = 0.56042.
  # With the data and w in other units (w is in the units of the data for
  # k = 1) the answer is the same wherever the first stage variance can be
  # represented, from a scale of 1e-154 (variance 4.9e-308; the smallest
  # normal double is 2.2e-308) to 1e153. At 1e-154 the square of b B / w
  # overflows although b B S / w does not, which used to stop the function.
  for (s in c(1e-154, 1e-150, 1e-100, 1, 1e100, 1e150, 1e153)) {
    out <- ts_br_norm_mean(x * s, A = 1, k = 1, w = 0.25 * s)
    expect_identical(attr(out, "criterion_n"), 61L)
    expect_equal(out$value[out$term == "risk"], 0.5604157 * s,
                 tolerance = 1e-6)
  }
  # Squared error (k = 2), A = 1, w = 0.5 in squared units:
  # N = (9 / 7) * 4.9333 / 0.5 = 12.7, that is 13.
  for (s in c(1e-150, 1, 1e150)) {
    expect_identical(attr(ts_br_norm_mean(x * s, A = 1, k = 2, w = 0.5 * s^2),
                          "criterion_n"), 13L)
  }
  # Beyond about 1e154 the squares overflow. The function used to stop with
  # "The implied total sample size (Inf) is too large to represent. Check
  # that the loss constant 'A' and the risk bound 'w' ...", although N is 61.
  expect_identical(msg(ts_br_norm_mean(x * 1e160, A = 1, k = 1, w = 0.25e160)),
                   paste0("The first stage variance", larger))
  expect_identical(msg(ts_br_norm_mean(x * 1e155, A = 1e-300, k = 2, w = 1)),
                   paste0("The first stage variance", larger))
  # Below about 1e-154 they underflow (the variance, 4.9e-320 here, is below
  # the smallest normal double and has lost its digits).
  expect_identical(msg(ts_br_norm_mean(x * 1e-160, A = 1, k = 1, w = 0.25e-160)),
                   paste0("The first stage variance is too small to ",
                          "represent (it underflows although the ",
                          "observations differ", smaller))
  expect_identical(msg(ts_br_norm_mean(x * 1e-200, A = 1, k = 1, w = 0.25e-200)),
                   paste0("The first stage variance is too small to ",
                          "represent (it underflows although the ",
                          "observations differ", smaller))
  # A second stage on such a scale: the rule was reported as met beside a
  # risk of Inf (for k = 1 too, where the risk itself, about 5e159, can be
  # represented but its computation squares the observations).
  for (k in 1:2) {
    expect_identical(msg(ts_br_norm_mean(c(1, 2, 3, 4, 5), c(1e160, 2e160),
                                         A = 1, k = k, w = 1)),
                     paste0("The variance of the observations supplied",
                            larger))
  }

  # The exponential mean. x = (1, 2, 3): m = 3, mean 2, B = 0.565 * 12 A =
  # 6.78 A, so N = 6.78 A * 4 s^2 / w for data in units of s. With
  # A = 1e100, s = 1e-150, and w = 1e-205 that is 27.12 * 1e5 = 2712000.
  expect_identical(attr(ts_br_exp_mean(c(1, 2, 3) * 1e-150, A = 1e100,
                                       w = 1e-205), "criterion_n"), 2712000L)
  # The same study with the data in units 1e20 times smaller (and w in the
  # square of those units) has the same N, but the square of the mean,
  # 4e-340, underflows to zero: the function returned N = 3, "met", and a
  # risk of 0.
  expect_identical(msg(ts_br_exp_mean(c(1, 2, 3) * 1e-170, A = 1e100,
                                      w = 1e-245)),
                   paste0("The square of the first stage mean is too small ",
                          "to represent (it underflows although the mean is ",
                          "positive", smaller))
  expect_identical(msg(ts_br_exp_mean(c(1, 2, 3) * 1e160, A = 1e-300, w = 1)),
                   paste0("The square of the first stage mean", larger))
  # A second stage on a scale whose squares overflow: was "met" beside a
  # risk of Inf.
  expect_identical(msg(ts_br_exp_mean(c(1, 2, 3, 4, 5), c(1e160, 2e160, 3e160),
                                      A = 1, w = 100)),
                   paste0("The square of the sample mean", larger))
  expect_identical(msg(ts_br_exp_mean(c(1, 2, 3, 4, 5) * 1e-150,
                                      rep(1e-170, 2e5), A = 1e100,
                                      w = 1e-205)),
                   paste0("The square of the sample mean is too small to ",
                          "represent (it underflows although the mean is ",
                          "positive", smaller))

  # A risk that overflows through the loss constant alone (the squares are
  # fine) is refused as a row that cannot be represented.
  # ts_br_norm_mean: variance 2.5e20, b = 2, N = 2 * 1e290 / 1e302 * 2.5e20
  # = 5e8; risk = 1e290 * 2.5e20 / 5 = 5e309.
  expect_error(ts_br_norm_mean(c(1, 2, 3, 4, 5) * 1e10, A = 1e290, k = 2,
                               w = 1e302),
               "The 'risk' row of the result cannot be represented (it is computed as Inf): the data or the targets are too large in magnitude. Express the data in larger units (divide them by a constant) and restate the loss constant 'A' for those units.",
               fixed = TRUE)
  # ts_br_exp_mean: N = 2.825e290 * 9 / 1e284 = 2.5e7 from the first stage;
  # the mean of all ten observations is about 5e11, so the risk is
  # 1e290 * 2.5e23 / 10.
  expect_error(ts_br_exp_mean(c(1, 2, 3, 4, 5), rep(1e12, 5), A = 1e290,
                              w = 1e284),
               "The 'risk' row of the result cannot be represented (it is computed as Inf)",
               fixed = TRUE)
})

test_that("no two stage check returns a row that is not finite (seeded hostile inputs)", {
  # Observations and targets on scales from 1e-300 to 1e307, the two stages
  # often on different scales. Whatever comes back as an object must have
  # finite rows, a definite stopping flag, and a total of at least the first
  # stage; everything else must be one of the package's own errors. Before
  # the rows were checked, about one call in a hundred of the bounded risk
  # functions in such a sweep returned a rule reported as met beside a risk
  # of Inf.
  withr::local_seed(113)
  scales <- 10^c(-300, -200, -160, -155, -100, -10, 0, 10, 100, 150, 154, 160,
                 200, 300, 307)
  pick <- function() scales[sample(length(scales), 1)]
  returned <- 0L; refused <- 0L; bad <- character(0)
  one <- function(label, m, expr) {
    out <- tryCatch(expr, error = function(e) e)
    if (inherits(out, "error")) {
      refused <<- refused + 1L
      # An error of the package's own (.smsd_stop() drops the call), not one
      # raised by R from inside a computation.
      if (!is.null(conditionCall(out))) bad <<- c(bad, paste(label, "R error"))
      return(invisible())
    }
    returned <<- returned + 1L
    ok <- all(is.finite(out$value)) &&
      isTRUE(attr(out, "is_satisfied") %in% c(TRUE, FALSE)) &&
      is.integer(attr(out, "criterion_n")) && attr(out, "criterion_n") >= m
    if (!ok) bad <<- c(bad, label)
  }
  for (i in 1:250) {
    m <- sample(4:8, 1)
    s1 <- pick(); s2 <- if (stats::runif(1) < 0.5) s1 else pick()
    first <- stats::rexp(m) * s1
    second <- if (stats::runif(1) < 0.4) NULL else stats::rexp(sample(1:6, 1)) * s2
    sgn <- sample(c(-1, 1), 1)
    d <- stats::runif(1, 0.1, 5) * pick()
    A <- stats::runif(1, 0.5, 2) * if (stats::runif(1) < 0.7) 1 else pick()
    w <- stats::runif(1, 0.1, 5) * pick()
    one("ts_fwd_norm_mean", m, ts_fwd_norm_mean(sgn * first, second, d = d))
    one("ts_fwd_dfr_mean", m, ts_fwd_dfr_mean(sgn * first, second, d = d))
    one("ts_fwd_nexp_mean", m, ts_fwd_nexp_mean(sgn * first, second, d = d))
    one("ts_br_norm_mean 1", m, ts_br_norm_mean(sgn * first, second, A = A, k = 1, w = w))
    one("ts_br_norm_mean 2", m, ts_br_norm_mean(sgn * first, second, A = A, k = 2, w = w))
    one("ts_br_exp_mean", m, ts_br_exp_mean(first, second, A = A, w = w))
  }
  expect_identical(bad, character(0))
  # The sweep exercised both outcomes.
  expect_gt(returned, 100L)
  expect_gt(refused, 100L)
})

# ---- ts_br_norm_mean: operating characteristics --------------------------

test_that("ts_br_norm_mean keeps the risk at or under the bound (seeded Monte Carlo)", {
  withr::local_seed(113)
  m <- 5; sigma <- 3; mu <- 10; A <- 2; k <- 2; w <- 0.5; reps <- 2000
  n_star <- (A / w) * sigma^2                      # B = A for k = 2: 36
  res <- replicate(reps, {
    first <- rnorm(m, mu, sigma)
    total <- attr(ts_br_norm_mean(first, A = A, k = k, w = w), "criterion_n")
    out <- ts_br_norm_mean(first, rnorm(total - m, mu, sigma), A = A, k = k,
                           w = w)
    stopifnot(isTRUE(attr(out, "is_satisfied")))
    c(n = total, loss = A * (out$value[out$term == "mean"] - mu)^2)
  })
  # Exact values: N = max(m, ceiling(b n* V / (m - 1))), V chi square on
  # m - 1 degrees of freedom, b = (m - 1) / (m - 3); risk / w = E[n* / N].
  dist <- ts_n_dist(m, m - 1, ((m - 1) / (m - 3)) * n_star / (m - 1))
  exact_risk <- sum(dist$p * n_star / dist$n)      # relative to w
  exact_n <- sum(dist$p * dist$n)
  expect_lte(exact_risk, 1)                        # the guarantee itself
  expect_equal(exact_risk, 0.9132, tolerance = 5e-4)
  expect_equal(exact_n / n_star, 2.014, tolerance = 1e-3)

  # Tolerances are four Monte Carlo standard errors.
  se_loss <- stats::sd(res["loss", ]) / sqrt(reps)
  expect_lte(mean(res["loss", ]), w + 4 * se_loss)
  expect_lt(abs(mean(res["loss", ]) - exact_risk * w), 4 * se_loss)
  expect_lt(abs(mean(res["n", ]) - exact_n),
            4 * stats::sd(res["n", ]) / sqrt(reps))
  # A sharper check with the same N: given N the risk is A sigma^2 / N.
  cond <- A * sigma^2 / res["n", ] / w
  se_cond <- stats::sd(cond) / sqrt(reps)
  expect_lt(abs(mean(cond) - exact_risk), 4 * se_cond)
  # At or under the bound, with the same four standard errors (the exact
  # value 0.913 is fewer than four of them below 1).
  expect_lte(mean(cond), 1 + 4 * se_cond)
})

test_that("ts_br_norm_mean with its own pilot holds the bound at scales 0.25, 1, and 3", {
  withr::local_seed(113)
  A <- 1; k <- 2; w <- 0.01; reps <- 1000
  m <- ts_br_norm_mean_pilot(A = A, k = k, w = w)
  expect_identical(m, 10L)                         # sqrt(100) = 10
  for (sigma in c(0.25, 1, 3)) {
    n_star <- (A / w) * sigma^2                    # 6.25, 100, 900
    total <- replicate(reps, attr(
      ts_br_norm_mean(rnorm(m, 0, sigma), A = A, k = k, w = w),
      "criterion_n"))
    dist <- ts_n_dist(m, m - 1, ((m - 1) / (m - 3)) * n_star / (m - 1))
    exact_risk <- sum(dist$p * n_star / dist$n)
    exact_n <- sum(dist$p * dist$n)
    expect_lte(exact_risk, 1)
    # Given N the risk relative to w is n* / N. Four standard errors.
    cond <- n_star / total
    expect_lt(abs(mean(cond) - exact_risk), 4 * stats::sd(cond) / sqrt(reps))
    expect_lt(abs(mean(total) - exact_n),
              4 * stats::sd(total) / sqrt(reps) + 1e-8)
    # The expected total stays within a small multiple of the optimum at
    # every scale; the pre-1.0.0 pilot (the optimum at sigma = 1, here 100)
    # could not give less than 100 / 6.25 = 16 times it at sigma = 0.25.
    expect_lt(exact_n / n_star, 2)
  }
})

test_that("ts_br_norm_mean_pilot follows the documented rule and minimum", {
  # k = 2: (B / w)^(2/k) = A / w; gamma = 1 gives ceiling(sqrt(A / w)).
  expect_identical(ts_br_norm_mean_pilot(A = 1, k = 2, w = 0.001), 32L)  # 31.6
  expect_identical(ts_br_norm_mean_pilot(A = 2, k = 2, w = 0.004), 23L)  # 22.4
  # k = 1: (B / w)^2 with B = A sqrt(2 / pi); gamma = 1 gives ceiling(B / w).
  # A = 2, w = 0.05: 2 * 0.7978846 / 0.05 = 31.92, so 32 (it was 1019).
  expect_identical(ts_br_norm_mean_pilot(A = 2, k = 1, w = 0.05), 32L)
  # Default minimum 10; never below k + 2.
  expect_identical(ts_br_norm_mean_pilot(A = 1, k = 2, w = 1), 10L)
  expect_identical(ts_br_norm_mean_pilot(A = 1, k = 2, w = 1, m_0 = 4), 4L)
  expect_identical(ts_br_norm_mean_pilot(A = 1, k = 1, w = 1, m_0 = 3), 3L)
  expect_error(ts_br_norm_mean_pilot(A = 1, k = 2, w = 1, m_0 = 3),
               "'m_0' must be a single whole number of at least 4")
  expect_error(ts_br_norm_mean_pilot(A = 1, k = 1, w = 1, m_0 = 2),
               "'m_0' must be a single whole number of at least 3")
  # gamma = 3: (1 / 2e-8)^(1/4) = 84.09, so 85; and at an exact power,
  # (1 / 1e-8)^(1/4) = 100, the answer is 100 and not 101.
  expect_identical(ts_br_norm_mean_pilot(A = 1, k = 2, w = 2e-8, gamma = 3),
                   85L)
  expect_identical(ts_br_norm_mean_pilot(A = 1, k = 2, w = 1e-8, gamma = 3),
                   100L)
  # Lower order than the unit scale optimum A / w.
  expect_lt(ts_br_norm_mean_pilot(A = 1, k = 2, w = 1e-6) / 1e6, 0.01)
  expect_error(ts_br_norm_mean_pilot(k = 1, w = 1), "'A' must be supplied")
  expect_error(ts_br_norm_mean_pilot(A = 1, k = 5, w = 1), "'k'")
  expect_error(ts_br_norm_mean_pilot(A = 1, k = 1, w = NA_real_), "'w'")
  expect_error(ts_br_norm_mean_pilot(A = 1, k = 1, w = 1, gamma = 0), "'gamma'")
  expect_error(ts_br_norm_mean_pilot(A = 1, k = 1, w = 1, m_0 = 10.5), "'m_0'")
  expect_error(ts_br_norm_mean_pilot(A = 1e300, k = 2, w = 1e-300),
               "implied pilot sample size .* is too large to represent")
  expect_error(ts_br_norm_mean_pilot(A = 1e300, k = 2, w = 1e-300),
               "Check that the loss constant 'A' and the risk bound 'w'")
  expect_identical(names(formals(ts_br_norm_mean_pilot)),
                   c("A", "k", "w", "m_0", "gamma"))
  # Every pilot it returns is accepted by the check function.
  withr::local_seed(113)
  for (k in 1:2) {
    m <- ts_br_norm_mean_pilot(A = 1, k = k, w = 1, m_0 = k + 2)
    expect_s3_class(ts_br_norm_mean(rnorm(m), A = 1, k = k, w = 1),
                    "smsd_check")
  }
})

test_that("a whole number pilot rule comes back as itself (k = 2, default gamma)", {
  # For squared error the documented rule is ceiling(sqrt(A / w)), the same
  # as for ts_br_exp_mean_pilot(). When A / w is a perfect square the answer
  # is its root, not the root plus one. (Before 1.0.0, and before B was
  # written in closed form, these gave 101, 14, and 26.)
  expect_identical(ts_br_norm_mean_pilot(A = 100, k = 2, w = 0.01), 100L)
  expect_identical(ts_br_norm_mean_pilot(A = 169, k = 2, w = 1), 13L)
  expect_identical(ts_br_norm_mean_pilot(A = 25, k = 2, w = 0.04), 25L)
  roots <- 4:400
  got <- vapply(roots, function(r) {
    ts_br_norm_mean_pilot(A = r^2, k = 2, w = 1, m_0 = 4)
  }, integer(1))
  expect_identical(got, roots)
  # The two bounded risk pilots document the same rule for k = 2 and agree.
  for (A in c(1, 2, 25, 100, 169, 1e4)) for (w in c(0.01, 0.04, 0.25, 1, 4)) {
    expect_identical(ts_br_norm_mean_pilot(A = A, k = 2, w = w, m_0 = 4),
                     ts_br_exp_mean_pilot(A = A, w = w, m_0 = 4))
  }
})

# ---- ts_br_exp_mean ------------------------------------------------------

test_that("ts_br_exp_mean criterion is Equations 8.3.9, 8.3.11, and 8.3.13", {
  x <- c(1, 2, 3, 4, 5)           # m = 5, mean 3
  # Equation 8.3.11: B_m(A) = 2 m (m + 1) A / ((m - 1)(m - 2)) = 60 A / 12
  # = 5 A. Equation 8.3.13: B* = 0.565 * 5 A = 2.825 A. Equation 8.3.9 with
  # A = 10, w = 2: N = ceiling(28.25 * 3^2 / 2) = ceiling(127.125) = 128.
  out <- ts_br_exp_mean(x, A = 10, w = 2)
  v <- stats::setNames(out$value, out$term)
  expect_equal(v[["optimal_n"]], 128)
  expect_equal(v[["mean"]], 3)
  expect_equal(v[["current_n"]], 5)
  expect_equal(v[["risk"]], 10 * 3^2 / 5)          # A xbar^2 / n = 18
  expect_false(attr(out, "is_satisfied"))
  expect_null(attr(out, "conf_level"))
  expect_equal(attr(out, "loss_A"), 10)
  expect_equal(attr(out, "risk_bound"), 2)
  expect_identical(out$term, c("mean", "risk", "optimal_n", "current_n"))
  expect_identical(attr(out, "method"),
                   "two stage bounded risk estimation (exponential mean)")
  expect_identical(attr(out, "current_n"), 5L)
  # m = 3, the smallest first stage: B* = 0.565 * 2 * 3 * 4 A / 2 = 6.78 A.
  # x = (1, 2, 3): mean 2; A = 1, w = 0.5: ceiling(6.78 * 4 / 0.5) = 55.
  expect_equal(attr(ts_br_exp_mean(c(1, 2, 3), A = 1, w = 0.5),
                    "criterion_n"), 55)
  # The estimate is the mean of both stages; N comes from the first alone.
  both <- ts_br_exp_mean(x, c(10, 20), A = 10, w = 2)
  expect_equal(both$value[both$term == "mean"], 45 / 7)
  expect_equal(attr(both, "criterion_n"), 128)
  # The risk row uses all seven observations and the current n:
  # A xbar^2 / n = 10 * (45 / 7)^2 / 7 = 59.04. (From the first stage mean
  # it would be 10 * 9 / 7 = 12.86; at the first stage size, 82.65.)
  expect_equal(both$value[both$term == "risk"], 10 * (45 / 7)^2 / 7)
  expect_equal(both$value[both$term == "risk"], 59.03790, tolerance = 1e-6)
  expect_equal(both$value[both$term == "current_n"], 7)
  expect_identical(attr(both, "current_n"), 7L)
  expect_true(attr(ts_br_exp_mean(x, rep(3, 123), A = 10, w = 2),
                   "is_satisfied"))
  expect_false(attr(ts_br_exp_mean(x, rep(3, 122), A = 10, w = 2),
                    "is_satisfied"))
  # A large bound: the first stage is already enough.
  expect_equal(attr(ts_br_exp_mean(x, A = 10, w = 1000), "criterion_n"), 5)
  # A constant positive first stage is valid here: N depends on its mean.
  expect_equal(attr(ts_br_exp_mean(rep(3, 5), A = 10, w = 2), "criterion_n"),
               128)
  # A criterion that is a whole number comes back as itself.
  # x = (4, 7, 8, 7, 6): m = 5, mean 32 / 5 = 6.4; B* = 2.825 A as above, so
  # with A = 1000, w = 2: N = 2825 * 6.4^2 / 2 = 2825 * 20.48 = 57856
  # exactly. (0.565 and 6.4 are not exact in floating point; the product
  # lands just above 57856 and a plain ceiling gave 57857.)
  expect_identical(attr(ts_br_exp_mean(c(4, 7, 8, 7, 6), A = 1000, w = 2),
                        "criterion_n"), 57856L)
  # x = (2, 2, 2, 2, 2): mean 2; A = 100, w = 1: 282.5 * 4 = 1130 exactly.
  expect_identical(attr(ts_br_exp_mean(rep(2, 5), A = 100, w = 1),
                        "criterion_n"), 1130L)
  # And one above a whole number still rounds up: A = 100.1 gives
  # 2.825 * 100.1 * 4 = 1131.13, so 1132.
  expect_identical(attr(ts_br_exp_mean(rep(2, 5), A = 100.1, w = 1),
                        "criterion_n"), 1132L)
})

test_that("whole number criteria come back as themselves (integer arithmetic)", {
  # Whole number first stages with round A and w; the exact criterion is
  # computed by integer arithmetic, with no rounding at all, and only the
  # cases in which it is a whole number are kept. Before the rounding
  # ignored representation error, about one such case in ten for
  # ts_br_norm_mean() came back one too high (the products and var() land a
  # few units in the last place above the whole number).
  withr::local_seed(113)
  pairs <- list(c(1, 1), c(2, 1), c(1, 2), c(3, 1), c(5, 1), c(10, 1),
                c(2, 5), c(4, 1), c(1, 4), c(25, 2))
  n_norm <- 0L; n_exp <- 0L
  for (i in 1:1500) {
    m <- sample(5:12, 1)
    x <- sample(1:7, m, replace = TRUE)
    if (max(x) == min(x)) next
    pq <- pairs[[sample(length(pairs), 1)]]
    A <- pq[1]; w <- pq[2]
    # ts_br_norm_mean, k = 2: ((m - 1) / (m - 3)) (A / w) S^2
    # = A (m sum x^2 - (sum x)^2) / (w m (m - 3)).
    num <- A * (m * sum(x^2) - sum(x)^2); den <- w * m * (m - 3)
    if (num %% den == 0) {
      n_norm <- n_norm + 1L
      expect_identical(attr(ts_br_norm_mean(x, A = A, k = 2, w = w),
                            "criterion_n"), as.integer(max(m, num / den)))
    }
    # ts_br_exp_mean: 0.565 * 2 m (m + 1) A / ((m - 1)(m - 2)) * xbar^2 / w
    # = 113 * 2 (m + 1) A (sum x)^2 / (200 (m - 1)(m - 2) m w).
    num <- 113 * 2 * (m + 1) * A * sum(x)^2
    den <- 200 * (m - 1) * (m - 2) * m * w
    if (num %% den == 0) {
      n_exp <- n_exp + 1L
      expect_identical(attr(ts_br_exp_mean(x, A = A, w = w), "criterion_n"),
                       as.integer(max(m, num / den)))
    }
  }
  # The sweep found whole number cases to check for both procedures.
  expect_gt(n_norm, 100L)
  expect_gt(n_exp, 2L)

  # Absolute error (k = 1) with an odd first stage, where (b B / A)^2 is a
  # rational number: with m = 2j + 1, Gamma(j - 1/2) = sqrt(pi) (2j - 2)! /
  # (4^(j - 1) (j - 1)!) and Gamma(j) = (j - 1)!, so
  # (b B / A)^2 = 2 j ((2j - 2)! / (4^(j - 1) ((j - 1)!)^2))^2, which is
  # 2, 1, 27 / 32, and 25 / 32 for m = 3, 5, 7, and 9. The criterion is
  # then (b B / A)^2 (A / w)^2 S^2
  # = p A^2 (m sum x^2 - (sum x)^2) / (q w^2 m (m - 1)) for (b B / A)^2 = p / q.
  # The function forms b B S / w (two square roots) before squaring it, so
  # the computed value is a few units in the last place off the whole
  # number either way.
  # (The list is looked up by sprintf("%d", m), which does not depend on the
  # session's display options as as.character() of a double does.)
  ratio <- list(`3` = c(2, 1), `5` = c(1, 1), `7` = c(27, 32), `9` = c(25, 32))
  for (m in c(3, 5, 7, 9)) {
    j <- (m - 1) / 2
    pq <- ratio[[sprintf("%d", m)]]
    expect_equal(2 * j * (factorial(2 * j - 2) /
                            (4^(j - 1) * factorial(j - 1)^2))^2,
                 pq[1] / pq[2])
  }
  n_abs <- 0L
  for (i in 1:1500) {
    m <- sample(c(3, 5, 7, 9), 1)
    x <- sample(1:9, m, replace = TRUE)
    if (max(x) == min(x)) next
    A <- sample(c(1, 2, 3, 4, 8, 16), 1); w <- sample(c(1, 2, 4), 1)
    pq <- ratio[[sprintf("%d", m)]]
    num <- pq[1] * A^2 * (m * sum(x^2) - sum(x)^2)
    den <- pq[2] * w^2 * m * (m - 1)
    if (num %% den == 0) {
      n_abs <- n_abs + 1L
      expect_identical(attr(ts_br_norm_mean(x, A = A, k = 1, w = w),
                            "criterion_n"), as.integer(max(m, num / den)))
    }
  }
  expect_gt(n_abs, 100L)
})

test_that("ts_br_exp_mean requires positive data and checks every argument", {
  x <- c(1, 2, 3, 4, 5)
  # Before 1.0.0 negative data were accepted.
  expect_error(ts_br_exp_mean(-x, A = 10, w = 2), "only positive values")
  expect_error(ts_br_exp_mean(c(x, 0), A = 10, w = 2), "only positive values")
  expect_error(ts_br_exp_mean(x, c(1, -1), A = 10, w = 2),
               "'additional_data' must contain only positive")
  expect_error(ts_br_exp_mean(c(1, 2), A = 10, w = 2),
               "'pilot_data' needs at least 3 observations; it has 2\\.")
  expect_error(ts_br_exp_mean(c(1, 2, NA), A = 10, w = 2),
               "'pilot_data' needs at least 3 observations; it has 2\\.")
  # The advice says what the minimum is for and is for a two stage design:
  # a larger first stage.
  expect_error(ts_br_exp_mean(c(1, 2), A = 10, w = 2),
               "it has 2\\. The procedure is defined for a first stage of at least three observations\\. Enlarge the first stage; ts_br_exp_mean_pilot\\(\\) returns a suitable size\\.")
  expect_false(grepl("check the stopping rule again",
                     tryCatch(ts_br_exp_mean(c(1, 2), A = 10, w = 2),
                              error = conditionMessage)))
  expect_error(ts_br_exp_mean(A = 10, w = 2), "'pilot_data' must be supplied")
  expect_error(ts_br_exp_mean(x, w = 2), "'A' must be supplied")
  expect_error(ts_br_exp_mean(x, A = 10), "'w' must be supplied")
  expect_error(ts_br_exp_mean(x, A = 0, w = 2), "'A'")
  expect_error(ts_br_exp_mean(x, A = NA_real_, w = 2), "'A'")
  expect_error(ts_br_exp_mean(x, A = c(1, 2), w = 2), "'A'")
  expect_error(ts_br_exp_mean(x, A = 10, w = -1), "'w'")
  expect_error(ts_br_exp_mean(x, A = 10, w = NA_real_), "'w'")
  expect_error(ts_br_exp_mean(x, A = 10, w = "a"), "'w'")
  expect_error(ts_br_exp_mean(letters, A = 10, w = 2), "'pilot_data' must be a numeric vector")
  expect_error(ts_br_exp_mean(c(x, NA), A = 10, w = 2, na.rm = FALSE), "missing values")
  expect_error(ts_br_exp_mean(x, A = 10, w = 2, na.rm = NA), "'na.rm'")
  expect_error(ts_br_exp_mean(x, A = 1e200, w = 1e-200),
               "implied total sample size .* is too large to represent")
  expect_error(ts_br_exp_mean(x, A = 1e200, w = 1e-200),
               "Check that the loss constant 'A' and the risk bound 'w'")
  expect_equal(ts_br_exp_mean(c(x, NA), A = 10, w = 2)$value,
               ts_br_exp_mean(x, A = 10, w = 2)$value)
})

test_that("ts_br_exp_mean keeps the risk at or under the bound (seeded Monte Carlo)", {
  withr::local_seed(113)
  m <- 10; theta <- 2; A <- 100; w <- 4; reps <- 2000
  n_star <- A * theta^2 / w                        # 100
  res <- replicate(reps, {
    first <- rexp(m, rate = 1 / theta)
    total <- attr(ts_br_exp_mean(first, A = A, w = w), "criterion_n")
    out <- ts_br_exp_mean(first, rexp(total - m, rate = 1 / theta), A = A,
                          w = w)
    stopifnot(isTRUE(attr(out, "is_satisfied")))
    c(n = total, loss = A * (out$value[out$term == "mean"] - theta)^2,
      # Risk given the first stage sum S: A ((S - m theta)^2 +
      # (N - m) theta^2) / N^2.
      cond = A * ((sum(first) - m * theta)^2 + (total - m) * theta^2) /
        total^2)
  })
  exact <- ts_br_exp_exact(m, n_star)
  expect_lte(exact[["risk_ratio"]], 1)             # the guarantee itself
  expect_equal(exact[["risk_ratio"]], 0.9152, tolerance = 5e-4)
  expect_equal(exact[["n_ratio"]], 1.904, tolerance = 5e-4)

  # Tolerances are four Monte Carlo standard errors.
  se_loss <- stats::sd(res["loss", ]) / sqrt(reps)
  expect_lte(mean(res["loss", ]), w + 4 * se_loss)
  expect_lt(abs(mean(res["loss", ]) / w - exact[["risk_ratio"]]),
            4 * se_loss / w)
  expect_lt(abs(mean(res["n", ]) / n_star - exact[["n_ratio"]]),
            4 * stats::sd(res["n", ]) / sqrt(reps) / n_star)
  se_cond <- stats::sd(res["cond", ]) / sqrt(reps)
  expect_lt(abs(mean(res["cond", ]) / w - exact[["risk_ratio"]]),
            4 * se_cond / w)
  # At or under the bound, with the same four standard errors (the exact
  # value 0.915 w is fewer than four of them below w).
  expect_lte(mean(res["cond", ]), w + 4 * se_cond)
})

test_that("the documented expected sample size and risk table for ts_br_exp_mean is right", {
  # Documented numbers. This block involves no package function: it checks
  # the table and the statements printed in ?ts_br_exp_mean_pilot and
  # ?ts_br_exp_mean against the exact calculation in helper-ts.R, so that
  # the help pages cannot drift from what the procedure does. The package
  # code is tied to the same exact calculation by the Monte Carlo test above.
  #
  # E[N] / n* ~ 0.565 * 2 (m + 1)^2 / ((m - 1)(m - 2)), the table in
  # ?ts_br_exp_mean_pilot, against exact calculation at a large n*.
  m <- c(3, 5, 10, 15, 20, 30, 50)
  approx <- 0.565 * 2 * (m + 1)^2 / ((m - 1) * (m - 2))
  expect_equal(round(approx, 2), c(9.04, 3.39, 1.90, 1.59, 1.46, 1.34, 1.25))
  exact <- vapply(m[-1], function(k) ts_br_exp_exact(k, 2000)[["n_ratio"]],
                  numeric(1))
  expect_equal(exact, approx[-1], tolerance = 2e-3)
  # The bound holds with little to spare where the risk is largest, and
  # with room at the extremes.
  for (worst in list(ts_br_exp_exact(14, 43), ts_br_exp_exact(15, 44.5))) {
    expect_equal(worst[["risk_ratio"]], 0.990, tolerance = 5e-4)
    expect_lt(worst[["risk_ratio"]], 1)
  }
  expect_lt(ts_br_exp_exact(3, 100)[["risk_ratio"]], 0.75)
  # Without the factor 0.565 (Theorem 8.3.1) the bound holds with a wide
  # margin at about twice the sample size.
  plain <- ts_br_exp_exact(20, 100, factor = 1)
  expect_lt(plain[["risk_ratio"]], 0.6)
  expect_gt(plain[["n_ratio"]], 2)
})

test_that("ts_br_exp_mean_pilot follows the documented rule and minimum", {
  # gamma = 1: ceiling(sqrt(A / w)).
  expect_identical(ts_br_exp_mean_pilot(A = 100, w = 0.01), 100L)
  expect_identical(ts_br_exp_mean_pilot(A = 100, w = 0.03), 58L)   # 57.7
  # Default minimum 20 (the old rule gave 3 here); never below 3.
  expect_identical(ts_br_exp_mean_pilot(A = 100, w = 4), 20L)
  expect_identical(ts_br_exp_mean_pilot(A = 100, w = 4, m_0 = 10), 10L)
  expect_identical(ts_br_exp_mean_pilot(A = 100, w = 4, m_0 = 3), 5L)
  expect_identical(ts_br_exp_mean_pilot(A = 1, w = 4, m_0 = 3), 3L)
  expect_error(ts_br_exp_mean_pilot(A = 100, w = 4, m_0 = 2),
               "'m_0' must be a single whole number of at least 3")
  # A perfect square comes back as its root (the rule uses sqrt() for the
  # default gamma, which is exact on every platform).
  expect_identical(vapply(3:300, function(r) {
    ts_br_exp_mean_pilot(A = r^2, w = 1, m_0 = 3)
  }, integer(1)), 3:300)
  # gamma = 0.25: 10000^(1 / 1.25) = 1584.9, so 1585.
  expect_identical(ts_br_exp_mean_pilot(A = 100, w = 0.01, gamma = 0.25),
                   1585L)
  # An exact power for a gamma other than 1 comes back as its root as well:
  # 32^(1 / 1.25) = 2^4 = 16 (computed as 16.000000000000004, which a plain
  # ceiling took to 17), and r^q to the power 1 / q is r.
  expect_identical(ts_br_exp_mean_pilot(A = 32, w = 1, m_0 = 3, gamma = 0.25),
                   16L)
  for (q in c(1.25, 1.5, 3, 4, 5)) {
    expect_identical(vapply(3:150, function(r) {
      ts_br_exp_mean_pilot(A = r^q, w = 1, m_0 = 3, gamma = q - 1)
    }, integer(1)), 3:150)
  }
  # Just above an exact power the rule still rounds up: 16.01^1.25 = 32.03.
  expect_identical(ts_br_exp_mean_pilot(A = 16.01^1.25, w = 1, m_0 = 3,
                                        gamma = 0.25), 17L)
  expect_lt(ts_br_exp_mean_pilot(A = 1e6, w = 1) / 1e6, 0.01)   # lower order
  expect_error(ts_br_exp_mean_pilot(w = 4), "'A' must be supplied")
  expect_error(ts_br_exp_mean_pilot(A = 100), "'w' must be supplied")
  expect_error(ts_br_exp_mean_pilot(A = -1, w = 4), "'A'")
  expect_error(ts_br_exp_mean_pilot(A = 100, w = NA_real_), "'w'")
  expect_error(ts_br_exp_mean_pilot(A = 100, w = 4, gamma = -1), "'gamma'")
  expect_error(ts_br_exp_mean_pilot(A = 100, w = 4, m_0 = NA), "'m_0'")
  expect_error(ts_br_exp_mean_pilot(A = 1e300, w = 1e-300),
               "implied pilot sample size .* is too large to represent")
  expect_error(ts_br_exp_mean_pilot(A = 1e300, w = 1e-300),
               "Check that the loss constant 'A' and the risk bound 'w'")
  expect_identical(names(formals(ts_br_exp_mean_pilot)),
                   c("A", "w", "m_0", "gamma"))
})
