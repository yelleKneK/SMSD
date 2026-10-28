# Residual hardening from the adversarial post-fix verification pass that is
# not covered in a family's own test file: argument-naming and formals pins
# for the renames, input validation of the Gini helpers, the linear-criterion
# nexp pilots, the max_n censoring row, and authority checks for the
# remaining variance estimators. (The ceiling pins for the classical
# criteria, the bounded risk row at k = 2, the negative exponential location
# row, and the two-sided pilot pins now live, against independent values, in
# the classical families' own test files.)

test_that("error messages name the caller's actual argument", {
  expect_error(seq_mr_cv("a", A = 400, c = 1), "'x' must be a numeric")
  expect_error(seq_mr_smd(1, rnorm(10), A = 100, c = 4),
               "'group_1' needs at least")
  expect_error(seq_mr_smd(rnorm(10), "b", A = 100, c = 4),
               "'group_2' must be a numeric")
})

test_that("seq_mr_smd formals are group_1 and group_2", {
  expect_identical(names(formals(seq_mr_smd))[1:2], c("group_1", "group_2"))
})

test_that("gini helpers validate input, including a sample of zeros", {
  expect_error(gini_index(rep(0, 20)), "Every value in 'x' is zero")
  expect_error(gini_index(letters[1:5]), "numeric")
  expect_error(gini_index(3), "at least 2")
  expect_error(gini_mean_diff(letters[1:5]), "numeric")
  expect_error(gini_mean_diff(3), "at least 2")
})

test_that("nexp pilots use the linear criterion constant -log(alpha)", {
  a <- -log(0.05)
  # sqrt(2.995732 / 0.2) = 3.87 -> 4 (a normal quantile with a quadratic
  # exponent would give 1.959964 / 0.2 = 9.80 -> 10).
  expect_identical(seq_fwd_nexp_mean_pilot(d = 0.2, conf_level = 0.95,
                                           m_0 = 2), 4L)
  # (m_0 = 2 so that the rule, not the default minimum of 10, is what is
  # pinned: log(1 / 0.05) = 2.995732 and sqrt(2.995732 / 0.2) = 3.87, so 4.
  # With the normal quantile in its place, sqrt(1.959964 / 0.2) = 3.13 would
  # also give 4, so a narrower target separates them: d = 0.02 gives
  # sqrt(149.79) = 12.24, so 13, against sqrt(98.0) = 9.90, so 10.)
  expect_identical(ts_fwd_nexp_mean_pilot(d = 0.2, conf_level = 0.95, m_0 = 2),
                   4L)
  expect_identical(ts_fwd_nexp_mean_pilot(d = 0.02, conf_level = 0.95,
                                          m_0 = 2), 13L)
})

test_that("pct_hit_max_n reports censoring by the cap", {
  capped <- seq_aipe_smd_sensitivity(true_value = 0.5, width = 0.2,
                                     G = 3, seed = 113, max_n = 25)
  expect_equal(capped$value[capped$term == "pct_hit_max_n"], 100)
  free <- seq_aipe_smd_sensitivity(true_value = 0.5, width = 0.8,
                                   G = 10, seed = 113)
  expect_equal(free$value[free$term == "pct_hit_max_n"], 0)
})

test_that(".var_cv and .var_sm match their normal theory values", {
  # Standard deviations of the estimates (400 samples each): 0.00072 for
  # .var_cv at n = 20,000, so its relative tolerance of 0.1 (0.0035 in
  # absolute terms) is 4.9 of them (at the n of 6,000 used before it was
  # 2.5); and 0.017 for .var_sm at n = 6,000, so its tolerance (0.1125) is
  # 6.6 of them.
  withr::local_seed(113)
  z <- rnorm(20000, 4, 1)         # kappa = 0.25
  expect_equal(.var_cv(z), 0.25^2 * (0.25^2 + 0.5), tolerance = 0.1)
  withr::local_seed(113)
  w <- rnorm(6000, 0.5, 1)        # theta = mu/sigma = 0.5
  expect_equal(.var_sm(w), 1 + 0.5^2 / 2, tolerance = 0.1)
})

test_that(".var_gini matches theory on fixed data and .var_c is the weighted sum", {
  # Exponential data: the asymptotic variance of the Gini index is 1/12
  # (Chattopadhyay and De, 2016, Equation 6, with Delta = 1, tau = 3/2,
  # sigma_1^2 = 1/3). The standard deviation of the estimate at n = 6,000
  # is 0.0029 (400 samples), so 0.014 is 4.9 of them. (This test used to
  # pin the function's own output, 0.0855750490.)
  withr::local_seed(113)
  x <- rgamma(6000, shape = 1, rate = 1)
  expect_lt(abs(.var_gini(x) - 1 / 12), 0.014)
  withr::local_seed(113)
  g <- cbind(rnorm(30, 1), rnorm(30, 2), rnorm(30, 3))
  wts <- c(1, -0.5, -0.5)
  expect_equal(.var_c(g, wts),
               sum(wts^2 * apply(g, 2, var)))
})
