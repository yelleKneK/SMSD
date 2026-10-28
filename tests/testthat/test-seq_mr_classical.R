# Tests for the classical purely sequential minimum risk point estimation
# procedures: seq_mr_norm_mean(), seq_mr_exp_mean(), seq_mr_nexp_mean(),
# seq_mr_dfr_mean_arb(), and their pilots. Each criterion is checked against
# the risk function it is meant to minimize (written out in the test, with
# expectations by quadrature where a distribution is involved), not against
# the package's closed form typed again.

# The criterion must be the integer just above the minimizer of a risk
# function: the risk is still falling one step below it and no longer
# falling at it. `risk` is a function of n.
expect_minimizer <- function(n_opt, risk) {
  slope <- function(n) (risk(n * (1 + 1e-6)) - risk(n * (1 - 1e-6))) / (2e-6 * n)
  expect_gte(slope(n_opt), 0)
  expect_lt(slope(n_opt - 1), 0)
}

# ---- pilots -----------------------------------------------------------------

test_that("the pilots reproduce the worked example of Chattopadhyay and Kelley (2016)", {
  # Their Equation 24, m = max{m_0, ceiling((A / c)^(1 / (2 + 2 gamma)))}, is
  # illustrated on p. 637 with A = 400,000, a sampling cost of 75, and
  # gamma = 0.49: the pilot is 18. Stated through price and epsilon
  # (1,000 and 0.05) it is again 18.
  expect_identical(seq_mr_norm_mean_pilot(A = 400000, c = 75, gamma = 0.49), 18L)
  expect_identical(seq_mr_norm_mean_pilot(c = 75, price = 1000, epsilon = 0.05,
                                          gamma = 0.49), 18L)
  expect_identical(seq_mr_exp_mean_pilot(A = 400000, c = 75, gamma = 0.49), 18L)
  expect_identical(seq_mr_dfr_mean_arb_pilot(A = 400000, c = 75, gamma = 0.49), 18L)
})

test_that("the pilots equal max(m_0, ceiling(K^(1 / (1 + gamma)))) by hand", {
  # Default gamma = 1: (A / c)^(1 / 4). A = 100, c = 0.004: 25000^(1/4)
  # = 12.57 -> 13. A = 100, c = 0.04: 2500^(1/4) = 7.07 -> 8.
  expect_identical(seq_mr_norm_mean_pilot(A = 100, c = 0.004), 13L)
  expect_identical(seq_mr_exp_mean_pilot(A = 100, c = 0.04), 8L)
  expect_identical(seq_mr_dfr_mean_arb_pilot(A = 100, c = 0.004), 13L)
  # A = 100, c = 0.5: 200^(1/4) = 3.76 -> 4 (the rule itself, not the
  # minimum; the minimum is tested below).
  expect_identical(seq_mr_norm_mean_pilot(A = 100, c = 0.5), 4L)
  expect_identical(seq_mr_norm_mean_pilot(A = 100, c = 0.5, m_0 = 9), 9L)
  # Negative exponential: K = (A k Gamma(k + 1) / (c t))^(1 / (t + k)).
  # k = 2, t = 1: (500 * 2 * 2 / 0.01)^(1/3) = 58.48, sqrt = 7.65 -> 8.
  # k = 1, t = 1: (500 / 0.01)^(1/2) = 223.61, sqrt = 14.95 -> 15.
  # k = 1, t = 2: (500 / 0.02)^(1/3) = 29.24, sqrt = 5.41 -> 6.
  expect_identical(seq_mr_nexp_mean_pilot(A = 500, c = 0.01, k = 2, t = 1), 8L)
  expect_identical(seq_mr_nexp_mean_pilot(A = 500, c = 0.01, k = 1, t = 1), 15L)
  expect_identical(seq_mr_nexp_mean_pilot(A = 500, c = 0.01, k = 1, t = 2), 6L)
  # gamma = 0.49 with k = 2, t = 1: 58.48^(1 / 1.49) = 15.34 -> 16.
  expect_identical(seq_mr_nexp_mean_pilot(A = 500, c = 0.01, k = 2, t = 1,
                                          gamma = 0.49), 16L)
})

test_that("the default minimum pilot of the minimum risk procedures is 4", {
  # Cases in which the minimum governs. A / c = 16: 16^(1/4) = 2, so the
  # rule alone gives 2 and the default m_0 = 4 gives 4. Negative exponential
  # with k = t = 1, A = 4, c = 1: K = sqrt(4) = 2 and sqrt(K) = 1.41 -> 2.
  expect_identical(seq_mr_norm_mean_pilot(A = 16, c = 1), 4L)
  expect_identical(seq_mr_exp_mean_pilot(A = 16, c = 1), 4L)
  expect_identical(seq_mr_dfr_mean_arb_pilot(A = 16, c = 1), 4L)
  expect_identical(seq_mr_nexp_mean_pilot(A = 4, c = 1, k = 1, t = 1), 4L)
  # With the minimum lowered the rule's own value shows through, and a
  # larger minimum is honored.
  expect_identical(seq_mr_norm_mean_pilot(A = 16, c = 1, m_0 = 2), 2L)
  expect_identical(seq_mr_exp_mean_pilot(A = 16, c = 1, m_0 = 2), 2L)
  expect_identical(seq_mr_dfr_mean_arb_pilot(A = 16, c = 1, m_0 = 2), 2L)
  expect_identical(seq_mr_nexp_mean_pilot(A = 4, c = 1, k = 1, t = 1, m_0 = 2), 2L)
  expect_identical(seq_mr_exp_mean_pilot(A = 16, c = 1, m_0 = 7), 7L)
  expect_identical(seq_mr_dfr_mean_arb_pilot(A = 16, c = 1, m_0 = 7), 7L)
  expect_identical(seq_mr_nexp_mean_pilot(A = 4, c = 1, k = 1, t = 1, m_0 = 7), 7L)
})

test_that("the pilots are of lower order than the risk minimizing sample size", {
  # With A = 100 and c = 0.01, n* = 100 sigma: 25 at sigma = 0.25. The old
  # pilot (n* at sigma = 1) was 100, four times too many, with risk 2.1 times
  # the minimum. The pilot must stay below 25.
  expect_lt(seq_mr_norm_mean_pilot(A = 100, c = 0.01), 25)
  expect_lt(seq_mr_exp_mean_pilot(A = 100, c = 0.01), 25)
  expect_lt(seq_mr_dfr_mean_arb_pilot(A = 100, c = 0.01), 25)
  # Negative exponential, k = 2, t = 1: n* = 58.48 sigma^(2/3) = 23.2 at
  # sigma = 0.25.
  expect_lt(seq_mr_nexp_mean_pilot(A = 500, c = 0.01, k = 2, t = 1), 23.2)
  cost <- c(1, 0.01, 0.0001)
  frac <- vapply(cost, function(ci) seq_mr_norm_mean_pilot(A = 100, c = ci) /
                   sqrt(100 / ci), numeric(1))
  expect_true(all(diff(frac) < 0))
})

test_that("the pilots check every argument", {
  expect_error(seq_mr_norm_mean_pilot(c = 1), "Supply either 'A'")
  expect_error(seq_mr_norm_mean_pilot(A = 1, c = 1, price = 2, epsilon = 1), "not both")
  expect_error(seq_mr_norm_mean_pilot(c = 1, price = 2), "Supply either 'A'")
  expect_error(seq_mr_norm_mean_pilot(c = 1, price = -2, epsilon = 1), "'price'")
  pilots <- list(
    function(...) seq_mr_norm_mean_pilot(...),
    function(...) seq_mr_exp_mean_pilot(...),
    function(...) seq_mr_dfr_mean_arb_pilot(...),
    function(...) seq_mr_nexp_mean_pilot(..., k = 1, t = 1))
  for (pilot in pilots) {
    expect_error(pilot(A = NA_real_, c = 1), "'A' must be")
    expect_error(pilot(A = -1, c = 1), "'A' must be")
    expect_error(pilot(A = c(1, 2), c = 1), "'A' must be")
    expect_error(pilot(A = 100, c = 0), "'c' must be")
    expect_error(pilot(A = 100, c = NA_real_), "'c' must be")
    expect_error(pilot(A = 100, c = 1, m_0 = 1), "'m_0'")
    expect_error(pilot(A = 100, c = 1, m_0 = 3.5), "'m_0'")
    expect_error(pilot(A = 100, c = 1, gamma = 0), "'gamma'")
    expect_error(pilot(A = 100, c = 1, gamma = NA_real_), "'gamma'")
    expect_error(pilot(A = 1e300, c = 1e-300),
                 "implied pilot sample size .* is too large to represent")
    expect_type(pilot(A = 100, c = 1), "integer")
  }
  # The advice that closes the too large message names what to check.
  expect_error(seq_mr_norm_mean_pilot(A = 1e300, c = 1e-300),
               "represent\\. Check the loss constant \\('A', or 'price' and 'epsilon'\\) and the cost 'c'\\.$")
  expect_error(seq_mr_exp_mean_pilot(A = 1e300, c = 1e-300),
               "represent\\. Check 'A' and the cost 'c'\\.$")
  expect_error(seq_mr_dfr_mean_arb_pilot(A = 1e300, c = 1e-300),
               "represent\\. Check 'A' and the cost 'c'\\.$")
  expect_error(seq_mr_nexp_mean_pilot(A = 1e300, c = 1e-300, k = 1, t = 1),
               "represent\\. Check 'A', 'c', 'k', and 't'\\.$")
  # Argument order: the method's own arguments, then m_0, then gamma.
  expect_identical(names(formals(seq_mr_norm_mean_pilot)),
                   c("A", "c", "price", "epsilon", "m_0", "gamma"))
  expect_identical(names(formals(seq_mr_exp_mean_pilot)),
                   c("A", "c", "m_0", "gamma"))
  expect_identical(names(formals(seq_mr_dfr_mean_arb_pilot)),
                   c("A", "c", "m_0", "gamma"))
  expect_identical(names(formals(seq_mr_nexp_mean_pilot)),
                   c("A", "c", "k", "t", "m_0", "gamma"))
  expect_error(seq_mr_nexp_mean_pilot(A = 500, c = 0.01, k = 0, t = 1), "'k' must be")
  expect_error(seq_mr_nexp_mean_pilot(A = 500, c = 0.01, k = 1, t = NA_real_), "'t' must be")
})

test_that("a required argument left out is named, in the package's own words", {
  set.seed(113)
  x <- rexp(30)
  expect_error(seq_mr_norm_mean_pilot(A = 100), "'c' must be supplied")
  expect_error(seq_mr_exp_mean_pilot(c = 1), "'A' must be supplied")
  expect_error(seq_mr_exp_mean_pilot(A = 100), "'c' must be supplied")
  expect_error(seq_mr_dfr_mean_arb_pilot(c = 1), "'A' must be supplied")
  expect_error(seq_mr_dfr_mean_arb_pilot(A = 100), "'c' must be supplied")
  expect_error(seq_mr_nexp_mean_pilot(c = 1, k = 1, t = 1), "'A' must be supplied")
  expect_error(seq_mr_nexp_mean_pilot(A = 100, k = 1, t = 1), "'c' must be supplied")
  expect_error(seq_mr_nexp_mean_pilot(A = 100, c = 1, t = 1), "'k' must be supplied")
  expect_error(seq_mr_nexp_mean_pilot(A = 100, c = 1, k = 1), "'t' must be supplied")
  expect_error(seq_mr_norm_mean(x, A = 100), "'c' must be supplied")
  expect_error(seq_mr_exp_mean(x, c = 1), "'A' must be supplied")
  expect_error(seq_mr_exp_mean(x, A = 100), "'c' must be supplied")
  expect_error(seq_mr_dfr_mean_arb(x, c = 1), "'A' must be supplied")
  expect_error(seq_mr_dfr_mean_arb(x, A = 100), "'c' must be supplied")
  expect_error(seq_mr_nexp_mean(x, c = 1, k = 1, t = 1), "'A' must be supplied")
  expect_error(seq_mr_nexp_mean(x, A = 100, k = 1, t = 1), "'c' must be supplied")
  expect_error(seq_mr_nexp_mean(x, A = 100, c = 1, t = 1), "'k' must be supplied")
  expect_error(seq_mr_nexp_mean(x, A = 100, c = 1, k = 1), "'t' must be supplied")
  for (f in list(seq_mr_norm_mean, seq_mr_exp_mean, seq_mr_dfr_mean_arb)) {
    expect_error(f(A = 100, c = 1), "'x' must be supplied")
  }
  expect_error(seq_mr_nexp_mean(A = 100, c = 1, k = 1, t = 1), "'x' must be supplied")
})

test_that("the data argument is x, as for every one sample procedure", {
  set.seed(113)
  obs <- rexp(30)
  expect_identical(names(formals(seq_mr_norm_mean)),
                   c("x", "A", "c", "price", "epsilon", "na.rm"))
  expect_identical(names(formals(seq_mr_exp_mean)), c("x", "A", "c", "na.rm"))
  expect_identical(names(formals(seq_mr_nexp_mean)),
                   c("x", "A", "c", "k", "t", "na.rm"))
  expect_identical(names(formals(seq_mr_dfr_mean_arb)),
                   c("x", "A", "c", "gamma", "na.rm"))
  # By name and by position give the same object; the name used before
  # 1.0.0 is not accepted (R's own "unused argument" error, whose wording is
  # not matched because R translates it).
  expect_identical(seq_mr_norm_mean(x = obs, A = 100, c = 1),
                   seq_mr_norm_mean(obs, 100, 1))
  expect_identical(seq_mr_exp_mean(x = obs, A = 100, c = 1),
                   seq_mr_exp_mean(obs, 100, 1))
  expect_identical(seq_mr_nexp_mean(x = obs, A = 100, c = 1, k = 1, t = 1),
                   seq_mr_nexp_mean(obs, 100, 1, 1, 1))
  expect_identical(seq_mr_dfr_mean_arb(x = obs, A = 100, c = 1),
                   seq_mr_dfr_mean_arb(obs, 100, 1))
  expect_error(seq_mr_norm_mean(data = obs, A = 100, c = 1))
  expect_error(seq_mr_exp_mean(data = obs, A = 100, c = 1))
  expect_error(seq_mr_nexp_mean(data = obs, A = 100, c = 1, k = 1, t = 1))
  expect_error(seq_mr_dfr_mean_arb(data = obs, A = 100, c = 1))
})

# ---- seq_mr_norm_mean -------------------------------------------------------

test_that("seq_mr_norm_mean's criterion minimizes A sigma^2 / n + c n", {
  # Mukhopadhyay and de Silva (2009), Equations 6.4.2 and 6.4.3, with the
  # sample variance in place of sigma^2.
  set.seed(113)
  for (cfg in list(c(A = 100, c = 0.5, v = 4), c(A = 2500, c = 4, v = 150),
                   c(A = 7, c = 0.03, v = 0.2))) {
    A <- cfg[["A"]]; cost <- cfg[["c"]]; v <- cfg[["v"]]
    x <- seqc_scale_var(rnorm(30, 70, 1), v)
    out <- seq_mr_norm_mean(x, A = A, c = cost)
    n_opt <- out$value[out$term == "optimal_n"]
    risk <- function(n) A * v / n + cost * n
    expect_minimizer(n_opt, risk)
    expect_equal(out$value[out$term == "risk"], risk(30))
    expect_identical(attr(out, "is_satisfied"), 30 >= n_opt)
  }
  # By hand: A = 100, c = 1, variance 1.1449 (sd 1.07): 10 * 1.07 = 10.7 -> 11.
  y <- seqc_scale_var(rnorm(10), 1.07^2)
  out <- seq_mr_norm_mean(y, A = 100, c = 1)
  expect_equal(out$value[out$term == "optimal_n"], 11)
  expect_false(attr(out, "is_satisfied"))           # 10 < 11: rounds up
})

test_that("seq_mr_norm_mean accepts A or price and epsilon, and returns the documented object", {
  set.seed(113)
  x <- rnorm(40, 70, 12)
  out <- seq_mr_norm_mean(x, A = 2500, c = 4)
  expect_s3_class(out, c("smsd_check", "data.frame"), exact = TRUE)
  expect_identical(out$term, c("mean", "risk", "optimal_n", "current_n"))
  expect_equal(out$value[c(1, 4)], c(mean(x), 40))
  expect_null(attr(out, "conf_level"))
  expect_null(attr(out, "width"))
  expect_identical(attr(out, "loss_A"), 2500)
  expect_identical(attr(out, "cost_c"), 4)
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_identical(attr(out, "method"), "minimum risk point estimation (normal mean)")
  # price / epsilon^2 = 10000 / 4 = 2500
  expect_equal(seq_mr_norm_mean(x, c = 4, price = 10000, epsilon = 2)$value, out$value)
  expect_error(seq_mr_norm_mean(x, A = 1, c = 4, price = 2, epsilon = 1), "not both")
  expect_error(seq_mr_norm_mean(x, c = 4), "Supply either 'A'")
})

test_that("a loss constant that price and epsilon cannot give is reported in their terms", {
  # A = price / epsilon^2 overflows for a tiny epsilon and underflows to zero
  # for a huge one. The user passed 'price' and 'epsilon', so the message
  # names them, with the value derived, and not 'A', which was not supplied.
  set.seed(113)
  x <- rnorm(20, 70, 12)
  expect_error(
    seq_mr_norm_mean(x, c = 1, price = 10, epsilon = 1e-200),
    paste0("^The loss constant implied by 'price' and 'epsilon', ",
           "price / epsilon\\^2, must be a positive finite number; it is ",
           "Inf for price = 10 and epsilon = 1e-200\\. Check both arguments ",
           "and their units\\.$")
  )
  expect_error(
    seq_mr_norm_mean_pilot(c = 1, price = 10, epsilon = 1e200),
    "price / epsilon\\^2, must be a positive finite number; it is 0 for price = 10 and epsilon = 1e\\+200"
  )
  for (call in list(
    function() seq_mr_norm_mean(x, c = 1, price = 10, epsilon = 1e-200),
    function() seq_mr_norm_mean(x, c = 1, price = 1e-200, epsilon = 1e100),
    function() seq_mr_norm_mean_pilot(c = 1, price = 1e200, epsilon = 1e-100)
  )) {
    msg <- tryCatch(call(), error = conditionMessage)
    expect_false(grepl("'A'", msg, fixed = TRUE))
  }
  # A supplied directly is still reported as 'A', and a bad price or epsilon
  # under its own name.
  expect_error(seq_mr_norm_mean(x, A = Inf, c = 1),
               "^'A' must be a single finite number greater than 0; got Inf")
  expect_error(seq_mr_norm_mean(x, c = 1, price = 10, epsilon = 0),
               "^'epsilon' must be a single finite number greater than 0")
  # The largest and smallest that can be represented still pass:
  # 1e300 / (1e-4)^2 = 1e308, and 1e-300 / (1e4)^2 = 1e-308.
  expect_identical(.mr_loss_constant(NULL, 1e300, 1e-4), 1e300 / 1e-4^2)
  expect_identical(.mr_loss_constant(NULL, 1e-300, 1e4), 1e-300 / 1e4^2)
  expect_identical(.mr_loss_constant(NULL, 10000, 2), 2500)
  expect_identical(.mr_loss_constant(2500, NULL, NULL), 2500)
})

# ---- seq_mr_exp_mean --------------------------------------------------------

test_that("seq_mr_exp_mean's criterion minimizes A lambda^2 / n + c n", {
  # The variance of an exponential variable with mean lambda is lambda^2
  # (checked here by quadrature), so the mean of n has variance lambda^2 / n.
  lambda <- 2.3
  v <- stats::integrate(function(x) (x - lambda)^2 * stats::dexp(x, 1 / lambda),
                        0, Inf)$value
  expect_equal(v, lambda^2, tolerance = 1e-6)
  set.seed(113)
  x <- rexp(40, rate = 0.5)
  for (cfg in list(c(A = 100, c = 0.04), c(A = 3, c = 0.5))) {
    A <- cfg[["A"]]; cost <- cfg[["c"]]
    out <- seq_mr_exp_mean(x, A = A, c = cost)
    n_opt <- out$value[out$term == "optimal_n"]
    risk <- function(n) A * mean(x)^2 / n + cost * n
    expect_minimizer(n_opt, risk)
    expect_equal(out$value[out$term == "risk"], risk(40))
    expect_identical(attr(out, "is_satisfied"), 40 >= n_opt)
  }
  out <- seq_mr_exp_mean(x, A = 100, c = 0.04)
  expect_identical(out$term, c("mean", "risk", "optimal_n", "current_n"))
  expect_equal(out$value[1], mean(x))
  expect_null(attr(out, "conf_level"))
  expect_identical(attr(out, "method"),
                   "minimum risk point estimation (exponential mean)")
  # By hand: three observations with mean 0.5, A = 100, c = 0.25:
  # sqrt(400) * 0.5 = 10; risk at n = 3 is 100 * 0.25 / 3 + 0.75.
  out <- seq_mr_exp_mean(c(0.25, 0.5, 0.75), A = 100, c = 0.25)
  expect_equal(out$value, c(0.5, 25 / 3 + 0.75, 10, 3))
  # A single observation is accepted (the rule needs only the mean).
  expect_false(attr(seq_mr_exp_mean(3, A = 4, c = 1), "is_satisfied"))
})

test_that("seq_mr_exp_mean requires positive data", {
  set.seed(113)
  x <- rexp(30)
  # Negative data used to give a negative optimal n and 'met'.
  expect_error(seq_mr_exp_mean(-x, A = 100, c = 0.04),
               "'x' must contain only positive values")
  expect_error(seq_mr_exp_mean(c(x, 0), A = 100, c = 0.04), "only positive values")
  expect_error(seq_mr_exp_mean(rep(0, 5), A = 100, c = 1), "only positive values")
})

# ---- seq_mr_nexp_mean -------------------------------------------------------

test_that("seq_mr_nexp_mean's criterion minimizes A E(X_(1) - theta)^k + c n^t", {
  # X_(1) - theta is exponential with scale sigma / n; its k-th moment is
  # taken by quadrature, with sigma replaced by U = sum(x - min x) / (n - 1).
  set.seed(113)
  x <- 3 + rexp(40, rate = 0.5)
  U <- sum(x - min(x)) / 39
  for (kt in list(c(1, 1), c(2, 1), c(1, 2), c(2, 2), c(2.5, 1.5), c(0.5, 0.7))) {
    k <- kt[1]; t <- kt[2]
    A <- 500; cost <- 0.01
    out <- seq_mr_nexp_mean(x, A = A, c = cost, k = k, t = t)
    n_opt <- out$value[out$term == "optimal_n"]
    risk <- function(n) {
      A * stats::integrate(function(y) y^k * stats::dexp(y, rate = n / U),
                           0, Inf, rel.tol = 1e-10)$value + cost * n^t
    }
    expect_minimizer(n_opt, risk)
    expect_equal(out$value[out$term == "risk"], risk(40), tolerance = 1e-6)
    expect_identical(attr(out, "is_satisfied"), 40 >= n_opt)
  }
})

test_that("seq_mr_nexp_mean reports the sample minimum and survives a large k", {
  set.seed(113)
  x <- 3 + rexp(40, rate = 1)
  out <- seq_mr_nexp_mean(x, A = 500, c = 0.01, k = 1, t = 2)
  expect_identical(out$term, c("location", "risk", "optimal_n", "current_n"))
  expect_equal(out$value[out$term == "location"], min(x))
  expect_null(attr(out, "conf_level"))
  expect_identical(attr(out, "k"), 1)
  expect_identical(attr(out, "t"), 2)
  expect_identical(attr(out, "method"),
                   "minimum risk point estimation (negative exponential location)")
  # A negative location is legitimate: shifting the data shifts only the
  # location row.
  shifted <- seq_mr_nexp_mean(x - 10, A = 500, c = 0.01, k = 1, t = 2)
  expect_equal(shifted$value, out$value - c(10, 0, 0, 0))

  # k = 150 with a scale estimate below 1: the closed form
  # (A k Gamma(k + 1) U^k / (c t))^(1 / (t + k)) is still computable directly.
  y <- rexp(30, rate = 2)
  U <- sum(y - min(y)) / 29
  direct <- (500 * 150 * gamma(151) * U^150 / (0.01 * 2))^(1 / 152)
  expect_true(is.finite(direct))
  out <- seq_mr_nexp_mean(y, A = 500, c = 0.01, k = 150, t = 2)
  expect_equal(out$value[out$term == "optimal_n"], ceiling(direct))
  # k = 200 overflows gamma(); the function works on the log scale, so the
  # result is still a definite answer (it used to be NA).
  out <- seq_mr_nexp_mean(y, A = 500, c = 0.01, k = 200, t = 2)
  expect_true(is.finite(out$value[out$term == "optimal_n"]))
  expect_gt(out$value[out$term == "optimal_n"], 0)
  expect_false(is.na(attr(out, "is_satisfied")))
})

# ---- seq_mr_dfr_mean_arb ----------------------------------------------------

test_that("seq_mr_dfr_mean_arb adds n^(-gamma) to the standard deviation", {
  set.seed(113)
  # By hand: A = 100, c = 1, sd 2, n = 30. gamma = 0.5: 10 * (2 + 30^-0.5)
  # = 21.83 -> 22. gamma = 1 (the default): 10 * (2 + 1 / 30) = 20.33 -> 21.
  x <- seqc_scale_var(rnorm(30, 50, 1), 4)
  expect_equal(seq_mr_dfr_mean_arb(x, A = 100, c = 1, gamma = 0.5)$value[3], 22)
  expect_equal(seq_mr_dfr_mean_arb(x, A = 100, c = 1, gamma = 1)$value[3], 21)
  expect_equal(seq_mr_dfr_mean_arb(x, A = 100, c = 1)$value[3], 21)
  # The normal theory rule has no correction term: 10 * 2 = 20 exactly at
  # most, so the distribution free rule asks for at least as many.
  expect_lte(seq_mr_norm_mean(x, A = 100, c = 1)$value[3], 21)
  # The criterion minimizes A s_eff^2 / n + c n with s_eff = S + n^(-gamma).
  for (g in c(0.25, 0.49, 1, 2)) {
    out <- seq_mr_dfr_mean_arb(x, A = 300, c = 0.7, gamma = g)
    s_eff <- 2 + 30^(-g)
    expect_minimizer(out$value[3], function(n) 300 * s_eff^2 / n + 0.7 * n)
    # The reported risk uses the sample variance itself.
    expect_equal(out$value[2], 300 * 4 / 30 + 0.7 * 30)
    expect_identical(attr(out, "gamma"), g)
    expect_identical(attr(out, "is_satisfied"), 30 >= out$value[3])
  }
  out <- seq_mr_dfr_mean_arb(x, A = 100, c = 1)
  expect_identical(out$term, c("mean", "risk", "optimal_n", "current_n"))
  expect_null(attr(out, "conf_level"))
  expect_identical(
    attr(out, "method"),
    "minimum risk point estimation (distribution free mean, arbitrary variance)")
})

# ---- degenerate and invalid input -------------------------------------------

test_that("constant data are 'not met' with a note, never an optimal n of 0", {
  checks <- list(
    function(x) seq_mr_norm_mean(x, A = 100, c = 1),
    function(x) seq_mr_nexp_mean(x, A = 500, c = 0.01, k = 1, t = 2),
    function(x) seq_mr_dfr_mean_arb(x, A = 100, c = 2, gamma = 0.5))
  for (check in checks) {
    for (x in list(rep(1, 2), rep(2, 5), rep(0.1, 30))) {
      out <- check(x)
      expect_identical(attr(out, "is_satisfied"), FALSE)
      expect_identical(attr(out, "variance_unusable"), TRUE)
      expect_true(all(is.na(out$value[out$term %in% c("risk", "optimal_n")])))
      # Exactly equal observations: the package's standard note, alone.
      expect_length(attr(out, "notes"), 1L)
      expect_match(attr(out, "notes"),
                   "not yet usable .*zero, negative, or not finite")
    }
  }
})

test_that("data that are constant up to rounding error are 'not met' as well", {
  # Ten gain scores of 0.3, each a posttest minus a pretest. On paper they
  # are equal; in floating point they differ by about 4e-15 and have a
  # "variance" near 2e-30, which satisfied every rule with a criterion of 1.
  pre <- c(101, 57, 12, 80, 33, 64, 90, 45, 71, 28)
  gain <- (pre + 0.3) - pre
  expect_gt(stats::var(gain), 0)
  expect_lt(max(gain) - min(gain), 1e-14)
  checks <- list(
    function(x) seq_mr_norm_mean(x, A = 2500, c = 4),
    function(x) seq_mr_nexp_mean(x, A = 500, c = 0.25, k = 1, t = 1),
    function(x) seq_mr_dfr_mean_arb(x, A = 400, c = 25))
  for (check in checks) {
    out <- check(gain)
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_identical(attr(out, "variance_unusable"), TRUE)
    expect_true(all(is.na(out$value[out$term %in% c("risk", "optimal_n")])))
    # Two notes: the standard one, and one that says why observations that
    # are not all equal were treated as constant and what to do about it.
    expect_length(attr(out, "notes"), 2L)
    expect_match(attr(out, "notes")[1], "not yet usable")
    expect_match(attr(out, "notes")[2], "agree to ten significant digits")
    expect_match(attr(out, "notes")[2], "subtract a common value")
  }
  # The exponential mean rule uses the sample mean, not a spread, so the
  # same data are an ordinary case there: sqrt(100 / 1) * 0.3 = 3 <= 10.
  out <- seq_mr_exp_mean(gain, A = 100, c = 1)
  expect_equal(out$value[out$term == "optimal_n"], 3)
  expect_true(attr(out, "is_satisfied"))
})

test_that("a variance or a risk that cannot be represented is an error", {
  x <- c(1.3, 0.2, 2.9, 0.7, 1.1)
  # Squares of data near 1e160 overflow; squares of data near 1e-200
  # underflow to zero although the observations differ.
  expect_error(seq_mr_norm_mean(x * 1e160, A = 1, c = 1e300),
               "sample variance of 'x' is too large to represent")
  expect_error(seq_mr_dfr_mean_arb(x * 1e160, A = 1, c = 1e300),
               "sample variance of 'x' is too large to represent")
  expect_error(seq_mr_norm_mean(x * 1e-200, A = 100, c = 1),
               "sample variance of 'x' is too small to represent")
  expect_error(seq_mr_dfr_mean_arb(x * 1e-200, A = 100, c = 1),
               "sample variance of 'x' is too small to represent")
  expect_error(seq_mr_nexp_mean(x * 6e307, A = 1, c = 1e300, k = 1, t = 1),
               "scale estimate of 'x' is too large to represent")
  # A = c = 1e308: the criterion is the sample standard deviation, about 1,
  # but the risk A s^2 / n + c n overflows. It used to be reported as Inf
  # beside a definite stopping decision. The message shows the arguments.
  expect_lt(stats::sd(x), 2)
  expect_error(seq_mr_norm_mean(x, A = 1e308, c = 1e308),
               "estimated risk is too large to represent \\(Inf\\) with A = 1e\\+308, c = 1e\\+308")
  expect_error(seq_mr_dfr_mean_arb(x, A = 1e308, c = 1e308),
               "estimated risk is too large to represent \\(Inf\\) with A = 1e\\+308, c = 1e\\+308")
  # Exponential mean: sqrt(A / c) * mean = 1.24, so the criterion is 2, but
  # c n = 5e308 overflows.
  expect_error(seq_mr_exp_mean(x, A = 1e308, c = 1e308),
               "estimated risk is too large to represent \\(Inf\\) with A = 1e\\+308, c = 1e\\+308")
  # The message ends by naming the data argument.
  expect_error(seq_mr_norm_mean(x, A = 1e308, c = 1e308),
               "Check these arguments and the units of 'x'\\.$")
})

test_that("no row is Inf or NaN at the edges of the double range", {
  # Every row holds a finite number or a documented NA.
  rows_ok <- function(out) {
    all(is.finite(out$value) | (is.na(out$value) & !is.nan(out$value)))
  }
  spread <- list(
    function(x) seq_mr_norm_mean(x, A = 100, c = 1),
    function(x) seq_mr_nexp_mean(x, A = 100, c = 1, k = 1, t = 1),
    function(x) seq_mr_dfr_mean_arb(x, A = 100, c = 1))
  # Samples with no usable spread whose squares or sums leave the range:
  # constant up to rounding near 2^1000 (the variance overflows) and near
  # 2^-1000 (it underflows), and exactly constant near the largest double
  # (the sum of the data overflows) and at the smallest.
  j <- c(0, 1, 3, 1)
  samples <- list(2^1000 * (1 + j * 2^-40), 2^-1000 * (1 + j * 2^-40),
                  rep(1.7e308, 2), rep(-1.7e308, 4), rep(2^-1074, 3))
  for (check in spread) {
    for (x in samples) {
      out <- check(x)
      expect_identical(attr(out, "is_satisfied"), FALSE)
      expect_identical(attr(out, "variance_unusable"), TRUE)
      expect_true(all(is.na(out$value[out$term %in% c("risk", "optimal_n")])))
      expect_true(rows_ok(out))
    }
  }
  # The first rows are the estimates themselves: the mean of two equal
  # observations and the minimum of four.
  expect_identical(seq_mr_norm_mean(rep(1.7e308, 2), A = 100, c = 1)$value[1],
                   1.7e308)
  expect_identical(seq_mr_nexp_mean(rep(-1.7e308, 4), A = 100, c = 1, k = 1,
                                    t = 1)$value[1], -1.7e308)
  # Three copies of the largest double itself: where mean() cannot form
  # their mean (it returns Inf on a platform without extended precision),
  # the function stops and says to change the units; no row is ever Inf.
  every <- c(spread, list(function(x) seq_mr_exp_mean(x, A = 100, c = 1)))
  for (check in every) {
    top <- tryCatch(check(rep(.Machine$double.xmax, 3)),
                    error = conditionMessage)
    if (is.character(top)) {
      expect_match(top, "too large to represent")
    } else {
      expect_true(rows_ok(top))
    }
  }
})

test_that("the checks refuse invalid arguments and data", {
  set.seed(113)
  x <- rexp(30)
  checks <- list(
    function(x, A = 100, c = 1, ...) seq_mr_norm_mean(x, A = A, c = c, ...),
    function(x, A = 100, c = 1, ...) seq_mr_exp_mean(x, A = A, c = c, ...),
    function(x, A = 100, c = 1, ...) seq_mr_dfr_mean_arb(x, A = A, c = c, ...),
    function(x, A = 100, c = 1, ...) seq_mr_nexp_mean(x, A = A, c = c, k = 1, t = 1, ...))
  for (check in checks) {
    expect_error(check(x, A = NA_real_), "'A' must be")
    expect_error(check(x, A = 0), "'A' must be")
    expect_error(check(x, A = "a"), "'A' must be")
    expect_error(check(x, c = NA_real_), "'c' must be")
    expect_error(check(x, c = -1), "'c' must be")
    expect_error(check(x, c = Inf), "'c' must be")
    expect_error(check("a"), "'x' must be a numeric vector")
    expect_error(check(matrix(x, 5)), "not a matrix")
    expect_error(check(c(x, Inf)), "infinite values")
    expect_error(check(c(x, NA), na.rm = FALSE), "missing values")
    expect_error(check(x, na.rm = "yes"), "'na.rm' must be TRUE or FALSE")
    expect_equal(check(c(x, NA))$value, check(x)$value)
    # A cost so small that the sample size cannot be represented is an error
    # rather than an NA in the stopping flag.
    expect_error(check(x, A = 1e12, c = 1e-12),
                 "implied sample size .* is too large to represent")
  }
  # The advice that closes the too large message names what to check.
  expect_error(seq_mr_norm_mean(x, A = 1e12, c = 1e-12),
               "represent\\. Check the loss constant \\('A', or 'price' and 'epsilon'\\), the cost 'c', and the units of 'x'\\.$")
  expect_error(seq_mr_exp_mean(x, A = 1e12, c = 1e-12),
               "represent\\. Check 'A', the cost 'c', and the units of 'x'\\.$")
  expect_error(seq_mr_dfr_mean_arb(x, A = 1e12, c = 1e-12),
               "represent\\. Check 'A', the cost 'c', and the units of 'x'\\.$")
  expect_error(seq_mr_nexp_mean(x, A = 1e12, c = 1e-12, k = 1, t = 1),
               "represent\\. Check 'A', 'c', 'k', and 't', and the units of 'x'\\.$")
  # Too few observations: the message names the argument and says what to
  # do in a sequential procedure.
  expect_error(seq_mr_norm_mean(1, A = 100, c = 1),
               "'x' needs at least 2 observations; it has 1\\. Collect more data, then check the stopping rule again\\.")
  expect_error(seq_mr_nexp_mean(1, A = 100, c = 1, k = 1, t = 1),
               "'x' needs at least 2 observations; it has 1")
  expect_error(seq_mr_dfr_mean_arb(1, A = 100, c = 1),
               "'x' needs at least 2 observations; it has 1")
  expect_error(seq_mr_exp_mean(numeric(0), A = 100, c = 1),
               "'x' needs at least 1 observation; it has 0")
  expect_error(seq_mr_nexp_mean(x, A = 500, c = 0.01, k = NA_real_, t = 2), "'k' must be")
  expect_error(seq_mr_nexp_mean(x, A = 500, c = 0.01, k = -1, t = 2), "'k' must be")
  expect_error(seq_mr_nexp_mean(x, A = 500, c = 0.01, k = 1, t = 0), "'t' must be")
  # n^t = 30^400 overflows: the message shows the arguments that enter the
  # risk with their values.
  expect_error(seq_mr_nexp_mean(x, A = 500, c = 0.01, k = 1, t = 400),
               "risk is too large to represent \\(Inf\\) with A = 500, c = 0.01, k = 1, t = 400")
  expect_error(seq_mr_dfr_mean_arb(x, A = 100, c = 2, gamma = 0), "'gamma' must be")
  expect_error(seq_mr_dfr_mean_arb(x, A = 100, c = 2, gamma = NA_real_), "'gamma' must be")
})

# ---- operating characteristics ----------------------------------------------
# Each procedure is run to stopping 1,500 times from its own pilot at scales
# 0.25, 1, and 3. Two things are checked against their targets: the mean
# final sample size against n*, and the risk (mean realized loss) against the
# minimum fixed sample size risk. The "expected" risk ratios are the results
# of an independent simulation with 400,000 replications (Monte Carlo SE
# 0.001; an engine that never calls the package), which measure the
# small sample bias of each procedure. The same simulation gives the
# standard deviation of one run's loss ratio, from which the standard error
# of the 1,500 run estimate follows. The loss is strongly right skewed (its
# estimation error part is close to a multiple of a chi square variable
# with one degree of freedom), so the mean of 1,500 losses has a longer
# right tail than a normal variable, and the tolerance is five of those
# standard errors, not four. The seed is fixed, so the tests cannot flake.

mr_oc <- function(reps, m, rgen, ok, n_star, min_risk, loss, check) {
  sim <- seqc_simulate(reps, m, rgen, ok, n_guess = n_star, measure = loss,
                       check = check, n_confirm = 30L, n_full = 2L)
  list(ratio_n = sim$mean_n / n_star, risk_eff = sim$mean_y / min_risk,
       agree = sim$agree)
}

test_that("seq_mr_norm_mean is near n* and the minimum risk at scales 0.25, 1, 3", {
  # A = 100, c = 0.01: n* = 100 sigma and minimum risk 2 c n*.
  A <- 100; cost <- 0.01
  expected_risk <- c(1.016, 1.002, 1.001)   # sd of loss ratio 0.76, 0.71, 0.71
  tol_risk <- 5 * c(0.76, 0.71, 0.71) / sqrt(1500)
  tol_n <- c(0.05, 0.02, 0.02)
  set.seed(113)
  for (i in 1:3) {
    sigma <- c(0.25, 1, 3)[i]
    n_star <- sqrt(A / cost) * sigma
    oc <- mr_oc(
      1500, seq_mr_norm_mean_pilot(A = A, c = cost),
      rgen = function(n) rnorm(n, 2, sigma),
      ok = function(x, n) n >= sqrt(A / cost) * sqrt(seqc_run_mv(x - 2, n)$var),
      n_star = n_star, min_risk = 2 * cost * n_star,
      loss = function(x, N) A * (mean(x[seq_len(N)]) - 2)^2 + cost * N,
      check = function(x) seq_mr_norm_mean(x, A = A, c = cost))
    expect_true(oc$agree)
    expect_lt(abs(oc$ratio_n - 1), tol_n[i])
    expect_lt(abs(oc$risk_eff - expected_risk[i]), tol_risk[i])
  }
})

test_that("seq_mr_exp_mean is near n* and the minimum risk at scales 0.25, 1, 3", {
  # A = 100, c = 0.01: n* = 100 lambda. At lambda = 0.25 (n* = 25) the risk
  # is about 8% above the minimum; the help page says so.
  A <- 100; cost <- 0.01
  expected_risk <- c(1.079, 1.017, 1.004)   # sd of loss ratio 0.84, 0.74, 0.71
  tol_risk <- 5 * c(0.84, 0.74, 0.71) / sqrt(1500)
  tol_n <- c(0.05, 0.02, 0.02)
  set.seed(113)
  for (i in 1:3) {
    lambda <- c(0.25, 1, 3)[i]
    n_star <- sqrt(A / cost) * lambda
    oc <- mr_oc(
      1500, seq_mr_exp_mean_pilot(A = A, c = cost),
      rgen = function(n) rexp(n, rate = 1 / lambda),
      ok = function(x, n) n >= sqrt(A / cost) * cumsum(x)[n] / n,
      n_star = n_star, min_risk = 2 * cost * n_star,
      loss = function(x, N) A * (mean(x[seq_len(N)]) - lambda)^2 + cost * N,
      check = function(x) seq_mr_exp_mean(x, A = A, c = cost))
    expect_true(oc$agree)
    expect_lt(abs(oc$ratio_n - 1), tol_n[i])
    expect_lt(abs(oc$risk_eff - expected_risk[i]), tol_risk[i])
  }
})

test_that("seq_mr_nexp_mean is near n* and the minimum risk at scales 0.25, 1, 3", {
  # Squared error, linear cost (k = 2, t = 1), A = 500, c = 0.01:
  # risk(n) = 2 A (sigma / n)^2 + c n, minimized at n* = (4 A sigma^2 / c)^(1/3).
  A <- 500; cost <- 0.01; theta <- 3
  expected_risk <- c(1.027, 1.006, 1.005)   # sd of loss ratio 0.88, 0.77, 0.77
  tol_risk <- 5 * c(0.88, 0.77, 0.77) / sqrt(1500)
  set.seed(113)
  for (i in 1:3) {
    sigma <- c(0.25, 1, 3)[i]
    n_star <- (4 * A * sigma^2 / cost)^(1 / 3)
    min_risk <- 2 * A * (sigma / n_star)^2 + cost * n_star
    oc <- mr_oc(
      1500, seq_mr_nexp_mean_pilot(A = A, c = cost, k = 2, t = 1),
      rgen = function(n) theta + rexp(n, rate = 1 / sigma),
      ok = function(x, n) n >= (4 * A * seqc_run_u(x, n)^2 / cost)^(1 / 3),
      n_star = n_star, min_risk = min_risk,
      loss = function(x, N) A * (min(x[seq_len(N)]) - theta)^2 + cost * N,
      check = function(x) seq_mr_nexp_mean(x, A = A, c = cost, k = 2, t = 1))
    expect_true(oc$agree)
    expect_lt(abs(oc$ratio_n - 1), 0.02)
    expect_lt(abs(oc$risk_eff - expected_risk[i]), tol_risk[i])
  }
})

test_that("seq_mr_dfr_mean_arb is near n* and the minimum risk at scales 0.25, 1, 3", {
  # Exponential data with standard deviation sigma, A = 100, c = 0.01,
  # gamma = 1 (the default): n* = 100 sigma. At sigma = 0.25 the correction
  # term 1 / n is not small beside the standard deviation and the mean final
  # sample size is about 1.10 n*; the help page says so.
  A <- 100; cost <- 0.01
  expected_n <- c(1.099, 0.988, 0.993)
  tol_n <- c(0.03, 0.02, 0.02)              # 4 SE is 0.022, 0.014, 0.008
  expected_risk <- c(1.019, 1.022, 1.007)   # sd of loss ratio 0.65, 0.76, 0.73
  tol_risk <- 5 * c(0.65, 0.76, 0.73) / sqrt(1500)
  set.seed(113)
  for (i in 1:3) {
    sigma <- c(0.25, 1, 3)[i]
    n_star <- sqrt(A / cost) * sigma
    oc <- mr_oc(
      1500, seq_mr_dfr_mean_arb_pilot(A = A, c = cost),
      rgen = function(n) rexp(n, rate = 1 / sigma),
      ok = function(x, n) n >= sqrt(A / cost) * (sqrt(seqc_run_mv(x, n)$var) + 1 / n),
      n_star = n_star, min_risk = 2 * cost * n_star,
      loss = function(x, N) A * (mean(x[seq_len(N)]) - sigma)^2 + cost * N,
      check = function(x) seq_mr_dfr_mean_arb(x, A = A, c = cost))
    expect_true(oc$agree)
    expect_lt(abs(oc$ratio_n - expected_n[i]), tol_n[i])
    expect_lt(abs(oc$risk_eff - expected_risk[i]), tol_risk[i])
  }
})
