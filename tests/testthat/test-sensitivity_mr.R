# Tests for the minimum risk point estimation sensitivity wrappers,
# seq_mr_cv_sensitivity() and seq_mr_smd_sensitivity(): the benchmark rows
# against closed forms and a numerical minimizer, each wrapper against the
# same procedure run by hand with the exported check functions, published
# operating characteristics from the two articles, and the argument checks.

sens <- function(out, term) out$value[out$term == term]

# The two helpers below run a procedure by hand until its rule is met. A
# defect in the code under test (rows that do not lead to the stopping
# point, a rule that is never met) must make a test fail, not make it run
# forever, so each loop gives up with an error after this many checks. No
# run these helpers make in this file needs more than 32 checks (the
# largest final sample size among them is 51), so the bound is about thirty
# times what a correct procedure uses.
by_hand_max_checks <- 1000
by_hand_give_up <- function(what) {
  stop("The ", what, " procedure run by hand had not stopped after ",
       by_hand_max_checks, " checks of the stopping rule.", call. = FALSE)
}

# The coefficient of variation procedure run by hand with the exported
# functions: `draw(k)` returns k new observations. One row per replication.
cv_by_hand <- function(G, seed, draw, A, cost, batch = 1, pilot_args = list(),
                       check_args = list()) {
  withr::with_seed(seed, {
    vapply(seq_len(G), function(g) {
      m <- do.call(seq_mr_cv_pilot, c(list(A = A, c = cost), pilot_args))
      x <- draw(m)
      checks <- 0
      repeat {
        check <- do.call(seq_mr_cv, c(list(x, A = A, c = cost), check_args))
        if (attr(check, "is_satisfied")) break
        checks <- checks + 1
        if (checks >= by_hand_max_checks) {
          by_hand_give_up("coefficient of variation")
        }
        x <- c(x, draw(batch))
      }
      c(n = length(x), estimate = check$value[check$term == "cv"])
    }, numeric(2))
  })
}

# The standardized mean difference procedure run by hand as ?seq_mr_smd
# describes it, with the exported functions only: while the rule is not
# met, the groups are brought up to the sizes in the rows next_n1 and
# next_n2. No formula of the article is written out here; the rows
# themselves are pinned to hand computed values in test-seq_mr.R.
#
# With batch > 1 and one cost, `batch` observations are added to each group.
# With batch > 1 and two costs, `batch` observations are added to the less
# costly group and the function is called; if the rows then ask only for
# what the costlier group is short of, that is collected and the function
# is called again (the batch recipe in ?seq_mr_smd).
#
# `draw` is a list of two functions; draw[[j]](k) returns k new observations
# for group j. The less costly group is drawn first at every stage, as in
# the wrapper, so the two use the random numbers in the same order.
smd_by_hand <- function(G, seed, draw, A, cost, batch = 1, pilot_args = list(),
                        check_args = list()) {
  check_now <- function(y) {
    do.call(seq_mr_smd, c(list(y[[1]], y[[2]], A = A, c = cost), check_args))
  }
  next_n <- function(check) c(sens(check, "next_n1"), sens(check, "next_n2"))
  grow <- function(y, j, to) c(y[[j]], draw[[j]](to - length(y[[j]])))
  two_costs <- length(cost) == 2L
  lo <- if (two_costs) which.min(cost) else 1L
  hi <- 3L - lo
  withr::with_seed(seed, {
    vapply(seq_len(G), function(g) {
      m <- do.call(seq_mr_smd_pilot, c(list(A = A, c = cost), pilot_args))
      y <- list(draw[[1]](m[[1]]), draw[[2]](m[[length(m)]]))
      check <- check_now(y)
      checks <- 0
      while (!attr(check, "is_satisfied")) {
        checks <- checks + 1
        if (checks >= by_hand_max_checks) {
          by_hand_give_up("standardized mean difference")
        }
        if (batch == 1) {
          target <- next_n(check)
          # Rows that ask for nothing, or for fewer observations than are in
          # hand, cannot lead anywhere.
          if (any(target < lengths(y)) || all(target == lengths(y))) {
            stop("The rule is not met at sizes ", lengths(y)[1], " and ",
                 lengths(y)[2], ", but next_n1 and next_n2 are ", target[1],
                 " and ", target[2], ".", call. = FALSE)
          }
        } else if (two_costs) {
          target <- lengths(y)
          target[lo] <- target[lo] + batch
        } else {
          target <- lengths(y) + batch
        }
        for (j in c(lo, hi)) y[[j]] <- grow(y, j, target[j])
        check <- check_now(y)
        if (batch > 1 && two_costs && !attr(check, "is_satisfied") &&
            next_n(check)[lo] == length(y[[lo]])) {
          y[[hi]] <- grow(y, hi, next_n(check)[hi])
          check <- check_now(y)
        }
      }
      c(n1 = length(y[[1]]), n2 = length(y[[2]]),
        estimate = check$value[check$term == "smd"])
    }, numeric(3))
  })
}

test_that("the wrappers return the documented rows, all finite", {
  cv <- seq_mr_cv_sensitivity(true_value = 0.5, A = 2000, c = 1, G = 20,
                              seed = 113)
  expect_s3_class(cv, c("smsd_sensitivity", "data.frame"), exact = TRUE)
  expect_identical(names(cv), c("term", "value"))
  expect_identical(cv$term,
                   c("true_value", "loss_A", "cost", "replications",
                     "optimal_n", "minimum_risk",
                     "mean_n", "median_n", "sd_n", "ratio_n",
                     "mean_estimate",
                     "mean_risk", "risk_efficiency", "se_risk_efficiency",
                     "regret", "pct_hit_max_n"))
  expect_true(all(is.finite(cv$value)))
  expect_equal(sens(cv, "replications"), 20)

  smd <- seq_mr_smd_sensitivity(true_value = 0.5, A = 500, c = c(2, 4),
                                G = 20, seed = 113)
  expect_s3_class(smd, c("smsd_sensitivity", "data.frame"), exact = TRUE)
  expect_identical(smd$term,
                   c("true_value", "loss_A", "cost_1", "cost_2",
                     "replications",
                     "optimal_n1", "optimal_n2", "minimum_risk",
                     "mean_n1", "mean_n2", "median_n1", "median_n2",
                     "sd_n1", "sd_n2", "ratio_n1", "ratio_n2",
                     "mean_estimate",
                     "mean_risk", "risk_efficiency", "se_risk_efficiency",
                     "regret", "pct_hit_max_n"))
  expect_true(all(is.finite(smd$value)))
  # Internal consistency of the derived rows.
  expect_equal(sens(smd, "ratio_n1"),
               sens(smd, "mean_n1") / sens(smd, "optimal_n1"))
  expect_equal(sens(smd, "risk_efficiency"),
               sens(smd, "mean_risk") / sens(smd, "minimum_risk"))
  expect_equal(sens(smd, "regret"),
               sens(smd, "mean_risk") - sens(smd, "minimum_risk"))
  # Both print through the shared sensitivity class: one line per row under
  # a header, in fixed notation (a plain data frame printed this table as
  # 5.000000e-01, 2.000000e+03, and so on).
  for (out in list(cv, smd)) {
    lines <- capture.output(print(out))
    expect_length(lines, nrow(out) + 1L)
    expect_false(any(grepl("e[+-][0-9]", lines)))
  }
  expect_match(capture.output(print(cv))[3], "^loss_A +2000$")
})

test_that("the wrappers take their arguments in the package wide order", {
  # Design arguments, then G, seed, m_0, gamma, max_n, batch.
  expect_identical(names(formals(seq_mr_cv_sensitivity)),
                   c("true_value", "A", "c", "distribution", "variance",
                     "G", "seed", "m_0", "gamma", "max_n", "batch"))
  expect_identical(names(formals(seq_mr_smd_sensitivity)),
                   c("true_value", "A", "c", "distribution",
                     "G", "seed", "m_0", "gamma", "max_n", "batch"))
  expect_identical(formals(seq_mr_cv_sensitivity)$G, 1000)
  expect_null(formals(seq_mr_smd_sensitivity)$seed)
})

test_that("the coefficient of variation benchmark is the closed form optimum", {
  # n_c = sqrt(A xi^2 / c) and the minimum risk 2 c n_c (Chattopadhyay &
  # Kelley, 2016, Equations 15 and 16), with xi^2 the asymptotic variance of
  # the population. Gamma with kappa^2 = 1 / shape: xi^2 = kappa^2 (1 +
  # kappa^2) / 2 (derived in test-seq_mr.R). Normal: kappa^4 + kappa^2 / 2.
  g <- seq_mr_cv_sensitivity(true_value = 0.4, A = 5e6, c = 10,
                             distribution = "gamma", G = 2, seed = 113,
                             max_n = 82)
  expect_equal(sens(g, "optimal_n"), sqrt(5e5 * 0.16 * 1.16 / 2))  # 215.41
  expect_equal(sens(g, "minimum_risk"), 20 * sqrt(5e5 * 0.16 * 1.16 / 2))
  n <- seq_mr_cv_sensitivity(true_value = 0.2, A = 5e6, c = 10,
                             distribution = "normal", G = 2, seed = 113,
                             max_n = 82)
  expect_equal(sens(n, "optimal_n"), sqrt(5e5 * (0.2^4 + 0.2^2 / 2)))
  # Lognormal: the population moments by numerical integration, put into
  # Kelley, Darku, and Chattopadhyay (2018), Equation 30.
  kappa <- 0.8
  sdlog <- sqrt(log(1 + kappa^2))
  moment <- function(f) {
    # X = exp(sdlog * Z); the standard normal mass beyond 12 is negligible.
    stats::integrate(function(z) f(exp(sdlog * z)) * stats::dnorm(z),
                     -12, 12, rel.tol = 1e-12, subdivisions = 1000L)$value
  }
  mu <- moment(identity)
  s2 <- moment(function(x) (x - mu)^2)
  mu3 <- moment(function(x) (x - mu)^3)
  mu4 <- moment(function(x) (x - mu)^4)
  xi2 <- mu4 / (4 * mu^2 * s2) - s2 / (4 * mu^2) - mu3 / mu^3 + s2^2 / mu^4
  l <- seq_mr_cv_sensitivity(true_value = kappa, A = 5e6, c = 10,
                             distribution = "lognormal", G = 2, seed = 113,
                             max_n = 82)
  expect_equal(sens(l, "optimal_n"), sqrt(5e5 * xi2), tolerance = 1e-6)
  # The benchmark does not depend on the variance estimator requested.
  p <- seq_mr_cv_sensitivity(true_value = 0.4, A = 5e6, c = 10,
                             variance = "published", G = 2, seed = 113,
                             max_n = 82)
  expect_identical(sens(p, "optimal_n"), sens(g, "optimal_n"))
})

test_that("the standardized mean difference benchmark minimizes the risk", {
  # One cost: n_c = sqrt(A / (2c)) sqrt(2 + delta^2 / 4), minimum risk
  # 4 c n_c (Chattopadhyay & Kelley, 2017, Equations 10 to 12).
  e <- seq_mr_smd_sensitivity(true_value = 0.8, A = 1e6, c = 500, G = 2,
                              seed = 113, max_n = 11)
  n_c <- sqrt(1e6 / 1000) * sqrt(2 + 0.8^2 / 4)             # 46.476
  expect_equal(sens(e, "optimal_n1"), n_c)
  expect_equal(sens(e, "optimal_n2"), n_c)
  expect_equal(sens(e, "minimum_risk"), 4 * 500 * n_c)
  # Two costs: a direct numerical minimizer of Equation 17.
  risk <- function(p, A, c1, c2, delta) {
    n <- exp(p)
    A * (1 / n[1] + 1 / n[2] + delta^2 / (2 * sum(n))) + c1 * n[1] + c2 * n[2]
  }
  fit <- stats::optim(log(c(150, 80)), risk, A = 50000, c1 = 8, c2 = 2,
                      delta = 2, method = "BFGS",
                      control = list(reltol = 1e-14, maxit = 2000))
  u <- seq_mr_smd_sensitivity(true_value = 2, A = 50000, c = c(8, 2), G = 2,
                              seed = 113, max_n = 30)
  expect_equal(c(sens(u, "optimal_n1"), sens(u, "optimal_n2")),
               exp(fit$par), tolerance = 1e-5)
  expect_equal(sens(u, "minimum_risk"), fit$value, tolerance = 1e-9)
  # Skewed population (standardized chi square on 3 df, kurtosis 7): the
  # large sample MSE of d is 1/n1 + 1/n2 + delta^2 (7 - 1) / (4 (n1 + n2)),
  # so with one cost the optimum is sqrt(A / (2c)) sqrt(2 + 3 delta^2 / 4).
  # (At n1 = n2 = 200 and delta = 1 the formula gives 0.01375 and normal
  # theory 0.01125; 20,000 simulated pairs of samples, drawn as
  # (rchisq(200, 3) - 3) / sqrt(6) after set.seed(1), give 0.01375 with a
  # standard error of 0.00014.)
  k <- seq_mr_smd_sensitivity(true_value = 1, A = 5000, c = 1,
                              distribution = "skewed", G = 2, seed = 113,
                              max_n = 15)
  expect_equal(sens(k, "optimal_n1"), sqrt(2500) * sqrt(2 + 3 / 4))
})

test_that("seq_mr_cv_sensitivity is the procedure run by hand", {
  A <- 3000; cost <- 2; kappa <- 0.5; G <- 12
  out <- seq_mr_cv_sensitivity(true_value = kappa, A = A, c = cost, G = G,
                               seed = 113)
  by_hand <- withr::with_seed(113, {
    vapply(seq_len(G), function(g) {
      x <- rgamma(seq_mr_cv_pilot(A, cost), shape = 1 / kappa^2, scale = 1)
      repeat {
        check <- seq_mr_cv(x, A = A, c = cost)
        if (attr(check, "is_satisfied")) break
        x <- c(x, rgamma(1, shape = 1 / kappa^2, scale = 1))
      }
      c(n = length(x), cv = check$value[check$term == "cv"])
    }, numeric(2))
  })
  expect_equal(sens(out, "mean_n"), mean(by_hand["n", ]))
  expect_equal(sens(out, "sd_n"), sd(by_hand["n", ]))
  expect_equal(sens(out, "mean_estimate"), mean(by_hand["cv", ]))
  expect_equal(sens(out, "mean_risk"),
               mean(A * (by_hand["cv", ] - kappa)^2 + cost * by_hand["n", ]))
})

test_that("seq_mr_smd_sensitivity is the procedure run by hand, in the article's order", {
  A <- 2000; delta <- 0.5; G <- 10
  # One cost: a pair of observations per stage.
  out <- seq_mr_smd_sensitivity(true_value = delta, A = A, c = 5, G = G,
                                seed = 113)
  by_hand <- withr::with_seed(113, {
    vapply(seq_len(G), function(g) {
      m <- seq_mr_smd_pilot(A, 5)
      g1 <- rnorm(m, delta); g2 <- rnorm(m)
      repeat {
        check <- seq_mr_smd(g1, g2, A = A, c = 5)
        if (attr(check, "is_satisfied")) break
        g1 <- c(g1, rnorm(1, delta)); g2 <- c(g2, rnorm(1))
      }
      c(n = length(g1), d = check$value[check$term == "smd"])
    }, numeric(2))
  })
  expect_equal(sens(out, "mean_n1"), mean(by_hand["n", ]))
  expect_equal(sens(out, "mean_n2"), mean(by_hand["n", ]))
  expect_equal(sens(out, "mean_estimate"), mean(by_hand["d", ]))

  # Two costs, group 2 the less costly: one observation to group 2, then
  # group 1 brought up to ceiling(n / sqrt(1 + (c_hi - c_lo) n^2 / A)). The
  # article's formula is written out here, independently of the next_n rows
  # the wrapper follows.
  cost <- c(8, 3)
  out <- seq_mr_smd_sensitivity(true_value = delta, A = A, c = cost, G = G,
                                seed = 113)
  by_hand <- withr::with_seed(113, {
    vapply(seq_len(G), function(g) {
      m <- seq_mr_smd_pilot(A, cost)
      g1 <- rnorm(m[["group_1"]], delta); g2 <- rnorm(m[["group_2"]])
      repeat {
        check <- seq_mr_smd(g1, g2, A = A, c = cost)
        if (attr(check, "is_satisfied")) break
        g2 <- c(g2, rnorm(1))
        paired <- ceiling(length(g2) / sqrt(1 + 5 * length(g2)^2 / A))
        g1 <- c(g1, rnorm(max(0, paired - length(g1)), delta))
      }
      c(n1 = length(g1), n2 = length(g2),
        d = check$value[check$term == "smd"])
    }, numeric(3))
  })
  expect_equal(sens(out, "mean_n1"), mean(by_hand["n1", ]))
  expect_equal(sens(out, "mean_n2"), mean(by_hand["n2", ]))
  expect_equal(sens(out, "mean_risk"),
               mean(A * (by_hand["d", ] - delta)^2 + 8 * by_hand["n1", ] +
                      3 * by_hand["n2", ]))
  expect_gt(sens(out, "mean_n2"), sens(out, "mean_n1"))   # cheaper: more
})

test_that("the equal cost procedure reproduces a row of the article's Table 3", {
  # Chattopadhyay and Kelley (2017), Table 3 (p. 104), delta = 0,
  # A = 1,000,000, c = 500, one observation per group per stage: mean final
  # sample size 50.0210 (standard error .0020 over 5,000 replications, so a
  # standard deviation of about 0.14). With G = 400 the standard error of
  # the mean here is 0.14 / 20 = 0.007; the tolerance of 0.04 is more than
  # five standard errors. The package's earlier rule gave 47.02.
  out <- seq_mr_smd_sensitivity(true_value = 0, A = 1e6, c = 500, G = 400,
                                seed = 113)
  expect_lt(abs(sens(out, "mean_n1") - 50.0210), 0.04)
  expect_equal(sens(out, "pct_hit_max_n"), 0)
  # Risk efficiency near 1 (the article reports 1.0023 for this row). Its
  # Monte Carlo standard error, reported by the function, is about 0.03
  # here; the tolerance is four of them. The reported standard error is
  # itself an estimate: over 60 other seeds its mean was 0.0318 and its
  # standard deviation 0.0031 (largest 0.0404), so the bound of 0.05, which
  # keeps the tolerance from being vacuous, is more than five standard
  # deviations above the mean.
  expect_lt(sens(out, "se_risk_efficiency"), 0.05)
  expect_lt(abs(sens(out, "risk_efficiency") - 1),
            4 * sens(out, "se_risk_efficiency"))
})

test_that("the two cost procedure attains the minimum risk and favors the cheaper group", {
  # A = 5,000, costs 1 and 4, delta = 0.5: the optimum is (72.67, 35.59)
  # (direct minimization of Equation 17).
  out <- seq_mr_smd_sensitivity(true_value = 0.5, A = 5000, c = c(1, 4),
                                G = 300, seed = 113)
  expect_equal(sens(out, "optimal_n1"), 72.675, tolerance = 1e-4)
  expect_equal(sens(out, "optimal_n2"), 35.593, tolerance = 1e-4)
  # Where the article's rule (with the squared term) stops, derived here and
  # not taken from the function. With r(n) = 1 / sqrt(1 + 3 n^2 / 5000), the
  # rule n >= sqrt(5000) [sqrt(1 + d^2 / (2 (1 + r(n))^2)) + n^-0.49] first
  # holds at n = 80.99 when d^2 is held at delta^2 = 0.25. The estimate d
  # varies around delta, so E(d^2) = delta^2 + MSE(d) = 0.25 + 0.041 (the
  # large sample mean squared error at these sizes), which moves the
  # crossing to 81.31; the stopping time is a whole number, which adds half
  # an observation on average. The mean of N for the less costly group
  # should therefore be near 81.8. Each replication's N has a standard
  # deviation of about 1.7, so the standard error of the mean with G = 300
  # is about 0.10; the tolerance of 0.5 is five of them. (Over 60 other
  # seeds the mean of N averaged 81.83, with a standard deviation of 0.085.
  # The same derivation with Equation 18 as printed, the first power, gives
  # 82.8, and with a correction term of sqrt(A / (2c)) n^-gamma it gives
  # 79.5; both are outside the tolerance.)
  r <- function(n) 1 / sqrt(1 + 3 * n^2 / 5000)
  crossing <- stats::uniroot(function(n) {
    n - sqrt(5000) * (sqrt(1 + (0.25 + 0.041) / (2 * (1 + r(n))^2)) + n^-0.49)
  }, c(40, 200), tol = 1e-9)$root
  expect_equal(crossing, 81.31, tolerance = 1e-3)
  expect_lt(abs(sens(out, "mean_n1") - (crossing + 0.5)), 0.5)
  # The costlier group has the size paired with N: n r(n) lies between 36
  # and 37 for every n from 77 to 87 (76.4 r(76.4) = 36.0 and 87.9 r(87.9)
  # = 37.0), so it rounds up to 37 in all but a few replications. Over 60
  # other seeds its mean was 37.006 with a standard deviation of 0.0045, so
  # 0.05 is more than nine standard errors.
  expect_gt(77 * r(77), 36)
  expect_lt(87 * r(87), 37)
  expect_lt(abs(sens(out, "mean_n2") - 37), 0.05)
  # The reported standard error is itself an estimate: over 60 other seeds
  # its mean was 0.0381 and its standard deviation 0.0049 (largest 0.0565),
  # so the bound of 0.065, which keeps the next tolerance from being
  # vacuous, is more than five standard deviations above the mean.
  expect_lt(sens(out, "se_risk_efficiency"), 0.065)
  expect_lt(abs(sens(out, "risk_efficiency") - 1),
            4 * sens(out, "se_risk_efficiency"))
})

test_that("the published variance reproduces a gamma row of the 2016 article's Table 1", {
  # Chattopadhyay and Kelley (2016), Table 1 (p. 642), gamma population
  # with kappa = 0.4, A = 5,000,000, c = 10: mean final sample size
  # 182.3674 (standard error 0.1598 over 5,000 replications, so a standard
  # deviation of about 11.3). With G = 300 the standard error of the mean
  # is about 0.65; the tolerance of 3 is more than four standard errors.
  # For this population the two variance estimators differ: the default one
  # gives a mean near 213 (the fixed sample size that minimizes the risk is
  # 215.4), far outside the tolerance, so this test fails if the wrapper
  # does not pass `variance` on.
  skip_on_cran()
  out <- seq_mr_cv_sensitivity(true_value = 0.4, A = 5e6, c = 10,
                               distribution = "gamma",
                               variance = "published", G = 300, seed = 113)
  expect_lt(abs(sens(out, "mean_n") - 182.3674), 3)
  expect_equal(sens(out, "optimal_n"), sqrt(5e5 * 0.0928))    # 215.41
})

test_that("the published variance reproduces a row of the 2016 article's Table 1", {
  # Chattopadhyay and Kelley (2016), Table 1 (p. 642), Normal(10, 2)
  # (kappa = 0.2), A = 5,000,000, c = 10: mean final sample size 123.8528
  # (standard error 0.0870 over 5,000 replications, so a standard deviation
  # of about 6.2). With G = 300 the standard error of the mean is about
  # 0.36; the tolerance of 1.5 is more than four standard errors.
  skip_on_cran()
  out <- seq_mr_cv_sensitivity(true_value = 0.2, A = 5e6, c = 10,
                               distribution = "normal",
                               variance = "published", G = 300, seed = 113)
  expect_lt(abs(sens(out, "mean_n") - 123.8528), 1.5)
})

test_that("the default coefficient of variation procedure is risk efficient for gamma data", {
  # Gamma, kappa = 0.4, A = 500,000, c = 10: n_c = sqrt(50000 * 0.0928)
  # = 68.12. With G = 5000 and seed = 1 the same call gives a mean final
  # sample size of 68.52 (standard deviation 11.1; ratio 1.006) and a risk
  # efficiency of 0.987 (standard error 0.0098). With G = 400 the standard
  # error of ratio_n is 11.1 / 20 / 68.12 = 0.008 and that of the risk
  # efficiency about 0.034; each tolerance is at least four standard
  # errors. With the published variance the ratio is 0.916, outside the
  # tolerance.
  out <- seq_mr_cv_sensitivity(true_value = 0.4, A = 5e5, c = 10, G = 400,
                               seed = 113)
  expect_equal(sens(out, "optimal_n"), 68.118, tolerance = 1e-4)
  expect_lt(abs(sens(out, "ratio_n") - 1.006), 0.04)
  expect_lt(sens(out, "se_risk_efficiency"), 0.05)
  expect_lt(abs(sens(out, "risk_efficiency") - 1),
            4 * sens(out, "se_risk_efficiency"))
  expect_equal(sens(out, "pct_hit_max_n"), 0)
})

test_that("a small coefficient of variation does not make the procedure sample past the pilot", {
  # Normal data with kappa = 0.02, A = 2,000, c = 1. The pilot size is
  # 2000^(1 / 2.98) = 12.8, so 13. The asymptotic variance,
  # kappa^2 / 2 + kappa^4 = 0.0002, is below 13^-3 = 0.00046, so the
  # estimate in use at the pilot is the bound, and Equation 23 of
  # Chattopadhyay and Kelley (2016) holds there:
  # 2000 (13^-3 + 13^-0.98) = 2000 * 0.081427 = 162.9 <= 169. The rule as
  # printed therefore stops at the pilot whenever the raw estimate is
  # positive and below 169 / 2000 - 13^-0.98 = 0.0035 (17 times its true
  # value). The package's rule does the same: an estimate below the bound is
  # small, not unusable, so it must not keep sampling past the pilot.
  expect_identical(seq_mr_cv_pilot(A = 2000, c = 1), 13L)
  expect_lt(2000 * (13^-3 + 13^-0.98), 13^2)
  # With G = 2000 and seed = 1 the same call stops at 13 in every
  # replication (mean_n 13, sd_n 0), so the chance that a replication goes
  # past the pilot is below 0.0015 (the 95% upper bound for 0 of 2,000).
  # Among 300 replications fewer than one is expected to, and the tolerance
  # of 0.1 on the mean allows 30 extra observations in all.
  out <- seq_mr_cv_sensitivity(true_value = 0.02, A = 2000, c = 1,
                               distribution = "normal", G = 300, seed = 113)
  expect_lt(abs(sens(out, "mean_n") - 13), 0.1)
  expect_equal(sens(out, "pct_hit_max_n"), 0)
  # The estimate at that size: for 13 normal observations the sample
  # coefficient of variation has mean kappa (1 - 1 / (4 n)) = 0.0196 to first
  # order and standard deviation about kappa / sqrt(2 n) = 0.0039 (0.01961
  # and 0.0041 in 20,000 samples of rnorm(13, 1, 0.02) after set.seed(1)),
  # so the mean of 300 has a standard error of 0.00024; the tolerance is
  # four of them.
  expect_lt(abs(sens(out, "mean_estimate") - 0.02 * (1 - 1 / 52)), 0.001)
})

test_that("seq_mr_cv_sensitivity passes every argument on to the procedure", {
  # Each call differs from the defaults in the arguments named, and is
  # compared with the procedure run by hand under the same seed. If the
  # wrapper ignored any of them the two would differ.
  A <- 3000; cost <- 2; kappa <- 0.6; G <- 12
  lognormal <- function(k) rlnorm(k, meanlog = 0, sdlog = sqrt(log(1 + kappa^2)))
  compare <- function(out, by_hand) {
    n <- by_hand["n", ]
    est <- by_hand["estimate", ]
    loss <- A * (est - kappa)^2 + cost * n
    expect_equal(sens(out, "mean_n"), mean(n))
    expect_equal(sens(out, "median_n"), stats::median(n))
    expect_equal(sens(out, "sd_n"), stats::sd(n))
    expect_equal(sens(out, "mean_estimate"), mean(est))
    expect_equal(sens(out, "mean_risk"), mean(loss))
    expect_equal(sens(out, "se_risk_efficiency"),
                 stats::sd(loss) / sqrt(G) / sens(out, "minimum_risk"))
  }
  # distribution, variance, gamma, and batch. The pilot is
  # 1500^(1 / 2.4) = 21.05, so 22 (12 with the default gamma).
  out <- seq_mr_cv_sensitivity(kappa, A = A, c = cost,
                               distribution = "lognormal",
                               variance = "published", G = G, seed = 113,
                               gamma = 0.2, batch = 3)
  by_hand <- cv_by_hand(G, 113, lognormal, A, cost, batch = 3,
                        pilot_args = list(gamma = 0.2),
                        check_args = list(gamma = 0.2,
                                          variance = "published"))
  expect_true(all((by_hand["n", ] - 22) %% 3 == 0))
  compare(out, by_hand)
  # The same data with the default variance stop elsewhere, so the
  # comparison above can tell the two estimators apart.
  other <- cv_by_hand(G, 113, lognormal, A, cost, batch = 3,
                      pilot_args = list(gamma = 0.2),
                      check_args = list(gamma = 0.2))
  expect_false(isTRUE(all.equal(mean(other["n", ]), mean(by_hand["n", ]))))
  # m_0 above the pilot the constants imply (12), with normal data.
  normal <- function(k) rnorm(k, mean = 1, sd = kappa)
  out <- seq_mr_cv_sensitivity(kappa, A = A, c = cost,
                               distribution = "normal", G = G, seed = 113,
                               m_0 = 30)
  by_hand <- cv_by_hand(G, 113, normal, A, cost,
                        pilot_args = list(m_0 = 30))
  expect_true(all(by_hand["n", ] >= 30))
  compare(out, by_hand)
})

test_that("seq_mr_smd_sensitivity passes every argument on to the procedure", {
  # As for the coefficient of variation: each call is compared with the
  # procedure run by hand under the same seed. A true difference of 1.5
  # makes the final sample sizes vary from replication to replication.
  A <- 2000; delta <- 1.5; G <- 10
  skewed_1 <- function(k) (rchisq(k, df = 3) - 3) / sqrt(6) + delta
  skewed_2 <- function(k) (rchisq(k, df = 3) - 3) / sqrt(6)
  normal_1 <- function(k) rnorm(k, mean = delta)
  normal_2 <- function(k) rnorm(k)
  compare <- function(out, by_hand, cost) {
    n1 <- by_hand["n1", ]
    n2 <- by_hand["n2", ]
    est <- by_hand["estimate", ]
    loss <- A * (est - delta)^2 + cost[1] * n1 + cost[length(cost)] * n2
    expect_equal(sens(out, "mean_n1"), mean(n1))
    expect_equal(sens(out, "mean_n2"), mean(n2))
    expect_equal(sens(out, "median_n1"), stats::median(n1))
    expect_equal(sens(out, "median_n2"), stats::median(n2))
    expect_equal(sens(out, "sd_n1"), stats::sd(n1))
    expect_equal(sens(out, "sd_n2"), stats::sd(n2))
    expect_equal(sens(out, "mean_estimate"), mean(est))
    expect_equal(sens(out, "mean_risk"), mean(loss))
    expect_equal(sens(out, "se_risk_efficiency"),
                 stats::sd(loss) / sqrt(G) / sens(out, "minimum_risk"))
  }
  # One cost: the skewed population, gamma, and batch. The pilot is
  # 200^(1 / 2.4) = 9.09, so 10 (6 with the default gamma).
  out <- seq_mr_smd_sensitivity(delta, A = A, c = 5, distribution = "skewed",
                                G = G, seed = 113, gamma = 0.2, batch = 5)
  by_hand <- smd_by_hand(G, 113, list(skewed_1, skewed_2), A, 5, batch = 5,
                         pilot_args = list(gamma = 0.2),
                         check_args = list(gamma = 0.2))
  expect_true(all((by_hand["n1", ] - 10) %% 5 == 0))
  expect_gt(stats::sd(by_hand["n1", ]), 0)
  compare(out, by_hand, 5)
  # One cost: m_0 above the pilot the constants imply
  # (200^(1 / 2.98) = 5.92, so 6).
  out <- seq_mr_smd_sensitivity(delta, A = A, c = 5, G = G, seed = 113,
                                m_0 = 15)
  by_hand <- smd_by_hand(G, 113, list(normal_1, normal_2), A, 5,
                         pilot_args = list(m_0 = 15))
  expect_true(all(by_hand["n1", ] >= 15))
  expect_gt(stats::sd(by_hand["n1", ]), 0)
  compare(out, by_hand, 5)
  # Two costs, the first group the less costly one: gamma and batch. The
  # pilots are (2000 / 3)^(1 / 2.4) = 15.02, so 16, and the paired size
  # 12.49, so 13.
  cost <- c(3, 8)
  out <- seq_mr_smd_sensitivity(delta, A = A, c = cost, G = G, seed = 113,
                                gamma = 0.2, batch = 5)
  by_hand <- smd_by_hand(G, 113, list(normal_1, normal_2), A, cost,
                         batch = 5, pilot_args = list(gamma = 0.2),
                         check_args = list(gamma = 0.2))
  expect_gt(stats::sd(by_hand["n1", ]), 0)
  compare(out, by_hand, cost)
  # Two costs, one observation per stage, in either order of the groups:
  # the wrapper is the loop of ?seq_mr_smd, which follows the next_n rows.
  for (cost in list(c(3, 8), c(8, 3))) {
    out <- seq_mr_smd_sensitivity(delta, A = A, c = cost, G = G, seed = 113)
    by_hand <- smd_by_hand(G, 113, list(normal_1, normal_2), A, cost)
    expect_gt(stats::sd(by_hand["n1", ]) + stats::sd(by_hand["n2", ]), 0)
    compare(out, by_hand, cost)
  }
  # Two costs, the second group the less costly one: m_0 above both pilots
  # (9 and 9 with the default m_0: (2000 / 3)^(1 / 2.98) = 8.86, paired size
  # 8.21), which holds the costlier group at 20.
  cost <- c(8, 3)
  out <- seq_mr_smd_sensitivity(delta, A = A, c = cost, G = G, seed = 113,
                                m_0 = 20, batch = 2)
  by_hand <- smd_by_hand(G, 113, list(normal_1, normal_2), A, cost,
                         batch = 2, pilot_args = list(m_0 = 20))
  expect_true(all(by_hand["n1", ] >= 20))
  expect_gt(stats::sd(by_hand["n2", ]), 0)
  compare(out, by_hand, cost)
})

test_that("with a skewed population the normal theory rule stops short of the optimum", {
  # A = 5,000, c = 1, delta = 2, standardized chi square (3 df) errors. The
  # fixed size that minimizes the risk for this population is
  # 50 sqrt(2 + 3 delta^2 / 4) = 111.80 per group. The rule uses the normal
  # theory variance, 2 + d^2 / 4: n >= 50 (sqrt(2 + d^2 / 4) + n^-0.49)
  # first holds at n = 92.05 for d = delta, and at 92.25 when d^2 is at its
  # expectation, delta^2 + MSE(d) = 4 + 0.054; the whole number stopping
  # time adds half an observation on average. So the ratio of the mean
  # final size to the optimum should be near 92.75 / 111.80 = 0.83. N has a
  # standard deviation of about 3.1, so with G = 150 the standard error of
  # the ratio is about 0.0023; the tolerance of 0.02 is more than eight of
  # them. The estimate itself has a standard deviation of about 0.23, so
  # the standard error of its mean is about 0.019 and 0.1 is five of them.
  crossing <- stats::uniroot(function(n) {
    n - 50 * (sqrt(2 + (4 + 0.054) / 4) + n^-0.49)
  }, c(40, 200), tol = 1e-9)$root
  expect_equal(crossing, 92.25, tolerance = 1e-3)
  out <- seq_mr_smd_sensitivity(true_value = 2, A = 5000, c = 1,
                                distribution = "skewed", G = 150, seed = 113)
  expect_equal(sens(out, "optimal_n1"), 50 * sqrt(5))
  expect_lt(abs(sens(out, "ratio_n1") - (crossing + 0.5) / (50 * sqrt(5))),
            0.02)
  expect_lt(abs(sens(out, "mean_estimate") - 2), 0.1)
  expect_equal(sens(out, "pct_hit_max_n"), 0)
})

test_that("a true difference that is tiny but not zero has the zero difference benchmark", {
  # With delta = 0 the risk separates and each group's optimum is
  # sqrt(A / c): sqrt(250) and sqrt(125), with minimum risk
  # 2 sqrt(500 * 2) + 2 sqrt(500 * 4). A difference that is tiny but not
  # zero (down to 1e-12) must give that benchmark too, and not an error.
  zero <- c(sqrt(250), sqrt(125), 2 * sqrt(1000) + 2 * sqrt(2000))
  for (tiny in c(0, 1e-12, 1e-9, 2e-8, -1e-9)) {
    out <- seq_mr_smd_sensitivity(tiny, A = 500, c = c(2, 4), G = 2,
                                  seed = 113, max_n = 20)
    expect_equal(c(sens(out, "optimal_n1"), sens(out, "optimal_n2"),
                   sens(out, "minimum_risk")), zero, tolerance = 1e-12)
  }
  for (delta2 in 10^seq(-20, -10, by = 0.25)) {
    expect_equal(.smd_mr_optimum(delta2, 500, 4, 2), zero[c(2, 1, 3)],
                 tolerance = 1e-9)
  }
})

test_that("a replication with no estimate is left out of the risk rows", {
  rows <- .mrpe_risk_rows(estimate = c(0.5, NA, 0.7), loss = c(10, NA, 14),
                          minimum_risk = 8)
  expect_equal(rows[["mean_estimate"]], 0.6)
  expect_equal(rows[["mean_risk"]], 12)
  expect_equal(rows[["risk_efficiency"]], 1.5)
  expect_equal(rows[["se_risk_efficiency"]], sqrt(8) / sqrt(2) / 8)
  expect_equal(rows[["regret"]], 4)
  none <- .mrpe_risk_rows(estimate = c(NA, NA), loss = c(NA_real_, NA),
                          minimum_risk = 8)
  expect_true(all(is.na(none)))
  # End to end: a normal population with a coefficient of variation of 3,
  # stopped at the pilot size of 22 by max_n. A sample mean of 22 such
  # observations is negative with probability pnorm(-sqrt(22) / 3) = 0.06,
  # and then there is no estimate. By hand, under the same seed:
  G <- 60
  by_hand <- withr::with_seed(113, {
    vapply(seq_len(G), function(g) {
      x <- rnorm(22, mean = 1, sd = 3)
      if (mean(x) > 0) sd(x) / mean(x) else NA_real_
    }, numeric(1))
  })
  expect_gt(sum(is.na(by_hand)), 0)
  out <- seq_mr_cv_sensitivity(3, A = 1e4, c = 1, distribution = "normal",
                               G = G, seed = 113, max_n = 22)
  expect_equal(sens(out, "mean_n"), 22)
  expect_equal(sens(out, "pct_hit_max_n"), 100)
  expect_equal(sens(out, "mean_estimate"), mean(by_hand, na.rm = TRUE))
  expect_equal(sens(out, "mean_risk"),
               mean(1e4 * (by_hand - 3)^2 + 22, na.rm = TRUE))
})

test_that("a criterion too large to represent counts as not met and does not stop the run", {
  # A sample mean barely above zero gives a huge coefficient of variation
  # and a criterion beyond the integer range. In the check function that is
  # an error; inside a simulation it only means "keep sampling". The first
  # four values have mean 0.00025 and SD 1.83.
  expect_error(seq_mr_cv(c(-2, -1, 1, 2.001), A = 1e6, c = 1),
               "too large to represent")
  stream <- c(-2, -1, 1, 2.001, 5, 6, 4, 5.5, 6.5, 4.5, 5, 6)
  position <- 0
  draw <- function(k) {
    out <- stream[position + seq_len(k)]
    position <<- position + k
    out
  }
  run <- .mr_cv_run_one(draw, pilot = 4, batch = 4, max_n = 12, A = 1e6,
                        c = 1, gamma = 0.49, variance = "corrected")
  expect_equal(run[["n"]], 12)
  expect_equal(run[["estimate"]], sd(stream) / mean(stream))
  expect_equal(run[["capped"]], 1)
})

test_that("a stage that would add no observation is an error, not an endless loop", {
  # A replication of the standardized mean difference procedure brings the
  # groups up to the sizes .smd_mr_next() gives until the rule is met. Those
  # sizes always include at least one new observation. If a defect ever made
  # them equal to the current sizes the loop would never end (and neither
  # would this test file), so the engine checks and stops with an error.
  # Here .smd_mr_next() is replaced by a version that returns the current
  # sizes on its first call and the right sizes after that: with the check
  # the first stage is an error, and without it the run would complete and
  # this test would fail. The rule cannot be met at the pilot (7 per group
  # at most, against a criterion of at least sqrt(500 / 2) = 15.8), so the
  # next sizes are always asked for.
  # with_mocked_bindings() is stable from testthat 3.2.0.
  skip_if_not_installed("testthat", "3.2.0")
  real_next <- .smd_mr_next
  calls <- 0
  stuck_once <- function(n1, n2, A, c1, c2, step = 1) {
    calls <<- calls + 1
    if (calls == 1) c(n1, n2) else real_next(n1, n2, A, c1, c2, step = step)
  }
  for (cost in list(2, c(2, 4))) {
    calls <- 0
    with_mocked_bindings(
      expect_error(
        seq_mr_smd_sensitivity(true_value = 0.5, A = 500, c = cost, G = 2,
                               seed = 113),
        "adds no observation", fixed = TRUE),
      .smd_mr_next = stuck_once
    )
    expect_equal(calls, 1)
  }
  # With the real function back, the same call runs.
  expect_s3_class(seq_mr_smd_sensitivity(true_value = 0.5, A = 500,
                                         c = c(2, 4), G = 2, seed = 113),
                  "smsd_sensitivity")
})

test_that("a seed makes the result reproducible and leaves the caller's stream alone", {
  a <- seq_mr_smd_sensitivity(true_value = 0.5, A = 500, c = 2, G = 10,
                              seed = 113)
  b <- seq_mr_smd_sensitivity(true_value = 0.5, A = 500, c = 2, G = 10,
                              seed = 113)
  expect_identical(a, b)
  other <- seq_mr_smd_sensitivity(true_value = 0.5, A = 500, c = 2, G = 10,
                                  seed = 114)
  expect_false(identical(a$value, other$value))
  # The caller's stream is restored after a seeded call. (The seeds 7, 21,
  # and 22 below only label distinct caller streams; the procedure itself is
  # seeded with 113.)
  with_call <- withr::with_seed(7, {
    first <- runif(1)
    seq_mr_cv_sensitivity(true_value = 0.5, A = 2000, c = 1, G = 5,
                          seed = 113)
    c(first, runif(1))
  })
  expect_identical(with_call, withr::with_seed(7, runif(2)))
  # seed = NULL draws from the caller's stream.
  u1 <- withr::with_seed(21, seq_mr_cv_sensitivity(0.5, A = 2000, c = 1,
                                                   G = 5))
  u2 <- withr::with_seed(21, seq_mr_cv_sensitivity(0.5, A = 2000, c = 1,
                                                   G = 5))
  u3 <- withr::with_seed(22, seq_mr_cv_sensitivity(0.5, A = 2000, c = 1,
                                                   G = 5))
  expect_identical(u1, u2)
  expect_false(identical(u1$value, u3$value))
})

test_that("batch and max_n do what the help says", {
  # Every final sample size is the pilot plus a whole number of batches.
  pilot <- seq_mr_cv_pilot(A = 2000, c = 1)                  # 13
  out <- seq_mr_cv_sensitivity(true_value = 0.5, A = 2000, c = 1, G = 11,
                               seed = 113, batch = 4)
  expect_equal((sens(out, "median_n") - pilot) %% 4, 0)
  # A larger batch cannot stop earlier on average than checking every time.
  one <- seq_mr_cv_sensitivity(true_value = 0.5, A = 2000, c = 1, G = 200,
                               seed = 113)
  ten <- seq_mr_cv_sensitivity(true_value = 0.5, A = 2000, c = 1, G = 200,
                               seed = 113, batch = 10)
  expect_gt(sens(ten, "mean_n"), sens(one, "mean_n"))
  # max_n is an exact cap, also when it is not a whole number of batches
  # past the pilot. A = 5e6 needs well over 85 observations.
  capped <- seq_mr_cv_sensitivity(true_value = 0.4, A = 5e6, c = 10, G = 5,
                                  seed = 113, batch = 2, max_n = 85)
  expect_equal(sens(capped, "mean_n"), 85)
  expect_equal(sens(capped, "sd_n"), 0)
  expect_equal(sens(capped, "pct_hit_max_n"), 100)
  capped2 <- seq_mr_smd_sensitivity(true_value = 0.5, A = 50000,
                                    c = c(2, 8), G = 5, seed = 113,
                                    batch = 4, max_n = 33)
  expect_equal(sens(capped2, "mean_n1"), 33)
  expect_equal(sens(capped2, "pct_hit_max_n"), 100)
  # The costlier group keeps the size paired with 33:
  # 33 / sqrt(1 + 6 * 33^2 / 50000) = 31.03, so 32.
  expect_equal(sens(capped2, "mean_n2"), 32)
  expect_error(seq_mr_cv_sensitivity(0.4, A = 5e6, c = 10, G = 5,
                                     max_n = 50),
               "'max_n' (50) is smaller than the pilot sample size (82)",
               fixed = TRUE)
})

test_that("the wrappers validate every argument", {
  f <- function(...) seq_mr_cv_sensitivity(..., seed = 113)
  expect_error(f(true_value = 0, A = 2000, c = 1, G = 5), "'true_value'")
  expect_error(f(true_value = -0.5, A = 2000, c = 1, G = 5), "'true_value'")
  expect_error(f(true_value = 0.5, A = NA, c = 1, G = 5), "'A'")
  expect_error(f(true_value = 0.5, A = 2000, c = c(1, 2), G = 5), "'c'")
  expect_error(f(true_value = 0.5, A = 2000, c = 1, G = 2.5),
               "'G' must be a single whole number of at least 2")
  expect_error(f(true_value = 0.5, A = 2000, c = 1, G = 1), "'G'")
  expect_error(f(true_value = 0.5, A = 2000, c = 1, G = 5, batch = 0),
               "'batch'")
  expect_error(f(true_value = 0.5, A = 2000, c = 1, G = 5, batch = 1.5),
               "'batch'")
  expect_error(f(true_value = 0.5, A = 2000, c = 1, G = 5, max_n = "a"),
               "'max_n'")
  expect_error(f(true_value = 0.5, A = 2000, c = 1, G = 5, gamma = 0.5),
               "'gamma'")
  expect_error(f(true_value = 0.5, A = 2000, c = 1, G = 5, m_0 = 3), "'m_0'")
  expect_error(f(true_value = 0.5, A = 2000, c = 1, G = 5,
                 distribution = "skewed"),
               "'distribution' must be one of \"gamma\", \"lognormal\", \"normal\"",
               fixed = TRUE)
  expect_error(f(true_value = 0.5, A = 2000, c = 1, G = 5,
                 variance = "none"), "'variance'")
  expect_error(seq_mr_cv_sensitivity(0.5, A = 2000, c = 1, G = 5,
                                     seed = 1.5), "'seed'")
  # A seed beyond the integer range is reported by name, with no coercion
  # warning.
  expect_warning(
    expect_error(seq_mr_cv_sensitivity(0.5, A = 2000, c = 1, G = 5,
                                       seed = 1e10),
                 "'seed' must be NULL or a single whole number .*; got 1e\\+10"),
    NA)
  expect_error(seq_mr_smd_sensitivity(0.5, A = 500, c = 2, G = 5,
                                      seed = "a"),
               paste("'seed' must be NULL or a single whole number that",
                     "set.seed() accepts; got \"a\"."),
               fixed = TRUE)
  expect_error(seq_mr_smd_sensitivity(0.5, A = 500, c = 2, G = 5,
                                      seed = c(1, 2)), "'seed'")
  h <- function(...) seq_mr_smd_sensitivity(..., seed = 113)
  expect_error(h(true_value = NA, A = 500, c = 2, G = 5), "'true_value'")
  expect_error(h(true_value = 0.5, A = 500, c = c(1, 2, 3), G = 5), "'c'")
  expect_error(h(true_value = 0.5, A = 500, c = 2, G = 5,
                 distribution = "gamma"),
               "'distribution' must be one of \"normal\", \"skewed\"",
               fixed = TRUE)
  expect_error(h(true_value = 0.5, A = 500, c = 2, G = 5, m_0 = 4.5),
               "'m_0'")
  expect_error(h(true_value = 0.5, A = 500, c = 2, G = 0), "'G'")
})

test_that("a population value too extreme to simulate is refused in the siblings' words", {
  # The accuracy in parameter estimation siblings refuse a value whose
  # generating distribution cannot be formed, and these two functions refuse
  # it in the same words. Unchecked, 1e155 with lognormal data would stop
  # with R's "missing value where TRUE/FALSE needed" and a warning from
  # rlnorm(), and 1e154 for the standardized mean difference would return
  # rows of Inf and NA.
  message_of <- function(expr) tryCatch(expr, error = conditionMessage)
  cv <- function(...) {
    seq_mr_cv_sensitivity(..., A = 20000, c = 75, G = 4, seed = 1,
                          max_n = 100)
  }
  cv_sibling <- function(...) {
    seq_aipe_cv_sensitivity(..., width = 0.2, G = 4, seed = 1, max_n = 100)
  }
  expect_error(
    cv(1e155, distribution = "lognormal"),
    paste0("^'true_value' is too extreme to simulate: the lognormal ",
           "distribution would need a log standard deviation, ",
           "sqrt\\(log\\(1 \\+ true_value\\^2\\)\\), of Inf, and that must ",
           "be a positive finite number; got 1e\\+155\\.$")
  )
  expect_identical(message_of(cv(1e155, distribution = "lognormal")),
                   message_of(cv_sibling(1e155, distribution = "lognormal")))
  # A gamma shape of 1 / (1e170)^2 underflows to zero.
  expect_error(
    cv(1e170, distribution = "gamma"),
    paste0("^'true_value' is too extreme to simulate: the gamma ",
           "distribution would need a shape, 1 / true_value\\^2, of 0, and ",
           "that must be a positive finite number; got 1e\\+170\\.$")
  )
  expect_identical(message_of(cv(1e170, distribution = "gamma")),
                   message_of(cv_sibling(1e170, distribution = "gamma")))
  # 1 + (1e-9)^2 is 1 in double precision, so the lognormal log standard
  # deviation would be zero: every observation equal.
  expect_error(cv(1e-9, distribution = "lognormal"),
               "would need a log standard deviation, .*, of 0, and that")

  smd <- function(...) {
    seq_mr_smd_sensitivity(..., A = 2000, c = 5, G = 4, seed = 1,
                           max_n = 100)
  }
  expect_error(
    smd(1e154),
    paste0("^'true_value' must be less than 1e\\+08 in absolute value; got ",
           "1e\\+154\\. The generated data have unit variance")
  )
  expect_identical(
    message_of(smd(1e154)),
    message_of(seq_aipe_smd_sensitivity(1e154, width = 1, G = 4, seed = 1,
                                        max_n = 100))
  )
  expect_error(smd(-1e8, distribution = "skewed"),
               "^'true_value' must be less than 1e\\+08 in absolute value")
  # Just inside the limit the study runs and every row is a finite number.
  inside <- seq_mr_smd_sensitivity(9.9e7, A = 2000, c = 5, G = 2, seed = 1,
                                   max_n = 20)
  expect_true(all(is.finite(inside$value)))
})

test_that("the benchmark of a study is a positive finite number, or the call stops", {
  # A coefficient of variation that every other check accepts can still put
  # the benchmark outside double precision: its fourth power overflows
  # beyond about 1e77 and underflows below about 1e-81, and the summary
  # would then hold optimal_n = Inf (or 0, with ratio_n = Inf) and rows of
  # NaN.
  cv <- function(...) {
    seq_mr_cv_sensitivity(..., A = 20000, c = 75, G = 2, seed = 1,
                          max_n = 20)
  }
  expect_error(
    cv(1e80, distribution = "normal"),
    paste0("^'true_value' is too extreme to simulate with this 'A' and ",
           "'c': the benchmark the procedure is compared with must be ",
           "positive finite numbers, and optimal_n is Inf; got ",
           "true_value = 1e\\+80, A = 20000, c = 75\\.$")
  )
  expect_error(cv(1e-170, distribution = "normal"),
               "positive finite numbers, and optimal_n is 0; got true_value = 1e-170")
  expect_error(cv(1e40, distribution = "lognormal"),
               "and optimal_n is Inf; got true_value = 1e\\+40")
  expect_error(cv(1e100, distribution = "gamma"), "and optimal_n is NaN")
  # The helper names the first entry that is not a positive finite number.
  expect_true(.mrpe_check_benchmark(c(optimal_n = 12.5, minimum_risk = 50),
                                    list(true_value = 0.5)))
  expect_error(.mrpe_check_benchmark(c(optimal_n1 = 12, optimal_n2 = Inf,
                                       minimum_risk = NaN),
                                     list(true_value = 3, A = 2, c1 = 1)),
               "and optimal_n2 is Inf; got true_value = 3, A = 2, c1 = 1\\.$")
  expect_error(.mrpe_check_benchmark(c(optimal_n = 12, minimum_risk = NA),
                                     list(true_value = 3)),
               "and minimum_risk is NA; got true_value = 3\\.$")
  # Over fifteen orders of magnitude on either side of ordinary values and
  # the three populations: a call either stops with one of these messages
  # or returns a summary with no Inf or NaN in it (NA is documented for
  # the estimate and risk rows when no replication gave an estimate).
  for (dist in c("gamma", "lognormal", "normal")) {
    for (kappa in 10^c(-200, -100, -15, -6, -2, 0, 2, 6, 15, 60, 100, 200)) {
      out <- tryCatch(cv(kappa, distribution = dist), error = conditionMessage)
      info <- paste(dist, kappa)
      if (is.character(out)) {
        expect_match(out, "^'true_value' is too extreme to simulate",
                     info = info)
      } else {
        expect_false(any(is.infinite(out$value) | is.nan(out$value)),
                     info = info)
        expect_false(anyNA(out$value[match(c("optimal_n", "minimum_risk",
                                             "mean_n", "ratio_n"),
                                           out$term)]), info = info)
      }
    }
  }
})
