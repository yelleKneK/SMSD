# Tests for the two stage fixed width procedures: ts_fwd_norm_mean() (Stein),
# ts_fwd_dfr_mean() (asymptotic), ts_fwd_nexp_mean() (negative exponential
# location), and their pilot functions.
#
# Every expected value is derived here, independently of the package code: a
# quantile from a printed table or a closed form, arithmetic shown in a
# comment, or an exact operating characteristic computed from the
# distribution of the first stage statistic that fixes N.
#
# A few expectations involve no package function at all (a quantile against
# its printed table value, the ratios and limits quoted in the help pages
# recomputed from R's distribution functions). They anchor the numbers the
# help pages print and the constants the other expectations rely on; they
# are not checks of the package code, and each is marked "documented
# number" where it stands.

# The helpers ts_n_dist() and ts_run_two_stage() are in helper-ts.R.

# ---- ts_fwd_norm_mean ----------------------------------------------------

test_that("ts_fwd_norm_mean criterion is Stein's rule with the two-sided t quantile", {
  x <- 1:10                       # mean 5.5, variance 55 / 6 = 9.1667
  out <- ts_fwd_norm_mean(x, d = 1, conf_level = 0.95)
  v <- stats::setNames(out$value, out$term)
  # t table, 9 degrees of freedom, upper 0.025 point: 2.262.
  # 2.262^2 * 9.1667 / 1^2 = 46.90, so N = 47.
  # (One-sided 1.833: 30.8 -> 31. Normal 1.960: 35.2 -> 36.)
  expect_equal(v[["optimal_n"]], 47)
  expect_equal(v[["optimal_n"]], ceiling(2.262^2 * (55 / 6)))
  expect_equal(attr(out, "criterion_n"), 47)
  expect_equal(v[["mean"]], 5.5)
  expect_equal(v[["lower_limit"]], 4.5)
  expect_equal(v[["upper_limit"]], 6.5)
  expect_equal(v[["current_n"]], 10)
  expect_equal(attr(out, "width"), 2)             # full width 2 d
  expect_equal(attr(out, "conf_level"), 0.95)
  expect_false(attr(out, "is_satisfied"))
  expect_identical(out$term, c("lower_limit", "mean", "upper_limit",
                               "optimal_n", "current_n"))
  expect_identical(attr(out, "method"),
                   "two stage fixed width interval (normal mean)")
})

test_that("ts_fwd_norm_mean fixes N from the first stage and stops exactly there", {
  withr::local_seed(113)
  x <- 1:10
  # N = 47 whatever the second stage holds: 36 more is one short, 37 is enough.
  short <- ts_fwd_norm_mean(x, rnorm(36, 100, 50), d = 1)
  enough <- ts_fwd_norm_mean(x, rnorm(37, 100, 50), d = 1)
  expect_equal(attr(short, "criterion_n"), 47)
  expect_equal(attr(enough, "criterion_n"), 47)
  expect_equal(attr(short, "current_n"), 46)
  expect_false(attr(short, "is_satisfied"))
  expect_true(attr(enough, "is_satisfied"))
  # The estimate is the mean of both stages together.
  second <- c(3, 9, 12)
  both <- ts_fwd_norm_mean(x, second, d = 1)
  expect_equal(both$value[both$term == "mean"], 79 / 13)   # (55 + 24) / 13
  expect_equal(both$value[both$term == "lower_limit"], 79 / 13 - 1)
  expect_equal(both$value[both$term == "upper_limit"], 79 / 13 + 1)
  expect_equal(both$value[both$term == "current_n"], 13)
  expect_identical(attr(both, "current_n"), 13L)
  # When the first stage is already enough, N is the first stage size.
  wide <- ts_fwd_norm_mean(x, d = 100)
  expect_equal(attr(wide, "criterion_n"), 10)
  expect_true(attr(wide, "is_satisfied"))
})

test_that("Stein's two stage interval has exact coverage (seeded Monte Carlo)", {
  withr::local_seed(113)
  m <- 5; sigma <- 2; d <- 0.5; mu <- 10; reps <- 2500
  res <- replicate(reps, ts_run_two_stage(
    ts_fwd_norm_mean, function(k) rnorm(k, mu, sigma), m, mu, d = d))
  coverage <- mean(res["covered", ])
  mean_n <- mean(res["n", ])

  # Exact values. N = max(m, ceiling(t^2 S^2 / d^2)) with (m - 1) S^2 /
  # sigma^2 chi square on m - 1 degrees of freedom, independent of the
  # sample means; given N the final mean is N(mu, sigma^2 / N), so the
  # coverage is E[2 pnorm(d sqrt(N) / sigma) - 1].
  dist <- ts_n_dist(m, m - 1,
                           stats::qt(0.975, m - 1)^2 * sigma^2 /
                             (d^2 * (m - 1)))
  exact_cov <- sum(dist$p * (2 * stats::pnorm(d * sqrt(dist$n) / sigma) - 1))
  exact_n <- sum(dist$p * dist$n)
  expect_gte(exact_cov, 0.95)                    # the guarantee itself
  expect_equal(exact_cov, 0.9513, tolerance = 5e-4)

  # Tolerances are four Monte Carlo standard errors.
  se_cov <- sqrt(exact_cov * (1 - exact_cov) / reps)      # about 0.0043
  se_n <- stats::sd(res["n", ]) / sqrt(reps)
  expect_lt(abs(coverage - exact_cov), 4 * se_cov)
  expect_gte(coverage, 0.95 - 4 * se_cov)
  expect_lt(abs(mean_n - exact_n), 4 * se_n)
  # Documented number (not a check of the code): known variance optimum
  # C = z^2 sigma^2 / d^2 = 61.46; Stein with m = 5 needs about
  # t^2 / z^2 = 2.0 times that.
  expect_equal(exact_n / (stats::qnorm(0.975)^2 * sigma^2 / d^2), 2.015,
               tolerance = 1e-3)
})

test_that("ts_fwd_norm_mean_pilot follows the documented rule and minimum", {
  # Rule with gamma = 1: ceiling(z / d). z_{0.025} = 1.959964.
  expect_identical(ts_fwd_norm_mean_pilot(d = 0.1), 20L)   # ceiling(19.6)
  expect_identical(ts_fwd_norm_mean_pilot(d = 0.01), 196L) # ceiling(195.996)
  # Below the minimum the default m_0 = 10 is returned.
  expect_identical(ts_fwd_norm_mean_pilot(d = 2), 10L)     # ceiling(0.98) = 1
  expect_identical(ts_fwd_norm_mean_pilot(d = 2, m_0 = 2), 2L)
  expect_identical(ts_fwd_norm_mean_pilot(d = 2, m_0 = 25), 25L)
  # gamma = 0.5: (z / d)^(2 / 1.5) = 19.59964^(4/3) = 52.8, so 53.
  expect_identical(ts_fwd_norm_mean_pilot(d = 0.1, gamma = 0.5), 53L)
  # conf_level enters through z: z_{0.005} = 2.575829, ceiling(25.76) = 26.
  expect_identical(ts_fwd_norm_mean_pilot(d = 0.1, conf_level = 0.99), 26L)
  # Lower order than the known variance optimum (z / d)^2 as d shrinks.
  expect_lt(ts_fwd_norm_mean_pilot(d = 0.01) / (1.959964 / 0.01)^2, 0.01)
  expect_identical(names(formals(ts_fwd_norm_mean_pilot)),
                   c("d", "conf_level", "m_0", "gamma"))
})

# ---- ts_fwd_dfr_mean -----------------------------------------------------

test_that("ts_fwd_dfr_mean criterion uses the two-sided normal quantile", {
  x <- 1:10                       # variance 55 / 6
  out <- ts_fwd_dfr_mean(x, d = 1, conf_level = 0.95)
  v <- stats::setNames(out$value, out$term)
  # 1.959964^2 * 9.1667 = 35.21, so N = 36 (Stein's t rule gives 47).
  expect_equal(v[["optimal_n"]], 36)
  expect_equal(v[["optimal_n"]], ceiling(1.959964^2 * (55 / 6)))
  expect_equal(v[["mean"]], 5.5)
  expect_equal(v[["lower_limit"]], 4.5)
  expect_equal(v[["upper_limit"]], 6.5)
  expect_equal(attr(out, "width"), 2)
  expect_identical(attr(out, "conf_level"), 0.95)
  expect_identical(attr(out, "method"),
                   "two stage fixed width interval (distribution free mean)")
  expect_identical(out$term, c("lower_limit", "mean", "upper_limit",
                               "optimal_n", "current_n"))
  expect_equal(v[["current_n"]], 10)
  expect_false(attr(out, "is_satisfied"))
  expect_true(attr(ts_fwd_dfr_mean(x, rep(5, 26), d = 1), "is_satisfied"))
  expect_false(attr(ts_fwd_dfr_mean(x, rep(5, 25), d = 1), "is_satisfied"))
  # With a second stage the rows describe both stages together: 13
  # observations with mean (55 + 3 + 9 + 12) / 13, and N still 36.
  both <- ts_fwd_dfr_mean(x, c(3, 9, 12), d = 1, conf_level = 0.90)
  vb <- stats::setNames(both$value, both$term)
  expect_equal(vb[["current_n"]], 13)
  expect_identical(attr(both, "current_n"), 13L)
  expect_equal(vb[["mean"]], 79 / 13)
  expect_equal(vb[["lower_limit"]], 79 / 13 - 1)
  expect_equal(vb[["upper_limit"]], 79 / 13 + 1)
  expect_identical(attr(both, "conf_level"), 0.90)
  # 90%: z = 1.644854; 2.705543 * 9.1667 = 24.80, so N = 25.
  expect_equal(attr(ts_fwd_dfr_mean(x, d = 1, conf_level = 0.90),
                    "criterion_n"), 25)
})

test_that("ts_fwd_dfr_mean undercovers by the documented amount (seeded Monte Carlo)", {
  withr::local_seed(113)
  sigma <- 2; d <- 0.5; mu <- 10; reps <- 2500
  m <- ts_fwd_dfr_mean_pilot(d = d)
  expect_identical(m, 30L)
  res <- replicate(reps, ts_run_two_stage(
    ts_fwd_dfr_mean, function(k) rnorm(k, mu, sigma), m, mu, d = d))
  coverage <- mean(res["covered", ])

  # Exact coverage for normal data, as for Stein's rule but with z in place
  # of t: the procedure is asymptotic, so this is below 0.95.
  dist <- ts_n_dist(m, m - 1,
                           stats::qnorm(0.975)^2 * sigma^2 / (d^2 * (m - 1)))
  exact_cov <- sum(dist$p * (2 * stats::pnorm(d * sqrt(dist$n) / sigma) - 1))
  expect_equal(exact_cov, 0.9417, tolerance = 5e-4)
  expect_lt(exact_cov, 0.95)
  # Tolerance: four Monte Carlo standard errors (about 0.019).
  expect_lt(abs(coverage - exact_cov),
            4 * sqrt(exact_cov * (1 - exact_cov) / reps))

  # Documented numbers (not checks of the code): the limit for a narrow
  # target quoted in the help pages, P(|t_{m-1}| <= z), and its agreement
  # with the exact coverage at a narrow target.
  limit <- function(m) 2 * stats::pt(stats::qnorm(0.975), m - 1) - 1
  expect_equal(limit(2), 0.700, tolerance = 1e-3)
  expect_equal(limit(10), 0.918, tolerance = 1e-3)
  expect_equal(limit(30), 0.940, tolerance = 1e-3)
  expect_equal(limit(100), 0.947, tolerance = 1e-3)
  narrow <- ts_n_dist(30, 29, stats::qnorm(0.975)^2 * 4 / (0.05^2 * 29))
  expect_equal(
    sum(narrow$p * (2 * stats::pnorm(0.05 * sqrt(narrow$n) / 2) - 1)),
    limit(30), tolerance = 1e-3)
  # The normal row of the table (standard deviation 2, first stage 30), to
  # the three decimals printed, and the statement that coverage falls below
  # the nominal level once the known variance sample size C is more than
  # about 1.2 times the first stage (C = 1.15 m and C = 1.25 m bracket it).
  exact_at_c <- function(C) {
    dist <- ts_n_dist(30, 29, C / 29)
    sum(dist$p * (2 * stats::pnorm(stats::qnorm(0.975) *
                                     sqrt(dist$n / C)) - 1))
  }
  big_c <- stats::qnorm(0.975)^2 * 4 / c(1, 0.5, 0.2, 0.1)^2
  expect_equal(round(vapply(big_c, exact_at_c, numeric(1)), 3),
               c(0.994, 0.942, 0.941, 0.940), tolerance = 1e-12)
  expect_gt(exact_at_c(1.15 * 30), 0.95)
  expect_lt(exact_at_c(1.25 * 30), 0.95)
})

test_that("ts_fwd_dfr_mean_pilot follows the documented rule and minimum", {
  expect_identical(ts_fwd_dfr_mean_pilot(d = 0.5), 30L)    # ceiling(3.92) = 4
  expect_identical(ts_fwd_dfr_mean_pilot(d = 0.05), 40L)   # ceiling(39.2)
  expect_identical(ts_fwd_dfr_mean_pilot(d = 0.5, m_0 = 2), 4L)
  # 90%: z = 1.644854, ceiling(1.644854 / 0.02) = ceiling(82.24) = 83.
  expect_identical(ts_fwd_dfr_mean_pilot(d = 0.02, conf_level = 0.90), 83L)
  expect_identical(names(formals(ts_fwd_dfr_mean_pilot)),
                   c("d", "conf_level", "m_0", "gamma"))
})

# ---- ts_fwd_nexp_mean ----------------------------------------------------

test_that("ts_fwd_nexp_mean criterion uses the one-sided F quantile and the first stage scale", {
  x <- c(3, 4, 5, 7, 11)          # minimum 3; U = (0 + 1 + 2 + 4 + 8) / 4 = 3.75
  out <- ts_fwd_nexp_mean(x, d = 0.5, conf_level = 0.95)
  v <- stats::setNames(out$value, out$term)
  # Upper alpha point of F(2, 2m - 2) in closed form:
  # (m - 1) (alpha^(-1 / (m - 1)) - 1) = 4 (0.05^(-1/4) - 1) = 4.4590
  # (F table, 2 and 8 degrees of freedom, 0.95: 4.46).
  f_upper <- 4 * (0.05^(-1 / 4) - 1)
  expect_equal(f_upper, 4.459, tolerance = 1e-4)   # documented number
  # 4.4590 * 3.75 / 0.5 = 33.44, so N = 34.
  expect_equal(v[["optimal_n"]], 34)
  expect_equal(v[["optimal_n"]], ceiling(f_upper * 3.75 / 0.5))
  expect_equal(v[["location"]], 3)
  expect_equal(v[["lower_limit"]], 2.5)
  expect_equal(v[["upper_limit"]], 3)
  expect_equal(attr(out, "width"), 0.5)           # d is the full width here
  expect_identical(attr(out, "conf_level"), 0.95)
  expect_identical(
    attr(out, "method"),
    "two stage fixed width interval (negative exponential location)")
  expect_equal(v[["current_n"]], 5)
  expect_false(attr(out, "is_satisfied"))
  expect_identical(out$term, c("lower_limit", "location", "upper_limit",
                               "optimal_n", "current_n"))
  # 99%: 4 (0.01^(-1/4) - 1) = 8.6491; 8.6491 * 3.75 / 0.5 = 64.87, so 65.
  out99 <- ts_fwd_nexp_mean(x, d = 0.5, conf_level = 0.99)
  expect_equal(attr(out99, "criterion_n"), 65)
  expect_identical(attr(out99, "conf_level"), 0.99)
  # The location estimate is the minimum over both stages; N is not moved by
  # the second stage.
  both <- ts_fwd_nexp_mean(x, c(2.2, 9, 6), d = 0.5)
  expect_equal(both$value[both$term == "location"], 2.2)
  expect_equal(both$value[both$term == "lower_limit"], 1.7)
  expect_equal(both$value[both$term == "upper_limit"], 2.2)
  expect_equal(both$value[both$term == "current_n"], 8)
  expect_identical(attr(both, "current_n"), 8L)
  expect_equal(attr(both, "criterion_n"), 34)
  expect_true(attr(ts_fwd_nexp_mean(x, rep(6, 29), d = 0.5), "is_satisfied"))
  expect_false(attr(ts_fwd_nexp_mean(x, rep(6, 28), d = 0.5), "is_satisfied"))
})

test_that("the negative exponential two stage interval has exact coverage (seeded Monte Carlo)", {
  withr::local_seed(113)
  m <- 5; sigma <- 2; theta <- 3; d <- 0.5; alpha <- 0.05; reps <- 2500
  res <- replicate(reps, ts_run_two_stage(
    ts_fwd_nexp_mean, function(k) theta + rexp(k, rate = 1 / sigma), m,
    theta, d = d))
  coverage <- mean(res["covered", ])
  mean_n <- mean(res["n", ])

  # Exact values. N = max(m, ceiling(F U / d)) with 2 (m - 1) U / sigma chi
  # square on 2m - 2 degrees of freedom, independent of the minimum; given N
  # the interval misses with probability exp(-N d / sigma).
  f_upper <- (m - 1) * (alpha^(-1 / (m - 1)) - 1)
  dist <- ts_n_dist(m, 2 * m - 2, f_upper * sigma / (d * 2 * (m - 1)))
  exact_cov <- sum(dist$p * (1 - exp(-dist$n * d / sigma)))
  exact_n <- sum(dist$p * dist$n)
  expect_gte(exact_cov, 0.95)                    # the guarantee itself
  # Tolerances are four Monte Carlo standard errors.
  se_cov <- sqrt(exact_cov * (1 - exact_cov) / reps)
  expect_lt(abs(coverage - exact_cov), 4 * se_cov)
  expect_gte(coverage, 0.95 - 4 * se_cov)
  expect_lt(abs(mean_n - exact_n), 4 * stats::sd(res["n", ]) / sqrt(reps))
  # Documented number (not a check of the code): the quantile makes the
  # bound tight, E[exp(-F U / sigma)] = alpha.
  expect_equal((1 + f_upper / (m - 1))^(-(m - 1)), alpha)
})

test_that("ts_fwd_nexp_mean_pilot follows the documented rule and minimum", {
  # Rule with gamma = 1: ceiling(sqrt(log(1 / alpha) / d)).
  # 99%, d = 0.01: sqrt(4.60517 / 0.01) = 21.46, so 22.
  expect_identical(ts_fwd_nexp_mean_pilot(d = 0.01, conf_level = 0.99), 22L)
  # 95%, d = 0.2: sqrt(2.995732 / 0.2) = 3.87, so 4; the default minimum wins.
  expect_identical(ts_fwd_nexp_mean_pilot(d = 0.2), 10L)
  expect_identical(ts_fwd_nexp_mean_pilot(d = 0.2, m_0 = 2), 4L)
  # gamma = 3: (2.995732 / 0.001)^(1/4) = 7.40, so 8.
  expect_identical(ts_fwd_nexp_mean_pilot(d = 0.001, m_0 = 2, gamma = 3), 8L)
  expect_identical(names(formals(ts_fwd_nexp_mean_pilot)),
                   c("d", "conf_level", "m_0", "gamma"))
})

# ---- argument checks and degenerate data ---------------------------------

test_that("the two stage fixed width checks have no gamma argument", {
  for (f in list(ts_fwd_norm_mean, ts_fwd_dfr_mean, ts_fwd_nexp_mean)) {
    expect_identical(names(formals(f)),
                     c("pilot_data", "additional_data", "d", "conf_level",
                       "na.rm"))
  }
  expect_error(ts_fwd_norm_mean(1:10, d = 1, gamma = 1), "unused argument")
  expect_error(ts_fwd_dfr_mean(1:10, d = 1, gamma = -5), "unused argument")
  expect_error(ts_fwd_nexp_mean(1:10, d = 1, gamma = "junk"),
               "unused argument")
})

test_that("two stage fixed width pilots check every argument", {
  for (f in list(ts_fwd_norm_mean_pilot, ts_fwd_dfr_mean_pilot,
                 ts_fwd_nexp_mean_pilot)) {
    expect_error(f(), "'d' must be supplied")
    expect_error(f(d = -1), "'d' must be a single finite number")
    expect_error(f(d = 0), "'d'")
    expect_error(f(d = NA_real_), "'d'")
    expect_error(f(d = Inf), "'d'")
    expect_error(f(d = c(0.1, 0.2)), "'d'")
    expect_error(f(d = "a"), "'d'")
    expect_error(f(d = 0.5, conf_level = 1), "'conf_level'")
    expect_error(f(d = 0.5, conf_level = 0), "'conf_level'")
    expect_error(f(d = 0.5, conf_level = NA_real_), "'conf_level'")
    expect_error(f(d = 0.5, m_0 = 1), "'m_0' must be a single whole number of at least 2")
    expect_error(f(d = 0.5, m_0 = 2.5), "'m_0'")
    expect_error(f(d = 0.5, m_0 = NA), "'m_0'")
    expect_error(f(d = 0.5, gamma = 0), "'gamma'")
    expect_error(f(d = 0.5, gamma = -1), "'gamma'")
    expect_error(f(d = 0.5, gamma = NA_real_), "'gamma'")
    expect_error(f(d = 0.5, gamma = "a"), "'gamma'")
    # A target too narrow to represent is an error, not NA, and the message
    # says what to look at (not the shared default about "the precision
    # target and the other arguments").
    expect_error(f(d = 1e-30), "implied pilot sample size .* is too large to represent")
    expect_error(f(d = 1e-30), "Check that 'd', the (half|full) width")
    expect_type(f(d = 0.5), "integer")
  }
  expect_error(ts_fwd_norm_mean_pilot(d = -1), "half width")
  expect_error(ts_fwd_nexp_mean_pilot(d = -1), "full width")
  expect_error(ts_fwd_norm_mean_pilot(d = 1e-30),
               "'d', the half width of the confidence interval, is in the units of the data")
  expect_error(ts_fwd_dfr_mean_pilot(d = 1e-30),
               "'d', the half width of the confidence interval, is in the units of the data")
  expect_error(ts_fwd_nexp_mean_pilot(d = 1e-30),
               "'d', the full width of the confidence interval, is in the units of the data")
})

test_that("a pilot rule that is an exact power comes back as itself for any gamma", {
  # The rule is ceiling(K^(1 / (1 + gamma))). When K is an exact power the
  # answer is its root, but K^(1 / (1 + gamma)) is computed a few units in
  # the last place off (4096^(1/4) need not be exactly 8), and a plain
  # ceiling of 8.000000000000002 is 9. The rounding ignores such
  # representation error.
  z <- stats::qnorm(0.975)
  # gamma = 1: K = (z / d)^2 with d = z / r, so the rule is r.
  for (r in c(7, 20, 50, 196)) {
    expect_identical(ts_fwd_norm_mean_pilot(d = z / r, m_0 = 2), as.integer(r))
    expect_identical(ts_fwd_dfr_mean_pilot(d = z / r, m_0 = 2), as.integer(r))
  }
  # gamma = 3: K = (z / d)^2 = r^4 when d = z / r^2, so the rule is r.
  # gamma = 0.25: K = r^(5/4) when d = z / r^(5/8).
  for (r in 3:60) {
    expect_identical(ts_fwd_norm_mean_pilot(d = z / r^2, m_0 = 2, gamma = 3),
                     as.integer(r))
    expect_identical(ts_fwd_dfr_mean_pilot(d = z / r^(5 / 8), m_0 = 2,
                                           gamma = 0.25), as.integer(r))
    # Negative exponential: K = log(1 / alpha) / d = r^4 when
    # d = log(20) / r^4 at 95%.
    expect_identical(ts_fwd_nexp_mean_pilot(d = log(20) / r^4, m_0 = 2,
                                            gamma = 3), as.integer(r))
  }
  # Just above an exact power the rule still rounds up: 8.01^4 = 4116.5.
  expect_identical(ts_fwd_nexp_mean_pilot(d = log(20) / 8.01^4, m_0 = 2,
                                          gamma = 3), 9L)
})

test_that("two stage fixed width checks check every argument", {
  x <- c(3, 4, 5, 7, 11, 6, 8, 5, 9, 4)
  for (f in list(ts_fwd_norm_mean, ts_fwd_dfr_mean, ts_fwd_nexp_mean)) {
    expect_error(f(d = 1), "'pilot_data' must be supplied")
    expect_error(f(x), "'d' must be supplied")
    expect_error(f(x, d = -1), "'d'")
    expect_error(f(x, d = NA_real_), "'d'")
    expect_error(f(x, d = Inf), "'d'")
    expect_error(f(x, d = 1, conf_level = 1.2), "'conf_level'")
    expect_error(f(x, d = 1, conf_level = NA_real_), "'conf_level'")
    expect_error(f(letters, d = 1), "'pilot_data' must be a numeric vector")
    expect_error(f(matrix(x, 5), d = 1), "'pilot_data' must be a numeric vector, not a matrix")
    expect_error(f(data.frame(x = x), d = 1), "'pilot_data'")
    expect_error(f(7, d = 1), "'pilot_data' needs at least 2 observations; it has 1")
    # The advice is for a two stage design: a larger first stage, not "check
    # the stopping rule again" (the first stage alone fixes N).
    expect_error(f(7, d = 1), "Enlarge the first stage; ts_fwd_(norm|dfr|nexp)_mean_pilot\\(\\) returns a suitable size")
    expect_false(grepl("check the stopping rule again",
                       tryCatch(f(7, d = 1), error = conditionMessage)))
    expect_error(f(c(NA, 7, NA), d = 1), "'pilot_data' needs at least 2 observations; it has 1")
    expect_error(f(c(x, Inf), d = 1), "infinite")
    expect_error(f(x, "b", d = 1), "'additional_data' must be a numeric vector")
    expect_error(f(x, c(1, Inf), d = 1), "'additional_data' contains infinite")
    expect_error(f(x, d = 1, na.rm = NA), "'na.rm' must be TRUE or FALSE")
    expect_error(f(c(x, NA), d = 1, na.rm = FALSE), "missing values")
    expect_error(f(x, c(5, NA), d = 1, na.rm = FALSE), "missing values")
    # NA removed by default: same answer as without it.
    expect_equal(f(c(x, NA), c(NA, 6, 7), d = 1)$value,
                 f(x, c(6, 7), d = 1)$value)
    # A criterion too large to represent is an error, never NA, and the
    # message says what to look at.
    expect_error(f(x, d = 1e-12), "implied total sample size .* is too large to represent")
    expect_error(f(x, d = 1e-12), "Check that 'd', the (half|full) width")
    # Integer input is handled as doubles.
    expect_equal(f(as.integer(x), d = 1)$value, f(as.numeric(x), d = 1)$value)
    out <- f(x, d = 1)
    expect_s3_class(out, "smsd_check")
    expect_true(is.logical(attr(out, "is_satisfied")) &&
                  !is.na(attr(out, "is_satisfied")))
    expect_true(all(is.finite(out$value)))
  }
  expect_error(ts_fwd_norm_mean(x, d = -1), "half width")
  expect_error(ts_fwd_nexp_mean(x, d = -1), "full width")
  # Each function names its own pilot function and its own measure of spread.
  expect_error(ts_fwd_norm_mean(7, d = 1),
               "A variance cannot be computed from fewer. Enlarge the first stage; ts_fwd_norm_mean_pilot\\(\\)")
  expect_error(ts_fwd_dfr_mean(7, d = 1),
               "A variance cannot be computed from fewer. Enlarge the first stage; ts_fwd_dfr_mean_pilot\\(\\)")
  expect_error(ts_fwd_nexp_mean(7, d = 1),
               "The scale cannot be estimated from fewer. Enlarge the first stage; ts_fwd_nexp_mean_pilot\\(\\)")
  expect_error(ts_fwd_norm_mean(x, d = 1e-12),
               "Check that 'd', the half width of the confidence interval, is in the units of the data; the target as given calls for more than two billion observations\\.")
  expect_error(ts_fwd_dfr_mean(x, d = 1e-12),
               "Check that 'd', the half width of the confidence interval, is in the units of the data; the target as given calls for more than two billion observations\\.")
  expect_error(ts_fwd_nexp_mean(x, d = 1e-12),
               "Check that 'd', the full width of the confidence interval, is in the units of the data; the target as given calls for more than two billion observations\\.")
})

test_that("a total too large to represent names the first stage when it is the cause", {
  msg <- function(expr) tryCatch(expr, error = conditionMessage)
  # Stein's rule with a first stage of two at a confidence level of 0.9999.
  # The t distribution on one degree of freedom is the Cauchy distribution,
  # so its upper 0.00005 point is 1 / tan(0.00005 pi) = 6366.2 (the normal
  # quantile is 3.8906). x = (1, 30) has variance 29^2 / 2 = 420.5, so with
  # d = 0.5 the rule gives 6366.2^2 * 420.5 / 0.25 = 6.8e10, too large to
  # represent, although d is in the units of the data: with the normal
  # quantile it is 3.8906^2 * 420.5 / 0.25 = 25459.9, so 25460. The message
  # used to say only to check the units of d.
  expect_equal(1 / tan(0.00005 * pi), 6366.198, tolerance = 1e-6)
  stein <- msg(ts_fwd_norm_mean(c(1, 30), d = 0.5, conf_level = 0.9999))
  expect_match(stein, "The implied total sample size (68168891233) is too large to represent.",
               fixed = TRUE)
  expect_match(stein, paste0(
    "With 2 first stage observations the t quantile on 1 degree of freedom ",
    "is 6366 at this confidence level, where the normal quantile is 3.891; ",
    "with the normal quantile the rule would call for 25460 observations. ",
    "Enlarge the first stage (ts_fwd_norm_mean_pilot() returns a suitable ",
    "size), or check that 'd', the half width of the confidence interval, ",
    "is in the units of the data."), fixed = TRUE)
  expect_false(grepl("two billion", stein))
  # On two degrees of freedom the quantile is (2p - 1) / sqrt(2p(1 - p)),
  # which is 99.99 at p = 0.99995. x = (1, 30, 12) has variance 214.33;
  # with d = 0.02 the rule gives 99.99^2 * 214.33 / 0.0004 = 5.4e9.
  expect_match(msg(ts_fwd_norm_mean(c(1, 30, 12), d = 0.02, conf_level = 0.9999)),
               "With 3 first stage observations the t quantile on 2 degrees of freedom is 99.99 at this confidence level",
               fixed = TRUE)
  # The sample size named for the normal quantile is the rule's, so never
  # below the first stage: x = (1, 1.1) has variance 0.005; at a confidence
  # level of 1 - 1e-12 the normal quantile is 7.13 and 7.13^2 * 0.005 =
  # 0.25, less than one observation, so the rule would stop at the first
  # stage of 2 (while the Cauchy quantile, about 6.4e11, gives 2e21).
  expect_match(msg(ts_fwd_norm_mean(c(1, 1.1), d = 1, conf_level = 1 - 1e-12)),
    "where the normal quantile is 7.131; with the normal quantile the rule would call for 2 observations.",
               fixed = TRUE)
  # At the default level the same first stage is no error, and N is what
  # the Cauchy quantile gives: 1 / tan(0.025 pi) = 12.7062, and
  # 12.7062^2 * 420.5 / 0.25 = 271554.9, so 271555.
  expect_identical(attr(ts_fwd_norm_mean(c(1, 30), d = 0.5), "criterion_n"),
                   271555L)
  expect_equal((1 / tan(0.025 * pi))^2 * 420.5 / 0.25, 271554.9,
               tolerance = 1e-6)

  # The negative exponential rule has the same feature. The F distribution
  # on 2 and 2 degrees of freedom has upper tail 1 / (1 + f), so its upper
  # 0.0001 point is 9999, where the constant for a known scale is
  # log(10000) = 9.2103. x = (1, 30) has U = 29, so with d = 0.0001 the rule
  # gives 9999 * 29 / 0.0001 = 2899710000 (the largest integer is
  # 2147483647), and with the known scale constant 9.2103 * 29 / 0.0001 =
  # 2670998.7, so 2670999.
  nexp <- msg(ts_fwd_nexp_mean(c(1, 30), d = 1e-4, conf_level = 0.9999))
  expect_match(nexp, "The implied total sample size (2899710000) is too large to represent.",
               fixed = TRUE)
  expect_match(nexp, paste0(
    "With 2 first stage observations the F quantile on 2 and 2 degrees of ",
    "freedom is 9999 at this confidence level, where the constant for a ",
    "known scale, log(1 / alpha), is 9.21; with that constant the rule ",
    "would call for 2670999 observations. Enlarge the first stage ",
    "(ts_fwd_nexp_mean_pilot() returns a suitable size), or check that ",
    "'d', the full width of the confidence interval, is in the units of ",
    "the data."), fixed = TRUE)

  # When the rule is too large with the known spread constant as well, the
  # first stage is not the cause and the advice is about d alone; and the
  # distribution free rule, whose constant does not depend on the first
  # stage, never names it.
  for (m in list(msg(ts_fwd_norm_mean(c(1, 30), d = 1e-4, conf_level = 0.9999)),
                 msg(ts_fwd_nexp_mean(c(1, 30), d = 1e-9, conf_level = 0.9999)),
                 msg(ts_fwd_dfr_mean(c(1, 30), d = 1e-3, conf_level = 0.9999)))) {
    expect_match(m, "is too large to represent. Check that 'd', the (half|full) width")
    expect_match(m, "calls for more than two billion observations")
    expect_false(grepl("first stage observations the", m))
  }
})

test_that("the two stage fixed width rules do not depend on the units of measurement", {
  # Each rule depends on the data only through the first stage spread in
  # units of d, so data and d expressed in other units give the same N. The
  # values are the hand computed ones of the criterion tests above: 47
  # (Stein, 2.262^2 * (55 / 6) = 46.9), 36 (normal quantile, 35.2), and 34
  # (negative exponential, 4.459 * 3.75 / 0.5 = 33.4). Computing S^2 and d^2
  # separately overflowed beyond a scale of 1e154 and underflowed below
  # 1e-162, which stopped the first two functions with "too large to
  # represent" although N is 47 and 36; at 1e-161, where the variance is a
  # subnormal number, they returned one too many (48 for these data with
  # Stein's rule, and 20 for the second sample below with the normal
  # quantile).
  for (s in c(1e-300, 1e-200, 1e-161, 1e-100, 1e-9, 1, 1e6, 1e100, 1e154,
              1e200, 1e300)) {
    expect_identical(attr(ts_fwd_norm_mean((1:10) * s, d = s), "criterion_n"),
                     47L)
    expect_identical(attr(ts_fwd_dfr_mean((1:10) * s, d = s), "criterion_n"),
                     36L)
    expect_identical(attr(ts_fwd_nexp_mean(c(3, 4, 5, 7, 11) * s, d = 0.5 * s),
                          "criterion_n"), 34L)
  }
  # The interval is reported in the units supplied.
  big <- ts_fwd_norm_mean((1:10) * 1e200, d = 1e200)
  expect_equal(big$value[big$term == "mean"], 5.5e200)
  expect_equal(big$value[big$term == "lower_limit"], 4.5e200)
  expect_equal(big$value[big$term == "upper_limit"], 6.5e200)
  expect_equal(attr(big, "width"), 2e200)
  # A second sample, x = (2, 4, 4, 6, 9, 5, 7, 3, 8, 6): variance 44.4 / 9;
  # 2.262^2 * 4.9333 = 25.2, so 26; 1.959964^2 * 4.9333 = 18.95, so 19.
  x <- c(2, 4, 4, 6, 9, 5, 7, 3, 8, 6)
  for (s in c(1e-250, 1e-161, 1, 1e155, 1e250)) {
    expect_identical(attr(ts_fwd_norm_mean(x * s, d = s), "criterion_n"), 26L)
    expect_identical(attr(ts_fwd_dfr_mean(x * s, d = s), "criterion_n"), 19L)
  }
  # Up to the largest numbers there are. At a scale of 1.5e307 the largest
  # observation is 1.65e308 (the largest double is 1.8e308) and the scale
  # estimate is 5.6e307; multiplying it by the F quantile before dividing
  # by d overflowed, which stopped the negative exponential rule with "too
  # large to represent" although N is 34.
  expect_identical(attr(ts_fwd_nexp_mean(c(3, 4, 5, 7, 11) * 1.5e307,
                                         d = 0.5 * 1.5e307), "criterion_n"),
                   34L)
  expect_identical(attr(ts_fwd_norm_mean((1:10) * 1.5e307, d = 1.5e307),
                        "criterion_n"), 47L)
})

test_that("a row that cannot be represented is an error, never Inf or NaN", {
  # The smallest observation is -1.7e308 and d is 1e308, so the lower limit
  # is below the most negative double; N itself is 3 (the scale estimate is
  # 4e307, 0.4 in units of d, and 6.944 * 0.4 = 2.8). The function used to
  # stop with a message about the sample size.
  expect_error(ts_fwd_nexp_mean(c(-1.7e308, -1.6e308, -1e308), d = 1e308),
               "The 'lower_limit' row of the result cannot be represented (it is computed as -Inf): the data or the targets are too large in magnitude. Express the data in larger units (divide them by a constant) and state 'd' in the same units.",
               fixed = TRUE)
  # A mean plus d beyond the largest double (N is 3 here too: the variance
  # in units of d is 0.01 and t_{2, 0.025}^2 = 18.5). Which row fails first
  # depends on how the platform sums, so only the form is pinned.
  for (f in list(ts_fwd_norm_mean, ts_fwd_dfr_mean)) {
    expect_error(f(c(1.7e308, 1.6e308, 1.5e308), d = 1e308),
                 "row of the result cannot be represented \\(it is computed as (Inf|NaN)\\).* state 'd' in the same units")
  }
  # Large numbers that can be represented are returned.
  ok <- ts_fwd_nexp_mean(c(-1.7e307, -1.6e307, -1e307), d = 1e307)
  expect_true(all(is.finite(ok$value)))
  expect_equal(ok$value[ok$term == "lower_limit"], -2.7e307)
})

test_that("a constant first stage is an error that says what to do", {
  # A two stage rule fixes N from the first stage spread alone; with zero
  # spread no second stage size exists and more second stage data cannot
  # change that. (Before 1.0.0 these returned 'satisfied' at the pilot size.)
  # It is an error, not the "not met" with variance_unusable = TRUE that the
  # purely sequential checks return: a loop on "not met" would never end.
  expect_error(ts_fwd_norm_mean(rep(5, 4), d = 0.5),
               "first stage variance is zero")
  expect_error(ts_fwd_dfr_mean(rep(5, 4), d = 0.5),
               "first stage variance is zero")
  expect_error(ts_fwd_nexp_mean(rep(5, 4), d = 0.5),
               "first stage scale estimate is zero")
  expect_error(ts_fwd_norm_mean(rep(5, 4), c(1, 9), d = 0.5),
               "Enlarge the first stage")
  # Constant apart from removed missing values is still constant.
  expect_error(ts_fwd_norm_mean(c(5, 5, NA, 5), d = 0.5), "All 3 first stage")
  # Real variation, however small next to the size of the data, is enough.
  expect_s3_class(ts_fwd_norm_mean(c(5, 5, 5, 5.001), d = 0.5), "smsd_check")
})

test_that("a first stage that is constant apart from rounding error is the same error", {
  # Ten respondents each scored 0.1, 0.2, and 0.3 on three items, entered in
  # different orders; the item mean is 0.2 for everyone, but floating point
  # addition gives two different doubles one unit in the last place apart.
  # The variance of such scores (about 1e-33) is cancellation noise; plugged
  # into the rule it gave N = m and "stopping rule met" at the first stage.
  items <- c(0.1, 0.2, 0.3)
  orders <- list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2),
                 c(3, 2, 1), c(1, 2, 3), c(3, 2, 1), c(2, 1, 3), c(1, 3, 2))
  score <- vapply(orders, function(o) {
    (items[o[1]] + items[o[2]] + items[o[3]]) / 3
  }, numeric(1))
  expect_gt(length(unique(score)), 1)              # not exactly constant
  expect_lt(max(score) - min(score), 1e-15)        # but only by rounding
  expect_error(ts_fwd_norm_mean(score, d = 0.001),
               "are equal up to rounding error")
  expect_error(ts_fwd_dfr_mean(score, d = 0.001),
               "are equal up to rounding error")
  expect_error(ts_fwd_nexp_mean(score, d = 0.001),
               "are equal up to rounding error")
  expect_error(ts_br_norm_mean(score, A = 1, k = 2, w = 1e-9),
               "are equal up to rounding error")
  expect_error(ts_fwd_norm_mean(score, d = 0.001),
               "first stage variance cannot be estimated")
  # For the negative exponential rule the message names the scale (it read
  # "the first stage scale estimate cannot be estimated").
  expect_error(ts_fwd_nexp_mean(score, d = 0.001),
               "\\), so the scale cannot be estimated from them and the total sample size cannot be determined\\.")
  expect_false(grepl("estimate cannot be estimated",
                     tryCatch(ts_fwd_nexp_mean(score, d = 0.001),
                              error = conditionMessage)))
  expect_error(ts_fwd_dfr_mean(score, d = 0.001),
               "\\), so the first stage variance cannot be estimated from them")
  expect_error(ts_fwd_norm_mean(c(0.1 + 0.2, 0.3, 0.3, 0.3), d = 0.01),
               "Enlarge the first stage")
  # The same at any magnitude, and for negative data.
  expect_error(ts_fwd_norm_mean(score * 1e9, d = 0.001),
               "are equal up to rounding error")
  expect_error(ts_fwd_norm_mean(-score * 1e-9, d = 0.001),
               "are equal up to rounding error")

  # The rule is the package's one rule for constant data (.is_constant()):
  # a range of no more than 1e-10 times the largest absolute value.
  # Agreement to eleven significant digits is treated as constant; agreement
  # to nine is data. (The two stage functions used 1e-12 before the rule was
  # made the same everywhere, so the 1e-11 cases below were accepted.)
  expect_error(ts_fwd_norm_mean(c(1, 1, 1, 1 + 1e-13), d = 0.5),
               "are equal up to rounding error")
  expect_error(ts_fwd_norm_mean(c(1, 1, 1, 1 + 1e-11), d = 0.5),
               "are equal up to rounding error")
  expect_error(ts_fwd_dfr_mean(1e6 + c(0, 0, 1e-7, 0), d = 0.5),
               "are equal up to rounding error")
  expect_error(ts_fwd_dfr_mean(1e6 + c(0, 0, 1e-5, 0), d = 0.5),
               "are equal up to rounding error")
  expect_error(ts_fwd_nexp_mean(c(1, 1, 1, 1 + 1e-11), d = 0.5),
               "are equal up to rounding error")
  expect_error(ts_br_norm_mean(c(1, 1, 1, 1, 1 + 1e-11), A = 1, k = 2, w = 1),
               "are equal up to rounding error")
  # The message states the rule with the package's tolerance.
  expect_error(ts_fwd_norm_mean(c(1, 1, 1, 1 + 1e-11), d = 0.5),
               "is no more than 1e-10 times their largest absolute value")
  expect_s3_class(ts_fwd_norm_mean(c(1, 1, 1, 1 + 1e-9), d = 0.5),
                  "smsd_check")
  expect_s3_class(ts_fwd_dfr_mean(1e6 + c(0, 0, 1e-3, 0), d = 0.5),
                  "smsd_check")
  expect_s3_class(ts_fwd_nexp_mean(c(1, 1, 1, 1 + 1e-9), d = 0.5),
                  "smsd_check")
  expect_s3_class(ts_br_norm_mean(c(1, 1, 1, 1, 1 + 1e-9), A = 1, k = 2,
                                  w = 1), "smsd_check")
  # The two stage functions and the shared rule agree case by case.
  for (r in c(1e-15, 1e-13, 1e-12, 1e-11, 3e-11, 3e-10, 1e-9, 1e-6)) {
    x <- c(2, 2, 2, 2 * (1 + r))
    refused <- grepl("are equal up to rounding error",
                     tryCatch({
                       ts_fwd_norm_mean(x, d = 0.5)
                       "accepted"
                     }, error = conditionMessage))
    expect_identical(refused, .is_constant(x))
    expect_identical(refused, r < 1e-10)
  }
  # Real variation that small beside the size of the data is refused too
  # (times in seconds since 1970 that differ by hundredths of a second have
  # a relative range of 5e-11), and the remedy the message gives works:
  # subtract a common value. x = (0.01, 0.04, 0.02, 0.09) has variance
  # 0.0038 / 3; t_{3, 0.025}^2 = 10.128; with d = 0.05,
  # 10.128 * 0.0012667 / 0.0025 = 5.13, so N = 6.
  stamps <- 1.7e9 + c(0.01, 0.04, 0.02, 0.09)
  expect_error(ts_fwd_norm_mean(stamps, d = 0.05),
               "If the variation is real, subtract a common value")
  expect_identical(attr(ts_fwd_norm_mean(stamps - 1.7e9, d = 0.05),
                        "criterion_n"), 6L)
  # It is relative, so small measurements with real variation pass, and
  # their N is what the rule gives: x = (1, 2, 3, 4) * 1e-9 has variance
  # (5 / 3) * 1e-18; t_{3, 0.025} = 3.182446; with d = 1e-9,
  # 3.182446^2 * 5 / 3 = 16.88, so N = 17.
  small <- ts_fwd_norm_mean(c(1, 2, 3, 4) * 1e-9, d = 1e-9)
  expect_equal(attr(small, "criterion_n"), 17)
  # The bounded risk procedure for an exponential mean does not depend on a
  # spread and accepts such a first stage.
  expect_s3_class(ts_br_exp_mean(score, A = 1, w = 1), "smsd_check")
})
