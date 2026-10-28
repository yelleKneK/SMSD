# Tests for the minimum risk point estimation (seq_mr_*) family after the
# migration to the modern single-call API (separate _pilot, tidy
# smsd_check return, conf_level = NULL since these are risk based).

test_that("seq_mr pilot functions return positive integers", {
  expect_gte(seq_mr_cv_pilot(A = 400, c = 1), 4)
  expect_gte(seq_mr_smd_pilot(A = 100, c = 1), 4)        # equal cost: scalar
  p <- seq_mr_smd_pilot(A = 100, c = c(4, 5))            # unequal: per group
  expect_length(p, 2)
  expect_true(all(p >= 4))
})

test_that("the pilot functions take m_0 before gamma", {
  # The package wide order for a pilot function: the method's own arguments,
  # then m_0, then gamma.
  expect_identical(names(formals(seq_mr_cv_pilot)), c("A", "c", "m_0", "gamma"))
  expect_identical(names(formals(seq_mr_smd_pilot)),
                   c("A", "c", "m_0", "gamma"))
  # So a third unnamed argument is m_0: 10^(1 / 2.98) = 2.17 and
  # (10 / 2)^(1 / 2.98) = 1.72 are both below it.
  expect_identical(seq_mr_cv_pilot(10, 1, 12), 12L)
  expect_identical(seq_mr_smd_pilot(10, 1, 9), 9L)
  expect_identical(seq_mr_smd_pilot(100, c(1e6, 1), 7),
                   c(group_1 = 7L, group_2 = 7L))
  # The check functions: data, the loss constants, then gamma, with na.rm
  # last.
  expect_identical(names(formals(seq_mr_cv)),
                   c("x", "A", "c", "gamma", "variance", "na.rm"))
  expect_identical(names(formals(seq_mr_smd)),
                   c("group_1", "group_2", "A", "c", "gamma", "na.rm"))
})

test_that("every seq_mr check returns an smsd_check with no confidence level", {
  set.seed(113)
  outs <- list(
    smd  = seq_mr_smd(rnorm(40, 0.5), rnorm(40, 0), A = 100, c = 4),
    cv   = seq_mr_cv(rgamma(40, shape = 4, scale = 2), A = 400, c = 1)
  )
  for (out in outs) {
    expect_s3_class(out, "smsd_check")
    expect_true(is.logical(attr(out, "is_satisfied")))
    expect_null(attr(out, "conf_level"))     # risk based: no confidence level
    expect_true(is.numeric(out$value))
  }
})

# ---------------------------------------------------------------------------
# Minimum risk point estimation of the standardized mean difference:
# Chattopadhyay and Kelley (2017). Published numbers are pinned from the
# article (journal page given with each); DMAR::mr_smd(), which implements
# the published equal cost rule from summary statistics, is the second
# oracle; a direct numerical minimizer of the article's risk is the oracle
# for the two cost procedure.
# ---------------------------------------------------------------------------

# Two groups with a standardized mean difference of exactly `d`: pooled SD 1,
# mean difference d. Deterministic, so no seed is involved.
smd_exact <- function(n1, n2, d) {
  z1 <- stats::qnorm(stats::ppoints(n1))
  z2 <- stats::qnorm(stats::ppoints(n2))
  s_p <- sqrt((sum((z1 - mean(z1))^2) + sum((z2 - mean(z2))^2)) /
                (n1 + n2 - 2))
  list(g1 = z1 / s_p + d, g2 = z2 / s_p)
}

# The article's risk for two costs (Equation 17), minimized directly on the
# log scale. Independent of the package's root finder.
smd_risk_minimizer <- function(A, c1, c2, delta) {
  risk <- function(p) {
    n <- exp(p)
    A * (1 / n[1] + 1 / n[2] + delta^2 / (2 * sum(n))) + c1 * n[1] + c2 * n[2]
  }
  fit <- stats::optim(log(c(sqrt(A / c1), sqrt(A / c2))), risk,
                      method = "BFGS",
                      control = list(reltol = 1e-14, maxit = 2000))
  list(n = exp(fit$par), risk = fit$value)
}

val <- function(out, term) out$value[out$term == term]

test_that("seq_mr_smd divides by the pooled SD, not the pooled variance", {
  set.seed(113)
  g1 <- rnorm(50, mean = 0.5, sd = 1)
  g2 <- rnorm(50, mean = 0.0, sd = 1)
  out <- seq_mr_smd(g1, g2, A = 100, c = 4)
  smd <- out$value[out$term == "smd"]
  s_p <- sqrt(((50 - 1) * var(g1) + (50 - 1) * var(g2)) / (50 + 50 - 2))
  expect_equal(smd, (mean(g1) - mean(g2)) / s_p)
})

test_that("seq_mr_smd reproduces the article's worked example: first met at n = 75", {
  # Chattopadhyay and Kelley (2017), p. 106: with d = 1.00, A = 10,000, and
  # c = 2.40 the rule (Equation 15) is TRUE at n = 75 "and FALSE for all
  # smaller sample sizes". By hand at n = 75:
  # sqrt(10000 / 4.8) * (sqrt(2.25) + 75^-0.49) = 45.6435 * 1.62056 = 73.97,
  # so the criterion is 74; at n = 74 it is 45.6435 * 1.62136 = 74.005 > 74.
  met <- vapply(13:76, function(n) {
    dat <- smd_exact(n, n, 1)
    attr(seq_mr_smd(dat$g1, dat$g2, A = 10000, c = 2.40), "is_satisfied")
  }, logical(1))
  expect_identical(which(met), c(63L, 64L))        # n = 75 and n = 76 only
  dat <- smd_exact(75, 75, 1)
  out <- seq_mr_smd(dat$g1, dat$g2, A = 10000, c = 2.40)
  expect_equal(val(out, "smd"), 1, tolerance = 1e-12)
  expect_equal(val(out, "optimal_n1"), 74)
  expect_equal(val(out, "optimal_n2"), 74)
  expect_equal(attr(out, "criterion_n"), 148)
})

test_that("with one cost and unequal group sizes the smaller group decides", {
  # The article states Equation 15 for n1 = n2 = n. With d = 1.00,
  # A = 10,000, c = 2.40 the criterion is 75 at n = 74 and 74 at n = 75 (see
  # the test above), so 80 and 74 observations do not meet the rule, and 80
  # and 75 do, whichever group is the larger.
  short <- smd_exact(80, 74, 1)
  out <- seq_mr_smd(short$g1, short$g2, A = 10000, c = 2.40)
  expect_equal(val(out, "optimal_n1"), 75)
  expect_equal(val(out, "optimal_n2"), 75)
  expect_false(attr(out, "is_satisfied"))
  enough <- smd_exact(80, 75, 1)
  expect_true(attr(seq_mr_smd(enough$g1, enough$g2, A = 10000, c = 2.40),
                   "is_satisfied"))
  expect_true(attr(seq_mr_smd(enough$g2, enough$g1, A = 10000, c = 2.40),
                   "is_satisfied"))
})

test_that("seq_mr_smd agrees with DMAR::mr_smd, including at large effects", {
  # DMAR::mr_smd() is the published equal cost rule from (d, n). The old
  # solver in this package failed from about |d| = 2.8.
  for (d in c(0, 0.5, 1, 2, 2.8, 3, 4, -3)) {
    for (n in c(20, 60, 150)) {
      dat <- smd_exact(n, n, d)
      out <- seq_mr_smd(dat$g1, dat$g2, A = 5000, c = 1)
      ref <- DMAR::mr_smd(d = d, n = n, A = 5000, sampling_cost = 1)
      expect_equal(val(out, "optimal_n1"), ceiling(attr(ref, "criterion")))
      expect_identical(attr(out, "is_satisfied"), attr(ref, "is_satisfied"))
      expect_equal(val(out, "risk"), ref$value[ref$term == "risk"])
    }
  }
})

test_that("seq_mr_smd reproduces the article's pilot data results", {
  # Chattopadhyay and Kelley (2017), pp. 105 to 106: d = 1.021484; the rule
  # is FALSE at n = 13 and again at n = 14.
  high <- c(11, 7, 22, 13, 6, 9, 11, 16, 12, 17, 14, 8, 16)
  low <- c(3, 6, 10, 8, 14, 5, 12, 10, 6, 8, 13, 5, 9)
  out <- seq_mr_smd(high, low, A = 10000, c = 2.40)
  expect_equal(val(out, "smd"), 1.021484, tolerance = 5e-7)
  expect_false(attr(out, "is_satisfied"))
  expect_false(attr(seq_mr_smd(c(high, 10), c(low, 8), A = 10000, c = 2.40),
                    "is_satisfied"))
})

test_that("seq_mr_smd_pilot gives the article's pilot sizes", {
  # Equation 14: 13 on p. 105 (A = 10,000, c = 2.40) and 30 on p. 102
  # (A = 50,000, c = 1).
  expect_identical(seq_mr_smd_pilot(A = 10000, c = 2.40), 13L)
  expect_identical(seq_mr_smd_pilot(A = 50000, c = 1), 30L)
  # Equations 20 and 21 for A = 50,000 and costs 2 and 8, by hand:
  # (50000 / 2)^(1 / 2.98) = 29.91, so 30; then
  # 30 / sqrt(1 + 6 * 900 / 50000) = 30 / 1.05262 = 28.50, so 29.
  expect_identical(seq_mr_smd_pilot(A = 50000, c = c(2, 8)),
                   c(group_1 = 30L, group_2 = 29L))
  # The less costly group gets the larger pilot whichever comes first.
  expect_identical(seq_mr_smd_pilot(A = 50000, c = c(8, 2)),
                   c(group_1 = 29L, group_2 = 30L))
  # Two equal costs give the one cost (Equation 14) size for both groups.
  expect_identical(seq_mr_smd_pilot(A = 50000, c = c(1, 1)),
                   c(group_1 = 30L, group_2 = 30L))
  # m_0 is a floor: (10 / 2)^(1 / 2.98) = 1.72.
  expect_identical(seq_mr_smd_pilot(A = 10, c = 1), 4L)
  expect_identical(seq_mr_smd_pilot(A = 10, c = 1, m_0 = 9), 9L)
  # m_0 is also a floor for the costlier group (Equation 21). By hand:
  # 1e6^(1 / 2.98) = 103.1, so 104, and 104 / sqrt(1 + (1e9 - 1) * 104^2 /
  # 1e6) = 0.03, which rounds up to 1 and is raised to 4. Likewise
  # 100^(1 / 2.98) = 4.69, so 5, and 5 / sqrt(1 + (1e6 - 1) * 25 / 100)
  # = 0.01.
  expect_identical(seq_mr_smd_pilot(A = 1e6, c = c(1, 1e9)),
                   c(group_1 = 104L, group_2 = 4L))
  expect_identical(seq_mr_smd_pilot(A = 100, c = c(1, 1e6)),
                   c(group_1 = 5L, group_2 = 4L))
  expect_identical(seq_mr_smd_pilot(A = 100, c = c(1e6, 1), m_0 = 7),
                   c(group_1 = 7L, group_2 = 7L))
  # gamma enters the exponent. With gamma = 0.25, by hand:
  # (10000 / 4.8)^(1 / 2.5) = 21.26, so 22 (13 with the default); and
  # 25000^(1 / 2.5) = 57.43, so 58, with 58 / sqrt(1 + 6 * 58^2 / 50000)
  # = 48.95, so 49 (30 and 29 with the default).
  expect_identical(seq_mr_smd_pilot(A = 10000, c = 2.40, gamma = 0.25), 22L)
  expect_identical(seq_mr_smd_pilot(A = 50000, c = c(2, 8), gamma = 0.25),
                   c(group_1 = 58L, group_2 = 49L))
})

test_that("the stopping rules use gamma", {
  # One cost, d = 1.00, n = 60, A = 10,000, c = 2.40, gamma = 0.25. By hand:
  # sqrt(10000 / 4.8) * (1.5 + 60^-0.25) = 45.6435 * 1.85930 = 84.87, so 85
  # (75 with the default gamma: 45.6435 * (1.5 + 60^-0.49) = 74.60).
  dat <- smd_exact(60, 60, 1)
  out <- seq_mr_smd(dat$g1, dat$g2, A = 10000, c = 2.40, gamma = 0.25)
  expect_equal(val(out, "optimal_n1"), 85)
  expect_identical(attr(out, "gamma"), 0.25)
  ref <- DMAR::mr_smd(d = 1, n = 60, A = 10000, sampling_cost = 2.40,
                      gamma = 0.25)
  expect_equal(val(out, "optimal_n1"), ceiling(attr(ref, "criterion")))
  expect_equal(val(seq_mr_smd(dat$g1, dat$g2, A = 10000, c = 2.40),
                   "optimal_n1"), 75)
  # Two costs (2 and 8), A = 50,000, d = 1.00, 60 observations in the less
  # costly group, gamma = 0.25. By hand: r = 1 / sqrt(1 + 6 * 3600 / 50000)
  # = 0.83566; sqrt(1 + 1 / (2 * 1.83566^2)) = 1.07163; 60^-0.25 = 0.35930;
  # sqrt(25000) * 1.43093 = 226.25, so 227 (191 with the default gamma:
  # 158.114 * (1.07163 + 60^-0.49) = 190.70). The size paired with 227 is
  # 227 / sqrt(1 + 6 * 227^2 / 50000) = 84.70, so 85.
  dat <- smd_exact(60, 51, 1)
  out <- seq_mr_smd(dat$g1, dat$g2, A = 50000, c = c(2, 8), gamma = 0.25)
  expect_equal(val(out, "optimal_n1"), 227)
  expect_equal(val(out, "optimal_n2"), 85)
  expect_equal(val(seq_mr_smd(dat$g1, dat$g2, A = 50000, c = c(2, 8)),
                   "optimal_n1"), 191)
})

test_that("two costs equal up to rounding error are one cost", {
  # 24 * 0.1 is not 2.4 in floating point, and the article's two procedures
  # differ at equal costs: the two cost pilot would be
  # (10000 / 2.4)^(1 / 2.98) = 16.39, so 17, where the article's example
  # (p. 105) has (10000 / 4.8)^(1 / 2.98) = 12.99, so 13.
  expect_false(24 * 0.1 == 2.4)
  expect_identical(seq_mr_smd_pilot(A = 10000, c = c(24 * 0.1, 2.4)),
                   c(group_1 = 13L, group_2 = 13L))
  expect_identical(seq_mr_smd_pilot(A = 10000, c = c(2.4, 24 * 0.1)),
                   c(group_1 = 13L, group_2 = 13L))
  # (50000 / 0.6)^(1 / 2.98) = 44.80, so 45 (two costs would give 57).
  expect_identical(seq_mr_smd_pilot(A = 50000, c = c(0.1 + 0.2, 0.3)),
                   c(group_1 = 45L, group_2 = 45L))
  # The article's worked example is still first met at n = 75, and the
  # whole object is the one cost object (the smaller of the two costs,
  # here 2.4 itself, is the cost used).
  for (cost in list(c(24 * 0.1, 2.4), c(2.4, 24 * 0.1))) {
    at_74 <- smd_exact(74, 74, 1)
    at_75 <- smd_exact(75, 75, 1)
    expect_false(attr(seq_mr_smd(at_74$g1, at_74$g2, A = 10000, c = cost),
                      "is_satisfied"))
    expect_identical(seq_mr_smd(at_75$g1, at_75$g2, A = 10000, c = cost),
                     seq_mr_smd(at_75$g1, at_75$g2, A = 10000, c = 2.4))
  }
  # The tolerance is a relative difference of sqrt(.Machine$double.eps),
  # about 1.5e-8, so a real difference, however small, still selects the
  # two cost procedure. By hand: one cost, 25000^(1 / 2.98) = 29.91, so 30;
  # two costs, 50000^(1 / 2.98) = 37.74, so 38 (and the paired size is 38).
  expect_identical(seq_mr_smd_pilot(A = 50000, c = c(1, 1 + 1e-9)),
                   c(group_1 = 30L, group_2 = 30L))
  expect_identical(seq_mr_smd_pilot(A = 50000, c = c(1, 1 + 1e-7)),
                   c(group_1 = 38L, group_2 = 38L))
  # c = c(2, 2.0001): two costs, 25000^(1 / 2.98) = 29.91, so 30 and 30;
  # c = 2: 12500^(1 / 2.98) = 23.70, so 24.
  expect_identical(seq_mr_smd_pilot(A = 50000, c = c(2, 2.0001)),
                   c(group_1 = 30L, group_2 = 30L))
  expect_identical(seq_mr_smd_pilot(A = 50000, c = 2), 24L)
})

test_that("two equal costs give exactly the one cost result", {
  for (d in c(0.4, 3)) {
    dat <- smd_exact(40, 40, d)
    expect_identical(seq_mr_smd(dat$g1, dat$g2, A = 5000, c = c(3, 3)),
                     seq_mr_smd(dat$g1, dat$g2, A = 5000, c = 3))
  }
})

test_that("the two cost optimum minimizes the article's risk (Equation 17)", {
  # Includes large |d| and cost ratios up to 8:1, where the old fixed point
  # iteration failed. Equation 18 as printed (first power) would give
  # n1 = 183.4 in the first case; the minimizer is 175.66.
  first <- .smd_mr_optimum(1, 50000, 2, 8)
  expect_equal(first[1:2], smd_risk_minimizer(50000, 2, 8, 1)$n,
               tolerance = 1e-6)
  expect_equal(first[1], 175.66, tolerance = 5e-5)
  for (A in c(500, 50000)) {
    for (ratio in c(2, 4, 8)) {
      for (delta in c(0, 0.5, 2, 3, 4)) {
        ref <- smd_risk_minimizer(A, 2, 2 * ratio, delta)
        opt <- .smd_mr_optimum(delta^2, A, 2, 2 * ratio)
        expect_equal(opt[1:2], ref$n, tolerance = 1e-5)
        expect_equal(opt[3], ref$risk, tolerance = 1e-9)
        # The same optimum with the groups in the other order.
        expect_equal(.smd_mr_optimum(delta^2, A, 2 * ratio, 2)[1:2],
                     rev(opt[1:2]), tolerance = 1e-12)
      }
    }
  }
  # With equal costs it is Equation 10: sqrt(A / (2c)) sqrt(2 + delta^2 / 4),
  # and the minimum risk is 4 c n_c (Equation 12).
  n_c <- sqrt(10000 / 4.8) * sqrt(2 + 1 / 4)
  expect_equal(.smd_mr_optimum(1, 10000, 2.4, 2.4), c(n_c, n_c, 4 * 2.4 * n_c))
})

test_that("the two cost rule stops just past the risk minimizing size", {
  # With d held fixed, the article's rule without its correction term is
  # first met at the risk minimizing n1 (the optimum is a fixed point), and
  # the correction term adds at most sqrt(A / c1) * n1^-gamma. So the first
  # n1 at which the check is met lies strictly above the minimizer's n1 and
  # no more than the correction (plus rounding) beyond it. Equation 18 as
  # printed fails this (it would stop near 196 in the first case below).
  for (cfg in list(c(50000, 2, 8, 1), c(50000, 2, 8, 3), c(5000, 1, 4, 2),
                   c(5000, 1, 2, 4))) {
    A <- cfg[1]; c1 <- cfg[2]; c2 <- cfg[3]; d <- cfg[4]
    ref <- smd_risk_minimizer(A, c1, c2, d)$n
    upper <- ref[1] + sqrt(A / c1) * ref[1]^(-0.49) + 2
    n1 <- floor(ref[1]) - 2
    repeat {
      n2 <- ceiling(n1 / sqrt(1 + (c2 - c1) * n1^2 / A))
      dat <- smd_exact(n1, n2, d)
      out <- seq_mr_smd(dat$g1, dat$g2, A = A, c = c(c1, c2))
      # The search ends just past the bound checked below, so that a rule
      # that is never met fails that check and does not run on.
      if (attr(out, "is_satisfied") || n1 > upper + 1) break
      n1 <- n1 + 1
    }
    expect_gt(n1, ref[1])
    expect_lte(n1, upper)
  }
})

test_that("with two costs the costlier group must reach its paired size", {
  # A = 50,000, costs 2 and 8, d = 0, 200 observations in the less costly
  # group. By hand: criterion = sqrt(25000) * (1 + 200^-0.49)
  # = 158.114 * 1.07456 = 169.9, so 170 <= 200; the paired size is
  # 200 / sqrt(1 + 6 * 200^2 / 50000) = 200 / sqrt(5.8) = 83.05, so 84.
  dat <- smd_exact(200, 84, 0)
  out <- seq_mr_smd(dat$g1, dat$g2, A = 50000, c = c(2, 8))
  expect_equal(val(out, "optimal_n1"), 170)
  expect_equal(val(out, "optimal_n2"), 84)
  expect_true(attr(out, "is_satisfied"))
  short <- seq_mr_smd(dat$g1, dat$g2[-1], A = 50000, c = c(2, 8))
  expect_equal(val(short, "optimal_n2"), 84)
  expect_false(attr(short, "is_satisfied"))
  # Before the less costly group reaches its criterion, the costlier
  # group's criterion is the size paired with that criterion:
  # sqrt(25000) * (1 + 100^-0.49) = 158.114 * 1.10471 = 174.7, so 175, and
  # 175 / sqrt(1 + 6 * 175^2 / 50000) = 175 / sqrt(4.675) = 80.94, so 81.
  early <- smd_exact(100, 60, 0)
  out <- seq_mr_smd(early$g1, early$g2, A = 50000, c = c(2, 8))
  expect_equal(val(out, "optimal_n1"), 175)
  expect_equal(val(out, "optimal_n2"), 81)
  expect_false(attr(out, "is_satisfied"))
  # The size is paired with the criterion as reported, that is, rounded up,
  # because that is the size the less costly group has to reach. With 90
  # observations: sqrt(25000) * (1 + 90^-0.49) = 158.114 * 1.11026 = 175.55,
  # so 176, and 176 / sqrt(1 + 6 * 176^2 / 50000) = 176 / sqrt(4.71712)
  # = 81.04, so 82. (Paired with the unrounded 175.55 it would be 80.99, so
  # 81: a costlier group of 81 would then meet its reported criterion while
  # the rule, once the less costly group reaches 176, still asks for 82.)
  dat <- smd_exact(90, 60, 0)
  out <- seq_mr_smd(dat$g1, dat$g2, A = 50000, c = c(2, 8))
  expect_equal(val(out, "optimal_n1"), 176)
  expect_equal(val(out, "optimal_n2"), 82)
  expect_equal(val(seq_mr_smd(dat$g2, dat$g1, A = 50000, c = c(8, 2)),
                   "optimal_n1"), 82)
  at_176 <- smd_exact(176, 81, 0)
  expect_false(attr(seq_mr_smd(at_176$g1, at_176$g2, A = 50000, c = c(2, 8)),
                    "is_satisfied"))
  at_176 <- smd_exact(176, 82, 0)
  expect_true(attr(seq_mr_smd(at_176$g1, at_176$g2, A = 50000, c = c(2, 8)),
                   "is_satisfied"))
})

test_that("the order of the groups does not matter with two costs", {
  set.seed(113)
  g1 <- rnorm(120, 0.7); g2 <- rnorm(70, 0)
  a <- seq_mr_smd(g1, g2, A = 20000, c = c(2, 8))
  b <- seq_mr_smd(g2, g1, A = 20000, c = c(8, 2))
  expect_equal(val(a, "smd"), -val(b, "smd"))
  expect_equal(val(a, "risk"), val(b, "risk"))
  expect_equal(val(a, "optimal_n1"), val(b, "optimal_n2"))
  expect_equal(val(a, "optimal_n2"), val(b, "optimal_n1"))
  expect_identical(attr(a, "is_satisfied"), attr(b, "is_satisfied"))
  expect_gt(val(a, "optimal_n1"), val(a, "optimal_n2"))   # cheaper: more
})

test_that("next_n1 and next_n2 give the sizes for the next check", {
  next_n <- function(out) c(val(out, "next_n1"), val(out, "next_n2"))
  # One cost, d = 1.00, A = 10,000, c = 2.40 (the article's example: not met
  # at 74 per group, met at 75). Not met: one more observation per group.
  dat <- smd_exact(74, 74, 1)
  expect_identical(next_n(seq_mr_smd(dat$g1, dat$g2, A = 10000, c = 2.40)),
                   c(75, 75))
  # Met: nothing more is called for.
  dat <- smd_exact(75, 75, 1)
  expect_identical(next_n(seq_mr_smd(dat$g1, dat$g2, A = 10000, c = 2.40)),
                   c(75, 75))
  # Unequal sizes: the smaller group is n, so it alone is brought up.
  dat <- smd_exact(80, 74, 1)
  expect_identical(next_n(seq_mr_smd(dat$g1, dat$g2, A = 10000, c = 2.40)),
                   c(80, 75))
  expect_identical(next_n(seq_mr_smd(dat$g2, dat$g1, A = 10000, c = 2.40)),
                   c(75, 80))

  # Two costs (2 and 8), A = 50,000, d = 0. The article pairs n observations
  # in the less costly group with ceiling(n / sqrt(1 + 6 n^2 / 50000)) in
  # the costlier one (Equation 21 and Stage II, p. 107). By hand:
  # 100 / sqrt(2.2) = 67.42, so 68; 101 / sqrt(2.22412) = 67.72, so 68;
  # 200 / sqrt(5.8) = 83.05, so 84.
  # The costlier group short of its paired size: that shortfall alone.
  dat <- smd_exact(100, 60, 0)
  expect_identical(next_n(seq_mr_smd(dat$g1, dat$g2, A = 50000, c = c(2, 8))),
                   c(100, 68))
  # Paired, and the rule not met (the criterion is 175): one more in the
  # less costly group, and the size paired with 101 in the costlier one.
  dat <- smd_exact(100, 68, 0)
  out <- seq_mr_smd(dat$g1, dat$g2, A = 50000, c = c(2, 8))
  expect_false(attr(out, "is_satisfied"))
  expect_identical(next_n(out), c(101, 68))
  # A costlier group already past the paired size is left as it is.
  dat <- smd_exact(100, 90, 0)
  expect_identical(next_n(seq_mr_smd(dat$g1, dat$g2, A = 50000, c = c(2, 8))),
                   c(101, 90))
  # The less costly group past its criterion (170) and the costlier group
  # one short of its paired size (84): only the costlier group grows.
  dat <- smd_exact(200, 83, 0)
  out <- seq_mr_smd(dat$g1, dat$g2, A = 50000, c = c(2, 8))
  expect_false(attr(out, "is_satisfied"))
  expect_identical(next_n(out), c(200, 84))
  # Met: the current sizes.
  dat <- smd_exact(200, 84, 0)
  out <- seq_mr_smd(dat$g1, dat$g2, A = 50000, c = c(2, 8))
  expect_true(attr(out, "is_satisfied"))
  expect_identical(next_n(out), c(200, 84))
  # The rows are in the order of the arguments, whichever group is the less
  # costly one.
  dat <- smd_exact(60, 100, 0)
  expect_identical(next_n(seq_mr_smd(dat$g1, dat$g2, A = 50000, c = c(8, 2))),
                   c(68, 100))
  dat <- smd_exact(68, 100, 0)
  expect_identical(next_n(seq_mr_smd(dat$g1, dat$g2, A = 50000, c = c(8, 2))),
                   c(68, 101))

  # The sizes do not depend on the estimate, so constant groups get them
  # too. One cost: 11 and 11. Costs 1 and 2 with A = 100: 10 / sqrt(2)
  # = 7.07, so 8, which the 12 observations in hand exceed; then one more
  # in the less costly group, and 11 / sqrt(2.21) = 7.40, so 8, again below
  # 12.
  out <- seq_mr_smd(rep(1, 10), rep(2, 10), A = 100, c = 1)
  expect_true(is.na(val(out, "smd")))
  expect_identical(next_n(out), c(11, 11))
  out <- seq_mr_smd(rep(0.3, 10), rep(0.1 + 0.2, 12), A = 100, c = c(1, 2))
  expect_identical(next_n(out), c(11, 12))
})

test_that("following next_n1 and next_n2 from the pilot is the article's sampling order", {
  # From the pilot sizes on, each step must add one observation to the less
  # costly group and hold the costlier group at
  # max(its size, ceiling(n / sqrt(1 + (c_hi - c_lo) n^2 / A))), with n the
  # less costly group's new size (Stage II, p. 107, with one observation per
  # stage). The formula is written out here, not taken from the package.
  # A = 2,000: the criterion for the less costly group is at least
  # sqrt(2000 / 3) = 25.8, and the walks below stop when it is reached or
  # after 25 steps.
  for (cost in list(c(3, 8), c(8, 3), c(3, 3.5))) {
    lo <- which.min(cost)
    hi <- which.max(cost)
    pilot <- seq_mr_smd_pilot(A = 2000, c = cost)
    size <- unname(pilot)
    steps <- 0
    repeat {
      dat <- smd_exact(size[1], size[2], 0.5)
      out <- seq_mr_smd(dat$g1, dat$g2, A = 2000, c = cost)
      if (attr(out, "is_satisfied") || steps == 25) break
      n_lo <- size[lo] + 1
      by_hand <- size
      by_hand[lo] <- n_lo
      by_hand[hi] <- max(size[hi],
                         ceiling(n_lo / sqrt(1 + (cost[hi] - cost[lo]) *
                                               n_lo^2 / 2000)))
      size <- c(val(out, "next_n1"), val(out, "next_n2"))
      expect_identical(size, by_hand)
      steps <- steps + 1
    }
    expect_gt(steps, 10)
  }
})

test_that("the decision and the reported criterion agree at a whole number", {
  # The criterion is rounded up by the shared .as_sample_size(), whose
  # ceiling ignores representation error, and the decision uses the same
  # rounding, so "met" is exactly "each current size is at least its
  # criterion". With d = 0 and one cost the criterion is
  # sqrt(A / (2c)) (sqrt(2) + n^-gamma); A0 makes it exactly 40 at n = 40.
  n <- 40
  A0 <- 2 * (n / (sqrt(2) + n^-0.49))^2
  dat <- smd_exact(n, n, 0)
  # A hair above 40 (by 2e-12, representation error): the criterion is 40
  # and the rule is met, not "not met" beside a reported criterion of 40.
  out <- seq_mr_smd(dat$g1, dat$g2, A = A0 * (1 + 1e-13), c = 1)
  expect_gt(sqrt(A0 * (1 + 1e-13) / 2) * (sqrt(2) + n^-0.49), 40)
  expect_equal(val(out, "optimal_n1"), 40)
  expect_true(attr(out, "is_satisfied"))
  # Really above 40 (by 2e-5): the criterion is 41 and the rule is not met.
  out <- seq_mr_smd(dat$g1, dat$g2, A = A0 * (1 + 1e-6), c = 1)
  expect_equal(val(out, "optimal_n1"), 41)
  expect_false(attr(out, "is_satisfied"))
  # The same with two costs (1 and 3) and d = 0: the less costly group's
  # criterion is sqrt(A / 1) (1 + n^-gamma), and B0 makes it exactly 40 at
  # n = 40. B0 = (40 / 1.16405)^2 = 1180.8, so the size paired with 40 is
  # 40 / sqrt(1 + 2 * 1600 / 1180.8) = 40 / sqrt(3.710) = 20.77, so 21 (and
  # the size paired with 41 is 41 / sqrt(3.847) = 20.90, also 21).
  B0 <- (n / (1 + n^-0.49))^2
  expect_equal(B0, 1180.8, tolerance = 1e-5)
  dat <- smd_exact(n, 21, 0)
  for (flip in c(FALSE, TRUE)) {
    check <- function(A) {
      if (flip) {
        seq_mr_smd(dat$g2, dat$g1, A = A, c = c(3, 1))
      } else {
        seq_mr_smd(dat$g1, dat$g2, A = A, c = c(1, 3))
      }
    }
    lo <- if (flip) "optimal_n2" else "optimal_n1"
    hi <- if (flip) "optimal_n1" else "optimal_n2"
    # A hair above 40 (by 2e-12): the criterion is 40 and the rule is met.
    out <- check(B0 * (1 + 1e-13))
    expect_gt(sqrt(B0 * (1 + 1e-13)) * (1 + n^-0.49), 40)
    expect_equal(val(out, lo), 40)
    expect_equal(val(out, hi), 21)
    expect_true(attr(out, "is_satisfied"))
    # Really above 40 (by 2e-5): the criterion is 41 and the rule is not met.
    out <- check(B0 * (1 + 1e-6))
    expect_equal(val(out, lo), 41)
    expect_equal(val(out, hi), 21)
    expect_false(attr(out, "is_satisfied"))
  }
  # The helper: rounds as .as_sample_size() does, and passes on what cannot
  # be represented without failing.
  expect_identical(.mrpe_ceiling(56.00000000000001), 56)
  expect_identical(.mrpe_ceiling(56.001), 57)
  expect_identical(.mrpe_ceiling(56), 56)
  expect_identical(.mrpe_ceiling(3e9), 3e9)
  expect_identical(.mrpe_ceiling(Inf), Inf)
  # In every result that has an estimate, is_satisfied is "each current size
  # is at least its criterion".
  withr::local_seed(113)
  for (i in 1:40) {
    n1 <- sample(8:60, 1)
    n2 <- sample(8:60, 1)
    cost <- list(1, c(1, 3), c(4, 2))[[sample(3, 1)]]
    out <- seq_mr_smd(rnorm(n1, 0.5), rnorm(n2), A = 1500, c = cost)
    expect_identical(attr(out, "is_satisfied"),
                     n1 >= val(out, "optimal_n1") &&
                       n2 >= val(out, "optimal_n2"))
  }
})

test_that("seq_mr_smd returns the documented rows and attributes", {
  set.seed(113)
  out <- seq_mr_smd(rnorm(30, 0.5), rnorm(30), A = 2000, c = c(2, 3))
  expect_identical(out$term, c("smd", "risk", "mean_1", "mean_2",
                               "optimal_n1", "optimal_n2",
                               "current_n1", "current_n2",
                               "next_n1", "next_n2"))
  expect_true(all(is.finite(out$value)))
  expect_identical(attr(out, "method"),
                   "minimum risk point estimation (standardized mean difference)")
  expect_identical(attr(out, "current_n"), 60L)
  expect_equal(attr(out, "criterion_n"),
               val(out, "optimal_n1") + val(out, "optimal_n2"))
  # A double on every path, as in the rest of the package.
  expect_type(attr(out, "criterion_n"), "double")
  expect_type(attr(seq_mr_smd(rep(1, 10), rep(1, 10), A = 100, c = 1),
                   "criterion_n"), "double")
  # The total is formed as a double: each criterion below is about 1.64e9
  # (sqrt(A / (2c)) = 1e9 times sqrt(2) + 20^-0.49 = 1.6446 at d = 0), which
  # is representable as an integer, while their sum is not.
  big <- smd_exact(20, 20, 0)
  expect_warning(out_big <- seq_mr_smd(big$g1, big$g2, A = 2e18, c = 1), NA)
  expect_equal(val(out_big, "optimal_n1"), 1.6446e9, tolerance = 1e-4)
  expect_gt(attr(out_big, "criterion_n"), .Machine$integer.max)
  expect_equal(attr(out_big, "criterion_n"), 2 * val(out_big, "optimal_n1"))
  expect_identical(attr(out, "loss_A"), 2000)
  expect_identical(c(attr(out, "cost_1"), attr(out, "cost_2")), c(2, 3))
  expect_identical(attr(out, "gamma"), 0.49)
  expect_false(attr(out, "variance_unusable"))
  expect_null(attr(out, "conf_level"))
  expect_identical(attr(out, "notes"), character(0))
  # The risk row is Equation 17 at the current sizes with d for delta.
  d <- val(out, "smd")
  expect_equal(val(out, "risk"),
               2000 * (1 / 30 + 1 / 30 + d^2 / 120) + 2 * 30 + 3 * 30)
})

test_that("constant groups are not estimable and never reported as met", {
  for (out in list(seq_mr_smd(rep(1, 10), rep(1, 10), A = 100, c = 1),
                   seq_mr_smd(rep(1, 10), rep(2, 10), A = 100, c = 1),
                   seq_mr_smd(rep(0.3, 10), rep(0.1 + 0.2, 12), A = 100,
                              c = c(1, 2)))) {
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_true(is.na(val(out, "smd")))
    expect_true(is.na(val(out, "risk")))
    expect_true(is.na(val(out, "optimal_n1")))
    expect_false(any(is.infinite(out$value) | is.nan(out$value)))
    expect_match(attr(out, "notes"), "not finite", fixed = TRUE)
    expect_true(is.na(attr(out, "criterion_n")))
    expect_false(attr(out, "variance_unusable"))
    expect_output(print(out), "Stopping rule met: FALSE", fixed = TRUE)
  }
  # One constant group is fine: the pooled SD is positive.
  set.seed(113)
  out <- seq_mr_smd(rep(1, 10), rnorm(10), A = 100, c = 1)
  expect_true(is.finite(val(out, "smd")))
  # Unless the group that varies does so on a scale so far below the other
  # group's values that its squares underflow beside them: the pooled SD is
  # then zero in floating point, and the result is "no estimate", not an
  # infinite one.
  out <- seq_mr_smd(c(1, 2, 3, 4) * 1e-200, rep(1, 4), A = 100, c = 1)
  expect_identical(attr(out, "is_satisfied"), FALSE)
  expect_true(is.na(val(out, "smd")))
  expect_false(any(is.infinite(out$value) | is.nan(out$value)))
  expect_identical(attr(out, "notes"), .note_not_estimable())
})

test_that("groups that vary only by rounding error are not estimable", {
  # The pooled SD here is positive (about 1.6e-15) but is rounding error in
  # values near 1 and 2, so a ratio would be noise of the order of 1e15.
  # The groups are not literally constant, so the note is the one that also
  # covers data with a large common offset (see the next test).
  g1 <- 1 + 1e-15 * c(-1, 1, 0, 2, -2)
  g2 <- 2 + 1e-15 * c(1, -1, 0, 2, -2)
  expect_gt(stats::sd(g1), 0)
  expect_gt(stats::sd(g2), 0)
  for (cost in list(1, c(1, 2))) {
    out <- seq_mr_smd(g1, g2, A = 100, c = cost)
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_true(is.na(val(out, "smd")))
    expect_true(is.na(val(out, "optimal_n1")))
    expect_identical(attr(out, "notes"), .note_smd_within_tolerance())
    expect_match(attr(out, "notes"), "agree to ten significant digits",
                 fixed = TRUE)
  }
  # "Constant" is the package wide test .is_constant(): a range of at most
  # 1e-10 times the largest absolute value. Both groups below have a
  # relative range of 4e-12, so there is no estimate.
  out <- seq_mr_smd(1 + 1e-12 * c(-1, 1, 0, 2, -2),
                    2 + 1e-12 * c(1, -1, 0, 2, -2), A = 100, c = 1)
  expect_identical(attr(out, "is_satisfied"), FALSE)
  expect_true(is.na(val(out, "smd")))
  expect_identical(attr(out, "notes"), .note_smd_within_tolerance())
  # One group literally constant and the other within the tolerance: the
  # same note.
  out <- seq_mr_smd(rep(1, 5), g2, A = 100, c = 1)
  expect_true(is.na(val(out, "smd")))
  expect_identical(attr(out, "notes"), .note_smd_within_tolerance())
  # A relative range of 4e-9 is ordinary data: the estimate is the mean
  # difference over the pooled SD, here (1e-9 * 0.2) / (1e-9 * 1.7607).
  h1 <- 1 + 1e-9 * c(-1, 1, 0, 2, -2)
  h2 <- 1 + 1e-9 * c(1, -1, 0, 2, -3)
  out <- seq_mr_smd(h1, h2, A = 100, c = 1)
  expect_equal(val(out, "smd"), 0.2 / sqrt((2.5 + 3.7) / 2), tolerance = 1e-6)
  # One group constant up to rounding error and one that varies: estimable.
  out <- seq_mr_smd(g1, c(1.5, 2, 2.5, 3, 1), A = 100, c = 1)
  expect_true(is.finite(val(out, "smd")))
})

test_that("data with a large common offset get a note that says what to do", {
  # The test for "constant" is relative to the size of the values, so two
  # groups that vary by about 1 around 1e11 (a relative range of 1.5e-11)
  # are constant to ten significant digits and there is no estimate. The
  # standardized mean difference does not depend on a common offset, so the
  # note says to subtract one; collecting more data would not help. The
  # values are multiples of 1/4, which are exact at this magnitude.
  g1 <- c(1.25, 0.5, 2, 0.75, 1.5, 1)
  g2 <- c(0.25, 1, 0.5, -0.25, 0.75, 0)
  expect_identical((g1 + 1e11) - 1e11, g1)
  for (cost in list(1, c(1, 2))) {
    out <- seq_mr_smd(g1 + 1e11, g2 + 1e11, A = 100, c = cost)
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_true(is.na(val(out, "smd")))
    expect_equal(val(out, "mean_1"), 1e11 + 7 / 6)
    expect_identical(attr(out, "notes"), .note_smd_within_tolerance())
    expect_match(attr(out, "notes"),
                 "subtract the same constant from both groups", fixed = TRUE)
    expect_output(print(out), "large common", fixed = TRUE)
  }
  # After the subtraction the estimate is that of the data without the
  # offset. By hand: the means are 7/6 and 3/8, the sums of squared
  # deviations 35/24 and 35/32, so the pooled variance is
  # (35/24 + 35/32) / 10 = 0.2552083 and d = (19/24) / sqrt(0.2552083)
  # = 1.567094.
  out <- seq_mr_smd((g1 + 1e11) - 1e11, (g2 + 1e11) - 1e11, A = 100, c = 1)
  expect_equal(val(out, "smd"), (19 / 24) / sqrt((35 / 24 + 35 / 32) / 10))
  expect_equal(val(out, "smd"), 1.567094, tolerance = 1e-6)
  expect_identical(attr(out, "notes"), character(0))
  # An offset that leaves more than ten significant digits of variation
  # needs no such step: at 1e6 the estimate is the same.
  out <- seq_mr_smd(g1 + 1e6, g2 + 1e6, A = 100, c = 1)
  expect_equal(val(out, "smd"), 1.567094, tolerance = 1e-6)
})

test_that("seq_mr_smd does not depend on the common units of the two groups", {
  # Two groups of 8 with means 4.5 and 3 and sums of squared deviations 18
  # and 12: the pooled variance is 30 / 14 and d = 1.5 / sqrt(15 / 7)
  # = 1.024695. Multiplying both groups by a power of two is exact, so every
  # row but the two means, and every attribute, must be identical, and the
  # means must be the scaled means. Computed on the data as given, the sums
  # of squares overflow from about 1e154 on (d comes out as 0; the test
  # uses 2^600, about 4e180) and underflow at 2^-700 (d comes out as Inf;
  # at 1e-160 it comes out as 1.02460, wrong in its fifth digit).
  g1 <- c(3, 5, 4, 6, 2, 7, 5, 4)
  g2 <- c(2, 4, 3, 3, 1, 5, 4, 2)
  for (cost in list(2, c(1, 3), c(3, 1))) {
    ref <- seq_mr_smd(g1, g2, A = 200, c = cost)
    expect_equal(val(ref, "smd"), 1.5 / sqrt(15 / 7))
    expect_identical(c(val(ref, "mean_1"), val(ref, "mean_2")), c(4.5, 3))
    for (scale in c(2^600, 2^-700, 2^1000, 2^-1000, 2^40)) {
      out <- seq_mr_smd(g1 * scale, g2 * scale, A = 200, c = cost)
      means <- out$term %in% c("mean_1", "mean_2")
      expect_identical(out$value[!means], ref$value[!means])
      expect_identical(out$value[means], c(4.5, 3) * scale)
      expect_identical(attributes(out)[c("is_satisfied", "criterion_n",
                                         "notes")],
                       attributes(ref)[c("is_satisfied", "criterion_n",
                                         "notes")])
    }
  }
  # Units that are not a power of two change the result only by rounding.
  ref <- seq_mr_smd(g1, g2, A = 200, c = 2)
  for (units in c(1e200, 1e-200, 1e-160, 2.54)) {
    out <- seq_mr_smd(g1 * units, g2 * units, A = 200, c = 2)
    expect_equal(val(out, "smd"), 1.5 / sqrt(15 / 7), tolerance = 1e-12)
    expect_identical(val(out, "optimal_n1"), val(ref, "optimal_n1"))
    expect_identical(attr(out, "notes"), character(0))
  }
})

test_that("seq_mr_smd needs four observations per group, as the article's pilot does", {
  # The article's stopping time is the smallest n at or above the pilot
  # size, and the pilot size is at least m_0 >= 4 (2017, p. 100). The rule
  # was reported as met for c(1, 2) against c(2, 4) with A = 1, c = 1.
  expect_error(seq_mr_smd(c(1, 2), c(2, 4), A = 1, c = 1),
               "'group_1' needs at least 4 observations")
  expect_error(seq_mr_smd(c(1, 2, 4, 3), c(2, 4, 3), A = 1, c = c(1, 2)),
               "'group_2' needs at least 4 observations")
  # Counted after missing values are removed.
  expect_error(seq_mr_smd(c(1, 2, 4, NA, 3), c(2, 4, 3, NA), A = 1, c = 1),
               paste("'group_2' needs at least 4 observations; it has 3.",
                     "Collect more data, then check the stopping rule again."),
               fixed = TRUE)
  expect_s3_class(seq_mr_smd(c(1, 2, 4, 3), c(2, 4, 3, 5), A = 1, c = 1),
                  "smsd_check")
})

test_that("a note explains a check that is not met although the totals suggest it is", {
  # One cost, d = 1.00, A = 10,000, c = 2.40: the criterion is 75 per group
  # at n = 74 (see the worked example above). With 80 and 74 observations
  # the total, 154, exceeds the total criterion, 150, but group 2 is short.
  dat <- smd_exact(80, 74, 1)
  out <- seq_mr_smd(dat$g1, dat$g2, A = 10000, c = 2.40)
  expect_false(attr(out, "is_satisfied"))
  expect_gte(attr(out, "current_n"), attr(out, "criterion_n"))
  expect_length(attr(out, "notes"), 1L)
  expect_match(attr(out, "notes"),
               "Group 2 has 74 observations and its criterion is 75",
               fixed = TRUE)
  expect_output(print(out), "Group 2 has 74 observations", fixed = TRUE)
  swapped <- seq_mr_smd(dat$g2, dat$g1, A = 10000, c = 2.40)
  expect_match(attr(swapped, "notes"),
               "Group 1 has 74 observations and its criterion is 75",
               fixed = TRUE)
  # Two costs (2 and 8), A = 50,000, d = 0, 200 and 83 observations: the
  # criteria are 170 and 84 (derived in the test of the paired size above).
  dat <- smd_exact(200, 83, 0)
  out <- seq_mr_smd(dat$g1, dat$g2, A = 50000, c = c(2, 8))
  expect_false(attr(out, "is_satisfied"))
  expect_match(attr(out, "notes"),
               "Group 2 has 83 observations and its criterion is 84",
               fixed = TRUE)
  # No note when the rule is met, or when the totals already show that it
  # is not (100 and 60 observations against criteria of 175 and 81).
  dat <- smd_exact(80, 75, 1)
  expect_identical(attr(seq_mr_smd(dat$g1, dat$g2, A = 10000, c = 2.40),
                        "notes"), character(0))
  dat <- smd_exact(100, 60, 0)
  expect_identical(attr(seq_mr_smd(dat$g1, dat$g2, A = 50000, c = c(2, 8)),
                        "notes"), character(0))
})

test_that("seq_mr_smd and its pilot validate every argument", {
  set.seed(113)
  g1 <- rnorm(10); g2 <- rnorm(10)
  expect_error(seq_mr_smd(g1, g2, A = 100, c = c(4, -1)),
               "'c[2]' must be a single finite number greater than 0",
               fixed = TRUE)
  expect_error(seq_mr_smd(g1, g2, A = 100, c = c(1, NA)), "'c[2]'",
               fixed = TRUE)
  expect_error(seq_mr_smd(g1, g2, A = 100, c = c(1, 2, 3)),
               "one positive number .* or two positive numbers")
  expect_error(seq_mr_smd(g1, g2, A = 100, c = "a"), "'c' must be one")
  expect_error(seq_mr_smd(g1, g2, A = 100, c = 0), "'c' must be a single")
  expect_error(seq_mr_smd(g1, g2, A = NA, c = 1), "'A' must be a single")
  expect_error(seq_mr_smd(g1, g2, A = c(1, 2), c = 1), "'A' must be a single")
  expect_error(seq_mr_smd(g1, g2, A = -5, c = 1), "'A' must be a single")
  # gamma is strictly between 0 and 1/2 (2017, p. 100).
  for (bad in list(0.6, 0.5, 0, -0.1, NA, "a")) {
    expect_error(seq_mr_smd(g1, g2, A = 100, c = 1, gamma = bad),
                 "'gamma' must be a single finite number greater than 0")
    expect_error(seq_mr_smd_pilot(A = 100, c = 1, gamma = bad), "'gamma'")
  }
  expect_error(seq_mr_smd(c(g1, NA), g2, A = 100, c = 1, na.rm = FALSE),
               "'group_1' contains missing values")
  expect_error(seq_mr_smd(g1, c(g2, NA), A = 100, c = 1, na.rm = FALSE),
               "'group_2' contains missing values")
  expect_error(seq_mr_smd(g1, g2, A = 100, c = 1, na.rm = NA), "'na.rm'")
  expect_error(seq_mr_smd(g1, c(g2, Inf), A = 100, c = 1), "infinite")
  expect_error(seq_mr_smd(cbind(g1, g1), g2, A = 100, c = 1),
               "'group_1' must be a numeric vector")
  # Missing values are dropped by default.
  expect_identical(seq_mr_smd(c(g1, NA), g2, A = 100, c = 1),
                   seq_mr_smd(g1, g2, A = 100, c = 1))
  # The pilot is not silently truncated or coerced to NA.
  expect_error(seq_mr_smd_pilot(A = 10, c = 1, m_0 = 4.5),
               "'m_0' must be a single whole number of at least 4")
  expect_error(seq_mr_smd_pilot(A = 10, c = 1, m_0 = 3), "'m_0'")
  expect_error(seq_mr_smd_pilot(A = 1e30, c = 1), "too large to represent")
  expect_error(seq_mr_smd_pilot(A = 1e30, c = c(1, 2)),
               "too large to represent")
  expect_error(seq_mr_smd_pilot(A = Inf, c = 1), "'A' must be a single")
  # The message names what to look at (a minimum risk procedure has no
  # precision target), whichever group is the less costly one.
  expect_error(seq_mr_smd_pilot(A = 1e30, c = 1),
               "pilot sample size .* is too large to represent. Check 'A' and 'c'")
  for (cost in list(1e-300, c(1e-300, 2e-300), c(2e-300, 1e-300))) {
    expect_error(seq_mr_smd(g1, g2, A = 1e300, c = cost),
                 paste("optimal sample size (Inf) is too large to represent.",
                       "Check 'A', 'c', and the data"),
                 fixed = TRUE)
  }
})

# ---------------------------------------------------------------------------
# Minimum risk point estimation of the coefficient of variation:
# Chattopadhyay and Kelley (2016). The published worked example is pinned
# with variance = "published" (the article's Equation 22), with
# DMAR::mr_cv() as a second oracle; the default variance is checked against
# the closed form population value for a gamma population, where the two
# estimators differ.
# ---------------------------------------------------------------------------

# The data of the worked example, Chattopadhyay and Kelley (2016), p. 638.
ck2016_books <- c(36, 53, 19, 11, 10, 24, 14, 65, 18, 48, 25, 35, 13, 18, 3,
                  41, 5, 3, 44, 26, 13, 39, 2, 3, 26, 22, 8, 15, 12, 22, 5,
                  21, 23, 40, 18)

# The two variance estimates written from the textbook h-statistics (the
# unbiased estimators of the third and fourth central moments in terms of
# the sample central moments m2, m3, m4), independently of the package's
# .mu3(), .mu4(), .var_cv(), and .var_cv_published(). Neither is bounded.
cv_h_statistics <- function(x) {
  n <- length(x)
  m2 <- mean((x - mean(x))^2)
  m3 <- mean((x - mean(x))^3)
  m4 <- mean((x - mean(x))^4)
  list(h3 = n^2 * m3 / ((n - 1) * (n - 2)),
       h4 = n * ((n^2 - 2 * n + 3) * m4 - 3 * (2 * n - 3) * m2^2) /
         ((n - 1) * (n - 2) * (n - 3)))
}
# Kelley, Darku, and Chattopadhyay (2018), Equation 32.
cv_eq32 <- function(x) {
  h <- cv_h_statistics(x)
  h$h4 / (4 * mean(x)^2 * var(x)) - var(x) / (4 * mean(x)^2) -
    h$h3 / mean(x)^3 + var(x)^2 / mean(x)^4
}
# Chattopadhyay and Kelley (2016), Equation 22.
cv_eq22 <- function(x) {
  h <- cv_h_statistics(x)
  var(x)^2 / (4 * mean(x)^4) + h$h4 / (4 * mean(x)^4) +
    var(x) / (2 * mean(x)^2) - h$h3 / mean(x)^3
}

test_that("seq_mr_cv computes the coefficient of variation as sd / mean", {
  set.seed(113)
  x <- rgamma(60, shape = 4, scale = 2)
  out <- seq_mr_cv(x, A = 400, c = 1)
  cv <- out$value[out$term == "cv"]
  expect_equal(cv, sd(x) / mean(x))
})

test_that("seq_mr_cv with the published variance reproduces the article's example", {
  # Chattopadhyay and Kelley (2016), p. 638 (A = 400,000, c = 75): the
  # printed output is risk 5964.345, cv 0.7391157, FALSE at n = 18; risk
  # 6224.861, cv 0.7113385, FALSE at n = 19; and risk 4891.284,
  # cv 0.7013904, TRUE at n = 35, "the first instance of the criterion
  # being satisfied". The printed risks have three decimals and the printed
  # estimates seven, which sets the (relative) tolerances.
  check <- function(n) {
    seq_mr_cv(ck2016_books[seq_len(n)], A = 400000, c = 75,
              variance = "published")
  }
  out <- check(18)
  expect_equal(val(out, "risk"), 5964.345, tolerance = 2e-7)
  expect_equal(val(out, "cv"), 0.7391157, tolerance = 1e-7)
  expect_false(attr(out, "is_satisfied"))
  out <- check(19)
  expect_equal(val(out, "risk"), 6224.861, tolerance = 2e-7)
  expect_equal(val(out, "cv"), 0.7113385, tolerance = 1e-7)
  expect_false(attr(out, "is_satisfied"))
  met <- vapply(18:35, function(n) attr(check(n), "is_satisfied"), logical(1))
  expect_identical(which(met), 18L)                 # n = 35 only
  out <- check(35)
  expect_equal(val(out, "risk"), 4891.284, tolerance = 2e-7)
  expect_equal(val(out, "cv"), 0.7013904, tolerance = 1e-7)
  expect_identical(attr(out, "variance"), "published")
  expect_match(attr(out, "notes"), "Equation 22", fixed = TRUE)
})

test_that("seq_mr_cv with the published variance agrees with DMAR::mr_cv", {
  for (n in c(18, 25, 34, 35)) {
    x <- ck2016_books[seq_len(n)]
    out <- seq_mr_cv(x, A = 400000, c = 75, variance = "published")
    ref <- DMAR::mr_cv(data = x, A = 400000, sampling_cost = 75)
    expect_equal(val(out, "risk"), ref$value[ref$term == "risk"])
    expect_equal(val(out, "optimal_n"), ceiling(sqrt(attr(ref, "criterion"))))
    expect_identical(attr(out, "is_satisfied"), attr(ref, "is_satisfied"))
  }
})

test_that("the default variance differs from the published one on the example", {
  # With the corrected variance the same 35 observations ask for one more.
  # Kelley, Darku, and Chattopadhyay (2018), Equation 32, computed here from
  # the textbook h-statistics (the unbiased estimators of the third and
  # fourth central moments written in the sample central moments m2, m3,
  # m4): the estimate is 0.2058226, so the criterion is
  # sqrt(400000 / 75 * (0.2058226 + 35^-0.98)) = 35.52, which rounds up to
  # 36. (Equation 22 gives 0.1982998 and 34.95, hence 35.)
  x <- ck2016_books
  n <- length(x)
  m2 <- mean((x - mean(x))^2)
  m3 <- mean((x - mean(x))^3)
  m4 <- mean((x - mean(x))^4)
  h3 <- n^2 * m3 / ((n - 1) * (n - 2))
  h4 <- n * ((n^2 - 2 * n + 3) * m4 - 3 * (2 * n - 3) * m2^2) /
    ((n - 1) * (n - 2) * (n - 3))
  eq32 <- h4 / (4 * mean(x)^2 * var(x)) - var(x) / (4 * mean(x)^2) -
    h3 / mean(x)^3 + var(x)^2 / mean(x)^4
  expect_equal(eq32, 0.2058226, tolerance = 1e-6)
  out <- seq_mr_cv(ck2016_books, A = 400000, c = 75)
  expect_identical(attr(out, "variance"), "corrected")
  expect_identical(attr(out, "notes"), character(0))
  expect_equal(val(out, "xi_squared"), eq32, tolerance = 1e-10)
  expect_equal(val(out, "optimal_n"), 36)
  expect_false(attr(out, "is_satisfied"))
  expect_identical(seq_mr_cv(ck2016_books, A = 400000, c = 75,
                             variance = "corrected"), out)
  expect_identical(seq_mr_cv(ck2016_books, A = 400000, c = 75,
                             variance = "pub"),
                   seq_mr_cv(ck2016_books, A = 400000, c = 75,
                             variance = "published"))
})

test_that("seq_mr_cv uses gamma", {
  # The article's printed risk at n = 35 (4891.284 = A V^2 / n + c n) gives
  # V^2 = (4891.284 - 2625) * 35 / 400000 = 0.1983. With gamma = 0.25 the
  # criterion is sqrt(400000 / 75 * (0.1983 + 35^-0.5)) = 44.26, so 45 and
  # not met; with the default it is 34.95, so 35 and met.
  out <- seq_mr_cv(ck2016_books, A = 400000, c = 75, variance = "published",
                   gamma = 0.25)
  expect_equal(val(out, "optimal_n"), 45)
  expect_false(attr(out, "is_satisfied"))
  expect_identical(attr(out, "gamma"), 0.25)
  ref <- DMAR::mr_cv(data = ck2016_books, A = 400000, sampling_cost = 75,
                     gamma = 0.25)
  expect_equal(val(out, "optimal_n"), ceiling(sqrt(attr(ref, "criterion"))))
  expect_identical(attr(out, "is_satisfied"), attr(ref, "is_satisfied"))
  expect_equal(val(seq_mr_cv(ck2016_books, A = 400000, c = 75,
                             variance = "published"), "optimal_n"), 35)
})

test_that("a gamma population distinguishes the two variance estimators", {
  # For a gamma population with shape a (scale 1): mu = a, sigma^2 = a,
  # mu3 = 2a, mu4 = 3a^2 + 6a, so kappa^2 = 1 / a. The asymptotic variance
  # of the sample coefficient of variation (Kelley, Darku, & Chattopadhyay,
  # 2018, Equation 30) is
  #   mu4 / (4 mu^2 sigma^2) - sigma^2 / (4 mu^2) - mu3 / mu^3
  #     + sigma^4 / mu^4 = 1 / (2a) + 1 / (2a^2),
  # and Chattopadhyay and Kelley (2016), Equation 9, is
  #   mu4 / (4 mu^4) + sigma^4 / (4 mu^4) - mu3 / mu^3 + sigma^2 / (2 mu^2)
  #     = 1 / (2a) - 1 / a^2 + 3 / (2a^3).
  # At a = 6.25 (kappa = 0.4) these are 0.0928 and 0.060544.
  a <- 6.25
  general <- 1 / (2 * a) + 1 / (2 * a^2)
  eq9 <- 1 / (2 * a) - 1 / a^2 + 3 / (2 * a^3)
  expect_equal(c(general, eq9), c(0.0928, 0.060544), tolerance = 1e-12)
  set.seed(113)
  x <- rgamma(200000, shape = a)
  corrected <- seq_mr_cv(x, A = 40000, c = 1)
  published <- seq_mr_cv(x, A = 40000, c = 1, variance = "published")
  # At n = 200,000 the relative standard deviation of the estimators is
  # 1.2% (corrected) and 0.4% (published) (60 samples of
  # rgamma(200000, shape = 6.25) after set.seed(1)), so 0.05 is more than
  # four standard deviations; the two population values differ by 53%.
  expect_equal(val(corrected, "xi_squared"), general, tolerance = 0.05)
  expect_equal(val(published, "xi_squared"), eq9, tolerance = 0.05)
  # Criterion sample sizes: sqrt(40000 * 0.0928) = 60.9 against
  # sqrt(40000 * 0.060544) = 49.2, a ratio of 1.24.
  expect_gt(val(corrected, "optimal_n"), 1.15 * val(published, "optimal_n"))
})

test_that("under normality the two variance estimators agree", {
  # Both expressions reduce to kappa^4 + kappa^2 / 2 when mu3 = 0 and
  # mu4 = 3 sigma^4. With kappa = 0.25 that is 0.03515625. At n = 200,000
  # the relative standard deviations are 0.7% and 0.5%, and 0.5% for the
  # difference (60 samples of rnorm(200000, 4, 1) after set.seed(1)); each
  # tolerance is at least four of those.
  set.seed(113)
  x <- rnorm(200000, mean = 4, sd = 1)
  corrected <- val(seq_mr_cv(x, A = 400, c = 1), "xi_squared")
  published <- val(seq_mr_cv(x, A = 400, c = 1, variance = "published"),
                   "xi_squared")
  expect_equal(corrected, 0.25^4 + 0.25^2 / 2, tolerance = 0.03)
  expect_equal(published, 0.25^4 + 0.25^2 / 2, tolerance = 0.03)
  expect_equal(corrected, published, tolerance = 0.02)
})

test_that("seq_mr_cv_pilot gives the article's pilot size", {
  # Chattopadhyay and Kelley (2016), p. 637: A = 400,000 and c = 75 give 18
  # (Equation 24: (400000 / 75)^(1 / 2.98) = 17.8).
  expect_identical(seq_mr_cv_pilot(A = 400000, c = 75), 18L)
  expect_identical(seq_mr_cv_pilot(A = 1000 / 0.05^2, c = 75, gamma = 0.49),
                   18L)
  # m_0 is a floor: 10^(1 / 2.98) = 2.17.
  expect_identical(seq_mr_cv_pilot(A = 10, c = 1), 4L)
  expect_identical(seq_mr_cv_pilot(A = 10, c = 1, m_0 = 12), 12L)
  # A smaller gamma gives a larger pilot: 400^(1 / 2.5) = 10.99.
  expect_identical(seq_mr_cv_pilot(A = 400, c = 1, gamma = 0.25), 11L)
})

test_that("seq_mr_cv returns the documented rows and attributes", {
  set.seed(113)
  x <- rgamma(80, shape = 4, scale = 2)
  out <- seq_mr_cv(x, A = 400, c = 2)
  expect_identical(out$term,
                   c("cv", "risk", "xi_squared", "optimal_n", "current_n"))
  expect_true(all(is.finite(out$value)))
  expect_identical(attr(out, "method"),
                   "minimum risk point estimation (coefficient of variation)")
  expect_identical(attr(out, "current_n"), 80L)
  expect_equal(attr(out, "criterion_n"), val(out, "optimal_n"))
  expect_type(attr(out, "criterion_n"), "double")
  expect_identical(attr(out, "loss_A"), 400)
  expect_identical(attr(out, "cost_c"), 2)
  expect_identical(attr(out, "gamma"), 0.49)
  expect_identical(attr(out, "variance"), "corrected")
  expect_false(attr(out, "variance_unusable"))
  expect_null(attr(out, "conf_level"))
  # The decision is the criterion on the scale of n (Equation 23), and the
  # risk row is Equation 14 with the estimate: A xi^2 / n + c n.
  expect_identical(attr(out, "is_satisfied"), 80 >= val(out, "optimal_n"))
  expect_equal(val(out, "risk"), 400 * val(out, "xi_squared") / 80 + 2 * 80)
})

test_that("a variance estimate that is zero, negative, or not finite means not met", {
  # Constant data: no variance can be estimated. Read literally, the
  # article's rule would be met at once, with cv = 0.
  for (v in c("corrected", "published")) {
    out <- seq_mr_cv(rep(5, 20), A = 400, c = 1, variance = v)
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_true(attr(out, "variance_unusable"))
    expect_equal(val(out, "cv"), 0)
    expect_equal(val(out, "xi_squared"), 20^-3)
    expect_identical(attr(out, "notes")[1], .note_variance_unusable())
    expect_match(attr(out, "notes")[1], "zero, negative, or not finite",
                 fixed = TRUE)
    # The criterion is at least one more observation, so that
    # current_n >= optimal_n says the same as is_satisfied. (The bound alone
    # would give sqrt(400 (20^-3 + 20^-0.98)) = 4.6, so 5.)
    expect_equal(val(out, "optimal_n"), 21)
    expect_equal(attr(out, "criterion_n"), 21)
  }
  expect_output(print(seq_mr_cv(rep(5, 20), A = 400, c = 1)),
                "Stopping rule met: FALSE", fixed = TRUE)
  # Constant up to rounding error (a relative range of 8e-16) is constant:
  # the moment ratios would be rounding noise.
  near <- 5 + 1e-15 * rep(c(-1, 1, 0, 2, -2), 4)
  expect_gt(stats::sd(near), 0)
  for (v in c("corrected", "published")) {
    out <- seq_mr_cv(near, A = 400, c = 1, variance = v)
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_true(attr(out, "variance_unusable"))
    expect_equal(val(out, "optimal_n"), 21)
  }

  # Data that vary, with a negative estimate. For c(9, 11, 9, 11, 9, 11)
  # the fourth h-statistic is zero, so Equation 32 is
  # -s^2 / (4 xbar^2) + s^4 / xbar^4 = -1.2 / 400 + 1.44 / 10000
  # = -0.002856, and Equation 22 is
  # s^4 / (4 xbar^4) + s^2 / (2 xbar^2) = 0.000036 + 0.006 = 0.006036.
  x <- c(9, 11, 9, 11, 9, 11)
  expect_equal(cv_eq32(x), -0.002856, tolerance = 1e-10)
  expect_equal(cv_eq22(x), 0.006036, tolerance = 1e-10)
  out <- seq_mr_cv(x, A = 100, c = 1)
  expect_identical(attr(out, "is_satisfied"), FALSE)
  expect_true(attr(out, "variance_unusable"))
  expect_equal(val(out, "xi_squared"), 6^-3)
  expect_equal(val(out, "optimal_n"), 7)
  expect_equal(val(out, "cv"), sqrt(1.2) / 10)
  expect_identical(attr(out, "notes"), .note_variance_unusable())
  # The published estimate is positive for the same data, so there the rule
  # is evaluated: sqrt(100 (0.006036 + 6^-0.98)) = 4.23, so 5 <= 6.
  out <- seq_mr_cv(x, A = 100, c = 1, variance = "published")
  expect_false(attr(out, "variance_unusable"))
  expect_equal(val(out, "xi_squared"), 0.006036, tolerance = 1e-10)
  expect_equal(val(out, "optimal_n"), 5)
  expect_true(attr(out, "is_satisfied"))

  # A negative Equation 22: c(2, 3, 2, 3, 2, 3, 2, 12) is strongly right
  # skewed, and the term -mu3 / xbar^3 outweighs the rest (-0.1198).
  y <- c(2, 3, 2, 3, 2, 3, 2, 12)
  expect_equal(cv_eq22(y), -0.119783, tolerance = 1e-5)
  out <- seq_mr_cv(y, A = 60, c = 1, variance = "published")
  expect_identical(attr(out, "is_satisfied"), FALSE)
  expect_true(attr(out, "variance_unusable"))
  expect_equal(val(out, "xi_squared"), 8^-3)
  expect_equal(val(out, "optimal_n"), 9)
  expect_identical(attr(out, "notes"),
                   c(.note_variance_unusable(), .note_cv_published()))
  # This is where the function departs from the rule as printed, and from
  # DMAR::mr_cv(), which applies Equation 22 with no bound and stops:
  # 60 (-0.1198 + 8^-0.98) = 0.63 <= 64.
  ref <- DMAR::mr_cv(data = y, A = 60, sampling_cost = 1)
  expect_equal(attr(ref, "criterion"), 60 * (cv_eq22(y) + 8^-0.98))
  expect_true(attr(ref, "is_satisfied"))
})

test_that("a small positive variance estimate is used, bounded below by n^-3", {
  # Heights in centimeters with a coefficient of variation near 0.02. Both
  # estimates are positive and below 8^-3 = 0.00195 (the true value,
  # kappa^2 / 2 + kappa^4, is about 0.0002). Such an estimate is small, not
  # unusable: it is used, with the bound of Kelley, Darku, and Chattopadhyay
  # (2018, Equation 31) in its place, and sampling does not go on because
  # of it.
  heights <- c(168.2, 171.9, 166.4, 173.5, 170.8, 169.1, 174.6, 167.3)
  expect_gt(cv_eq32(heights), 0)
  expect_lt(cv_eq32(heights), 8^-3)
  expect_gt(cv_eq22(heights), 0)
  expect_lt(cv_eq22(heights), 8^-3)
  for (v in c("corrected", "published")) {
    # A = 400: sqrt(400 (8^-3 + 8^-0.98)) = sqrt(400 * 0.132264) = 7.27, so
    # the criterion is 8 and the rule is met at n = 8.
    out <- seq_mr_cv(heights, A = 400, c = 1, variance = v)
    expect_false(attr(out, "variance_unusable"))
    expect_identical(val(out, "xi_squared"), 8^-3)
    expect_equal(val(out, "optimal_n"), 8)
    expect_true(attr(out, "is_satisfied"))
    expect_equal(val(out, "risk"), 400 * 8^-3 / 8 + 8)
    expect_length(attr(out, "notes"), if (v == "published") 1L else 0L)
    # A = 500: sqrt(500 * 0.132264) = 8.13, so 9, and the rule is not met;
    # the estimate is still usable.
    out <- seq_mr_cv(heights, A = 500, c = 1, variance = v)
    expect_false(attr(out, "variance_unusable"))
    expect_equal(val(out, "optimal_n"), 9)
    expect_false(attr(out, "is_satisfied"))
  }
  # The decision is that of DMAR::mr_cv() at both values of A. The risk
  # differs, because DMAR uses Equation 22 itself where this function uses
  # the bound.
  for (A in c(400, 500)) {
    out <- seq_mr_cv(heights, A = A, c = 1, variance = "published")
    ref <- DMAR::mr_cv(data = heights, A = A, sampling_cost = 1)
    expect_identical(attr(out, "is_satisfied"), attr(ref, "is_satisfied"))
    expect_equal(ref$value[ref$term == "risk"],
                 A * cv_eq22(heights) / 8 + 8)
    expect_gt(val(out, "risk"), ref$value[ref$term == "risk"])
  }
  # Between those two values of A the decisions can differ, and the help
  # says so. The bound is larger than Equation 22 (0.00195 against 0.00015),
  # so this function's criterion is the larger one, and n^2 = 64 lies
  # between the two for A from 64 / (8^-3 + 8^-0.98) = 483.9 to
  # 64 / (0.00015 + 8^-0.98) = 490.6. At A = 485:
  # 485 (8^-3 + 8^-0.98) = 485 * 0.132262 = 64.15 > 64, so the criterion
  # here is sqrt(64.15) = 8.009, reported as 9, and the rule is not met;
  # 485 (0.00015 + 8^-0.98) = 485 * 0.130460 = 63.27 <= 64, so mr_cv()
  # stops. The estimate is usable in both.
  expect_equal(cv_eq22(heights), 0.000150514, tolerance = 1e-5)
  expect_gt(485 * (8^-3 + 8^-0.98), 64)
  expect_lt(485 * (cv_eq22(heights) + 8^-0.98), 64)
  for (v in c("corrected", "published")) {
    out <- seq_mr_cv(heights, A = 485, c = 1, variance = v)
    expect_false(attr(out, "variance_unusable"))
    expect_identical(val(out, "xi_squared"), 8^-3)
    expect_equal(val(out, "optimal_n"), 9)
    expect_false(attr(out, "is_satisfied"))
  }
  ref <- DMAR::mr_cv(data = heights, A = 485, sampling_cost = 1)
  expect_equal(attr(ref, "criterion"), 485 * (cv_eq22(heights) + 8^-0.98))
  expect_true(attr(ref, "is_satisfied"))
  # The result does not depend on the units: the same heights in meters and
  # in millimeters.
  for (v in c("corrected", "published")) {
    cm <- seq_mr_cv(heights, A = 400, c = 1, variance = v)
    for (units in c(1 / 100, 10)) {
      other <- seq_mr_cv(heights * units, A = 400, c = 1, variance = v)
      expect_identical(attr(other, "is_satisfied"), attr(cm, "is_satisfied"))
      expect_identical(val(other, "optimal_n"), val(cm, "optimal_n"))
      expect_false(attr(other, "variance_unusable"))
    }
  }
})

test_that("seq_mr_cv decides on the criterion it reports", {
  # The criterion is sqrt((A / c) (xi^2 + n^(-2 gamma))), rounded up by the
  # shared .as_sample_size(), whose ceiling ignores representation error;
  # the decision uses the same rounding. For the 35 observations of the
  # article's example, A0 below makes the criterion exactly 35 (with c = 1
  # and the default variance, Equation 32 written from the h-statistics).
  n <- 35
  A0 <- n^2 / (cv_eq32(ck2016_books) + n^-0.98)
  # Above 35 by a relative 5e-13 (representation error): the criterion is
  # 35 and the rule is met, not "not met" beside a reported criterion of 35.
  out <- seq_mr_cv(ck2016_books, A = A0 * (1 + 1e-12), c = 1)
  expect_equal(val(out, "optimal_n"), 35)
  expect_true(attr(out, "is_satisfied"))
  # Really above 35 (by 2e-5): the criterion is 36 and the rule is not met.
  out <- seq_mr_cv(ck2016_books, A = A0 * (1 + 1e-6), c = 1)
  expect_equal(val(out, "optimal_n"), 36)
  expect_false(attr(out, "is_satisfied"))
  # On every path with an estimate, is_satisfied is current_n >= optimal_n,
  # including when the variance cannot be estimated.
  withr::local_seed(113)
  for (i in 1:30) {
    k <- sample(6:60, 1)
    x <- if (i %% 5 == 0) rep(3, k) else rgamma(k, shape = 4)
    out <- seq_mr_cv(x, A = sample(c(50, 400, 3000), 1), c = 1,
                     variance = c("corrected", "published")[1 + i %% 2])
    expect_identical(attr(out, "is_satisfied"), k >= val(out, "optimal_n"))
  }
})

test_that("neither the estimate nor either variance estimate depends on the scale of the data", {
  # Multiplying by a power of two is exact, so the whole object must be
  # identical. Equation 22 evaluated on the data as given overflows at
  # 2^270 (about 1.9e81: the fourth power of the mean is not representable)
  # and underflows at 2^-280; the published estimate was then reported as
  # unusable while the default one was not. The standard deviation of the
  # data as given overflows at 2^600 (about 4e180: from 1e154 on there used
  # to be no estimate) and underflows at 2^-700 (about 2e-211: the cv row
  # was 0 beside a rule reported as met).
  set.seed(113)
  x <- rgamma(30, shape = 4, rate = 1)
  for (v in c("corrected", "published")) {
    ref <- seq_mr_cv(x, A = 400, c = 1, variance = v)
    expect_false(attr(ref, "variance_unusable"))
    expect_equal(val(ref, "cv"), sd(x) / mean(x))
    for (scale in c(2^270, 2^-280, 2^40, 2^600, 2^-700, 2^1000, 2^-1000)) {
      expect_identical(seq_mr_cv(x * scale, A = 400, c = 1, variance = v),
                       ref)
    }
  }
  # And the value is still Equation 22 (Equation 32 for the default),
  # written from the h-statistics.
  expect_equal(val(seq_mr_cv(x, A = 400, c = 1, variance = "published"),
                   "xi_squared"), cv_eq22(x), tolerance = 1e-12)
  expect_equal(val(seq_mr_cv(x, A = 400, c = 1), "xi_squared"), cv_eq32(x),
               tolerance = 1e-12)
  # The helper for Equation 22 is free of the units on its own account (it
  # works on the data divided by their mean), whatever its caller does.
  expect_equal(.var_cv_published(x, floor = FALSE), cv_eq22(x),
               tolerance = 1e-12)
  for (scale in c(2^270, 2^-280)) {
    expect_identical(.var_cv_published(x * scale, floor = FALSE),
                     .var_cv_published(x, floor = FALSE))
  }
  # Units that are not a power of two change the result only by rounding.
  # For c(1, 2, 3, 5) the mean is 2.75 and the sum of squared deviations
  # 8.75, so the coefficient of variation is sqrt(8.75 / 3) / 2.75
  # = 0.6210273. Computed on the data as given, sd() / mean() is Inf at
  # 1e200, 0.62106 at 1e-160, and 0 at 1e-200.
  y <- c(1, 2, 3, 5)
  ref <- seq_mr_cv(y, A = 400, c = 1)
  expect_equal(val(ref, "cv"), sqrt(8.75 / 3) / 2.75)
  for (units in c(1e200, 1e-160, 1e-200, 1e300, 1e-300, 2.54)) {
    out <- seq_mr_cv(y * units, A = 400, c = 1)
    expect_equal(val(out, "cv"), sqrt(8.75 / 3) / 2.75, tolerance = 1e-12)
    expect_equal(val(out, "xi_squared"), val(ref, "xi_squared"),
                 tolerance = 1e-10)
    expect_identical(val(out, "optimal_n"), val(ref, "optimal_n"))
    expect_identical(attr(out, "is_satisfied"), attr(ref, "is_satisfied"))
    expect_identical(attr(out, "notes"), character(0))
  }
})

test_that("a sample mean that is not positive is not estimable", {
  set.seed(113)
  x <- rgamma(20, shape = 4, rate = 1)
  negative <- seq_mr_cv(-x, A = 400, c = 1)          # old code: met, cv < 0
  zero <- seq_mr_cv(rep(0, 10), A = 400, c = 1)      # old code: an error
  balanced <- seq_mr_cv(c(-1, 1, -2, 2), A = 400, c = 1)
  for (out in list(negative, zero, balanced)) {
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_true(all(is.na(out$value[1:4])))
    expect_false(any(is.infinite(out$value) | is.nan(out$value)))
    expect_true(is.na(attr(out, "criterion_n")))
    expect_length(attr(out, "notes"), 1L)
  }
  # The same note whether the mean is negative, zero, or balanced at zero.
  for (out in list(negative, zero, balanced)) {
    expect_match(attr(out, "notes"), "sample mean is not positive",
                 fixed = TRUE)
  }
  expect_output(print(negative), "Stopping rule met: FALSE", fixed = TRUE)
})

test_that("a mean of zero up to rounding error is not positive, whatever its sign", {
  # Standardized scores have a mean of about 1e-17 with either sign. The
  # sign used to decide between a note and an error about an optimal sample
  # size near 1e34. Seeds 1 and 2 give one mean of each sign.
  signs <- numeric(0)
  for (seed in 1:6) {
    x <- withr::with_seed(seed, as.numeric(scale(rnorm(30))))
    signs <- c(signs, sign(mean(x)))
    out <- seq_mr_cv(x, A = 400, c = 1)
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_true(is.na(val(out, "cv")))
    expect_match(attr(out, "notes"), "sample mean is not positive",
                 fixed = TRUE)
  }
  expect_true(any(signs > 0) && any(signs < 0))
  out <- seq_mr_cv(c(-2, -1, 1, 2 + 4e-16), A = 400, c = 1)
  expect_gt(mean(c(-2, -1, 1, 2 + 4e-16)), 0)
  expect_match(attr(out, "notes"), "sample mean is not positive",
               fixed = TRUE)
  # A mean that is small but real is estimable: mean 0.01, SD 1.83.
  out <- seq_mr_cv(c(-2, -1, 1, 2.04), A = 400, c = 1)
  expect_equal(val(out, "cv"), sd(c(-2, -1, 1, 2.04)) / 0.01)
  # The test of the mean is relative to the size of the data, so it gives
  # the same answer in any units: standardized scores multiplied by 2^600
  # or 2^-700 still have no estimate, with the same note.
  z <- withr::with_seed(1, as.numeric(scale(rnorm(30))))
  for (scale in c(2^600, 2^-700)) {
    out <- seq_mr_cv(z * scale, A = 400, c = 1)
    expect_true(is.na(val(out, "cv")))
    expect_identical(attr(out, "notes"), .note_cv_mean_not_positive())
  }
})

test_that("for two equally likely values the rule is met at the pilot with probability one half", {
  # ?seq_mr_cv: observations 9 or 11 with equal probability, A = 400, c = 1,
  # pilot 8. The check depends on the data only through the number k of
  # elevens among the 8, so its probability of stopping at the pilot is a
  # sum of binomial probabilities. The rule is written out here from
  # Equation 32 of Kelley, Darku, and Chattopadhyay (2018) (cv_eq32() above)
  # and the stopping rule n^2 >= (A / c) (max(V, n^-3) + n^(-2 gamma)), with
  # a raw estimate that is not positive counted as not met (the package's
  # rule for an unusable variance).
  expect_identical(seq_mr_cv_pilot(A = 400, c = 1), 8L)
  met_by_hand <- vapply(0:8, function(k) {
    x <- c(rep(11, k), rep(9, 8 - k))
    v <- if (k %in% c(0, 8)) NA_real_ else cv_eq32(x)
    isTRUE(v > 0) && 8^2 >= 400 * (max(v, 8^-3) + 8^-0.98)
  }, logical(1))
  met <- vapply(0:8, function(k) {
    attr(seq_mr_cv(c(rep(11, k), rep(9, 8 - k)), A = 400, c = 1),
         "is_satisfied")
  }, logical(1))
  expect_identical(met, met_by_hand)
  expect_identical((0:8)[met], c(1L, 2L, 5L, 6L, 7L))
  # Of the 2^8 = 256 equally likely samples, 8 + 28 + 56 + 28 + 8 = 128 meet
  # the rule: exactly one half.
  expect_identical(sum(choose(8, (0:8)[met])), 128)
})

test_that("seq_mr_cv and its pilot validate every argument", {
  set.seed(113)
  x <- rgamma(10, shape = 4, rate = 1)
  expect_error(seq_mr_cv(x, A = -1, c = 1), "'A' must be a single")
  expect_error(seq_mr_cv(x, A = NA, c = 1), "'A' must be a single")
  expect_error(seq_mr_cv(x, A = Inf, c = 1), "'A' must be a single")
  expect_error(seq_mr_cv(x, A = 400, c = 0), "'c' must be a single")
  expect_error(seq_mr_cv(x, A = 400, c = Inf), "'c' must be a single")
  expect_error(seq_mr_cv(x, A = 400, c = c(1, 2)), "'c' must be a single")
  # gamma is strictly between 0 and 1/2 (2016, p. 636).
  for (bad in list(0.6, 0.5, 0, NA)) {
    expect_error(seq_mr_cv(x, A = 400, c = 1, gamma = bad),
                 "'gamma' must be a single finite number greater than 0")
    expect_error(seq_mr_cv_pilot(A = 400, c = 1, gamma = bad), "'gamma'")
  }
  expect_error(seq_mr_cv(x, A = 400, c = 1, variance = "eq22"),
               "'variance' must be one of")
  expect_error(seq_mr_cv(x, A = 400, c = 1, variance = NA),
               "'variance' must be one of")
  expect_error(seq_mr_cv(c(x, NA), A = 400, c = 1, na.rm = FALSE),
               "'x' contains missing values")
  expect_error(seq_mr_cv(c(x, Inf), A = 400, c = 1), "infinite")
  expect_error(seq_mr_cv(x[1:3], A = 400, c = 1), "at least 4")
  expect_error(seq_mr_cv(matrix(x, 2), A = 400, c = 1),
               "'x' must be a numeric vector")
  expect_identical(seq_mr_cv(c(x, NA), A = 400, c = 1),
                   seq_mr_cv(x, A = 400, c = 1))
  expect_error(seq_mr_cv_pilot(A = 10, c = 1, m_0 = 4.5),
               "'m_0' must be a single whole number of at least 4")
  expect_error(seq_mr_cv_pilot(A = 10, c = 1, m_0 = 3), "'m_0'")
  expect_error(seq_mr_cv_pilot(A = 1e30, c = 1), "too large to represent")
  expect_error(seq_mr_cv_pilot(A = 400, c = NA), "'c' must be a single")
  # A criterion too large to represent is an error that names what to look
  # at. Mean 0.00025 and SD 1.83 give a coefficient of variation near 7300,
  # whose variance estimate is of the order of 1e15.
  expect_error(seq_mr_cv_pilot(A = 1e30, c = 1),
               "pilot sample size .* is too large to represent. Check 'A' and 'c'")
  expect_error(seq_mr_cv(c(-2, -1, 1, 2.001), A = 1e6, c = 1),
               "is too large to represent. Check 'A', 'c', and the data")
  expect_error(seq_mr_cv(x, A = 1e300, c = 1e-300),
               "optimal sample size (Inf) is too large to represent",
               fixed = TRUE)
})

test_that("integer data give the same answer as the same values as doubles", {
  out_int <- seq_mr_cv(as.integer(ck2016_books), A = 400000, c = 75)
  expect_identical(out_int, seq_mr_cv(ck2016_books, A = 400000, c = 75))
})
