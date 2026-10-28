# Tests for the budget constrained K-population procedure of Chattopadhyay,
# Bandyopadhyay, Kelley, and Padalunkal (2025, Psychological Methods, 30(2),
# 425-439): seq_budget_optimal(), seq_budget_pilot(), seq_budget(),
# seq_budget_next(), the decide_seq_budget_*() functions, and
# seq_budget_sensitivity().
#
# Expected values are numbers printed in the article, or are worked out by
# hand in the comments from data built to have exact sample means and
# standard deviations. None is copied from the functions' own output.
# Where a comment cites a simulation ("independent simulation", "long run"),
# it gives the design and the number of replications, so the figure can be
# regenerated with seq_budget_sensitivity() or the loop the comment
# describes.

# A group of n observations whose sample mean and sample standard deviation
# are exactly `mean` and `sd` (up to floating point error).
exact_group <- function(n, mean, sd) {
  z <- seq_len(n) - (n + 1) / 2
  mean + sd * z / stats::sd(z)
}
exact_groups <- function(n, mean, sd) {
  Map(exact_group, n, mean, sd)
}
term <- function(check, name) check$value[check$term == name]

# The hand-computed design used throughout:
#   weights c = (1, -2, 1), unit costs a = (4, 9, 16), sample SDs s = (3, 1, 2).
#   sum |c| s sqrt(a) = 1*3*2 + 2*1*3 + 1*2*4 = 20.
#   Equation 18 targets for a budget A0: A0 |c_i| s_i / (sqrt(a_i) * 20)
#     = A0 * (3/40, 2/60, 2/80) = A0 * (0.075, 0.0333..., 0.025).
hand_w <- c(1, -2, 1)
hand_a <- c(4, 9, 16)
hand_s <- c(3, 1, 2)

# --- seq_budget_optimal: published values --------------------------------

test_that("seq_budget_optimal reproduces the article's p. 430 example", {
  # Three populations, 30 per observation, a budget of 10,000: "the optimal
  # sample size will be n_i0 = 111.11", 111 per group costs 9,990, and "a
  # ceiling function ... cannot be taken" (112 per group costs 10,080).
  # 111.11 = 10000 / (3 * 30) needs equal products |c_i| sigma_i.
  out <- seq_budget_optimal(c(-1, 0.5, 0.5), c(1, 2, 2), c(30, 30, 30),
                            10000)
  expect_equal(round(out$n_per_group, 2), rep(111.11, 3))
  expect_equal(sum(30 * floor(out$n_per_group)), 9990)
  expect_gt(sum(30 * ceiling(out$n_per_group)), 10000)
  # Minimized variance by hand: (3 * 1 * sqrt(30))^2 / 10000 = 270 / 10000.
  expect_equal(out$minimum_variance, 0.027)
})

test_that("seq_budget_optimal reproduces Table 1 of the article", {
  # Table 1 (p. 431): contrast weights (1, -2, 1), budget 5,000. Columns
  # n_i0 and V_o as printed (four decimals). Standard deviations follow
  # from the stated distributions.
  w <- c(1, -2, 1)
  mixture_sd <- function(m1, s1, m2, s2) {
    mu <- 0.95 * m1 + 0.05 * m2
    sqrt(0.95 * (s1^2 + m1^2) + 0.05 * (s2^2 + m2^2) - mu^2)
  }
  rows <- list(
    gamma = list(
      sd = c(0.2 * sqrt(16), 0.125 * sqrt(50), 1.50001 * sqrt(6)),
      cost = c(10, 12, 36),
      n = c(41.2035, 83.1149, 99.7385), V = 0.1885),
    logistic = list(
      sd = c(2, 1, 7) * pi / sqrt(3),
      cost = c(16, 14, 39),
      n = c(42.2309, 45.1468, 94.6731), V = 2.3058),
    normal = list(
      sd = c(4.5, 0.9, 1),
      cost = c(38, 15, 10),
      n = c(96.373, 61.3566, 41.7479), V = 0.2869),
    # The mixture row prints V_o as 0.8486, a misprint: with the standard
    # deviations that reproduce its three n_i0 to four decimals, Equation 13
    # gives 8.4788 (and the same row's simulated mean variance is 8.228).
    mixture = list(
      sd = c(mixture_sd(1, 19, 20, 15), mixture_sd(2, 5, 20, 8),
             mixture_sd(3, 11, 20, 11)),
      cost = c(38, 15, 10),
      n = c(75.913, 81.5976, 89.1344), V = 8.4788)
  )
  for (nm in names(rows)) {
    r <- rows[[nm]]
    out <- seq_budget_optimal(w, r$sd, r$cost, 5000)
    expect_lt(max(abs(out$n_per_group - r$n)), 5.5e-5, label = nm)
    expect_lt(abs(out$minimum_variance - r$V), 5.5e-5, label = nm)
    expect_equal(sum(r$cost * out$n_per_group), 5000, label = nm)
  }
})

test_that("seq_budget_optimal agrees with a constrained numerical optimizer", {
  # Minimize sum(c^2 sigma^2 / n) subject to sum(a n) = A0 directly, with
  # the budget shares parametrized on the simplex.
  withr::local_seed(113)
  for (i in 1:5) {
    K <- sample(2:5, 1)
    w <- stats::rnorm(K)
    sigma <- stats::runif(K, 0.2, 5)
    a <- stats::runif(K, 1, 50)
    A0 <- stats::runif(1, 500, 50000)
    objective <- function(theta) {
      p <- exp(c(0, theta))
      p <- p / sum(p)
      sum(w^2 * sigma^2 / (A0 * p / a))
    }
    fit <- stats::optim(rep(0, K - 1), objective, method = "BFGS",
                        control = list(reltol = 1e-14, maxit = 5000))
    p <- exp(c(0, fit$par))
    p <- p / sum(p)
    out <- seq_budget_optimal(w, sigma, a, A0)
    expect_equal(out$n_per_group, A0 * p / a, tolerance = 1e-4)
    expect_equal(out$minimum_variance, fit$value, tolerance = 1e-8)
  }
})

test_that("seq_budget_optimal gives no observations to a zero weight group", {
  # Hand: c = (1, -1, 0), sigma = 1, a = 10, A0 = 5000.
  # sum |c| sigma sqrt(a) = 2 sqrt(10); n = 5000 / (sqrt(10) * 2 sqrt(10)) = 250.
  out <- seq_budget_optimal(c(1, -1, 0), c(1, 1, 1), c(10, 10, 10), 5000)
  expect_equal(out$n_per_group, c(250, 250, 0))
})

test_that("seq_budget_optimal validates its arguments", {
  expect_error(seq_budget_optimal(1, 1, 1, 100), "at least 2 populations")
  expect_error(seq_budget_optimal(c(0, 0), c(1, 1), c(1, 1), 100),
               "at least one nonzero")
  # The vector arguments go through the package's shared check, which names
  # the argument, says what its entries are, and points at the bad entry.
  expect_error(seq_budget_optimal(c(1, NA), c(1, 1), c(1, 1), 100),
               "'contrast_weights' must contain only finite numbers \\(one contrast weight per group\\); entry 2 is NA")
  expect_error(seq_budget_optimal(c(1, -1), c(1, 1, 1), c(10, 10), 1000),
               "'sigmas' must have 2 entries \\(one per entry of 'contrast_weights'\\); it has 3")
  expect_error(seq_budget_optimal(c(1, -1), c(-1, 1), c(10, 10), 1000),
               "'sigmas' must contain only positive numbers \\(one standard deviation per population\\); entry 1 is -1")
  expect_error(seq_budget_optimal(c(1, -1), c(1, NA), c(10, 10), 1000),
               "'sigmas' must contain only finite numbers.*entry 2 is NA")
  expect_error(seq_budget_optimal(c(1, -1), c(1, 1), c(10, -10), 1000),
               "'cost_per_group' must contain only positive numbers.*entry 2 is -10")
  expect_error(seq_budget_optimal(c(1, -1), c(1, 1), c(10, 10, 10), 1000),
               "'cost_per_group' must have 2 entries.*it has 3")
  expect_error(seq_budget_optimal(numeric(0), c(1, 1), c(10, 10), 1000),
               "'contrast_weights' must be a numeric vector.*a vector of length 0")
  expect_error(seq_budget_optimal(matrix(c(1, -1), 1), c(1, 1), c(10, 10), 100),
               "'contrast_weights' must be a numeric vector")
  expect_error(seq_budget_optimal(c(1, -1), c(1, 1), c(10, 10), -1000),
               "'budget'.*greater than 0")
  expect_error(seq_budget_optimal(c(1, -1), c(1, 1), c(10, 10), NA),
               "'budget'")
  expect_error(seq_budget_optimal(c(1, -1), c(1, 1), c(10, 10), c(1, 2)),
               "'budget'")
  expect_error(seq_budget_optimal("a", c(1, 1), c(10, 10), 100),
               "'contrast_weights' must be a numeric vector")
})

test_that("seq_budget_optimal refuses vectors named in different orders", {
  # Entries are matched by position. Hand: c = (1, -1), sigma = (2, 1),
  # a = (4, 9), A0 = 100: sum |c| sigma sqrt(a) = 4 + 3 = 7;
  # n = 100 * (2 / (2 * 7), 1 / (3 * 7)) = (14.2857, 4.7619).
  w <- c(new = 1, old = -1)
  ok <- seq_budget_optimal(w, c(new = 2, old = 1), c(new = 4, old = 9), 100)
  expect_equal(ok$n_per_group, c(100 / 7, 100 / 21))
  # Names on one argument only are no ground for complaint.
  expect_equal(seq_budget_optimal(w, c(2, 1), c(4, 9), 100), ok)
  expect_error(seq_budget_optimal(w, c(old = 1, new = 2), c(4, 9), 100),
               "'sigmas' is named \\(old, new\\) but 'contrast_weights' is named \\(new, old\\)")
  expect_error(seq_budget_optimal(w, c(2, 1), c(old = 9, new = 4), 100),
               "'cost_per_group' is named \\(old, new\\).*matched to groups by position")
  expect_error(seq_budget_optimal(c(1, -1), c(a = 2, b = 1), c(x = 4, y = 9),
                                  100),
               "'cost_per_group' is named \\(x, y\\) but 'sigmas' is named \\(a, b\\)")
})

test_that("seq_budget_optimal stops when the result cannot be represented", {
  # Products that underflow to zero or overflow stop with a message; they
  # must not come back as NaN sample sizes with a variance of 0, or sizes of
  # 0 with a variance of Inf.
  expect_error(seq_budget_optimal(c(1e-300, -1e-300), c(1e-300, 1e-300),
                                  c(1e-300, 1e-300), 1e-300),
               "optimal allocation cannot be computed: the values in 'contrast_weights', 'sigmas', 'cost_per_group', and 'budget' are too large or too small.*Rescale")
  expect_error(seq_budget_optimal(c(1, -1), c(1e200, 1e200), c(1e200, 1e200),
                                  1),
               "optimal allocation cannot be computed.*Rescale")
  # A very large or very small unit of money alone is fine: the allocation
  # does not depend on it. Hand: equal groups, so n = A0 / (2 a) = 500 each,
  # and V_o = (2 * 3 * sqrt(a))^2 / A0 = 36 a / A0 = 0.036. (Multiplying the
  # budget by the shares before dividing by the costs would overflow at
  # 1e300 and underflow at 1e-300.)
  for (unit in c(1e-300, 1e300)) {
    out <- seq_budget_optimal(c(1, -1), c(3, 3), c(unit, unit), 1000 * unit)
    expect_equal(out$n_per_group, c(500, 500))
    expect_equal(out$minimum_variance, 0.036)
  }
})


# --- seq_budget_pilot ----------------------------------------------------

test_that("seq_budget_pilot returns the chosen size as an integer", {
  expect_identical(seq_budget_pilot(), 5L)
  expect_identical(seq_budget_pilot(m_0 = 10), 10L)
  # The article's minimum: two observations per population (p. 429).
  expect_identical(seq_budget_pilot(m_0 = 2), 2L)
})

test_that("seq_budget_pilot refuses what it used to truncate or pass on", {
  expect_error(seq_budget_pilot(m_0 = 1), "at least 2")
  expect_error(seq_budget_pilot(m_0 = 2.7), "whole number")
  expect_error(seq_budget_pilot(m_0 = Inf), "whole number")
  expect_error(seq_budget_pilot(m_0 = NA), "whole number")
  expect_error(seq_budget_pilot(m_0 = c(5, 6)), "whole number")
})


# --- seq_budget: the estimate, the interval, and the criterion -----------

test_that("seq_budget returns the documented object", {
  groups <- exact_groups(c(30, 20, 12), c(4, 4, 1), hand_s)
  out <- seq_budget(groups, hand_w, hand_a, budget = 1000)
  expect_s3_class(out, c("smsd_check", "data.frame"), exact = TRUE)
  expect_identical(out$term,
                   c("lower_limit", "contrast", "upper_limit", "variance",
                     "current_total_n", "current_total_cost"))
  expect_identical(attr(out, "method"),
                   "budget constrained K-population contrast")
  # Every attribute the help page lists, and no attribute it does not list
  # (names, row.names, and class belong to the data frame).
  documented <- c("method", "is_satisfied", "rule_met", "budget_exhausted",
                  "budget_exceeded", "variance_unusable", "current_n",
                  "criterion_n", "conf_level", "notes", "n_per_group",
                  "mean_per_group", "sd_per_group", "target_per_group",
                  "criterion_per_group", "met_per_group",
                  "not_estimable_per_group", "cost_per_group",
                  "contrast_weights", "budget", "remaining_budget")
  expect_setequal(setdiff(names(attributes(out)),
                          c("names", "row.names", "class")),
                  documented)
  expect_true(all(is.finite(out$value)))
  expect_silent(seq_budget(groups, hand_w, hand_a, budget = 1000))
})

test_that("seq_budget matches Equation 18 on a hand-computed case", {
  # Budget 1000: targets 1000 * (0.075, 0.0333..., 0.025) = (75, 33.33, 25);
  # the smallest whole numbers satisfying N >= target are (75, 34, 25).
  # n = (30, 20, 12): cost 4*30 + 9*20 + 16*12 = 492, so 508 remains.
  # Means (4, 4, 1): contrast 4 - 8 + 1 = -3.
  # Variance (Equation 5 with sample variances):
  #   1*9/30 + 4*1/20 + 1*4/12 = 0.3 + 0.2 + 0.3333 = 5/6.
  # 95% limits: -3 -/+ 1.959964 * sqrt(5/6) = -3 -/+ 1.789194.
  groups <- exact_groups(c(30, 20, 12), c(4, 4, 1), hand_s)
  out <- seq_budget(groups, hand_w, hand_a, budget = 1000)
  expect_equal(attr(out, "target_per_group"), c(75, 100 / 3, 25))
  expect_identical(attr(out, "criterion_per_group"), c(75L, 34L, 25L))
  expect_identical(attr(out, "n_per_group"), c(30L, 20L, 12L))
  expect_equal(attr(out, "mean_per_group"), c(4, 4, 1))
  expect_equal(attr(out, "sd_per_group"), hand_s)
  expect_equal(term(out, "contrast"), -3)
  expect_equal(term(out, "variance"), 5 / 6)
  expect_equal(term(out, "lower_limit"), -4.789194, tolerance = 1e-6)
  expect_equal(term(out, "upper_limit"), -1.210806, tolerance = 1e-6)
  expect_equal(term(out, "current_total_n"), 62)
  expect_equal(term(out, "current_total_cost"), 492)
  expect_equal(attr(out, "remaining_budget"), 508)
  expect_identical(attr(out, "current_n"), 62L)
  expect_equal(attr(out, "criterion_n"), 134)
  expect_identical(attr(out, "met_per_group"), c(FALSE, FALSE, FALSE))
  expect_false(attr(out, "rule_met"))
  expect_false(attr(out, "budget_exhausted"))
  expect_false(attr(out, "budget_exceeded"))
  expect_false(attr(out, "variance_unusable"))
  expect_false(attr(out, "is_satisfied"))
  expect_length(attr(out, "notes"), 0)
})

test_that("conf_level changes the limits and nothing else", {
  # 90% limits by hand: -3 -/+ 1.644854 * sqrt(5/6) = -3 -/+ 1.501539.
  groups <- exact_groups(c(30, 20, 12), c(4, 4, 1), hand_s)
  out95 <- seq_budget(groups, hand_w, hand_a, budget = 1000)
  out90 <- seq_budget(groups, hand_w, hand_a, budget = 1000,
                      conf_level = 0.90)
  expect_equal(term(out90, "lower_limit"), -4.501539, tolerance = 1e-6)
  expect_equal(term(out90, "upper_limit"), -1.498461, tolerance = 1e-6)
  expect_identical(attr(out90, "conf_level"), 0.90)
  expect_equal(out90$value[-c(1, 3)], out95$value[-c(1, 3)])
  expect_identical(attr(out90, "criterion_per_group"),
                   attr(out95, "criterion_per_group"))
})

test_that("the criteria of Equation 18 always cost the budget in total", {
  # Multiply each target by its unit cost and add: the sum of |c| s sqrt(a)
  # cancels and A0 is left. This is why, judged on the current data, every
  # group can meet its criterion only when the budget is fully spent.
  withr::local_seed(113)
  groups <- list(stats::rnorm(12, 0, 3), stats::rexp(9), stats::runif(20),
                 stats::rnorm(7, 5, 0.4))
  a <- c(38, 15, 10, 4.5)
  out <- seq_budget(groups, c(1, -0.5, -0.25, -0.25), a, budget = 7321)
  expect_equal(sum(a * attr(out, "target_per_group")), 7321)
})


# --- seq_budget: when sampling ends --------------------------------------

test_that("the rule is met only at exactly the budget (whole-number targets)", {
  # Budget 3000: targets 3000 * (0.075, 0.0333..., 0.025) = (225, 100, 75),
  # whole numbers, costing 900 + 900 + 1200 = 3000. A target that is a whole
  # number up to floating point error must not be rounded up to the next one.
  groups <- exact_groups(c(225, 100, 75), c(4, 4, 1), hand_s)
  out <- seq_budget(groups, hand_w, hand_a, budget = 3000)
  expect_identical(attr(out, "criterion_per_group"), c(225L, 100L, 75L))
  expect_identical(attr(out, "met_per_group"), c(TRUE, TRUE, TRUE))
  expect_true(attr(out, "rule_met"))
  expect_false(attr(out, "budget_exhausted"))
  expect_false(attr(out, "budget_exceeded"))
  expect_true(attr(out, "is_satisfied"))
  expect_equal(attr(out, "remaining_budget"), 0)
  expect_length(attr(out, "notes"), 0)
  expect_identical(seq_budget_next(out), NA_integer_)
})

test_that("a criterion is rounded up by the package's rule for a sample size", {
  # Two groups with equal weights, costs, and standard deviations, so each
  # target is half the budget. The package rounds a criterion up after
  # forgiving one part in 1e10 of representation error, and no more:
  #   budget 200 + 2e-13: target 100 + 1e-13, one part in 1e15 above 100
  #     -> criterion 100; 100 per group meets the rule.
  #   budget 200.0000001: target 100.00000005, five parts in 1e10 above 100
  #     -> criterion 101; 100 per group is below it, and the 1e-7 that is
  #     left cannot pay for an observation, so the budget ends the sampling.
  # (The family's own earlier allowance, 1e-9 * (1 + target), gave 100 in
  # both cases.)
  g <- exact_groups(c(100, 100), c(0, 0), c(2, 2))
  out <- seq_budget(g, c(1, -1), c(1, 1), budget = 200 + 2e-13)
  expect_gt(attr(out, "target_per_group")[1], 100)
  expect_identical(attr(out, "criterion_per_group"), c(100L, 100L))
  expect_true(attr(out, "rule_met"))
  expect_false(attr(out, "budget_exhausted"))
  out <- seq_budget(g, c(1, -1), c(1, 1), budget = 200.0000001)
  expect_equal(attr(out, "target_per_group"), rep(100.00000005, 2),
               tolerance = 1e-13)
  expect_identical(attr(out, "criterion_per_group"), c(101L, 101L))
  expect_equal(attr(out, "criterion_n"), 202)
  expect_false(attr(out, "rule_met"))
  expect_true(attr(out, "budget_exhausted"))
  expect_true(attr(out, "is_satisfied"))
  # The allowance is the one part in 1e10 that the help states, no more and
  # no less: a target half a part in 1e10 above 100 is taken to be 100 (the
  # group has met its criterion with a size a hair below its target), and a
  # target two parts in 1e10 above 100 needs 101.
  out <- seq_budget(g, c(1, -1), c(1, 1), budget = 200 * (1 + 5e-11))
  expect_gt(attr(out, "target_per_group")[1], 100 * (1 + 4e-11))
  expect_identical(attr(out, "criterion_per_group"), c(100L, 100L))
  expect_identical(attr(out, "met_per_group"), c(TRUE, TRUE))
  out <- seq_budget(g, c(1, -1), c(1, 1), budget = 200 * (1 + 2e-10))
  expect_lt(attr(out, "target_per_group")[1], 100 * (1 + 3e-10))
  expect_identical(attr(out, "criterion_per_group"), c(101L, 101L))
  expect_identical(attr(out, "met_per_group"), c(FALSE, FALSE))
})

test_that("one observation short and affordable: sampling continues", {
  # As above with group 3 at 74: cost 2984, 16 left, and group 3 costs 16.
  groups <- exact_groups(c(225, 100, 74), c(4, 4, 1), hand_s)
  out <- seq_budget(groups, hand_w, hand_a, budget = 3000)
  expect_identical(attr(out, "met_per_group"), c(TRUE, TRUE, FALSE))
  expect_false(attr(out, "rule_met"))
  expect_false(attr(out, "budget_exhausted"))
  expect_false(attr(out, "is_satisfied"))
  expect_identical(seq_budget_next(out), 3L)
})

test_that("one observation short and not affordable: stop by exhaustion", {
  # Budget 2999: targets (224.925, 99.967, 74.975), criteria (225, 100, 75).
  # n = (225, 100, 74) costs 2984; 15 is left and group 3 costs 16.
  groups <- exact_groups(c(225, 100, 74), c(4, 4, 1), hand_s)
  out <- seq_budget(groups, hand_w, hand_a, budget = 2999)
  expect_identical(attr(out, "criterion_per_group"), c(225L, 100L, 75L))
  expect_false(attr(out, "rule_met"))
  expect_true(attr(out, "budget_exhausted"))
  expect_false(attr(out, "budget_exceeded"))
  expect_true(attr(out, "is_satisfied"))
  expect_identical(seq_budget_next(out), NA_integer_)
  # The printed object shows current n 399 below criterion n 400. Its
  # stopping line does not read "TRUE" there: the rule itself (every group
  # at its criterion) was not met, and the line says what ended the sampling
  # (is_satisfied, above, stays TRUE). The note gives the amounts.
  expect_length(attr(out, "notes"), 1)
  expect_match(attr(out, "notes"), "budget cannot pay for another")
  expect_match(attr(out, "notes"), "15 of 2999 is left")
  expect_match(attr(out, "notes"), "\\(group 3\\) costs at least 16")
  lines <- format(out)
  expect_true("  Current n       : 399" %in% lines)
  expect_true("  Criterion n     : 400" %in% lines)
  expect_true(paste("  Stopping rule met: FALSE (sampling has ended: the",
                    "budget is exhausted)") %in% lines)
  expect_false(any(grepl("Stopping rule met: TRUE", lines, fixed = TRUE)))
  expect_match(paste(lines, collapse = "\n"),
               "Note: Sampling has ended because the budget")
  # One more observation in group 3 (which a budget of 3000 pays for) and
  # every group is at its criterion: that ending prints as TRUE.
  met <- seq_budget(exact_groups(c(225, 100, 75), c(4, 4, 1), hand_s), hand_w,
                    hand_a, budget = 3000)
  expect_true(attr(met, "rule_met"))
  expect_false(attr(met, "budget_exhausted"))
  expect_true("  Stopping rule met: TRUE" %in% format(met))
})

test_that("budget_exhausted looks only at groups that still need data", {
  # The article's Table 1 costs (38, 15, 10), budget 5000, weights (1, -2, 1),
  # sample SDs (4.8, 0.8, 0.8). By hand:
  #   |c| s sqrt(a) = 4.8*6.164414, 1.6*3.872983, 0.8*3.162278
  #                 = 29.589187, 6.196773, 2.529822; sum 38.315782.
  #   targets = 5000 * share / (a * sum) = 101.611, 53.910, 33.013;
  #   criteria (102, 54, 34).
  # n = (96, 62, 40) costs 3648 + 930 + 400 = 4978, so 22 is left. Only
  # group 1 is short, and it costs 38. The cheapest group overall (10) has
  # met its criterion and must not count; if it did, the procedure would be
  # left at "not satisfied" with no affordable observation.
  w <- c(1, -2, 1)
  a <- c(38, 15, 10)
  groups <- exact_groups(c(96, 62, 40), c(1, 2, 3), c(4.8, 0.8, 0.8))
  out <- seq_budget(groups, w, a, budget = 5000)
  expect_equal(attr(out, "target_per_group"), c(101.611, 53.910, 33.013),
               tolerance = 1e-5)
  expect_identical(attr(out, "criterion_per_group"), c(102L, 54L, 34L))
  expect_identical(attr(out, "met_per_group"), c(FALSE, TRUE, TRUE))
  expect_equal(attr(out, "remaining_budget"), 22)
  expect_false(attr(out, "rule_met"))
  expect_true(attr(out, "budget_exhausted"))
  expect_true(attr(out, "is_satisfied"))
  expect_identical(seq_budget_next(out), NA_integer_)
  expect_match(attr(out, "notes"), "\\(group 1\\) costs at least 38")
})

test_that("a dear group that cannot be paid for does not stop a cheaper one", {
  # Same design, n = (101, 49, 40): cost 3838 + 735 + 400 = 4973, 27 left.
  # Group 1 is short (101 < 102) but costs 38; group 2 is short (49 < 54)
  # and costs 15, so sampling continues with group 2.
  w <- c(1, -2, 1)
  a <- c(38, 15, 10)
  sds <- c(4.8, 0.8, 0.8)
  out <- seq_budget(exact_groups(c(101, 49, 40), c(1, 2, 3), sds), w, a,
                    budget = 5000)
  expect_identical(attr(out, "met_per_group"), c(FALSE, FALSE, TRUE))
  expect_equal(attr(out, "remaining_budget"), 27)
  expect_false(attr(out, "budget_exhausted"))
  expect_false(attr(out, "is_satisfied"))
  expect_identical(seq_budget_next(out), 2L)
  expect_identical(seq_budget_next(out, order = "rounds"), 2L)
  # One more in group 2: cost 4988, 12 left, group 2 still short (50 < 54)
  # but 15 > 12: exhausted.
  out <- seq_budget(exact_groups(c(101, 50, 40), c(1, 2, 3), sds), w, a,
                    budget = 5000)
  expect_true(attr(out, "budget_exhausted"))
  expect_true(attr(out, "is_satisfied"))
  expect_identical(seq_budget_next(out), NA_integer_)
  expect_match(attr(out, "notes"), "groups 1 and 2")
})

test_that("data that already cost more than the budget are flagged", {
  # 30 observations at 10 each cost 300 against a budget of 100 (a pilot the
  # budget cannot pay for). The result is flagged with a note, not reported
  # as "met" with no comment.
  groups <- exact_groups(c(10, 10, 10), c(0, 0, 0), c(1, 1, 1))
  out <- seq_budget(groups, c(1, -2, 1), c(10, 10, 10), budget = 100)
  expect_true(attr(out, "budget_exceeded"))
  expect_true(attr(out, "is_satisfied"))
  expect_equal(attr(out, "remaining_budget"), -200)
  expect_identical(seq_budget_next(out), NA_integer_)
  expect_match(attr(out, "notes")[1],
               "already cost 300, which is more than the budget of 100")
  expect_match(paste(format(out), collapse = "\n"),
               "budget constraint does not hold")
  # A pilot the budget cannot pay for is not cured by seq_budget_next(); the
  # note says what is.
  expect_match(attr(out, "notes")[1], "Raise the budget or use a smaller pilot")
  # Within budget: not flagged.
  out <- seq_budget(groups, c(1, -2, 1), c(10, 10, 10), budget = 300)
  expect_false(attr(out, "budget_exceeded"))
})

test_that("costs that are not whole numbers do not trip the budget checks", {
  # Six observations at 0.1 cost 0.6000000000000001 in floating point; that
  # is the budget, not more than the budget.
  groups <- exact_groups(c(3, 3), c(0, 0), c(1, 1))
  out <- seq_budget(groups, c(1, -1), c(0.1, 0.1), budget = 0.6)
  expect_false(attr(out, "budget_exceeded"))
  expect_true(attr(out, "rule_met"))
  # 0.5 spent of 0.6: one more at 0.1 is affordable.
  out <- seq_budget(exact_groups(c(3, 2), c(0, 0), c(1, 1)), c(1, -1),
                    c(0.1, 0.1), budget = 0.6)
  expect_false(attr(out, "is_satisfied"))
  expect_identical(seq_budget_next(out), 2L)
})

test_that("affordability carries its floating point tolerance", {
  # 3 + 3 observations at 0.1 cost 0.6000000000000001 in floating point, and
  # one more makes 0.7000000000000001: that is the budget of 0.7, not more
  # than it, so a seventh observation can be paid for. With equal costs and
  # weights the targets are 7 * s_i / (s_1 + s_2) = 2.56 and 4.44 for the
  # SDs 1.53 and 2.65, so group 2 (n = 3) is short and is named.
  expect_gt(sum(c(0.1, 0.1) * c(3, 3)) + 0.1, 0.7)
  out <- seq_budget(list(c(1, 2, 4), c(2, 3, 7)), c(1, -1), c(0.1, 0.1), 0.7)
  expect_false(attr(out, "budget_exhausted"))
  expect_false(attr(out, "is_satisfied"))
  expect_identical(seq_budget_next(out), 2L)
  # With the seventh observation the budget is spent: what is left is
  # -1.1e-16, which the note shows as 0, not as a negative amount.
  out <- seq_budget(list(c(1, 2, 4), c(2, 3, 7, 9)), c(1, -1), c(0.1, 0.1),
                    0.7)
  expect_lt(abs(attr(out, "remaining_budget")), 1e-12)
  expect_false(attr(out, "budget_exceeded"))
  expect_match(attr(out, "notes"), "0 of 0.7 is left")
  expect_match(attr(out, "notes"), "costs at least 0.1\\.")
})

test_that("amounts in the notes keep the digits of the unit costs", {
  # Budget ten million, costs 1200.4 and 2500, SDs 50 and 1, n = (6501, 878).
  # Cost 6501 * 1200.4 + 878 * 2500 = 7803800.4 + 2195000 = 9998800.4, so
  # 1199.6 is left, which is less than the 1200.4 that group 1 costs.
  # Group 1 is short: its target is (1e7 / 1200.4) * 1732.3 / 1782.3 = 8097.
  # Seven significant digits relative to the budget printed "1200 is left ...
  # costs at least 1200".
  out <- seq_budget(exact_groups(c(6501, 878), c(0, 0), c(50, 1)), c(1, -1),
                    c(1200.4, 2500), budget = 1e7)
  expect_true(attr(out, "budget_exhausted"))
  expect_equal(attr(out, "remaining_budget"), 1199.6, tolerance = 1e-9)
  expect_match(attr(out, "notes"), "1199.6 of 10000000 is left")
  expect_match(attr(out, "notes"), "\\(group 1\\) costs at least 1200.4\\.")
  # A budget that is not a whole number is shown in full, and so is a small
  # remainder. n = (4165, 2000) with SDs (1, 1) costs 4999666 + 5000000 =
  # 9999666, which leaves 0.3 of 9999666.3; group 2 is short (its target is
  # (9999666.3 / 2500) * 50 / (34.647 + 50) = 2363) and costs 2500.
  out <- seq_budget(exact_groups(c(4165, 2000), c(0, 0), c(1, 1)), c(1, -1),
                    c(1200.4, 2500), budget = 9999666.3)
  expect_true(attr(out, "budget_exhausted"))
  expect_match(attr(out, "notes"), "0.3 of 9999666.3 is left")
  expect_match(attr(out, "notes"), "\\(group 2\\) costs at least 2500\\.")
  # The display helper itself: small amounts against a large budget, and
  # floating point dust.
  expect_identical(.budget_money(0.25, 5e7), "0.25")
  expect_identical(.budget_money(-1.110223e-16, 0.7), "0")
  expect_identical(.budget_money(3000, 3000), "3000")
})


# --- seq_budget: a group whose variance is unusable ----------------------
# The package's rule: a variance estimate that cannot be used is not an
# estimate, so the stopping rule is reported as not met, the attribute
# variance_unusable is TRUE, and a note says so. Here that is a sample
# variance of zero, or one computed from observations that agree to ten
# significant digits. (A sample variance that is not finite never reaches
# that state: the scores have overflowed and the function stops, which
# "extreme magnitudes give a defined answer or a clear message" tests.)

test_that("a constant group is 'not yet estimable', never a silent 'met'", {
  # Group 1 is constant. Equation 18 as printed would give it a criterion of
  # zero and call it finished; the check instead reports the rule as not
  # met, asks for one more observation from that group, and says why.
  withr::local_seed(113)
  groups <- list(rep(3, 5), stats::rnorm(5), stats::rnorm(5))
  out <- seq_budget(groups, c(1, -0.5, -0.5), c(1, 1, 1), budget = 100)
  expect_true(attr(out, "variance_unusable"))
  expect_identical(attr(out, "not_estimable_per_group"), c(TRUE, FALSE, FALSE))
  expect_false(attr(out, "rule_met"))
  expect_false(attr(out, "is_satisfied"))
  expect_identical(attr(out, "criterion_per_group")[1], 6L)
  expect_false(attr(out, "met_per_group")[1])
  expect_true(all(is.finite(out$value)))
  expect_true(all(is.finite(attr(out, "target_per_group"))))
  # The constant group has no share; the other two share the whole budget.
  expect_identical(attr(out, "target_per_group")[1], 0)
  expect_equal(sum(attr(out, "target_per_group")[2:3]), 100)
  expect_identical(attr(out, "notes")[1], .note_variance_unusable())
  expect_match(attr(out, "notes")[2],
               "observations in group 1 are all equal")
  # With equal sizes no other group has fewer observations, so the constant
  # group is next under either order.
  expect_identical(seq_budget_next(out), 1L)
  expect_identical(seq_budget_next(out, order = "rounds"), 1L)
})

test_that("a constant group is never sampled ahead of the others", {
  # Sizes (6, 5, 5): the constant group already has more observations than
  # the others, so the next observation goes to one of them, and not to the
  # constant group every time.
  withr::local_seed(113)
  groups <- list(rep(3, 6), stats::rnorm(5), stats::rnorm(5))
  out <- seq_budget(groups, c(1, -0.5, -0.5), c(1, 1, 1), budget = 100)
  expect_true(attr(out, "variance_unusable"))
  expect_true(seq_budget_next(out) %in% 2:3)
  expect_true(seq_budget_next(out, order = "rounds") %in% 2:3)
  # Sizes (6, 6, 9): level with the smallest of the others, so it is next.
  groups <- list(rep(3, 6), stats::rnorm(6), stats::rnorm(9))
  out <- seq_budget(groups, c(1, -0.5, -0.5), c(1, 1, 1), budget = 100)
  expect_identical(seq_budget_next(out), 1L)
  expect_identical(seq_budget_next(out, order = "rounds"), 1L)
})

test_that("a group that stays constant cannot absorb the budget (hand case)", {
  # Two groups, weights (1, -1). Group 1 is constant at 5 and costs 3; group
  # 2 alternates between 1 and -1 and costs 2. Budget 50, pilot 2 each.
  # Group 1 has no share, so group 2's target is the whole budget, 50 / 2 =
  # 25, and it stays below it. The rule gives the constant group an
  # observation only when group 2 has no fewer observations than it has:
  #   (2, 2) cost 10 -> group 1 -> (3, 2) 13 -> group 2 -> (3, 3) 15 -> ...
  #   each pair of steps adds one to both and costs 5: (9, 9) costs 45;
  #   5 is left: group 1 (cost 3) -> (10, 9), 48; 2 is left: only group 2
  #   can be paid for -> (10, 10), cost 50, nothing left.
  # Sending every observation to the constant group, as before, ended at
  # (15, 2) with a cost of 49.
  loop <- function(costs, budget) {
    x <- list(c(5, 5), c(1, -1))
    path <- integer(0)
    repeat {
      out <- seq_budget(x, c(1, -1), costs, budget)
      g <- seq_budget_next(out)
      if (is.na(g)) break
      path <- c(path, g)
      x[[g]] <- c(x[[g]], if (g == 1L) 5 else (-1)^length(x[[2]]))
    }
    list(out = out, path = path)
  }
  res <- loop(c(3, 2), 50)
  expect_identical(attr(res$out, "n_per_group"), c(10L, 10L))
  expect_identical(res$path, rep(1:2, 8))
  expect_equal(term(res$out, "current_total_cost"), 50)
  expect_true(attr(res$out, "variance_unusable"))
  expect_true(attr(res$out, "budget_exhausted"))
  expect_false(attr(res$out, "is_satisfied"))
  # The constant group is taken out of turn only as a last resort. Costs
  # (2, 3), budget 49: (9, 9) costs 45; 4 is left -> group 1 -> (10, 9), 47;
  # 2 is left, group 2 (cost 3) cannot be paid for and group 1 (cost 2) can
  # -> (11, 9), cost 49.
  res <- loop(c(2, 3), 49)
  expect_identical(attr(res$out, "n_per_group"), c(11L, 9L))
  expect_equal(term(res$out, "current_total_cost"), 49)
})

test_that("a constant group with unequal costs keeps pace, no more", {
  # The article's Table 1 costs (38, 15, 10) with weights (1, -0.5, -0.5):
  # the dear group scores 10 every time on a 0 to 10 quiz. The package's
  # loop is compared, step for step, with the documented rule written out
  # here without package code: the constant group takes its turn when no
  # other group has fewer observations; otherwise the group furthest below
  # its Equation 18 target (the constant group having no share) is next.
  w <- c(1, -0.5, -0.5)
  a <- c(38, 15, 10)
  A0 <- 1500
  gen <- function(g, k) {
    if (g == 1L) rep(10, k) else pmin(10, pmax(0, round(stats::rnorm(k, 7, 1.5))))
  }
  by_package <- withr::with_seed(113, {
    x <- lapply(1:3, gen, k = 5)
    repeat {
      out <- seq_budget(x, w, a, A0)
      g <- seq_budget_next(out)
      expect_identical(is.na(g),
                       attr(out, "rule_met") || attr(out, "budget_exhausted"))
      if (is.na(g)) break
      n <- lengths(x)
      if (g == 1L) expect_lte(n[1], min(n[2:3]))
      x[[g]] <- c(x[[g]], gen(g, 1))
    }
    out
  })
  by_rule <- withr::with_seed(113, {
    x <- lapply(1:3, gen, k = 5)
    repeat {
      n <- lengths(x)
      s <- vapply(x, stats::sd, numeric(1))
      target <- A0 * abs(w) * s / (sqrt(a) * sum(abs(w) * s * sqrt(a)))
      can_pay <- sum(a * n) + a <= A0
      open <- n < target
      g <- if (can_pay[1] && n[1] <= min(n[2:3])) {
        1L
      } else if (any(open[2:3] & can_pay[2:3])) {
        ratio <- ifelse(open & can_pay & s > 0, n / target, Inf)
        which.min(ratio)
      } else if (can_pay[1]) {
        1L
      } else {
        NA_integer_
      }
      if (is.na(g)) break
      x[[g]] <- c(x[[g]], gen(g, 1))
    }
    lengths(x)
  })
  final_n <- attr(by_package, "n_per_group")
  expect_identical(final_n, by_rule)
  # 24, 23, and 24 observations at 38, 15, and 10 cost 1497. The earlier
  # rule ended at (36, 5, 5): 1368 of the 1500 on the constant group.
  expect_identical(final_n, c(24L, 23L, 24L))
  expect_lte(final_n[1], min(final_n[2:3]) + 1L)
  expect_equal(term(by_package, "current_total_cost"), 1497)
  # The study ends by the budget with the variance still not estimable:
  # sampling has ended, the rule is reported as not met, and no decision is
  # made from it.
  expect_true(attr(by_package, "budget_exhausted"))
  expect_true(attr(by_package, "variance_unusable"))
  expect_false(attr(by_package, "is_satisfied"))
  expect_error(decide_seq_budget_noninferiority(by_package, margin = 1),
               "No decision can be made")
})

test_that("a group constant up to floating point error is not estimable", {
  # 0.1 + 0.2 is not 0.3 in floating point, so this group's standard
  # deviation is about 2.5e-17 rather than 0. It carries no information
  # about the variance all the same.
  x <- c(0.3, 0.1 + 0.2, 0.3, 0.3, 0.3)
  expect_gt(stats::sd(x), 0)
  out <- seq_budget(list(x, c(1, 2, 4, 3, 5)), c(1, -1), c(1, 1), 100)
  expect_true(attr(out, "variance_unusable"))
  expect_identical(attr(out, "not_estimable_per_group"), c(TRUE, FALSE))
  expect_false(attr(out, "is_satisfied"))
  expect_identical(seq_budget_next(out), 1L)
  expect_match(attr(out, "notes")[2],
               "differ too little relative to their size")
})

test_that("'no variation' is the package's one tolerance: ten digits", {
  # The package calls a sample constant when its range is no more than 1e-10
  # of its largest absolute value (agreement to ten significant digits).
  # Scores recorded near one million:
  #   differing by 1e-5: range 4e-5, 4e-11 of a million -> constant;
  #   differing by 1e-3: range 4e-3, 4e-9 of a million  -> not constant.
  # (The family's own earlier threshold, a standard deviation below
  # 100 * .Machine$double.eps of the largest value, called both usable.)
  other <- c(1, 2, 4, 3, 5)
  eleventh_digit <- 1e6 + (0:4) * 1e-5
  ninth_digit <- 1e6 + (0:4) * 1e-3
  expect_gt(stats::sd(eleventh_digit), 1e-5)
  out <- seq_budget(list(eleventh_digit, other), c(1, -1), c(1, 1), 100)
  expect_true(attr(out, "variance_unusable"))
  expect_identical(attr(out, "not_estimable_per_group"), c(TRUE, FALSE))
  expect_identical(attr(out, "criterion_per_group")[1], 6L)
  expect_false(attr(out, "is_satisfied"))
  expect_identical(attr(out, "notes")[1], .note_variance_unusable())
  out <- seq_budget(list(ninth_digit, other), c(1, -1), c(1, 1), 100)
  expect_false(attr(out, "variance_unusable"))
  expect_identical(attr(out, "not_estimable_per_group"), c(FALSE, FALSE))
  # Subtracting the same constant from every score in every group, as the
  # help advises, restores a usable variance and leaves a contrast whose
  # weights sum to zero where it was: (1e6 + 2e-5) - 3 = 999997.00002.
  flat <- seq_budget(list(eleventh_digit, other), c(1, -1), c(1, 1), 100)
  out <- seq_budget(list(eleventh_digit - 1e6, other - 1e6), c(1, -1),
                    c(1, 1), 100)
  expect_false(attr(out, "variance_unusable"))
  expect_equal(term(out, "contrast"), 999997.00002, tolerance = 1e-12)
  expect_equal(term(out, "contrast"), term(flat, "contrast"),
               tolerance = 1e-12)
})

test_that("a group that varies still varies when an observation is added", {
  # The simulation engine relies on this to look again only at a group that
  # has been constant so far. Random groups whose range sits just above the
  # tolerance (1 to 50 times 1e-10 of the largest value), about any location
  # and of either sign, with a new observation of any size.
  withr::local_seed(113)
  broke <- 0L
  tried <- 0L
  for (i in 1:400) {
    location <- sample(c(-1, 1), 1) * 10^stats::runif(1, -3, 6)
    x <- location * (1 + stats::runif(4, 0, 50) * 1e-10)
    if (.is_constant(x)) next
    tried <- tried + 1L
    new <- switch(sample(4, 1),
                  location,                                # inside the range
                  location * (1 + stats::rnorm(1) * 1e-9), # just outside it
                  location * 10^stats::runif(1, 1, 8),     # far larger
                  -location * 10^stats::runif(1, -2, 8))   # other sign
    broke <- broke + .is_constant(c(x, new))
  }
  expect_gt(tried, 300L)
  expect_identical(broke, 0L)
})

test_that("all groups constant: defined values, not NaN, and not 'met'", {
  groups <- list(rep(3, 5), rep(2, 5), rep(1, 5))
  out <- seq_budget(groups, c(1, -0.5, -0.5), c(1, 1, 1), budget = 100)
  expect_true(attr(out, "variance_unusable"))
  expect_false(attr(out, "is_satisfied"))
  expect_true(all(is.finite(out$value)))
  expect_equal(term(out, "contrast"), 3 - 1 - 0.5)
  expect_identical(attr(out, "criterion_per_group"), c(6L, 6L, 6L))
  expect_identical(seq_budget_next(out), 1L)
  # One more in group 1: the fewest observations come first, in rounds.
  out <- seq_budget(list(rep(3, 6), rep(2, 5), rep(1, 5)), c(1, -0.5, -0.5),
                    c(1, 1, 1), budget = 100)
  expect_identical(seq_budget_next(out), 2L)
})

test_that("no budget left and a variance not estimable: ended, not satisfied", {
  # The rule for an unusable variance holds to the end. The budget (15) is
  # spent, so sampling has ended: budget_exhausted is TRUE and there is no
  # group to name. But the variance is unusable, so the rule is not reported
  # as met; the one-row summaries, which have no column for the flag, show
  # FALSE too.
  groups <- list(rep(3, 5), rep(2, 5), rep(1, 5))
  out <- seq_budget(groups, c(1, -0.5, -0.5), c(1, 1, 1), budget = 15)
  expect_true(attr(out, "variance_unusable"))
  expect_false(attr(out, "rule_met"))
  expect_true(attr(out, "budget_exhausted"))
  expect_false(attr(out, "is_satisfied"))
  expect_identical(seq_budget_next(out), NA_integer_)
  expect_false(generics::glance(out)$is_satisfied)
  expect_false(summary(out)$is_satisfied)
  printed <- paste(format(out), collapse = "\n")
  expect_match(printed, "Stopping rule met: FALSE")
  expect_match(attr(out, "notes")[1], "Sampling has ended because the budget")
  # The second note says what the package's standard note says (not usable,
  # reported as not met), without its advice to collect more data, which the
  # budget no longer allows.
  expect_match(attr(out, "notes")[2], "variance estimate is not usable")
  expect_match(attr(out, "notes")[2], "reported as not met")
  expect_match(attr(out, "notes")[2], "too narrow to rely on")
  expect_match(attr(out, "notes")[2], "would take a larger budget")
  expect_length(attr(out, "notes"), 2)
  expect_false(.note_variance_unusable() %in% attr(out, "notes"))
  # The decision functions refuse it, all four.
  refusal <- "variance estimate is not usable.*variance of the contrast is not estimable"
  expect_error(decide_seq_budget_noninferiority(out, margin = 1), refusal)
  expect_error(decide_seq_budget_superiority(out, margin = 0), refusal)
  expect_error(decide_seq_budget_inferiority(out, margin = 0), refusal)
  expect_error(decide_seq_budget_equivalence(out, -1, 1), refusal)
  # The same when the data already cost more than the budget.
  out <- seq_budget(groups, c(1, -0.5, -0.5), c(1, 1, 1), budget = 10)
  expect_true(attr(out, "budget_exceeded"))
  expect_false(attr(out, "is_satisfied"))
  expect_identical(seq_budget_next(out), NA_integer_)
  # Two groups, one constant.
  out <- seq_budget(list(rep(0, 5), c(1, 2, 3, 4, 5)), c(1, -1), c(1, 1), 10)
  expect_true(attr(out, "budget_exhausted"))
  expect_false(attr(out, "is_satisfied"))
  expect_match(attr(out, "notes")[2], "observations in group 1 are all equal")
})

test_that("a constant group with a zero weight does not matter", {
  withr::local_seed(113)
  groups <- list(rep(3, 5), stats::rnorm(5), stats::rnorm(5))
  out <- seq_budget(groups, c(0, 1, -1), c(1, 1, 1), budget = 100)
  expect_false(attr(out, "variance_unusable"))
  expect_identical(attr(out, "criterion_per_group")[1], 0L)
  expect_true(attr(out, "met_per_group")[1])
})


# --- seq_budget: the first passage rule (closed) -------------------------

test_that("closed keeps a group closed, as in the article's Stage 1 and 2", {
  # Hand design at budget 1000, n = (30, 20, 12): every group is below its
  # criterion (75, 34, 25) on the current data.
  groups <- exact_groups(c(30, 20, 12), c(4, 4, 1), hand_s)
  out <- seq_budget(groups, hand_w, hand_a, budget = 1000,
                    closed = c(TRUE, FALSE, FALSE))
  expect_identical(attr(out, "met_per_group"), c(TRUE, FALSE, FALSE))
  expect_false(attr(out, "is_satisfied"))
  # Open groups 2 and 3: n / target = 20 / 33.33 = 0.60 and 12 / 25 = 0.48.
  expect_identical(seq_budget_next(out), 3L)
  # All closed: the rule is met with 508 of the budget unspent, and the
  # printed object says why the sizes are below the criteria.
  out <- seq_budget(groups, hand_w, hand_a, budget = 1000,
                    closed = c(TRUE, TRUE, TRUE))
  expect_true(attr(out, "rule_met"))
  expect_true(attr(out, "is_satisfied"))
  expect_false(attr(out, "budget_exhausted"))
  expect_equal(attr(out, "remaining_budget"), 508)
  expect_match(attr(out, "notes"), "first passage rule")
  expect_identical(seq_budget_next(out), NA_integer_)
})

test_that("a closed group whose variance is not estimable is reopened", {
  # A constant group cannot have met a criterion, whatever 'closed' says.
  out <- seq_budget(list(rep(3, 5), c(1, 2, 4, 3, 5)), c(1, -1), c(1, 1), 100,
                    closed = c(TRUE, FALSE))
  expect_false(attr(out, "met_per_group")[1])
  expect_true(attr(out, "variance_unusable"))
  expect_identical(seq_budget_next(out, order = "rounds"), 1L)
  # A closed group that can be estimated stays closed (n = 5 against a
  # target of 100 * s_1 / (s_1 + s_2), far above 5).
  out <- seq_budget(list(c(1, 3, 2, 5, 4), c(1, 2, 4, 3, 5)), c(1, -1),
                    c(1, 1), 100, closed = c(TRUE, FALSE))
  expect_true(attr(out, "met_per_group")[1])
  expect_identical(seq_budget_next(out, order = "rounds"), 2L)
})

test_that("closed is validated", {
  groups <- exact_groups(c(30, 20, 12), c(4, 4, 1), hand_s)
  expect_error(seq_budget(groups, hand_w, hand_a, 1000, closed = c(TRUE, NA, FALSE)),
               "'closed' must be NULL or a logical vector")
  expect_error(seq_budget(groups, hand_w, hand_a, 1000, closed = c(TRUE, FALSE)),
               "'closed' must be NULL or a logical vector")
  expect_error(seq_budget(groups, hand_w, hand_a, 1000, closed = c(1, 0, 0)),
               "'closed' must be NULL or a logical vector")
})


# --- seq_budget: argument checks -----------------------------------------

test_that("seq_budget rejects bad data with a message that names the group", {
  good <- exact_group(10, 0, 1)
  expect_error(seq_budget(list(c(good, NA), good), c(1, -1), c(1, 1), 100),
               "'groups\\[\\[1\\]\\]' contains missing values")
  expect_error(seq_budget(list(good, c(good, NaN)), c(1, -1), c(1, 1), 100),
               "'groups\\[\\[2\\]\\]' contains missing values")
  expect_error(seq_budget(list(good, c(good, Inf)), c(1, -1), c(1, 1), 100),
               "'groups\\[\\[2\\]\\]' contains infinite values")
  expect_error(seq_budget(list(good, 1), c(1, -1), c(1, 1), 100),
               "'groups\\[\\[2\\]\\]' needs at least 2 observations")
  expect_error(seq_budget(list(good, "a"), c(1, -1), c(1, 1), 100),
               "'groups\\[\\[2\\]\\]' must be a numeric vector")
  expect_error(seq_budget(list(good, matrix(1:6, 2)), c(1, -1), c(1, 1), 100),
               "'groups\\[\\[2\\]\\]' must be a numeric vector, not a matrix")
  expect_error(seq_budget(good, c(1, -1), c(1, 1), 100),
               "'groups' must be a list")
  expect_error(seq_budget(data.frame(a = good, b = good), c(1, -1), c(1, 1),
                          100),
               "'groups' must be a list")
  expect_error(seq_budget(list(good), 1, 1, 100),
               "at least 2 populations; 'groups' describes 1")
})

test_that("seq_budget rejects bad weights, costs, budget, and conf_level", {
  g <- exact_groups(c(10, 10, 10), c(0, 0, 0), c(1, 1, 1))
  w <- c(1, -2, 1)
  a <- c(1, 1, 1)
  expect_error(seq_budget(g, c(0, 0, 0), a, 100), "at least one nonzero")
  expect_error(seq_budget(g, c(1, -1), a, 100),
               "'contrast_weights' must have 3 entries \\(one per group in 'groups'\\); it has 2")
  expect_error(seq_budget(g, c(1, NA, -1), a, 100),
               "'contrast_weights' must contain only finite numbers.*entry 2 is NA")
  expect_error(seq_budget(g, c(1, Inf, -1), a, 100),
               "'contrast_weights' must contain only finite numbers.*entry 2 is Inf")
  expect_error(seq_budget(g, "a", a, 100),
               "'contrast_weights' must be a numeric vector")
  expect_error(seq_budget(g, w, c(1, 1), 100),
               "'cost_per_group' must have 3 entries \\(one per group in 'groups'\\); it has 2")
  expect_error(seq_budget(g, w, c(1, NA, 1), 100),
               "'cost_per_group' must contain only finite numbers.*entry 2 is NA")
  expect_error(seq_budget(g, w, c(1, 0, 1), 100),
               "'cost_per_group' must contain only positive numbers \\(the cost of one observation in each group\\); entry 2 is 0")
  expect_error(seq_budget(g, w, c(1, -3, 1), 100),
               "'cost_per_group' must contain only positive numbers.*entry 2 is -3")
  expect_error(seq_budget(g, w, a, NA), "'budget'")
  expect_error(seq_budget(g, w, a, 0), "'budget'.*greater than 0")
  expect_error(seq_budget(g, w, a, -5), "'budget'.*greater than 0")
  expect_error(seq_budget(g, w, a, Inf), "'budget'")
  expect_error(seq_budget(g, w, a, c(100, 200)), "'budget'")
  # The package's rule for a sample size that cannot be represented, with
  # advice that names the two arguments that set it here.
  expect_error(seq_budget(g, w, a, 1e300),
               "number of observations that the budget pays for in the cheapest group \\(1e\\+300\\) is too large to represent. Check 'budget' and 'cost_per_group'")
  expect_error(seq_budget(g, w, c(1e-9, 1, 1), 100), "too large to represent")
  expect_error(seq_budget(g, w, a, 100, conf_level = NA), "'conf_level'")
  expect_error(seq_budget(g, w, a, 100, conf_level = 1), "'conf_level'")
  expect_error(seq_budget(g, w, a, 100, conf_level = 95), "'conf_level'")
})

test_that("seq_budget refuses names that disagree with the groups", {
  # A list from split() is in alphabetical order; weights and costs typed in
  # the order of the design are not. Matched by position they would give
  # another contrast (here 0.5 * 7 - 3 + 0.5 * 4 = 2.5 instead of -3.5)
  # without a word.
  groups <- list(face = c(6, 7, 8), hybrid = c(2, 3, 4), video = c(3, 4, 5))
  w <- c(video = 0.5, face = -1, hybrid = 0.5)
  costs <- c(video = 150, face = 500, hybrid = 250)
  expect_error(seq_budget(groups, w, costs, 50000),
               "'contrast_weights' is named \\(video, face, hybrid\\) but 'groups' is named \\(face, hybrid, video\\)")
  expect_error(seq_budget(groups, unname(w), costs, 50000),
               "'cost_per_group' is named \\(video, face, hybrid\\) but 'groups' is named")
  expect_error(seq_budget(unname(groups), w, costs[c(2, 3, 1)], 50000),
               "'cost_per_group' is named \\(face, hybrid, video\\) but 'contrast_weights' is named")
  # In the order of the groups the names are accepted, and the contrast is
  # the intended one: -7 + 0.5 * 3 + 0.5 * 4 = -3.5.
  out <- seq_budget(groups, w[names(groups)], costs[names(groups)], 50000)
  expect_equal(term(out, "contrast"), -3.5)
  expect_equal(term(out, "current_total_cost"), 3 * (500 + 250 + 150))
  # Names on one argument only, or no names, are matched by position as
  # documented.
  expect_equal(term(seq_budget(groups, c(-1, 0.5, 0.5), c(500, 250, 150),
                               50000), "contrast"), -3.5)
  expect_equal(term(seq_budget(unname(groups), w[names(groups)],
                               c(500, 250, 150), 50000), "contrast"), -3.5)
})

test_that("extreme magnitudes give a defined answer or a clear message", {
  withr::local_seed(113)
  g <- list(stats::rnorm(10), stats::rnorm(10), stats::rnorm(10))
  w <- c(1, -0.5, -0.5)
  # The unit of money does not matter: costs of 1e-300 or 1e210 with a budget
  # of 100 observations' worth give the same targets as costs of 1, and a
  # definite stopping flag (never NA, which the check object refuses).
  ref <- seq_budget(g, w, rep(1, 3), 100)
  for (unit in c(1e-300, 1e210)) {
    out <- seq_budget(g, w, rep(unit, 3), 100 * unit)
    expect_false(attr(out, "is_satisfied"))
    expect_equal(attr(out, "target_per_group"), attr(ref, "target_per_group"))
    expect_identical(seq_budget_next(out), seq_budget_next(ref))
  }
  # Values too large to combine: the contrast overflows.
  expect_error(seq_budget(list(c(1e308, -1e308, 1e308), c(1, 2, 3)), c(1, -1),
                          c(1, 1), 100),
               "contrast estimate or its variance cannot be computed.*'groups' and 'contrast_weights'.*Rescale")
  # A sample variance that is not finite is this case, not the "variance
  # unusable" state of a check: scores of the order of 1e160 have a finite
  # mean, but their squared deviations overflow, so the standard deviation
  # is Inf and no variance of the contrast exists to report.
  huge <- c(1e160, -1e160, 3e159, 2e159)
  expect_true(is.finite(mean(huge)))
  expect_identical(stats::sd(huge), Inf)
  expect_error(seq_budget(list(huge, c(1, 2, 3, 4)), c(1, -1), c(1, 1), 100),
               "contrast estimate or its variance cannot be computed.*Rescale")
  # The same with that group's weight at zero: nothing is returned for data
  # the function cannot summarize.
  expect_error(seq_budget(list(huge, c(1, 2, 3, 4), c(2, 4, 5, 9)),
                          c(0, 1, -1), c(1, 1, 1), 100),
               "contrast estimate or its variance cannot be computed.*Rescale")
  # Weights so small that their squares underflow: the variance would be
  # reported as exactly 0 with no flag.
  expect_error(seq_budget(g, w * 1e-170, c(38, 15, 10), 1500),
               "contrast estimate or its variance cannot be computed")
  # Weights and data so small that the products |c| s sqrt(a) underflow to
  # zero (1e-170 * 1e-155 is below the smallest double), or so large that
  # they overflow: no share can be computed.
  expect_error(seq_budget(lapply(g, function(x) x * 1e-155), w * 1e-170,
                          c(38, 15, 10), 1500),
               "share of the budget.*cannot be computed.*Rescale")
  expect_error(seq_budget(lapply(g, function(x) x * 1e150), w * 1e170,
                          c(38, 15, 10), 1500),
               "share of the budget.*cannot be computed.*Rescale")
  # Data so small that the sample variances themselves underflow to zero are
  # in the "not yet estimable" state, never "met".
  out <- seq_budget(lapply(g, function(x) x * 1e-200), w, c(38, 15, 10), 1500)
  expect_true(attr(out, "variance_unusable"))
  expect_false(attr(out, "is_satisfied"))
  # A confidence level within rounding of 1 still gives finite limits (the
  # upper tail quantile is taken directly, not as 1 minus a tiny number).
  out <- seq_budget(g, w, c(38, 15, 10), 1500, conf_level = 1 - 1e-16)
  expect_true(all(is.finite(out$value)))
  # qnorm(5.55e-17, lower.tail = FALSE) is 8.29; the limits are symmetric.
  half <- (term(out, "upper_limit") - term(out, "lower_limit")) / 2
  expect_equal(half / sqrt(term(out, "variance")), 8.29, tolerance = 1e-3)
})

test_that("integer data are handled as doubles", {
  g_int <- list(c(1L, 5L, 9L, 3L), c(2L, 2L, 7L, 4L))
  g_dbl <- lapply(g_int, as.numeric)
  expect_equal(seq_budget(g_int, c(1, -1), c(2, 3), 500)$value,
               seq_budget(g_dbl, c(1, -1), c(2, 3), 500)$value)
})


# --- seq_budget_next -----------------------------------------------------

test_that("seq_budget_next follows the documented order", {
  # Hand design, budget 1000, n = (30, 20, 12), targets (75, 33.33, 25):
  # n / target = 0.40, 0.60, 0.48, so the group furthest below its estimated
  # optimal size is group 1; the group with the fewest observations is 3.
  groups <- exact_groups(c(30, 20, 12), c(4, 4, 1), hand_s)
  out <- seq_budget(groups, hand_w, hand_a, budget = 1000)
  expect_identical(seq_budget_next(out), 1L)
  expect_identical(seq_budget_next(out, order = "shortfall"), 1L)
  expect_identical(seq_budget_next(out, order = "rounds"), 3L)
  # The choice goes through the package's shared check: the message names
  # the argument, and an unambiguous abbreviation is accepted, as for every
  # such argument in the package.
  expect_identical(seq_budget_next(out, order = "round"), 3L)
  expect_identical(seq_budget_next(out, order = "s"), 1L)
  expect_error(seq_budget_next(out, order = "largest"),
               "'order' must be one of \"shortfall\", \"rounds\"; got \"largest\"")
  expect_error(seq_budget_next(out, order = 2), "'order' must be one of")
  expect_error(seq_budget_next(out, order = c("rounds", "shortfall")),
               "'order' must be one of")
  expect_error(seq_budget_next(out, order = NA_character_),
               "'order' must be one of")
})

test_that("seq_budget_next breaks ties by the larger SD, then by position", {
  # Equal sizes, so "rounds" ties on n: the article's order, the larger
  # standard deviation first (group 1, s = 3). Budget 1000 as above.
  groups <- exact_groups(c(12, 12, 12), c(4, 4, 1), hand_s)
  out <- seq_budget(groups, hand_w, hand_a, budget = 1000)
  expect_identical(seq_budget_next(out, order = "rounds"), 1L)
  # The same with the larger standard deviation in the LAST position, so
  # that the tie is seen to go to the larger SD and not to the first group.
  groups <- exact_groups(c(12, 12, 12), c(4, 4, 1), c(1, 2, 3))
  out <- seq_budget(groups, hand_w, hand_a, budget = 1000)
  expect_identical(seq_budget_next(out, order = "rounds"), 3L)
  # Fully symmetric groups tie on everything; the first is returned.
  groups <- exact_groups(c(12, 12), c(0, 0), c(2, 2))
  out <- seq_budget(groups, c(1, -1), c(5, 5), budget = 1000)
  expect_identical(seq_budget_next(out), 1L)
  expect_identical(seq_budget_next(out, order = "rounds"), 1L)
})

test_that("seq_budget_next accepts only a seq_budget() check", {
  withr::local_seed(113)
  smd_check <- seq_aipe_smd(stats::rnorm(40, 0.5), stats::rnorm(40),
                            width = 0.5)
  expect_error(seq_budget_next(smd_check),
               "must be an object returned by seq_budget\\(\\)")
  expect_error(seq_budget_next(data.frame(term = "contrast", value = 1)),
               "must be an object returned by seq_budget\\(\\)")
  expect_error(seq_budget_next(3), "must be an object returned by seq_budget")
})

test_that("a loop on seq_budget_next never passes the budget (unequal costs)", {
  # The article's Table 1 costs (38, 15, 10) and its normal row. At every
  # step: the group named is below its criterion and affordable; NA is
  # returned exactly when is_satisfied is TRUE; the cost stays in budget.
  withr::local_seed(113)
  mu <- c(1, 2, 3)
  sigma <- c(4.5, 0.9, 1)
  w <- c(1, -2, 1)
  a <- c(38, 15, 10)
  for (order in c("shortfall", "rounds")) {
    groups <- lapply(1:3, function(i) stats::rnorm(10, mu[i], sigma[i]))
    steps <- 0
    repeat {
      out <- seq_budget(groups, w, a, budget = 1500)
      g <- seq_budget_next(out, order = order)
      # With every variance estimable, "no group to name" and is_satisfied
      # are the same thing.
      expect_identical(is.na(g), attr(out, "is_satisfied"))
      expect_identical(attr(out, "is_satisfied"),
                       attr(out, "rule_met") || attr(out, "budget_exhausted"))
      expect_lte(term(out, "current_total_cost"), 1500)
      if (is.na(g)) break
      expect_false(attr(out, "met_per_group")[g])
      expect_lte(a[g], attr(out, "remaining_budget"))
      groups[[g]] <- c(groups[[g]], stats::rnorm(1, mu[g], sigma[g]))
      steps <- steps + 1
    }
    expect_gt(steps, 20)
    expect_true(attr(out, "budget_exhausted"))
    # Nothing that is still needed is affordable.
    short <- !attr(out, "met_per_group")
    expect_true(all(a[short] > attr(out, "remaining_budget")))
  }
})


test_that("coarse data: NA exactly when sampling has ended, never a silent 'met'", {
  # Ratings on a 1 to 3 scale and a pilot of 2, so that groups are often
  # constant at first. At every step of every study: seq_budget_next()
  # returns NA exactly when rule_met or budget_exhausted is TRUE; is_satisfied
  # is that and a usable variance; rule_met is never TRUE with an unusable
  # one; the named group is short and affordable;
  # the cost stays within the budget. Violations are counted over all steps
  # (about 500) and each count must be zero.
  withr::local_seed(113)
  a <- c(7, 3, 5)
  w <- c(1, -0.5, -0.5)
  steps <- 0L
  unusable_steps <- 0L
  bad <- c(na_iff_ended = 0L, satisfied = 0L, rule_met_unusable = 0L,
           over_budget = 0L, named_met = 0L, named_unaffordable = 0L)
  for (study in 1:15) {
    order <- if (study %% 2L) "shortfall" else "rounds"
    groups <- lapply(1:3, function(i) sample(1:3, 2, replace = TRUE))
    repeat {
      out <- seq_budget(groups, w, a, budget = 150)
      g <- seq_budget_next(out, order = order)
      ended <- attr(out, "rule_met") || attr(out, "budget_exhausted")
      unusable <- attr(out, "variance_unusable")
      steps <- steps + 1L
      unusable_steps <- unusable_steps + unusable
      bad["na_iff_ended"] <- bad["na_iff_ended"] + (is.na(g) != ended)
      bad["satisfied"] <- bad["satisfied"] +
        (attr(out, "is_satisfied") != (ended && !unusable))
      bad["rule_met_unusable"] <- bad["rule_met_unusable"] +
        (unusable && attr(out, "rule_met"))
      bad["over_budget"] <- bad["over_budget"] +
        (term(out, "current_total_cost") > 150)
      if (is.na(g)) break
      bad["named_met"] <- bad["named_met"] + attr(out, "met_per_group")[g]
      bad["named_unaffordable"] <- bad["named_unaffordable"] +
        (a[g] > attr(out, "remaining_budget"))
      groups[[g]] <- c(groups[[g]], sample(1:3, 1))
    }
  }
  expect_identical(unname(bad), rep(0L, 6), label = paste(names(bad)[bad > 0],
                                                         collapse = ", "))
  # The design does visit the state it is meant to test.
  expect_gt(steps, 300L)
  expect_gt(unusable_steps, 5L)
})


# --- the decision functions ----------------------------------------------

# A completed check with hand-computable numbers: the whole-number design at
# budget 3000, n = (225, 100, 75), means (4, 4, 1).
#   contrast = 4 - 8 + 1 = -3
#   variance = 9/225 + 4*1/100 + 4/75 = 0.04 + 0.04 + 0.05333 = 2/15
#   standard error = sqrt(2/15) = 0.3651484
#   one-sided 0.05 : 1.644854 * SE = 0.600616, limits (-3.600616, -2.399384)
#   one-sided 0.025: 1.959964 * SE = 0.715678, limits (-3.715678, -2.284322)
#   one-sided 0.10 : 1.281552 * SE = 0.467957, limits (-3.467957, -2.532043)
done_check <- function(means = c(4, 4, 1), conf_level = 0.95,
                       weights = hand_w) {
  seq_budget(exact_groups(c(225, 100, 75), means, hand_s), weights, hand_a,
             budget = 3000, conf_level = conf_level)
}
decisions <- list(
  noninferiority = decide_seq_budget_noninferiority,
  superiority = decide_seq_budget_superiority,
  inferiority = decide_seq_budget_inferiority,
  equivalence = decide_seq_budget_equivalence
)

test_that("the four decisions are one family: names and default level", {
  expect_identical(names(formals(decide_seq_budget_noninferiority)),
                   c("check", "margin", "alpha_level"))
  expect_identical(names(formals(decide_seq_budget_superiority)),
                   c("check", "margin", "alpha_level"))
  expect_identical(names(formals(decide_seq_budget_inferiority)),
                   c("check", "margin", "alpha_level"))
  expect_identical(names(formals(decide_seq_budget_equivalence)),
                   c("check", "margin_lower", "margin_upper", "alpha_level"))
  # One default level, the 0.05 of the article's tables.
  for (f in decisions) {
    expect_identical(formals(f)$alpha_level, 0.05)
  }
})

test_that("noninferiority is the article's Equation C3, in its direction", {
  # H0 (inferiority): psi <= -margin. Reject at one-sided level alpha when
  # psi_hat > -margin + z_alpha sqrt(V), z_alpha the upper alpha point.
  out <- done_check()
  expect_true(attr(out, "is_satisfied"))
  # alpha = 0.05: reject iff -3.600616 > -margin, that is, margin > 3.600616.
  expect_true(decide_seq_budget_noninferiority(out, margin = 3.61))
  expect_false(decide_seq_budget_noninferiority(out, margin = 3.59))
  # alpha = 0.025: reject iff margin > 3.715678.
  expect_false(decide_seq_budget_noninferiority(out, margin = 3.71,
                                                alpha_level = 0.025))
  expect_true(decide_seq_budget_noninferiority(out, margin = 3.72,
                                               alpha_level = 0.025))
  # Direction: a contrast far above -margin is noninferior, one far below
  # is not, whatever the upper limit does.
  expect_true(decide_seq_budget_noninferiority(
    done_check(means = c(10, 4, 1)), margin = 0.5))      # contrast +3
  expect_false(decide_seq_budget_noninferiority(
    done_check(means = c(-3, 4, 1)), margin = 0.5))      # contrast -10
  expect_type(decide_seq_budget_noninferiority(out, margin = 3.61), "logical")
})

test_that("each decision is a one-sided test at alpha_level, 0.05 by default", {
  # By hand from the estimate and its standard error (see done_check above).
  # Superiority rejects when the lower limit is above the margin (the
  # article's Equation C13 with the margin as psi_1, and C4 when the margin
  # is |psi_0|), inferiority when the upper limit is below it (Equation C14
  # with the margin as psi_2), equivalence when both hold for its two
  # margins (C15).
  out <- done_check()
  # Default level 0.05: limits (-3.600616, -2.399384).
  expect_true(decide_seq_budget_superiority(out, margin = -3.61))
  expect_false(decide_seq_budget_superiority(out, margin = -3.60))
  expect_true(decide_seq_budget_inferiority(out, margin = -2.39))
  expect_false(decide_seq_budget_inferiority(out, margin = -2.40))
  expect_true(decide_seq_budget_equivalence(out, margin_lower = -3.61,
                                            margin_upper = -2.39))
  expect_false(decide_seq_budget_equivalence(out, margin_lower = -3.60,
                                             margin_upper = -2.39))
  expect_false(decide_seq_budget_equivalence(out, margin_lower = -3.61,
                                             margin_upper = -2.40))
  # alpha_level = 0.025: limits (-3.715678, -2.284322).
  expect_true(decide_seq_budget_superiority(out, margin = -3.72,
                                            alpha_level = 0.025))
  expect_false(decide_seq_budget_superiority(out, margin = -3.71,
                                             alpha_level = 0.025))
  expect_true(decide_seq_budget_inferiority(out, margin = -2.28,
                                            alpha_level = 0.025))
  expect_false(decide_seq_budget_inferiority(out, margin = -2.29,
                                             alpha_level = 0.025))
  expect_true(decide_seq_budget_equivalence(out, -3.72, -2.28,
                                            alpha_level = 0.025))
  expect_false(decide_seq_budget_equivalence(out, -3.71, -2.28,
                                             alpha_level = 0.025))
  expect_false(decide_seq_budget_equivalence(out, -3.72, -2.29,
                                             alpha_level = 0.025))
  # alpha_level = 0.10: limits (-3.467957, -2.532043).
  expect_true(decide_seq_budget_superiority(out, margin = -3.47,
                                            alpha_level = 0.10))
  expect_false(decide_seq_budget_superiority(out, margin = -3.46,
                                             alpha_level = 0.10))
  expect_true(decide_seq_budget_inferiority(out, margin = -2.53,
                                            alpha_level = 0.10))
  expect_false(decide_seq_budget_inferiority(out, margin = -2.54,
                                             alpha_level = 0.10))
  # Superiority in the article's own sense, a positive margin (Equations 8
  # and 9): contrast +3 with the same standard error, lower limit 2.399384.
  better <- done_check(means = c(10, 4, 1))
  expect_true(decide_seq_budget_superiority(better, margin = 2.39))
  expect_false(decide_seq_budget_superiority(better, margin = 2.40))
  # The article's footnote 1 (a smaller response is better): H_A is
  # psi < -|psi_0|. Contrast -3, upper limit -2.399384.
  expect_true(decide_seq_budget_inferiority(out, margin = -2.39))
  for (f in decisions[2:3]) expect_type(f(out, margin = 0), "logical")
})

test_that("no decision reads the interval stored on the check", {
  # Checks of the same data at four confidence levels store four different
  # intervals; every decision is the same, the one computed at alpha_level
  # from the estimate and its standard error.
  for (cl in c(0.80, 0.90, 0.95, 0.99)) {
    out <- done_check(conf_level = cl)
    expect_true(decide_seq_budget_noninferiority(out, margin = 3.61))
    expect_false(decide_seq_budget_noninferiority(out, margin = 3.59))
    expect_true(decide_seq_budget_superiority(out, margin = -3.61))
    expect_false(decide_seq_budget_superiority(out, margin = -3.60))
    expect_true(decide_seq_budget_inferiority(out, margin = -2.39))
    expect_false(decide_seq_budget_inferiority(out, margin = -2.40))
    expect_true(decide_seq_budget_equivalence(out, -3.61, -2.39))
    expect_false(decide_seq_budget_equivalence(out, -3.60, -2.40))
  }
  # The 80 percent check stores (-3.467957, -2.532043). Read from that
  # interval, as three of the functions once did, each of these would be
  # TRUE; at one-sided level 0.05 each is FALSE.
  out80 <- done_check(conf_level = 0.80)
  expect_equal(term(out80, "lower_limit"), -3.467957, tolerance = 1e-6)
  expect_false(decide_seq_budget_superiority(out80, margin = -3.5))
  expect_false(decide_seq_budget_inferiority(out80, margin = -2.5))
  expect_false(decide_seq_budget_equivalence(out80, -3.5, -2.5))
  # Stored limits that are plainly wrong change nothing.
  spoiled <- done_check()
  spoiled$value[spoiled$term == "lower_limit"] <- 100
  spoiled$value[spoiled$term == "upper_limit"] <- -100
  expect_true(decide_seq_budget_superiority(spoiled, margin = -3.61))
  expect_false(decide_seq_budget_superiority(spoiled, margin = -3.60))
  expect_true(decide_seq_budget_inferiority(spoiled, margin = -2.39))
  expect_true(decide_seq_budget_equivalence(spoiled, -3.61, -2.39))
  expect_true(decide_seq_budget_noninferiority(spoiled, margin = 3.61))
  # The limits the tests use are those a check stores when it is computed
  # with conf_level = 1 - 2 * alpha_level, as the help says.
  out90 <- done_check(conf_level = 1 - 2 * 0.05)
  expect_equal(term(out90, "lower_limit"), -3.600616, tolerance = 1e-6)
  expect_equal(term(out90, "upper_limit"), -2.399384, tolerance = 1e-6)
  # NULL is not a level: a level is always stated.
  for (f in decisions[1:3]) {
    expect_error(f(done_check(), 1, alpha_level = NULL),
                 "'alpha_level' must be a single finite number")
  }
  expect_error(decide_seq_budget_equivalence(done_check(), -1, 1,
                                             alpha_level = NULL),
               "'alpha_level' must be a single finite number")
})

test_that("the four decisions agree with direct implementations and with DMAR", {
  # Direct implementations of the rejection rules from the estimate and its
  # standard error: Equation C3 (noninferiority), C4 and C13 (superiority,
  # the margin being |psi_0| or psi_1), C14 (inferiority, the margin being
  # psi_2), and C15 (equivalence). And DMAR::equivalence_c(), the fixed
  # sample size cousin,
  # on its direct interface with degrees of freedom so large that its t
  # quantile is the normal one; it states the same hypotheses about the
  # bounds (-delta_lower, delta_upper) and reports the four decisions as
  # flags.
  withr::local_seed(113)
  times_true <- c(noninferior = 0L, superior = 0L, inferior = 0L,
                  equivalent = 0L)
  cases <- 60L
  for (i in seq_len(cases)) {
    psi_hat <- stats::runif(1, -2, 2)
    out <- done_check(means = c(psi_hat + 7, 4, 1))
    lower <- stats::runif(1, 0.1, 3)
    upper <- stats::runif(1, 0.1, 3)
    alpha <- stats::runif(1, 0.01, 0.2)
    c_hat <- term(out, "contrast")
    half <- stats::qnorm(alpha, lower.tail = FALSE) * sqrt(term(out, "variance"))
    direct <- c(noninferior = c_hat > -lower + half,
                superior = c_hat > upper + half,
                inferior = c_hat < -lower - half,
                equivalent = -lower + half < c_hat && c_hat < upper - half)
    dmar <- DMAR::equivalence_c(psi_hat = c_hat,
                                se = sqrt(term(out, "variance")),
                                df_error = 1e9, delta_lower = lower,
                                delta_upper = upper, alpha_level = alpha)
    got <- c(
      noninferior = decide_seq_budget_noninferiority(out, margin = lower,
                                                     alpha_level = alpha),
      superior = decide_seq_budget_superiority(out, margin = upper,
                                               alpha_level = alpha),
      inferior = decide_seq_budget_inferiority(out, margin = -lower,
                                               alpha_level = alpha),
      equivalent = decide_seq_budget_equivalence(out, margin_lower = -lower,
                                                 margin_upper = upper,
                                                 alpha_level = alpha)
    )
    expect_identical(got, direct)
    expect_identical(unname(got),
                     unname(dmar$value[match(names(got), dmar$term)] == 1))
    times_true <- times_true + got
  }
  # Each decision came out both ways, so no comparison above was empty.
  expect_true(all(times_true > 0L & times_true < cases))
})

test_that("noninferiority is the case margin = -abs(delta) of superiority", {
  # The identity the help states, over margins on both sides of each
  # threshold (3.600616 at 0.05, 3.715678 at 0.025, 3.467957 at 0.10) and
  # for checks stored at two confidence levels.
  outcomes <- logical(0)
  for (cl in c(0.95, 0.80)) {
    out <- done_check(conf_level = cl)
    for (alpha in c(0.05, 0.025, 0.10)) {
      for (delta in c(3.3, 3.59, 3.61, 3.9)) {
        ni <- decide_seq_budget_noninferiority(out, margin = delta,
                                               alpha_level = alpha)
        expect_identical(
          decide_seq_budget_superiority(out, margin = -abs(delta),
                                        alpha_level = alpha),
          ni)
        outcomes <- c(outcomes, ni)
      }
    }
  }
  # By hand: at 0.05 the margins 3.61 and 3.9 reject; at 0.025 only 3.9; at
  # 0.10 the margins 3.59, 3.61, and 3.9. Six of twelve, at either level of
  # confidence of the check.
  expect_identical(sum(outcomes), 12L)
  expect_identical(outcomes[1:12],
                   c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE,
                     FALSE, TRUE, TRUE, TRUE))
  # At the default level too.
  expect_identical(decide_seq_budget_noninferiority(done_check(), 3.61),
                   decide_seq_budget_superiority(done_check(), -abs(3.61)))
})

test_that("inferiority is the mirror image of superiority", {
  # Reversing the sign of every weight reverses the contrast and leaves its
  # variance alone, so 'below m' for one is 'above -m' for the other.
  out <- done_check()
  mirrored <- done_check(weights = -hand_w)
  expect_equal(term(mirrored, "contrast"), 3)
  expect_equal(term(mirrored, "variance"), 2 / 15)
  outcomes <- logical(0)
  for (alpha in c(0.05, 0.025)) {
    for (m in c(-2.5, -2.39, -2.35, -2.28, -2)) {
      res <- decide_seq_budget_inferiority(out, margin = m, alpha_level = alpha)
      expect_identical(
        decide_seq_budget_superiority(mirrored, margin = -m,
                                      alpha_level = alpha),
        res)
      outcomes <- c(outcomes, res)
    }
  }
  # By hand: the upper limit is -2.399384 at 0.05 (margins -2.39, -2.35,
  # -2.28, and -2 reject) and -2.284322 at 0.025 (-2.28 and -2 reject).
  expect_identical(outcomes, c(FALSE, TRUE, TRUE, TRUE, TRUE,
                               FALSE, FALSE, FALSE, TRUE, TRUE))
})

test_that("equivalence is the two one-sided tests, each at alpha_level", {
  out <- done_check()
  outcomes <- logical(0)
  for (alpha in c(0.05, 0.025)) {
    for (lo in c(-3.75, -3.65, -3.55)) {
      for (hi in c(-2.45, -2.35, -2.25)) {
        res <- decide_seq_budget_equivalence(out, margin_lower = lo,
                                             margin_upper = hi,
                                             alpha_level = alpha)
        expect_identical(
          decide_seq_budget_superiority(out, margin = lo,
                                        alpha_level = alpha) &&
            decide_seq_budget_inferiority(out, margin = hi,
                                          alpha_level = alpha),
          res)
        outcomes <- c(outcomes, res)
      }
    }
  }
  # By hand. At 0.05, limits (-3.600616, -2.399384): lower margins -3.75
  # and -3.65 with upper margins -2.35 and -2.25, four of nine. At 0.025,
  # limits (-3.715678, -2.284322): only (-3.75, -2.25).
  expect_identical(sum(outcomes[1:9]), 4L)
  expect_identical(which(outcomes[10:18]), 3L)
})

test_that("decide_* accept only a seq_budget() check", {
  withr::local_seed(113)
  smd_check <- seq_aipe_smd(stats::rnorm(40, 0.5), stats::rnorm(40),
                            width = 0.5)
  cv_check <- seq_mr_cv(stats::rgamma(40, 4), A = 400, c = 1)
  msg <- "must be an object returned by seq_budget\\(\\)"
  expect_error(decide_seq_budget_noninferiority(smd_check, 1), msg)
  expect_error(decide_seq_budget_superiority(smd_check, 0), msg)
  expect_error(decide_seq_budget_inferiority(cv_check, 0), msg)
  expect_error(decide_seq_budget_equivalence(smd_check, -1, 1), msg)
  expect_error(decide_seq_budget_superiority(data.frame(), 0), msg)
  # A partial object is no longer a check.
  expect_error(decide_seq_budget_superiority(done_check()[1:2, ], 0), msg)
  # The method label is what identifies the procedure: the same rows and
  # attributes under another label are refused.
  relabeled <- done_check()
  attr(relabeled, "method") <- "contrast"
  expect_error(decide_seq_budget_noninferiority(relabeled, 4),
               "a check from another procedure \\(\"contrast\"\\)")
  expect_error(seq_budget_next(relabeled), msg)
})

test_that("decide_* validate their margins and the level", {
  out <- done_check()
  # Noninferiority takes the size of the margin, a positive number.
  expect_error(decide_seq_budget_noninferiority(out),
               "'margin' must be supplied: the noninferiority margin")
  expect_error(decide_seq_budget_noninferiority(out, margin = -1),
               "'margin' must be greater than 0.*at or below -margin")
  expect_error(decide_seq_budget_noninferiority(out, margin = 0),
               "'margin' must be greater than 0")
  expect_error(decide_seq_budget_noninferiority(out, margin = NA), "'margin'")
  expect_error(decide_seq_budget_noninferiority(out, margin = c(1, 2)),
               "'margin'")
  expect_error(decide_seq_budget_noninferiority(out, margin = "1"), "'margin'")
  # Superiority and inferiority take the value itself, of either sign.
  expect_error(decide_seq_budget_superiority(out),
               "'margin' must be supplied: the value the contrast must exceed")
  expect_error(decide_seq_budget_inferiority(out),
               "'margin' must be supplied: the value the contrast must fall below")
  expect_error(decide_seq_budget_superiority(out, margin = NA),
               "'margin' must be a single finite number")
  expect_error(decide_seq_budget_superiority(out, margin = c(0, 1)),
               "'margin'")
  expect_error(decide_seq_budget_inferiority(out, margin = "0"), "'margin'")
  expect_error(decide_seq_budget_inferiority(out, margin = Inf), "'margin'")
  expect_silent(decide_seq_budget_superiority(out, margin = -1))
  expect_silent(decide_seq_budget_inferiority(out, margin = 1))
  # Equivalence takes its two margins with their signs, in order.
  expect_error(decide_seq_budget_equivalence(out, margin_upper = 1),
               "'margin_lower' must be supplied")
  expect_error(decide_seq_budget_equivalence(out, margin_lower = -1),
               "'margin_upper' must be supplied")
  expect_error(decide_seq_budget_equivalence(out, NA, 1), "'margin_lower'")
  expect_error(decide_seq_budget_equivalence(out, -1, NA), "'margin_upper'")
  expect_error(decide_seq_budget_equivalence(out, 1, 0),
               "'margin_lower' must be less than 'margin_upper'; got margin_lower = 1 and margin_upper = 0")
  expect_error(decide_seq_budget_equivalence(out, 1, 1),
               "'margin_lower' must be less than 'margin_upper'")
  # The level, for all four.
  for (f in decisions[1:3]) {
    expect_error(f(out, 1, alpha_level = 0), "'alpha_level'")
    expect_error(f(out, 1, alpha_level = 0.5), "'alpha_level'")
    expect_error(f(out, 1, alpha_level = 0.7), "'alpha_level'")
    expect_error(f(out, 1, alpha_level = NA), "'alpha_level'")
    expect_error(f(out, 1, alpha_level = c(0.05, 0.1)), "'alpha_level'")
  }
  expect_error(decide_seq_budget_equivalence(out, -1, 1, alpha_level = 0.5),
               "'alpha_level'")
  # The earlier argument names are gone.
  expect_error(decide_seq_budget_superiority(out, c_0 = 0), "unused argument")
  expect_error(decide_seq_budget_equivalence(out, c_lower = -1, c_upper = 1),
               "unused argument")
})

test_that("a decision before sampling has ended comes with a warning", {
  # Hand design at budget 1000 with n = (30, 20, 12): not ended.
  # Variance 5/6; 1.644854 * sqrt(5/6) = 1.501539, so the limits at
  # one-sided 0.05 are -3 -/+ 1.501539 = (-4.501539, -1.498461).
  out <- seq_budget(exact_groups(c(30, 20, 12), c(4, 4, 1), hand_s), hand_w,
                    hand_a, budget = 1000)
  expect_false(attr(out, "is_satisfied"))
  expect_warning(res <- decide_seq_budget_noninferiority(out, margin = 4.6),
                 "decision is provisional")
  expect_true(res)
  expect_warning(res <- decide_seq_budget_noninferiority(out, margin = 4.4),
                 "Sampling has not ended")
  expect_false(res)
  expect_warning(res <- decide_seq_budget_superiority(out, -4.6), "provisional")
  expect_true(res)
  expect_warning(res <- decide_seq_budget_inferiority(out, -1.4), "provisional")
  expect_true(res)
  expect_warning(res <- decide_seq_budget_equivalence(out, -4.6, -1.4),
                 "provisional")
  expect_true(res)
  expect_warning(res <- decide_seq_budget_equivalence(out, -4.4, -1.4),
                 "provisional")
  expect_false(res)
  # A completed check decides silently.
  for (f in decisions[1:3]) expect_silent(f(done_check(), margin = 4))
  expect_silent(decide_seq_budget_equivalence(done_check(), -4, 4))
})


# --- seq_budget_sensitivity ----------------------------------------------

test_that("seq_budget_sensitivity returns the documented rows", {
  out <- seq_budget_sensitivity(
    true_means = c(50, 50, 50), true_sds = c(9, 6, 6),
    contrast_weights = c(1, -0.5, -0.5), cost_per_group = c(38, 15, 10),
    budget = 1000, margin = 5, G = 20, seed = 113, m_0 = 5
  )
  # The summary carries the class every *_sensitivity() function returns,
  # and is otherwise a term/value data frame.
  expect_s3_class(out, c("smsd_sensitivity", "data.frame"), exact = TRUE)
  expect_identical(names(out), c("term", "value"))
  expect_type(out$term, "character")
  expect_type(out$value, "double")
  expect_identical(
    out$term,
    c("true_contrast", "budget", "conf_level", "replications", "coverage",
      "mean_total_cost", "sd_total_cost", "pct_within_budget",
      "pct_rule_met", "mean_n_1", "mean_n_2", "mean_n_3", "optimal_n_1",
      "optimal_n_2", "optimal_n_3", "mean_estimated_variance",
      "mean_true_variance", "empirical_variance", "optimal_variance",
      "margin", "alpha_level", "noninferiority_rejection_rate"))
  expect_true(all(is.finite(out$value)))
  expect_equal(term(out, "true_contrast"), 0)
  expect_equal(term(out, "replications"), 20)
  expect_equal(term(out, "margin"), 5)
  expect_equal(term(out, "alpha_level"), 0.05)
  # The optimal rows are the closed form for the true standard deviations.
  # Hand: |c| sigma sqrt(a) = 9*6.164414, 3*3.872983, 3*3.162278
  #       = 55.479726, 11.618950, 9.486833; sum 76.585509.
  #       V_o = 76.585509^2 / 1000 = 5.865340.
  expect_equal(term(out, "optimal_variance"), 5.865340, tolerance = 1e-6)
  expect_equal(sum(c(38, 15, 10) * out$value[13:15]), 1000)
  # Without a margin the test rows are NA and everything else is unchanged.
  out2 <- seq_budget_sensitivity(c(50, 50, 50), c(9, 6, 6), c(1, -0.5, -0.5),
                                 c(38, 15, 10), 1000, G = 20, seed = 113)
  expect_identical(out2$term, out$term)
  expect_true(all(is.na(out2$value[20:22])))
  expect_equal(out2$value[1:19], out$value[1:19])
  # It prints one line per row in fixed notation, with NA for the rows that
  # do not apply: the budget of 1000 and the 20 replications as written, not
  # as 1e+03 and 2e+01.
  printed <- utils::capture.output(print(out2))
  expect_length(printed, 23)
  expect_match(printed[1], "^term +value$")
  expect_match(printed[3], "^budget +1000$")
  expect_match(printed[5], "^replications +20$")
  expect_match(printed[23], "^noninferiority_rejection_rate +NA$")
  expect_false(any(grepl("e[+-][0-9]", printed)))
})

test_that("the argument order is that of the package's sensitivity functions", {
  # The design (populations, contrast, costs, budget, confidence level,
  # distribution, and the test reported), then G, seed, and m_0, then the
  # option that belongs to this procedure alone.
  expect_identical(
    names(formals(seq_budget_sensitivity)),
    c("true_means", "true_sds", "contrast_weights", "cost_per_group",
      "budget", "conf_level", "distribution", "margin", "alpha_level",
      "G", "seed", "m_0", "procedure"))
  expect_identical(formals(seq_budget_sensitivity)$G, 1000)
  expect_null(formals(seq_budget_sensitivity)$seed)
})

test_that("seq_budget_sensitivity is reproducible and leaves the RNG alone", {
  run <- function() {
    seq_budget_sensitivity(c(0, 0), c(1, 2), c(1, -1), c(3, 5), 300, G = 10,
                           seed = 113)
  }
  expect_identical(run(), run())
  # The session's stream is seeded with 1 here, not 113, on purpose: it must
  # differ from the seed the function is given, so that a function that left
  # its own stream behind would be caught.
  withr::local_seed(1)
  before <- stats::runif(1)
  withr::local_seed(1)
  invisible(run())
  expect_identical(stats::runif(1), before)
})

test_that("the engine is the loop a user writes with the exported functions", {
  # Same random number stream, same draws: final sizes, cost, estimate, and
  # variance must be identical, for the package's reading and the article's.
  mu <- c(1, 2, 3)
  sigma <- c(4.5, 0.9, 1)
  w <- c(1, -2, 1)
  a <- c(38, 15, 10)
  user_loop <- function(first_passage) {
    groups <- lapply(1:3, function(i) .budget_draw(6L, mu[i], sigma[i], "skewed"))
    closed <- NULL
    repeat {
      out <- seq_budget(groups, w, a, budget = 1500, closed = closed)
      if (first_passage) closed <- attr(out, "met_per_group")
      g <- seq_budget_next(out, order = if (first_passage) "rounds" else "shortfall")
      if (is.na(g)) break
      groups[[g]] <- c(groups[[g]], .budget_draw(1L, mu[g], sigma[g], "skewed"))
    }
    c(lengths(groups), term(out, "current_total_cost"), term(out, "contrast"),
      term(out, "variance"), as.numeric(attr(out, "rule_met")))
  }
  # Three seeds (113 and the next two) so that three different sampling paths
  # are compared, not one.
  for (seed in c(113, 114, 115)) {
    for (first_passage in c(FALSE, TRUE)) {
      by_hand <- withr::with_seed(seed, user_loop(first_passage))
      engine <- withr::with_seed(seed, .budget_sim_one(
        mu, sigma, w, a, 1500, 6L, "skewed", first_passage,
        if (first_passage) "rounds" else "shortfall"))
      expect_identical(engine, by_hand)
    }
  }
})

test_that("the engine is that loop for groups that show no variation too", {
  # The engine keeps each group's 'does it vary' flag from step to step and
  # looks again only at a group that has been constant so far; seq_budget()
  # works it out from the data at every check. With coarse data the two must
  # still agree to the last digit. Two data generating schemes stand in for
  # the engine's own (which is continuous and never produces a tie):
  #   ratings: whole numbers around 3, so a pilot of 2 is often constant;
  #   near one million: steps of 4e-5, so a group counts as varying only
  #     once its range reaches three steps (1.2e-10 of a million), and a
  #     group can cross the package's tolerance as observations arrive.
  w <- c(1, -0.5, -0.5)
  a <- c(7, 3, 5)
  schemes <- list(
    ratings = function(k, mean, sd, distribution) {
      round(mean + sd * stats::rnorm(k))
    },
    near_million = function(k, mean, sd, distribution) {
      1e6 + round(2 * sd * stats::rnorm(k)) * 4e-5
    }
  )
  unusable_steps <- c(ratings = 0L, near_million = 0L)
  became_usable <- c(ratings = 0L, near_million = 0L)
  for (scheme in names(schemes)) {
    local_mocked_bindings(.budget_draw = schemes[[scheme]])
    user_loop <- function(first_passage) {
      groups <- lapply(1:3, function(i) .budget_draw(2L, 3, 0.7, "normal"))
      closed <- NULL
      was_unusable <- FALSE
      repeat {
        out <- seq_budget(groups, w, a, budget = 150, closed = closed)
        unusable <- attr(out, "variance_unusable")
        unusable_steps[scheme] <<- unusable_steps[scheme] + unusable
        became_usable[scheme] <<- became_usable[scheme] +
          (was_unusable && !unusable)
        was_unusable <- unusable
        if (first_passage) closed <- attr(out, "met_per_group")
        g <- seq_budget_next(out, order = if (first_passage) "rounds" else "shortfall")
        if (is.na(g)) break
        groups[[g]] <- c(groups[[g]], .budget_draw(1L, 3, 0.7, "normal"))
      }
      c(lengths(groups), term(out, "current_total_cost"),
        term(out, "contrast"), term(out, "variance"),
        as.numeric(attr(out, "rule_met")))
    }
    for (seed in 113:118) {
      for (first_passage in c(FALSE, TRUE)) {
        by_hand <- withr::with_seed(seed, user_loop(first_passage))
        engine <- withr::with_seed(seed, .budget_sim_one(
          rep(3, 3), rep(0.7, 3), w, a, 150, 2L, "normal", first_passage,
          if (first_passage) "rounds" else "shortfall"))
        expect_identical(engine, by_hand,
                         label = paste(scheme, seed, first_passage))
      }
    }
  }
  # Both schemes visit the state under test, and in both a group goes from
  # unusable to usable along the way (110 and 160 such steps and 9 such
  # changes each under these seeds).
  expect_gt(unusable_steps[["ratings"]], 10L)
  expect_gt(unusable_steps[["near_million"]], 10L)
  expect_gt(became_usable[["ratings"]], 2L)
  expect_gt(became_usable[["near_million"]], 2L)
})

# A small valid call of seq_budget_sensitivity() in which one argument at a
# time is replaced.
sens_call <- function(...) {
  args <- utils::modifyList(
    list(true_means = c(0, 0), true_sds = c(1, 1),
         contrast_weights = c(1, -1), cost_per_group = c(10, 10),
         budget = 500, G = 5, seed = 113),
    list(...))
  do.call(seq_budget_sensitivity, args)
}

test_that("seq_budget_sensitivity validates its arguments", {
  f <- sens_call
  expect_error(f(budget = 50), "pilot alone \\(m_0 = 5 per group\\) costs 100")
  expect_error(f(G = 2.5), "'G' must be a single whole number")
  expect_error(f(G = 1), "'G' must be a single whole number of at least 2")
  expect_error(f(m_0 = 1), "'m_0'")
  expect_error(f(m_0 = 4.5), "'m_0'")
  # The seed goes through the package's shared check, before any study is
  # run: a number too large for set.seed() is reported in the package's own
  # words, not base R's.
  seed_msg <- "'seed' must be NULL or a single whole number that set.seed\\(\\) accepts"
  expect_error(f(seed = 1.5), seed_msg)
  expect_error(f(seed = 1e10), paste0(seed_msg, "; got 1e\\+10"))
  expect_error(f(seed = "a"), seed_msg)
  expect_error(f(seed = c(1, 2)), seed_msg)
  # Every argument is checked before the design is judged: a bad seed is
  # reported even when the budget could not pay for the pilot.
  expect_error(f(seed = 1.5, budget = 50), seed_msg)
  expect_error(f(true_sds = c(1, 0)),
               "'true_sds' must contain only positive numbers \\(one standard deviation per population\\); entry 2 is 0")
  expect_error(f(true_means = c(0, NA)),
               "'true_means' must contain only finite numbers \\(one mean per population\\); entry 2 is NA")
  expect_error(f(true_means = c(0, 0, 0)),
               "'true_means' must have 2 entries \\(one per entry of 'contrast_weights'\\); it has 3")
  expect_error(f(contrast_weights = c(0, 0)), "at least one nonzero")
  expect_error(f(cost_per_group = c(10, -1)),
               "'cost_per_group' must contain only positive numbers.*entry 2 is -1")
  expect_error(f(conf_level = 1.5), "'conf_level'")
  expect_error(f(margin = -1), "'margin'")
  expect_error(f(margin = 1, alpha_level = 0.9), "'alpha_level'")
  # alpha_level is checked whether or not a margin is given, so a bad value
  # is not accepted silently when margin is NULL.
  expect_error(f(alpha_level = 5), "'alpha_level' must be a single finite number")
  expect_error(f(alpha_level = NA), "'alpha_level'")
  expect_error(f(distribution = "cauchy"),
               "'distribution' must be one of \"normal\", \"skewed\"; got \"cauchy\"")
  expect_error(f(procedure = "greedy"),
               "'procedure' must be one of \"reallocate\", \"first_passage\"; got \"greedy\"")
  expect_error(f(distribution = c("skewed", "normal")),
               "'distribution' must be one of")
  # An unambiguous abbreviation selects the choice, as for every such
  # argument in the package; the result is that of the full name.
  expect_identical(f(distribution = "skew"), f(distribution = "skewed"))
  expect_identical(f(procedure = "first"), f(procedure = "first_passage"))
  expect_false(identical(f(distribution = "skew"), f()))
  expect_error(f(true_means = c(a = 0, b = 0), true_sds = c(b = 1, a = 1)),
               "'true_sds' is named \\(b, a\\) but 'true_means' is named \\(a, b\\)")
})

test_that("a design that cannot be run is refused before any study is simulated", {
  # The engine is replaced by one that fails if it is reached, so each of
  # these is seen to be refused by the argument checks, in the words of this
  # function's own arguments, and not after G studies have been run.
  local_mocked_bindings(
    .budget_sim_one = function(...) stop("a study was simulated")
  )
  f <- sens_call
  expect_error(f(), "a study was simulated")
  expect_error(f(true_sds = c(1, 1, 1)),
               "'true_sds' must have 2 entries \\(one per entry of 'contrast_weights'\\); it has 3")
  expect_error(f(cost_per_group = c(10, 10, 10)),
               "'cost_per_group' must have 2 entries \\(one per entry of 'contrast_weights'\\); it has 3")
  # A budget that pays for more observations than can be counted: the
  # package's message for a sample size that cannot be represented, with the
  # advice that names the two arguments that set it (1e300 / 10 = 1e299).
  expect_error(f(budget = 1e300),
               "number of observations that the budget pays for in the cheapest group \\(1e\\+299\\) is too large to represent. Check 'budget' and 'cost_per_group'")
  # A budget that pays for more observations than a simulated study may
  # have. A study is simulated one observation at a time in a time that
  # grows with the square of its size, so a budget of 3e9 at a cost of 10
  # (300 million observations) would be a run that does not finish; it is
  # refused before anything is simulated. The limit is one million whole
  # observations in the cheapest group: 1e7 + 9 pays for 1,000,000 and
  # reaches the engine, 1e7 + 10 pays for 1,000,001.
  expect_error(
    f(budget = 3e9),
    paste0("^The budget pays for 300000000 observations in the cheapest ",
           "group, and seq_budget_sensitivity\\(\\) simulates studies of at ",
           "most 1000000: a study is simulated one observation at a time, ",
           "and the time for one replication grows with the square of its ",
           "size\\. Simulate a smaller budget to see how the procedure ",
           "behaves, and use seq_budget_optimal\\(\\) for the allocation at ",
           "the full budget\\.$")
  )
  expect_error(f(budget = 1e7 + 10),
               "^The budget pays for 1000001 observations in the cheapest group")
  expect_error(f(budget = 1e7 + 9), "a study was simulated")
  expect_error(f(budget = 1e7), "a study was simulated")
  # The cheapest group sets the count, whichever position it is in.
  expect_error(f(cost_per_group = c(40, 0.004), budget = 5000),
               "^The budget pays for 1250000 observations")
  expect_error(f(cost_per_group = c(40, 0.005), budget = 5000),
               "a study was simulated")
  expect_identical(.budget_sens_max_n, 1e6)
  # Standard deviations so large that the optimal variance overflows
  # ((2 * 1e200 * sqrt(10))^2 / 500 is 8e398): the message names true_sds,
  # not the 'sigmas' of seq_budget_optimal().
  expect_error(f(true_sds = c(1e200, 1e200)),
               "optimal allocation cannot be computed: the values in 'contrast_weights', 'true_sds', 'cost_per_group', and 'budget'.*Rescale")
})

test_that("a simulated study that overflows is reported, not summarized", {
  # Weights of 0.001 and standard deviations of 1e155: the optimal variance,
  # (2 * 0.001 * 1e155 * sqrt(10))^2 / 500 = 8e302, can be represented, so
  # the design passes the checks. But the scores are of the order of 1e155,
  # so a sample variance is of the order of 1e310, beyond the largest
  # double (about 1.8e308) however the sum of squares is accumulated: on
  # every platform, with or without extended precision, it is Inf.
  # seq_budget() stops on such data, and so does the simulation, rather than
  # return a summary built on estimated variances of Inf.
  expect_error(
    seq_budget_sensitivity(c(0, 0), c(1e155, 1e155), c(0.001, -0.001),
                           c(10, 10), 500, G = 5, seed = 113),
    "contrast estimate or its variance in a simulated study cannot be computed: the values in 'true_means', 'true_sds', and 'contrast_weights'.*Rescale")
  # A hundred times smaller and everything is finite: V_o = 8e296.
  out <- seq_budget_sensitivity(c(0, 0), c(1e152, 1e152), c(0.001, -0.001),
                                c(10, 10), 500, G = 5, seed = 113)
  expect_true(all(is.finite(out$value[1:17])))
  expect_equal(out$value[out$term == "optimal_variance"], 8e296)
})

test_that("means so far from zero that the location would matter are refused", {
  # Scores are mean + sd * Z, and the package treats scores that agree to ten
  # significant digits as showing no variation. Beside a mean of M standard
  # deviations that is agreement within M * 1e-10 standard deviations. A mean
  # of 1e6 standard deviations or more (a width of 1e-4 of a standard
  # deviation or more) is refused in a group with a nonzero weight: beyond
  # it, pilots whose scores happen to fall close together begin to be
  # treated as constant, and the summary depends on where the scores are
  # located. (In 10,000 studies of this design with a pilot of 2, each
  # compared with the same study at means of zero, 1 took another path at
  # means of 1e6, 158 at 1e8, and 1,618 at 1e9, where the coverage had moved
  # from 0.901 to 0.942; from 1e11 every study was sampled by rounds, 43, 43,
  # 42 against an optimal 55.1, 33.7, 33.7.)
  f <- function(means, sds = c(1, 1, 1), w = c(1, -0.5, -0.5), m_0 = 5,
                G = 40) {
    seq_budget_sensitivity(means, sds, w, c(3, 2, 2), 300, G = G,
                           seed = 113, m_0 = m_0)
  }
  expect_error(f(c(1e12, 0, 0)),
               "'true_means' must each be less than 1e\\+06 standard deviations \\('true_sds'\\) from zero in a group with a nonzero weight; entry 1 is 1e\\+12 with a standard deviation of 1\\.")
  # The limit itself, on either side of zero, and just inside it.
  expect_error(f(c(1e6, 0, 0)), "entry 1 is 1e\\+06 with")
  expect_error(f(c(0, 0, -1e6)), "entry 3 is -1e\\+06 with")
  expect_error(f(c(1e12, 0, 1e13)), "entry 1 is 1e\\+12.*2 entries in all")
  expect_s3_class(f(c(999999, 0, 0), G = 5), "smsd_sensitivity")
  # It is the mean beside its own standard deviation that matters:
  # 50 / 1e-5 = 5e6 standard deviations is refused, 1e8 / 1000 = 1e5 is not.
  expect_error(f(c(50, 50, 50), sds = c(1e-5, 1, 1)),
               "entry 1 is 50 with a standard deviation of 1e-05")
  expect_s3_class(f(c(1e8, 0, 0), sds = c(1000, 1, 1), G = 5),
                  "smsd_sensitivity")
  # A group with a weight of zero takes no part in the contrast and receives
  # no observation after the pilot, so its mean may be anything: the summary
  # is the one obtained with that mean at zero, row for row.
  expect_identical(f(c(5, 6, 1e12), w = c(1, -1, 0), G = 20),
                   f(c(5, 6, 0), w = c(1, -1, 0), G = 20))
  # Inside the limit the location does not matter, even with the smallest
  # pilot, where it would matter first: moving every mean by 900,000
  # standard deviations leaves a contrast whose weights sum to zero, and
  # every row of the summary, where it was (up to the rounding error of
  # scores carried with six fewer digits).
  at_zero <- f(c(0, 0, 0), m_0 = 2)
  far <- f(rep(9e5, 3), m_0 = 2)
  expect_equal(far$value[1:19], at_zero$value[1:19], tolerance = 1e-6)
  # The optimal sizes there are the closed form: with equal standard
  # deviations, 300 * (1, 0.5, 0.5) / (sqrt(c(3, 2, 2)) * (sqrt(3) +
  # sqrt(2))), that is, 55.05, 33.71, 33.71.
  expect_equal(far$value[far$term %in% paste0("optimal_n_", 1:3)],
               c(55.05, 33.71, 33.71), tolerance = 1e-3)
})

test_that("the two data generating shapes are the documented ones", {
  # mean + sd * Z with Z standard normal, or Z a chi square variable with 3
  # degrees of freedom standardized to mean 0 and variance 1.
  expect_identical(withr::with_seed(113, .budget_draw(6L, 2, 3, "normal")),
                   withr::with_seed(113, 2 + 3 * stats::rnorm(6)))
  expect_identical(
    withr::with_seed(113, .budget_draw(6L, 2, 3, "skewed")),
    withr::with_seed(113, 2 + 3 * ((stats::rchisq(6, df = 3) - 3) / sqrt(6))))
})

test_that("every summary row is what a user computes from the same studies", {
  # The studies are rebuilt with the exported functions only, under the same
  # seed, and every row of the summary is recomputed from them. Nothing here
  # is at its default: conf_level 0.80, alpha_level 0.10, a pilot of 7, both
  # distributions, and both procedures. The margin of 2.5 puts the rejection
  # rate in the middle of its range (the true contrast is -1), and under the
  # first passage rule some studies end by the rule and some by the budget.
  mu <- c(0, 2, 3)
  sigma <- c(4.5, 0.9, 1)
  w <- c(1, -2, 1)
  a <- c(38, 15, 10)
  A0 <- 1500
  G <- 12
  m <- 7L
  cl <- 0.80
  alpha <- 0.10
  margin <- 2.5
  psi <- sum(w * mu)
  draw <- function(k, i, distribution) {
    z <- if (distribution == "normal") {
      stats::rnorm(k)
    } else {
      (stats::rchisq(k, df = 3) - 3) / sqrt(6)
    }
    mu[i] + sigma[i] * z
  }
  results <- list()
  for (procedure in c("reallocate", "first_passage")) {
    for (distribution in c("normal", "skewed")) {
      fp <- procedure == "first_passage"
      by_hand <- withr::with_seed(113, replicate(G, {
        groups <- lapply(1:3, function(i) draw(m, i, distribution))
        closed <- NULL
        repeat {
          out <- seq_budget(groups, w, a, A0, conf_level = cl, closed = closed)
          if (fp) closed <- attr(out, "met_per_group")
          g <- seq_budget_next(out, order = if (fp) "rounds" else "shortfall")
          if (is.na(g)) break
          groups[[g]] <- c(groups[[g]], draw(1L, g, distribution))
        }
        c(lengths(groups), term(out, "current_total_cost"),
          term(out, "lower_limit") <= psi && psi <= term(out, "upper_limit"),
          decide_seq_budget_noninferiority(out, margin, alpha),
          attr(out, "rule_met"), term(out, "contrast"), term(out, "variance"))
      }))
      got <- seq_budget_sensitivity(mu, sigma, w, a, A0, conf_level = cl,
                                    distribution = distribution,
                                    margin = margin, alpha_level = alpha,
                                    G = G, seed = 113, m_0 = m,
                                    procedure = procedure)
      v <- stats::setNames(got$value, got$term)
      label <- paste(procedure, distribution)
      expect_equal(unname(v[c("mean_n_1", "mean_n_2", "mean_n_3")]),
                   rowMeans(by_hand[1:3, ]), label = label)
      expect_equal(v[["mean_total_cost"]], mean(by_hand[4, ]), label = label)
      expect_equal(v[["sd_total_cost"]], stats::sd(by_hand[4, ]),
                   label = label)
      expect_equal(v[["pct_within_budget"]], 100 * mean(by_hand[4, ] <= A0),
                   label = label)
      expect_equal(v[["coverage"]], mean(by_hand[5, ]), label = label)
      expect_equal(v[["noninferiority_rejection_rate"]], mean(by_hand[6, ]),
                   label = label)
      expect_equal(v[["pct_rule_met"]], 100 * mean(by_hand[7, ]),
                   label = label)
      expect_equal(v[["empirical_variance"]], mean((by_hand[8, ] - psi)^2),
                   label = label)
      expect_equal(v[["mean_estimated_variance"]], mean(by_hand[9, ]),
                   label = label)
      expect_equal(v[["mean_true_variance"]],
                   mean(colSums(w^2 * sigma^2 / by_hand[1:3, ])),
                   label = label)
      expect_equal(v[["true_contrast"]], psi)
      expect_equal(v[["conf_level"]], cl)
      expect_equal(v[["alpha_level"]], alpha)
      expect_equal(v[["margin"]], margin)
      results[[label]] <- by_hand
    }
  }
  # The comparison would be empty if the four configurations coincided, if
  # no test rejected or all did, or if no first passage study ended by the
  # rule or all did.
  expect_false(isTRUE(all.equal(results[["reallocate normal"]],
                                results[["first_passage normal"]])))
  expect_false(isTRUE(all.equal(results[["reallocate normal"]],
                                results[["reallocate skewed"]])))
  rejected <- mean(results[["reallocate normal"]][6, ])
  expect_gt(rejected, 0)
  expect_lt(rejected, 1)
  ended_by_rule <- mean(results[["first_passage skewed"]][7, ])
  expect_gt(ended_by_rule, 0)
  expect_lt(ended_by_rule, 1)
  expect_identical(mean(results[["reallocate normal"]][7, ]), 0)
})

test_that("the summary of a set of studies matches a hand computation", {
  # Four studies of two groups, written out by hand, with weights (1, -1),
  # true means (0, 0), true SDs (1, 1), unit costs (10, 10), budget 100:
  #   study:      1      2      3      4
  #   n_1:        5      4      6      5
  #   n_2:        5      5      5      4
  #   cost:     100     90    110     90
  #   estimate: 0.5   -1.2      0    2.1
  #   variance: 0.25   0.16     1      1
  #   rule met:   1      0      0      1
  # The third costs more than the budget. The engine cannot produce such a
  # study; it is here to show that pct_within_budget is computed from the
  # costs and would report an engine that overspent.
  runs <- rbind(c(5, 4, 6, 5), c(5, 5, 5, 4), c(100, 90, 110, 90),
                c(0.5, -1.2, 0, 2.1), c(0.25, 0.16, 1, 1), c(1, 0, 0, 1))
  out <- .budget_summarize(runs, means = c(0, 0), sds = c(1, 1), w = c(1, -1),
                           a = c(10, 10), budget = 100, conf_level = 0.95,
                           margin = 1, alpha_level = 0.05)
  v <- stats::setNames(out$value, out$term)
  expect_equal(v[["true_contrast"]], 0)
  expect_equal(v[["replications"]], 4)
  # 95% half widths 1.96 * sqrt(variance) = 0.98, 0.78, 1.96, 1.96: the
  # true contrast 0 is inside for studies 1 and 3.
  expect_equal(v[["coverage"]], 0.5)
  # Costs 100, 90, 110, 90: mean 97.5; squared deviations 6.25, 56.25,
  # 156.25, 56.25 sum to 275, so the SD is sqrt(275 / 3).
  expect_equal(v[["mean_total_cost"]], 97.5)
  expect_equal(v[["sd_total_cost"]], sqrt(275 / 3))
  expect_equal(v[["pct_within_budget"]], 75)
  expect_equal(v[["pct_rule_met"]], 50)
  expect_equal(unname(v[c("mean_n_1", "mean_n_2")]), c(5, 4.75))
  # Equal groups: n_i0 = 100 / (2 * 10) = 5; V_o = (2 sqrt(10))^2 / 100.
  expect_equal(unname(v[c("optimal_n_1", "optimal_n_2")]), c(5, 5))
  expect_equal(v[["optimal_variance"]], 0.4)
  expect_equal(v[["mean_estimated_variance"]], (0.25 + 0.16 + 1 + 1) / 4)
  # True variance 1 / n_1 + 1 / n_2: 0.4, 0.45, 0.3667, 0.45.
  expect_equal(v[["mean_true_variance"]],
               (2 / 5 + 9 / 20 + 11 / 30 + 9 / 20) / 4)
  expect_equal(v[["empirical_variance"]], (0.25 + 1.44 + 0 + 4.41) / 4)
  # Noninferiority, margin 1, alpha 0.05: lower limits estimate - 1.645 *
  # sqrt(variance) = -0.32, -1.86, -1.64, 0.46; above -1 for studies 1, 4.
  expect_equal(v[["noninferiority_rejection_rate"]], 0.5)
  expect_equal(v[["margin"]], 1)
  expect_equal(v[["alpha_level"]], 0.05)
  # At alpha 0.25 the critical value is 0.674: lower limits 0.16, -1.47,
  # -0.67, 1.43, so study 3 rejects as well.
  out <- .budget_summarize(runs, c(0, 0), c(1, 1), c(1, -1), c(10, 10), 100,
                           conf_level = 0.95, margin = 1, alpha_level = 0.25)
  expect_equal(out$value[out$term == "noninferiority_rejection_rate"], 0.75)
  # At 50 percent confidence the half widths are 0.674 * sqrt(variance) =
  # 0.34, 0.27, 0.67, 0.67: only study 3 covers.
  out <- .budget_summarize(runs, c(0, 0), c(1, 1), c(1, -1), c(10, 10), 100,
                           conf_level = 0.50, margin = NULL, alpha_level = 0.05)
  expect_equal(out$value[out$term == "coverage"], 0.25)
  expect_true(all(is.na(out$value[out$term %in%
                                    c("margin", "alpha_level",
                                      "noninferiority_rejection_rate")])))
})

test_that("operating characteristics with the article's unequal costs", {
  # The article's Table 1 normal row: means (1, 2, 3), SDs (4.5, 0.9, 1),
  # weights (1, -2, 1), costs (38, 15, 10), budget 5000, pilot 10.
  # G = 625 replications. Every tolerance leaves at least four Monte Carlo
  # standard errors on each side of the procedure's own long-run value
  # (6,000 replications of this design under another seed), so the test
  # would not fail under another random stream either:
  #   coverage: long run 0.943 against the nominal 0.95.
  #     SE = sqrt(0.95 * 0.05 / 625) = 0.0087. The tolerance of 0.044 around
  #     the nominal level is five SEs and leaves 4.2 SEs below the long-run
  #     value (0.943 - 0.906 = 0.037).
  #   mean n_i: the SD of N_i is about 2.7, 5.8, 5.1, so the SE is about
  #     0.11, 0.23, 0.21; tolerances 1.5, 2, 2.
  #   mean cost: long run 4989.0 with an SD of 8.8, so the SE is 0.35.
  out <- seq_budget_sensitivity(
    true_means = c(1, 2, 3), true_sds = c(4.5, 0.9, 1),
    contrast_weights = c(1, -2, 1), cost_per_group = c(38, 15, 10),
    budget = 5000, G = 625, seed = 113, m_0 = 10
  )
  v <- stats::setNames(out$value, out$term)
  # Never over budget, by construction.
  expect_identical(v[["pct_within_budget"]], 100)
  expect_lte(v[["mean_total_cost"]], 5000)
  # Judged on the current data the rule cannot be met short of the budget,
  # so every run ends by exhaustion, within one unit cost of the budget.
  expect_identical(v[["pct_rule_met"]], 0)
  expect_gt(v[["mean_total_cost"]], 5000 - 15)
  expect_lt(abs(v[["coverage"]] - 0.95), 0.044)
  # Optimal allocation printed in Table 1: 96.373, 61.3566, 41.7479.
  expect_equal(unname(v[c("optimal_n_1", "optimal_n_2", "optimal_n_3")]),
               c(96.373, 61.3566, 41.7479), tolerance = 1e-5)
  expect_lt(abs(v[["mean_n_1"]] - 96.373), 1.5)
  expect_lt(abs(v[["mean_n_2"]] - 61.3566), 2)
  expect_lt(abs(v[["mean_n_3"]] - 41.7479), 2)
  # The variance is close to the optimum V_o = 0.2869 (Table 1). No
  # allocation within the budget can beat V_o, so the true variance at the
  # final sizes is above it, by less than 2 percent here (long run: 0.6
  # percent above). The empirical variance is 2.8 percent above V_o in the
  # long run, and one squared deviation has an SD of 1.45 times V_o, so at
  # G = 625 its SE is 0.058 of V_o: the tolerance of 0.28 leaves 4.4 SEs
  # above the long-run value.
  expect_equal(v[["optimal_variance"]], 0.2869, tolerance = 1e-4)
  expect_gt(v[["mean_true_variance"]], v[["optimal_variance"]])
  expect_lt(v[["mean_true_variance"]] / 0.2869 - 1, 0.02)
  expect_lt(abs(v[["mean_estimated_variance"]] / 0.2869 - 1), 0.05)
  expect_lt(abs(v[["empirical_variance"]] / 0.2869 - 1), 0.28)
})

test_that("the noninferiority test holds its level at the boundary", {
  skip_on_cran()
  # True contrast exactly at -margin (the boundary of H0): means (0, 2, 3)
  # give 0 - 4 + 3 = -1, margin 1. The rejection rate estimates the Type I
  # error rate, nominally 0.05; the long-run rate of this design (budget
  # 3000, pilot 10) is 0.052 (20,000 replications). G = 800:
  # SE = sqrt(0.05 * 0.95 / 800) = 0.0077. The tolerance of 0.039 around
  # the nominal level is five SEs and leaves 4.7 SEs above the long-run
  # rate. (The direction and the critical value of the rejection rule are
  # tested without simulation error, and on CRAN too, in "every summary row
  # is what a user computes from the same studies".)
  out <- seq_budget_sensitivity(
    true_means = c(0, 2, 3), true_sds = c(4.5, 0.9, 1),
    contrast_weights = c(1, -2, 1), cost_per_group = c(38, 15, 10),
    budget = 3000, margin = 1, alpha_level = 0.05, G = 800, seed = 113,
    m_0 = 10
  )
  v <- stats::setNames(out$value, out$term)
  expect_equal(v[["true_contrast"]], -1)
  expect_lt(abs(v[["noninferiority_rejection_rate"]] - 0.05), 0.039)
  expect_identical(v[["pct_within_budget"]], 100)
})

test_that("the first passage reading reproduces the article's mean cost", {
  skip_on_cran()
  # Under the first passage rule every group can have met its criterion with
  # money left, so a good share of runs end by the rule (about 31 percent in
  # the long run) and the mean cost is lower than under the package's
  # default (4989.0, tested above). The mean final sample sizes printed in
  # the normal row of the article's Table 1 imply a mean cost of
  #   38 * 95.868 + 15 * 61.3886 + 10 * 41.3042 = 4976.9,
  # which the default reading does not reproduce and this one does: its
  # long-run mean cost is 4976.0 (6,000 replications) with an SD of 18.8.
  # G = 400, so the SE of the mean is 0.94; the tolerance of 6 around the
  # article's figure leaves 5 SEs beyond the long-run value. (That the
  # procedure argument reaches the engine at all is tested on CRAN too, in
  # "every summary row is what a user computes from the same studies".)
  out <- seq_budget_sensitivity(
    true_means = c(1, 2, 3), true_sds = c(4.5, 0.9, 1),
    contrast_weights = c(1, -2, 1), cost_per_group = c(38, 15, 10),
    budget = 5000, G = 400, seed = 113, m_0 = 10,
    procedure = "first_passage"
  )
  v <- stats::setNames(out$value, out$term)
  expect_identical(v[["pct_within_budget"]], 100)
  expect_gt(v[["pct_rule_met"]], 10)
  expect_lt(v[["pct_rule_met"]], 60)
  # The default reading's long-run mean cost, 4989.0, is 12 away and would
  # fail this.
  expect_lt(abs(v[["mean_total_cost"]] - 4976.9), 6)
})
