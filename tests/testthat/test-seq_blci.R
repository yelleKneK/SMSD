# Tests for the purely sequential bounded length confidence interval for a
# median, seq_blci_dfr_median(), and its pilot. The order statistic interval
# [X_(b), X_(n - b + 1)] covers the median of a continuous distribution
# exactly when a Binomial(n, 1/2) count lies between b and n - b, so the
# expected values here come from the binomial distribution, not from the
# package's normal approximation.

exact_coverage <- function(b, n) sum(stats::dbinom(b:(n - b), n, 0.5))

# ---- minimum sample size ----------------------------------------------------

test_that("the minimum n is the smallest with 1 - 2^(1 - n) >= conf_level", {
  # By hand: 1 - 2^(1 - n) is 0.5, 0.75, 0.875, 0.9375, 0.96875, 0.984375,
  # 0.9921875 for n = 2..8, and 0.99805, 0.99902 for n = 10, 11.
  expected <- c("0.5" = 2L, "0.8" = 4L, "0.9" = 5L, "0.95" = 6L,
                "0.99" = 8L, "0.999" = 11L)
  for (cf in names(expected)) {
    conf <- as.numeric(cf)
    n_min <- expected[[cf]]
    out <- seq_blci_dfr_median(seq_len(n_min), d = 1e6, conf_level = conf)
    expect_identical(attr(out, "min_n"), n_min)
    # At the minimum n the interval is [min, max] and, with a bound this
    # loose, the rule is met.
    expect_true(attr(out, "is_satisfied"))
    expect_equal(out$value[c(1, 3)], c(1, n_min))
    expect_identical(attr(out, "ranks"), c(lower = 1L, upper = n_min))
    # The pilot never returns less than the minimum n.
    expect_identical(seq_blci_dfr_median_pilot(d = 1e6, conf_level = conf,
                                               m_0 = 2), n_min)
  }
  # A level that the widest interval attains exactly: at n = 5 the coverage
  # of [X_(1), X_(5)] is exactly 0.9375.
  out <- seq_blci_dfr_median(1:5, d = 100, conf_level = 0.9375)
  expect_identical(attr(out, "min_n"), 5L)
  expect_true(attr(out, "is_satisfied"))
})

test_that("below the minimum n the rule is 'not met' with a note, not an error", {
  # At 95% at least 6 observations are needed. With 2 to 5 the widest
  # interval has coverage 0.50, 0.75, 0.875, 0.9375, so no interval is
  # reported however loose the length bound.
  for (n in 1:5) {
    out <- seq_blci_dfr_median(c(1, 2, 3, 4, 100)[seq_len(n)], d = 1000)
    expect_identical(attr(out, "is_satisfied"), FALSE)
    expect_true(all(is.na(out$value[out$term %in% c("lower_limit",
                                                    "upper_limit",
                                                    "interval_length")])))
    expect_equal(out$value[out$term == "median"],
                 stats::median(c(1, 2, 3, 4, 100)[seq_len(n)]))
    expect_match(attr(out, "notes"), "at least 6 are needed")
    expect_match(attr(out, "notes"), "95% confidence")
    # The count is in the singular for one observation.
    expect_match(attr(out, "notes"),
                 if (n == 1) "^With 1 observation no order" else
                   sprintf("^With %d observations no order", n))
    expect_output(print(out), "Stopping rule met: FALSE")
  }
  expect_true(attr(seq_blci_dfr_median(c(1, 2, 3, 4, 100, 5), d = 1000),
                   "is_satisfied"))
  # At 99% the minimum is 8.
  expect_false(attr(seq_blci_dfr_median(1:7, d = 1000, conf_level = 0.99),
                    "is_satisfied"))
  expect_true(attr(seq_blci_dfr_median(1:8, d = 1000, conf_level = 0.99),
                   "is_satisfied"))
})

# ---- the interval -----------------------------------------------------------

test_that("the interval is the right pair of order statistics (n = 20, 95%)", {
  # With n = 20 the pair of ranks (6, 15) has exact coverage 0.9586 and the
  # next pair in, (7, 14), has 0.8847, so (6, 15) is the narrowest admissible
  # interval. On a shuffle of 1..20 the limits are therefore 6 and 15.
  expect_gte(exact_coverage(6, 20), 0.95)
  expect_lt(exact_coverage(7, 20), 0.95)
  set.seed(113)
  x <- sample(1:20)
  out <- seq_blci_dfr_median(x, d = 4.5)
  expect_s3_class(out, c("smsd_check", "data.frame"), exact = TRUE)
  expect_identical(out$term, c("lower_limit", "median", "upper_limit",
                               "interval_length", "current_n"))
  expect_equal(out$value, c(6, 10.5, 15, 9, 20))
  expect_identical(attr(out, "ranks"), c(lower = 6L, upper = 15L))
  # The rule is "length at most 2 d": met at d = 4.5 (9 <= 9), not at 4.4.
  expect_true(attr(out, "is_satisfied"))
  expect_false(attr(seq_blci_dfr_median(x, d = 4.4), "is_satisfied"))
  expect_true(is.na(attr(out, "criterion_n")))
  expect_equal(attr(out, "width"), 9)               # the bound 2 d
  expect_equal(attr(out, "conf_level"), 0.95)
  expect_identical(attr(out, "method"),
                   "bounded length confidence interval (distribution free)")
})

test_that("a met rule carries the note about sequential coverage; an unmet one does not", {
  # At a fixed n the interval has at least the stated coverage, but a
  # sequential procedure that stops when the interval is short covers less
  # often (0.92 to 0.94 at a nominal 0.95 in simulation; see the operating
  # characteristic test below). The object says so when the rule is met, and
  # print() shows it.
  set.seed(113)
  x <- sample(1:20)
  met <- seq_blci_dfr_median(x, d = 4.5)
  expect_true(attr(met, "is_satisfied"))
  expect_length(attr(met, "notes"), 1L)
  expect_match(attr(met, "notes"), "less often than stated")
  # The ranges are those measured for normal, exponential, and Cauchy data
  # at four to six target widths (40,000 runs each; 0.842 to 0.879, 0.918
  # to 0.947, and 0.983 to 0.989, the figures in the help page).
  expect_match(attr(met, "notes"), "0.84 to 0.88 for a nominal 0.90")
  expect_match(attr(met, "notes"), "0.92 to 0.94 for 0.95")
  expect_match(attr(met, "notes"), "0.98 to 0.99 for 0.99")
  expect_output(print(met), "Note: The confidence level holds for a sample size fixed")
  # Enough observations but an interval still too long: no note.
  unmet <- seq_blci_dfr_median(x, d = 4.4)
  expect_false(attr(unmet, "is_satisfied"))
  expect_identical(attr(unmet, "notes"), character(0))
  # Too few observations: only the note that says how many are needed.
  few <- seq_blci_dfr_median(x[1:4], d = 100)
  expect_length(attr(few, "notes"), 1L)
  expect_match(attr(few, "notes"), "at least 6 are needed")
})

test_that("exact coverage is at least conf_level at every n, and not wastefully more", {
  # The whole range of n from the minimum to 3000, for six levels including
  # low ones where the normal approximation alone falls short.
  for (conf in c(0.50, 0.55, 0.80, 0.90, 0.95, 0.99)) {
    n_min <- .seq_blci_min_n(conf)
    n <- n_min:3000
    b <- vapply(n, .seq_blci_rank, integer(1), conf_level = conf)
    expect_true(all(b >= 1L))
    cover <- 1 - 2 * stats::pbinom(b - 1, n, 0.5)
    expect_true(all(cover >= conf))
    # Never more than one rank wider than the narrowest admissible interval:
    # moving both limits in by two ranks always loses the confidence level.
    cover_in2 <- 1 - 2 * stats::pbinom(b + 1, n, 0.5)
    expect_true(all(cover_in2 < conf))
    # Below the minimum n there is no admissible rank.
    expect_identical(.seq_blci_rank(n_min - 1L, conf), 0L)
  }
  # The exported function reports those ranks.
  for (n in c(6L, 9L, 25L, 100L, 1000L)) {
    out <- seq_blci_dfr_median(seq_len(n), d = 1)
    b <- attr(out, "ranks")[["lower"]]
    expect_identical(b, .seq_blci_rank(n, 0.95))
    expect_equal(out$value[c(1, 3)], c(b, n - b + 1))
    expect_gte(exact_coverage(b, n), 0.95)
  }
})

test_that("the rank is the exact binomial rank or one below it, never above", {
  # The exact binomial rule takes the largest rank b whose fixed n coverage
  # P(b <= Bin(n, 1/2) <= n - b) is at least conf_level, that is, the largest
  # b with P(Bin <= b - 1) <= alpha / 2. With q the alpha / 2 quantile of the
  # binomial (the smallest k with P(Bin <= k) >= alpha / 2), that b is q + 1
  # when P(Bin <= q) equals alpha / 2 exactly and q otherwise. The package's
  # rank (a continuity corrected normal approximation, lowered if short) is
  # that rank at all but a few sample sizes and one below it at the rest:
  # 7 of 1996 sample sizes at 90%, 17 of 1995 at 95%, 48 of 1993 at 99%
  # (counts from a search over the ranks that does not use qbinom). A rank
  # above the exact one would have fixed n coverage below conf_level.
  one_below <- c("0.9" = 7L, "0.95" = 17L, "0.99" = 48L)
  for (cf in names(one_below)) {
    conf <- as.numeric(cf)
    a2 <- (1 - conf) / 2
    n <- .seq_blci_min_n(conf):2000
    q <- stats::qbinom(a2, n, 0.5)
    exact <- ifelse(stats::pbinom(q, n, 0.5) <= a2, q + 1, q)
    b <- vapply(n, .seq_blci_rank, integer(1), conf_level = conf)
    expect_true(all(b <= exact))
    expect_true(all(b >= exact - 1))
    expect_identical(sum(b == exact - 1), one_below[[cf]])
    # The smallest fixed n coverage over the whole range, by summing binomial
    # probabilities: at or above the level, and within 0.0001 of it.
    cover <- vapply(seq_along(n), function(i) exact_coverage(b[i], n[i]),
                    numeric(1))
    expect_gte(min(cover), conf)
    expect_lt(min(cover), conf + 1e-4)
  }
})

test_that("the rank is lowered where the normal approximation is anti-conservative", {
  # n = 7 at 55%: the approximation gives ceiling((7 - 0.7554 * sqrt(7) - 1) / 2)
  # = 3, but P(3 <= Bin(7, 1/2) <= 4) = 70 / 128 = 0.5469 < 0.55. Rank 2 has
  # coverage 112 / 128 = 0.875.
  expect_equal(ceiling((7 - stats::qnorm(0.775) * sqrt(7) - 1) / 2), 3)
  expect_equal(exact_coverage(3, 7), 70 / 128)
  out <- seq_blci_dfr_median(1:7, d = 10, conf_level = 0.55)
  expect_identical(attr(out, "ranks"), c(lower = 2L, upper = 6L))
  expect_equal(out$value[c(1, 3)], c(2, 6))
})

test_that("ties are allowed and can give an interval of length zero", {
  out <- seq_blci_dfr_median(rep(3, 12), d = 0.5)
  expect_equal(out$value[out$term == "interval_length"], 0)
  expect_true(attr(out, "is_satisfied"))
})

# ---- pilot ------------------------------------------------------------------

test_that("the pilot is max(minimum n, m_0, ceiling((z / d)^(2 / (1 + gamma))))", {
  # 95%, d = 0.15: z / d = 13.07 -> 14.
  expect_identical(seq_blci_dfr_median_pilot(d = 0.15), 14L)
  # d = 2: z / d = 0.98 -> 1, so the default m_0 = 10 governs; with m_0 = 2
  # the minimum n for 95% (6) governs.
  expect_identical(seq_blci_dfr_median_pilot(d = 2), 10L)
  expect_identical(seq_blci_dfr_median_pilot(d = 2, m_0 = 2), 6L)
  # 99%, d = 0.05: z_.005 / d = 2.575829 / 0.05 = 51.52 -> 52.
  expect_identical(seq_blci_dfr_median_pilot(d = 0.05, conf_level = 0.99), 52L)
  # gamma = 0.49 at d = 0.15: (13.0664^2)^(1 / 1.49) = 31.50 -> 32.
  expect_identical(seq_blci_dfr_median_pilot(d = 0.15, gamma = 0.49), 32L)
  # It no longer equals the required sample size at unit scale: for normal
  # data n* = (pi / 2) (z / d)^2 sigma^2, which is 16.8 at d = 0.15 and
  # sigma = 0.25; the pilot must not exceed it.
  z <- stats::qnorm(0.975)
  expect_lt(seq_blci_dfr_median_pilot(d = 0.15), (pi / 2) * (z / 0.15)^2 * 0.25^2)
})

# ---- guards -----------------------------------------------------------------

test_that("an interval length that cannot be represented is an error, never Inf", {
  # Eight observations at a 95% level use the smallest and the largest
  # (b = 1: P(1 <= Bin(8, 1/2) <= 7) = 1 - 2 / 256 = 0.992, and b = 2 gives
  # 1 - 18 / 256 = 0.930). With a smallest value of -1.7e308 and a largest
  # of 1.7e308 the length is 3.4e308, beyond the largest double (1.797e308);
  # it was returned as Inf beside a rule that could then never be met.
  x <- c(-1.7e308, -1.6e308, 1.5e308, 1.7e308, 1.6e308, 1.5e308, 1.7e308,
         1.6e308)
  expect_error(seq_blci_dfr_median(x, d = 1),
               "'interval_length' row of the result is too large to represent \\(Inf\\)\\. Express the data in larger units \\(divide them by a constant\\) and state the precision target in the same units\\.$")
  # Half the size, the same data are an ordinary case: the length is
  # 0.85e308 - (-0.85e308) = 1.7e308, which is a double.
  out <- seq_blci_dfr_median(x / 2, d = 1)
  expect_identical(out$value[out$term == "interval_length"], 1.7e308)
  expect_identical(attr(out, "ranks"), c(lower = 1L, upper = 8L))
  expect_false(attr(out, "is_satisfied"))
  # A bound 2 d that is not a double (d beyond half of the largest double)
  # is refused as well; the width attribute was Inf and every interval met
  # the rule.
  expect_error(seq_blci_dfr_median(1:8, d = 1e308),
               "full width of the interval, twice 'd', is too large to represent \\(Inf\\)")
  expect_identical(attr(seq_blci_dfr_median(1:8, d = 8.9e307), "width"),
                   1.78e308)
})

test_that("seq_blci_dfr_median and its pilot check every argument", {
  set.seed(113)
  x <- rnorm(30)
  expect_error(seq_blci_dfr_median(x), "'d' must be supplied")
  expect_error(seq_blci_dfr_median(x, d = 0), "largest confidence interval length")
  expect_error(seq_blci_dfr_median(x, d = NA_real_), "'d' must be a single finite number")
  expect_error(seq_blci_dfr_median(x, d = Inf), "'d' must be a single finite number")
  expect_error(seq_blci_dfr_median(x, d = 0.5, conf_level = 1), "'conf_level'")
  expect_error(seq_blci_dfr_median(x, d = 0.5, conf_level = NA_real_), "'conf_level'")
  expect_error(seq_blci_dfr_median(d = 0.5), "'x' must be supplied")
  expect_error(seq_blci_dfr_median("a", d = 0.5), "'x' must be a numeric vector")
  expect_error(seq_blci_dfr_median(matrix(x, 6), d = 0.5), "not a matrix")
  expect_error(seq_blci_dfr_median(c(x, Inf), d = 0.5), "infinite values")
  expect_error(seq_blci_dfr_median(c(x, NA), d = 0.5, na.rm = FALSE), "missing values")
  expect_error(seq_blci_dfr_median(numeric(0), d = 0.5),
               "'x' needs at least 1 observation; it has 0\\. Collect more data")
  # The data argument is x, as for every one sample procedure.
  expect_identical(names(formals(seq_blci_dfr_median)),
                   c("x", "d", "conf_level", "na.rm"))
  expect_identical(seq_blci_dfr_median(x = x, d = 0.5),
                   seq_blci_dfr_median(x, 0.5))
  # The name used before 1.0.0 is R's "unused argument" error (its wording
  # is not matched because R translates it).
  expect_error(seq_blci_dfr_median(data = x, d = 0.5))
  expect_equal(seq_blci_dfr_median(c(x, NA), d = 0.5)$value,
               seq_blci_dfr_median(x, d = 0.5)$value)

  expect_error(seq_blci_dfr_median_pilot(), "'d' must be supplied")
  expect_error(seq_blci_dfr_median_pilot(d = -1), "'d' must be a single finite number")
  expect_error(seq_blci_dfr_median_pilot(d = 0), "largest confidence interval length")
  expect_error(seq_blci_dfr_median_pilot(d = 0.3, conf_level = 0), "'conf_level'")
  expect_error(seq_blci_dfr_median_pilot(d = 0.3, m_0 = 1), "'m_0'")
  expect_error(seq_blci_dfr_median_pilot(d = 0.3, m_0 = NA), "'m_0'")
  expect_error(seq_blci_dfr_median_pilot(d = 0.3, gamma = -1), "'gamma'")
  expect_error(seq_blci_dfr_median_pilot(d = 1e-20),
               "implied pilot sample size .* is too large to represent\\. Check 'd' \\(half of the largest confidence interval length allowed\\)\\.$")
  expect_identical(names(formals(seq_blci_dfr_median_pilot)),
                   c("d", "conf_level", "m_0", "gamma"))
  expect_type(seq_blci_dfr_median_pilot(d = 0.3), "integer")
})

# ---- operating characteristics ----------------------------------------------

# One complete sequential run at 95% from a pilot of m observations, with the
# stopping rule coded independently of the package: the sample is kept sorted
# by insertion, and the rank is the continuity corrected normal approximation
# (which a test above shows has exact coverage of at least 0.95 at every n).
blci_run <- function(m, d, rgen) {
  z <- stats::qnorm(0.975)
  x <- rgen(m)
  s <- sort(x)
  n <- m
  repeat {
    b <- max(1, ceiling((n - z * sqrt(n) - 1) / 2))
    if (s[n - b + 1] - s[b] <= 2 * d) break
    x_new <- rgen(1)
    x <- c(x, x_new)
    s <- append(s, x_new, after = findInterval(x_new, s))
    n <- n + 1L
  }
  list(x = x, N = n, lower = s[b], upper = s[n - b + 1])
}

# Does the package agree with a run? Met at the stopping point with the same
# limits, and not met one observation earlier (when beyond the pilot).
blci_agrees <- function(run, m, d) {
  out <- seq_blci_dfr_median(run$x, d = d)
  isTRUE(attr(out, "is_satisfied")) &&
    isTRUE(all.equal(out$value[c(1, 3)], c(run$lower, run$upper))) &&
    (run$N == m || identical(attr(seq_blci_dfr_median(run$x[-run$N], d = d),
                                  "is_satisfied"), FALSE))
}

test_that("the sequential median interval performs as documented at a moderate target", {
  # N(0, 1) data, 95%, d chosen so that n* = (pi / 2) (z / d)^2 = 100. The
  # complete procedure is run 1,500 times from the package's pilot. The
  # package function confirms the first 40 runs: met at the stopping point
  # with the same limits and a reported length within the bound, and not met
  # one observation earlier.
  set.seed(113)
  z <- stats::qnorm(0.975)
  d <- z * sqrt(pi / 2) / 10
  m <- seq_blci_dfr_median_pilot(d)
  expect_identical(m, 10L)
  reps <- 1500L
  N <- integer(reps)
  cover <- logical(reps)
  pkg_len <- numeric(40L)
  agree <- TRUE
  for (i in seq_len(reps)) {
    run <- blci_run(m, d, rnorm)
    N[i] <- run$N
    cover[i] <- run$lower <= 0 && 0 <= run$upper
    if (i <= 40L) {
      out <- seq_blci_dfr_median(run$x, d = d)
      pkg_len[i] <- out$value[out$term == "interval_length"]
      agree <- agree && blci_agrees(run, m, d)
    }
  }
  expect_true(agree)
  # The interval the package reports at the stopping point is within the
  # length bound.
  expect_true(all(pkg_len > 0 & pkg_len <= 2 * d))
  # Coverage of the sequential procedure is below nominal: the help page
  # documents roughly 0.92 to 0.94 (three independent simulations of 20,000
  # runs gave 0.921, 0.923, and 0.928 at this target, and one of 40,000 runs
  # gave 0.926). SE from 1,500 runs is 0.0067; the band is the documented
  # range widened by 4 SE.
  expect_gt(mean(cover), 0.92 - 4 * 0.0067)
  expect_lt(mean(cover), 0.94 + 4 * 0.0067)
  # Mean final sample size: independent simulations gave 0.968 n* and
  # 0.971 n* with sd(N) of 39, so the SE of the mean is 1.0 (0.010 n*) and
  # the tolerance is 5 SE. A one-sided quantile would give about 0.70 n*
  # and a rank one lower 1.19 n*.
  expect_lt(abs(mean(N) / 100 - 0.968), 0.05)
})

test_that("the median interval's pilot does not presume a scale: sd 0.25, 1, and 3", {
  # N(0, sd^2) data, 95%, with d = 0.1535 chosen so that
  # n* = (pi / 2) (z / d)^2 sd^2 is 16 at sd 0.25, hence 256 at sd 1 and 2304
  # at sd 3. The pilot is 13 whatever the scale (the pilot this replaced,
  # (z / d)^2 + 1 = 164, was ten times the 16 needed at sd 0.25). An
  # independent simulation (20,000 runs; 3,000 at sd 3) gave
  # E[N] / n* = 1.32, 0.95, and 0.96; the first is above 1 because the pilot
  # is most of the 16 needed. sd(N) there was 7.3, 83, and 444, so with 400,
  # 200, and 40 runs four standard errors of E[N] / n* are 0.091, 0.092, and
  # 0.122; the tolerances are 0.10, 0.10, and 0.13.
  z <- stats::qnorm(0.975)
  d <- z * sqrt(pi / 2) * 0.25 / sqrt(16)
  m <- seq_blci_dfr_median_pilot(d)
  expect_identical(m, 13L)
  scale_leg <- function(sd, reps, expected, tol) {
    n_star <- (pi / 2) * (z / d)^2 * sd^2
    N <- integer(reps)
    agree <- TRUE
    for (i in seq_len(reps)) {
      run <- blci_run(m, d, function(k) rnorm(k, 0, sd))
      N[i] <- run$N
      if (i <= 5L) agree <- agree && blci_agrees(run, m, d)
    }
    expect_true(agree)
    expect_lt(abs(mean(N) / n_star - expected), tol)
  }
  set.seed(113)
  scale_leg(0.25, 400L, 1.32, 0.10)
  scale_leg(1, 200L, 0.95, 0.10)
  # The third scale needs about 2,200 observations a run; it is skipped on
  # CRAN for time only.
  skip_on_cran()
  scale_leg(3, 40L, 0.96, 0.13)
})
