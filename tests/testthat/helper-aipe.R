# Helpers for the tests of the sequential accuracy in parameter estimation
# (seq_aipe_*) family.

# Standardized chi square with 3 degrees of freedom: mean 0, variance 1,
# third central moment sqrt(8/3), fourth central moment 7.
aipe_skewed <- function(k) (stats::rchisq(k, df = 3) - 3) / sqrt(6)

# Bivariate normal pairs with correlation rho, as a list(x, y).
aipe_pairs <- function(k, rho) {
  x <- stats::rnorm(k)
  list(x = x, y = rho * x + sqrt(1 - rho^2) * stats::rnorm(k))
}

# The stopping rule as printed in the source articles (Kelley, Darku, and
# Chattopadhyay, 2018, Equation 22; 2019, Equation 10), written here
# independently of the package: is n at least
# ceiling((2 z / width)^2 (xi2 + 1 / n))?
aipe_rule_met <- function(n, xi2, width, conf_level = 0.95) {
  # The quantile is taken from the upper tail, which stays accurate as
  # conf_level approaches 1 (qnorm(1 - alpha / 2) does not).
  z <- stats::qnorm((1 - conf_level) / 2, lower.tail = FALSE)
  n >= ceiling((2 * z / width)^2 * (xi2 + 1 / n))
}

# The bounded variance estimate of the articles, max{raw, n^-3} (for example
# Kelley, Darku, and Chattopadhyay, 2018, Equation 27), written here
# independently of the package; the bound alone when the raw value is not a
# number.
aipe_bounded <- function(raw, n) {
  if (is.finite(raw)) max(raw, n^(-3)) else n^(-3)
}

# How the windows of the seeded operating characteristic tests are set. Each
# test compares a mean over a few replications with a reference value that
# does not come from the package (a published table row, or a theoretical
# sample size or the nominal coverage). The half width of the window is
#
#   4 standard errors of the tested mean + |reference - long run value|,
#
# where the long run value is this package's own mean over thousands of
# replications, so that the window reaches at least four standard errors on
# both sides of what the procedure actually delivers, not only of the
# reference. The standard error is the standard deviation of the tested
# mean over 200 other seeds of the same test (or the binomial standard
# error for a coverage), or the value implied by a larger run when that is
# bigger. Each test states its numbers. The long runs were made when a
# positive variance estimate below n^-3 also kept the rule from being met.
# In the conditions tested against a window that event did not occur in 300
# to 4,000 further replications of each, and the means of those
# replications agree with the long run values quoted in the tests, so the
# values stand.

# Run one sequential study to its stopping point, adding one observation
# (per group) at each stage, as in the simulation studies of the source
# articles.
#
#   data:  the whole stream of observations (longer than will be needed);
#   take:  function(data, n) returning the first n observations (per group);
#   stat:  function(d) returning c(estimate, raw), the point estimate and
#          the RAW variance estimate (floor = FALSE) from the package's
#          internal estimators;
#   check: function(d) calling the exported check function.
#
# Calling the exported function at every stage costs several times as much
# as the estimator alone, so the stopping time is located with `stat` and
# the rule above, and then confirmed with the exported function: it must
# report the rule as met at the stopping point, with the same estimate and
# variance, and as not met one stage earlier. (That the exported function
# agrees with this search at every stage is tested separately on whole
# paths.)
#
# Returns c(n, lower, upper, n_printed): the stopping time and interval of
# the procedure as shipped, which does not stop while the raw variance
# estimate is zero, negative, or not finite, or the estimate is not finite,
# and the stopping time of the rule as printed in the articles, which does
# (n_printed <= n, and the two are equal unless the printed rule stopped on
# such an estimate).
aipe_run <- function(data, take, stat, check, width, pilot,
                     conf_level = 0.95, max_n = 5000L) {
  n <- pilot
  n_printed <- NA_integer_
  repeat {
    s <- stat(take(data, n))
    xi2 <- aipe_bounded(s[2], n)
    met <- aipe_rule_met(n, xi2, width, conf_level)
    if (met && is.na(n_printed)) n_printed <- n
    usable <- is.finite(s[1]) && is.finite(s[2]) && s[2] > 0
    if (usable && met) break
    n <- n + 1L
    if (n > max_n) stop("the stopping rule was not met by max_n")
  }
  out <- check(take(data, n))
  v <- stats::setNames(out$value, out$term)
  ok <- isTRUE(attr(out, "is_satisfied")) &&
    isTRUE(all.equal(unname(v[2]), s[1])) &&
    isTRUE(all.equal(v[["xi_squared"]], xi2))
  if (ok && n > pilot) {
    ok <- !isTRUE(attr(check(take(data, n - 1L)), "is_satisfied"))
  }
  if (!ok) stop("the exported check function disagrees with the search")
  c(n = n, lower = v[["lower_limit"]], upper = v[["upper_limit"]],
    n_printed = n_printed)
}

# Summaries of a matrix of aipe_run() results (one column per replication).
aipe_summary <- function(runs, truth) {
  c(mean_n = mean(runs["n", ]),
    se_n = stats::sd(runs["n", ]) / sqrt(ncol(runs)),
    coverage = mean(runs["lower", ] <= truth & truth <= runs["upper", ]),
    mean_width = mean(runs["upper", ] - runs["lower", ]),
    max_width = max(runs["upper", ] - runs["lower", ]),
    mean_n_printed = mean(runs["n_printed", ]))
}

# Collect the warnings raised by an expression (and its value).
aipe_warnings <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

# The unbiased fourth central moment estimator as printed: Kelley, Darku,
# and Chattopadhyay (2018, Equations 76 and 79; 2019, Equations 71 and 72).
mu4_as_printed <- function(x) {
  n <- length(x)
  (n^2 * sum((x - mean(x))^4) - (2 * n - 3) * sum(x^4) +
     (8 * n - 12) * mean(x) * sum(x^3) - (6 * n - 9) / n * sum(x^2)^2) /
    ((n - 1) * (n - 2) * (n - 3))
}

# The bivariate k-statistics as printed in Kelley, Darku, and Chattopadhyay
# (2019, Equation 74), in the raw power sums W_pq = sum(x^p y^q). The
# coefficient of W_12 W_10 in k22 is 2 (n + 1) / n, as symmetry in x and y
# requires (the article prints (n + 1) / n there).
k_as_printed <- function(x, y) {
  n <- length(x)
  W <- function(p, q) sum(x^p * y^q)
  cst <- n / ((n - 1) * (n - 2) * (n - 3))
  k13 <- cst * ((n + 1) * W(1, 3) - (n + 1) / n * W(0, 3) * W(1, 0) -
                  3 * (n - 1) / n * W(1, 1) * W(0, 2) -
                  3 * (n + 1) / n * W(1, 2) * W(0, 1) +
                  6 / n * W(1, 1) * W(0, 1)^2 +
                  6 / n * W(0, 2) * W(0, 1) * W(1, 0) -
                  6 / n^2 * W(1, 0) * W(0, 1)^3)
  k22 <- cst * ((n + 1) * W(2, 2) - 2 * (n + 1) / n * W(2, 1) * W(0, 1) -
                  2 * (n + 1) / n * W(1, 2) * W(1, 0) -
                  (n - 1) / n * W(2, 0) * W(0, 2) -
                  2 * (n - 1) / n * W(1, 1)^2 +
                  8 / n * W(1, 1) * W(0, 1) * W(1, 0) +
                  2 / n * W(0, 2) * W(1, 0)^2 + 2 / n * W(2, 0) * W(0, 1)^2 -
                  6 / n^2 * W(1, 0)^2 * W(0, 1)^2)
  c(k13 = k13, k22 = k22)
}

# Walk the exported check function over every stage from the pilot until it
# reports the rule as met, and return that sample size. Used to confirm
# that aipe_run() and the exported function agree on whole paths.
aipe_walk <- function(data, take, check, pilot, max_n = 5000L) {
  n <- pilot
  repeat {
    out <- check(take(data, n))
    if (!identical(attr(out, "is_satisfied"),
                   attr(out, "current_n") >= attr(out, "criterion_n"))) {
      stop("is_satisfied disagrees with current_n >= criterion_n")
    }
    if (attr(out, "is_satisfied")) return(n)
    n <- n + 1L
    if (n > max_n) stop("the stopping rule was not met by max_n")
  }
}

# Expectations shared by the argument checks of every pilot function.
expect_pilot_checks <- function(pilot) {
  expect_error(pilot(), "'width' must be supplied")
  expect_error(pilot(width = 0), "'width' must be a single finite number")
  expect_error(pilot(width = -1), "'width'")
  expect_error(pilot(width = NA_real_), "'width'.*got NA")
  expect_error(pilot(width = Inf), "'width'")
  expect_error(pilot(width = c(0.1, 0.2)), "'width'.*vector of length 2")
  expect_error(pilot(width = "0.2"), "'width'")
  expect_error(pilot(width = 0.3, conf_level = 1), "'conf_level'")
  expect_error(pilot(width = 0.3, conf_level = 0), "'conf_level'")
  expect_error(pilot(width = 0.3, conf_level = 1.5), "'conf_level'")
  expect_error(pilot(width = 0.3, conf_level = NA_real_),
               "'conf_level'.*got NA")
  expect_error(pilot(width = 0.3, m_0 = 2), "'m_0'.*at least 4")
  expect_error(pilot(width = 0.3, m_0 = 4.5), "'m_0'.*whole number")
  expect_error(pilot(width = 0.3, m_0 = NA_real_), "'m_0'")
  expect_error(pilot(width = 0.3, m_0 = Inf), "'m_0'")
  # A width so small that the pilot cannot be represented is an error, not
  # an NA with a coercion warning, and the message points at the argument
  # that caused it.
  expect_error(pilot(width = 1e-12),
               paste0("pilot sample size .* is too large to represent\\. ",
                      "Check 'width'"))
}

# Expectations shared by the scalar argument checks of every check function;
# `call` takes width and conf_level and supplies valid data.
expect_scalar_checks <- function(call) {
  expect_error(call(), "'width' must be supplied")
  expect_error(call(width = NA_real_), "'width'.*got NA")
  expect_error(call(width = 0), "'width'")
  expect_error(call(width = Inf), "'width'")
  expect_error(call(width = c(0.2, 0.3)), "'width'")
  expect_error(call(width = "a"), "'width'")
  expect_error(call(width = 0.3, conf_level = NA_real_), "'conf_level'")
  expect_error(call(width = 0.3, conf_level = 1), "'conf_level'")
  expect_error(call(width = 0.3, conf_level = 0), "'conf_level'")
  expect_error(call(width = 0.3, conf_level = c(0.9, 0.95)), "'conf_level'")
}

# Expectations shared by the checks of one data vector; `call` takes the
# vector and supplies everything else. `name` is the argument's name. For
# paired arguments a shortened vector is a length mismatch first, so the too
# few observations check is tested separately (short = FALSE).
expect_data_checks <- function(call, name, good, short = TRUE) {
  q <- function(msg) paste0("'", name, "'.*", msg)
  expect_error(call(letters[1:10]), q("must be a numeric vector"))
  expect_error(call(rep(c(TRUE, FALSE), 5)), q("must be a numeric vector"))
  expect_error(call(matrix(good, ncol = 2)), q("not a matrix"))
  expect_error(call(data.frame(v = good)), q("not a data frame"))
  expect_error(call(replace(good, 2, NA)), q("missing values"))
  expect_error(call(replace(good, 2, NaN)), q("missing values"))
  expect_error(call(replace(good, 2, Inf)), q("infinite"))
  if (short) expect_error(call(good[1:3]), q("at least 4"))
}
