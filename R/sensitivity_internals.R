# Internal helpers shared by the seq_aipe_*_sensitivity Monte Carlo
# functions. Not exported. Each function supplies a data generating closure
# (draw) with a known true value and a check closure (evaluate); the engine
# runs the sequential procedure to stopping G times and summarizes the
# realized sample size, interval width, and coverage.

# Standardized chi square with 3 degrees of freedom: mean 0, variance 1,
# right skewed (skewness sqrt(8 / 3), about 1.63). The "skewed" data
# generating mechanism.
#' @noRd
.skewed_std <- function(k) (stats::rchisq(k, df = 3) - 3) / sqrt(6)

# Evaluate a thunk under a temporarily set seed. The caller's random number
# generator state (including the case of a session in which no seed has yet
# been created) is restored afterward, on error too. Saving and restoring is
# delegated to withr; the package itself never reads or writes the
# generator's state. A NULL seed is a no-op: the draws come from the user's
# current stream.
#' @noRd
.with_seed <- function(seed, thunk) {
  seed <- .check_seed(seed)
  if (is.null(seed)) return(thunk())
  withr::with_seed(seed, thunk())
}

# A required argument of a sensitivity function that was not supplied: the
# package's own message, which says what the argument is, in place of R's
# 'argument "true_value" is missing, with no default'. Call it first, as
# .sensitivity_need(true_value, "true_value", "..."); `value` is only
# tested for being missing, never evaluated.
#' @noRd
.sensitivity_need <- function(value, name, what) {
  if (missing(value)) {
    .smsd_stop("'", name, "' must be supplied: ", what, ".")
  }
  invisible(TRUE)
}

# Check the scalar arguments every sensitivity function shares and return
# them cleaned (doubles for true_value, width, and conf_level; integers for
# G, m_0, max_n, and batch; an integer or NULL for seed). A function with a
# restricted true value (a correlation, a coefficient of variation, a Gini
# index) checks that range itself before calling this. `width_range` is the
# length of the range of an effect size confined to one (2 for a
# correlation, 1 for the Gini index): a target width that large is refused,
# as the check and pilot functions refuse it.
#' @noRd
.sensitivity_check_args <- function(true_value, width, conf_level, G, m_0,
                                    max_n = 1e5L, batch = 1L, seed = NULL,
                                    width_range = Inf) {
  list(
    true_value = .check_number(true_value, "true_value"),
    width = .aipe_check_width_bounded(width, width_range),
    conf_level = .check_conf_level(conf_level),
    G = .check_count(G, "G"),
    m_0 = .check_count(m_0, "m_0", min = 4L),
    max_n = .check_count(max_n, "max_n"),
    batch = .check_count(batch, "batch"),
    seed = .check_seed(seed)
  )
}

# The data generating distribution must be one that the procedure can be run
# on. Without the checks below an extreme population value gives data that
# are infinite, missing, or without variation that a check function can use,
# and the user sees an error from a check function about a data argument
# they never supplied, a study in which no replication can stop, or a
# summary that describes the package's tolerance for constant data and not
# the population asked for.
#
# A parameter of the generating distribution that is computed from the
# requested population value (a gamma shape, a lognormal log standard
# deviation) must be a positive finite number. `distribution` names the
# distribution, `what` the parameter and its formula.
#' @noRd
.sensitivity_check_parameter <- function(parameter, distribution, what,
                                         true_value) {
  if (is.finite(parameter) && parameter > 0) return(parameter)
  .smsd_stop("'true_value' is too extreme to simulate: the ", distribution,
             " distribution would need a ", what, ", of ", format(parameter),
             ", and that must be a positive finite number; got ",
             .smsd_describe(true_value), ".")
}

# Every generated variable must vary by the package's own standard. The
# check functions take a sample whose values agree to ten significant digits
# to have no variation (.is_constant(): a range of no more than
# .smsd_constant_tol times the largest absolute value). Depending on the
# procedure, such a sample gives no estimate, or its group adds nothing to
# the variance of a contrast. A population whose standard deviation is that
# small beside its mean would therefore give a study in which no replication
# can stop, or one whose intervals leave a group's variance out, and the
# summary would describe the tolerance, not the population.
#
# The sensitivity functions refuse such a population, with a margin: the
# absolute value of the mean of a generated variable must be less than
# .sensitivity_mean_limit() times its standard deviation, where the limit is
# a hundredth of the reciprocal of the tolerance (1e8 for the tolerance of
# 1e-10). At that limit the tolerance is a range of 0.01 standard
# deviations, and a sample of 4, the smallest pilot sample, has so small a
# range with probability 1.3e-7 for normal data and 3.7e-7 for the skewed
# mechanism (numerical integration of the distribution of the range; a
# sample of 5 is at 6e-10), so the tolerance plays no part in a study of any
# practical size. Without the margin (a mean of 1e10 standard deviations)
# the probability is 0.11.
#' @noRd
.sensitivity_mean_limit <- function() 0.01 / .smsd_constant_tol

# Which generated variables have too little variation: an absolute mean of
# .sensitivity_mean_limit() standard deviations or more.
#' @noRd
.sensitivity_too_constant <- function(mean, sd) {
  abs(mean) >= .sensitivity_mean_limit() * sd
}

# Why such a population is refused: the closing sentence of the messages.
#' @noRd
.sensitivity_why_too_constant <- function() {
  paste0("The check functions take a sample whose values agree to ten ",
         "significant digits to have no variation, and the generated data ",
         "must stay well clear of that.")
}

# A mechanism that adds noise with standard deviation `scale` to a location
# (the standardized mean difference, the standardized mean, the contrast):
# the mean of the generated variable is `location`. `name` is "true_value"
# (unit variance noise) or "true_means" (one location and one scale per
# group).
#' @noRd
.sensitivity_check_location <- function(location, scale = 1,
                                        name = "true_value") {
  bad <- .sensitivity_too_constant(location, scale)
  if (!any(bad)) return(invisible(location))
  limit <- format(.sensitivity_mean_limit())
  if (name == "true_value") {
    .smsd_stop("'true_value' must be less than ", limit, " in absolute ",
               "value; got ", .smsd_describe(location), ". The generated ",
               "data have unit variance, which is too little beside a ",
               "larger value. ", .sensitivity_why_too_constant())
  }
  .smsd_stop("'", name, "' must each be less than ", limit, " standard ",
             "deviations ('true_sds') in absolute value; ", sum(bad), " of ",
             length(bad), if (sum(bad) == 1L) " is" else " are", " not. ",
             .sensitivity_why_too_constant())
}

# A mechanism that generates positive data whose standard deviation is
# `rel_sd` times their mean (the coefficient of variation, the Gini index):
# `rel_sd` must be more than 1 / .sensitivity_mean_limit(). `least` is the
# smallest population value that gives such data, in words.
#' @noRd
.sensitivity_check_rel_sd <- function(rel_sd, least, true_value) {
  if (!isTRUE(.sensitivity_too_constant(1, rel_sd))) {
    return(invisible(true_value))
  }
  .smsd_stop("'true_value' must be greater than ", least, "; got ",
             .smsd_describe(true_value), ". ",
             .sensitivity_why_too_constant())
}

# The regression mechanism y = beta x + e has variables with mean zero, so
# both vary whatever the slope. What can be lost is the error e beside
# beta x: from |beta| = 1 / .Machine$double.eps on, adjacent doubles near
# beta x are at least one error standard deviation apart.
#' @noRd
.sensitivity_check_slope <- function(beta) {
  if (abs(beta) < 1 / .Machine$double.eps) return(invisible(beta))
  .smsd_stop("'true_value' must be less than 1 / .Machine$double.eps ",
             "(about ", format(1 / .Machine$double.eps, digits = 2),
             ") in absolute value; got ", .smsd_describe(beta), ". Beside ",
             "a larger slope the unit variance of the errors cannot be ",
             "represented.")
}

# Standardized (mean 0, variance 1) draws for the "normal" and "skewed"
# mechanisms: a function of the number of observations.
#' @noRd
.sensitivity_std <- function(distribution) {
  if (distribution == "normal") stats::rnorm else .skewed_std
}

# Parameter maps from a requested population value to the generating
# distribution. Each is pinned in the tests against the population value
# obtained by numerical integration of the generating density.
#
# Gamma with shape k and any scale has coefficient of variation 1 / sqrt(k).
#' @noRd
.sensitivity_gamma_shape_cv <- function(cv) 1 / cv^2

# Lognormal with log standard deviation s has coefficient of variation
# sqrt(exp(s^2) - 1), so s = sqrt(log(1 + cv^2)). It is computed with
# log1p(), which is that logarithm without the rounding of 1 + cv^2: with
# log(1 + cv^2) the generated data have a coefficient of variation that is
# off by 4e-5 of its value at 1e-6 and by a third just above 1e-8.
#' @noRd
.sensitivity_sdlog_cv <- function(cv) sqrt(log1p(cv^2))

# Lognormal with log standard deviation s has Gini index
# 2 * pnorm(s / sqrt(2)) - 1.
#' @noRd
.sensitivity_gini_lognormal <- function(sdlog) {
  2 * stats::pnorm(sdlog / sqrt(2)) - 1
}

# The log standard deviation that gives a requested Gini index: the inverse
# of .sensitivity_gini_lognormal().
#' @noRd
.sensitivity_sdlog_gini <- function(gini) sqrt(2) * stats::qnorm((gini + 1) / 2)

# One replication: draw the pilot sample, then check the stopping rule,
# adding `batch` observations (per group, where the design has groups)
# between checks until the rule is met or the sample reaches max_n. The last
# batch is shortened so that max_n is an exact cap for any batch size.
#
# batch = 1 is the purely sequential procedure of the articles. Jumping ahead
# to the current criterion instead would overshoot the first crossing time,
# because the criterion moves as data accrue.
#
# A check counts as a stop only when attr(res, "is_satisfied") is TRUE. A
# check at which no variance has been estimated (the raw estimate is zero,
# negative, or not finite; the check functions report it as "not met", with
# attr(res, "variance_unusable") = TRUE) therefore leads to more sampling,
# never to a stop.
#
# Returns the realized sample size, the full width of the final interval,
# and three 0/1 indicators: the interval covers true_value; the replication
# was stopped by the cap rather than by the rule; at least one check, the
# last included, reported an unusable variance estimate.
#
# The final check can be one that gives no interval: its limits are NA
# because the point estimate is not finite or the variance estimate is
# unusable. Such a check is never reported as met, so it can be the final
# one only when the cap ends the replication. The replication then has no
# width (NA) and does not cover; .sensitivity_run() says what the summary
# makes of it.
#' @noRd
.sensitivity_one <- function(pilot_n, draw, evaluate, true_value, max_n,
                             batch = 1L) {
  data <- draw(NULL, pilot_n)
  n <- pilot_n
  unusable <- FALSE
  repeat {
    res <- evaluate(data)
    if (isTRUE(attr(res, "variance_unusable"))) unusable <- TRUE
    stopped <- isTRUE(attr(res, "is_satisfied"))
    if (stopped || n >= max_n) break
    k <- min(batch, max_n - n)
    data <- draw(data, k)
    n <- n + k
  }
  lower <- res$value[match("lower_limit", res$term)]
  upper <- res$value[match("upper_limit", res$term)]
  full_width <- upper - lower
  c(n = n,
    width = if (is.finite(full_width)) full_width else NA_real_,
    covered = as.numeric(isTRUE(lower <= true_value && true_value <= upper)),
    capped = as.numeric(!stopped),
    unusable = as.numeric(unusable))
}

# Run the Monte Carlo study and return the term/value summary, an
# smsd_sensitivity object.
#
# draw(data, k) returns the data with k more observations (per group)
# appended, or a fresh sample of size k when data is NULL; evaluate(data)
# returns the smsd_check for the data so far. pilot_n is the sample size of
# the first check. G, max_n, and batch are checked here as well as in
# .sensitivity_check_args(), and the seed by .with_seed(), so the engine is
# safe for any caller.
#
# What the summary reports when a replication ends without an interval (see
# .sensitivity_one(); it can happen only at the cap). Every row stays a
# statement about all G replications or says which ones it covers:
#   - mean_n, median_n, sd_n, pct_hit_max_n, and pct_variance_unusable count
#     the replication like any other;
#   - coverage counts it as not covering, and pct_width_within as not within
#     the target, so both remain proportions of all G replications (an
#     interval that was never obtained neither contains the true value nor
#     meets the width);
#   - mean_width and median_width are over the replications that ended with
#     an interval, and are NA when none did. No other row can be NA on this
#     account.
#' @noRd
.sensitivity_run <- function(G, seed, pilot_n, draw, evaluate, true_value,
                             width, conf_level, max_n, batch = 1L) {
  G <- .check_count(G, "G")
  max_n <- .check_count(max_n, "max_n")
  batch <- .check_count(batch, "batch")
  pilot_n <- as.integer(pilot_n)
  if (is.na(pilot_n)) {
    .smsd_stop("The pilot sample size implied by 'width', 'conf_level', and ",
               "'m_0' is too large to represent. Use a wider 'width'.")
  }
  if (max_n < pilot_n) {
    .smsd_stop("'max_n' (", max_n, ") is below the pilot sample size (",
               pilot_n, ") implied by 'width', 'conf_level', and 'm_0'. ",
               "Every replication begins with the pilot sample, so the cap ",
               "cannot be smaller than it. Raise 'max_n' to at least ",
               pilot_n, ", or use a wider 'width'.")
  }
  mat <- .with_seed(seed, function() {
    vapply(seq_len(G), function(i) {
      .sensitivity_one(pilot_n, draw, evaluate, true_value, max_n, batch)
    }, numeric(5L))
  })
  n <- mat["n", ]
  w <- mat["width", ]
  has_interval <- !is.na(w)
  final_width <- w[has_interval]
  # Proportions and percentages are formed from counts, so that, for
  # example, 950 covering intervals in 1000 give exactly 0.95.
  new_smsd_sensitivity(
    term = c("true_value", "target_width", "conf_level", "replications",
             "mean_n", "median_n", "sd_n",
             "mean_width", "median_width", "pct_width_within", "coverage",
             "pct_hit_max_n", "pct_variance_unusable"),
    value = c(true_value, width, conf_level, G,
              mean(n), stats::median(n), stats::sd(n),
              if (any(has_interval)) mean(final_width) else NA_real_,
              if (any(has_interval)) stats::median(final_width) else NA_real_,
              100 * sum(final_width <= width) / G,
              sum(mat["covered", ]) / G,
              # Replications stopped by the max_n cap rather than by the
              # rule: when this is above zero the other summaries describe
              # a truncated procedure and should be read with care.
              100 * sum(mat["capped", ]) / G,
              # Replications in which at least one check reported that no
              # variance had been estimated (and was therefore not a stop).
              100 * sum(mat["unusable", ]) / G)
  )
}
