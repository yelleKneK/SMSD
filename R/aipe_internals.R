# Internal helpers shared by the sequential accuracy in parameter estimation
# (seq_aipe_*) pilot and check functions. Not exported. Scalar and data
# argument checking is delegated to the shared checks in R/checks.R; the
# helpers here add only what is specific to this family (paired data,
# balanced groups, contrast weights) and assemble the returned object in one
# place, so that the stopping rule, the interval, and the treatment of an
# unusable variance estimate cannot drift apart across the seven effect
# sizes. A sample with no variation is recognized with the package wide
# .is_constant() (values that agree to ten significant digits), never with
# an exact comparison: arithmetic on equal values leaves noise, and under
# the rule in R/var_floor.R a variance "estimate" of 1e-30 computed from
# that noise would count as positive and be used.

# Pilot sample size shared by the family: max{m_0, ceiling(2 z / width)},
# Kelley, Darku, and Chattopadhyay (2018, Equation 21) and (2019,
# Equation 9); for the Gini index, Chattopadhyay and De (2016), whose half
# width d is width / 2. `range` is the length of the range of an effect size
# that is confined to one (see .aipe_check_width_bounded()).
#' @noRd
.aipe_pilot <- function(width, conf_level, m_0, range = Inf) {
  width <- .aipe_check_width_bounded(width, range)
  conf_level <- .check_conf_level(conf_level)
  m_0 <- .check_count(m_0, "m_0", min = 4L)
  z <- .aipe_z(conf_level)
  max(m_0, .as_sample_size(2 * z / width, "pilot sample size",
                           advice = .aipe_advice_width()))
}

# The closing sentence of the message when a pilot or criterion sample size
# cannot be represented: in this family that can only come from the width.
#' @noRd
.aipe_advice_width <- function() {
  paste("Check 'width': it is the full width of the confidence interval,",
        "on the scale of the effect size.")
}

# The upper alpha / 2 quantile of the standard normal distribution, with
# alpha = 1 - conf_level. It is taken from the upper tail at alpha / 2
# rather than as qnorm(1 - alpha / 2): the two agree at ordinary confidence
# levels, but 1 - alpha / 2 rounds to 1 as conf_level approaches 1, and the
# quantile then loses accuracy and finally becomes infinite.
#' @noRd
.aipe_z <- function(conf_level) {
  stats::qnorm((1 - conf_level) / 2, lower.tail = FALSE)
}

# One data vector of a check function (which has no na.rm argument). Returns
# the vector as doubles.
#' @noRd
.aipe_check_data <- function(x, name, min_n = 4L, ...) {
  if (missing(x)) {
    .smsd_stop("'", name, "' must be supplied: the data collected so far.")
  }
  .check_data(x, name, min_n = min_n, na.rm = NULL, ...)
}

# Two paired vectors (observed on the same units). Returns list(x, y) as
# doubles. With na.rm = TRUE a pair is dropped when either member is
# missing; with na.rm = FALSE or NULL a missing value is an error (the
# wording is that of .check_data() for the same na.rm). The too few
# observations message closes with the advice of .check_data() for a
# sequential check; a plain estimator (slr_slope()) sets
# `plain_estimator = TRUE` and gets the message without it.
#' @noRd
.aipe_check_pair <- function(x, y, x_name, y_name, min_n = 4L, na.rm = NULL,
                             plain_estimator = FALSE) {
  if (missing(x)) {
    .smsd_stop("'", x_name, "' must be supplied.")
  }
  if (missing(y)) {
    .smsd_stop("'", y_name, "' must be supplied.")
  }
  drop_pairs <- isTRUE(na.rm)
  # Type, shape, infinite values, and (unless pairs are to be dropped)
  # missing values, one vector at a time.
  .check_data(x, x_name, min_n = 0L, na.rm = if (drop_pairs) TRUE else na.rm)
  .check_data(y, y_name, min_n = 0L, na.rm = if (drop_pairs) TRUE else na.rm)
  if (length(x) != length(y)) {
    .smsd_stop("'", x_name, "' and '", y_name, "' must have the same ",
               "length, because they are paired; they have ", length(x),
               " and ", length(y), " values.")
  }
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (drop_pairs) {
    keep <- !(is.na(x) | is.na(y))
    x <- x[keep]
    y <- y[keep]
  }
  enough <- function(v, name) {
    if (plain_estimator) {
      .check_data(v, name, min_n = min_n, na.rm = NULL, advice = NULL)
    } else {
      .check_data(v, name, min_n = min_n, na.rm = NULL)
    }
  }
  list(x = enough(x, x_name), y = enough(y, y_name))
}

# Contrast weights: one finite number per group (the shared .check_vector()),
# and, which is specific to a contrast, not all zero. Returns the weights as
# doubles.
#' @noRd
.aipe_check_weights <- function(contrast_weights, n_groups) {
  if (missing(contrast_weights)) {
    .smsd_stop("'contrast_weights' must be supplied: one weight per group.")
  }
  contrast_weights <- .check_vector(contrast_weights, "contrast_weights",
                                    what = "one weight per group",
                                    len = n_groups,
                                    len_what = "one per group in 'groups'")
  if (all(contrast_weights == 0)) {
    .smsd_stop("'contrast_weights' must not all be zero: a contrast with ",
               "every weight zero estimates nothing.")
  }
  contrast_weights
}

# Deviations from the sample mean, computed in two passes. One pass leaves
# the deviations with a common offset of the order of the rounding error of
# the mean, which matters when the mean is many orders of magnitude larger
# than the standard deviation (a timestamp in seconds, say); subtracting the
# mean of the deviations removes it. The moment estimators are built on
# these deviations, which is what makes them location invariant to working
# precision.
#' @noRd
.aipe_center <- function(x) {
  d <- x - mean(x)
  d - mean(d)
}

# A power of two of the order of the largest absolute value in the data
# (1 when every value is zero). The estimators divide the data by it before
# forming sums of powers, which keeps those sums inside the range of a
# double whatever the units of the data. Without it the sixth power of the
# pooled standard deviation (standardized mean difference) and the eighth
# power of the predictor's standard deviation (slope) overflow for data of
# magnitude 1e52 and 1e39 and underflow at 1e-54 and 1e-41; the estimate is
# then wrong with no sign of trouble (the variance of the standardized mean
# difference came out 2.6% low at 1e52) or sits at its bound, so that
# ordinary data are reported as not yet usable. Because the divisor is a
# power of two the division is exact, so for data on ordinary scales a
# statistic that does not depend on the scale changes at most by rounding
# in the last place (a relative change of at most 5e-14 in 3,000 samples).
#' @noRd
.aipe_scale <- function(...) {
  largest <- max(abs(c(...)))
  if (!is.finite(largest) || largest <= 0) {
    return(1)
  }
  .aipe_power_of_two(largest)
}

# The power of two 2^e with 2^e <= value < 2^(e + 1), for a positive finite
# value. log2() is rounded, so its integer part can be one off for a value
# within rounding of a power of two (log2 of the largest double is returned
# as 1024, and 2^1024 is infinite); the exponent is therefore checked
# against the value itself and moved by one if it is off in either
# direction. `e` is the first guess of the exponent, an argument only so
# that both corrections can be tested.
#' @noRd
.aipe_power_of_two <- function(value, e = floor(log2(value))) {
  if (2^e > value) {
    e <- e - 1
  } else if (2^(e + 1) <= value) {
    e <- e + 1
  }
  2^e
}

# A raw variance estimate that is nothing but cancellation noise, as zero.
#
# Several of the estimators are a sum of moment terms of both signs, and for
# some data that vary the terms cancel exactly on paper (for the Pearson
# correlation and the slope: collinear pairs whose sample kurtosis is
# (n^2 - 3) / (n - 1)^2, which discrete data can produce). What is computed
# then is rounding noise of the order of 1e-16 times the terms, of either
# sign. Under the rule in R/var_floor.R the sign would decide between "used
# through the bound n^-3" (the stopping rule is evaluated, and at a small
# sample met at once) and "not yet usable", so the same data would get
# different answers after a change of units (y = x was met where y = 0.3 x
# was not). A value that is no more than the package wide tolerance
# (.smsd_constant_tol, ten significant digits) times the sum of the absolute
# terms is therefore returned as exactly zero, which the check functions
# report as not yet usable. Every other value is returned untouched, so the
# estimate for data that do not cancel is the same number as without this
# step.
#
# value: the estimate, computed by the caller as the sum of `terms`.
# terms: the signed terms whose sum it is, on the same scale.
#' @noRd
.aipe_zero_if_cancelled <- function(value, terms) {
  if (isTRUE(is.finite(value) &&
             abs(value) <= .smsd_constant_tol * sum(abs(terms)))) {
    return(0)
  }
  value
}

# Is the sample mean positive by more than rounding error? Centered or
# standardized scores have a mean of about 1e-17 of either sign; a positive
# one taken as a real mean would give a coefficient of variation near 1e16.
# The tolerance has the form stats::t.test() uses to call data essentially
# constant, and is the one seq_mr_cv() applies to the same question, so the
# two procedures for the coefficient of variation agree on which samples
# have no estimate.
#' @noRd
.aipe_mean_positive <- function(x) {
  mean(x) > 10 * .Machine$double.eps * max(abs(x))
}

# The note attached when the sample mean is not positive in a procedure
# whose effect size is defined only for a positive mean.
#' @noRd
.note_mean_not_positive <- function(what) {
  paste0("The sample mean is not positive (it is zero, negative, or not ",
         "distinguishable from zero at machine precision), so the ", what,
         " (defined for a positive mean) cannot be estimated from the ",
         "current data and the stopping rule is reported as not met. ",
         "Collect more data and check again.")
}

# The note attached, for an effect size that has units (the slope, the
# contrast), when the raw variance estimate is positive but below the bound
# n^-3, so that the bound was used in its place. The decision is the
# article's (the rule is evaluated with the bound and the interval is built
# on it); the note only tells the user that this happened and how to avoid
# it, because the bound is a pure number and whether an estimate falls below
# it depends on the units of the data. Without it the user would see a
# required sample size in the millions and no reason for it.
#
# effect: the effect size and its units, for example "the slope is in units
#   of y per unit of x".
# cause: what the situation means, in the user's terms.
# remedy: what to rescale.
# topic: the help topic whose Units paragraph explains it.
#' @noRd
.note_bound_applied <- function(effect, cause, remedy, topic) {
  paste0("The variance estimate is positive but below its lower bound, ",
         "n^-3, so the bound was used in its place, as the source article ",
         "prescribes: the interval is wider, and the required sample size ",
         "larger, than the estimate itself would give. The bound is a pure ",
         "number, whereas ", effect, ": ", cause, ". ", remedy,
         " to avoid it (see 'Units' in the help of ", topic, "()).")
}

# Is a raw (unbounded) variance estimate a genuine positive number that is
# nevertheless below the bound n^-3, so that the bound replaces it?
#' @noRd
.aipe_below_bound <- function(raw, n) {
  is.finite(raw) && raw > 0 && raw < n^(-3)
}

# The note attached, by the procedures with groups, when a group has no
# variation while another has: the variance of that group has not been
# estimated, so neither has the variance of the effect size (the rule of
# R/var_floor.R, as seq_budget() applies it to a contrast). `which_groups`
# says which, in words ("group_1", "group 2", "groups 1 and 3").
#' @noRd
.note_group_constant <- function(which_groups, several = FALSE) {
  paste0("The observations in ", which_groups, " are all equal, or differ ",
         "too little relative to their size for a variance to be computed, ",
         "so the variance of ", if (several) "those groups" else "that group",
         " has not been estimated. If a group shows no variation well past ",
         "the pilot sample, reconsider the measure (a ceiling or a floor).")
}

# The target full width of the interval for an effect size confined to a
# range of length `range` (2 for a correlation, 1 for the Gini index): an
# interval as wide as the range says nothing, so such a target is refused,
# as seq_aipe_R2() refuses a width of 1 or more. Without this a width typed
# as a percentage, or ten times too large, gives a pilot of 4 and a rule
# met almost at once.
#' @noRd
.aipe_check_width_bounded <- function(width, range) {
  .check_width(width)
  .check_number(width, "width", lower = 0, upper = range,
                lower_open = TRUE, upper_open = TRUE)
}

# Assemble the smsd_check for one evaluation of the stopping rule
#
#   n >= ceiling((2 z / width)^2 (xi2 + 1 / n))
#
# of Kelley, Darku, and Chattopadhyay (2018, Equation 22; 2019, Equation 10),
# with the interval estimate -/+ z sqrt(xi2 / n) (2018, Equation 23; 2019,
# Equation 11), where xi2 = max{raw, n^-3} is the bounded variance estimate
# of the articles.
#
# estimate: the point estimate; a non-finite value means "cannot be
#   estimated from the current data".
# raw: the RAW variance estimate, as returned by the matching .var_* helper
#   with floor = FALSE (before the bound n^-3 is applied). It may be zero,
#   negative, or not finite.
# n: the sample size the stopping rule counts (per group for the two group
#   and contrast procedures), which is also the n of the bound.
# not_estimable_note: the note to attach when `estimate` is not finite
#   (default: the standard note).
# bound_note: an informational note to attach when the raw estimate is
#   positive but below the bound, so that the bound was used (default:
#   none). Only the slope and the contrast, whose variance has units, pass
#   one. It never changes the decision.
# unusable_note: a note that says why the variance estimate is unusable,
#   attached after the standard note when it is (default: none). The
#   procedures with groups pass one that names a group with no variation.
#
# The rule of R/var_floor.R, as this family applies it:
#   - raw positive: xi2 = max{raw, n^-3} is used as the articles prescribe;
#     the rule is evaluated and the interval is built.
#   - raw zero, negative, or not finite (.var_unusable()), or the point
#     estimate not finite: nothing has been estimated yet. The rule is
#     reported as not met, the confidence limits are NA (no interval is
#     built on an unusable variance), and criterion_n is at least n + 1, so
#     that is_satisfied is identical to current_n >= criterion_n in every
#     case and a loop that collects "criterion_n - current_n" more
#     observations always moves forward. xi_squared then holds the bound,
#     which is what the articles would have used. This is a deliberate
#     departure from the rule as printed, which would stop at once with a
#     near zero width interval.
#
# A criterion that is not finite (a width so small that (2 z / width)^2
# overflows) is an error, as it is in the pilot functions, so that a
# non-finite number never reaches the returned object. A criterion that is
# finite but beyond the range of an integer is reported as it is (a double),
# which is why the ceiling is taken here and not by .as_sample_size().
#' @noRd
.aipe_result <- function(estimate, raw, n, width, conf_level, term, method,
                         extra = NULL, not_estimable_note = NULL,
                         bound_note = NULL, unusable_note = NULL) {
  z <- .aipe_z(conf_level)
  estimable <- is.finite(estimate)
  unusable <- .var_unusable(raw)
  xi2 <- .var_floor(raw, n)
  usable <- estimable && !unusable

  criterion <- ceiling((2 * z / width)^2 * (xi2 + 1 / n))
  if (!is.finite(criterion)) {
    .as_sample_size(criterion, "sample size required by the stopping rule",
                    advice = .aipe_advice_width())
  }
  if (usable) {
    half_width <- z * sqrt(xi2 / n)
  } else {
    half_width <- NA_real_
    criterion <- max(criterion, n + 1)
    if (!estimable) estimate <- NA_real_
  }

  notes <- c(
    if (!estimable) {
      if (is.null(not_estimable_note)) .note_not_estimable() else
        not_estimable_note
    },
    if (unusable) c(.note_variance_unusable(), unusable_note),
    if (usable && .aipe_below_bound(raw, n)) bound_note
  )

  new_smsd_check(
    term = c("lower_limit", term, "upper_limit", "xi_squared",
             "criterion_n", "current_n"),
    value = c(estimate - half_width, estimate, estimate + half_width, xi2,
              criterion, n),
    method = method,
    is_satisfied = usable && n >= criterion,
    current_n = n,
    criterion_n = criterion,
    conf_level = conf_level,
    width = width,
    extra = c(list(variance_unusable = unusable,
                   coverage_guarantee = "asymptotic"), extra),
    notes = notes
  )
}
