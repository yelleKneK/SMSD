# Internal helpers for the two stage (ts_*) procedures, beside those in
# R/ts_internals.R. Not exported.

# The rule for a first stage (pilot) sample size.
#
# `constant` is the known-parameter optimal sample size evaluated at a unit
# scale (for example (z / d)^2 for a fixed width interval for a normal mean,
# whose optimal size is (z / d)^2 * sigma^2). Raising it to the power
# 1 / (1 + gamma), gamma > 0, gives a first stage that grows with the size of
# the problem but more slowly than the optimal sample size itself. This is
# the form of pilot that Chattopadhyay and Kelley (2016, Equation 24)
# recommend for their sequential procedure, where they give it as following
# Mukhopadhyay and de Silva (2009); applying it to the two stage procedures
# here is this package's choice. `m_0` is the smallest first stage accepted.
# `advice` closes the message when the result cannot be represented; each
# pilot function supplies the sentence that fits its arguments. Returns an
# integer.
#
# For gamma = 1 (the default everywhere) the power is a square root, and
# sqrt() is used because it is correctly rounded on every platform, so a
# constant that is a perfect square gives exactly its root; `^` carries no
# such guarantee. The rounding up is left to .as_sample_size(), which
# ignores representation error, so that an exact power for any other gamma
# (32^(1 / 1.25) = 16, computed as 16.000000000000004) also comes back as
# itself.
#' @noRd
.ts_pilot_size <- function(constant, m_0, gamma, advice) {
  root <- if (gamma == 1) sqrt(constant) else constant^(1 / (1 + gamma))
  max(m_0, .as_sample_size(root, what = "pilot sample size", advice = advice))
}

# Validate the two pilot rule arguments shared by every ts_*_pilot function.
# `min_m_0` is the smallest first stage for which the procedure is defined.
# Returns list(m_0 = <integer>, gamma = <double>).
#' @noRd
.ts_check_pilot_args <- function(m_0, gamma, min_m_0 = 2L) {
  list(m_0 = .check_count(m_0, "m_0", min = min_m_0),
       gamma = .check_positive(gamma, "gamma"))
}

# The sentences that close the shared messages, worded for the two stage
# procedures.
#
# .ts_advice_first_stage(): for a first stage with too few observations
# (.check_data()). The shared default, "Collect more data, then check the
# stopping rule again", suits a purely sequential check; here the first
# stage fixes the total sample size, so what is needed is a larger first
# stage, and the pilot function gives one. `needs` says what the minimum is
# for.
#
# .ts_advice_width() and .ts_advice_risk(): for a total or pilot sample size
# too large to represent (.as_sample_size()), that is, above
# .Machine$integer.max, a little over two billion. That happens when the
# precision target is many orders of magnitude finer than the data support,
# which in practice is a target stated in the wrong units. `width` is "half"
# or "full", the role of d in the calling procedure.
#
# For the two rules whose quantile depends on the first stage size there is
# a second cause, which `first_stage` (from .ts_first_stage_cause()) carries
# when it applies: see there.
#' @noRd
.ts_advice_first_stage <- function(needs, pilot_fun) {
  paste0(needs, " Enlarge the first stage; ", pilot_fun,
         "() returns a suitable size.")
}

#' @noRd
.ts_advice_width <- function(width, first_stage = NULL) {
  check_d <- paste0("heck that 'd', the ", width, " width of the confidence ",
                    "interval, is in the units of the data")
  if (is.null(first_stage)) {
    return(paste0("C", check_d, "; the target as given calls for more than ",
                  "two billion observations."))
  }
  paste0(first_stage, ", or c", check_d, ".")
}

# Stein's rule and the negative exponential rule pay for estimating the
# spread from the first stage through a quantile (t on m - 1 degrees of
# freedom; F on 2 and 2m - 2) that has no upper limit as the confidence
# level approaches 1 when the first stage is very small: on one degree of
# freedom the t quantile at a confidence level of 0.9999 is 6366, where the
# normal quantile is 3.89. A total sample size too large to represent can
# then be the work of the first stage size and not of d. That is the case
# when the same rule with the constant for a known spread (`known_n`, on the
# sample size scale) is representable, and the message then leads with it:
# the two constants, the sample size the known spread constant would give,
# and the remedy. `used` and `known` are the two clauses that name the
# constants with their values; `m` is the first stage size. Returns NULL
# when the rule is not representable with the known spread constant either,
# so that .ts_advice_width() gives the advice about d alone.
#' @noRd
.ts_first_stage_cause <- function(m, known_n, used, known, pilot_fun) {
  if (!is.finite(known_n) || known_n > .Machine$integer.max) {
    return(NULL)
  }
  paste0("With ", m, " first stage observations ", used, " at this ",
         "confidence level, where ", known, " the rule would call for ",
         max(m, .as_sample_size(known_n)), " observations. Enlarge the ",
         "first stage (", pilot_fun, "() returns a suitable size)")
}

#' @noRd
.ts_advice_risk <- function() {
  paste0("Check that the loss constant 'A' and the risk bound 'w' are in ",
         "the units intended; the bound as given calls for more than two ",
         "billion observations.")
}

# A two stage procedure fixes the total sample size from a measure of spread
# in the first stage alone. When every first stage observation is equal that
# measure is zero, later observations cannot change it, and no second stage
# size can be defined; stop and say what to do. The same holds when the
# observations are equal up to rounding error, judged by the package wide
# rule of .is_constant() (a range of no more than .smsd_constant_tol times
# the largest absolute value): the spread is then numerically zero, and
# plugging it into the rule would report the stopping rule as met at the
# first stage size. Data whose real variation is that small beside their
# size are refused as well, and the message gives the remedy (subtract a
# common value; N does not depend on it).
#
# This is an error and not the "not met" result with variance_unusable = TRUE
# that the purely sequential checks return for the same data. There, more
# data can make the variance estimable and the rule is evaluated again; here
# N depends on the first stage alone, so a rule reported as not met could
# never be met by collecting more.
#
# `what` names the first stage estimate that is zero ("variance", "scale
# estimate"); `unknown` names what cannot be estimated when the observations
# are equal only up to rounding error (the "first stage variance"; the
# "scale", because "the scale estimate cannot be estimated" is not a
# sentence).
#' @noRd
.ts_stop_if_constant <- function(pilot, what,
                                 unknown = paste("first stage", what)) {
  if (!.is_constant(pilot)) {
    return(invisible(TRUE))
  }
  spread <- max(pilot) - min(pilot)
  advice <- paste0("Enlarge the first stage until it shows variation ",
                   "(recording the outcome more finely if it is heavily ",
                   "rounded), and pass every first stage observation as ",
                   "'pilot_data'.")
  if (spread == 0) {
    .smsd_stop("All ", length(pilot), " first stage observations in ",
               "'pilot_data' are equal, so the first stage ", what,
               " is zero and the total sample size cannot be determined ",
               "from it. ", advice)
  }
  .smsd_stop("The ", length(pilot), " first stage observations in ",
             "'pilot_data' are equal up to rounding error (their range, ",
             format(spread), ", is no more than ",
             format(.smsd_constant_tol), " times their largest absolute ",
             "value, ", format(max(abs(pilot))), "), so the ", unknown,
             " cannot be estimated from them and the total sample size ",
             "cannot be determined. ", advice, " If the variation is ",
             "real, subtract a common value from every observation in both ",
             "stages and add it back to the estimate; the total sample ",
             "size does not depend on it.")
}

# The bounded risk rules square the observations: a variance, or the square
# of a mean. When that square is not a usable number although it is positive
# in exact arithmetic, the scale of the data is outside what double
# precision arithmetic can carry (squares overflow beyond about 1e154 and
# lose their digits below about 1e-154, where they fall under
# .Machine$double.xmin). Plugged into the rule, an overflow gives a message
# about the loss constant and the risk bound although the data are the
# cause, and an underflow gives a total sample size equal to the first
# stage and a risk of zero. More data cannot change either, so it is an
# error that says what to do, in the words the purely sequential checks use
# for the same case (.seq_zero_spread() in R/seq_internals.R). `x` is the
# computed square; `what` names it; `positive` says why it cannot be zero.
# Returns `x`.
#' @noRd
.ts_usable_square <- function(x, what, positive) {
  restate <- " and restate the loss constant 'A' for those units."
  if (!is.finite(x)) {
    .smsd_stop("The ", what, " is too large to represent. Express the data ",
               "in larger units (divide them by a constant)", restate)
  }
  if (x < .Machine$double.xmin) {
    .smsd_stop("The ", what, " is too small to represent (it underflows ",
               "although ", positive, "). Express the data in smaller ",
               "units (multiply them by a constant)", restate)
  }
  x
}

# Every row of a returned check is a finite number. Data or targets of a
# magnitude that double precision cannot carry would otherwise give a rule
# reported as met beside an estimate, a confidence limit, or a risk of Inf
# or NaN (the purely sequential checks stop in the same case, see
# .seq_finite_risk() in R/seq_internals.R). `term` and `value` are the rows
# about to be returned; `restate` closes the message with what to do about
# the targets of the calling procedure.
#' @noRd
.ts_finite_rows <- function(term, value, restate) {
  bad <- which(!is.finite(value))
  if (length(bad)) {
    .smsd_stop("The '", term[bad[1L]], "' row of the result cannot be ",
               "represented (it is computed as ", format(value[bad[1L]]),
               "): the data or the targets are too large in magnitude. ",
               "Express the data in larger units (divide them by a ",
               "constant) and ", restate)
  }
  invisible(TRUE)
}

# The stopping constant of the two stage bounded risk procedure for a normal
# mean: b = E[(S_m^2 / sigma^2)^(-k/2)]
#         = ((m - 1) / 2)^(k/2) * Gamma((m - k - 1) / 2) / Gamma((m - 1) / 2),
# with S_m^2 the variance of m normal observations. The expectation exists
# only for m >= k + 2. For k = 2 it reduces to (m - 1) / (m - 3), which is
# used as it stands. For k = 1 it is computed on the log scale so that it
# neither overflows (gamma() does beyond m = 344) nor loses its sign.
#' @noRd
.ts_br_norm_b <- function(m, k) {
  if (k == 2) {
    return((m - 1) / (m - 3))
  }
  exp((k / 2) * log((m - 1) / 2) + lgamma((m - k - 1) / 2) -
        lgamma((m - 1) / 2))
}
