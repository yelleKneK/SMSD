# Internal helpers shared by the purely sequential classical procedures
# (seq_fwd_*, seq_blci_*, seq_br_*, and the classical seq_mr_*_mean
# functions). Not exported. Argument checking is delegated to the shared
# checks in R/checks.R.

# The package's pilot rule for the purely sequential classical procedures.
#
# Each of these procedures has an optimal fixed sample size of the form
# n* = K * g(scale): K (`constant`, the criterion constant on the sample size
# scale) is known from the precision target, and the unknown scale of the
# data enters only through g. A pilot equal to n* at unit scale is of the
# same order as the answer, so it oversamples whenever the scale is small in
# the units of the data and offers no protection when the scale is large.
# The rule here is of lower order,
#
#   m = max(m_0, ceiling(K^(1 / (1 + gamma)))),
#
# the form of the pilot in Chattopadhyay and Kelley (2016, Equation 24),
# where K = (A / c)^(1 / 2), and, with gamma = 1, of Kelley, Darku, and
# Chattopadhyay (2018, Equation 21), where K = (2 z / omega)^2. It is the
# package's rule for these procedures; it is not presented as the pilot size
# prescribed by Mukhopadhyay and de Silva (2009).
#
# `min_m_0` is the smallest pilot the procedure can use at all (two
# observations for a sample variance or a scale estimate). `advice` closes
# the message when the pilot is too large to represent: it names the
# arguments of the calling procedure that set the constant.
#' @noRd
.seq_pilot <- function(constant, m_0, gamma, advice, min_m_0 = 2L) {
  m_0 <- .check_count(m_0, "m_0", min = min_m_0)
  gamma <- .check_positive(gamma, "gamma")
  max(m_0, .as_sample_size(constant^(1 / (1 + gamma)), "pilot sample size",
                           advice = advice))
}

# Round a criterion sample size up to an integer of at least 1.
#
# .as_sample_size() refuses a value too large to represent; `advice` closes
# that message and names the arguments of the calling procedure to check. At
# the other extreme, a precision target that is astronomically loose for the
# spread of the data can make the raw criterion underflow to exactly zero
# even though the variance estimate is positive; a criterion of zero
# observations means nothing, so the smallest criterion reported is 1.
#' @noRd
.seq_criterion <- function(value, advice) {
  max(1L, .as_sample_size(value, advice = advice))
}

# Do the data carry no usable information about their spread?
#
# A spread of zero (every observation equal, as can happen early with tied or
# coarsely recorded data) carries no information about the scale. Plugged
# into a stopping rule it gives a criterion sample size of zero, so the rule
# would be met at once. Following the package wide rule for a variance
# estimate that cannot be used (R/var_floor.R), the check functions report
# the stopping rule as not met, the object carries
# attr(., "variance_unusable") = TRUE, the quantities that depend on the
# estimate are reported as NA, and the standard note is attached.
#
# The n^(-3) lower bound of the effect size estimators does not apply here:
# these variances are in the units of the data, so no fixed threshold on the
# variance is meaningful. The data are judged directly instead, with the
# package wide test .is_constant() (R/checks.R): a range of no more than
# 1e-10 times the largest absolute value counts as no spread. Exact equality
# is not enough, because observations that are equal on paper are often not
# equal in floating point: ten gain scores of 0.3, each computed as a
# posttest minus a pretest of up to about 100, differ from one another by as
# much as 4e-15 and have a "variance" near 2e-30, which would pass as a
# usable estimate and satisfy any rule with a criterion of one observation.
# The price is that data whose real spread is that small beside their size
# (measurements near 1e15 that differ by a few units) are also reported as
# not yet usable; subtracting a common value first removes the problem, and
# the note attached in that case (.note_seq_rounding()) says so.
#
# When the data do vary but the estimate is zero or not finite, the scale of
# the data is outside what double precision arithmetic can carry (the
# squares overflow near 1e154 or underflow near 1e-162). More data cannot
# change that, so it is an error that says what to do, not a "not yet"
# result.
#
# `v` is the computed sample variance or scale estimate; `x` is the sample
# it came from; `what` names the estimate in the error messages.
#' @noRd
.seq_zero_spread <- function(v, x, what = "sample variance") {
  if (.is_constant(x)) {
    return(TRUE)
  }
  if (.var_unusable(v)) {
    if (!is.finite(v)) {
      .smsd_stop("The ", what, " of 'x' is too large to represent. ",
                 "Express the data in larger units (divide them by a ",
                 "constant) and state the precision target in the same ",
                 "units.")
    }
    .smsd_stop("The ", what, " of 'x' is too small to represent (it ",
               "underflows to zero although the observations differ). ",
               "Express the data in smaller units (multiply them by a ",
               "constant) and state the precision target in the same units.")
  }
  FALSE
}

# The notes attached to a check whose data show no usable spread: the
# package's standard note, and, when the observations are not exactly equal,
# a second note that says why they were treated as constant and what to do
# if the variation is real.
#' @noRd
.seq_notes_unusable <- function(x) {
  c(.note_variance_unusable(), if (max(x) > min(x)) .note_seq_rounding())
}

#' @noRd
.note_seq_rounding <- function() {
  paste0("The observations are not all equal, but they agree to ten ",
         "significant digits (their range is no more than ",
         format(.smsd_constant_tol), " times their largest absolute ",
         "value), which cannot be told apart from rounding error, so their ",
         "spread is treated as zero. If the variation is real, subtract a ",
         "common value from every observation and add it back to the ",
         "estimate.")
}

# The sample standard deviation for a reported row: sqrt(s2) whenever the
# sample variance s2 is a positive finite number, which is every case in
# which the stopping rule is evaluated. The squares inside a variance
# overflow for deviations beyond about 1e154 and underflow below about
# 1e-162, while the standard deviation itself is still an ordinary number:
# for 1e300 * (1 + c(0, 1, 3, 1) * 1e-13) the variance is Inf and the
# standard deviation is about 1.26e287. Data like these are constant up to
# rounding and are reported as not met, but the row must still hold a
# number, so the standard deviation is then computed on the data divided by
# a power of two near their largest absolute value (an exact operation) and
# scaled back. The exponent is kept within the range of the powers of two
# that are doubles, -1074 to 1023: log2() of the largest double rounds to
# 1024, and log2(0) is -Inf (a sample of zeros is then divided by 2^-1074
# and stays zero). For exactly equal observations the result is exactly
# zero at any magnitude (stats::var() can give Inf for three observations
# of 1.7e308, when their sum overflows).
#' @noRd
.seq_sd <- function(x, s2) {
  if (is.finite(s2) && s2 > 0) {
    return(sqrt(s2))
  }
  unit <- 2^min(1023, max(-1074, floor(log2(max(abs(x))))))
  unit * stats::sd(x / unit)
}

# The check object of this family: new_smsd_check(), which refuses a row
# that is Inf or NaN, behind one more protection. The full width 2 d that
# these procedures report overflows when d is at the very top of the double
# precision range; that is an error, stated in terms of 'd', that says to
# change the units, not a result.
#' @noRd
.seq_check <- function(term, value, ..., width = NULL) {
  if (!is.null(width) && !is.finite(width)) {
    .smsd_stop("The full width of the interval, twice 'd', is too large to ",
               "represent (", format(width), "). Express the data in larger ",
               "units (divide them by a constant) and state 'd' in the same ",
               "units.")
  }
  new_smsd_check(term = term, value = value, ..., width = width)
}

# A required argument with no default. Stops with the package's own message
# when the caller did not supply it, in place of R's 'argument "A" is
# missing, with no default'. Call as .seq_need(missing(A), "A", "...").
#' @noRd
.seq_need <- function(is_missing, name, what) {
  if (is_missing) {
    .smsd_stop("'", name, "' must be supplied (", what, ").")
  }
  invisible(TRUE)
}

# An estimated risk must be a finite number to be reported. Loss and cost
# arguments of absurd size can make it overflow while the criterion sample
# size is still representable; the result would be a 'met' rule beside a
# risk of Inf. `...` are the named arguments that enter the risk, shown in
# the message with their values.
#' @noRd
.seq_finite_risk <- function(risk, ...) {
  if (!is.finite(risk)) {
    args <- list(...)
    shown <- paste0(names(args), " = ", vapply(args, format, character(1)),
                    collapse = ", ")
    .smsd_stop("The estimated risk is too large to represent (", format(risk),
               ") with ", shown, ". Check these arguments and the units of ",
               "'x'.")
  }
  risk
}

# The smallest sample size at which an order statistic interval for the
# median can reach the confidence level. The widest such interval,
# [X_(1), X_(n)], misses the median only when all n observations fall on the
# same side of it, so its coverage is 1 - 2^(1 - n); the minimum is the
# smallest n with 1 - 2^(1 - n) >= conf_level.
#' @noRd
.seq_blci_min_n <- function(conf_level) {
  n <- max(2, ceiling(1 - log2(1 - conf_level)))
  while (1 - 2^(1 - n) < conf_level) n <- n + 1
  while (n > 2 && 1 - 2^(2 - n) >= conf_level) n <- n - 1
  as.integer(n)
}

# Exact coverage of [X_(b), X_(n - b + 1)] for the median of a continuous
# distribution at a fixed sample size n: the interval covers when the number
# of observations below the median, a Binomial(n, 1/2) count, lies between
# b and n - b. For b = 1 this is 1 - 2^(1 - n), written out so that it agrees
# exactly with .seq_blci_min_n().
#' @noRd
.seq_blci_coverage <- function(b, n) {
  if (b <= 1) {
    return(1 - 2^(1 - n))
  }
  1 - 2 * stats::pbinom(b - 1, n, 0.5)
}

# The lower rank b of the order statistic interval [X_(b), X_(n - b + 1)].
#
# The starting value is ceiling((n - z sqrt(n) - 1) / 2), the continuity
# corrected normal approximation to the binomial quantile. Its exact
# coverage is at least conf_level for every n up to 3000 at confidence
# levels of 0.80 and above, but at lower levels it occasionally falls short
# (by 0.003 at a 55% level with n = 7, the largest shortfall found), so b is
# lowered until the exact coverage is at least conf_level. Returns 0 when
# even the widest interval cannot reach conf_level (n below
# .seq_blci_min_n()).
#
# How this rank compares with the exact binomial rule (the largest b whose
# exact coverage is at least conf_level, that is, the narrowest admissible
# interval), for every n from the minimum to 2,000: it is the same rank at
# 1989 of 1996 sample sizes at 90%, 1978 of 1995 at 95%, and 1945 of 1993 at
# 99%, and one rank lower (a slightly wider interval) at the rest. It is
# never higher, so its fixed n coverage is never below conf_level (the
# smallest over that range is 0.9000357, 0.9500466, and 0.9900257, by
# summing binomial probabilities; tests/testthat/test-seq_blci.R repeats the
# comparison).
#' @noRd
.seq_blci_rank <- function(n, conf_level) {
  z <- .aipe_z(conf_level)
  b <- max(1, ceiling((n - z * sqrt(n) - 1) / 2))
  while (b >= 1 && .seq_blci_coverage(b, n) < conf_level) b <- b - 1
  as.integer(b)
}

# The note attached to a bounded length interval check with too few
# observations for the confidence level.
#' @noRd
.note_blci_too_few <- function(n, min_n, conf_level) {
  paste0("With ", n, if (n == 1L) " observation" else " observations",
         " no order statistic interval for the ",
         "median can reach ", format(100 * conf_level), "% confidence (at ",
         "least ", min_n, " are needed), so the stopping rule is reported ",
         "as not met. Collect more data and check again.")
}

# The note attached to a bounded length interval check whose stopping rule is
# met. At a fixed sample size the interval has at least the stated coverage,
# but a sequential procedure stops when the interval happens to be short, and
# short intervals miss the median more often. The figures are those of the
# simulations reported in the help page of seq_blci_dfr_median(): 0.842 to
# 0.879, 0.918 to 0.947, and 0.983 to 0.989 over normal, exponential, and
# Cauchy data and four to six target widths, 40,000 replications each.
#' @noRd
.note_blci_sequential <- function() {
  paste0("The confidence level holds for a sample size fixed in advance. A ",
         "sequential procedure that stops at this check covers the median ",
         "less often than stated: in simulation, about 0.84 to 0.88 for a ",
         "nominal 0.90, 0.92 to 0.94 for 0.95, and 0.98 to 0.99 for 0.99. ",
         "See the Details of ?seq_blci_dfr_median.")
}

# The criterion constant of the bounded risk procedure, (B / w)^(2 / k), so
# that the smallest sample size with risk at most w is this constant times
# sigma^2. Written out for the two powers so that no step is left to pow().
#' @noRd
.seq_br_scale <- function(A, k, w) {
  ratio <- .br_constant(A, k) / w
  if (k == 2) ratio else ratio^2
}
