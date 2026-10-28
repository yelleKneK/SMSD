# Sequential accuracy in parameter estimation for the squared multiple
# correlation coefficient: Kelley, Darku, and Chattopadhyay (2019), section
# "An Extension: Squared Multiple Correlation Coefficient" (pp. 503--507) and
# Appendix C (p. 513). The pilot, the check, the Monte Carlo sensitivity
# wrapper, and their internal helpers are together in this file.
#
# Notation: the article writes k for the number of predictors; this file and
# the user facing arguments write p, as the R2 functions of the sibling
# package DMAR do.

# ---------------------------------------------------------------------------
# Internal helpers (not exported)
# ---------------------------------------------------------------------------

# The target full width for an interval on a squared multiple correlation. A
# parameter confined to [0, 1] cannot have an informative interval of width 1
# or more, so such a target is refused rather than silently accepted.
#' @noRd
.R2_check_width <- function(width) {
  .check_width(width)
  .check_number(width, "width", lower = 0, upper = 1,
                lower_open = TRUE, upper_open = TRUE)
}

# The number of predictors. The smallest usable sample size is p + 2, so p
# may not be so large that p + 2 cannot be represented as an integer.
#' @noRd
.R2_check_p <- function(p) {
  if (missing(p)) {
    .smsd_stop("'p' must be supplied: the number of predictors in the ",
               "regression.")
  }
  p <- .check_count(p, "p", min = 1L)
  if (p > .Machine$integer.max - 2L) {
    .smsd_stop("'p' must be a single whole number of at most ",
               .Machine$integer.max - 2L, " (the procedure needs p + 2 ",
               "observations); got ", p, ".")
  }
  p
}

# The two confidence intervals of the article. The argument is checked with
# the shared .check_choice(), so the untouched default selects the first and
# a single string may abbreviate a choice.
#' @noRd
.R2_ci_methods <- c("wald", "log_transformed")

# What to do when a width is so narrow that the sample size it implies cannot
# be represented: the closing sentence of the .as_sample_size() message.
#' @noRd
.R2_width_advice <- paste0("Check 'width': a target this narrow calls for ",
                           "more observations than a study can collect.")

# What is said, in an error or a note, about a variable whose values differ
# only beyond the tenth significant digit. The package treats such a variable
# as constant (.is_constant(), one tolerance for the whole package). That is
# right for the rounding noise of equal values; a large constant plus a small
# spread needs a word of explanation, because "collect more data" is not the
# remedy for it.
#' @noRd
.R2_ten_digits <- function(subject) {
  paste0(subject, " only beyond the tenth significant digit, which counts ",
         "as not varying. If it is a large constant plus a small spread, ",
         "subtract the constant before calling; that does not change the ",
         "squared multiple correlation.")
}

# The predictors, as a numeric matrix with one row per observation and one
# column per predictor. Accepts a numeric matrix, a data frame whose columns
# are all numeric, or a numeric vector (a single predictor).
#' @noRd
.R2_check_predictors <- function(x) {
  x <- .drop_one_dim(x)
  if (is.data.frame(x)) {
    if (ncol(x) < 1L) {
      .smsd_stop("'x' must have at least one predictor column; got ",
                 .smsd_describe(x), ".")
    }
    not_numeric <- !vapply(x, is.numeric, logical(1L))
    if (any(not_numeric)) {
      .smsd_stop("Every column of 'x' must be numeric; not numeric: ",
                 paste(names(x)[not_numeric], collapse = ", "),
                 ". Code a categorical predictor as numeric indicator ",
                 "columns before calling.")
    }
    x <- as.matrix(x)
  } else if (is.matrix(x)) {
    if (!is.numeric(x)) {
      .smsd_stop("'x' must be a numeric matrix or data frame of predictors; ",
                 "got a ", typeof(x), " matrix.")
    }
  } else if (is.numeric(x) && is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  } else {
    .smsd_stop("'x' must be a numeric matrix or data frame of predictors ",
               "(one row per observation, one column per predictor), or a ",
               "numeric vector for a single predictor; got ",
               .smsd_describe(x), ".")
  }
  storage.mode(x) <- "double"
  if (ncol(x) < 1L) {
    .smsd_stop("'x' must have at least one predictor column; got ",
               .smsd_describe(x), ".")
  }
  if (anyNA(x)) {
    .smsd_stop("'x' contains missing values (NA or NaN). Remove the ",
               "incomplete rows from 'x' and 'y', or impute, before calling.")
  }
  if (any(!is.finite(x))) {
    .smsd_stop("'x' contains infinite values. Remove them before calling.")
  }
  x
}

# Names of predictor columns for messages: the column names when there are
# any, otherwise the column positions. A matrix can have names for some of
# its columns only (cbind(1, x) with a named x, for example); a column
# without a name is then given by its position.
#' @noRd
.R2_column_labels <- function(x, which) {
  which <- as.integer(which)
  nms <- colnames(x)
  if (is.null(nms)) {
    return(paste("column", paste(which, collapse = ", ")))
  }
  labels <- nms[which]
  unnamed <- is.na(labels) | !nzchar(labels)
  labels[unnamed] <- paste("column", which[unnamed])
  paste(labels, collapse = ", ")
}

# Validate the predictors and the response together and return the sample R2.
#
# Returns a list with n, p, R2, one_minus_R2 (computed directly as the
# residual share of the sum of squares, so it stays accurate when R2 is close
# to one), resolution, and response_differs. R2 is NaN when the response does
# not vary. The fit is by a QR decomposition of the centered predictors, and
# R2 is the model sum of squares over the model plus residual sums of
# squares, the quantity that summary(stats::lm()) reports as r.squared.
# Every variable is first divided by a power of two (.R2_binary_unit()), so
# that the result does not depend on the units of the data.
#
# "Does not vary" is the package's one definition, .is_constant(): a range of
# at most 1e-10 times the largest magnitude, that is, values that agree to
# ten significant digits. response_differs says whether a response that does
# not vary in that sense nevertheless has unequal values (a large constant
# plus a small spread, or rounding noise), so that the check function can
# say so. It is FALSE for a response that varies.
#
# resolution is the smallest R2, or 1 - R2, that these data can tell apart
# from zero. A stored value is rounded to within eps * max|v| / 2, which is
# eps * kappa / 2 standard deviations of the variable when
# kappa = max|v| / sd(v). A response that is an exact linear function of the
# predictors therefore has a residual share of up to the order of
# (eps * kappa)^2, possibly exactly zero, when a variable is a large
# constant plus a small spread. In 4,518 simulated data sets of that kind
# that pass the test for variation, the residual share was at most
# 24 (eps * kappa)^2. (n of p + 2, 12, 40, or 200; p of 1, 2, 5, or 10
# standard normal predictors, each multiplied by 10^u with u uniform on
# (-2, 2) and shifted by an offset of 1e3 to 1e12 of either sign; a
# response that is the predictors times coefficients uniform on (-3, 3)
# plus 0, the offset, minus the offset, or 100 times it.)
# resolution is (100 eps * kappa)^2 with kappa the largest over the
# predictors and the response, and never less than eps, the value for data
# without such an offset. Because a variable must vary in its first ten
# significant digits to be accepted at all, kappa stays near 1e10 or below
# (5.6e10 at most in those data sets) and resolution near 1e-6 or below, so
# a response with real error is not taken for a perfect fit: in data sets
# drawn the same way with a residual standard deviation of 1% or 10% of
# that of the fitted values, none was.
#
# A predictor that does not vary, or predictors that are collinear, stop
# with an error that says what to do. The sensitivity engine instead passes
# stop_on_degenerate = FALSE and receives an R2 of NaN: in a simulated study
# such a sample is a check at which the rule is not met, and sampling goes
# on (the recourse the error message gives the user).
#' @noRd
.R2_estimate <- function(x, y, stop_on_degenerate = TRUE) {
  if (missing(x) || missing(y)) {
    .smsd_stop("Both 'x' (the predictors) and 'y' (the response) must be ",
               "supplied.")
  }
  x <- .R2_check_predictors(x)
  y <- .check_data(y, "y", min_n = 0L, na.rm = NULL)
  n <- nrow(x)
  p <- ncol(x)
  if (length(y) != n) {
    .smsd_stop("'x' and 'y' must describe the same observations: 'x' has ",
               n, " row", if (n == 1L) "" else "s", " and 'y' has ",
               length(y), " value", if (length(y) == 1L) "" else "s", ".")
  }
  if (n < p + 2L) {
    .smsd_stop("With ", p, " predictor", if (p == 1L) "" else "s",
               " this procedure needs at least ", p + 2L,
               " observations (with p + 1 or fewer the sample squared ",
               "multiple correlation is identically 1); there are ", n,
               ". Collect more data, then check the stopping rule again.")
  }

  not_estimable <- list(n = n, p = p, R2 = NaN, one_minus_R2 = NaN,
                        resolution = .Machine$double.eps,
                        response_differs = FALSE)
  constant <- apply(x, 2L, .is_constant)
  if (any(constant)) {
    if (!stop_on_degenerate) return(not_estimable)
    differs <- apply(x[, constant, drop = FALSE], 2L,
                     function(v) min(v) != max(v))
    .smsd_stop("A predictor does not vary in the current data (",
               .R2_column_labels(x, which(constant)),
               " of 'x'), so the squared multiple correlation for ", p,
               " predictor", if (p == 1L) "" else "s", " is not defined. ",
               "Remove the column, or collect more data if it is expected ",
               "to vary.",
               if (any(apply(x[, constant, drop = FALSE], 2L,
                             function(v) all(v == 1)))) {
                 paste0(" (The function fits the intercept itself, so a ",
                        "column of ones is not needed.)")
               },
               if (any(differs)) {
                 paste0(" ", .R2_ten_digits(paste0(
                   "The values of ",
                   .R2_column_labels(x, which(constant)[differs]),
                   " differ")))
               })
  }
  x <- sweep(x, 2L, apply(x, 2L, .R2_binary_unit), `/`)
  x_centered <- sweep(x, 2L, colMeans(x))
  x_qr <- qr(x_centered)
  if (x_qr$rank < p) {
    if (!stop_on_degenerate) return(not_estimable)
    .smsd_stop("The predictors in 'x' are collinear: ", p,
               " columns were supplied but only ", x_qr$rank,
               " are linearly independent in the current data. Remove the ",
               "redundant predictor",
               if (p - x_qr$rank == 1L) "" else "s",
               " (for example an indicator column for every category), or ",
               "collect more data if the dependence is an accident of a ",
               "small sample.")
  }

  if (.is_constant(y)) {
    not_estimable$response_differs <- min(y) != max(y)
    return(not_estimable)
  }
  y <- y / .R2_binary_unit(y)
  y_centered <- y - mean(y)
  ss_model <- sum(qr.fitted(x_qr, y_centered)^2)
  ss_resid <- sum(qr.resid(x_qr, y_centered)^2)
  ss_total <- ss_model + ss_resid

  # kappa = max|v| / sd(v) for every variable, from values divided by their
  # largest magnitude.
  largest <- c(apply(abs(x), 2L, max), max(abs(y)))
  relative <- sweep(cbind(x_centered, y_centered), 2L, largest, `/`)
  kappa <- max(1 / sqrt(colSums(relative^2) / (n - 1L)))
  resolution <- max(.Machine$double.eps,
                    (100 * .Machine$double.eps * kappa)^2)

  # response_differs is about a response that does not vary; this one does.
  list(n = n, p = p, R2 = ss_model / ss_total,
       one_minus_R2 = ss_resid / ss_total, resolution = resolution,
       response_differs = FALSE)
}

# The power of two that brings the largest magnitude of a variable into
# [1, 2). R2 does not depend on the units of a variable, and its computation
# must not either: before anything is squared or summed, .R2_estimate()
# divides every variable by this number. Sums of squares of a response
# recorded in units that make its values of the order of 1e-160 or 1e155
# would otherwise underflow or overflow (and lose digits on the way there),
# and the mean of a predictor near the largest representable number would
# overflow. Dividing by a power of two is exact in binary floating point, so
# for data on any ordinary scale every result is the same to the last bit as
# without it. The variable must not be all zero (the caller has by then
# established that it varies).
#' @noRd
.R2_binary_unit <- function(v) {
  2^floor(log2(max(abs(v))))
}

# The stopping rule and the confidence limits, vectorized over R2 and n so
# that the check function (one value) and the sensitivity engine (many sample
# sizes at once) run exactly the same arithmetic.
#
# Equation numbers refer to Kelley, Darku, and Chattopadhyay (2019).
#
# "wald": interval R2 +/- z * 2 R (1 - R2) / sqrt(n) (Equation 59, with the
#   Olkin and Finn, 1995, asymptotic variance 4 P2 (1 - P2)^2 evaluated at
#   R2). Stopping rule: n >= (16 z^2 / omega^2) (R2 + 1/n) (1 - R2 + 1/n)^2.
#   Equation 63 as printed has (1 - (R2 + 1/n))^2 for the last factor; the
#   form used here, with each factor increased by 1/n, is the one that
#   reproduces the article's Tables 13 and 14 and that keeps the interval at
#   the stopping point no wider than omega. See the help page of
#   seq_aipe_R2().
# "log_transformed": interval 1 - exp(ln(1 - R2) +/- z * 2 R / sqrt(n - p - 2))
#   (Equation 58; Bonett and Wright, 2011). Stopping rule (Equation 62):
#   n >= p + 2 + 4 R2 z^2 / L^2 with
#   L = ln[(1/2)(a + sqrt(a^2 + 4))] = asinh(a / 2), a = omega / (1 - R2).
#   The interval is not defined at n = p + 2.
#
# xi2 is the estimated asymptotic variance of R2, 4 R2 (1 - R2)^2. It is
# never negative. It is unusable, in the sense the whole package gives the
# word (zero, negative, or not finite), when R2 cannot be computed or when
# the sample R2 is 0 or 1 to the precision of the data, which makes the
# estimate zero or rounding noise: the printed rules would then be met at
# once with an interval of zero width. `resolution` is that precision: R2 or
# 1 - R2 at or below it cannot be told apart from zero. It comes from
# .R2_estimate() and is .Machine$double.eps unless a variable is a large
# constant plus a small spread. Exactly orthogonal or exactly fitted data
# give, through the QR decomposition of .R2_estimate(), an R2 within about
# 1e-28 of 0 or 1, far inside it; continuous data essentially never come
# that close. (The running sums of the sensitivity engine resolve 1 - R2
# only to about 1e-15, so the engine recomputes such a sample through
# .R2_estimate() before it applies this rule; see .R2_sensitivity_one().)
# A positive estimate, however small, is used as it stands: the lower bound
# of n^(-3) that the distribution free estimators carry is not part of this
# procedure (it would change the published operating characteristics of the
# log transformed interval at small n).
#
# When no decision can be made (the variance estimate is unusable, or the
# log transformed interval is not yet defined) the limits are NA, the rule
# is not met, and the criterion is at least n + 1, as in the other seq_aipe_*
# check functions: at least one more observation is needed, and
# `satisfied` is TRUE exactly when n reaches the criterion as
# .as_sample_size() rounds it (see .R2_reached()).
#' @noRd
.R2_rule <- function(R2, one_minus_R2, n, p, width, z, ci_method,
                     resolution = .Machine$double.eps) {
  xi2 <- 4 * R2 * one_minus_R2^2
  estimable <- is.finite(R2)
  unusable <- !estimable | R2 <= resolution | one_minus_R2 <= resolution
  if (ci_method == "wald") {
    criterion <- (4 * z / width)^2 * (R2 + 1 / n) * (one_minus_R2 + 1 / n)^2
    half_width <- z * sqrt(xi2 / n)
    lower <- R2 - half_width
    upper <- R2 + half_width
    defined <- !unusable
  } else {
    criterion <- p + 2 + 4 * R2 * z^2 / asinh(width / (2 * one_minus_R2))^2
    defined <- !unusable & n > p + 2
    se_log <- 2 * sqrt(R2) / sqrt(pmax(n - p - 2, 0))
    lower <- 1 - one_minus_R2 * exp(z * se_log)
    upper <- 1 - one_minus_R2 * exp(-z * se_log)
  }
  lower[!defined] <- NA_real_
  upper[!defined] <- NA_real_
  criterion[!defined] <- pmax(criterion, n + 1, na.rm = TRUE)[!defined]
  satisfied <- defined & .R2_reached(n, criterion)
  satisfied[is.na(satisfied)] <- FALSE
  list(xi2 = xi2, criterion = criterion, lower = lower, upper = upper,
       estimable = estimable, unusable = unusable, defined = defined,
       satisfied = satisfied)
}

# Has the sample size reached a criterion? TRUE exactly when n is at least
# the criterion as .as_sample_size() rounds it up, so that the decision and
# the criterion sample size reported beside it cannot disagree: a criterion
# that is a whole number in exact arithmetic and a few units in the last
# place above it in floating point is that whole number. Vectorized, for the
# sensitivity engine. (A test pins the agreement with .as_sample_size().)
#' @noRd
.R2_reached <- function(n, criterion) {
  n >= criterion * (1 - 1e-10)
}

# Is the sample small relative to the number of predictors? The sample R2
# overestimates the population value by roughly (1 - P2) p / n, about
# p / (2 sqrt(P2) sqrt(n)) standard errors, so the interval built around it
# covers less often than the confidence level until n is large relative to
# p (for a few predictors) and to p^2 (for many). The threshold is
# n < p * max(20, 2 p): fewer than 20 observations per predictor, or fewer
# than 2 p per predictor with more than 10 predictors. It is where the
# coverage of the Wald interval at a fixed sample size and a nominal .95
# reaches about .90 for population values of .30 or more (the test "at the
# threshold the fixed sample coverage of the Wald interval is about .90" in
# tests/testthat/test-seq_aipe_R2.R draws nine such coverages from the
# distribution of the sample R2); below it coverage is lower, and lower
# still for smaller population values. Both numbers are doubles, so that a
# large p cannot overflow.
#' @noRd
.R2_small_sample_per_predictor <- function(p) {
  max(20, 2 * as.numeric(p))
}

#' @noRd
.R2_small_sample_below <- function(p) {
  as.numeric(p) * .R2_small_sample_per_predictor(p)
}

# A whole number as the prose of a note writes it: in full, with a comma
# between thousands, whatever the options of the session say. Left to the
# options, a round number of 100,000 or more (and any number under a
# negative `scipen`) is written in scientific notation, and with
# options(OutDec = ",") format() warns that the thousands mark and the
# decimal mark are the same character. A whole number has no decimal mark,
# so the mark is fixed here and nothing depends on the option.
#' @noRd
.R2_whole <- function(x) {
  format(as.numeric(x), big.mark = ",", decimal.mark = ".",
         scientific = FALSE, trim = TRUE)
}

# The note attached to a check whose sample is small relative to the number
# of predictors. It informs; it does not change the decision.
#' @noRd
.note_R2_small_sample <- function(n, p) {
  paste0("With ", .R2_whole(n), " observations and ", .R2_whole(p),
         " predictor", if (p == 1L) "" else "s",
         " the sample is small relative to the number of predictors ",
         "(fewer than ", .R2_whole(.R2_small_sample_below(p)),
         " observations, ", .R2_whole(.R2_small_sample_per_predictor(p)),
         " per predictor). The sample ",
         "squared multiple correlation then overestimates the population ",
         "value, and in simulations the interval of this procedure ",
         "contained the population value less often than the confidence ",
         "level, far less often the smaller the sample. This note does not ",
         "change the stopping rule. See 'Small Samples Relative to the ",
         "Number of Predictors' in the help of seq_aipe_R2().")
}

# The article's pilot sample size, raised to m_0. For the Wald interval the
# article gives max{k, 4 z / omega} (p. 505 and Appendix C, Equation 97); for
# the log transformed interval k + 2 (Appendix C, Equation 95). m_0 is at
# least p + 2, the smallest sample size at which R2 is not identically one.
# The rounding up is the shared .as_sample_size(), as for every sample size
# the package returns.
#' @noRd
.R2_pilot <- function(width, p, z, ci_method, m_0) {
  article <- if (ci_method == "wald") {
    max(p, 4 * z / width)
  } else {
    p + 2
  }
  .as_sample_size(max(m_0, article), "pilot sample size",
                  advice = .R2_width_advice)
}

# Sample R2 at many sample sizes at once, from running sums.
#
# sums: one row per sample size, the column sums of the p predictors and then
#   the response. cross: the matching sums of products for every pair of
#   columns (l <= m), in the column order of `position` (a (p + 1) square
#   matrix whose upper triangle holds the column of `cross` for each pair).
# The centered sums of squares and products are reduced by Gaussian
# elimination on the predictors, done as arithmetic on whole vectors; what is
# left for the response is its residual sum of squares.
#
# Running sums are fast but less exact than the QR decomposition of
# .R2_estimate(). `reliable` is FALSE for a sample size at which they cannot
# be trusted to give the same answer: a variable whose centered sum of
# squares is less than 1e-4 of its raw sum of squares (four or more digits
# lost to cancellation: a mean far from zero relative to the spread, or a
# variable that does not vary), or an R2 that is not finite or is within
# 1e-8 of 0 or 1. The caller recomputes those through .R2_estimate().
#' @noRd
.R2_from_sums <- function(sums, cross, n, position) {
  d <- ncol(sums)
  centered <- vector("list", ncol(cross))
  for (l in seq_len(d)) {
    for (m in l:d) {
      centered[[position[l, m]]] <-
        cross[, position[l, m]] - sums[, l] * sums[, m] / n
    }
  }
  reliable <- rep(TRUE, length(n))
  for (l in seq_len(d)) {
    reliable <- reliable &
      centered[[position[l, l]]] > 1e-4 * cross[, position[l, l]]
  }
  ss_total <- centered[[position[d, d]]]
  for (j in seq_len(d - 1L)) {
    pivot <- centered[[position[j, j]]]
    for (l in (j + 1L):d) {
      multiplier <- centered[[position[j, l]]] / pivot
      for (m in l:d) {
        centered[[position[l, m]]] <-
          centered[[position[l, m]]] - multiplier * centered[[position[j, m]]]
      }
    }
  }
  one_minus_R2 <- pmin(pmax(centered[[position[d, d]]] / ss_total, 0), 1)
  reliable <- reliable & is.finite(one_minus_R2) &
    one_minus_R2 > 1e-8 & one_minus_R2 < 1 - 1e-8
  reliable[is.na(reliable)] <- FALSE
  list(R2 = 1 - one_minus_R2, one_minus_R2 = one_minus_R2,
       reliable = reliable)
}

# One replication of the sequential procedure, run to its stopping point.
#
# draw(m) returns the next m observations as an m by (p + 1) numeric matrix:
# the p predictors, then the response. The stopping rule is checked at the
# pilot sample size and then after every `batch` further observations (and at
# max_n, where sampling stops whatever the rule says). Observations are drawn
# several at a time and the rule is evaluated at every checkpoint among them
# through .R2_from_sums() and .R2_rule(). Wherever the running sums are not
# reliable (see .R2_from_sums(); this does not happen with continuous data
# centered near zero, such as seq_aipe_R2_sensitivity() generates) the
# estimate is recomputed from the observations themselves through
# .R2_estimate(), exactly as seq_aipe_R2() computes it. The result is
# therefore the one a loop of seq_aipe_R2() calls on the same observations
# gives, degenerate samples included; observations drawn beyond the stopping
# point are discarded. A sample on which seq_aipe_R2() would stop with an
# error (a predictor that does not vary, collinear predictors) is a check at
# which the rule is not met, and sampling goes on.
#
# Returns c(n, lower, upper, capped, unusable): the realized sample size, the
# limits there (NA when the cap stops a replication at a check that has no
# interval), whether the cap rather than the rule stopped the replication,
# and whether the variance estimate was unusable at any check up to the
# stopping point.
#' @noRd
.R2_sensitivity_one <- function(draw, pilot_n, p, width, z, ci_method,
                                max_n, batch) {
  d <- p + 1L
  position <- matrix(0L, d, d)
  position[upper.tri(position, diag = TRUE)] <- seq_len(d * (d + 1L) / 2L)
  pairs <- which(upper.tri(position, diag = TRUE), arr.ind = TRUE)
  # Rows handled at a time, so that the matrix of running sums of products
  # stays near 2e6 cells whatever the number of predictors and the pilot
  # sample size.
  chunk_cap <- max(1, floor(2e6 / nrow(pairs)))

  sums <- numeric(d)
  cross <- numeric(nrow(pairs))
  kept <- list()
  n_have <- 0
  n_to <- min(pilot_n, chunk_cap)
  unusable_any <- FALSE
  last_criterion <- NA_real_

  repeat {
    obs <- draw(n_to - n_have)
    kept[[length(kept) + 1L]] <- obs
    rows <- seq.int(n_have + 1, n_to)
    sums_at <- matrix(0, length(rows), d)
    for (l in seq_len(d)) {
      sums_at[, l] <- sums[l] + cumsum(obs[, l])
    }
    cross_at <- matrix(0, length(rows), nrow(pairs))
    for (i in seq_len(nrow(pairs))) {
      cross_at[, i] <- cross[i] +
        cumsum(obs[, pairs[i, 1L]] * obs[, pairs[i, 2L]])
    }
    sums <- sums_at[length(rows), ]
    cross <- cross_at[length(rows), ]
    n_have <- n_to

    # Still collecting the pilot sample: nothing to check yet.
    if (n_have < pilot_n) {
      n_to <- min(pilot_n, n_have + chunk_cap)
      next
    }

    check <- rows >= pilot_n &
      ((rows - pilot_n) %% batch == 0 | rows == max_n)
    if (any(check)) {
      n_check <- rows[check]
      est <- .R2_from_sums(sums_at[check, , drop = FALSE],
                           cross_at[check, , drop = FALSE], n_check,
                           position)
      resolution <- rep(.Machine$double.eps, length(n_check))
      redo <- which(!est$reliable)
      if (length(redo) > 0L) {
        if (length(kept) > 1L) kept <- list(do.call(rbind, kept))
        for (i in redo) {
          upto <- seq_len(n_check[i])
          exact <- .R2_estimate(kept[[1L]][upto, seq_len(p), drop = FALSE],
                                kept[[1L]][upto, d],
                                stop_on_degenerate = FALSE)
          est$R2[i] <- exact$R2
          est$one_minus_R2[i] <- exact$one_minus_R2
          resolution[i] <- exact$resolution
        }
      }
      rule <- .R2_rule(est$R2, est$one_minus_R2, n_check, p, width, z,
                       ci_method, resolution)
      stop_at <- which(rule$satisfied)[1L]
      capped <- FALSE
      if (is.na(stop_at) && n_have >= max_n) {
        stop_at <- length(n_check)
        capped <- TRUE
      }
      if (!is.na(stop_at)) {
        unusable_any <- unusable_any || any(rule$unusable[seq_len(stop_at)])
        return(c(n = n_check[stop_at], lower = rule$lower[stop_at],
                 upper = rule$upper[stop_at], capped = as.numeric(capped),
                 unusable = as.numeric(unusable_any)))
      }
      unusable_any <- unusable_any || any(rule$unusable)
      last_criterion <- rule$criterion[length(n_check)]
    }

    # Draw far enough ahead to reach the current criterion, plus a margin.
    ahead <- if (is.finite(last_criterion)) {
      ceiling(1.1 * last_criterion) - n_have
    } else {
      0
    }
    step <- min(max(ahead, batch, 50), chunk_cap, max_n - n_have)
    n_to <- n_have + step
  }
}

# The term/value summary of a set of replications, as an smsd_sensitivity
# object.
#
# runs: one column per replication, with the rows that .R2_sensitivity_one()
# returns (n, lower, upper, capped, unusable). The rows of the result and
# their order are those of the other seq_aipe_*_sensitivity functions.
#
# A replication stopped by the cap can end without an interval (limits of
# NA: no estimate, an unusable variance estimate, or the log transformed
# interval at n = p + 2). It then has no width and does not cover; the width
# summaries are over the replications that have an interval. This cannot
# happen with the data that seq_aipe_R2_sensitivity() generates, but the
# summaries are formed so that it could not turn them into NA. Proportions
# and percentages are formed from counts, so that, for example, 950 covering
# intervals in 1000 give exactly 0.95.
#' @noRd
.R2_sensitivity_summary <- function(runs, true_value, width, conf_level) {
  G <- ncol(runs)
  realized_n <- runs["n", ]
  realized_width <- runs["upper", ] - runs["lower", ]
  has_interval <- is.finite(realized_width)
  covered <- has_interval & runs["lower", ] <= true_value &
    true_value <= runs["upper", ]
  usable_width <- realized_width[has_interval]

  new_smsd_sensitivity(
    term = c("true_value", "target_width", "conf_level", "replications",
             "mean_n", "median_n", "sd_n",
             "mean_width", "median_width", "pct_width_within", "coverage",
             "pct_hit_max_n", "pct_variance_unusable"),
    value = c(true_value, width, conf_level, G,
              mean(realized_n), stats::median(realized_n),
              stats::sd(realized_n),
              if (any(has_interval)) mean(usable_width) else NA_real_,
              if (any(has_interval)) stats::median(usable_width) else NA_real_,
              100 * sum(has_interval & realized_width <= width) / G,
              sum(covered) / G,
              100 * sum(runs["capped", ]) / G,
              100 * sum(runs["unusable", ]) / G)
  )
}

# ---------------------------------------------------------------------------
# Exported functions
# ---------------------------------------------------------------------------

#' Pilot Sample Size for Sequential AIPE: Squared Multiple Correlation Coefficient
#'
#' Returns the number of observations to collect before the sequential
#' accuracy in parameter estimation (AIPE) stopping rule for the population
#' squared multiple correlation coefficient is evaluated for the first time
#' with \code{\link{seq_aipe_R2}}.
#'
#' Kelley, Darku, and Chattopadhyay (2019, p. 505, and Appendix C, Equations
#' 95 and 97) give a pilot sample size for each of the two confidence
#' intervals they study. With \eqn{p} predictors (\eqn{k} in the article),
#' target full width \eqn{\omega}, and standard normal quantile
#' \eqn{z_{1-\alpha/2}}:
#' \itemize{
#'   \item \code{ci_method = "wald"}:
#'     \eqn{\max\{p, 4 z_{1-\alpha/2}/\omega\}}, rounded up to a whole
#'     number;
#'   \item \code{ci_method = "log_transformed"}: \eqn{p + 2}.
#' }
#' The value returned is the larger of that size and \code{m_0}. The article
#' prints \eqn{4 z_{1-\alpha/2}/\omega} without saying how it is made a whole
#' number; rounding up is what reproduces the mean final sample sizes of its
#' Tables 13 and 14 in the conditions where most replications stop at the
#' pilot sample size.
#'
#' With \eqn{p + 1} or fewer observations the sample squared multiple
#' correlation is identically 1, so \code{m_0} is at least \eqn{p + 2} and
#' the returned pilot sample size is never smaller than that. This departs
#' from the article's expression for the Wald interval only when
#' \eqn{\lceil 4 z_{1-\alpha/2}/\omega \rceil < p + 2}, that is, for a wide
#' target with many predictors; it never happens in the conditions of the
#' article's tables.
#'
#' For the log transformed interval the article's pilot sample size of
#' \eqn{p + 2} lets the procedure stop within a few observations of the
#' pilot when the sample value happens to be near 0 or near 1, and the
#' article's Tables 11 and 12 show the loss of coverage that follows. With
#' either interval, a pilot sample that is not much larger than the number
#' of predictors gives a sample value that is far too high, and the
#' procedure can then stop at once with an interval that misses the
#' population value. Raising \code{m_0} guards against both.
#' \code{\link{seq_aipe_R2}} attaches a note to an interval reported from
#' fewer than \eqn{p \max(20, 2p)}{p * max(20, 2p)} observations (20 per
#' predictor, or
#' \eqn{2p} per predictor with more than 10 predictors); an \code{m_0} of
#' that size keeps every check out of that range, at the price of a larger
#' sample whenever the stopping rule would have been met earlier. In
#' simulations with 5, 10, and 20 predictors, population values of .30 to
#' .90, widths of .20 and .30, and a 95\% level, coverage with the Wald
#' interval was then between .896 and .932, where the default pilot sample
#' size gave between .019 and .932. See the measured coverage in the help
#' for \code{\link{seq_aipe_R2}}, and use
#' \code{\link{seq_aipe_R2_sensitivity}} to check a planned design.
#'
#' @param width The target full width of the confidence interval for the
#'   population squared multiple correlation coefficient (\eqn{\omega} in the
#'   article). A single number greater than 0 and less than 1.
#' @param p The number of predictors in the regression (\eqn{k} in the
#'   article). A whole number of at least 1.
#' @param conf_level The confidence level (default 0.95). A single number
#'   strictly between 0 and 1.
#' @param ci_method Which confidence interval the procedure is built on:
#'   \code{"wald"} (the default) or \code{"log_transformed"}. See
#'   \code{\link{seq_aipe_R2}}.
#' @param m_0 The smallest pilot sample size allowed. A whole number of at
#'   least \code{p + 2}, which is also the default. (The other pilot sample
#'   size functions in this package default to 4; here the smallest usable
#'   sample depends on the number of predictors.)
#'
#' @return A length one integer: the pilot sample size.
#'
#' @author
#' Ken Kelley \email{kkelley@@nd.edu},
#' Francis Bilson Darku \email{fbilsond@@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @references
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
#' accuracy in parameter estimation for population correlation
#' coefficients. \emph{Psychological Methods}, \emph{24}(4), 492--515.
#' \doi{10.1037/met0000203}
#'
#' Kelley, K. (2008). Sample size planning for the squared multiple
#' correlation coefficient: Accuracy in parameter estimation via narrow
#' confidence intervals. \emph{Multivariate Behavioral Research},
#' \emph{43}(4), 524--555. \doi{10.1080/00273170802490632}
#'
#' @examples
#' # Three predictors of job satisfaction; a 95% confidence interval for the
#' # population squared multiple correlation no wider than 0.20.
#' seq_aipe_R2_pilot(width = 0.20, p = 3)
#'
#' # The log transformed interval starts from p + 2 unless m_0 raises it.
#' seq_aipe_R2_pilot(width = 0.20, p = 3, ci_method = "log_transformed")
#' seq_aipe_R2_pilot(width = 0.20, p = 3, ci_method = "log_transformed",
#'                   m_0 = 30)
#'
#' @seealso \code{\link{seq_aipe_R2}}, \code{\link{seq_aipe_R2_sensitivity}};
#'   in \pkg{DMAR}, the fixed sample size planner
#'   \code{\link[DMAR]{ss_aipe_R2}}.
#' @family seq_aipe_R2
#' @export
seq_aipe_R2_pilot <- function(width, p, conf_level = 0.95,
                              ci_method = c("wald", "log_transformed"),
                              m_0 = p + 2) {
  ci_method <- .check_choice(ci_method, .R2_ci_methods, "ci_method")
  width <- .R2_check_width(width)
  p <- .R2_check_p(p)
  conf_level <- .check_conf_level(conf_level)
  m_0 <- .check_count(m_0, "m_0", min = p + 2L)
  z <- .aipe_z(conf_level)
  .R2_pilot(width, p, z, ci_method, m_0)
}


#' Sequential AIPE Check for the Squared Multiple Correlation Coefficient
#'
#' Given the predictors and the response collected so far, evaluates whether
#' the sequential accuracy in parameter estimation (AIPE) stopping rule of
#' Kelley, Darku, and Chattopadhyay (2019) has been satisfied for the
#' population squared multiple correlation coefficient, and returns the
#' current estimate with its confidence interval. Use it when the goal is a
#' confidence interval for the population squared multiple correlation that
#' is no wider than a chosen width and the population value needed to plan a
#' fixed sample size is unknown.
#'
#' \strong{This procedure assumes multivariate normality.} Unlike the
#' correlation procedures of the same article (\code{\link{seq_aipe_r}}) and
#' the other \code{seq_aipe_*} functions, it is not distribution free: the
#' article develops it for multivariate normal data only (p. 503), because
#' no distribution free interval for the squared multiple correlation was
#' available. See the section on assumptions below.
#'
#' \strong{The coverage guarantee is asymptotic and can fail badly in small
#' samples, above all with many predictors.} Coverage approaches the
#' confidence level as the target width narrows and the sample grows for
#' interior population values, \eqn{0 < \mathrm{P}^2 < 1}, under the
#' stated assumptions. At the null boundary \eqn{\mathrm{P}^2 = 0}, the
#' normal approximation is degenerate; increasing the sample size does not
#' restore the usual Wald calibration. See the section on the null boundary.
#' At a wide target it falls short even with few predictors: with five
#' predictors, a 95\% level, and a target width of .30 it was .851 for a
#' population value of .50 and .797 for a population value of .90. It falls
#' much further when the number of predictors is not small relative to the
#' sample size that the width calls for. The sample \eqn{R^2} overestimates
#' the population value by roughly \eqn{(1 - \mathrm{P}^2) p / n}, both
#' intervals are built around it, and the stopping rule does not correct for
#' it, so the rule is met with an interval that tends to lie above the
#' population value. In
#' simulations of the full procedure with multivariate normal data and a
#' 95\% level (10,000 replications per condition), coverage with the default
#' interval for a population value of .50 was .865 with 10 predictors and
#' .790 with 15 at a target width of .20, and .668 with 10 and .361 with 15
#' at a target width of .30. With 20 predictors unrelated to the response (a
#' population value of 0) and a target width of .30, the final interval
#' contained 0 in 4\% of replications; in 36\% of data sets the rule was met
#' at the very first check, with 27 observations, a sample value near .88,
#' and a lower limit near .80. The article offers no remedy, and the function
#' does not alter the stopping rule to correct for it. What it does is say
#' so: whenever an
#' interval is reported from a sample that is small relative to the number
#' of predictors (fewer than \eqn{p \max(20, 2p)}{p * max(20, 2p)}
#' observations), the object
#' returned carries a note (see the section on small samples below). The
#' article's own simulations cover 2, 5, and 10 predictors at widths of .05
#' and .10. Before relying on the procedure, run
#' \code{\link{seq_aipe_R2_sensitivity}} at the planned width and number of
#' predictors; a larger pilot sample (\code{m_0} in
#' \code{\link{seq_aipe_R2_pilot}}) helps. See the sections on the
#' asymptotic guarantee and on small samples below.
#'
#' The function is stateless. After collecting the pilot sample, whose
#' size \code{\link{seq_aipe_R2_pilot}} gives, the user calls
#' \code{seq_aipe_R2()} and, while the rule is not met, adds observations
#' and calls it again on all the data collected so far.
#'
#' \code{\link{seq_path}} applies the function to data already in hand at
#' every sample size from the pilot sample size on. For most procedures it
#' works out the pilot sample size itself. For this one it cannot, because
#' the pilot sample size depends on the number of predictors, which the
#' check function reads from the data and \code{seq_path()} has no way to
#' hand to the pilot function; called without \code{start} it stops with a
#' message. Supply the start, as in the last example:
#' \preformatted{seq_path(seq_aipe_R2, list(x = x, y = y), width = width,
#'          start = seq_aipe_R2_pilot(width, p = ncol(x)))}
#' The printed path gives the stopping point and the interval there but not
#' the notes of the check at that point, so a path that stops while the
#' sample is small relative to the number of predictors (see the section on
#' small samples below) prints without that note. The check at the stopping
#' point, with its notes, is \code{attr(path, "final_check")}; look at it
#' before the interval is reported.
#'
#' @section Estimator and Confidence Intervals:
#' With \eqn{n} observations on \eqn{p} predictors (\eqn{k} in the article)
#' and a response, the estimate is the sample squared multiple correlation
#' coefficient \eqn{R^2} (Equation 57 of the article), the value
#' \code{summary(lm(y ~ x))$r.squared} reports. Write \eqn{R} for its
#' positive square root, \eqn{\omega} for \code{width}, and
#' \eqn{z = z_{1-\alpha/2}} for the standard normal quantile.
#'
#' The article studies two intervals, both with \eqn{R} in place of the
#' population value in the standard error:
#' \itemize{
#'   \item \code{ci_method = "wald"} (Equation 59), built on the asymptotic
#'     variance \eqn{4 \mathrm{P}^2 (1 - \mathrm{P}^2)^2 / n} of Olkin and
#'     Finn (1995, p. 161, Equation 10):
#'     \deqn{R^2 \pm z \, \frac{2 R (1 - R^2)}{\sqrt{n}}.}
#'   \item \code{ci_method = "log_transformed"} (Equation 58), the interval
#'     of Bonett and Wright (2011), formed for \eqn{\ln(1 - R^2)} and
#'     transformed back:
#'     \deqn{1 - \exp\left\{\ln(1 - R^2) \pm z \, \frac{2 R}{\sqrt{n - p - 2}}\right\}.}
#' }
#' Neither interval is truncated to \eqn{[0, 1]} in the article, and the
#' limits are returned as computed: with either interval the lower limit can
#' be negative, and with the Wald interval the upper limit can exceed 1 when
#' the sample is very small (\eqn{n < 4 z^2 R^2}). Once the stopping rule is
#' met the interval is no wider than \eqn{\omega}. The log transformed
#' interval is not defined at \eqn{n = p + 2}, where its limits are returned
#' as \code{NA}.
#'
#' The Wald interval is the default because the article concludes from its
#' simulations (p. 506) that the procedure built on it did better than the
#' one built on the log transformed interval in terms of the sample size
#' required; its coverage is also closer to the confidence level in 59 of
#' the 60 conditions of the article's Tables 11 to 14 (see the section on
#' the asymptotic guarantee).
#'
#' @section Stopping Rule:
#' The rule is met when \eqn{n} is at least \code{criterion_n}, the smallest
#' whole number not below
#' \itemize{
#'   \item for \code{"wald"}:
#'     \deqn{\frac{16 z^2}{\omega^2} \left(R^2 + \frac{1}{n}\right) \left(1 - R^2 + \frac{1}{n}\right)^2;}
#'   \item for \code{"log_transformed"} (Equation 62):
#'     \deqn{p + 2 + \frac{4 R^2 z^2}{L^2}, \quad L = \ln\left[\frac{1}{2}\left(\frac{\omega}{1 - R^2} + \sqrt{\frac{\omega^2}{(1 - R^2)^2} + 4}\right)\right].}
#' }
#' These estimate the sample sizes that would be needed if the population
#' value \eqn{\mathrm{P}^2} were known (Equations 61 and 60):
#' \eqn{16 z^2 \mathrm{P}^2 (1 - \mathrm{P}^2)^2 / \omega^2} for the Wald
#' interval, and the log transformed expression with \eqn{\mathrm{P}^2} in
#' place of \eqn{R^2}. The terms \eqn{1/n} in the Wald rule are the
#' correction that keeps sampling from stopping too early; the log
#' transformed rule has no such correction in the article.
#'
#' \strong{The Wald rule differs from Equation 63 as printed.} The printed
#' equation has \eqn{(1 - (R^2 + 1/n))^2} as its last factor, which
#' \emph{lowers} the criterion whenever \eqn{R^2} is above about one third,
#' allows the interval at the stopping point to be wider than \eqn{\omega},
#' and does not reproduce the article's Tables 13 and 14. In simulations of
#' all 60 conditions of those tables the printed form gave a mean final
#' sample size below the tabled value in every condition, by more than three
#' standard errors in 49 of them (and by as much as 39.7), and a mean
#' interval width above \eqn{\omega} in 27. The form above, with both factors
#' increased by \eqn{1/n}, reproduced all 60 conditions within Monte Carlo
#' error; it never lowers the criterion, and it guarantees a width of at
#' most \eqn{\omega} at the stopping point. The package implements the form
#' that reproduces the published results.
#'
#' \strong{When no decision can be made the rule is "not met".} As in the
#' other \code{seq_aipe_*} check functions, when the variance estimate is
#' unusable or the estimate cannot be computed, the rule is reported as not
#' met, the limits are \code{NA}, \code{criterion_n} is at least one more
#' than \code{current_n}, the attribute \code{variance_unusable} is
#' \code{TRUE}, and a note says why (see \code{\link{smsd_check_methods}}).
#' This is a deliberate departure from the rules as printed, which would be
#' met at once with an interval of zero width.
#'
#' For this procedure the variance estimate, \eqn{4 R^2 (1 - R^2)^2}, is
#' unusable when the sample \eqn{R^2} is 0 or 1 to within the precision of
#' the data (for example because the response is an exact linear function
#' of the predictors), which makes the estimate zero, and it cannot be
#' computed when the response does not vary. A positive estimate, however
#' small, is used as it stands: the lower bound of \eqn{n^{-3}} that the
#' distribution free procedures in this package place on their variance
#' estimates is not part of this procedure. The rule is also reported as
#' not met, with limits of \code{NA} and a note, while the log transformed
#' interval is not yet defined (\eqn{n = p + 2}).
#'
#' A variable whose values agree to ten significant digits counts as not
#' varying, as everywhere in this package: such a predictor is an error, and
#' such a response is the case in which the estimate cannot be computed.
#' Whether \eqn{R^2} is 0 or 1 is likewise judged to the precision of the
#' data. For ordinary data that is the precision of the arithmetic. When a
#' variable is a large constant plus a small spread (a mean more than about
#' 700,000 times the standard deviation, as with long identifiers or time
#' stamps), the last digits of its stored values are rounding error, and a
#' share of the sum of squares small enough to be that error alone is
#' treated as zero. Subtracting a constant from such a variable before
#' calling avoids both and does not change \eqn{R^2}.
#'
#' @section Assumptions:
#' The predictors and the response are assumed to be independent draws from
#' a multivariate normal distribution (so the predictors are random, not
#' fixed by design). Both standard errors come from the normal theory
#' variance of \eqn{R^2}, which is too small for heavy tailed data, so the
#' intervals are then too narrow however large the sample. In simulations of
#' the full procedure with three predictors, a population value of .50, a
#' target width of .10, and a 95\% level (40,000 replications per
#' condition), coverage with the Wald interval was .946 for multivariate
#' normal data, .902 for a multivariate \emph{t} distribution with 10
#' degrees of freedom, .788 for a multivariate \emph{t} distribution with 5
#' degrees of freedom, and .860 when the predictors and the error were
#' independent right skewed variables (standardized chi square with 3
#' degrees of freedom). A narrower target does not help: at a width of .05
#' coverage was .948 for normal data and .782 for the \emph{t} distribution
#' with 5 degrees of freedom. The log transformed interval behaved the same
#' way (.938, .896, .773, and .851 at a width of .10). The article offers no
#' remedy, and this function should not be used for data that are far from
#' multivariate normal.
#'
#' @section The Null Boundary:
#' When all \eqn{p} predictors and the response are independent normal
#' variables and the population squared multiple correlation is zero,
#' the fixed sample distribution is
#' \eqn{R^2 \sim \mathrm{Beta}(p/2, (n-p-1)/2)}. In particular,
#' \eqn{n R^2} approaches a chi squared distribution with \eqn{p} degrees
#' of freedom, rather than the normal approximation with positive variance
#' used by the Wald interval. Its lower limit includes zero only when
#' \eqn{n R^2 \leq 4 z^2 (1-R^2)^2}. Its limiting fixed sample coverage
#' is therefore \eqn{\Pr(\chi^2_p \leq 4 z^2)}, which is about 24.49\%
#' for \eqn{p=20} and a nominal 95\% level. This is a fixed sample limit,
#' not the coverage of the sequential stopping rule. It shows why a larger
#' pilot alone cannot establish nominal coverage at the null boundary.
#' Moreover, when \eqn{p=20}, no legal fixed sample size gives 95\%
#' coverage with this Wald interval. Put \eqn{c=4z^2}. The inclusion event
#' implies \eqn{R^2 \leq c/n}. The beta distribution's mean and variance,
#' together with a one sided variance inequality, bound its probability by
#' \eqn{40/(40+(20-c)^2)}, or about 65.07\%, for every \eqn{n>21}.
#' This loose bound rules out 95\% coverage at any fixed sample size;
#' the limit above gives the sharper large sample result.
#'
#' @section Asymptotic Guarantee:
#' For fixed \eqn{p} and an interior population value
#' \eqn{0 < \mathrm{P}^2 < 1}, the coverage of the interval at the stopping
#' point approaches the confidence level as the target width narrows;
#' it is not guaranteed at any
#' particular width, and the article notes (p. 508) that the procedure works
#' well as the sample size increases but not for small sample sizes. With
#' five predictors, multivariate normal data, and a 95\% level (10,000
#' replications per condition), coverage with the Wald interval for a
#' population value of .50 was .851 at a target width of .30, .912 at .20,
#' .936 at .10, and .949 at .05; with the log transformed interval it was
#' .831, .896, .930, and .939. Among these interior conditions, Wald
#' coverage is lowest when the population value is large: for .90 it was
#' .797 at a width of .30, .848 at .20, .901 at .10, and .913 at .05. The log
#' transformed interval, which can stop within a few observations of its
#' pilot sample size, does not improve steadily as the target narrows: for a
#' population value of .90 it gave .925, .906, .774, and .847 at the same
#' four widths. The article's Tables 11 to 14 show the same pattern at
#' widths of .05 and .10: at the 95\% level its coverage ranges from .830 to
#' .954 for the Wald interval and from .646 to .949 for the log transformed
#' interval.
#'
#' Coverage also falls, much further, when the sample is small relative to
#' the number of predictors; see the next section. (The Monte Carlo standard
#' error of each coverage from a simulation of the full procedure reported
#' on this page for a single condition is at most .011, and at most .005
#' where 10,000 replications are stated; a condition that appears in two
#' places comes from two separate runs.)
#'
#' @section Small Samples Relative to the Number of Predictors:
#' The sample \eqn{R^2} overestimates the population value by roughly
#' \eqn{(1 - \mathrm{P}^2) p / n}. In standard errors that is about
#' \eqn{p / (2 \mathrm{P} \sqrt{n})}, and both intervals are built around
#' the sample value. Their coverage is therefore below the confidence level
#' until \eqn{n} is large relative to \eqn{p} and, with many predictors,
#' relative to \eqn{p^2}. Olkin and Finn (1995, p. 161) give their variance
#' as an approximation for samples in which, in the notation used here,
#' \eqn{2p + 5} is small relative to \eqn{n}. The stopping rule takes no
#' account of this: with a wide target or many predictors it can be met at
#' the pilot sample size.
#'
#' \strong{The note.} Whenever an interval is reported from fewer than
#' \eqn{p \max(20, 2p)}{p * max(20, 2p)} observations, the object returned
#' carries a note
#' that says so. That is fewer than 20 observations per predictor or, with
#' more than 10 predictors, fewer than \eqn{2p} per predictor: 40
#' observations with 2 predictors, 100 with 5, 200 with 10, 800 with 20. The
#' note depends on \eqn{n} and \eqn{p} alone. It does not change
#' \code{is_satisfied}: the article has no remedy for the shortfall, and
#' none is invented here. It marks the sample sizes at which the interval
#' should not be relied on without checking the design with
#' \code{\link{seq_aipe_R2_sensitivity}}.
#'
#' The threshold is where the coverage of the Wald interval at a
#' \emph{fixed} sample size and a 95\% level reaches about .90 for
#' population values of .30 or more. With 2 to 40 predictors and population
#' values of .30, .50, .70, and .90 (400,000 samples per condition from the
#' exact distribution of \eqn{R^2} under multivariate normality, which
#' Gatsonis and Sampson, 1989, p. 517, give), that coverage was between .902
#' and .928 at exactly \eqn{p \max(20, 2p)}{p * max(20, 2p)} observations,
#' between .847 and .904 at half as many, and between .714 and .855 at a
#' quarter as many; at twice as many it was between .927 and .938. (With
#' one predictor it was between .864 and .891 at 20 observations. At the
#' 90\% and 99\% levels it was between .834 and .870 and between .953 and
#' .981 at the threshold.) A smaller population value needs a larger sample:
#' for a population value of .10 the coverage at the threshold was between
#' .834 and .934. The log transformed interval does better at a fixed sample
#' size (between .912 and .950 at the threshold and between .799 and .957 at
#' a quarter of it), but the sequential procedure built on it stops earlier.
#'
#' In simulations of the full sequential procedure at a 95\% level
#' (population values of .10, .30, .50, .70, and .90; 500 to 2,000
#' replications per condition; every condition given the same weight),
#' coverage among the replications that stopped with the note depended
#' above all on the number of predictors. With the Wald interval it was
#' .793 with 2 predictors, .817 with 5, .663 with 10, .428 with 20, and
#' .180 with 30. Among the replications that stopped without the note it
#' was .915 with 1 predictor (with which the note was never attached), .933
#' with 2, .942 with 5, .917 with 10, .927 with 20, and .892 with 30. With
#' the log transformed interval coverage with the note was .466, .732, .807,
#' .755, .615, and .422 with 1, 2, 5, 10, 20, and 30 predictors, and
#' without it between .923 and .951. The conditions are not fully crossed:
#' the target widths were .05, .07, .10, .15, .20, .30, and .40 with 1, 2,
#' and 5 predictors, .10 to .40 with 10 predictors, .07 to .40 with 20, and
#' .05, .20, and .30 with 30 (175 conditions for each interval). A single
#' figure for all replications that stopped with the note therefore depends
#' on which conditions are pooled, and none is offered: for the Wald
#' interval it is .530 over all 175 conditions, .811 over those with 1, 2,
#' and 5 predictors, and .439 over those with widths of .20 and .30, which
#' every number of predictors has. The figure without the note depends on
#' that far less (.926, .928, and .914).
#'
#' By the sample size at the stopping point relative to the threshold, for
#' population values of .30 or more and over all of those conditions,
#' coverage with the Wald interval was .780 from a quarter to a half of the
#' threshold, .877 from a half to the threshold, .901 from the threshold to
#' twice it, .928 from twice to four times it, and .949 beyond; with the log
#' transformed interval it was .815, .911, .938, .946, and .951. Below a
#' quarter of the threshold it again depended on the number of predictors:
#' .601 with 5, .395 with 10, .243 with 20, and .115 with 30 for the Wald
#' interval (with 1 and 2 predictors no replication stopped there), and
#' .239, .587, .643, .607, .506, and .399 with 1, 2, 5, 10, 20, and 30 for
#' the log transformed interval.
#'
#' The note is a statement about \eqn{n} and \eqn{p} and nothing else, so
#' it is neither a verdict nor a guarantee. Coverage can be at the
#' confidence level with it: with five predictors, a population value of
#' .10, and a width of .30, 90\% of replications stopped with the note and
#' coverage was .941. And coverage can be below the confidence level without
#' it, for the reasons given in the previous section: at a wide target and
#' for a small population value. Among the replications with 1, 2, and 5
#' predictors that stopped without the note, Wald coverage for a population
#' value of .50 was .949 at a width of .10 and .918 at .30, and for a
#' population value of .10 it was between .864 and .942 across the seven
#' widths. With a single predictor, for which the pilot sample size at
#' these widths is 20 or more, so that the note is never attached, coverage
#' for a population value of .10 was between .80 and .89 at widths of .10
#' to .40.
#'
#' Some single conditions, each from a run of its own. For a population
#' value of .50 at the 95\% level (10,000 replications per condition; 4,000
#' for 40 predictors at a width of .10), coverage with the Wald interval at
#' a target width of .30 was .900 with 2 predictors, .845 with 5, .668 with
#' 10, .361 with 15, .071 with 20, and 0 with 40 (where every replication
#' stopped at the pilot sample size of 42); at a width of .20 it was .865
#' with 10 predictors and .790 with 15; at a width of .10 it was .948, .942,
#' .934, .903, and .778 with 2, 5, 10, 20, and 40. The log transformed
#' interval gave .885, .837, .709, .563, .427, and .076 at a width of .30,
#' .829 and .748 at .20, and .936, .926, .906, .851, and .671 at .10. When
#' the predictors are unrelated to the response (a population value of 0),
#' coverage at a width of .30 was .971 with 5 predictors, .711 with 10, and
#' .040 with 20 for the Wald interval, and .967, .847, and .273 for the log
#' transformed interval. This case does not improve with the sample size,
#' because the variance \eqn{4 \mathrm{P}^2 (1 - \mathrm{P}^2)^2} on which
#' both intervals rest is zero there: the Wald interval contains 0 when
#' \eqn{n R^2 \le 4 z^2 (1 - R^2)^2}, so as \eqn{n} grows its coverage
#' tends to the probability that a chi square variable with \eqn{p} degrees
#' of freedom is at most \eqn{4 z^2}, which at the 95\% level is .991, .881,
#' .575, and .245 for 5, 10, 15, and 20 predictors (at a fixed sample size
#' of 400 observations per predictor, 100,000 samples each: .990, .877,
#' .569, and .240).
#'
#' \strong{A larger pilot sample.} The one thing in the user's hands is
#' \code{m_0} in \code{\link{seq_aipe_R2_pilot}}. With 20 predictors, a
#' population value of .50, and a width of .30, pilot samples of 100, 200,
#' and 400 gave coverage of .601, .789, and .877 with the Wald interval.
#' With the pilot sample size at the threshold itself,
#' \code{m_0 = p * max(20, 2 * p)}, so that no interval is ever reported
#' with the note, coverage with the Wald interval was between .896 and .932
#' for 5, 10, and 20 predictors, population values of .30 to .90, and widths
#' of .20 and .30 (2,000 replications per condition), where the default
#' pilot sample size gave between .019 and .932; for a population value of
#' .10 it was between .853 and .953. The log transformed interval gave
#' between .920 and .955 (default pilot sample size: .373 to .930), and
#' between .875 and .966 for a population value of .10. The price is the
#' sample size: with 20 predictors, a population value of .50, and a width
#' of .30 that pilot sample is 800 observations, where the rule by itself
#' stops near 33. This is not part of the article's procedure and is not a
#' default; the article's own simulations cover 2, 5, and 10 predictors at
#' widths of .05 and .10. Use \code{\link{seq_aipe_R2_sensitivity}} to see
#' what to expect at the width, number of predictors, and pilot sample size
#' of a planned study before relying on the procedure.
#'
#' @param x The predictors collected so far: a numeric matrix or a data
#'   frame of numeric columns with one row per observation and one column
#'   per predictor, or a numeric vector for a single predictor. Do not
#'   include a column of ones; the intercept is part of the model.
#' @param y Numeric vector. The response collected so far, one value per
#'   row of \code{x}.
#' @param width The target full width of the confidence interval for the
#'   population squared multiple correlation coefficient (\eqn{\omega} in the
#'   article). A single number greater than 0 and less than 1.
#' @param conf_level The confidence level (default 0.95). A single number
#'   strictly between 0 and 1.
#' @param ci_method Which confidence interval to report and to build the
#'   stopping rule on: \code{"wald"} (the default; the symmetric interval
#'   with the Olkin and Finn, 1995, variance) or \code{"log_transformed"}
#'   (the Bonett and Wright, 2011, interval).
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} with
#'   columns \code{term} and \code{value} and six rows:
#'   \describe{
#'     \item{\code{lower_limit}}{The lower confidence limit.}
#'     \item{\code{R2}}{The sample squared multiple correlation coefficient.}
#'     \item{\code{upper_limit}}{The upper confidence limit.}
#'     \item{\code{xi_squared}}{The estimated asymptotic variance of
#'       \eqn{R^2} under multivariate normality,
#'       \eqn{4 R^2 (1 - R^2)^2} (Olkin and Finn, 1995). The Wald standard
#'       error is its square root over \eqn{\sqrt{n}}. It is reported as
#'       computed; no lower bound is applied to it.}
#'     \item{\code{criterion_n}}{The criterion sample size of the stopping
#'       rule, given the current data. The rule is met when
#'       \code{current_n} is at least this number. When no decision can be
#'       made it is at least \code{current_n + 1}.}
#'     \item{\code{current_n}}{The number of observations supplied.}
#'   }
#'   A quantity that is not defined for the current data is \code{NA}: the
#'   two limits whenever no decision can be made (the variance estimate is
#'   unusable, the estimate cannot be computed, or
#'   \code{ci_method = "log_transformed"} and \eqn{n = p + 2}), and
#'   \code{R2} and \code{xi_squared} when the response does not vary. The
#'   attributes are:
#'   \describe{
#'     \item{\code{method}}{The label of the method, which names the
#'       interval: \code{"squared multiple correlation, Wald interval"} or
#'       \code{"squared multiple correlation, log transformed interval"}.}
#'     \item{\code{is_satisfied}}{\code{TRUE} when the stopping rule is
#'       met, otherwise \code{FALSE}; never \code{NA}. It is \code{TRUE}
#'       exactly when \code{current_n} is at least \code{criterion_n}.}
#'     \item{\code{current_n}, \code{criterion_n}}{The two sample sizes, as
#'       integers.}
#'     \item{\code{conf_level}, \code{width}}{As supplied.}
#'     \item{\code{notes}}{A character vector, of length zero when there is
#'       nothing to say, that \code{print()} shows under the stopping line.
#'       The notes that say why no decision can be made (the estimate cannot
#'       be computed; the variance estimate is unusable; the log transformed
#'       interval is not yet defined) accompany a rule that is not met. The
#'       note that the sample is small relative to the number of predictors
#'       can accompany either decision and does not change it.}
#'     \item{\code{ci_method}}{The interval used, \code{"wald"} or
#'       \code{"log_transformed"}.}
#'     \item{\code{p}}{The number of predictors.}
#'     \item{\code{variance_unusable}}{\code{TRUE} when the variance
#'       estimate is unusable (zero, or not finite because the estimate
#'       cannot be computed), otherwise \code{FALSE}.}
#'   }
#'
#' @author
#' Ken Kelley \email{kkelley@@nd.edu},
#' Francis Bilson Darku \email{fbilsond@@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @references
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
#' accuracy in parameter estimation for population correlation
#' coefficients. \emph{Psychological Methods}, \emph{24}(4), 492--515.
#' \doi{10.1037/met0000203}
#'
#' Kelley, K. (2008). Sample size planning for the squared multiple
#' correlation coefficient: Accuracy in parameter estimation via narrow
#' confidence intervals. \emph{Multivariate Behavioral Research},
#' \emph{43}(4), 524--555. \doi{10.1080/00273170802490632}
#'
#' Bonett, D. G., & Wright, T. A. (2011). Sample size requirements for
#' multiple regression interval estimation. \emph{Journal of Organizational
#' Behavior}, \emph{32}(6), 822--830. \doi{10.1002/job.717}
#'
#' Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact power
#' and sample size calculations. \emph{Psychological Bulletin},
#' \emph{106}(3), 516--524.
#'
#' Olkin, I., & Finn, J. D. (1995). Correlations redux. \emph{Psychological
#' Bulletin}, \emph{118}(1), 155--164. \doi{10.1037/0033-2909.118.1.155}
#'
#' @examples
#' # How well do autonomy, supervisor support, and workload together predict
#' # job satisfaction? The goal is a 95% confidence interval for the
#' # population squared multiple correlation no wider than 0.25.
#' set.seed(113)
#' survey <- function(n) {
#'   predictors <- cbind(autonomy = rnorm(n), support = rnorm(n),
#'                       workload = rnorm(n))
#'   satisfaction <- drop(predictors %*% c(0.45, 0.40, -0.30)) + rnorm(n)
#'   list(x = predictors, y = satisfaction)
#' }
#' m <- seq_aipe_R2_pilot(width = 0.25, p = 3)
#' employees <- survey(m)
#' out <- seq_aipe_R2(employees$x, employees$y, width = 0.25)
#' # Not met, with a note: 32 observations are few for three predictors.
#' out
#'
#' # Survey ten more employees at a time until the stopping rule is met.
#' while (!attr(out, "is_satisfied")) {
#'   more <- survey(10)
#'   employees <- list(x = rbind(employees$x, more$x),
#'                     y = c(employees$y, more$y))
#'   out <- seq_aipe_R2(employees$x, employees$y, width = 0.25)
#' }
#' out
#'
#' # The same data one employee at a time, as a path. seq_path() cannot work
#' # out the pilot sample size for this procedure (it depends on the number
#' # of predictors), so the start is supplied.
#' path <- seq_path(seq_aipe_R2, employees, width = 0.25,
#'                  start = seq_aipe_R2_pilot(width = 0.25, p = 3))
#' attr(path, "stop_n")
#' # The check at the stopping point. Its notes, if it has any, are shown
#' # here and not when the path itself is printed.
#' attr(path, "final_check")
#'
#' @seealso \code{\link{seq_aipe_R2_pilot}},
#'   \code{\link{seq_aipe_R2_sensitivity}}; \code{\link{seq_path}} for the
#'   whole stopping path of data in hand; \code{\link{smsd_check_methods}}
#'   for the methods of the returned object; \code{\link{seq_aipe_r}} for the
#'   distribution free procedures for a single correlation; in \pkg{DMAR},
#'   the fixed sample size planner \code{\link[DMAR]{ss_aipe_R2}} and the
#'   confidence interval \code{\link[DMAR]{ci_R2}} (Kelley, 2008), which
#'   invert the exact sampling distribution of \eqn{R^2} rather than use an
#'   approximate standard error.
#' @family seq_aipe_R2
#' @export
seq_aipe_R2 <- function(x, y, width, conf_level = 0.95,
                        ci_method = c("wald", "log_transformed")) {
  ci_method <- .check_choice(ci_method, .R2_ci_methods, "ci_method")
  est <- .R2_estimate(x, y)
  width <- .R2_check_width(width)
  conf_level <- .check_conf_level(conf_level)
  z <- .aipe_z(conf_level)
  n <- est$n
  p <- est$p

  rule <- .R2_rule(est$R2, est$one_minus_R2, n, p, width, z, ci_method,
                   est$resolution)
  criterion_n <- .as_sample_size(rule$criterion, "criterion sample size",
                                 advice = .R2_width_advice)

  # Why no decision can be made, when none can (the notes of the whole
  # package for an estimate that is not finite and for a variance estimate
  # that is unusable), and, when there is an interval, whether the sample is
  # small relative to the number of predictors.
  notes <- c(
    if (!rule$estimable) .note_not_estimable(),
    if (!rule$estimable && est$response_differs) {
      .R2_ten_digits("The values of the response differ")
    },
    if (rule$unusable) .note_variance_unusable(),
    if (!rule$unusable && !rule$defined) {
      paste0("The log transformed interval is not defined with p + 2 ",
             "observations (its standard error divides by n - p - 2), so ",
             "the limits are NA and the stopping rule is reported as not ",
             "met. Collect more data and check again.")
    },
    if (rule$defined && n < .R2_small_sample_below(p)) {
      .note_R2_small_sample(n, p)
    }
  )
  if (is.null(notes)) notes <- character(0)

  new_smsd_check(
    term = c("lower_limit", "R2", "upper_limit", "xi_squared",
             "criterion_n", "current_n"),
    value = c(rule$lower, if (rule$estimable) est$R2 else NA_real_,
              rule$upper, if (rule$estimable) rule$xi2 else NA_real_,
              criterion_n, n),
    method = paste0("squared multiple correlation, ",
                    if (ci_method == "wald") "Wald" else "log transformed",
                    " interval"),
    is_satisfied = rule$satisfied,
    current_n = n,
    criterion_n = criterion_n,
    conf_level = conf_level,
    width = width,
    extra = list(ci_method = ci_method, p = p,
                 coverage_guarantee = "asymptotic",
                 variance_unusable = rule$unusable),
    notes = notes
  )
}


#' Monte Carlo Sensitivity Analysis for Sequential AIPE: Squared Multiple Correlation Coefficient
#'
#' Runs the sequential accuracy in parameter estimation (AIPE) procedure for
#' the squared multiple correlation coefficient
#' (\code{\link{seq_aipe_R2_pilot}}, then \code{\link{seq_aipe_R2}} until
#' the stopping rule is met) \code{G} times on multivariate normal data with
#' a known population value, and summarizes the realized sample size, the
#' realized confidence interval width, and the empirical coverage. Use it
#' before a study to see what sample size to expect and how close the
#' coverage is to the confidence level at the planned width and number of
#' predictors: the procedure's guarantee is asymptotic, and its coverage
#' falls, in some conditions far below the confidence level, as the number
#' of predictors grows relative to the sample size the width calls for (see
#' \code{\link{seq_aipe_R2}}, which attaches a note to an interval
#' reported from fewer than \eqn{p \max(20, 2p)}{p * max(20, 2p)}
#' observations; a design
#' whose \code{mean_n} and \code{median_n} here are below that number is one
#' in which a study will usually end with that note).
#'
#' Each replication draws \eqn{p} independent standard normal predictors and
#' a response equal to \eqn{\sqrt{\mathrm{P}^2 / p}} times their sum plus
#' independent normal error with variance \eqn{1 - \mathrm{P}^2}, where
#' \eqn{\mathrm{P}^2} is \code{true_value}. This is the design of the
#' simulations in Kelley, Darku, and Chattopadhyay (2019, p. 506). Under
#' multivariate normality the distribution of the sample squared multiple
#' correlation depends only on \eqn{\mathrm{P}^2}, \eqn{p}, and the sample
#' size, so the results apply to any means, variances, and correlations
#' among the predictors with the same population value.
#'
#' The stopping rule is checked at the pilot sample size and then after
#' every \code{batch} further observations. With \code{batch = 1} this is
#' the purely sequential procedure of the article; a larger value shows the
#' overshoot to expect from a study that checks the rule only after every
#' \code{batch} observations. The function gives the results that a loop of
#' \code{\link{seq_aipe_R2}} calls on the same observations would give, but
#' it evaluates the rule from running sums, so a narrow target that needs
#' thousands of observations per replication remains practical.
#'
#' The function generates multivariate normal data only, the assumption
#' under which the procedure is derived; see \code{\link{seq_aipe_R2}} for
#' measured coverage when that assumption fails.
#'
#' @param true_value The population squared multiple correlation
#'   coefficient under the data generating mechanism. A single number of at
#'   least 0 and less than 1.
#' @param width The target full width of the confidence interval
#'   (\eqn{\omega} in the article). A single number greater than 0 and less
#'   than 1.
#' @param p The number of predictors (\eqn{k} in the article). A whole
#'   number of at least 1.
#' @param conf_level The confidence level (default 0.95).
#' @param ci_method \code{"wald"} (the default) or \code{"log_transformed"},
#'   as in \code{\link{seq_aipe_R2}}.
#' @param G The number of Monte Carlo replications (default 1000). A whole
#'   number of at least 1.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single whole number that \code{set.seed()} accepts
#'   to make the result reproducible; the caller's random number stream is
#'   then left as it was found.
#' @param m_0 The smallest pilot sample size allowed, passed to
#'   \code{\link{seq_aipe_R2_pilot}}. A whole number of at least
#'   \code{p + 2}, which is also the default.
#' @param max_n The largest sample size a replication may reach: a
#'   replication in which the rule is still not met with \code{max_n}
#'   observations stops there. It must be at least the pilot sample size,
#'   and at least \code{p + 3} with the log transformed interval (which is
#'   not defined with \code{p + 2} observations). The cap is exact for any
#'   \code{batch} (the last batch is shortened).
#' @param batch The number of observations added between checks of the
#'   stopping rule after the pilot sample (default 1).
#'
#' @return An object of class \code{c("smsd_sensitivity", "data.frame")}:
#'   a data frame with columns \code{term} and \code{value}, which can be
#'   used as any data frame and which prints its values in fixed notation.
#'   See \code{\link{smsd_sensitivity_methods}} for the print method. It
#'   has one row per summary, in this order:
#'   \describe{
#'     \item{\code{true_value}}{The population squared multiple correlation
#'       supplied.}
#'     \item{\code{target_width}}{The target width supplied.}
#'     \item{\code{conf_level}}{The confidence level supplied.}
#'     \item{\code{replications}}{The number of replications run.}
#'     \item{\code{mean_n}, \code{median_n}, \code{sd_n}}{The mean, median,
#'       and standard deviation of the realized sample size. \code{sd_n} is
#'       \code{NA} when \code{G} is 1.}
#'     \item{\code{mean_width}, \code{median_width}}{The mean and median of
#'       the realized confidence interval width, over the replications that
#'       ended with an interval.}
#'     \item{\code{pct_width_within}}{The percentage (0 to 100) of
#'       replications whose realized width is at most the target.}
#'     \item{\code{coverage}}{The proportion (0 to 1) of replications whose
#'       interval contains \code{true_value}. A replication stopped by
#'       \code{max_n} at a check that has no interval (see
#'       \code{\link{seq_aipe_R2}}) counts as not covering and as not
#'       within the target width; this cannot happen with the data
#'       generated here.}
#'     \item{\code{pct_hit_max_n}}{The percentage (0 to 100) of replications
#'       stopped by \code{max_n} before the stopping rule was met (a
#'       replication in which the rule is met with exactly \code{max_n}
#'       observations is not counted). When it is above zero the other
#'       summaries are censored and should be read with care.}
#'     \item{\code{pct_variance_unusable}}{The percentage (0 to 100) of
#'       replications in which, at one or more checks up to the stopping
#'       point, the variance estimate was unusable: zero (a sample value of
#'       0 or 1 to within the precision of the data) or not finite (see
#'       \code{\link{seq_aipe_R2}}). For the data generated here this is 0
#'       in practice.}
#'   }
#'
#' @author
#' Ken Kelley \email{kkelley@@nd.edu},
#' Francis Bilson Darku \email{fbilsond@@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @references
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
#' accuracy in parameter estimation for population correlation
#' coefficients. \emph{Psychological Methods}, \emph{24}(4), 492--515.
#' \doi{10.1037/met0000203}
#'
#' Kelley, K. (2008). Sample size planning for the squared multiple
#' correlation coefficient: Accuracy in parameter estimation via narrow
#' confidence intervals. \emph{Multivariate Behavioral Research},
#' \emph{43}(4), 524--555. \doi{10.1080/00273170802490632}
#'
#' @examples
#' # Three predictors, a population squared multiple correlation of .40, and
#' # a 95% confidence interval no wider than 0.25. G is small to keep the
#' # example quick; a reported analysis deserves G = 10000.
#' seq_aipe_R2_sensitivity(true_value = 0.40, width = 0.25, p = 3, G = 100,
#'                         seed = 113)
#'
#' # The same study checking the stopping rule after every 20 participants.
#' seq_aipe_R2_sensitivity(true_value = 0.40, width = 0.25, p = 3, G = 100,
#'                         seed = 113, batch = 20)
#'
#' @seealso \code{\link{seq_aipe_R2}}, \code{\link{seq_aipe_R2_pilot}};
#'   \code{\link{smsd_sensitivity_methods}} for how the result prints; in
#'   \pkg{DMAR}, \code{\link[DMAR]{ss_aipe_R2_sensitivity}} for the fixed
#'   sample size analog.
#' @family seq_aipe_R2
#' @export
seq_aipe_R2_sensitivity <- function(true_value, width, p, conf_level = 0.95,
                                    ci_method = c("wald", "log_transformed"),
                                    G = 1000, seed = NULL, m_0 = p + 2,
                                    max_n = 1e5L, batch = 1L) {
  ci_method <- .check_choice(ci_method, .R2_ci_methods, "ci_method")
  if (missing(true_value)) {
    .smsd_stop("'true_value' must be supplied: the population squared ",
               "multiple correlation coefficient.")
  }
  true_value <- .check_number(true_value, "true_value", lower = 0, upper = 1,
                              upper_open = TRUE)
  width <- .R2_check_width(width)
  p <- .R2_check_p(p)
  conf_level <- .check_conf_level(conf_level)
  G <- .check_count(G, "G")
  m_0 <- .check_count(m_0, "m_0", min = p + 2L)
  batch <- .check_count(batch, "batch")
  max_n <- .check_count(max_n, "max_n")
  seed <- .check_seed(seed)
  z <- .aipe_z(conf_level)
  pilot_n <- .R2_pilot(width, p, z, ci_method, m_0)
  if (max_n < pilot_n) {
    .smsd_stop("'max_n' (", max_n, ") is below the pilot sample size (",
               pilot_n, "). Raise 'max_n' to at least the pilot sample ",
               "size.")
  }
  if (ci_method == "log_transformed" && max_n < p + 3) {
    .smsd_stop("'max_n' must be at least p + 3 (", p + 3, ") for the log ",
               "transformed interval, which is not defined with p + 2 ",
               "observations; got ", max_n, ".")
  }

  slope <- sqrt(true_value / p)
  error_sd <- sqrt(1 - true_value)
  draw <- function(m) {
    obs <- matrix(stats::rnorm(m * (p + 1L)), nrow = m, ncol = p + 1L,
                  byrow = TRUE)
    obs[, p + 1L] <- slope * rowSums(obs[, seq_len(p), drop = FALSE]) +
      error_sd * obs[, p + 1L]
    obs
  }

  runs <- .with_seed(seed, function() {
    vapply(seq_len(G), function(g) {
      .R2_sensitivity_one(draw, pilot_n, p, width, z, ci_method, max_n, batch)
    }, numeric(5L))
  })
  .R2_sensitivity_summary(runs, true_value, width, conf_level)
}
