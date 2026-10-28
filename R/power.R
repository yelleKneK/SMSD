#' Type II Error and Sample Size for a One Sample, One Sided z Test
#'
#' For the one sample, upper tailed \emph{z} test of
#' \eqn{H_0\!: \mu = \mu_0} against \eqn{H_1\!: \mu = \mu_1 > \mu_0} with a
#' known population standard deviation, returns either the sample size
#' needed for a chosen Type II error rate, the Type II error rate at a
#' chosen sample size, or a table of Type II error rates over a range of
#' sample sizes. Power is one minus the Type II error rate. The function
#' is named \code{power_z} so as not to mask \code{stats::power}.
#'
#' The test rejects \eqn{H_0} when
#' \eqn{\sqrt{n}\,(\bar{X} - \mu_0)/\sigma > z_{1-\alpha}}, with one-sided
#' Type I error rate \eqn{\alpha = 1 - } \code{conf_level}. Its Type II
#' error rate at \eqn{\mu = \mu_1} with \eqn{n} observations is
#' \deqn{\beta(n) = \Phi\!\left(z_{1-\alpha} -
#'   \sqrt{n}\,\frac{\mu_1 - \mu_0}{\sigma}\right),}
#' and the sample size at which it equals a target \eqn{\beta} is
#' \deqn{n_\beta = \left(\frac{(z_{1-\alpha} + z_{1-\beta})\,\sigma}
#'   {\mu_1 - \mu_0}\right)^2.}
#' \eqn{n_\beta} is rarely a whole number. \code{method = "sample_size"}
#' returns the smallest whole sample size whose Type II error rate
#' \eqn{\beta(n)} is at most \code{beta}, which is \eqn{n_\beta} rounded
#' up. (When \eqn{n_\beta} is within rounding error of a whole number,
#' \eqn{\beta(n)} decides between that number and the next, so asking for
#' the Type II error rate that this function reports at \eqn{n} returns
#' \eqn{n}.) When \code{beta} is at least \code{conf_level} (a power no
#' greater than the Type I error rate) every sample size meets the target
#' and 1 is returned.
#'
#' Only the upper tailed alternative is computed, so \code{mu_1} must
#' exceed \code{mu_0}. For a lower tailed alternative the results are the
#' same with both means reflected: pass \code{-mu_0} and \code{-mu_1}. For
#' an unknown standard deviation, two groups, or a two-sided test, use the
#' planning functions in \pkg{DMAR}, for example
#' \code{\link[DMAR]{ss_power_smd}}.
#'
#' @param mu_0 The population mean under the null hypothesis. A single
#'   number.
#' @param mu_1 The population mean under the alternative hypothesis. A
#'   single number greater than \code{mu_0}.
#' @param sigma The known population standard deviation. A single positive
#'   number.
#' @param conf_level One minus the one-sided Type I error rate, strictly
#'   between 0 and 1. Default 0.95, a test at the 0.05 level.
#' @param n For \code{method = "type_2"}: the sample size at which to
#'   compute the Type II error rate. A whole number of at least 1.
#' @param n_min,n_max For \code{method = "type_2_seq"}: the smallest and
#'   largest sample sizes of the table. Whole numbers with
#'   \code{n_max >= n_min >= 1}.
#' @param n_values For \code{method = "type_2_seq"}: the number of sample
#'   sizes in the table. They are equally spaced from \code{n_min} to
#'   \code{n_max}, and each is taken to the nearest whole number (a value
#'   halfway between two whole numbers goes to the larger); duplicates are
#'   dropped, so the table can have fewer rows than \code{n_values}. When
#'   \code{n_values} is at least the number of whole sample sizes from
#'   \code{n_min} to \code{n_max}, the table has every one of them.
#' @param beta For \code{method = "sample_size"}: the target Type II error
#'   rate (one minus the desired power), strictly between 0 and 1.
#' @param method What to compute: \code{"sample_size"}, which needs
#'   \code{beta}; \code{"type_2"}, which needs \code{n}; or
#'   \code{"type_2_seq"}, which needs \code{n_min}, \code{n_max}, and
#'   \code{n_values}. When \code{method} is not given, it is the method
#'   whose arguments are supplied (\code{beta} alone, \code{n} alone, or
#'   the three grid arguments); supplying arguments for more than one result
#'   without naming \code{method} is an error. When \code{method} is
#'   given, arguments it does not use are ignored. A unique abbreviation
#'   of a method name is accepted; \code{NULL} is the same as not giving
#'   \code{method}.
#'
#' @return Depends on \code{method}:
#'   \describe{
#'     \item{\code{"sample_size"}}{A length one integer: the smallest
#'       sample size whose Type II error rate is at most \code{beta}
#'       (\eqn{n_\beta} of Details, rounded up).}
#'     \item{\code{"type_2"}}{A length one number: the Type II error rate
#'       \eqn{\beta(n)}.}
#'     \item{\code{"type_2_seq"}}{A \code{data.frame} with one row per
#'       sample size and two columns: \code{n} (integer, the sample size)
#'       and \code{type_2_error} (the Type II error rate \eqn{\beta(n)}).}
#'   }
#'
#' @references
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#' Francis Bilson Darku \email{fbilsond@@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # A training program is expected to raise the mean score on a
#' # standardized assessment from 100 to 104; the standard deviation of
#' # the assessment is known to be 15.
#'
#' # Sample size for a Type II error rate of 0.20 (power 0.80)
#' power_z(mu_0 = 100, mu_1 = 104, sigma = 15, beta = 0.20)
#'
#' # Type II error rate with 60 trainees
#' power_z(mu_0 = 100, mu_1 = 104, sigma = 15, n = 60)
#'
#' # Type II error rates for 20 to 140 trainees
#' power_z(mu_0 = 100, mu_1 = 104, sigma = 15, n_min = 20, n_max = 140,
#'         n_values = 7)
#'
#' @seealso \code{\link[DMAR]{ss_power_smd}} in \pkg{DMAR} for two groups
#'   with an unknown standard deviation.
#' @export
power_z <- function(mu_0, mu_1, sigma, conf_level = 0.95, n = NULL,
                    n_min = NULL, n_max = NULL, n_values = NULL, beta = NULL,
                    method = c("sample_size", "type_2", "type_2_seq")) {
  if (missing(mu_0) || missing(mu_1) || missing(sigma)) {
    .smsd_stop("'mu_0', 'mu_1', and 'sigma' must all be supplied: the mean ",
               "under the null hypothesis, the mean under the alternative, ",
               "and the known population standard deviation.")
  }
  mu_0 <- .check_number(mu_0, "mu_0")
  mu_1 <- .check_number(mu_1, "mu_1")
  if (mu_1 <= mu_0) {
    .smsd_stop("'mu_1' must be greater than 'mu_0' (got mu_0 = ",
               format(mu_0), " and mu_1 = ", format(mu_1), "): power_z() ",
               "computes the upper tailed test. For a lower tailed ",
               "alternative pass -mu_0 and -mu_1; the results are the same.")
  }
  sigma <- .check_positive(sigma, "sigma")
  conf_level <- .check_conf_level(conf_level)
  method <- .power_z_method(method, missing(method), n, n_min, n_max,
                            n_values, beta)
  difference <- mu_1 - mu_0
  # Opposite signed finite means can overflow on subtraction even when
  # their difference in standard deviation units is representable.
  effect <- if (is.finite(difference)) difference / sigma else
    mu_1 / sigma - mu_0 / sigma

  if (method == "sample_size") {
    if (is.null(beta)) {
      .smsd_stop("method = \"sample_size\" needs 'beta', the target Type II ",
                 "error rate (one minus the desired power). For the Type II ",
                 "error rate at a given sample size, supply 'n' and use ",
                 "method = \"type_2\".")
    }
    beta <- .check_number(beta, "beta", lower = 0, upper = 1,
                          lower_open = TRUE, upper_open = TRUE)
    return(.power_z_sample_size(effect, conf_level, beta))
  }
  if (method == "type_2") {
    if (is.null(n)) {
      .smsd_stop("method = \"type_2\" needs 'n', the sample size. For the ",
                 "sample size that achieves a given Type II error rate, ",
                 "supply 'beta' and use method = \"sample_size\".")
    }
    n <- .check_count(n, "n", min = 1L)
    return(.power_z_type_2(effect, conf_level, n))
  }
  # method == "type_2_seq"
  if (is.null(n_min) || is.null(n_max) || is.null(n_values)) {
    .smsd_stop("method = \"type_2_seq\" needs 'n_min', 'n_max', and ",
               "'n_values': the smallest and largest sample sizes and how ",
               "many sample sizes to tabulate.")
  }
  n_min <- .check_count(n_min, "n_min", min = 1L)
  n_max <- .check_count(n_max, "n_max", min = n_min)
  n_values <- .check_count(n_values, "n_values", min = 1L)
  n_grid <- .power_z_grid(n_min, n_max, n_values)
  data.frame(n = n_grid,
             type_2_error = .power_z_type_2(effect, conf_level, n_grid))
}

# The three results power_z() can return.
#' @noRd
.power_z_methods <- c("sample_size", "type_2", "type_2_seq")

# The sample sizes of the table: `n_values` points equally spaced from n_min
# to n_max, each taken to the nearest whole number with a half going up, and
# duplicates dropped. This builds a grid of whole sample sizes; it is not
# display rounding. Point i of 0, ..., n_values - 1 is
# n_min + i * (n_max - n_min) / (n_values - 1), and its nearest whole number
# is found by whole number division, so that a half is decided by the rule
# and never by rounding error (and not by round(), which sends a half to the
# even neighbor).
#
# The table cannot have more rows than there are whole sample sizes from
# n_min to n_max. When n_values is at least that many, the points are no
# more than 1 apart, so their nearest whole numbers are every sample size in
# the range, which is also what exactly that many points give. n_values is
# therefore capped there: the result is the same, and the work is in
# proportion to the table and not to n_values (which may be as large as
# .Machine$integer.max, where the uncapped grid alone would be 2^31 - 1
# doubles of 8 bytes each, 16 GiB).
#' @noRd
.power_z_grid <- function(n_min, n_max, n_values) {
  n_values <- min(as.numeric(n_values),
                  as.numeric(n_max) - as.numeric(n_min) + 1)
  steps <- n_values - 1
  if (steps == 0) {
    return(n_min)
  }
  i <- seq(from = 0, to = steps)
  span <- as.numeric(n_max) - as.numeric(n_min)
  unique(as.integer(n_min + (2 * i * span + steps) %/% (2 * steps)))
}

# Settle the method. A `method` that is not given (missing, NULL, or the
# whole vector of choices that is its default, as a wrapper passing its own
# default along would supply) is chosen from the arguments supplied.
# Otherwise it must be a single string that names, or uniquely abbreviates,
# one of the three choices, which the shared .check_choice() settles.
#' @noRd
.power_z_method <- function(method, not_given, n, n_min, n_max, n_values,
                            beta) {
  if (!(not_given || is.null(method) ||
          identical(method, .power_z_methods))) {
    return(.check_choice(method, .power_z_methods, "method"))
  }
  given <- c(sample_size = !is.null(beta),
             type_2 = !is.null(n),
             type_2_seq = !is.null(n_min) || !is.null(n_max) ||
               !is.null(n_values))
  if (sum(given) == 1L) {
    return(names(given)[given])
  }
  what <- paste0("'beta' (for the sample size that achieves that Type II ",
                 "error rate), 'n' (for the Type II error rate at that ",
                 "sample size), or 'n_min', 'n_max', and 'n_values' (for a ",
                 "table of Type II error rates)")
  if (sum(given) == 0L) {
    .smsd_stop("Nothing to compute. Supply ", what, ".")
  }
  .smsd_stop("More than one result was requested. Supply only one of ", what,
             ", or name the result wanted with 'method'.")
}

# Type II error rate of the upper tailed z test at standardized effect
# (mu_1 - mu_0) / sigma, one-sided level 1 - conf_level, and sample size n.
# The critical value z_{1 - alpha} is the conf_level quantile.
#' @noRd
.power_z_type_2 <- function(effect, conf_level, n) {
  stats::pnorm(stats::qnorm(conf_level) - sqrt(n) * effect)
}

# Smallest whole sample size whose Type II error rate is at most beta.
#
# The closed form is ((z_{1 - alpha} + z_{1 - beta}) / effect)^2, rounded
# up. z_{1 - beta} is taken from the upper tail so that a very small beta
# keeps its quantile (1 - beta would round to 1). When the sum of quantiles
# is not positive (beta >= conf_level) every sample size meets the target,
# and squaring the negative sum would return a spurious positive number, so
# 1 is returned.
#
# The closed form and the Type II error rate it inverts can disagree by one
# when the unrounded value is within rounding error of a whole number. The
# Type II error rate decides, so the result is the smallest n for which
# .power_z_type_2() itself is at most beta.
#' @noRd
.power_z_sample_size <- function(effect, conf_level, beta) {
  z_sum <- stats::qnorm(conf_level) + stats::qnorm(beta, lower.tail = FALSE)
  if (z_sum <= 0) {
    return(1L)
  }
  advice <- paste0("Check 'mu_0', 'mu_1', and 'sigma'; a difference this ",
                   "small next to the standard deviation calls for more ",
                   "than two billion observations.")
  n <- max(1L, .as_sample_size((z_sum / effect)^2, advice = advice))
  if (n > 1L && .power_z_type_2(effect, conf_level, n - 1L) <= beta) {
    n <- n - 1L
  } else if (.power_z_type_2(effect, conf_level, n) > beta) {
    n <- .as_sample_size(n + 1, advice = advice)
  }
  n
}
