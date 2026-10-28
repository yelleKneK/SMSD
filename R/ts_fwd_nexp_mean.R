#' Pilot Sample Size for a Two Stage Fixed Width Interval for a Negative
#' Exponential Location
#'
#' Returns the first stage (pilot) sample size for the two stage fixed
#' width confidence interval for the location (the threshold, or minimum
#' possible value) of a negative exponential distribution, evaluated by
#' \code{\link{ts_fwd_nexp_mean}}. The interval is
#' \eqn{[X_{(1)} - d,\ X_{(1)}]}, with \eqn{X_{(1)}} the smallest
#' observation: \code{d} is the full width of the interval.
#'
#' The pilot size is
#' \deqn{m = \max\{m_0,\ \lceil (a / d)^{1/(1 + \gamma)} \rceil\},
#'   \qquad a = \log(1/\alpha),}
#' with \eqn{\alpha = 1 - } \code{conf_level}. The sample size that would
#' be used if the scale \eqn{\sigma} were known is \eqn{C = a \sigma / d};
#' the rule raises the constant \eqn{a / d} to the power
#' \eqn{1/(1 + \gamma)}, so the first stage grows as the target narrows but
#' more slowly than \eqn{C} does. This is the form of pilot that
#' Chattopadhyay and Kelley (2016, Equation 24) recommend, following
#' Mukhopadhyay and de Silva (2009), for a purely sequential procedure.
#' Using it here, with the known scale sample size for \eqn{\sigma = 1} as
#' the constant, is this package's choice and not part of the two stage
#' procedure, which is valid for any first stage of at least two
#' observations.
#'
#' The minimum \code{m_0} matters more than the rule. The procedure pays
#' for not knowing the scale through an \emph{F} quantile on \eqn{2} and
#' \eqn{2m - 2} degrees of freedom, so the expected total sample size is
#' about \eqn{F_{2,2m-2,\alpha} / a} times \eqn{C}. At a confidence level
#' of 0.95 that ratio is 6.34 for \eqn{m = 2}, 2.32 for \eqn{m = 3}, 1.49
#' for \eqn{m = 5}, 1.19 for \eqn{m = 10}, and 1.08 for \eqn{m = 20}. The
#' default \code{m_0 = 10} holds the expected total to about 1.2 times
#' \eqn{C} at 0.95 (about 1.3 times \eqn{C} at 0.99). A first stage
#' larger than \eqn{C} is wasted, so lower \code{m_0} only when even ten
#' observations may be more than the study needs.
#'
#' @param d The full width of the confidence interval, in the units of the
#'   data. A single positive number.
#' @param conf_level The confidence level, strictly between 0 and 1.
#'   Default 0.95.
#' @param m_0 The smallest first stage sample size to return. A whole
#'   number of at least 2 (the scale estimate needs two observations).
#'   Default 10; see Details.
#' @param gamma Positive tuning constant of the pilot rule. Larger values
#'   give a smaller first stage. Default 1, for which the rule is
#'   \eqn{\lceil \sqrt{a / d} \rceil}.
#'
#' @return A length one integer: the first stage sample size.
#'
#' @references
#' Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient
#' of variation with minimum risk: A sequential method for minimizing
#' sampling error and study cost. \emph{Multivariate Behavioral Research},
#' \emph{51}(5), 627--648. \doi{10.1080/00273171.2016.1203279}
#'
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#'   Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # Shortest possible task completion time, to within half a second.
#' ts_fwd_nexp_mean_pilot(d = 0.5, conf_level = 0.95)
#'
#' # A narrow target: the rule exceeds the minimum.
#' ts_fwd_nexp_mean_pilot(d = 0.01, conf_level = 0.99)
#'
#' @seealso \code{\link{ts_fwd_nexp_mean}}
#' @family ts_fwd_nexp_mean
#' @export
ts_fwd_nexp_mean_pilot <- function(d, conf_level = 0.95, m_0 = 10, gamma = 1) {
  .ts_check_width(d, what = "the confidence interval full width")
  alpha <- .ts_alpha(conf_level)
  rule <- .ts_check_pilot_args(m_0, gamma, min_m_0 = 2L)
  # The criterion is linear in 1 / d with constant log(1 / alpha), so the
  # constant on the sample size scale is log(1 / alpha) / d.
  .ts_pilot_size(-log(alpha) / d, rule$m_0, rule$gamma,
                 advice = .ts_advice_width("full"))
}


#' Two Stage Fixed Width Interval for a Negative Exponential Location
#'
#' Evaluates the two stage procedure for a confidence interval of fixed
#' width for the location \eqn{\theta} of a negative exponential
#' distribution: the model \eqn{X = \theta + \sigma E} with \eqn{E} a
#' standard exponential variable, so that \eqn{\theta} is the smallest
#' value the outcome can take (a threshold, a guarantee time, a minimum
#' response time) and \eqn{\sigma} is the scale. Despite the function's
#' name, the quantity estimated is the location \eqn{\theta}, not the mean
#' \eqn{\theta + \sigma}. The first stage (\code{pilot_data}) fixes the
#' total sample size; the second stage (\code{additional_data}) completes
#' it. The interval is \eqn{[X_{(1)} - d,\ X_{(1)}]}, where \eqn{X_{(1)}}
#' is the smallest of all observations from both stages and \code{d} is
#' the full width.
#'
#' With \eqn{m} first stage observations, their minimum \eqn{X_{(1:m)}},
#' and the scale estimate
#' \eqn{U_m = \sum_{i=1}^{m} (X_i - X_{(1:m)}) / (m - 1)}, the total
#' sample size is
#' \deqn{N = \max\{m,\ \lceil F_{2,2m-2,\alpha} \, U_m / d \rceil\},}
#' where \eqn{F_{2,2m-2,\alpha}} is the upper \eqn{\alpha} point of the
#' \emph{F} distribution with 2 and \eqn{2m - 2} degrees of freedom and
#' \eqn{\alpha = 1 - } \code{conf_level}. If \eqn{N = m} the first stage is
#' enough; otherwise \eqn{N - m} further observations are collected, once.
#' \eqn{N} depends on the first stage alone. The row \code{optimal_n} is
#' this total \eqn{N}, not the number still to collect.
#'
#' The quantile is one-sided because the interval is: the smallest
#' observation can never fall below \eqn{\theta}, so the interval misses
#' only when \eqn{X_{(1)} - \theta > d}.
#'
#' \strong{The guarantee is exact.} For negative exponential data the
#' interval covers \eqn{\theta} with probability at least
#' \code{conf_level} for every \eqn{\theta}, every \eqn{\sigma}, every
#' \code{d}, and every first stage size \eqn{m \ge 2}. (Given \eqn{N}, the
#' interval misses with probability \eqn{\exp(-N d / \sigma)}; \eqn{N}
#' depends only on \eqn{U_m}, with \eqn{2(m - 1) U_m / \sigma} distributed
#' as chi square on \eqn{2m - 2} degrees of freedom; and the expectation
#' of \eqn{\exp(-F U_m / \sigma)} equals \eqn{\alpha} exactly at the
#' \emph{F} quantile used.) The price is sample size: the expected \eqn{N}
#' is about \eqn{F_{2,2m-2,\alpha} / \log(1/\alpha)} times the known scale
#' sample size \eqn{C = \log(1/\alpha) \sigma / d}, which is why
#' \code{\link{ts_fwd_nexp_mean_pilot}} does not return a very small first
#' stage. The guarantee rests on the negative exponential model.
#'
#' The interval rows are reported at every call so that progress can be
#' followed, but the interval carries the stated confidence level only
#' when the stopping rule is met, that is, when the current sample size
#' has reached \code{optimal_n}.
#'
#' \strong{A first stage with no variation is an error.} If every first
#' stage observation is equal, or equal up to rounding error (a range of no
#' more than \eqn{10^{-10}} times the largest absolute value, the rule used
#' throughout the package), the scale cannot be estimated and the function
#' stops with an error that says what to do. The purely sequential
#' functions report such data as a stopping rule that is not yet met, with
#' the attribute \code{variance_unusable}; here that would mislead, because
#' the total sample size is fixed by the first stage alone, so no later
#' observation could turn a rule reported as not met into one that is met.
#' Data whose real variation is that small beside their size (times
#' recorded in seconds since 1970 that differ by hundredths of a second,
#' for example) should have a common value subtracted from both stages
#' first; the total sample size does not depend on it.
#'
#' @param pilot_data Numeric vector of the first stage observations, at
#'   least two. The whole first stage must be supplied at every call,
#'   because it alone determines the total sample size.
#' @param additional_data Numeric vector of the second stage observations
#'   collected so far, or \code{NULL} (the default) when none have been
#'   collected.
#' @param d The full width of the confidence interval, in the units of the
#'   data. A single positive number.
#' @param conf_level The confidence level, strictly between 0 and 1.
#'   Default 0.95.
#' @param na.rm Logical. If \code{TRUE} (the default), missing values are
#'   removed from both stages before anything is computed; if
#'   \code{FALSE}, a missing value is an error.
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} (see
#'   \code{\link{smsd_check_methods}} for its \code{print()},
#'   \code{summary()}, and \code{tidy()} methods) with columns \code{term}
#'   and \code{value} and five rows:
#'   \describe{
#'     \item{\code{lower_limit}}{\eqn{X_{(1)} - d}, with \eqn{X_{(1)}} the
#'       smallest of all observations supplied.}
#'     \item{\code{location}}{\eqn{X_{(1)}}, the estimate of
#'       \eqn{\theta}.}
#'     \item{\code{upper_limit}}{\eqn{X_{(1)}}; the interval is one-sided,
#'       so its upper limit is the estimate itself.}
#'     \item{\code{optimal_n}}{The total sample size \eqn{N} fixed by the
#'       first stage.}
#'     \item{\code{current_n}}{The number of observations supplied, both
#'       stages together.}
#'   }
#'   Attributes: \code{is_satisfied} (\code{TRUE} once \code{current_n} is at
#'   least \code{optimal_n}), \code{current_n}, \code{criterion_n} (equal to
#'   \code{optimal_n}), \code{conf_level}, \code{width} (the full width
#'   \code{d}), \code{method}, and \code{notes} (always empty for this
#'   function). There is no \code{variance_unusable} attribute: a first stage
#'   whose scale cannot be estimated is an error (see Details).
#'
#' @references
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#'   Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # A usability researcher models the time to complete a task as a
#' # shortest possible time plus an exponential delay, and wants the
#' # shortest possible time to within half a second, with 95% confidence.
#' set.seed(113)
#' m <- ts_fwd_nexp_mean_pilot(d = 0.5, conf_level = 0.95)
#' first_stage <- 3 + rexp(m, rate = 1 / 4)
#' check <- ts_fwd_nexp_mean(first_stage, d = 0.5)
#' check
#'
#' # The first stage fixed the total sample size; collect the remainder.
#' n_more <- attr(check, "criterion_n") - attr(check, "current_n")
#' second_stage <- 3 + rexp(n_more, rate = 1 / 4)
#' ts_fwd_nexp_mean(first_stage, second_stage, d = 0.5)
#'
#' @seealso \code{\link{ts_fwd_nexp_mean_pilot}} for the first stage size;
#'   \code{\link{seq_fwd_nexp_mean}} for the purely sequential procedure.
#' @family ts_fwd_nexp_mean
#' @export
ts_fwd_nexp_mean <- function(pilot_data, additional_data = NULL, d,
                             conf_level = 0.95, na.rm = TRUE) {
  .ts_check_width(d, what = "the confidence interval full width")
  alpha <- .ts_alpha(conf_level)
  cc <- .ts_combine(pilot_data, additional_data, na.rm,
                    advice = .ts_advice_first_stage(
                      "The scale cannot be estimated from fewer.",
                      "ts_fwd_nexp_mean_pilot"))
  .ts_stop_if_constant(cc$pilot, "scale estimate", unknown = "scale")
  loc <- min(cc$all)

  # The interval [min - d, min] is one-sided (the minimum never falls below
  # the location), so the quantile is the upper alpha point of F(2, 2m - 2).
  #
  # When N is too large to represent, the cause is d unless the same rule
  # with log(1 / alpha), the constant for a known scale, is representable;
  # then it is the F quantile of a very small first stage, and the message
  # says so (the advice is evaluated only when it is needed).
  #
  # The rule depends on the data only through U_m / d, the scale estimate in
  # units of d, which is computed as such (as the variance is in
  # ts_fwd_norm_mean()): the product of the quantile and a scale estimate
  # near the largest representable number would overflow although the
  # ratio, and N, are small.
  scale_in_d <- sum((cc$pilot - min(cc$pilot)) / d) / (cc$m - 1)
  f_quantile <- stats::qf(conf_level, 2, 2 * cc$m - 2)
  optimal_n <- max(cc$m, .as_sample_size(
    f_quantile * scale_in_d, what = "total sample size",
    advice = .ts_advice_width("full", first_stage = .ts_first_stage_cause(
      cc$m, -log(alpha) * scale_in_d,
      used = paste0("the F quantile on 2 and ", 2 * cc$m - 2,
                    " degrees of freedom is ", format(f_quantile, digits = 4)),
      known = paste0("the constant for a known scale, log(1 / alpha), is ",
                     format(-log(alpha), digits = 4), "; with that constant"),
      pilot_fun = "ts_fwd_nexp_mean_pilot"))))
  ci <- c(loc - d, loc)

  term <- c("lower_limit", "location", "upper_limit", "optimal_n",
            "current_n")
  value <- c(ci[1], loc, ci[2], optimal_n, cc$n)
  .ts_finite_rows(term, value, "state 'd' in the same units.")
  new_smsd_check(
    term = term,
    value = value,
    method = "two stage fixed width interval (negative exponential location)",
    is_satisfied = cc$n >= optimal_n,
    current_n = cc$n,
    criterion_n = optimal_n,
    conf_level = conf_level,
    width = d
  )
}
