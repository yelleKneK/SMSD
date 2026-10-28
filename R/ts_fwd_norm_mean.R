#' Pilot Sample Size for a Two Stage Fixed Width Interval for a Normal Mean
#'
#' Returns the first stage (pilot) sample size for Stein's two stage fixed
#' width confidence interval for the mean of a normal distribution,
#' evaluated by \code{\link{ts_fwd_norm_mean}}. The interval is
#' \eqn{\bar{X} \pm d}: \code{d} is the half width (the margin of error),
#' so the full width is \eqn{2d}.
#'
#' The pilot size is
#' \deqn{m = \max\{m_0,\ \lceil (z_{\alpha/2} / d)^{2/(1 + \gamma)} \rceil\},}
#' where \eqn{z_{\alpha/2}} is the standard normal quantile for the
#' confidence level. The sample size that would be used if the variance
#' \eqn{\sigma^2} were known is \eqn{C = z_{\alpha/2}^2 \sigma^2 / d^2}
#' (see \code{\link{seq_fwd_norm_mean_optimal}}); the rule raises the
#' constant \eqn{(z_{\alpha/2} / d)^2} to the power \eqn{1/(1 + \gamma)},
#' so the first stage grows as the target narrows but more slowly than
#' \eqn{C} does. This is the form of pilot that Chattopadhyay and Kelley
#' (2016, Equation 24) recommend, following Mukhopadhyay and de Silva
#' (2009), for a purely sequential procedure. Using it here, with the known
#' variance sample size for \eqn{\sigma = 1} as the constant, is this
#' package's choice and not part of Stein's procedure, which is valid for
#' any first stage of at least two observations.
#'
#' The minimum \code{m_0} matters more than the rule. Stein's procedure
#' pays for not knowing the variance through the \emph{t} quantile on
#' \eqn{m - 1} degrees of freedom, so the expected total sample size is
#' about \eqn{t_{m-1,\alpha/2}^2 / z_{\alpha/2}^2} times \eqn{C}. At a
#' confidence level of 0.95 that ratio is 42.0 for \eqn{m = 2}, 4.82 for
#' \eqn{m = 3}, 2.01 for \eqn{m = 5}, 1.33 for \eqn{m = 10}, 1.14 for
#' \eqn{m = 20}, and 1.09 for \eqn{m = 30}. The default \code{m_0 = 10}
#' holds the expected total to about 1.33 times \eqn{C} at 0.95 (about
#' 1.6 times \eqn{C} at 0.99). A first stage larger than \eqn{C} is
#' wasted, so lower \code{m_0} only when even ten observations may be
#' more than the study needs.
#'
#' @param d The half width of the confidence interval (the margin of
#'   error), in the units of the data. A single positive number.
#' @param conf_level The confidence level, strictly between 0 and 1.
#'   Default 0.95.
#' @param m_0 The smallest first stage sample size to return. A whole
#'   number of at least 2 (a variance needs two observations). Default 10;
#'   see Details.
#' @param gamma Positive tuning constant of the pilot rule. Larger values
#'   give a smaller first stage. Default 1, for which the rule is
#'   \eqn{\lceil z_{\alpha/2} / d \rceil}.
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
#' Stein, C. (1945). A two-sample test for a linear hypothesis whose power
#' is independent of the variance. \emph{The Annals of Mathematical
#' Statistics}, \emph{16}(3), 243--258. \doi{10.1214/aoms/1177731088}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#'   Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # Mean score on a 100 point skills assessment, to within 2 points.
#' ts_fwd_norm_mean_pilot(d = 2, conf_level = 0.95)
#'
#' # A narrow target on a standardized scale: the rule exceeds the minimum.
#' ts_fwd_norm_mean_pilot(d = 0.1, conf_level = 0.95)
#'
#' @seealso \code{\link{ts_fwd_norm_mean}};
#'   \code{\link{seq_fwd_norm_mean_optimal}} for the known variance sample
#'   size.
#' @family ts_fwd_norm_mean
#' @export
ts_fwd_norm_mean_pilot <- function(d, conf_level = 0.95, m_0 = 10, gamma = 1) {
  .ts_check_width(d)
  alpha <- .ts_alpha(conf_level)
  rule <- .ts_check_pilot_args(m_0, gamma, min_m_0 = 2L)
  .ts_pilot_size((.aipe_z(conf_level) / d)^2, rule$m_0, rule$gamma,
                 advice = .ts_advice_width("half"))
}


#' Two Stage Fixed Width Confidence Interval for a Normal Mean
#'
#' Evaluates Stein's two stage procedure for a confidence interval of
#' fixed width for the mean \eqn{\mu} of a normal distribution whose
#' variance is unknown. The first stage (\code{pilot_data}) fixes the
#' total sample size; the second stage (\code{additional_data}) completes
#' it. The interval is \eqn{\bar{X}_N \pm d}, where \code{d} is the half
#' width (the margin of error) and \eqn{\bar{X}_N} is the mean of all
#' observations from both stages.
#'
#' With \eqn{m} first stage observations and their sample variance
#' \eqn{S_m^2}, the total sample size is
#' \deqn{N = \max\{m,\ \lceil t_{m-1,\alpha/2}^2 \, S_m^2 / d^2 \rceil\},}
#' where \eqn{t_{m-1,\alpha/2}} is the upper \eqn{\alpha/2} point of the
#' \emph{t} distribution with \eqn{m - 1} degrees of freedom and
#' \eqn{\alpha = 1 - } \code{conf_level}. If \eqn{N = m} the first stage
#' is enough; otherwise \eqn{N - m} further observations are collected,
#' once. \eqn{N} depends on the first stage alone and does not change as
#' second stage observations arrive. The row \code{optimal_n} is this
#' total \eqn{N}, not the number still to collect.
#'
#' \strong{The guarantee is exact.} For normal data the interval
#' \eqn{\bar{X}_N \pm d} covers \eqn{\mu} with probability at least
#' \code{conf_level} for every \eqn{\mu}, every \eqn{\sigma}, every
#' \code{d}, and every first stage size \eqn{m \ge 2} (Stein, 1945;
#' Mukhopadhyay and de Silva, 2009, discuss the procedure in their Section
#' 6.2.1). The price is sample size: the expected \eqn{N} is about
#' \eqn{t_{m-1,\alpha/2}^2 / z_{\alpha/2}^2} times the known variance
#' sample size \eqn{C = z_{\alpha/2}^2 \sigma^2 / d^2}, which is why
#' \code{\link{ts_fwd_norm_mean_pilot}} does not return a very small first
#' stage. The guarantee rests on normality; for data that are not normal
#' see \code{\link{ts_fwd_dfr_mean}}, whose level is only asymptotic.
#'
#' The interval rows are reported at every call so that progress can be
#' followed, but the interval carries the stated confidence level only
#' when the stopping rule is met, that is, when the current sample size
#' has reached \code{optimal_n}.
#'
#' \strong{A first stage with no variation is an error.} If every first
#' stage observation is equal, or equal up to rounding error (a range of no
#' more than \eqn{10^{-10}} times the largest absolute value, the rule used
#' throughout the package), the first stage variance cannot be estimated
#' and the function stops with an error that says what to do. The purely
#' sequential functions report such data as a stopping rule that is not yet
#' met, with the attribute \code{variance_unusable}; here that would
#' mislead, because the total sample size is fixed by the first stage
#' alone, so no later observation could turn a rule reported as not met
#' into one that is met. Data whose real variation is that small beside
#' their size (times recorded in seconds since 1970 that differ by
#' hundredths of a second, for example) should have a common value
#' subtracted from both stages first; the total sample size does not depend
#' on it.
#'
#' @param pilot_data Numeric vector of the first stage observations, at
#'   least two. The whole first stage must be supplied at every call,
#'   because it alone determines the total sample size.
#' @param additional_data Numeric vector of the second stage observations
#'   collected so far, or \code{NULL} (the default) when none have been
#'   collected.
#' @param d The half width of the confidence interval (the margin of
#'   error), in the units of the data. A single positive number. The full
#'   width is \eqn{2d}.
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
#'     \item{\code{lower_limit}}{\eqn{\bar{X} - d}, with \eqn{\bar{X}} the
#'       mean of all observations supplied.}
#'     \item{\code{mean}}{\eqn{\bar{X}}, the estimate of \eqn{\mu}.}
#'     \item{\code{upper_limit}}{\eqn{\bar{X} + d}.}
#'     \item{\code{optimal_n}}{The total sample size \eqn{N} fixed by the
#'       first stage.}
#'     \item{\code{current_n}}{The number of observations supplied, both
#'       stages together.}
#'   }
#'   Attributes: \code{is_satisfied} (\code{TRUE} once \code{current_n} is at
#'   least \code{optimal_n}), \code{current_n}, \code{criterion_n} (equal to
#'   \code{optimal_n}), \code{conf_level}, \code{width} (the full width
#'   \eqn{2d}), \code{method}, and \code{notes} (always empty for this
#'   function). There is no \code{variance_unusable} attribute: a first stage
#'   whose variance cannot be estimated is an error (see Details).
#'
#' @references
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' Stein, C. (1945). A two-sample test for a linear hypothesis whose power
#' is independent of the variance. \emph{The Annals of Mathematical
#' Statistics}, \emph{16}(3), 243--258. \doi{10.1214/aoms/1177731088}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#'   Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # A training director wants the mean score on a 100 point skills
#' # assessment to within 2 points, with 95% confidence.
#' set.seed(113)
#' m <- ts_fwd_norm_mean_pilot(d = 2, conf_level = 0.95)
#' first_stage <- rnorm(m, mean = 70, sd = 8)
#' check <- ts_fwd_norm_mean(first_stage, d = 2)
#' check
#'
#' # The first stage fixed the total sample size; collect the remainder.
#' n_more <- attr(check, "criterion_n") - attr(check, "current_n")
#' second_stage <- rnorm(n_more, mean = 70, sd = 8)
#' ts_fwd_norm_mean(first_stage, second_stage, d = 2)
#'
#' @seealso \code{\link{ts_fwd_norm_mean_pilot}} for the first stage size;
#'   \code{\link{seq_fwd_norm_mean}} for the purely sequential procedure;
#'   \code{\link{seq_fwd_norm_mean_optimal}} for the sample size when the
#'   variance is known; \code{\link{ts_fwd_dfr_mean}} when normality is not
#'   assumed.
#' @family ts_fwd_norm_mean
#' @export
ts_fwd_norm_mean <- function(pilot_data, additional_data = NULL, d,
                             conf_level = 0.95, na.rm = TRUE) {
  .ts_check_width(d)
  alpha <- .ts_alpha(conf_level)
  cc <- .ts_combine(pilot_data, additional_data, na.rm,
                    advice = .ts_advice_first_stage(
                      "A variance cannot be computed from fewer.",
                      "ts_fwd_norm_mean_pilot"))
  .ts_stop_if_constant(cc$pilot, "variance")
  xbar <- mean(cc$all)

  # Stein's rule. The interval xbar +/- d is two-sided, so the quantile is
  # the upper alpha / 2 point of t on m - 1 degrees of freedom; the variance
  # is that of the first stage alone, which is what makes N independent of
  # the first stage mean and the coverage exact. S_m^2 / d^2 is computed as
  # the variance of the first stage in units of d: the rule depends on the
  # data only through that ratio, and in this form neither squares of very
  # large observations nor of a very small d overflow or underflow on the
  # way, so N is the same in any units of measurement.
  #
  # When N is too large to represent, the cause is d unless the same rule
  # with the normal quantile is representable; then it is the t quantile of
  # a very small first stage, and the message says so (the advice is
  # evaluated only when it is needed).
  t_quantile <- stats::qt(alpha / 2, cc$m - 1, lower.tail = FALSE)
  spread <- stats::var(cc$pilot / d)
  first_stage_cause <- function() {
    z <- .aipe_z(conf_level)
    .ts_first_stage_cause(
      cc$m, z^2 * spread,
      used = paste0("the t quantile on ", cc$m - 1,
                    if (cc$m == 2L) " degree" else " degrees",
                    " of freedom is ", format(t_quantile, digits = 4)),
      known = paste0("the normal quantile is ", format(z, digits = 4),
                     "; with the normal quantile"),
      pilot_fun = "ts_fwd_norm_mean_pilot")
  }
  optimal_n <- max(cc$m, .as_sample_size(
    t_quantile^2 * spread, what = "total sample size",
    advice = .ts_advice_width("half", first_stage = first_stage_cause())))
  ci <- xbar + c(-1, 1) * d

  term <- c("lower_limit", "mean", "upper_limit", "optimal_n", "current_n")
  value <- c(ci[1], xbar, ci[2], optimal_n, cc$n)
  .ts_finite_rows(term, value, "state 'd' in the same units.")
  new_smsd_check(
    term = term,
    value = value,
    method = "two stage fixed width interval (normal mean)",
    is_satisfied = cc$n >= optimal_n,
    current_n = cc$n,
    criterion_n = optimal_n,
    conf_level = conf_level,
    width = 2 * d
  )
}
