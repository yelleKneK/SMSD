#' Pilot Sample Size for Two Stage Bounded Risk Estimation of a Normal Mean
#'
#' Returns the first stage (pilot) sample size for the two stage bounded
#' risk point estimation procedure for the mean of a normal distribution,
#' evaluated by \code{\link{ts_br_norm_mean}}. The loss in estimating the
#' mean \eqn{\mu} by the sample mean is \eqn{A\,|\bar{X} - \mu|^k}, and
#' the procedure keeps its expectation, the risk, at or below \code{w}.
#'
#' The pilot size is
#' \deqn{m = \max\{m_0,\ \lceil ((B / w)^{2/k})^{1/(1 + \gamma)} \rceil\},
#'   \qquad B = A \, \frac{2^{k/2}}{\sqrt{\pi}} \,
#'   \Gamma\!\left(\frac{k + 1}{2}\right),}
#' where \eqn{B = A\,E|Z|^k} for a standard normal \eqn{Z}, which is
#' \eqn{A} for \eqn{k = 2} and \eqn{A\sqrt{2/\pi}} for \eqn{k = 1}. (For
#' \eqn{k = 2} and the default \eqn{\gamma = 1} the rule is therefore
#' \eqn{\lceil \sqrt{A / w} \rceil}.) The sample size
#' that would be used if the variance \eqn{\sigma^2} were known is
#' \eqn{n^* = (B / w)^{2/k} \sigma^2}; the rule raises the constant
#' \eqn{(B / w)^{2/k}} to the power \eqn{1/(1 + \gamma)}, so the first
#' stage grows with the size of the problem but more slowly than \eqn{n^*}
#' does. This is the form of pilot that Chattopadhyay and Kelley (2016,
#' Equation 24) recommend, following Mukhopadhyay and de Silva (2009), for
#' a purely sequential procedure. Using it here, with \eqn{n^*} for
#' \eqn{\sigma = 1} as the constant, is this package's choice and not part
#' of the two stage procedure. A pilot equal to \eqn{(B / w)^{2/k}}
#' itself, the optimal sample size for \eqn{\sigma = 1}, would exceed
#' \eqn{n^*} whenever \eqn{\sigma} is below 1. The rule's first stage
#' exceeds \eqn{n^*} only when \eqn{\sigma^2} is below
#' \eqn{((B / w)^{2/k})^{-\gamma/(1 + \gamma)}}, a threshold that shrinks
#' as the problem grows.
#'
#' The minimum \code{m_0} matters more than the rule. The procedure is
#' defined only for a first stage of at least \eqn{k + 2} observations,
#' and it pays for not knowing the variance through a constant that is
#' large when the first stage is small: the expected total sample size is
#' about \eqn{b_m^{2/k}} times \eqn{n^*} (see \code{\link{ts_br_norm_mean}}
#' for \eqn{b_m}). For squared error (\eqn{k = 2}) that ratio is
#' \eqn{(m - 1)/(m - 3)}: 3 for \eqn{m = 4}, 2 for \eqn{m = 5}, 1.29 for
#' \eqn{m = 10}, and 1.12 for \eqn{m = 20}. For absolute error
#' (\eqn{k = 1}) it is 3.14 for \eqn{m = 3}, 1.57 for \eqn{m = 5}, 1.20
#' for \eqn{m = 10}, and 1.09 for \eqn{m = 20}. The default
#' \code{m_0 = 10} holds the expected total to about 1.3 times
#' \eqn{n^*}. A first stage larger than \eqn{n^*} is wasted, so lower
#' \code{m_0} only when even ten observations may be more than the study
#' needs.
#'
#' @param A The loss constant: the loss is \eqn{A\,|\bar{X} - \mu|^k}. A
#'   single positive number, in the units of \code{w} per unit of error to
#'   the power \code{k}.
#' @param k The error index: 1 for absolute error loss, 2 for squared
#'   error loss.
#' @param w The risk bound: the largest acceptable expected loss. A single
#'   positive number.
#' @param m_0 The smallest first stage sample size to return. A whole
#'   number of at least \eqn{k + 2}. Default 10; see Details.
#' @param gamma Positive tuning constant of the pilot rule. Larger values
#'   give a smaller first stage. Default 1.
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
#' # Mean waiting time, squared error loss, mean squared error at most 4.
#' ts_br_norm_mean_pilot(A = 1, k = 2, w = 4)
#'
#' # A demanding bound: the rule exceeds the minimum.
#' ts_br_norm_mean_pilot(A = 1, k = 2, w = 0.001)
#'
#' @seealso \code{\link{ts_br_norm_mean}}
#' @family ts_br_norm_mean
#' @export
ts_br_norm_mean_pilot <- function(A, k, w, m_0 = 10, gamma = 1) {
  .br_check_Akw(A, k, w)
  rule <- .ts_check_pilot_args(m_0, gamma, min_m_0 = as.integer(k) + 2L)
  .ts_pilot_size((.br_constant(A, k) / w)^(2 / k), rule$m_0, rule$gamma,
                 advice = .ts_advice_risk())
}


#' Two Stage Bounded Risk Point Estimation of a Normal Mean
#'
#' Evaluates the two stage bounded risk point estimation procedure for the
#' mean \eqn{\mu} of a normal distribution whose variance is unknown. The
#' loss in estimating \eqn{\mu} by the sample mean is
#' \eqn{A\,|\bar{X} - \mu|^k} (absolute error for \eqn{k = 1}, squared
#' error for \eqn{k = 2}), and the goal is a sample size at which the
#' expected loss, the risk, is at most \code{w}. The first stage
#' (\code{pilot_data}) fixes the total sample size; the second stage
#' (\code{additional_data}) completes it. The estimate is the mean of all
#' observations from both stages. The procedure is risk based, so there is
#' no confidence level.
#'
#' If the standard deviation \eqn{\sigma} were known, the risk of the mean
#' of \eqn{n} observations would be \eqn{B\,(\sigma^2 / n)^{k/2}} with
#' \eqn{B = A\,2^{k/2}\,\Gamma((k + 1)/2)/\sqrt{\pi}} (which is \eqn{A}
#' for \eqn{k = 2} and \eqn{A\sqrt{2/\pi}} for \eqn{k = 1}), and the
#' smallest sample size meeting the bound would be
#' \eqn{n^* = (B / w)^{2/k} \sigma^2}. With \eqn{m} first stage
#' observations and their sample variance \eqn{S_m^2}, the total sample
#' size is
#' \deqn{N = \max\{m,\ \lceil (b_m B / w)^{2/k} \, S_m^2 \rceil\},
#'   \qquad b_m = \left(\frac{m - 1}{2}\right)^{k/2}
#'   \frac{\Gamma((m - k - 1)/2)}{\Gamma((m - 1)/2)}.}
#' If \eqn{N = m} the first stage is enough; otherwise \eqn{N - m} further
#' observations are collected, once. \eqn{N} depends on the first stage
#' alone. The row \code{optimal_n} is this total \eqn{N}, not the number
#' still to collect.
#'
#' The constant \eqn{b_m} is the expectation of
#' \eqn{(S_m^2 / \sigma^2)^{-k/2}}, which exists only when
#' \eqn{m \ge k + 2}; for \eqn{k = 2} it is \eqn{(m - 1)/(m - 3)}. A
#' smaller first stage is an error.
#'
#' \strong{The guarantee is exact.} For normal data the risk of the final
#' mean is at most \code{w} for every \eqn{\mu}, every \eqn{\sigma}, and
#' every first stage size \eqn{m \ge k + 2}: \eqn{N} depends only on
#' \eqn{S_m^2}, which is independent of the sample means, so the risk is
#' \eqn{B \sigma^k E[N^{-k/2}]}, and
#' \eqn{N \ge (b_m B / w)^{2/k} S_m^2} bounds that by \code{w}. The price
#' is sample size: the expected \eqn{N} is about \eqn{b_m^{2/k}} times
#' \eqn{n^*} (for \eqn{k = 2}, \eqn{(m - 1)/(m - 3)}), which is why
#' \code{\link{ts_br_norm_mean_pilot}} does not return a very small first
#' stage. The guarantee rests on normality.
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
#' Data on a scale so large or so small that a variance cannot be
#' represented in double precision (observations beyond about
#' \eqn{10^{154}}, or differing by less than about \eqn{10^{-154}}) are an
#' error that says to express them in other units; no row of the result is
#' ever infinite.
#'
#' @param pilot_data Numeric vector of the first stage observations, at
#'   least \eqn{k + 2}. The whole first stage must be supplied at every
#'   call, because it alone determines the total sample size.
#' @param additional_data Numeric vector of the second stage observations
#'   collected so far, or \code{NULL} (the default) when none have been
#'   collected.
#' @param A The loss constant: the loss is \eqn{A\,|\bar{X} - \mu|^k}. A
#'   single positive number, in the units of \code{w} per unit of error to
#'   the power \code{k}.
#' @param k The error index: 1 for absolute error loss, 2 for squared
#'   error loss.
#' @param w The risk bound: the largest acceptable expected loss. A single
#'   positive number.
#' @param na.rm Logical. If \code{TRUE} (the default), missing values are
#'   removed from both stages before anything is computed; if
#'   \code{FALSE}, a missing value is an error.
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} (see
#'   \code{\link{smsd_check_methods}} for its \code{print()},
#'   \code{summary()}, and \code{tidy()} methods) with columns \code{term}
#'   and \code{value} and four rows:
#'   \describe{
#'     \item{\code{mean}}{The mean of all observations supplied, the
#'       estimate of \eqn{\mu}.}
#'     \item{\code{risk}}{The estimated risk of that mean at the current
#'       sample size, \eqn{B\,(S^2 / n)^{k/2}}, with \eqn{S^2} the sample
#'       variance of all observations supplied. A descriptive estimate;
#'       the guarantee is the one stated in Details.}
#'     \item{\code{optimal_n}}{The total sample size \eqn{N} fixed by the
#'       first stage.}
#'     \item{\code{current_n}}{The number of observations supplied, both
#'       stages together.}
#'   }
#'   Attributes: \code{is_satisfied} (\code{TRUE} once \code{current_n} is at
#'   least \code{optimal_n}), \code{current_n}, \code{criterion_n} (equal to
#'   \code{optimal_n}), \code{loss_A}, \code{k}, \code{risk_bound},
#'   \code{method}, and \code{notes} (always empty for this function). There
#'   is no \code{conf_level} or \code{width}, and no \code{variance_unusable}
#'   attribute: a first stage whose variance cannot be estimated is an error
#'   (see Details).
#'
#' @references
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#'   Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # A clinic manager wants the mean waiting time (minutes) estimated with
#' # a mean squared error of at most 4, that is, squared error loss with
#' # A = 1 and risk bound w = 4.
#' set.seed(113)
#' m <- ts_br_norm_mean_pilot(A = 1, k = 2, w = 4)
#' first_stage <- rnorm(m, mean = 45, sd = 12)
#' check <- ts_br_norm_mean(first_stage, A = 1, k = 2, w = 4)
#' check
#'
#' # The first stage fixed the total sample size; collect the remainder.
#' n_more <- attr(check, "criterion_n") - attr(check, "current_n")
#' second_stage <- rnorm(n_more, mean = 45, sd = 12)
#' ts_br_norm_mean(first_stage, second_stage, A = 1, k = 2, w = 4)
#'
#' @seealso \code{\link{ts_br_norm_mean_pilot}} for the first stage size;
#'   \code{\link{seq_br_norm_mean}} for the purely sequential procedure.
#' @family ts_br_norm_mean
#' @export
ts_br_norm_mean <- function(pilot_data, additional_data = NULL, A, k, w,
                            na.rm = TRUE) {
  .br_check_Akw(A, k, w)
  # The stopping constant is the expectation of the first stage standard
  # deviation to the power -k, which exists only for m >= k + 2.
  cc <- .ts_combine(
    pilot_data, additional_data, na.rm, min_n = as.integer(k) + 2L,
    advice = .ts_advice_first_stage(
      paste0("With k = ", k, " the stopping constant of the two stage ",
             "bounded risk procedure (the expectation of the first stage ",
             "standard deviation to the power -k) exists only for a first ",
             "stage of at least k + 2 = ", k + 2, " observations."),
      "ts_br_norm_mean_pilot"))
  .ts_stop_if_constant(cc$pilot, "variance")
  # The variances must be numbers the rule can use: observations whose
  # squares overflow or underflow are an error that says to change the
  # units, not a message about A and w or a risk of Inf or zero.
  s2_pilot <- .ts_usable_square(stats::var(cc$pilot),
                                "first stage variance",
                                "the observations differ")

  B <- .br_constant(A, k)
  b <- .ts_br_norm_b(cc$m, k)
  # N = (b B S_m^k / w)^(2 / k). For k = 2 that is (b B / w) S_m^2 as it
  # stands. For k = 1 the risk bound is in the units of the data, and the
  # ratio b B S_m / w is formed before it is squared: the square of b B / w
  # alone overflows when the data and w are both on a very small scale,
  # although the ratio, and N, are ordinary numbers.
  ratio <- (b * B) / w
  criterion <- if (k == 2) ratio * s2_pilot else (ratio * sqrt(s2_pilot))^2
  optimal_n <- max(cc$m, .as_sample_size(
    criterion, what = "total sample size", advice = .ts_advice_risk()))
  # Estimated risk of the mean at the current n: B * (S^2 / n)^(k/2).
  s2_all <- .ts_usable_square(stats::var(cc$all),
                              "variance of the observations supplied",
                              "the observations differ")
  risk <- B * (s2_all / cc$n)^(k / 2)

  term <- c("mean", "risk", "optimal_n", "current_n")
  value <- c(mean(cc$all), risk, optimal_n, cc$n)
  .ts_finite_rows(term, value, "restate the loss constant 'A' for those units.")
  new_smsd_check(
    term = term,
    value = value,
    method = "two stage bounded risk estimation (normal mean)",
    is_satisfied = cc$n >= optimal_n,
    current_n = cc$n,
    criterion_n = optimal_n,
    extra = list(loss_A = A, k = k, risk_bound = w)
  )
}
