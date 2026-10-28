#' Pilot Sample Size for Two Stage Bounded Risk Estimation of an
#' Exponential Mean
#'
#' Returns the first stage (pilot) sample size for the two stage bounded
#' risk point estimation procedure for the mean of an exponential
#' distribution, evaluated by \code{\link{ts_br_exp_mean}}. The loss in
#' estimating the mean \eqn{\theta} by the sample mean is
#' \eqn{A\,(\bar{X} - \theta)^2}, and the procedure keeps its expectation,
#' the risk, at or below \code{w}.
#'
#' The source of the procedure (Mukhopadhyay and de Silva, 2009, Section
#' 8.3.1) requires only that the first stage have at least three
#' observations and gives no rule for choosing its size. The choice does
#' not affect the risk bound, but it decides how many observations the
#' procedure uses. Writing \eqn{n^* = A \theta^2 / w} for the sample size
#' at which the risk \eqn{A \theta^2 / n} of a fixed sample equals
#' \code{w}, the expected total sample size for a first stage of \eqn{m}
#' observations is, when \eqn{n^*} is well above \eqn{m},
#' \deqn{E[N] / n^* \approx 0.565 \times
#'   \frac{2 (m + 1)^2}{(m - 1)(m - 2)},}
#' where 0.565 is the constant of Zacks and Mukhopadhyay (2006; the book's
#' Equation 8.3.13) and the rest is the expectation of the stopping
#' variable (Equations 8.3.9 and 8.3.11).
#' \tabular{rrrrrrrr}{
#'   \eqn{m} \tab 3 \tab 5 \tab 10 \tab 15 \tab 20 \tab 30 \tab 50 \cr
#'   \eqn{E[N] / n^*} \tab 9.04 \tab 3.39 \tab 1.90 \tab 1.59 \tab 1.46
#'     \tab 1.34 \tab 1.25
#' }
#' The ratio falls toward 1.13 as \eqn{m} grows. (Exact calculation agrees
#' with the approximation to within 1 percent once \eqn{n^*} is at least
#' five times \eqn{m}, and to within 6 percent when \eqn{n^*} is as small
#' as \eqn{m}.) A first stage of three observations therefore uses about
#' nine times the optimal sample size. On the other side, a first stage
#' larger than \eqn{n^*} is itself the waste: the total can never be
#' smaller than \eqn{m}. The first stage should be as large as one can be
#' confident does not exceed \eqn{n^*}.
#'
#' The pilot size returned is
#' \deqn{m = \max\{m_0,\ \lceil (A / w)^{1/(1 + \gamma)} \rceil\}.}
#' With the default \code{m_0 = 20}, a first stage of 20 keeps the
#' expected total within about 1.5 times \eqn{n^*} whenever \eqn{n^*} is
#' at least 20 (1.54 at \eqn{n^* = 20}, 1.47 at 50). The second term
#' raises \eqn{A / w}, which is \eqn{n^*} for \eqn{\theta = 1}, to the
#' power \eqn{1/(1 + \gamma)}, so that the first stage grows with the size
#' of the problem but more slowly than \eqn{n^*} does. This is the form of
#' pilot that Chattopadhyay and Kelley (2016, Equation 24) recommend,
#' following Mukhopadhyay and de Silva (2009), for a purely sequential
#' procedure. Using it here is this package's choice; as said above, the
#' source of this two stage procedure gives no rule. If the mean is
#' expected to be so small that \eqn{n^* = A \theta^2 / w} may be below 20,
#' lower \code{m_0}.
#'
#' @param A The loss constant: the loss is \eqn{A\,(\bar{X} - \theta)^2}.
#'   A single positive number, in the units of \code{w} per squared unit
#'   of the data.
#' @param w The risk bound: the largest acceptable expected loss. A single
#'   positive number.
#' @param m_0 The smallest first stage sample size to return. A whole
#'   number of at least 3. Default 20; see Details.
#' @param gamma Positive tuning constant of the pilot rule. Larger values
#'   give a smaller first stage. Default 1, for which the rule is
#'   \eqn{\lceil \sqrt{A / w} \rceil}.
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
#' Zacks, S., & Mukhopadhyay, N. (2006). Bounded risk estimation of the
#' exponential parameter in a two-stage sampling. \emph{Sequential
#' Analysis}, \emph{25}(4), 437--452. \doi{10.1080/07474940600934896} (As
#' described in Mukhopadhyay and de Silva, 2009, Section 8.3.2, Equation
#' 8.3.13.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#'   Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # Mean time between service calls, squared error loss with A = 100,
#' # risk at most 4.
#' ts_br_exp_mean_pilot(A = 100, w = 4)
#'
#' # A demanding bound: the rule exceeds the minimum.
#' ts_br_exp_mean_pilot(A = 100, w = 0.01)
#'
#' # A smaller first stage, accepted knowingly (see the table in Details).
#' ts_br_exp_mean_pilot(A = 100, w = 4, m_0 = 10)
#'
#' @seealso \code{\link{ts_br_exp_mean}}
#' @family ts_br_exp_mean
#' @export
ts_br_exp_mean_pilot <- function(A, w, m_0 = 20, gamma = 1) {
  .ts_br_exp_check_Aw(A, w)
  rule <- .ts_check_pilot_args(m_0, gamma, min_m_0 = 3L)
  .ts_pilot_size(A / w, rule$m_0, rule$gamma, advice = .ts_advice_risk())
}

#' @noRd
.ts_br_exp_check_Aw <- function(A, w) {
  if (missing(A)) .smsd_stop("'A' must be supplied (the loss constant).")
  if (missing(w)) .smsd_stop("'w' must be supplied (the risk bound).")
  .check_positive(A, "A")
  .check_positive(w, "w")
  invisible(TRUE)
}


#' Two Stage Bounded Risk Point Estimation of an Exponential Mean
#'
#' Evaluates the two stage bounded risk point estimation procedure for the
#' mean \eqn{\theta} of an exponential distribution (a mean time to
#' failure, a mean waiting time, a mean time between events). The loss in
#' estimating \eqn{\theta} by the sample mean is
#' \eqn{A\,(\bar{X} - \theta)^2}, and the goal is a sample size at which
#' the expected loss, the risk, is at most \code{w} whatever the value of
#' \eqn{\theta}. The first stage (\code{pilot_data}) fixes the total
#' sample size; the second stage (\code{additional_data}) completes it.
#' The estimate is the mean of all observations from both stages. The
#' procedure is risk based, so there is no confidence level.
#'
#' The procedure is the two stage methodology that Mukhopadhyay and de
#' Silva (2009, Section 8.3.1) attribute to Mukhopadhyay and Pepe (2006),
#' with the improved constant that their Section 8.3.2 attributes to Zacks
#' and Mukhopadhyay (2006). With \eqn{m \ge 3} first stage observations
#' and their mean \eqn{\bar{X}_m}, the total sample size is (the book's
#' Equation 8.3.9)
#' \deqn{N = \max\{m,\ \lceil B \, \bar{X}_m^2 / w \rceil\}.}
#' If \eqn{N = m} the first stage is enough; otherwise \eqn{N - m} further
#' observations are collected, once, and \eqn{\theta} is estimated by the
#' mean \eqn{\bar{X}_N} of all \eqn{N} observations. \eqn{N} depends on the
#' first stage alone. The row \code{optimal_n} is this total \eqn{N}, not
#' the number still to collect. The constant is (Equations 8.3.11 and
#' 8.3.13)
#' \deqn{B = 0.565 \times \frac{2 m (m + 1) A}{(m - 1)(m - 2)}.}
#' (The book writes the rounding as one more than the largest integer
#' below \eqn{B \bar{X}_m^2 / w}, which is the ceiling.)
#'
#' \strong{What is guaranteed.} For the constant without the factor 0.565,
#' Theorem 8.3.1 of Mukhopadhyay and de Silva (2009) states that the risk
#' is at most \code{w} for every \eqn{\theta}, every \eqn{m \ge 3}, and
#' every \code{A}; but that procedure uses about twice the sample size
#' \eqn{n^* = A \theta^2 / w} at which a fixed sample has risk \code{w}.
#' The factor 0.565 (Equation 8.3.13) is the improvement of Zacks and
#' Mukhopadhyay (2006), found by exact computation of the risk; the book
#' reports that it keeps the risk "only a shade under" the bound. That is
#' a computed result, not a theorem, and it was recomputed for this
#' package: the risk relative to \code{w} depends only on \eqn{m} and
#' \eqn{n^*}, and exact calculation over first stage sizes from 3 to 100
#' and values of \eqn{n^*} from 1 to 2,000 gives a largest value of 0.990
#' (for a first stage of 14 or 15 observations with \eqn{n^*} near 44).
#' The bound holds throughout that range, with little to spare by design.
#' The guarantee rests on the exponential model.
#'
#' \strong{What it costs.} The expected total sample size exceeds \eqn{n^*}
#' by a factor that depends on the first stage size: about
#' \eqn{0.565 \times 2 (m + 1)^2 / ((m - 1)(m - 2))}, which is 9.04 for
#' \eqn{m = 3}, 1.90 for \eqn{m = 10}, 1.46 for \eqn{m = 20}, and tends to
#' 1.13. See \code{\link{ts_br_exp_mean_pilot}} for the table and for the
#' choice of first stage size, which the source leaves open.
#'
#' The data must be positive, as exponential observations are. Data on a
#' scale so large or so small that the square of a mean cannot be
#' represented in double precision (a mean beyond about \eqn{10^{154}} or
#' below about \eqn{10^{-154}}) are an error that says to express them in
#' other units; no row of the result is ever infinite.
#'
#' @param pilot_data Numeric vector of the first stage observations, at
#'   least three, all positive. The whole first stage must be supplied at
#'   every call, because it alone determines the total sample size.
#' @param additional_data Numeric vector of the second stage observations
#'   collected so far, all positive, or \code{NULL} (the default) when
#'   none have been collected.
#' @param A The loss constant: the loss is \eqn{A\,(\bar{X} - \theta)^2}.
#'   A single positive number, in the units of \code{w} per squared unit
#'   of the data.
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
#'       estimate of \eqn{\theta}.}
#'     \item{\code{risk}}{\eqn{A \bar{X}^2 / n}: the risk
#'       \eqn{A \theta^2 / n} of the mean of a fixed sample of the current
#'       size, with \eqn{\theta} replaced by its estimate. A descriptive
#'       estimate; the guarantee is the one stated in Details.}
#'     \item{\code{optimal_n}}{The total sample size \eqn{N} fixed by the
#'       first stage.}
#'     \item{\code{current_n}}{The number of observations supplied, both
#'       stages together.}
#'   }
#'   Attributes: \code{is_satisfied} (\code{TRUE} once \code{current_n} is at
#'   least \code{optimal_n}), \code{current_n}, \code{criterion_n} (equal to
#'   \code{optimal_n}), \code{loss_A}, \code{risk_bound}, \code{method}, and
#'   \code{notes} (always empty for this function). There is no
#'   \code{conf_level} or \code{width}, and no \code{variance_unusable}
#'   attribute: the total sample size depends on the first stage mean, not on
#'   a variance, so a first stage of equal observations is accepted.
#'
#' @references
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' Mukhopadhyay, N., & Pepe, W. (2006). Exact bounded risk estimation when
#' the terminal sample size and estimator are dependent: The exponential
#' case. \emph{Sequential Analysis}, \emph{25}(1), 85--101.
#' \doi{10.1080/07474940500452254} (As described in Mukhopadhyay and de
#' Silva, 2009, Section 8.3.1, Equations 8.3.9 and 8.3.11.)
#'
#' Zacks, S., & Mukhopadhyay, N. (2006). Bounded risk estimation of the
#' exponential parameter in a two-stage sampling. \emph{Sequential
#' Analysis}, \emph{25}(4), 437--452. \doi{10.1080/07474940600934896} (As
#' described in Mukhopadhyay and de Silva, 2009, Section 8.3.2, Equation
#' 8.3.13.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#'   Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # An operations analyst wants the mean time between service calls
#' # (hours) estimated under squared error loss with A = 100, keeping the
#' # risk at or below 4.
#' set.seed(113)
#' m <- ts_br_exp_mean_pilot(A = 100, w = 4)
#' first_stage <- rexp(m, rate = 1 / 2)
#' check <- ts_br_exp_mean(first_stage, A = 100, w = 4)
#' check
#'
#' # The first stage fixed the total sample size; collect the remainder.
#' n_more <- attr(check, "criterion_n") - attr(check, "current_n")
#' second_stage <- rexp(n_more, rate = 1 / 2)
#' ts_br_exp_mean(first_stage, second_stage, A = 100, w = 4)
#'
#' @seealso \code{\link{ts_br_exp_mean_pilot}} for the first stage size;
#'   \code{\link{seq_mr_exp_mean}} for purely sequential minimum risk
#'   estimation of an exponential mean.
#' @family ts_br_exp_mean
#' @export
ts_br_exp_mean <- function(pilot_data, additional_data = NULL, A, w,
                           na.rm = TRUE) {
  .ts_br_exp_check_Aw(A, w)
  cc <- .ts_combine(pilot_data, additional_data, na.rm, min_n = 3L,
                    positive = TRUE,
                    advice = .ts_advice_first_stage(
                      paste0("The procedure is defined for a first stage ",
                             "of at least three observations."),
                      "ts_br_exp_mean_pilot"))
  xbar <- mean(cc$all)

  # Mukhopadhyay and de Silva (2009): the stopping variable is Equation
  # 8.3.9 (Mukhopadhyay and Pepe, 2006), B_m(A) is Equation 8.3.11, and the
  # factor 0.565 is Equation 8.3.13 (Zacks and Mukhopadhyay, 2006).
  #
  # The squared means must be numbers the rule can use: a mean whose square
  # overflows or underflows is an error that says to change the units, not
  # a message about A and w, a total sample size equal to the first stage,
  # or a risk of Inf or zero.
  B <- 0.565 * (2 * cc$m * (cc$m + 1) * A) / ((cc$m - 1) * (cc$m - 2))
  mean_sq_pilot <- .ts_usable_square(mean(cc$pilot)^2,
                                     "square of the first stage mean",
                                     "the mean is positive")
  optimal_n <- max(cc$m, .as_sample_size(
    B * mean_sq_pilot / w,
    what = "total sample size", advice = .ts_advice_risk()))
  # Plug in risk of the mean of a fixed sample of the current size.
  mean_sq_all <- .ts_usable_square(xbar^2, "square of the sample mean",
                                   "the mean is positive")
  risk <- A * mean_sq_all / cc$n

  term <- c("mean", "risk", "optimal_n", "current_n")
  value <- c(xbar, risk, optimal_n, cc$n)
  .ts_finite_rows(term, value, "restate the loss constant 'A' for those units.")
  new_smsd_check(
    term = term,
    value = value,
    method = "two stage bounded risk estimation (exponential mean)",
    is_satisfied = cc$n >= optimal_n,
    current_n = cc$n,
    criterion_n = optimal_n,
    extra = list(loss_A = A, risk_bound = w)
  )
}
