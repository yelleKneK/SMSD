#' Pilot Sample Size for Minimum Risk Point Estimation of an Exponential
#' Mean
#'
#' Returns the number of observations to collect before the stopping rule of
#' \code{\link{seq_mr_exp_mean}} is checked for the first time. That
#' procedure is the purely sequential minimum risk point estimation (MRPE)
#' procedure for the mean of an exponential population (Mukhopadhyay and de
#' Silva, 2009).
#'
#' @details
#' If the population mean \eqn{\lambda} were known, the sample size that
#' minimizes the risk would be \eqn{n^* = (A/c)^{1/2}\lambda} (the standard
#' deviation of an exponential distribution equals its mean). The pilot must
#' be chosen before anything is known about \eqn{\lambda}, so it uses only
#' the known constant \eqn{(A/c)^{1/2}}:
#' \deqn{m = \max\{m_0,\ \lceil (A/c)^{1/(2 + 2\gamma)} \rceil\}.}
#' With the default \eqn{\gamma = 1} the second term is \eqn{(A/c)^{1/4}},
#' which grows more slowly than \eqn{n^*} as the cost per observation
#' falls. It therefore does not presume a value of \eqn{\lambda}: a pilot
#' equal to \eqn{n^*} evaluated at \eqn{\lambda = 1} already exceeds the
#' best sample size whenever the mean is below 1 in the units of the data.
#'
#' This is the package's pilot rule. It has the form of the pilot sample
#' size of Chattopadhyay and Kelley (2016, Equation 24). It should not be
#' read as the pilot size prescribed by Mukhopadhyay and de Silva (2009).
#'
#' @param A The loss constant: the price of one unit of squared estimation
#'   error, in the same currency as \code{c}, so that the loss is
#'   \eqn{A(\bar{X}_n - \lambda)^2 + c\,n}. A positive number.
#' @param c The cost of sampling one observation. A positive number.
#' @param m_0 The smallest pilot sample size allowed. A whole number of at
#'   least 2; default 4.
#' @param gamma Positive constant that sets how fast the pilot grows as the
#'   cost per observation falls; see Details. Default 1.
#'
#' @return A length one integer giving the pilot sample size.
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
#' @author Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in},
#'   Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # An error of one minute (squared) in the mean time to resolve a support
#' # ticket is valued at 100 dollars; timing one more ticket costs 0.04.
#' seq_mr_exp_mean_pilot(A = 100, c = 0.04)
#'
#' @seealso \code{\link{seq_mr_exp_mean}}
#' @family seq_mr_exp_mean
#' @export
seq_mr_exp_mean_pilot <- function(A, c, m_0 = 4, gamma = 1) {
  .seq_need(missing(A), "A", "the loss constant")
  .seq_need(missing(c), "c", "the cost of sampling one observation")
  A <- .check_positive(A, "A")
  .mr_check_cost(c)
  .seq_pilot(sqrt(A / c), m_0, gamma,
             advice = "Check 'A' and the cost 'c'.")
}


#' Minimum Risk Point Estimation of an Exponential Mean
#'
#' Given the data collected so far, checks the stopping rule of the purely
#' sequential minimum risk point estimation (MRPE) procedure for the mean of
#' an exponential population, such as a mean waiting time or time to failure
#' (Mukhopadhyay and de Silva, 2009). Estimation error and sampling both
#' cost something; the procedure looks for the sample size at which their
#' expected total, the risk, is smallest. That sample size depends on the
#' unknown mean, so observations are collected until the sample size reaches
#' the size the current estimate calls for. The procedure is risk based, so
#' there is no confidence level.
#'
#' The function does not collect data and keeps no state. Call it after each
#' new observation (or batch), starting once
#' \code{\link{seq_mr_exp_mean_pilot}} observations are in hand, and stop
#' sampling the first time \code{attr(out, "is_satisfied")} is \code{TRUE}.
#'
#' @details
#' The loss from estimating the mean \eqn{\lambda} by the sample mean of
#' \eqn{n} observations is \eqn{A(\bar{X}_n - \lambda)^2 + c\,n}: \eqn{A}
#' prices a unit of squared estimation error and \eqn{c} is the cost of one
#' observation. For exponential data the variance of the sample mean is
#' \eqn{\lambda^2/n}, so the risk for a fixed \eqn{n} is
#' \eqn{A\lambda^2/n + c\,n}, which is smallest at
#' \eqn{n^* = (A/c)^{1/2}\lambda}, where it equals \eqn{2 c\,n^*}. The
#' sequential rule replaces \eqn{\lambda} by the sample mean: sampling
#' stops at the first \eqn{n} (at or beyond the pilot size) with
#' \deqn{n \ge (A/c)^{1/2}\,\bar{X}_n .}
#' The \code{optimal_n} row is the right hand side rounded up, and the
#' \code{risk} row is the estimated risk at the current sample size,
#' \eqn{A\bar{X}_n^2/n + c\,n}; both are evaluated at the current data.
#'
#' \strong{The optimality is asymptotic.} As the cost per observation falls
#' (so that \eqn{n^*} grows), the expected final sample size approaches
#' \eqn{n^*} and the risk of the procedure approaches the minimum
#' \eqn{2 c\,n^*}. When \eqn{n^*} is small the risk is noticeably above the
#' minimum, more so than for a normal mean, because the sample mean of
#' skewed data is itself the estimate of spread and a low early mean stops
#' the sampling. In simulation (exponential data, the pilot of
#' \code{seq_mr_exp_mean_pilot()} with its defaults, 40,000 replications,
#' Monte Carlo standard error at most 0.006), with a mean of 1 the ratio of
#' the risk to the minimum was 1.08, 1.16, 1.15, 1.05, and 1.02 when
#' \eqn{n^*} was 10, 15, 20, 50, and 100 (that is, when \eqn{A/c} was 100,
#' 225, 400, 2,500, and 10,000, for which the pilot is 4, 4, 5, 8, and 10).
#' The ratio is not a function of \eqn{n^*} alone: the pilot depends only
#' on \eqn{A/c}, so the same \eqn{n^*} reached with another mean starts
#' from a different pilot. With a mean of 0.25 (pilots of 7 to 20) the five
#' ratios were 0.95, 1.04, 1.07, 1.04, and 1.01, and with a mean of 4
#' (pilots of 4 or 5) they were 1.08, 1.16, 1.16, 1.06, and 1.02. (A ratio
#' below 1 is possible because the final sample size depends on the data.)
#'
#' \strong{Data.} The exponential model is for positive data. Zero or
#' negative values are an error: they are not possible under the model, and
#' a negative sample mean would give a negative sample size.
#'
#' @param x Numeric vector of the observations collected so far (at least
#'   1). All values must be positive.
#' @param A The loss constant: the price of one unit of squared estimation
#'   error, in the same currency as \code{c}. A positive number.
#' @param c The cost of sampling one observation. A positive number.
#' @param na.rm Logical; should \code{NA} values be removed first? Default
#'   \code{TRUE}. With \code{FALSE}, missing values are an error.
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} with
#'   columns \code{term} and \code{value} and these rows:
#'   \describe{
#'     \item{\code{mean}}{The sample mean, the point estimate.}
#'     \item{\code{risk}}{The estimated risk at the current sample size,
#'       \eqn{A\bar{X}_n^2/n + c\,n}.}
#'     \item{\code{optimal_n}}{The sample size the current estimate of the
#'       mean calls for, \eqn{\lceil (A/c)^{1/2}\,\bar{X}_n \rceil}.}
#'     \item{\code{current_n}}{The number of observations used.}
#'   }
#'   Attributes: \code{is_satisfied} (\code{TRUE} once \code{current_n} is at
#'   least \code{optimal_n}), \code{current_n}, \code{criterion_n} (equal to
#'   \code{optimal_n}), \code{method}, \code{loss_A}, \code{cost_c}, and
#'   \code{notes}. \code{conf_level} and \code{width} are \code{NULL}.
#'
#' @references
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @author Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in},
#'   Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # Mean time to resolve a support ticket, in minutes. An error of one
#' # minute (squared) is valued at 100 dollars and timing one more ticket
#' # costs 0.04 dollars. Times are simulated here with a mean of 6
#' # minutes, so the best sample size is sqrt(100 / 0.04) * 6 = 300.
#' set.seed(113)
#' m <- seq_mr_exp_mean_pilot(A = 100, c = 0.04)
#' minutes <- rexp(m, rate = 1 / 6)
#' out <- seq_mr_exp_mean(minutes, A = 100, c = 0.04)
#' out
#'
#' # Add one ticket at a time until the stopping rule is met.
#' while (!attr(out, "is_satisfied")) {
#'   minutes <- c(minutes, rexp(1, rate = 1 / 6))
#'   out <- seq_mr_exp_mean(minutes, A = 100, c = 0.04)
#' }
#' out
#'
#' @seealso \code{\link{seq_mr_exp_mean_pilot}};
#'   \code{\link{ts_br_exp_mean}} for two stage bounded risk estimation of
#'   the same mean; \code{\link{seq_mr_dfr_mean_arb}}, which makes no
#'   distributional assumption.
#' @family seq_mr_exp_mean
#' @export
seq_mr_exp_mean <- function(x, A, c, na.rm = TRUE) {
  .seq_need(missing(A), "A", "the loss constant")
  .seq_need(missing(c), "c", "the cost of sampling one observation")
  A <- .check_positive(A, "A")
  .mr_check_cost(c)
  x <- .mr_check_data(x, na.rm, min_n = 1L, arg = "x", positive = TRUE)
  n <- length(x)
  xbar <- mean(x)

  optimal_n <- .seq_criterion(
    sqrt(A / c) * xbar,
    advice = "Check 'A', the cost 'c', and the units of 'x'.")
  risk <- .seq_finite_risk(A * xbar^2 / n + c * n, A = A, c = c)

  .seq_check(
    term = c("mean", "risk", "optimal_n", "current_n"),
    value = c(xbar, risk, optimal_n, n),
    method = "minimum risk point estimation (exponential mean)",
    is_satisfied = n >= optimal_n,
    current_n = n,
    criterion_n = optimal_n,
    extra = list(loss_A = A, cost_c = c)
  )
}
