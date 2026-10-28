# The criterion constant of the minimum risk procedure for a negative
# exponential location, K = (A k Gamma(k + 1) / (c t))^(1 / (t + k)), so
# that the risk minimizing sample size is K * sigma^(k / (t + k)). Computed
# on the log scale so that a large k cannot overflow the gamma function.
#' @noRd
.seq_mr_nexp_log_constant <- function(A, c, k, t) {
  (log(A) + log(k) + lgamma(k + 1) - log(c) - log(t)) / (t + k)
}


#' Pilot Sample Size for Minimum Risk Point Estimation of the Location of a
#' Negative Exponential Distribution
#'
#' Returns the number of observations to collect before the stopping rule of
#' \code{\link{seq_mr_nexp_mean}} is checked for the first time. That
#' procedure is the purely sequential minimum risk point estimation (MRPE)
#' procedure for the location parameter (the lower endpoint) of a negative
#' exponential distribution (Mukhopadhyay and de Silva, 2009).
#'
#' @details
#' If the scale \eqn{\sigma} of the distribution were known, the sample size
#' that minimizes the risk would be
#' \eqn{n^* = K\,\sigma^{k/(t+k)}} with
#' \eqn{K = \{A\,k\,\Gamma(k+1)/(c\,t)\}^{1/(t+k)}} (see
#' \code{\link{seq_mr_nexp_mean}}). The pilot must be chosen before anything
#' is known about \eqn{\sigma}, so it uses only the known constant \eqn{K}:
#' \deqn{m = \max\{m_0,\ \lceil K^{1/(1 + \gamma)} \rceil\}.}
#' With the default \eqn{\gamma = 1} the second term is the square root of
#' \eqn{K}, which grows more slowly than \eqn{n^*} as the cost per
#' observation falls. It therefore does not presume a value of
#' \eqn{\sigma}: a pilot equal to \eqn{n^*} evaluated at \eqn{\sigma = 1}
#' already exceeds the best sample size whenever the scale is below 1 in
#' the units of the data.
#'
#' This is the package's pilot rule. It has the form of the pilot sample
#' size of Chattopadhyay and Kelley (2016, Equation 24). It should not be
#' read as the pilot size prescribed by Mukhopadhyay and de Silva (2009).
#'
#' @param A The loss constant: the price of one unit of estimation error
#'   raised to the power \code{k}, in the same currency as \code{c}, so that
#'   the loss is \eqn{A\,(X_{(1)} - \theta)^k + c\,n^t}. A positive number.
#' @param c The cost constant: sampling \eqn{n} observations costs
#'   \eqn{c\,n^t}. A positive number.
#' @param k The power applied to the estimation error in the loss: 1 for
#'   absolute error, 2 for squared error. Any positive number.
#' @param t The power applied to the sample size in the cost: 1 when each
#'   observation costs the same, more than 1 when later observations cost
#'   more. Any positive number.
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
#' # Shortest possible delivery time (days): an error of one day is valued
#' # at 500 dollars and tracking one more order costs 0.25 dollars.
#' seq_mr_nexp_mean_pilot(A = 500, c = 0.25, k = 1, t = 1)
#'
#' @seealso \code{\link{seq_mr_nexp_mean}}
#' @family seq_mr_nexp_mean
#' @export
seq_mr_nexp_mean_pilot <- function(A, c, k, t, m_0 = 4, gamma = 1) {
  .seq_need(missing(A), "A", "the loss constant")
  .seq_need(missing(c), "c", "the cost constant")
  .seq_need(missing(k), "k", "the power of the estimation error in the loss")
  .seq_need(missing(t), "t", "the power of the sample size in the cost")
  A <- .check_positive(A, "A")
  .mr_check_cost(c)
  k <- .check_positive(k, "k")
  t <- .check_positive(t, "t")
  .seq_pilot(exp(.seq_mr_nexp_log_constant(A, c, k, t)), m_0, gamma,
             advice = "Check 'A', 'c', 'k', and 't'.")
}


#' Minimum Risk Point Estimation of the Location of a Negative Exponential
#' Distribution
#'
#' Given the data collected so far, checks the stopping rule of the purely
#' sequential minimum risk point estimation (MRPE) procedure for the
#' location parameter of a negative exponential distribution (Mukhopadhyay
#' and de Silva, 2009). The negative exponential (two parameter exponential)
#' distribution has density
#' \eqn{\sigma^{-1}\exp\{-(x - \theta)/\sigma\}} for \eqn{x > \theta}:
#' \eqn{\theta} is the location, the smallest value the variable can take (a
#' minimum response time, a guarantee period), and \eqn{\sigma} is the
#' scale. Despite the function's name, the quantity estimated is the
#' location \eqn{\theta}, not the mean (which is \eqn{\theta + \sigma}). The
#' procedure is risk based, so there is no confidence level.
#'
#' The function does not collect data and keeps no state. Call it after each
#' new observation (or batch), starting once
#' \code{\link{seq_mr_nexp_mean_pilot}} observations are in hand, and stop
#' sampling the first time \code{attr(out, "is_satisfied")} is \code{TRUE}.
#'
#' @details
#' The location is estimated by the sample minimum \eqn{X_{(1)}}, and the
#' loss from \eqn{n} observations is
#' \deqn{A\,(X_{(1)} - \theta)^k + c\,n^t .}
#' Here \eqn{A} prices a unit of estimation error raised to the power
#' \eqn{k}, and \eqn{c\,n^t} is the cost of sampling. Because
#' \eqn{X_{(1)} - \theta} is exponential with scale \eqn{\sigma/n}, the
#' risk for a fixed \eqn{n} is
#' \deqn{A\,\Gamma(k+1)\,(\sigma/n)^k + c\,n^t ,}
#' which is smallest at
#' \eqn{n^* = \{A\,k\,\Gamma(k+1)\,\sigma^k/(c\,t)\}^{1/(t+k)}}. The
#' sequential rule replaces \eqn{\sigma} by
#' \eqn{U_n = \sum (X_i - X_{(1)})/(n - 1)}: sampling stops at the first
#' \eqn{n} (at or beyond the pilot size) with
#' \deqn{n \ge \{A\,k\,\Gamma(k+1)\,U_n^k/(c\,t)\}^{1/(t+k)} .}
#' The \code{optimal_n} row is the right hand side rounded up, and the
#' \code{risk} row is the estimated risk at the current sample size,
#' \eqn{A\,\Gamma(k+1)\,(U_n/n)^k + c\,n^t}; both are evaluated at the
#' current data.
#'
#' \strong{The optimality is asymptotic.} As the cost constant falls (so
#' that \eqn{n^*} grows), the expected final sample size approaches
#' \eqn{n^*} and the risk of the procedure approaches the minimum. In
#' simulation (\code{A = 500}, \code{c = 0.01}, each combination of
#' \code{k} and \code{t} equal to 1 or 2, scales 0.25, 1, and 3, the pilot
#' of \code{seq_mr_nexp_mean_pilot()} with its defaults, 40,000
#' replications, Monte Carlo standard error at most 0.007), the ratio of
#' the risk to the minimum was between 1.00 and 1.03 in ten of the twelve
#' conditions, in which \eqn{n^*} ran from 18 to 387. The two exceptions
#' both had \code{k = 2} and \code{t = 2}: the ratio was 1.04 at scale 1,
#' where \eqn{n^*} is 18, and 1.09 at scale 0.25, where \eqn{n^*} is 9
#' (400,000 replications each, standard error 0.002).
#'
#' \strong{Data.} Any finite values are accepted, because the location may
#' be negative. The model is right skewed with a sharp lower endpoint; the
#' function cannot check that.
#'
#' \strong{Constant data.} If every observation is equal, or equal up to
#' rounding error (a range of no more than \eqn{10^{-10}} times the largest
#' absolute value), \eqn{U_n} is zero, or nothing but rounding error, and
#' the rule as printed would be met at once. The function departs from the
#' printed rule deliberately: no scale has been estimated, so the rule is
#' reported as not met, \code{risk} and \code{optimal_n} are \code{NA},
#' \code{attr(out, "variance_unusable")} is \code{TRUE}, and a note says so
#' (see \code{\link{smsd_check_methods}}); data whose real spread is that
#' small beside their size should have a common value subtracted first.
#'
#' \strong{Extreme units.} Data that do vary, but whose spread is at the
#' very edge of double precision, so that the scale estimate cannot be
#' represented (differences among the observations near \eqn{10^{308}} or
#' near \eqn{10^{-323}}), are an error that says to change the units. That
#' case is not reported as a rule not yet met, because collecting more data
#' cannot help.
#'
#' @param x Numeric vector of the observations collected so far (at
#'   least 2).
#' @param A The loss constant: the price of one unit of estimation error
#'   raised to the power \code{k}, in the same currency as \code{c}. A
#'   positive number.
#' @param c The cost constant: sampling \eqn{n} observations costs
#'   \eqn{c\,n^t}. A positive number.
#' @param k The power applied to the estimation error in the loss: 1 for
#'   absolute error, 2 for squared error. Any positive number.
#' @param t The power applied to the sample size in the cost: 1 when each
#'   observation costs the same, more than 1 when later observations cost
#'   more. Any positive number.
#' @param na.rm Logical; should \code{NA} values be removed first? Default
#'   \code{TRUE}. With \code{FALSE}, missing values are an error.
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} with
#'   columns \code{term} and \code{value} and these rows:
#'   \describe{
#'     \item{\code{location}}{The sample minimum, the point estimate of the
#'       location \eqn{\theta}.}
#'     \item{\code{risk}}{The estimated risk at the current sample size,
#'       \eqn{A\,\Gamma(k+1)\,(U_n/n)^k + c\,n^t}.}
#'     \item{\code{optimal_n}}{The sample size the current scale estimate
#'       calls for,
#'       \eqn{\lceil \{A\,k\,\Gamma(k+1)\,U_n^k/(c\,t)\}^{1/(t+k)} \rceil}.}
#'     \item{\code{current_n}}{The number of observations used.}
#'   }
#'   \code{risk} and \code{optimal_n} are \code{NA} when the data are
#'   constant. Attributes: \code{is_satisfied} (\code{TRUE} when
#'   \code{current_n} has reached \code{optimal_n}), \code{current_n},
#'   \code{criterion_n} (equal to \code{optimal_n}), \code{method},
#'   \code{loss_A}, \code{cost_c}, \code{k}, \code{t},
#'   \code{variance_unusable} (\code{TRUE} when the data are constant, so
#'   that no scale has been estimated; the attribute has the same name in
#'   every procedure of the package), and \code{notes}. \code{conf_level}
#'   and \code{width} are \code{NULL}.
#'
#' @references
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @author Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in},
#'   Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # Shortest possible delivery time for an online retailer, in days. An
#' # error of one day is valued at 500 dollars and tracking one more order
#' # costs 0.25 dollars. Times are simulated here as a minimum of 2 days
#' # plus an exponential delay averaging 1.5 days, so the best sample size
#' # is sqrt(500 * 1.5 / 0.25), about 55.
#' set.seed(113)
#' m <- seq_mr_nexp_mean_pilot(A = 500, c = 0.25, k = 1, t = 1)
#' days <- 2 + rexp(m, rate = 1 / 1.5)
#' out <- seq_mr_nexp_mean(days, A = 500, c = 0.25, k = 1, t = 1)
#' out
#'
#' # Add one order at a time until the stopping rule is met.
#' while (!attr(out, "is_satisfied")) {
#'   days <- c(days, 2 + rexp(1, rate = 1 / 1.5))
#'   out <- seq_mr_nexp_mean(days, A = 500, c = 0.25, k = 1, t = 1)
#' }
#' out
#'
#' @seealso \code{\link{seq_mr_nexp_mean_pilot}};
#'   \code{\link{seq_fwd_nexp_mean}} for a fixed width confidence interval
#'   for the same location.
#' @family seq_mr_nexp_mean
#' @export
seq_mr_nexp_mean <- function(x, A, c, k, t, na.rm = TRUE) {
  .seq_need(missing(A), "A", "the loss constant")
  .seq_need(missing(c), "c", "the cost constant")
  .seq_need(missing(k), "k", "the power of the estimation error in the loss")
  .seq_need(missing(t), "t", "the power of the sample size in the cost")
  A <- .check_positive(A, "A")
  .mr_check_cost(c)
  k <- .check_positive(k, "k")
  t <- .check_positive(t, "t")
  x <- .mr_check_data(x, na.rm, arg = "x")
  n <- length(x)
  loc <- min(x)
  U <- sum(x - loc) / (n - 1)
  unusable <- .seq_zero_spread(U, x, "scale estimate")

  if (unusable) {
    optimal_n <- NA_real_
    risk <- NA_real_
  } else {
    # n* = K * U^(k / (t + k)) and risk = A Gamma(k + 1) (U / n)^k + c n^t,
    # with the gamma function handled on the log scale.
    optimal_n <- .seq_criterion(
      exp(.seq_mr_nexp_log_constant(A, c, k, t) + k * log(U) / (t + k)),
      advice = "Check 'A', 'c', 'k', and 't', and the units of 'x'.")
    risk <- .seq_finite_risk(
      exp(log(A) + lgamma(k + 1) + k * (log(U) - log(n))) + c * n^t,
      A = A, c = c, k = k, t = t)
  }

  .seq_check(
    term = c("location", "risk", "optimal_n", "current_n"),
    value = c(loc, risk, optimal_n, n),
    method = "minimum risk point estimation (negative exponential location)",
    is_satisfied = !unusable && n >= optimal_n,
    current_n = n,
    criterion_n = optimal_n,
    extra = list(loss_A = A, cost_c = c, k = k, t = t,
                 variance_unusable = unusable),
    notes = if (unusable) .seq_notes_unusable(x)
  )
}
