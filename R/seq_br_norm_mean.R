#' Pilot Sample Size for Bounded Risk Point Estimation of a Normal Mean
#'
#' Returns the number of observations to collect before the stopping rule of
#' \code{\link{seq_br_norm_mean}} is checked for the first time. That
#' procedure is the purely sequential bounded risk point estimator of the
#' mean of a normal population with unknown variance (Mukhopadhyay and de
#' Silva, 2009).
#'
#' @details
#' If the population standard deviation \eqn{\sigma} were known, the
#' smallest sample size with risk at most \code{w} would be
#' \eqn{n^* = (B/w)^{2/k}\sigma^2}, where
#' \eqn{B = A\,2^{k/2}\,\Gamma\{(k+1)/2\}/\sqrt{\pi}} (see
#' \code{\link{seq_br_norm_mean}}). The pilot must be chosen before anything
#' is known about \eqn{\sigma}, so it uses only the known constant
#' \eqn{(B/w)^{2/k}}:
#' \deqn{m = \max\{m_0,\ \lceil \{(B/w)^{2/k}\}^{1/(1 + \gamma)} \rceil\}.}
#' With the default \eqn{\gamma = 1} the second term is the square root of
#' that constant, which grows more slowly than \eqn{n^*} as \code{w}
#' shrinks. It therefore does not presume a value of \eqn{\sigma}: a pilot
#' equal to \eqn{n^*} evaluated at \eqn{\sigma = 1} is far too large when
#' the standard deviation is well below 1 in the units of the data, and far
#' too small to protect against early stopping when it is well above 1.
#'
#' This is the package's pilot rule, of the same form as the pilot sample
#' sizes of the package's other sequential procedures (Chattopadhyay and
#' Kelley, 2016, Equation 24; Kelley, Darku, and Chattopadhyay, 2018,
#' Equation 21). It should not be read as the pilot size prescribed by
#' Mukhopadhyay and de Silva (2009).
#'
#' @param A The loss constant: the loss from estimating the mean
#'   \eqn{\mu} by the sample mean is \eqn{A\,|\bar{X}_n - \mu|^k}. A
#'   positive number; it sets the units in which the risk bound \code{w} is
#'   expressed (use \code{A = 1} to bound the mean absolute or mean squared
#'   error itself).
#' @param k The power in the loss: 1 for absolute error, 2 for squared
#'   error.
#' @param w The risk bound: the largest expected loss that is acceptable. A
#'   positive number, in the units of \code{A} times the data's units to the
#'   power \code{k}.
#' @param m_0 The smallest pilot sample size allowed. A whole number of at
#'   least 2; default 10. A larger value reduces the amount by which the
#'   risk can exceed the bound when the bound is loose (see the Details of
#'   \code{\link{seq_br_norm_mean}}); it costs at most \code{m_0}
#'   observations.
#' @param gamma Positive constant that sets how fast the pilot grows as the
#'   bound tightens; see Details. Default 1.
#'
#' @return A length one integer giving the pilot sample size.
#'
#' @references
#' Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient
#' of variation with minimum risk: A sequential method for minimizing
#' sampling error and study cost. \emph{Multivariate Behavioral Research},
#' \emph{51}(5), 627--648. \doi{10.1080/00273171.2016.1203279}
#'
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#' parameter estimation for a general class of effect sizes: A sequential
#' approach. \emph{Psychological Methods}, \emph{23}(2), 226--243.
#' \doi{10.1037/met0000127}
#'
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#'   Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # Mean absolute error of at most 0.25 points on a satisfaction scale
#' seq_br_norm_mean_pilot(A = 1, k = 1, w = 0.25)
#' # Mean squared error of at most 0.01
#' seq_br_norm_mean_pilot(A = 1, k = 2, w = 0.01)
#'
#' @seealso \code{\link{seq_br_norm_mean}}
#' @family seq_br_norm_mean
#' @export
seq_br_norm_mean_pilot <- function(A, k, w, m_0 = 10, gamma = 1) {
  .br_check_Akw(A, k, w)
  .seq_pilot(.seq_br_scale(A, k, w), m_0, gamma,
             advice = "Check 'A' and the risk bound 'w'.")
}


#' Bounded Risk Point Estimation of a Normal Mean
#'
#' Given the data collected so far, checks the stopping rule of the purely
#' sequential bounded risk point estimator of the mean of a normal
#' population with unknown variance (Mukhopadhyay and de Silva, 2009). The
#' aim is a sample mean whose expected loss (its risk) is no larger than a
#' bound \code{w} chosen in advance. No fixed sample size can promise that,
#' because the risk depends on the unknown variance, so observations are
#' collected until the sample size reaches the size the current variance
#' estimate calls for. The procedure is risk based, so there is no
#' confidence level.
#'
#' The function does not collect data and keeps no state. Call it after each
#' new observation (or batch), starting once
#' \code{\link{seq_br_norm_mean_pilot}} observations are in hand, and stop
#' sampling the first time \code{attr(out, "is_satisfied")} is \code{TRUE}.
#'
#' @details
#' The loss from estimating the mean \eqn{\mu} by the sample mean is
#' \eqn{A\,|\bar{X}_n - \mu|^k} with \eqn{k} equal to 1 (absolute error) or
#' 2 (squared error). For normal data and a fixed \eqn{n} its expectation is
#' \deqn{B\,(\sigma^2/n)^{k/2}, \qquad
#'   B = A\,2^{k/2}\,\Gamma\{(k+1)/2\}/\sqrt{\pi},}
#' that is, \eqn{B = A\sqrt{2/\pi}} for \eqn{k = 1} and \eqn{B = A} for
#' \eqn{k = 2}. The risk is at most \code{w} when
#' \eqn{n \ge n^* = (B/w)^{2/k}\sigma^2}. The sequential rule replaces
#' \eqn{\sigma^2} by the sample variance \eqn{S_n^2}: sampling stops at the
#' first \eqn{n} (at or beyond the pilot size) with
#' \deqn{n \ge (B/w)^{2/k} S_n^2 .}
#' The \code{optimal_n} row is the right hand side rounded up, and the
#' \code{risk} row is the estimated risk at the current sample size,
#' \eqn{B\,(S_n^2/n)^{k/2}}; both are evaluated at the current data.
#'
#' \strong{The risk bound is asymptotic.} The risk of the sequential
#' procedure approaches \code{w} as \code{w} shrinks (so that the required
#' sample size grows), and it approaches it from above: for squared error
#' and normal data the risk is \eqn{A\sigma^2 E(1/N)}, where \eqn{N} is the
#' final sample size, which exceeds \code{w} whenever the expected final
#' sample size is at or below \eqn{n^*}. With a loose bound the excess is
#' not small. In simulation (normal data with standard deviation 3,
#' \code{A = 2}, the pilot of \code{seq_br_norm_mean_pilot()} with its
#' defaults, which was 10 at every bound used, 40,000 replications), the
#' ratio of the risk to \code{w} for squared error (\code{k = 2}) was 1.11
#' when the bound called for 15 observations, highest, about 1.22, when it
#' called for 25 to 30, and then 1.20, 1.05, and 1.01 when it called for
#' 36, 100, and 400 (\code{w} of 0.5, 0.18, and 0.045). These were computed
#' from the identity above, with a Monte Carlo standard error of at most
#' 0.003. For absolute error (\code{k = 1}) the ratio was 1.04 at 15
#' observations, highest, between 1.08 and 1.09, at 25 to 30, and 1.07,
#' 1.02, and 1.00 at 36, 100, and 400 (standard error 0.004).
#' A larger \code{m_0} in the pilot reduces the excess for loose bounds.
#' \code{\link{ts_br_norm_mean}} is the two stage alternative, whose risk
#' stays at or under the bound for normal data at the price of a larger
#' expected sample size.
#'
#' \strong{Constant data.} If every observation is equal, or equal up to
#' rounding error (a range of no more than \eqn{10^{-10}} times the largest
#' absolute value, as with gain scores of 0.3 computed as posttest minus
#' pretest), the sample variance is zero, or nothing but rounding error, and
#' the rule as printed would be met at once. The function departs from the
#' printed rule deliberately: no variance has been estimated, so the rule is
#' reported as not met, \code{risk} and \code{optimal_n} are \code{NA},
#' \code{attr(out, "variance_unusable")} is \code{TRUE}, and a note says so
#' (see \code{\link{smsd_check_methods}}); data whose real spread is that
#' small beside their size should have a common value subtracted first.
#'
#' \strong{Extreme units.} Data that do vary, but whose spread is so large
#' or so small that the sample variance cannot be represented in double
#' precision (differences among the observations beyond about
#' \eqn{10^{154}} or below about \eqn{10^{-162}}), are an error that says
#' to change the units. That case is not reported as a rule not yet met,
#' because collecting more data cannot help.
#'
#' @param x Numeric vector of the observations collected so far (at
#'   least 2).
#' @param A The loss constant: the loss from estimating the mean
#'   \eqn{\mu} by the sample mean is \eqn{A\,|\bar{X}_n - \mu|^k}. A
#'   positive number; it sets the units in which the risk bound \code{w} is
#'   expressed (use \code{A = 1} to bound the mean absolute or mean squared
#'   error itself).
#' @param k The power in the loss: 1 for absolute error, 2 for squared
#'   error.
#' @param w The risk bound: the largest expected loss that is acceptable. A
#'   positive number, in the units of \code{A} times the data's units to the
#'   power \code{k}.
#' @param na.rm Logical; should \code{NA} values be removed first? Default
#'   \code{TRUE}. With \code{FALSE}, missing values are an error.
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} with
#'   columns \code{term} and \code{value} and these rows:
#'   \describe{
#'     \item{\code{mean}}{The sample mean, the point estimate.}
#'     \item{\code{sd}}{The sample standard deviation.}
#'     \item{\code{risk}}{The estimated risk at the current sample size,
#'       \eqn{B\,(S_n^2/n)^{k/2}}.}
#'     \item{\code{optimal_n}}{The sample size the current variance estimate
#'       calls for, \eqn{\lceil (B/w)^{2/k} S_n^2 \rceil}.}
#'     \item{\code{current_n}}{The number of observations used.}
#'   }
#'   \code{risk} and \code{optimal_n} are \code{NA} when the data are
#'   constant. Attributes: \code{is_satisfied} (\code{TRUE} when
#'   \code{current_n} has reached \code{optimal_n}, which is when the
#'   estimated risk is at most \code{w}; the rounding up of
#'   \code{optimal_n} allows for rounding error of one part in
#'   \eqn{10^{10}}, so when the estimated risk and \code{w} agree to about
#'   ten significant digits the two statements can differ, and
#'   \code{is_satisfied} is the authority), \code{current_n},
#'   \code{criterion_n} (equal to \code{optimal_n}), \code{method},
#'   \code{loss_A}, \code{k}, \code{risk_bound}, \code{variance_unusable}
#'   (\code{TRUE} when the data are constant, so that no variance has been
#'   estimated), and \code{notes}. \code{conf_level} and \code{width} are
#'   \code{NULL}.
#'
#' @references
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#'   Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # Mean rating on a 7 point satisfaction scale, with a mean absolute
#' # error of at most 0.25 points. Ratings are simulated here as normal
#' # with a standard deviation of 1.2.
#' set.seed(113)
#' m <- seq_br_norm_mean_pilot(A = 1, k = 1, w = 0.25)
#' ratings <- rnorm(m, mean = 4.8, sd = 1.2)
#' out <- seq_br_norm_mean(ratings, A = 1, k = 1, w = 0.25)
#' out
#'
#' # Add one respondent at a time until the stopping rule is met.
#' while (!attr(out, "is_satisfied")) {
#'   ratings <- c(ratings, rnorm(1, mean = 4.8, sd = 1.2))
#'   out <- seq_br_norm_mean(ratings, A = 1, k = 1, w = 0.25)
#' }
#' out
#'
#' @seealso \code{\link{seq_br_norm_mean_pilot}};
#'   \code{\link{ts_br_norm_mean}} for the two stage procedure;
#'   \code{\link{seq_mr_norm_mean}}, which trades the same loss against the
#'   cost of sampling instead of bounding it.
#' @family seq_br_norm_mean
#' @export
seq_br_norm_mean <- function(x, A, k, w, na.rm = TRUE) {
  .br_check_Akw(A, k, w)
  x <- .mr_check_data(x, na.rm, arg = "x")
  n <- length(x)
  s2 <- stats::var(x)
  unusable <- .seq_zero_spread(s2, x)

  if (unusable) {
    optimal_n <- NA_real_
    risk <- NA_real_
  } else {
    optimal_n <- .seq_criterion(
      .seq_br_scale(A, k, w) * s2,
      advice = "Check 'A', the risk bound 'w', and the units of 'x'.")
    # B (S^2 / n)^(k / 2), with the power written out for k = 1 and k = 2.
    risk <- .seq_finite_risk(
      .br_constant(A, k) * if (k == 2) s2 / n else sqrt(s2 / n),
      A = A, k = k)
  }

  .seq_check(
    term = c("mean", "sd", "risk", "optimal_n", "current_n"),
    value = c(mean(x), .seq_sd(x, s2), risk, optimal_n, n),
    method = "bounded risk estimation (normal mean)",
    is_satisfied = !unusable && n >= optimal_n,
    current_n = n,
    criterion_n = optimal_n,
    extra = list(loss_A = A, k = k, risk_bound = w,
                 variance_unusable = unusable),
    notes = if (unusable) .seq_notes_unusable(x)
  )
}
