#' Pilot Sample Size for a Fixed Width Confidence Interval for the Location
#' of a Negative Exponential Distribution
#'
#' Returns the number of observations to collect before the stopping rule of
#' \code{\link{seq_fwd_nexp_mean}} is checked for the first time. That
#' procedure is the purely sequential fixed width confidence interval for
#' the location parameter (the lower endpoint, or minimum guaranteed value)
#' of a negative exponential distribution (Mukhopadhyay and de Silva, 2009).
#'
#' @details
#' If the scale \eqn{\sigma} of the distribution were known, the sample size
#' that gives the interval \eqn{[X_{(1)} - d,\ X_{(1)}]} the stated
#' confidence level would be \eqn{n^* = \log(1/\alpha)\,\sigma / d}. The
#' pilot must be chosen before anything is known about \eqn{\sigma}, so it
#' uses only the known constant \eqn{\log(1/\alpha)/d}:
#' \deqn{m = \max\{m_0,\ \lceil (\log(1/\alpha)/d)^{1/(1 + \gamma)} \rceil\}.}
#' With the default \eqn{\gamma = 1} this is the square root of that
#' constant, which grows more slowly than \eqn{n^*} as \code{d} shrinks and
#' so does not presume a value of \eqn{\sigma}.
#'
#' This is the package's pilot rule, of the same form as the pilot sample
#' sizes of the package's other sequential procedures (Chattopadhyay and
#' Kelley, 2016, Equation 24). It should not be read as the pilot size
#' prescribed by Mukhopadhyay and de Silva (2009).
#'
#' @param d The full width of the confidence interval, in the units of the
#'   data: the interval reported runs from the sample minimum minus \code{d}
#'   up to the sample minimum. A positive number. (In the fixed width
#'   procedures for a mean, \code{d} is a half width; here the interval is
#'   one-sided, so \code{d} is its whole width.)
#' @param conf_level The confidence level (default 0.95).
#' @param m_0 The smallest pilot sample size allowed. A whole number of at
#'   least 2; default 10. A larger value protects against stopping too early
#'   when the target is wide (see the Details of
#'   \code{\link{seq_fwd_nexp_mean}}); it costs at most \code{m_0}
#'   observations.
#' @param gamma Positive constant that sets how fast the pilot grows as the
#'   target narrows; see Details. Default 1.
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
#' @author Ken Kelley \email{kkelley@@nd.edu},
#'   Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # Shortest possible response time (seconds), to within 0.05 seconds
#' seq_fwd_nexp_mean_pilot(d = 0.05, conf_level = 0.95)
#' seq_fwd_nexp_mean_pilot(d = 0.005, conf_level = 0.95)
#'
#' @seealso \code{\link{seq_fwd_nexp_mean}}
#' @family seq_fwd_nexp_mean
#' @export
seq_fwd_nexp_mean_pilot <- function(d, conf_level = 0.95, m_0 = 10,
                                    gamma = 1) {
  .ts_check_width(d, "the confidence interval full width")
  alpha <- .ts_alpha(conf_level)
  # The criterion is log(1 / alpha) * U / d, linear in 1 / d, so the pilot's
  # constant is log(1 / alpha) / d. (The fixed width procedures for a mean
  # have criteria quadratic in 1 / d and use (z / d)^2.)
  .seq_pilot(-log(alpha) / d, m_0, gamma,
             advice = "Check 'd' (the full width of the confidence interval).")
}


#' Fixed Width Confidence Interval for the Location of a Negative
#' Exponential Distribution
#'
#' Given the data collected so far, checks the stopping rule of the purely
#' sequential fixed width confidence interval for the location parameter of
#' a negative exponential distribution (Mukhopadhyay and de Silva, 2009).
#' The negative exponential (two parameter exponential) distribution has
#' density \eqn{\sigma^{-1}\exp\{-(x - \theta)/\sigma\}} for
#' \eqn{x > \theta}: \eqn{\theta} is the location, the smallest value the
#' variable can take (a minimum response time, a guarantee period), and
#' \eqn{\sigma} is the scale. Despite the function's name, the quantity
#' estimated is the location \eqn{\theta}, not the mean (which is
#' \eqn{\theta + \sigma}).
#'
#' The function does not collect data and keeps no state. Call it after each
#' new observation (or batch), starting once
#' \code{\link{seq_fwd_nexp_mean_pilot}} observations are in hand, and stop
#' sampling the first time \code{attr(out, "is_satisfied")} is \code{TRUE}.
#'
#' @details
#' The location is estimated by the sample minimum \eqn{X_{(1)}}, which can
#' only overestimate it, so the interval is
#' \eqn{[X_{(1)} - d,\ X_{(1)}]}, of width \code{d}. Because
#' \eqn{X_{(1)} - \theta} is exponential with scale \eqn{\sigma/n}, the
#' interval covers \eqn{\theta} with probability
#' \eqn{1 - \exp(-n d/\sigma)}, which is at least \eqn{1 - \alpha} when
#' \eqn{n \ge n^* = \log(1/\alpha)\,\sigma/d}. The sequential rule replaces
#' \eqn{\sigma} by \eqn{U_n = \sum (X_i - X_{(1)})/(n - 1)}: sampling stops
#' at the first \eqn{n} (at or beyond the pilot size) with
#' \deqn{n \ge \log(1/\alpha)\, U_n / d .}
#' The \code{optimal_n} row is the right hand side rounded up, evaluated at
#' the current data.
#'
#' \strong{The confidence level is asymptotic.} It is attained in the limit
#' as \code{d} shrinks, not at every \code{d}; for wide targets sampling
#' tends to stop early on an underestimate of the scale. In simulation
#' (scale 2, nominal 0.95, the pilot of \code{seq_fwd_nexp_mean_pilot()}
#' with its defaults, 40,000 replications, Monte Carlo standard error
#' 0.001), coverage was 0.938 at \code{d = 0.4} (where the required sample
#' size \eqn{n^*} is about 15). It was lowest, about 0.930, when \eqn{n^*}
#' was 25, and rose from there: 0.941 at \eqn{n^*} of 40, 0.944 at
#' \code{d = 0.1} (\eqn{n^*} about 60), and 0.948 at \code{d = 0.025}
#' (\eqn{n^*} about 240; 400,000 replications).
#' \code{\link{ts_fwd_nexp_mean}} is the two stage alternative, whose
#' confidence level is guaranteed for negative exponential data.
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
#' reported as not met, \code{optimal_n} is \code{NA},
#' \code{attr(out, "variance_unusable")} is \code{TRUE}, and a note says so
#' (see \code{\link{smsd_check_methods}}); data whose real spread is that
#' small beside their size should have a common value subtracted first.
#' The limits are still reported (the sample minimum minus \code{d}, and
#' the sample minimum, which involve no estimate of the scale), but they are
#' then not an interval to report: \code{summary()}, \code{glance()}, and
#' \code{tidy()} return them without the note, so test
#' \code{attr(out, "variance_unusable")} before using them.
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
#' @param d The full width of the confidence interval, in the units of the
#'   data: the interval reported runs from the sample minimum minus \code{d}
#'   up to the sample minimum. A positive number. (In the fixed width
#'   procedures for a mean, \code{d} is a half width; here the interval is
#'   one-sided, so \code{d} is its whole width.)
#' @param conf_level The confidence level (default 0.95).
#' @param na.rm Logical; should \code{NA} values be removed first? Default
#'   \code{TRUE}. With \code{FALSE}, missing values are an error.
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} with
#'   columns \code{term} and \code{value} and these rows:
#'   \describe{
#'     \item{\code{lower_limit}}{The sample minimum minus \code{d}.}
#'     \item{\code{location}}{The sample minimum, the point estimate of the
#'       location \eqn{\theta}.}
#'     \item{\code{upper_limit}}{The sample minimum.}
#'     \item{\code{optimal_n}}{The sample size the current scale estimate
#'       calls for, \eqn{\lceil \log(1/\alpha)\, U_n / d \rceil} (\code{NA}
#'       when the data are constant).}
#'     \item{\code{current_n}}{The number of observations used.}
#'   }
#'   Attributes: \code{is_satisfied} (\code{TRUE} once \code{current_n} is at
#'   least \code{optimal_n}), \code{current_n}, \code{criterion_n} (equal to
#'   \code{optimal_n}), \code{conf_level}, \code{width} (the full width
#'   \code{d}), \code{method}, \code{variance_unusable} (\code{TRUE} when the
#'   data are constant, so that no scale has been estimated; the attribute has
#'   the same name in every procedure of the package), and \code{notes}.
#'
#' @references
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#'   Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # Shortest possible response time on a reaction task (seconds), to
#' # within 0.05 seconds with 95% confidence. Times are simulated as a
#' # minimum of 0.20 seconds plus an exponential delay averaging 0.30.
#' set.seed(113)
#' m <- seq_fwd_nexp_mean_pilot(d = 0.05, conf_level = 0.95)
#' times <- 0.20 + rexp(m, rate = 1 / 0.30)
#' out <- seq_fwd_nexp_mean(times, d = 0.05, conf_level = 0.95)
#' out
#'
#' # Add one trial at a time until the stopping rule is met.
#' while (!attr(out, "is_satisfied")) {
#'   times <- c(times, 0.20 + rexp(1, rate = 1 / 0.30))
#'   out <- seq_fwd_nexp_mean(times, d = 0.05, conf_level = 0.95)
#' }
#' out
#'
#' @seealso \code{\link{seq_fwd_nexp_mean_pilot}};
#'   \code{\link{ts_fwd_nexp_mean}} for the two stage procedure;
#'   \code{\link{seq_mr_nexp_mean}} for minimum risk point estimation of the
#'   same location.
#' @family seq_fwd_nexp_mean
#' @export
seq_fwd_nexp_mean <- function(x, d, conf_level = 0.95, na.rm = TRUE) {
  .ts_check_width(d, "the confidence interval full width")
  alpha <- .ts_alpha(conf_level)
  x <- .mr_check_data(x, na.rm, arg = "x")
  n <- length(x)
  loc <- min(x)
  U <- sum(x - loc) / (n - 1)
  unusable <- .seq_zero_spread(U, x, "scale estimate")

  optimal_n <- if (unusable) {
    NA_real_
  } else {
    .seq_criterion(
      -log(alpha) * U / d,
      advice = paste("Check 'd' (the full width of the confidence",
                     "interval) and that it is in the units of 'x'."))
  }

  .seq_check(
    term = c("lower_limit", "location", "upper_limit", "optimal_n",
             "current_n"),
    value = c(loc - d, loc, loc, optimal_n, n),
    method = "fixed width interval (negative exponential location)",
    is_satisfied = !unusable && n >= optimal_n,
    current_n = n,
    criterion_n = optimal_n,
    conf_level = conf_level,
    width = d,
    extra = list(variance_unusable = unusable),
    notes = if (unusable) .seq_notes_unusable(x)
  )
}
