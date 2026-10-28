#' Pilot Sample Size for a Fixed Width Confidence Interval for a Mean,
#' Distribution Free
#'
#' Returns the number of observations to collect before the stopping rule of
#' \code{\link{seq_fwd_dfr_mean}} is checked for the first time. That
#' procedure is the purely sequential fixed width confidence interval for
#' the mean of a population whose distribution is not specified
#' (Mukhopadhyay and de Silva, 2009).
#'
#' @details
#' If the population standard deviation \eqn{\sigma} were known, the sample
#' size that gives the interval \eqn{\bar{X}_n \pm d} approximately the
#' stated confidence level would be
#' \eqn{n^* = z_{\alpha/2}^2 \sigma^2 / d^2}. The pilot must be chosen
#' before anything is known about \eqn{\sigma}, so it uses only the known
#' constant \eqn{(z_{\alpha/2}/d)^2}:
#' \deqn{m = \max\{m_0,\ \lceil (z_{\alpha/2}/d)^{2/(1 + \gamma)} \rceil\}.}
#' With the default \eqn{\gamma = 1} this is
#' \eqn{\max\{m_0, \lceil z_{\alpha/2}/d \rceil\}}. Because of the
#' \eqn{1/n} term in the stopping rule of \code{seq_fwd_dfr_mean()}, that
#' rule can never be met with fewer than \eqn{z_{\alpha/2}/d} observations,
#' so nothing is lost by collecting that many at the outset.
#'
#' With \eqn{\gamma = 1} this is the pilot sample size of Kelley, Darku, and
#' Chattopadhyay (2018, Equation 21), whose full width \eqn{\omega} is
#' \eqn{2d}; those authors present it as following Mukhopadhyay and de Silva
#' (2009). The argument \code{gamma} and the default value of \code{m_0}
#' are the package's.
#'
#' @param d The half width of the confidence interval, in the units of the
#'   data: the interval reported is the sample mean plus and minus \code{d},
#'   so its full width is \eqn{2d}. A positive number.
#' @param conf_level The confidence level (default 0.95).
#' @param m_0 The smallest pilot sample size allowed. A whole number of at
#'   least 2; default 10. A larger value protects against stopping too early
#'   when the target is wide (see the Details of
#'   \code{\link{seq_fwd_dfr_mean}}); it costs at most \code{m_0}
#'   observations.
#' @param gamma Positive constant that sets how fast the pilot grows as the
#'   target narrows; see Details. Default 1.
#'
#' @return A length one integer giving the pilot sample size.
#'
#' @references
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
#' # Mean weekly hours of remote work, to within half an hour, 90% confidence
#' seq_fwd_dfr_mean_pilot(d = 0.5, conf_level = 0.90)
#' seq_fwd_dfr_mean_pilot(d = 0.1, conf_level = 0.90)
#'
#' @seealso \code{\link{seq_fwd_dfr_mean}}
#' @family seq_fwd_dfr_mean
#' @export
seq_fwd_dfr_mean_pilot <- function(d, conf_level = 0.95, m_0 = 10,
                                   gamma = 1) {
  .ts_check_width(d, "the confidence interval half width")
  alpha <- .ts_alpha(conf_level)
  z <- .aipe_z(conf_level)
  .seq_pilot((z / d)^2, m_0, gamma,
             advice = "Check 'd' (the half width of the confidence interval).")
}


#' Fixed Width Confidence Interval for a Mean, Distribution Free
#'
#' Given the data collected so far, checks the stopping rule of the purely
#' sequential fixed width confidence interval for a population mean when the
#' distribution is not specified (Mukhopadhyay and de Silva, 2009). The aim
#' is an interval of half width \code{d} (full width \eqn{2d}) around the
#' sample mean with approximately the stated confidence level. The sample
#' size needed depends on the unknown variance, so observations are
#' collected until the sample size reaches the size the current variance
#' estimate calls for.
#'
#' The function does not collect data and keeps no state. Call it after each
#' new observation (or batch), starting once
#' \code{\link{seq_fwd_dfr_mean_pilot}} observations are in hand, and stop
#' sampling the first time \code{attr(out, "is_satisfied")} is \code{TRUE}.
#'
#' @details
#' Sampling stops at the first \eqn{n} (at or beyond the pilot size) with
#' \deqn{n \ge \frac{z_{\alpha/2}^2}{d^2}\left(S_n^2 + \frac{1}{n}\right),}
#' where \eqn{S_n^2} is the sample variance. The \eqn{1/n} term keeps the
#' rule from being met on a variance estimate that happens to be very small
#' early on; it vanishes as the sample grows. The \code{optimal_n} row is
#' the right hand side rounded up, evaluated at the current data. The
#' reported interval is always \eqn{\bar{X}_n \pm d}; it has the stated
#' confidence level, approximately, only once the rule is met. The rule has
#' the form of the stopping rule of Kelley, Darku, and Chattopadhyay (2018,
#' Equation 22), with the sample variance as the estimate of the asymptotic
#' variance and their full width \eqn{\omega} equal to \eqn{2d}.
#'
#' \strong{The confidence level is asymptotic.} It rests on the central
#' limit theorem and on the consistency of \eqn{S_n^2}, and is attained in
#' the limit as \code{d} shrinks. For wide targets, and for skewed data even
#' at fairly narrow ones, coverage is below the stated level. In simulation
#' (standard deviation 2, nominal 0.95, the pilot of
#' \code{seq_fwd_dfr_mean_pilot()} with its defaults, 40,000 replications,
#' Monte Carlo standard error at most 0.002), coverage for normal data was
#' 0.935 at \code{d = 1} (where the required sample size
#' \eqn{n^* = z_{\alpha/2}^2\sigma^2/d^2} is about 15), 0.942 at
#' \code{d = 0.5} (\eqn{n^*} about 61), and 0.948 at \code{d = 0.25}
#' (\eqn{n^*} about 246). For exponential data with the same standard
#' deviation it was 0.960 at \code{d = 1} (where the pilot of 10 is already
#' most of the 15 required) and then far lower: 0.879 at \eqn{n^*} of 30,
#' and lowest, about 0.855, when \eqn{n^*} was between 40 and 61
#' (\code{d = 0.5}). It recovered slowly: 0.886 at \eqn{n^*} of 100, 0.934
#' at \code{d = 0.25} (\eqn{n^*} about 246), 0.946 at \code{d = 0.125}
#' (\eqn{n^*} about 983; 400,000 replications), and still 0.946 at
#' \code{d = 0.0625} (\eqn{n^*} about 3,934; 10,000 replications, standard
#' error 0.002). Skewness slows the approach to the stated level
#' considerably.
#'
#' \strong{The units of the data.} The \eqn{1/n} term is added to a
#' variance, so it is in the squared units of the data, and the procedure
#' depends on the units in which the data are recorded. It adds up to about
#' \eqn{1/\sigma^2} observations to the sample size that would otherwise be
#' needed: a negligible number when the standard deviation is 1 or more in
#' the units of the data, and a substantial one when it is much smaller. In
#' simulation with normal data, \eqn{\sigma = 0.25}, and \code{d = 0.1}
#' (\eqn{n^*} about 24), the average final sample size was 35 and coverage
#' was 0.98. The interval is then conservative, not wrong.
#'
#' \strong{Constant data.} If every observation is equal, or equal up to
#' rounding error (a range of no more than \eqn{10^{-10}} times the largest
#' absolute value, as with gain scores of 0.3 computed as posttest minus
#' pretest), the sample variance is zero, or nothing but rounding error, and
#' under the rule as printed the \eqn{1/n} term alone would eventually stop
#' the sampling. The function departs from the printed rule deliberately: no
#' variance has been estimated, so the rule is reported as not met,
#' \code{optimal_n} is \code{NA}, \code{attr(out, "variance_unusable")} is
#' \code{TRUE}, and a note says so (see \code{\link{smsd_check_methods}});
#' data whose real spread is that small beside their size should have a
#' common value subtracted first.
#' The limits are still reported (the sample mean minus and plus
#' \code{d}, which involve no variance), but they are then not an interval
#' to report: \code{summary()}, \code{glance()}, and \code{tidy()} return
#' them without the note, so test \code{attr(out, "variance_unusable")}
#' before using them.
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
#' @param d The half width of the confidence interval, in the units of the
#'   data: the interval reported is the sample mean plus and minus \code{d},
#'   so its full width is \eqn{2d}. A positive number.
#' @param conf_level The confidence level (default 0.95).
#' @param na.rm Logical; should \code{NA} values be removed first? Default
#'   \code{TRUE}. With \code{FALSE}, missing values are an error.
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} with
#'   columns \code{term} and \code{value} and these rows:
#'   \describe{
#'     \item{\code{lower_limit}}{The sample mean minus \code{d}.}
#'     \item{\code{mean}}{The sample mean, the point estimate.}
#'     \item{\code{upper_limit}}{The sample mean plus \code{d}.}
#'     \item{\code{optimal_n}}{The sample size the current data call for,
#'       \eqn{\lceil z_{\alpha/2}^2 (S_n^2 + 1/n) / d^2 \rceil} (\code{NA}
#'       when the data are constant).}
#'     \item{\code{current_n}}{The number of observations used.}
#'   }
#'   Attributes: \code{is_satisfied} (\code{TRUE} once \code{current_n} is at
#'   least \code{optimal_n}), \code{current_n}, \code{criterion_n} (equal to
#'   \code{optimal_n}), \code{conf_level}, \code{width} (the full width
#'   \eqn{2d}), \code{method}, \code{variance_unusable} (\code{TRUE} when the
#'   data are constant, so that no variance has been estimated), and
#'   \code{notes}.
#'
#' @references
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
#' # Mean weekly hours of remote work in a firm, to within half an hour
#' # with 90% confidence. Hours are right skewed; simulated here.
#' set.seed(113)
#' m <- seq_fwd_dfr_mean_pilot(d = 0.5, conf_level = 0.90)
#' hours <- rgamma(m, shape = 4, scale = 2)
#' out <- seq_fwd_dfr_mean(hours, d = 0.5, conf_level = 0.90)
#' out
#'
#' # Add one employee at a time until the stopping rule is met.
#' while (!attr(out, "is_satisfied")) {
#'   hours <- c(hours, rgamma(1, shape = 4, scale = 2))
#'   out <- seq_fwd_dfr_mean(hours, d = 0.5, conf_level = 0.90)
#' }
#' out
#'
#' @seealso \code{\link{seq_fwd_dfr_mean_pilot}};
#'   \code{\link{seq_fwd_norm_mean}} for the normal theory rule without the
#'   \eqn{1/n} term; \code{\link{ts_fwd_dfr_mean}} for the two stage
#'   procedure; \code{\link{seq_aipe_sm}} for a mean standardized by its
#'   standard deviation.
#' @family seq_fwd_dfr_mean
#' @export
seq_fwd_dfr_mean <- function(x, d, conf_level = 0.95, na.rm = TRUE) {
  .ts_check_width(d, "the confidence interval half width")
  alpha <- .ts_alpha(conf_level)
  x <- .mr_check_data(x, na.rm, arg = "x")
  n <- length(x)
  xbar <- mean(x)
  s2 <- stats::var(x)
  unusable <- .seq_zero_spread(s2, x)

  optimal_n <- if (unusable) {
    NA_real_
  } else {
    .seq_criterion(
      .aipe_z(conf_level)^2 * (s2 + 1 / n) / d^2,
      advice = paste("Check 'd' (the half width of the confidence",
                     "interval) and that it is in the units of 'x'."))
  }

  .seq_check(
    term = c("lower_limit", "mean", "upper_limit", "optimal_n", "current_n"),
    value = c(xbar - d, xbar, xbar + d, optimal_n, n),
    method = "fixed width interval (distribution free mean)",
    is_satisfied = !unusable && n >= optimal_n,
    current_n = n,
    criterion_n = optimal_n,
    conf_level = conf_level,
    width = 2 * d,
    extra = list(variance_unusable = unusable),
    notes = if (unusable) .seq_notes_unusable(x)
  )
}
