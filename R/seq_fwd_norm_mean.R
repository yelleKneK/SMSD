#' Pilot Sample Size for a Fixed Width Confidence Interval for a Normal Mean
#'
#' Returns the number of observations to collect before the stopping rule of
#' \code{\link{seq_fwd_norm_mean}} is checked for the first time. That
#' procedure is the purely sequential fixed width confidence interval for
#' the mean of a normal population with unknown variance (Mukhopadhyay and
#' de Silva, 2009).
#'
#' @details
#' If the population standard deviation \eqn{\sigma} were known, the sample
#' size that gives the interval \eqn{\bar{X}_n \pm d} the stated confidence
#' level would be \eqn{n^* = z_{\alpha/2}^2 \sigma^2 / d^2}
#' (\code{\link{seq_fwd_norm_mean_optimal}} returns it, rounded up to a
#' whole number, for a given variance). The pilot must be chosen before
#' anything is known about \eqn{\sigma}, so it uses only the known constant
#' \eqn{(z_{\alpha/2}/d)^2}:
#' \deqn{m = \max\{m_0,\ \lceil (z_{\alpha/2}/d)^{2/(1 + \gamma)} \rceil\}.}
#' With the default \eqn{\gamma = 1} this is
#' \eqn{\max\{m_0, \lceil z_{\alpha/2}/d \rceil\}}, which grows more slowly
#' than \eqn{n^*} as \eqn{d} shrinks. It therefore does not presume a value
#' of \eqn{\sigma}: a pilot equal to \eqn{n^*} evaluated at \eqn{\sigma = 1}
#' would be far too large whenever the standard deviation is well below 1 in
#' the units of the data.
#'
#' This is the package's pilot rule. It has the form of the pilot sample
#' size of Kelley, Darku, and Chattopadhyay (2018, Equation 21), whose full
#' width \eqn{\omega} is \eqn{2d}. It should not be read as the pilot size
#' prescribed by Mukhopadhyay and de Silva (2009).
#'
#' @param d The half width of the confidence interval, in the units of the
#'   data: the interval reported is the sample mean plus and minus \code{d},
#'   so its full width is \eqn{2d}. A positive number.
#' @param conf_level The confidence level (default 0.95).
#' @param m_0 The smallest pilot sample size allowed. A whole number of at
#'   least 2; default 10. A larger value protects against stopping too early
#'   on an unluckily small variance estimate when the target is wide (see
#'   the Details of \code{\link{seq_fwd_norm_mean}}); it costs at most
#'   \code{m_0} observations.
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
#' # Mean score on a 100 point exam, to within 2 points with 95% confidence
#' seq_fwd_norm_mean_pilot(d = 2, conf_level = 0.95)
#' # A narrower target calls for a larger pilot
#' seq_fwd_norm_mean_pilot(d = 0.1, conf_level = 0.95)
#'
#' @seealso \code{\link{seq_fwd_norm_mean}};
#'   \code{\link{seq_fwd_norm_mean_optimal}} for the sample size when the
#'   variance is known.
#' @family seq_fwd_norm_mean
#' @export
seq_fwd_norm_mean_pilot <- function(d, conf_level = 0.95, m_0 = 10,
                                    gamma = 1) {
  .ts_check_width(d, "the confidence interval half width")
  alpha <- .ts_alpha(conf_level)
  z <- .aipe_z(conf_level)
  .seq_pilot((z / d)^2, m_0, gamma,
             advice = "Check 'd' (the half width of the confidence interval).")
}


#' Fixed Width Confidence Interval for a Normal Mean
#'
#' Given the data collected so far, checks the stopping rule of the purely
#' sequential fixed width confidence interval for the mean of a normal
#' population with unknown variance (Mukhopadhyay and de Silva, 2009). The
#' aim is an interval of half width \code{d} (full width \eqn{2d}) around
#' the sample mean with the stated confidence level. No fixed sample size
#' can promise that, because the size needed depends on the unknown
#' variance, so observations are collected until the sample size reaches the
#' size the current variance estimate calls for.
#'
#' The function does not collect data and keeps no state. Call it after each
#' new observation (or batch), starting once
#' \code{\link{seq_fwd_norm_mean_pilot}} observations are in hand, and stop
#' sampling the first time \code{attr(out, "is_satisfied")} is \code{TRUE}.
#'
#' @details
#' With the population standard deviation \eqn{\sigma} known, the interval
#' \eqn{\bar{X}_n \pm d} has confidence level \eqn{1 - \alpha} when
#' \eqn{n \ge n^* = z_{\alpha/2}^2 \sigma^2 / d^2}
#' (\code{\link{seq_fwd_norm_mean_optimal}} returns that sample size,
#' rounded up to a whole number, for a given variance). The sequential rule
#' replaces \eqn{\sigma^2} by the sample variance \eqn{S_n^2}: sampling
#' stops at the first \eqn{n} (at or beyond the pilot size) with
#' \deqn{n \ge z_{\alpha/2}^2 S_n^2 / d^2 .}
#' The \code{optimal_n} row is the right hand side rounded up, evaluated at
#' the current data. The reported interval is always
#' \eqn{\bar{X}_n \pm d}; it has the stated confidence level only once the
#' rule is met.
#'
#' \strong{The confidence level is asymptotic.} It is attained in the limit
#' as \code{d} shrinks (so that the required sample size grows), not at
#' every \code{d}. When the target is wide relative to the standard
#' deviation, sampling tends to stop early on an underestimate of the
#' variance, and coverage falls below the stated level. In simulation
#' (normal data with standard deviation 2, nominal 0.95, the pilot of
#' \code{seq_fwd_norm_mean_pilot()} with its defaults, which was 10 at
#' every target used, 40,000 replications, Monte Carlo standard error
#' 0.001), coverage was 0.933 at \code{d = 1} (where the required sample
#' size \eqn{n^*} is about 15). It was lowest, about 0.922, when
#' \eqn{n^*} was between 24 and 28, and rose from there: 0.930 at
#' \eqn{n^*} of 40, 0.938 at \code{d = 0.5} (\eqn{n^*} about 61), 0.943
#' at \eqn{n^*} of 100, and 0.949 at \code{d = 0.25} (\eqn{n^*} about
#' 246; 400,000 replications). The pilot matters: with \code{m_0 = 4} in
#' place of the default 10, coverage at \code{d = 1} was 0.880.
#' \code{\link{ts_fwd_norm_mean}} is the two stage alternative, whose
#' confidence level is guaranteed for normal data at the price of a larger
#' expected sample size.
#'
#' \strong{Constant data.} If every observation is equal, or equal up to
#' rounding error (a range of no more than \eqn{10^{-10}} times the largest
#' absolute value, as with gain scores of 0.3 computed as posttest minus
#' pretest), the sample variance is zero, or nothing but rounding error, and
#' the rule as printed would be met at once. The function departs from the
#' printed rule deliberately: no variance has been estimated, so the rule is
#' reported as not met, \code{optimal_n} is \code{NA},
#' \code{attr(out, "variance_unusable")} is \code{TRUE}, and a note says so
#' (see \code{\link{smsd_check_methods}}); data whose real spread is that
#' small beside their size should have a common value subtracted first.
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
#'     \item{\code{optimal_n}}{The sample size the current variance estimate
#'       calls for, \eqn{\lceil z_{\alpha/2}^2 S_n^2 / d^2 \rceil}
#'       (\code{NA} when the data are constant).}
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
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#'   Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # Mean score on a 100 point exam, to within 2 points with 95% confidence.
#' # Scores are simulated here with a standard deviation of 12.
#' set.seed(113)
#' m <- seq_fwd_norm_mean_pilot(d = 2, conf_level = 0.95)
#' scores <- rnorm(m, mean = 70, sd = 12)
#' out <- seq_fwd_norm_mean(scores, d = 2, conf_level = 0.95)
#' out
#'
#' # Add one student at a time until the stopping rule is met.
#' while (!attr(out, "is_satisfied")) {
#'   scores <- c(scores, rnorm(1, mean = 70, sd = 12))
#'   out <- seq_fwd_norm_mean(scores, d = 2, conf_level = 0.95)
#' }
#' out
#'
#' @seealso \code{\link{seq_fwd_norm_mean_pilot}};
#'   \code{\link{seq_fwd_norm_mean_optimal}} for the sample size when the
#'   variance is known, the benchmark this procedure aims for;
#'   \code{\link{seq_fwd_dfr_mean}}, which makes no normality assumption;
#'   \code{\link{ts_fwd_norm_mean}} for the two stage procedure.
#' @family seq_fwd_norm_mean
#' @export
seq_fwd_norm_mean <- function(x, d, conf_level = 0.95, na.rm = TRUE) {
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
      .aipe_z(conf_level)^2 * s2 / d^2,
      advice = paste("Check 'd' (the half width of the confidence",
                     "interval) and that it is in the units of 'x'."))
  }

  .seq_check(
    term = c("lower_limit", "mean", "upper_limit", "optimal_n", "current_n"),
    value = c(xbar - d, xbar, xbar + d, optimal_n, n),
    method = "fixed width interval (normal mean)",
    is_satisfied = !unusable && n >= optimal_n,
    current_n = n,
    criterion_n = optimal_n,
    conf_level = conf_level,
    width = 2 * d,
    extra = list(variance_unusable = unusable),
    notes = if (unusable) .seq_notes_unusable(x)
  )
}
