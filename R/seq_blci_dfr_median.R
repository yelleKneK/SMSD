#' Pilot Sample Size for a Bounded Length Confidence Interval for a Median,
#' Distribution Free
#'
#' Returns the number of observations to collect before the stopping rule of
#' \code{\link{seq_blci_dfr_median}} is checked for the first time. That
#' procedure is the purely sequential bounded length confidence interval for
#' the median of a population whose distribution is not specified
#' (Mukhopadhyay and de Silva, 2009).
#'
#' @details
#' The pilot is the largest of three numbers:
#' \deqn{m = \max\{n_{\min},\ m_0,\
#'   \lceil (z_{\alpha/2}/d)^{2/(1 + \gamma)} \rceil\}.}
#' The first, \eqn{n_{\min}}, is the smallest sample size at which an
#' interval formed from order statistics can reach the confidence level at
#' all. The widest such interval, from the smallest to the largest
#' observation, misses the median only when every observation falls on the
#' same side of it, so its coverage is \eqn{1 - 2^{1-n}}; \eqn{n_{\min}} is
#' the smallest \eqn{n} for which this is at least \code{conf_level} (5, 6,
#' and 8 observations at 90\%, 95\%, and 99\%).
#'
#' The third is of lower order than the sample size the procedure will need,
#' which is approximately
#' \eqn{n^* = (z_{\alpha/2}/d)^2 / \{4 f(\theta)^2\}} for a population with
#' density \eqn{f} at its median \eqn{\theta}. The density is unknown before
#' sampling, so the pilot uses only the known constant
#' \eqn{(z_{\alpha/2}/d)^2}; with the default \eqn{\gamma = 1} the third
#' term is \eqn{\lceil z_{\alpha/2}/d \rceil}. A pilot equal to
#' \eqn{(z_{\alpha/2}/d)^2} itself would be several times too large for a
#' population that is concentrated near its median in the units of the
#' data.
#'
#' This is the package's pilot rule, of the same form as the pilot sample
#' size of Kelley, Darku, and Chattopadhyay (2018, Equation 21). It should
#' not be read as the pilot size prescribed by Mukhopadhyay and de Silva
#' (2009).
#'
#' @param d Half of the largest interval length allowed, in the units of the
#'   data: sampling continues until the length of the interval is at most
#'   \eqn{2d}. A positive number.
#' @param conf_level The confidence level (default 0.95).
#' @param m_0 The smallest pilot sample size allowed, in addition to the
#'   minimum the confidence level itself requires. A whole number of at
#'   least 2; default 10.
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
#' # Median household income (thousands of dollars), interval no longer
#' # than 10 thousand
#' seq_blci_dfr_median_pilot(d = 5, conf_level = 0.95)
#' # A tight target and a high confidence level
#' seq_blci_dfr_median_pilot(d = 0.1, conf_level = 0.99)
#'
#' @seealso \code{\link{seq_blci_dfr_median}}
#' @family seq_blci_dfr_median
#' @export
seq_blci_dfr_median_pilot <- function(d, conf_level = 0.95, m_0 = 10,
                                      gamma = 1) {
  .ts_check_width(d, "half of the largest confidence interval length allowed")
  alpha <- .ts_alpha(conf_level)
  z <- .aipe_z(conf_level)
  max(.seq_blci_min_n(conf_level),
      .seq_pilot((z / d)^2, m_0, gamma,
                 advice = paste("Check 'd' (half of the largest confidence",
                                "interval length allowed).")))
}


#' Bounded Length Confidence Interval for a Median, Distribution Free
#'
#' Given the data collected so far, checks the stopping rule of the purely
#' sequential bounded length confidence interval for a population median
#' when the distribution is not specified (Mukhopadhyay and de Silva, 2009).
#' The interval runs between two order statistics of the sample, and
#' observations are collected until its length is at most \eqn{2d}. The
#' quantity estimated is the population median, not the mean.
#'
#' The function does not collect data and keeps no state. Call it after each
#' new observation (or batch), starting once
#' \code{\link{seq_blci_dfr_median_pilot}} observations are in hand, and
#' stop sampling the first time \code{attr(out, "is_satisfied")} is
#' \code{TRUE}.
#'
#' @details
#' With the \eqn{n} observations sorted as
#' \eqn{X_{(1)} \le \dots \le X_{(n)}}, the interval is
#' \eqn{[X_{(b)},\ X_{(n-b+1)}]}. For a continuous distribution it covers
#' the median exactly when the number of observations below the median, a
#' binomial count with \eqn{n} trials and probability 1/2, lies between
#' \eqn{b} and \eqn{n - b}. The rank \eqn{b} starts from the continuity
#' corrected normal approximation to that binomial,
#' \eqn{b = \max\{1, \lceil (n - z_{\alpha/2}\sqrt{n} - 1)/2 \rceil\}},
#' and is lowered if necessary until the exact binomial coverage is at least
#' \code{conf_level}. (At confidence levels of 0.80 and above the
#' approximation needed no lowering for any \eqn{n} up to 3,000; at lower
#' levels it is occasionally short by a small amount.) At the 90\%, 95\%,
#' and 99\% levels this is the largest rank with that property, that is,
#' the narrowest such interval, at more than 97\% of the sample sizes up to
#' 2,000, and one rank lower at the rest. Sampling stops at the first
#' \eqn{n} (at or beyond the pilot size) at which
#' \eqn{X_{(n-b+1)} - X_{(b)} \le 2d}. Unlike the fixed width intervals for
#' a mean, this interval is not symmetric about the point estimate and its
#' length is not fixed in advance; \eqn{2d} is an upper bound on it.
#'
#' \strong{Minimum sample size.} Below a certain sample size no pair of
#' order statistics can reach the confidence level: even the interval from
#' the smallest to the largest observation has coverage \eqn{1 - 2^{1-n}}
#' (0.50, 0.75, 0.875, and 0.9375 for \eqn{n} = 2 to 5). With fewer than
#' that many observations (6 at 95\%), the rule is reported as not met, the
#' limits and the interval length are \code{NA}, and a note says how many
#' observations are needed. \code{attr(out, "min_n")} holds that number.
#'
#' \strong{The confidence level is asymptotic, and the shortfall is not
#' small.} At any fixed \eqn{n} the interval has at least the stated
#' coverage. The sequential procedure, however, stops when the interval
#' happens to be short, and short intervals miss the median more often, so
#' its coverage is below the stated level and approaches it only slowly as
#' \code{d} shrinks. In simulation (the pilot of
#' \code{seq_blci_dfr_median_pilot()} with its defaults, 40,000
#' replications in each condition, Monte Carlo standard error at most
#' 0.002), with \code{d} chosen so that the sample size needed (\eqn{n^*}
#' in the Details of \code{\link{seq_blci_dfr_median_pilot}}) was 25, 50,
#' 100, 200, 400, and 800, a nominal 0.95 interval covered the median with
#' probability 0.936, 0.928, 0.926, 0.929, 0.935, and 0.938 for standard
#' normal data, 0.947, 0.925, 0.918, 0.925, 0.929, and 0.937 for unit
#' exponential data, and 0.947, 0.937, 0.936, 0.936, 0.938, and 0.937 for
#' standard Cauchy data. The average final sample size was between 0.94 and
#' 1.05 times \eqn{n^*} once \eqn{n^*} was 100 or more. Expect coverage of
#' roughly 0.92 to 0.94 from a nominal 0.95 interval, unless the pilot is
#' already a large part of the sample size needed. The shortfall is larger
#' at a lower level and smaller at a higher one: for the same three
#' distributions with \eqn{n^*} of 50, 100, 200, and 800, a nominal 0.90
#' gave between 0.842 and 0.879 and a nominal 0.99 between 0.983 and 0.989.
#' Because the shortfall is not negligible, a check whose rule is met
#' carries a note saying so, which \code{print()} shows.
#'
#' \strong{If the stated level must be met,} ask for more than is needed.
#' With \code{d} held at the value for which the 0.95 procedure needs
#' about 100 and about 400 observations (normal data; 20,000 replications
#' at the first target and 6,000 at the second, standard error at most
#' 0.003), a nominal 0.96 was not enough (coverage 0.941 and 0.948). A
#' nominal 0.97 gave 0.955 and 0.962, with an average final sample size of
#' 1.20 and 1.17 times those 100 and 400 observations, and a nominal 0.975
#' gave 0.961 and 0.970, at 1.29 and 1.25 times as many. For exponential
#' data at the first of the two targets, a nominal 0.97 gave 0.952.
#'
#' \strong{Ties and constant data.} Ties are permitted. With heavily tied
#' data the two order statistics can coincide, giving an interval of length
#' zero; that is a valid (if uninformative about precision) interval for
#' the median of a discrete population, and the rule is then met. (This
#' procedure estimates no variance, so the package's rule for a variance
#' estimate that cannot be used does not arise here.)
#'
#' @param x Numeric vector of the observations collected so far (at least
#'   1; an interval is reported once there are enough for the confidence
#'   level, see Details).
#' @param d Half of the largest interval length allowed, in the units of the
#'   data: sampling continues until the length of the interval is at most
#'   \eqn{2d}. A positive number.
#' @param conf_level The confidence level (default 0.95).
#' @param na.rm Logical; should \code{NA} values be removed first? Default
#'   \code{TRUE}. With \code{FALSE}, missing values are an error.
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} with
#'   columns \code{term} and \code{value} and these rows:
#'   \describe{
#'     \item{\code{lower_limit}}{The order statistic \eqn{X_{(b)}}.}
#'     \item{\code{median}}{The sample median, the point estimate.}
#'     \item{\code{upper_limit}}{The order statistic \eqn{X_{(n-b+1)}}.}
#'     \item{\code{interval_length}}{The upper limit minus the lower limit.}
#'     \item{\code{current_n}}{The number of observations used.}
#'   }
#'   The two limits and the length are \code{NA} when there are too few
#'   observations for the confidence level. Attributes: \code{is_satisfied}
#'   (\code{TRUE} when the interval length is at most \eqn{2d});
#'   \code{current_n}; \code{criterion_n}, which is \code{NA} because the
#'   rule is stated in terms of the interval length and not a target sample
#'   size; \code{conf_level}; \code{width} (the bound \eqn{2d});
#'   \code{method}; \code{min_n} (the smallest sample size at which the
#'   confidence level can be reached); \code{ranks} (the ranks \eqn{b} and
#'   \eqn{n - b + 1} of the two order statistics); and \code{notes} (when
#'   the rule is met, the reminder that a sequential procedure covers the
#'   median less often than stated; when there are too few observations,
#'   how many are needed; otherwise empty).
#'
#' @references
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#'   Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @examples
#' # Median household income in a district (thousands of dollars), with an
#' # interval no longer than 10 thousand, 95% confidence. Incomes are right
#' # skewed; simulated here from a lognormal distribution.
#' set.seed(113)
#' m <- seq_blci_dfr_median_pilot(d = 5, conf_level = 0.95)
#' income <- rlnorm(m, meanlog = log(60), sdlog = 0.5)
#' out <- seq_blci_dfr_median(income, d = 5, conf_level = 0.95)
#' out
#'
#' # Add one household at a time until the stopping rule is met. The note
#' # printed with the result is a reminder that a completed sequential
#' # procedure covers the median somewhat less often than the stated 95%
#' # (see Details).
#' while (!attr(out, "is_satisfied")) {
#'   income <- c(income, rlnorm(1, meanlog = log(60), sdlog = 0.5))
#'   out <- seq_blci_dfr_median(income, d = 5, conf_level = 0.95)
#' }
#' out
#'
#' # Too few observations for 95% confidence: not met, with a note.
#' seq_blci_dfr_median(income[1:4], d = 5, conf_level = 0.95)
#'
#' @seealso \code{\link{seq_blci_dfr_median_pilot}};
#'   \code{\link{seq_fwd_dfr_mean}} for a fixed width interval for a mean.
#' @family seq_blci_dfr_median
#' @export
seq_blci_dfr_median <- function(x, d, conf_level = 0.95, na.rm = TRUE) {
  .ts_check_width(d, "half of the largest confidence interval length allowed")
  .ts_alpha(conf_level)
  x <- .mr_check_data(x, na.rm, min_n = 1L, arg = "x")
  n <- length(x)
  min_n <- .seq_blci_min_n(conf_level)
  b <- .seq_blci_rank(n, conf_level)
  too_few <- b < 1L

  if (too_few) {
    lower <- upper <- NA_real_
    ranks <- c(lower = NA_integer_, upper = NA_integer_)
  } else {
    sorted <- sort(x)
    lower <- sorted[b]
    upper <- sorted[n - b + 1L]
    ranks <- c(lower = b, upper = n - b + 1L)
  }
  interval_length <- upper - lower
  is_satisfied <- !too_few && interval_length <= 2 * d

  .seq_check(
    term = c("lower_limit", "median", "upper_limit", "interval_length",
             "current_n"),
    value = c(lower, stats::median(x), upper, interval_length, n),
    method = "bounded length confidence interval (distribution free)",
    is_satisfied = is_satisfied,
    current_n = n,
    criterion_n = NA_real_,
    conf_level = conf_level,
    width = 2 * d,
    extra = list(min_n = min_n, ranks = ranks),
    notes = if (too_few) {
      .note_blci_too_few(n, min_n, conf_level)
    } else if (is_satisfied) {
      .note_blci_sequential()
    }
  )
}
