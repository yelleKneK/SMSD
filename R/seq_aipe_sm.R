#' Pilot Sample Size for Sequential AIPE: Standardized Mean
#'
#' Returns the number of observations to collect before the sequential
#' accuracy in parameter estimation (AIPE) stopping rule for the population
#' standardized mean \eqn{\mu/\sigma} can be evaluated for the first time
#' with \code{\link{seq_aipe_sm}}.
#'
#' The pilot size is
#' \deqn{m = \max\{m_0, \lceil 2 z_{\alpha/2} / \omega \rceil\},}
#' where \eqn{\omega} is \code{width}, \eqn{z_{\alpha/2}} is the upper
#' \eqn{\alpha/2} quantile of the standard normal distribution with
#' \eqn{\alpha} equal to one minus \code{conf_level}, and \eqn{m_0} is
#' \code{m_0} (Kelley, Darku, and Chattopadhyay, 2018, Equation 21).
#' Whatever the data, the stopping rule cannot be met with fewer than
#' \eqn{\lceil 2 z_{\alpha/2} / \omega \rceil} observations, which is why
#' the pilot is at least that large.
#'
#' @param width The target full width of the confidence interval for the
#'   population standardized mean (\eqn{\omega} in the cited article): the
#'   upper limit minus the lower limit, in standard deviation units. A
#'   single positive number.
#' @param conf_level The confidence level of the interval, a single number
#'   strictly between 0 and 1. Default 0.95.
#' @param m_0 The smallest pilot size allowed, a whole number of at least 4
#'   (the variance estimator needs four observations). Default 4. Raising it
#'   makes an unusable early variance estimate less likely; see
#'   \code{\link{seq_aipe_sm}}.
#'
#' @return A single integer: the pilot sample size.
#'
#' @author
#' Ken Kelley \email{kkelley@@nd.edu},
#' Francis Bilson Darku \email{fbilsond@@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @references
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#' parameter estimation for a general class of effect sizes: A sequential
#' approach. \emph{Psychological Methods}, \emph{23}(2), 226--243.
#' \doi{10.1037/met0000127}
#'
#' @examples
#' # A 95% interval for a standardized mean no wider than 0.20
#' seq_aipe_sm_pilot(width = 0.20)
#'
#' # A wider target at 90% confidence, with at least 15 in the pilot
#' seq_aipe_sm_pilot(width = 0.50, conf_level = 0.90, m_0 = 15)
#'
#' @seealso \code{\link{seq_aipe_sm}}
#' @family seq_aipe_sm
#' @export
seq_aipe_sm_pilot <- function(width, conf_level = 0.95, m_0 = 4) {
  .aipe_pilot(width, conf_level, m_0)
}


#' Sequential AIPE Check for the Standardized Mean
#'
#' Checks, with the data collected so far, whether enough observations are
#' in hand for a confidence interval for the population standardized mean
#' \eqn{\mu/\sigma} (the mean in standard deviation units; for difference
#' scores, the standardized mean change) that is no wider than a chosen
#' width, without assuming normality and without a planning value. The
#' standardized mean is one of the effect sizes covered by the sequential
#' accuracy in parameter estimation (AIPE) procedure of Kelley, Darku, and
#' Chattopadhyay (2018, Example 3).
#'
#' The function evaluates the stopping rule once and returns. The user
#' collects the pilot sample (see \code{\link{seq_aipe_sm_pilot}}), calls
#' \code{seq_aipe_sm()}, and, while the rule is not met, adds one or more
#' observations and calls it again. The function prints nothing and keeps no
#' state.
#'
#' \strong{Estimate.} The standardized mean is estimated by
#' \eqn{\bar{X}_n / s_n}. Any real mean is allowed (the estimate is negative
#' when the sample mean is); only a sample with no variation at all cannot
#' be used.
#'
#' \strong{Stopping rule.} Sampling stops at the first \eqn{n}, at or after
#' the pilot size, with
#' \deqn{n \ge \left\lceil \frac{4 z_{\alpha/2}^2}{\omega^2}
#'   \left(\hat{\xi}_n^2 + \frac{1}{n}\right) \right\rceil}
#' (their Equation 22), and the reported interval is
#' \eqn{\bar{X}_n / s_n \pm z_{\alpha/2} \hat{\xi}_n / \sqrt{n}} (their
#' Equation 23). Whenever the rule is met, the width of that interval is at
#' most \code{width}. The article names the standardized mean as a member of
#' its class of effect sizes but does not print its asymptotic variance. It
#' follows from the article's Theorem 1 by the delta method:
#' \deqn{\xi^2 = 1 - \frac{\mu \mu_3}{\sigma^4}
#'   + \frac{\mu^2 (\mu_4 - \sigma^4)}{4 \sigma^6},}
#' where \eqn{\mu_3} and \eqn{\mu_4} are the third and fourth central
#' moments (for normal data this is \eqn{1 + (\mu/\sigma)^2/2}). The
#' estimate \eqn{\hat{\xi}_n^2} substitutes the sample mean and variance and
#' the unbiased central moment estimators of the article's Equations 78 and
#' 79.
#'
#' \strong{What the guarantee is.} The width guarantee is exact. The
#' coverage guarantee is asymptotic: coverage approaches \code{conf_level}
#' as \code{width} shrinks (their Theorem 2). In simulated studies with a
#' 95\% interval and one observation added at each stage, coverage for
#' normal data with \eqn{\mu/\sigma = 0.4} was 0.949 at \code{width = 0.4}
#' (6,000 studies, about 105 observations at the end) and 0.956 at
#' \code{width = 0.2} (2,000 studies, about 415 observations). For right
#' skewed data (a standardized chi square distribution with 3 degrees of
#' freedom, shifted so that \eqn{\mu/\sigma = 1}) it was 0.893 at
#' \code{width = 0.4} and 0.940 at \code{width = 0.2}: with skewed data the
#' variance estimate runs low in small samples, so the procedure stops
#' early on average (at 62 observations where the known variance would
#' require 84, and at 289 where it would require 334). The variance
#' estimate was zero or negative in fewer than 1\% of the pilot samples, and
#' positive but below \eqn{n^{-3}}, so that the bound was used and the rule
#' met with an interval far narrower than requested, in 2 of the 6,000
#' studies with skewed data at \code{width = 0.4} and in no other. Use
#' \code{\link{seq_aipe_sm_sensitivity}} to study a planned study before
#' collecting data.
#'
#' @inheritSection seq_aipe_smd Estimates That Are Not Yet Usable
#'
#' @param x A numeric vector: the observations collected so far, at least 4
#'   of them. Missing or infinite values are an error.
#' @param width The target full width of the confidence interval for
#'   \eqn{\mu/\sigma} (\eqn{\omega} in the cited article): the upper limit
#'   minus the lower limit, in standard deviation units. A single positive
#'   number.
#' @param conf_level The confidence level of the interval, a single number
#'   strictly between 0 and 1. Default 0.95.
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} (see
#'   \code{\link{smsd_check_methods}} for its \code{print()},
#'   \code{summary()}, and \code{tidy()} methods) with one row per quantity
#'   and columns \code{term} and \code{value}. The rows are, in order:
#'   \describe{
#'     \item{\code{lower_limit}}{the lower confidence limit (\code{NA}
#'       when the variance estimate is not yet usable or the estimate
#'       cannot be computed; see the section below);}
#'     \item{\code{sm}}{the estimated standardized mean (\code{NA} when the
#'       sample has no variation);}
#'     \item{\code{upper_limit}}{the upper confidence limit (\code{NA}
#'       whenever \code{lower_limit} is);}
#'     \item{\code{xi_squared}}{the estimated asymptotic variance
#'       \eqn{\hat{\xi}_n^2}, at least \eqn{n^{-3}} (the bound itself
#'       when the estimate is not yet usable);}
#'     \item{\code{criterion_n}}{the sample size the stopping rule currently
#'       requires (at least \code{current_n + 1} when the variance estimate
#'       is not yet usable or the estimate cannot be computed);}
#'     \item{\code{current_n}}{the current sample size.}
#'   }
#'   The object also carries these attributes:
#'   \describe{
#'     \item{\code{is_satisfied}}{\code{TRUE} when the stopping rule is
#'       met, otherwise \code{FALSE}; never \code{NA}, and always
#'       identical to \code{current_n >= criterion_n};}
#'     \item{\code{method}}{the string \code{"standardized mean"};}
#'     \item{\code{current_n}, \code{criterion_n}}{the values in the
#'       rows of those names;}
#'     \item{\code{conf_level}, \code{width}}{the values supplied;}
#'     \item{\code{variance_unusable}}{\code{TRUE} when the raw variance
#'       estimate is zero, negative, or not finite (see the section below),
#'       otherwise \code{FALSE};}
#'     \item{\code{notes}}{a character vector, empty unless something
#'       needs explaining.}
#'   }
#'
#' @author
#' Ken Kelley \email{kkelley@@nd.edu},
#' Francis Bilson Darku \email{fbilsond@@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @references
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#' parameter estimation for a general class of effect sizes: A sequential
#' approach. \emph{Psychological Methods}, \emph{23}(2), 226--243.
#' \doi{10.1037/met0000127}
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#' for statistical power and accuracy in parameter estimation. \emph{Annual
#' Review of Psychology}, \emph{59}, 537--563.
#' \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' @examples
#' # By how many standard deviations does a reading score change from fall
#' # to spring? The data are the students' change scores, and the aim is a
#' # 95% confidence interval for the standardized mean change that is no
#' # wider than 0.50.
#' set.seed(113)
#' m <- seq_aipe_sm_pilot(width = 0.50)
#' change <- rnorm(m, mean = 4, sd = 10)
#' out <- seq_aipe_sm(change, width = 0.50)
#' out
#'
#' # The rule is not met at the pilot. Add one student at a time and check
#' # again, until it is.
#' while (!attr(out, "is_satisfied")) {
#'   change <- c(change, rnorm(1, mean = 4, sd = 10))
#'   out <- seq_aipe_sm(change, width = 0.50)
#' }
#' out
#'
#' @seealso \code{\link{seq_aipe_sm_pilot}} for the pilot size and
#'   \code{\link{seq_aipe_sm_sensitivity}} to simulate the procedure before
#'   collecting data. In \pkg{DMAR}: \code{\link[DMAR]{ss_aipe_sm}}, the
#'   fixed sample size analog for when a planning value is available, and
#'   \code{\link[DMAR]{ci_sm}}.
#' @family seq_aipe_sm
#' @export
seq_aipe_sm <- function(x, width, conf_level = 0.95) {
  x <- .aipe_check_data(x, "x")
  width <- .check_width(width)
  conf_level <- .check_conf_level(conf_level)
  n <- length(x)
  # The standardized mean does not depend on the scale of the data, which
  # are divided by a power of two of their magnitude so that no sum of
  # squares can overflow or underflow; see .aipe_scale().
  x <- x / .aipe_scale(x)

  # A sample with no variation has no standardized mean (a division by
  # zero); it is reported as not estimable.
  std_mean <- if (.is_constant(x)) NA_real_ else mean(x) / stats::sd(x)
  raw <- .var_sm(x, floor = FALSE)

  .aipe_result(estimate = std_mean, raw = raw, n = n, width = width,
               conf_level = conf_level, term = "sm",
               method = "standardized mean")
}
