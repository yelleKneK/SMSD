#' Pilot Sample Size for Sequential AIPE: Simple Linear Regression Slope
#'
#' Returns the number of paired observations to collect before the
#' sequential accuracy in parameter estimation (AIPE) stopping rule for the
#' slope of a simple linear regression,
#' \eqn{\beta_1 = \sigma_{XY}/\sigma_X^2}, can be evaluated for the first
#' time with \code{\link{seq_aipe_reg_coef}}.
#'
#' The pilot size is
#' \deqn{m = \max\{m_0, \lceil 2 z_{\alpha/2} / \omega \rceil\},}
#' where \eqn{\omega} is \code{width}, \eqn{z_{\alpha/2}} is the upper
#' \eqn{\alpha/2} quantile of the standard normal distribution with
#' \eqn{\alpha} equal to one minus \code{conf_level}, and \eqn{m_0} is
#' \code{m_0} (Kelley, Darku, and Chattopadhyay, 2018, Equation 21).
#' Whatever the data, the stopping rule cannot be met with fewer than
#' \eqn{\lceil 2 z_{\alpha/2} / \omega \rceil} pairs, which is why the pilot
#' is at least that large. Because the slope is not scale free, the pilot
#' size depends on the units of \code{width}: see
#' \code{\link{seq_aipe_reg_coef}}.
#'
#' @param width The target full width of the confidence interval for the
#'   slope (\eqn{\omega} in the cited article): the upper limit minus the
#'   lower limit, in units of the outcome per unit of the predictor. A
#'   single positive number.
#' @param conf_level The confidence level of the interval, a single number
#'   strictly between 0 and 1. Default 0.95.
#' @param m_0 The smallest pilot size allowed, a whole number of at least 4
#'   (the variance estimator needs four pairs). Default 4. Raising it makes
#'   an unusable early variance estimate less likely; see
#'   \code{\link{seq_aipe_reg_coef}}.
#'
#' @return A single integer: the pilot number of pairs.
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
#' # A 95% interval for a slope no wider than 0.20
#' seq_aipe_reg_coef_pilot(width = 0.20)
#'
#' # The same target with at least 30 pairs in the pilot
#' seq_aipe_reg_coef_pilot(width = 0.20, m_0 = 30)
#'
#' @seealso \code{\link{seq_aipe_reg_coef}}
#' @family seq_aipe_reg_coef
#' @export
seq_aipe_reg_coef_pilot <- function(width, conf_level = 0.95, m_0 = 4) {
  .aipe_pilot(width, conf_level, m_0)
}


#' Sequential AIPE Check for the Simple Linear Regression Slope
#'
#' Checks, with the paired data collected so far, whether enough pairs are
#' in hand for a confidence interval for the slope of the simple linear
#' regression of \code{y} on \code{x},
#' \eqn{\beta_1 = \sigma_{XY} / \sigma_X^2}, that is no wider than a chosen
#' width, without assuming normal or constant variance errors and without
#' planning values. This is the sequential accuracy in parameter estimation
#' (AIPE) procedure of Kelley, Darku, and Chattopadhyay (2018).
#'
#' The function evaluates the stopping rule once and returns. The user
#' collects the pilot sample (see \code{\link{seq_aipe_reg_coef_pilot}}),
#' calls \code{seq_aipe_reg_coef()}, and, while the rule is not met, adds
#' one or more pairs and calls it again. The function prints nothing and
#' keeps no state.
#'
#' \strong{Estimate.} The slope is the ordinary least squares estimate
#' \eqn{b_{1n} = s_{XYn} / s_{Xn}^2} (their Equation 35), as returned by
#' \code{\link{slr_slope}}.
#'
#' \strong{Stopping rule.} Sampling stops at the first \eqn{n}, at or after
#' the pilot size, with
#' \deqn{n \ge \left\lceil \frac{4 z_{\alpha/2}^2}{\omega^2}
#'   \left(\hat{\xi}_n^2 + \frac{1}{n}\right) \right\rceil}
#' (their Equation 22), where \eqn{\hat{\xi}_n^2} is the estimator in their
#' Equations 84 and 85 of the asymptotic variance of
#' \eqn{\sqrt{n}(b_{1n} - \beta_1)} given in their Equation 37, built from
#' estimates of the joint central moments of \code{x} and \code{y}. The
#' reported interval is \eqn{b_{1n} \pm z_{\alpha/2} \hat{\xi}_n / \sqrt{n}}
#' (their Equation 38). Whenever the rule is met, the width of that interval
#' is at most \code{width}.
#'
#' \strong{Units.} Unlike a standardized effect size, the slope and
#' \code{width} are in units of \code{y} per unit of \code{x}, while the
#' correction term \eqn{1/n} in the stopping rule and the lower bound
#' \eqn{n^{-3}} on the variance estimate are pure numbers, as printed in the
#' article. They play their intended roles (a small correction, and a bound
#' far below any plausible estimate) when \eqn{\xi^2} is not many orders of
#' magnitude from one, that is, when the residual standard deviation of
#' \code{y} is of the same order as the standard deviation of \code{x}. If
#' the two variables are recorded in very different units, rescale one of
#' them (and \code{width} with it) first. When the variance estimate is
#' positive but below \eqn{n^{-3}}, the bound is used in its place, as the
#' article prescribes, so the interval is wider and the required sample size
#' larger than the estimate itself would give; a note on the returned object
#' says that this happened and advises rescaling. On ordinary scales that is
#' a rare chance event of a small sample: with \eqn{\xi^2 = 1} it happened
#' in 0.05\% of samples of 10 pairs and 0.002\% of samples of 20 pairs
#' (200,000 simulated samples of normal data each), against 16\% and 1\% in
#' which the estimate was zero or negative.
#'
#' \strong{What the guarantee is.} The width guarantee is exact. The
#' coverage guarantee is asymptotic: coverage approaches \code{conf_level}
#' as \code{width} shrinks (their Theorem 7). At wide targets, and so small
#' final sample sizes, it falls short. In simulated studies (a 95\%
#' interval, the default \code{m_0}, one pair added at each stage, a slope
#' of 0.5, and \code{x} and the errors with unit variance, so that
#' \eqn{\xi^2 = 1}), coverage for normal \code{x} and errors was 0.94 at
#' \code{width = 0.2} (5,000 studies) and 0.83 at \code{width = 0.4}
#' (10,000 studies); for right skewed \code{x} and errors (a standardized
#' chi square distribution with 3 degrees of freedom) it was 0.86 (2,000
#' studies) and 0.71 (10,000 studies). The
#' reason is the variance estimate at small sample sizes. It was zero or
#' negative in 16.1\% (normal) and 26.5\% (skewed) of pilot samples of 10
#' pairs, the pilot at \code{width = 0.4}, and in 1.3\% and 9.9\% of pilot
#' samples of 20, the pilot at \code{width = 0.2} (200,000 simulated pilot
#' samples each); and when it is usable it tends to run low. When it is
#' positive but below \eqn{n^{-3}} the bound is used and the rule is then
#' met with an interval far narrower than requested: in other runs, 7
#' (normal) and 15 (skewed) of 6,000 studies ended that way at
#' \code{width = 0.4}, and none of 2,000 at \code{width = 0.2}. A larger
#' minimum pilot helps:
#' with \code{m_0 = 50} the coverages were 0.95 and 0.91 at
#' \code{width = 0.2} and 0.92 and 0.89 at \code{width = 0.4}. Use
#' \code{\link{seq_aipe_reg_coef_sensitivity}} to study a planned study
#' before collecting data, and prefer a narrow target or a larger
#' \code{m_0}.
#'
#' @inheritSection seq_aipe_smd Estimates That Are Not Yet Usable
#'
#' @param x A numeric vector: the predictor values observed so far.
#' @param y A numeric vector: the outcome values observed so far, paired by
#'   position with \code{x}. The two must have the same length, with at
#'   least 4 pairs; missing or infinite values are an error.
#' @param width The target full width of the confidence interval for
#'   \eqn{\beta_1} (\eqn{\omega} in the cited article): the upper limit
#'   minus the lower limit, in units of \code{y} per unit of \code{x}. A
#'   single positive number.
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
#'     \item{\code{reg_coef}}{the estimated slope (\code{NA} when \code{x}
#'       has no variation);}
#'     \item{\code{upper_limit}}{the upper confidence limit (\code{NA}
#'       whenever \code{lower_limit} is);}
#'     \item{\code{xi_squared}}{the estimated asymptotic variance
#'       \eqn{\hat{\xi}_n^2}, at least \eqn{n^{-3}} (the bound itself
#'       when the estimate is not yet usable);}
#'     \item{\code{criterion_n}}{the number of pairs the stopping rule
#'       currently requires (at least \code{current_n + 1} when the
#'       variance estimate is not yet usable or the estimate cannot be
#'       computed);}
#'     \item{\code{current_n}}{the current number of pairs.}
#'   }
#'   The object also carries these attributes:
#'   \describe{
#'     \item{\code{is_satisfied}}{\code{TRUE} when the stopping rule is
#'       met, otherwise \code{FALSE}; never \code{NA}, and always
#'       identical to \code{current_n >= criterion_n};}
#'     \item{\code{method}}{the string
#'       \code{"simple linear regression slope"};}
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
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#' regression: Obtaining regression coefficients that are accurate, not
#' simply significant. \emph{Psychological Methods}, \emph{8}(3), 305--321.
#' \doi{10.1037/1082-989X.8.3.305}
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#' for statistical power and accuracy in parameter estimation. \emph{Annual
#' Review of Psychology}, \emph{59}, 537--563.
#' \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' @examples
#' # How much does job satisfaction change per standard deviation of
#' # perceived autonomy (both standardized)? The aim is a 95% confidence
#' # interval for the slope that is no wider than 0.40.
#' set.seed(113)
#' m <- seq_aipe_reg_coef_pilot(width = 0.40, m_0 = 20)
#' autonomy <- rnorm(m)
#' satisfaction <- 0.5 * autonomy + rnorm(m, sd = sqrt(0.75))
#' out <- seq_aipe_reg_coef(autonomy, satisfaction, width = 0.40)
#' out
#'
#' # The rule is not met at the pilot. Add one employee at a time and check
#' # again, until it is.
#' while (!attr(out, "is_satisfied")) {
#'   new_autonomy <- rnorm(1)
#'   autonomy <- c(autonomy, new_autonomy)
#'   satisfaction <- c(satisfaction,
#'                     0.5 * new_autonomy + rnorm(1, sd = sqrt(0.75)))
#'   out <- seq_aipe_reg_coef(autonomy, satisfaction, width = 0.40)
#' }
#' out
#'
#' @seealso \code{\link{seq_aipe_reg_coef_pilot}} for the pilot size,
#'   \code{\link{seq_aipe_reg_coef_sensitivity}} to simulate the procedure
#'   before collecting data, and \code{\link{slr_slope}} for the point
#'   estimate.
#'
#'   \code{\link[DMAR]{ss_aipe_reg_coef}} in \pkg{DMAR} is the fixed sample
#'   size analog for when planning values are available, and
#'   \code{\link[DMAR]{ci_reg_coef}} gives a fixed sample size interval.
#' @family seq_aipe_reg_coef
#' @export
seq_aipe_reg_coef <- function(x, y, width, conf_level = 0.95) {
  pair <- .aipe_check_pair(x, y, "x", "y")
  width <- .check_width(width)
  conf_level <- .check_conf_level(conf_level)
  x <- pair$x
  y <- pair$y
  n <- length(x)

  # A predictor with no variation has no slope. That can happen by chance
  # in a small early sample (a discrete predictor, say), so it is reported
  # as not estimable rather than as an error. A response with no variation
  # has a slope of zero (and a variance estimate of zero: see
  # .var_reg_coef()); the slope formula would return the rounding noise of
  # values that are equal on paper.
  beta_1 <- if (.is_constant(x)) {
    NA_real_
  } else if (.is_constant(y)) {
    0
  } else {
    .slr_slope_value(x, y)
  }
  raw <- .var_reg_coef(x, y, floor = FALSE)

  # The slope has units, the bound n^-3 does not. When the variance
  # estimate is positive but below the bound, the bound is used, as for
  # every effect size, and the rule is evaluated; the note says that this
  # happened and that it is usually a matter of units, because the user
  # would otherwise have no way to see why the required sample size is so
  # large.
  bound_note <- .note_bound_applied(
    effect = "the slope is in units of y per unit of x",
    cause = paste0("this is what happens when x and y are recorded in ",
                   "units that make the variance of the slope very small ",
                   "(and, rarely, by chance in a small sample)"),
    remedy = paste0("Rescale x or y, and width with it, so that the ",
                    "residual standard deviation of y and the standard ",
                    "deviation of x are of the same order,"),
    topic = "seq_aipe_reg_coef"
  )

  .aipe_result(estimate = beta_1, raw = raw, n = n, width = width,
               conf_level = conf_level, term = "reg_coef",
               method = "simple linear regression slope",
               bound_note = bound_note)
}
