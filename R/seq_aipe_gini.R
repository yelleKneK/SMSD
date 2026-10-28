#' Pilot Sample Size for Sequential AIPE: Gini Index
#'
#' Returns the number of observations to collect before the sequential
#' accuracy in parameter estimation (AIPE) stopping rule for the population
#' Gini index can be evaluated for the first time with
#' \code{\link{seq_aipe_gini}}.
#'
#' The pilot size is
#' \deqn{m = \max\{m_0, \lceil 2 z_{\alpha/2} / \omega \rceil\},}
#' where \eqn{\omega} is \code{width}, \eqn{z_{\alpha/2}} is the upper
#' \eqn{\alpha/2} quantile of the standard normal distribution with
#' \eqn{\alpha} equal to one minus \code{conf_level}, and \eqn{m_0} is
#' \code{m_0}. This is the pilot size of Chattopadhyay and De (2016),
#' \eqn{\max\{4, \lceil z_{\alpha/2} / d \rceil\}}, written with the full
#' width \eqn{\omega = 2d} in place of their half width \eqn{d}. Whatever
#' the data, the stopping rule cannot be met with fewer than
#' \eqn{\lceil 2 z_{\alpha/2} / \omega \rceil} observations, which is why
#' the pilot is at least that large.
#'
#' @param width The target full width of the confidence interval for the
#'   population Gini index: the upper limit minus the lower limit (twice the
#'   half width \eqn{d} of Chattopadhyay and De, 2016). A single positive
#'   number less than 1, the length of the range of the Gini index.
#' @param conf_level The confidence level of the interval, a single number
#'   strictly between 0 and 1. Default 0.95.
#' @param m_0 The smallest pilot size allowed, a whole number of at least 4.
#'   Default 4.
#'
#' @return A single integer: the pilot sample size.
#'
#' @author
#' Ken Kelley \email{kkelley@@nd.edu},
#' Francis Bilson Darku \email{fbilsond@@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @references
#' Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
#' pre-specified error bound. \emph{Econometrics}, \emph{4}(3), 30.
#' \doi{10.3390/econometrics4030030}
#'
#' @examples
#' # A 95% interval for a Gini index no wider than 0.10
#' seq_aipe_gini_pilot(width = 0.10)
#'
#' # A 90% interval no wider than 0.04, a half width of 0.02
#' seq_aipe_gini_pilot(width = 0.04, conf_level = 0.90)
#'
#' @seealso \code{\link{seq_aipe_gini}}
#' @family seq_aipe_gini
#' @export
seq_aipe_gini_pilot <- function(width, conf_level = 0.95, m_0 = 4) {
  .aipe_pilot(width, conf_level, m_0, range = 1)
}


#' Sequential AIPE Check for the Gini Index
#'
#' Checks, with the data collected so far, whether enough observations are
#' in hand for a confidence interval for the population Gini index of
#' inequality (of incomes, wealth, sales, and the like) that is no wider
#' than a chosen width, without assuming a distribution and without a
#' planning value. This is the sequential fixed width interval procedure of
#' Chattopadhyay and De (2016), expressed in the package's accuracy in
#' parameter estimation (AIPE) form with the full width of the interval as
#' the target.
#'
#' The function evaluates the stopping rule once and returns. The user
#' collects the pilot sample (see \code{\link{seq_aipe_gini_pilot}}), calls
#' \code{seq_aipe_gini()}, and, while the rule is not met, adds one or more
#' observations and calls it again. The function prints nothing and keeps no
#' state.
#'
#' \strong{Estimate.} The Gini index \eqn{G = \Delta / (2\mu)}, with
#' \eqn{\Delta = E|X_1 - X_2|} Gini's mean difference, is estimated by
#' \eqn{G_n = \hat{\Delta}_n / (2 \bar{X}_n)}, as returned by
#' \code{\link{gini_index}}. It is defined for nonnegative data with a
#' positive mean, so a negative value in \code{x} is an error. A sample in
#' which every value is zero cannot estimate the index; because that can
#' happen by chance in a small early sample, it is reported as "not yet
#' estimable" rather than as an error (see the section below).
#'
#' \strong{Stopping rule.} Sampling stops at the first \eqn{n}, at or after
#' the pilot size, with
#' \deqn{n \ge \left\lceil \frac{4 z_{\alpha/2}^2}{\omega^2}
#'   \left(\hat{\xi}_n^2 + \frac{1}{n}\right) \right\rceil,}
#' which is the stopping rule of Chattopadhyay and De (2016, Equation 12),
#' \eqn{n \ge (z_{\alpha/2}/d)^2 (V_n^2 + 1/n)}, with \eqn{\omega = 2d}.
#' Here \eqn{\hat{\xi}_n^2} is their estimator \eqn{V_n^2} (Equation 11) of
#' the asymptotic variance of \eqn{\sqrt{n}(G_n - G)} (Equation 6), built
#' from the sample mean and variance, the Gini mean difference, a second
#' U-statistic, and a jackknife variance. The reported interval is
#' \eqn{G_n \pm z_{\alpha/2} \hat{\xi}_n / \sqrt{n}}. Whenever the rule is
#' met, the width of that interval is at most \code{width}. The article
#' bounds \eqn{V_n^2} below by zero; the package bounds a positive estimate
#' below by \eqn{n^{-3}}, as for the other effect sizes, and treats an
#' estimate that is zero or negative as not yet usable (see the section
#' below).
#'
#' \strong{Which interval.} The interval reported is built as for every
#' effect size in the family, from the estimated standard error. It is not
#' the interval of Chattopadhyay and De (2016), which is the fixed width
#' interval \eqn{(G_n - d, G_n + d)}, always exactly \eqn{\omega} wide,
#' whose coverage their stopping rule is designed to secure. The interval
#' reported is never wider than that one once the rule is met, and it
#' covers less often: in simulated studies with lognormal data, a Gini
#' index of 0.3, and a 95\% level (10,000 studies each, the default
#' \code{m_0}, one observation added at each stage), the interval reported
#' covered the Gini index in 0.890 of the studies at \code{width = 0.2} and
#' 0.916 at \code{width = 0.1}, and the article's interval in 0.986 and
#' 0.960. For the article's interval, take the estimate (the row
#' \code{gini}) plus and minus half of \code{width}.
#'
#' Kelley, Darku, and Chattopadhyay (2018) cite this procedure and give the
#' general sequential AIPE framework that the rest of the \code{seq_aipe_*}
#' family implements; their article does not itself contain a Gini index
#' procedure.
#'
#' \strong{What the guarantee is.} The width guarantee is exact. The
#' coverage guarantee is asymptotic: coverage approaches \code{conf_level}
#' as \code{width} shrinks. In simulated studies with a 95\% interval and
#' one observation added at each stage, coverage for exponential data (a
#' Gini index of 0.5) was 0.951 at \code{width = 0.1} (3,000 studies, about
#' 142 observations at the end), 0.937 at \code{width = 0.2} (2,000 studies,
#' about 44 observations), and 0.945 at \code{width = 0.4} (6,000 studies,
#' about 17 observations). For lognormal data with a Gini index of 0.3 it
#' was 0.890 at \code{width = 0.2} (10,000 studies, about 30 observations)
#' and 0.870 at \code{width = 0.4} (about 13 observations). The variance
#' estimate was never zero, negative, or below its lower bound in these
#' studies; the shortfall comes from the small final samples that wide
#' targets allow.
#' Use \code{\link{seq_aipe_gini_sensitivity}} to study a planned study
#' before collecting data.
#'
#' @inheritSection seq_aipe_smd Estimates That Are Not Yet Usable
#'
#' @param x A numeric vector of nonnegative values: the observations
#'   collected so far, at least 4 of them. Missing, infinite, or negative
#'   values are an error.
#' @param width The target full width of the confidence interval for the
#'   Gini index: the upper limit minus the lower limit (twice the half
#'   width \eqn{d} of Chattopadhyay and De, 2016). A single positive number
#'   less than 1: the Gini index lies between 0 and 1, so an interval 1 or
#'   more wide would say nothing, and such a target is an error.
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
#'     \item{\code{gini}}{the estimated Gini index (\code{NA} when every
#'       value is zero);}
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
#'     \item{\code{method}}{the string \code{"Gini index"};}
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
#' Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
#' pre-specified error bound. \emph{Econometrics}, \emph{4}(3), 30.
#' \doi{10.3390/econometrics4030030}
#'
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#' parameter estimation for a general class of effect sizes: A sequential
#' approach. \emph{Psychological Methods}, \emph{23}(2), 226--243.
#' \doi{10.1037/met0000127}
#'
#' @examples
#' # How unequal are household incomes (in thousands of dollars) in a
#' # region? The aim is a 95% confidence interval for the Gini index that
#' # is no wider than 0.15.
#' set.seed(113)
#' m <- seq_aipe_gini_pilot(width = 0.15)
#' income <- rlnorm(m, meanlog = 4, sdlog = 0.6)
#' out <- seq_aipe_gini(income, width = 0.15)
#' out
#'
#' # The rule is not met at the pilot. Add one household at a time and check
#' # again, until it is.
#' while (!attr(out, "is_satisfied")) {
#'   income <- c(income, rlnorm(1, meanlog = 4, sdlog = 0.6))
#'   out <- seq_aipe_gini(income, width = 0.15)
#' }
#' out
#'
#' @seealso \code{\link{seq_aipe_gini_pilot}} for the pilot size,
#'   \code{\link{seq_aipe_gini_sensitivity}} to simulate the procedure
#'   before collecting data, and \code{\link{gini_index}} and
#'   \code{\link{gini_mean_diff}} for the point estimates. (\pkg{DMAR} has
#'   no fixed sample size analog for the Gini index.)
#' @family seq_aipe_gini
#' @export
seq_aipe_gini <- function(x, width, conf_level = 0.95) {
  x <- .aipe_check_data(x, "x", nonnegative = TRUE)
  width <- .aipe_check_width_bounded(width, range = 1)
  conf_level <- .check_conf_level(conf_level)
  n <- length(x)

  if (max(x) > 0) {
    # A sample with no variation has a Gini index of zero (and a variance
    # estimate of zero: see .var_gini()); gini_index() would return the
    # rounding noise of values that are equal on paper.
    gi <- if (.is_constant(x)) 0 else gini_index(x)
    raw <- .var_gini(x, floor = FALSE)
  } else {
    # Every value is zero: the mean is zero and the index is undefined.
    gi <- NA_real_
    raw <- NaN
  }

  .aipe_result(estimate = gi, raw = raw, n = n, width = width,
               conf_level = conf_level, term = "gini", method = "Gini index")
}
