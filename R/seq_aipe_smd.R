#' Pilot Sample Size for Sequential AIPE: Standardized Mean Difference
#'
#' Returns the number of observations to collect from each of the two groups
#' before the sequential accuracy in parameter estimation (AIPE) stopping
#' rule for the standardized mean difference can be evaluated for the first
#' time with \code{\link{seq_aipe_smd}}.
#'
#' The pilot size is
#' \deqn{m = \max\{m_0, \lceil 2 z_{\alpha/2} / \omega \rceil\},}
#' where \eqn{\omega} is \code{width}, \eqn{z_{\alpha/2}} is the upper
#' \eqn{\alpha/2} quantile of the standard normal distribution with
#' \eqn{\alpha} equal to one minus \code{conf_level}, and \eqn{m_0} is
#' \code{m_0} (Kelley, Darku, and Chattopadhyay, 2018, Equation 21).
#' Whatever the data, the stopping rule cannot be met with fewer than
#' \eqn{\lceil 2 z_{\alpha/2} / \omega \rceil} observations per group, which
#' is why the pilot is at least that large.
#'
#' @param width The target full width of the confidence interval for the
#'   population standardized mean difference (\eqn{\omega} in the cited
#'   article): the upper limit minus the lower limit. A single positive
#'   number.
#' @param conf_level The confidence level of the interval, a single number
#'   strictly between 0 and 1. Default 0.95.
#' @param m_0 The smallest pilot size allowed, a whole number of at least 4
#'   (the variance estimator needs four observations per group). Default 4.
#'   Raising it makes an unusable early variance estimate less likely; see
#'   \code{\link{seq_aipe_smd}}.
#'
#' @return A single integer: the pilot sample size per group.
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
#' # A 95% interval for a standardized mean difference no wider than 0.30
#' seq_aipe_smd_pilot(width = 0.30)
#'
#' # The same target at 99% confidence, with at least 20 per group
#' seq_aipe_smd_pilot(width = 0.30, conf_level = 0.99, m_0 = 20)
#'
#' @seealso \code{\link{seq_aipe_smd}}
#' @family seq_aipe_smd
#' @export
seq_aipe_smd_pilot <- function(width, conf_level = 0.95, m_0 = 4) {
  .aipe_pilot(width, conf_level, m_0)
}


#' Sequential AIPE Check for the Standardized Mean Difference
#'
#' Checks, with the data collected so far from two independent groups,
#' whether enough observations are in hand for a confidence interval for the
#' population standardized mean difference
#' \eqn{\delta = (\mu_1 - \mu_2)/\sigma} that is no wider than a chosen
#' width, without assuming normality and without a planning value for
#' \eqn{\delta}. This is the sequential accuracy in parameter estimation
#' (AIPE) procedure of Kelley, Darku, and Chattopadhyay (2018).
#'
#' The function evaluates the stopping rule once and returns. The user
#' collects the pilot sample (see \code{\link{seq_aipe_smd_pilot}}), calls
#' \code{seq_aipe_smd()}, and, while the rule is not met, adds one or more
#' observations to \emph{each} group and calls it again. The function prints
#' nothing and keeps no state.
#'
#' \strong{Estimate.} The standardized mean difference is
#' \eqn{d_n = (\bar{X}_{1n} - \bar{X}_{2n}) / s_{pn}}, with
#' \eqn{s_{pn}^2 = (s_{1n}^2 + s_{2n}^2)/2} the pooled variance of the two
#' groups of common size \eqn{n} (their Equation 24), computed by
#' \code{\link[DMAR]{smd}}.
#'
#' \strong{Stopping rule.} Sampling stops at the first \eqn{n} (per group),
#' at or after the pilot size, with
#' \deqn{n \ge \left\lceil \frac{4 z_{\alpha/2}^2}{\omega^2}
#'   \left(\hat{\xi}_n^2 + \frac{1}{n}\right) \right\rceil}
#' (their Equation 22), where \eqn{\hat{\xi}_n^2} estimates the asymptotic
#' variance of \eqn{\sqrt{n}(d_n - \delta)} from the sample means, the pooled
#' variance, and unbiased estimates of the third and fourth central moments
#' of each group (their Equation 27). The reported interval is
#' \eqn{d_n \pm z_{\alpha/2} \hat{\xi}_n / \sqrt{n}} (their Equation 28).
#' Whenever the rule is met, the width of that interval is at most
#' \code{width}.
#'
#' \strong{Correction to Equation 26.} Equation 26 of the article, and the
#' sample version printed under its Equation 27, omit a factor of
#' \eqn{1/2} in the term that carries the third central moments. The package
#' implements the corrected variance,
#' \deqn{\xi^2 = 2 - \frac{(\mu_1 - \mu_2)(\mu_{13} - \mu_{23})}{2 \sigma^4}
#'   + \frac{(\mu_1 - \mu_2)^2}{4 \sigma^6}
#'   \left(\frac{\mu_{14} + \mu_{24}}{4} - \frac{\sigma^4}{2}\right),}
#' where \eqn{\mu_{kj}} is the \eqn{j}th central moment of group \eqn{k}.
#' The two forms differ only when the means differ and the third central
#' moments differ; they agree when the groups have the same third central
#' moment (in particular for symmetric distributions) and when the two means
#' are equal. The authors' correction is posted at
#' \url{https://kenkelley.org/wp-content/uploads/articles/Kelley_Darku_Chattopadhyay_Equation_26.pdf}
#' and is installed with the package, at the path returned by
#' \preformatted{system.file("eq26_correction.pdf", package = "SMSD")}
#'
#' \strong{What the guarantee is.} The width guarantee is exact. The
#' coverage guarantee is asymptotic: the coverage of the interval approaches
#' \code{conf_level} as \code{width} shrinks, and so as the final sample size
#' grows (their Theorem 4). In simulated studies with a 95\% interval, a
#' true \eqn{\delta} of 0.5, and one observation per group added at each
#' stage, coverage was 0.952 at \code{width = 0.4} (6,000 studies, about 200
#' per group at the end) and 0.946 at \code{width = 0.2} (2,000 studies,
#' about 790 per group) for normal data, and 0.944 and 0.949 for right skewed
#' data (a standardized chi square distribution with 3 degrees of freedom);
#' the Monte Carlo standard errors are 0.003 and 0.005. The variance estimate
#' was never zero, negative, or below its lower bound in these studies. At a
#' wide target the final sample is small, and coverage falls short for
#' skewed data and, to a lesser degree, for a large \eqn{\delta}: at
#' \code{width = 0.8} (10,000 studies each, about 50 to 105 per group at the
#' end) coverage was 0.950, 0.946, and 0.943 for normal data and 0.926,
#' 0.926, and 0.919 for the skewed data at \eqn{\delta} of 0.5, 1, and 2.
#' Use \code{\link{seq_aipe_smd_sensitivity}} to study other conditions
#' before collecting data.
#'
#' \strong{A group with no variation.} When every observation in one group
#' is the same (values that agree to ten significant digits count as the
#' same) while the other group varies, the estimate is still reported, but
#' the variance estimate is not yet usable: the variance of that group has
#' not been estimated, and the pooled standard deviation rests on the other
#' group alone, so that it is too small and the estimate too large in size.
#' The check reports the rule as not met and a note names the group (see
#' the section below), as \code{\link{seq_aipe_c}} and
#' \code{\link{seq_budget}} do with a constant group. In simulated studies
#' of a 5 point item with 70\% of one group's answers at the top of the
#' scale (\eqn{\delta = 1.08}, \code{width = 1.5}, 4,000 studies), 2.4\%
#' of the studies would have stopped with that group constant, and 41\% of
#' their intervals covered \eqn{\delta}; with the group treated as not yet
#' estimated, coverage over all studies was 0.927 against 0.917. The price
#' is paid where a group stays constant: for a binary outcome with
#' proportions 0.30 and 0.02 (\code{width = 0.8}, 2,000 studies) the mean
#' final sample size was 64 per group against 33, and coverage 0.967 against
#' 0.943. If a group shows no variation well past the pilot sample,
#' reconsider the measure.
#'
#' @section Estimates That Are Not Yet Usable:
#' When the raw variance estimate is zero (exactly, or to within rounding
#' error), negative, or not finite, which is common at small sample sizes
#' for some effect sizes, or when the point estimate cannot be computed (a
#' sample with no variation, where values that agree to ten significant
#' digits count as equal), the check reports the stopping rule as not met,
#' sets \code{attr(out, "variance_unusable")} to \code{TRUE}, and attaches
#' a note that \code{print()} shows; see \code{\link{smsd_check_methods}}.
#' This is a deliberate departure from the stopping rule as printed in the
#' source articles, under which the lower bound \eqn{n^{-3}} would stand in
#' for such an estimate and the rule would be met at once with an interval
#' far narrower than requested; a positive estimate, however small, is used
#' as the articles prescribe (the larger of it and \eqn{n^{-3}}).
#'
#' @param group_1,group_2 Numeric vectors: the observations collected so far
#'   from the first and the second group. Each needs at least 4
#'   observations, the two must have the same length (the procedure adds the
#'   same number of observations to each group at every stage), and missing
#'   or infinite values are an error.
#' @param width The target full width of the confidence interval for
#'   \eqn{\delta} (\eqn{\omega} in the cited article): the upper limit minus
#'   the lower limit. A single positive number.
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
#'     \item{\code{smd}}{the estimated standardized mean difference
#'       (\code{group_1} minus \code{group_2}, in pooled standard deviation
#'       units; \code{NA} when neither group has any variation);}
#'     \item{\code{upper_limit}}{the upper confidence limit (\code{NA}
#'       whenever \code{lower_limit} is);}
#'     \item{\code{xi_squared}}{the estimated asymptotic variance
#'       \eqn{\hat{\xi}_n^2}, at least \eqn{n^{-3}} (the bound itself
#'       when the estimate is not yet usable);}
#'     \item{\code{criterion_n}}{the per group sample size the stopping rule
#'       currently requires (at least \code{current_n + 1} when the
#'       variance estimate is not yet usable or the estimate cannot be
#'       computed);}
#'     \item{\code{current_n}}{the current per group sample size.}
#'   }
#'   The object also carries these attributes:
#'   \describe{
#'     \item{\code{is_satisfied}}{\code{TRUE} when the stopping rule is
#'       met, otherwise \code{FALSE}; never \code{NA}, and always
#'       identical to \code{current_n >= criterion_n};}
#'     \item{\code{method}}{the string \code{"standardized mean difference"};}
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
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#' standardized mean difference: Accuracy in parameter estimation via narrow
#' confidence intervals. \emph{Psychological Methods}, \emph{11}(4),
#' 363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#' for statistical power and accuracy in parameter estimation. \emph{Annual
#' Review of Psychology}, \emph{59}, 537--563.
#' \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' @examples
#' # Does a training program raise employee engagement (rated 1 to 7)? The
#' # aim is a 95% confidence interval for the standardized mean difference
#' # that is no wider than 0.60.
#' set.seed(113)
#' m <- seq_aipe_smd_pilot(width = 0.60)
#' trained <- rnorm(m, mean = 4.9, sd = 1)
#' control <- rnorm(m, mean = 4.5, sd = 1)
#' out <- seq_aipe_smd(trained, control, width = 0.60)
#' out
#'
#' # The rule is not met at the pilot. Add one employee to each group and
#' # check again, until it is.
#' while (!attr(out, "is_satisfied")) {
#'   trained <- c(trained, rnorm(1, mean = 4.9, sd = 1))
#'   control <- c(control, rnorm(1, mean = 4.5, sd = 1))
#'   out <- seq_aipe_smd(trained, control, width = 0.60)
#' }
#' out
#'
#' @seealso \code{\link{seq_aipe_smd_pilot}} for the pilot size and
#'   \code{\link{seq_aipe_smd_sensitivity}} to simulate the procedure before
#'   collecting data. In \pkg{DMAR}: \code{\link[DMAR]{ss_aipe_smd}}, the
#'   fixed sample size analog for when a planning value of \eqn{\delta} is
#'   available, and \code{\link[DMAR]{smd}}, the point estimator used here.
#' @family seq_aipe_smd
#' @export
seq_aipe_smd <- function(group_1, group_2, width, conf_level = 0.95) {
  group_1 <- .aipe_check_data(group_1, "group_1")
  group_2 <- .aipe_check_data(group_2, "group_2")
  if (length(group_1) != length(group_2)) {
    .smsd_stop("'group_1' and 'group_2' must have the same length (the ",
               "procedure adds the same number of observations to each ",
               "group at every stage); they have ", length(group_1),
               " and ", length(group_2), " observations.")
  }
  width <- .check_width(width)
  conf_level <- .check_conf_level(conf_level)
  n <- length(group_1)

  # Point estimate reused from DMAR. With the equal group sizes this
  # function requires, DMAR's pooled standard deviation is
  # sqrt((s_1^2 + s_2^2) / 2), the s_pn of the article, so the dependency
  # keeps the two packages in step on the definition of the standardized
  # mean difference. When neither group has any variation the pooled
  # standard deviation is zero and there is no standardized mean
  # difference, which .aipe_result() reports as not estimable. The
  # standardized mean difference does not depend on the common scale of the
  # two groups, which are divided by a power of two of their magnitude so
  # that no sum of squares can overflow or underflow; see .aipe_scale().
  scale <- .aipe_scale(group_1, group_2)
  group_1 <- group_1 / scale
  group_2 <- group_2 / scale
  if (.is_constant(group_1) && .is_constant(group_2)) {
    smd <- NA_real_
  } else {
    smd_out <- DMAR::smd(group_1, group_2)
    smd <- smd_out$value[smd_out$term == "smd"]
  }
  raw <- .var_smd(group_1, group_2, floor = FALSE)
  # One group with no variation while the other varies: the estimate is
  # still reported, the variance estimate is not usable (.var_smd() returns
  # NaN), and the note says which group it is. With none in either group
  # the estimate itself is missing and the standard notes say so.
  flat <- c(.is_constant(group_1), .is_constant(group_2))
  flat_note <- if (xor(flat[1L], flat[2L])) {
    .note_group_constant(if (flat[1L]) "'group_1'" else "'group_2'")
  }

  .aipe_result(estimate = smd, raw = raw, n = n, width = width,
               conf_level = conf_level, term = "smd",
               method = "standardized mean difference",
               unusable_note = flat_note)
}
