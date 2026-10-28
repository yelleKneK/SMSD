#' Pilot Sample Size for Sequential AIPE: Contrast
#'
#' Returns the number of observations to collect from each group before the
#' sequential accuracy in parameter estimation (AIPE) stopping rule for a
#' contrast of group means, \eqn{\psi = \sum_{k=1}^{K} c_k \mu_k}, can be
#' evaluated for the first time with \code{\link{seq_aipe_c}}.
#'
#' The pilot size is
#' \deqn{m = \max\{m_0, \lceil 2 z_{\alpha/2} / \omega \rceil\},}
#' where \eqn{\omega} is \code{width}, \eqn{z_{\alpha/2}} is the upper
#' \eqn{\alpha/2} quantile of the standard normal distribution with
#' \eqn{\alpha} equal to one minus \code{conf_level}, and \eqn{m_0} is
#' \code{m_0} (Kelley, Darku, and Chattopadhyay, 2018, Equation 21).
#' Whatever the data, the stopping rule cannot be met with fewer than
#' \eqn{\lceil 2 z_{\alpha/2} / \omega \rceil} observations per group, which
#' is why the pilot is at least that large. Because a contrast is not scale
#' free, the pilot size depends on the units of \code{width}: see
#' \code{\link{seq_aipe_c}}.
#'
#' @param width The target full width of the confidence interval for the
#'   contrast (\eqn{\omega} in the cited article): the upper limit minus the
#'   lower limit, in the units of the outcome. A single positive number.
#' @param conf_level The confidence level of the interval, a single number
#'   strictly between 0 and 1. Default 0.95.
#' @param m_0 The smallest pilot size allowed per group, a whole number of
#'   at least 4. Default 4.
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
#' # A 95% interval for a contrast no wider than 0.40 outcome units
#' seq_aipe_c_pilot(width = 0.40)
#'
#' # The same target with at least 20 per group in the pilot
#' seq_aipe_c_pilot(width = 0.40, m_0 = 20)
#'
#' @seealso \code{\link{seq_aipe_c}}
#' @family seq_aipe_c
#' @export
seq_aipe_c_pilot <- function(width, conf_level = 0.95, m_0 = 4) {
  .aipe_pilot(width, conf_level, m_0)
}


#' Sequential AIPE Check for a Contrast
#'
#' Checks, with the data collected so far from K independent groups, whether
#' enough observations are in hand for a confidence interval for a contrast
#' of the group means, \eqn{\psi = \sum_{k=1}^{K} c_k \mu_k} (for example
#' one condition against the average of two others), that is no wider than
#' a chosen width, without assuming normality or equal variances and without
#' planning values. The contrast is one of the effect sizes covered by the
#' sequential accuracy in parameter estimation (AIPE) procedure of Kelley,
#' Darku, and Chattopadhyay (2018, Example 6).
#'
#' The function evaluates the stopping rule once and returns. The user
#' collects the pilot sample from every group (see
#' \code{\link{seq_aipe_c_pilot}}), calls \code{seq_aipe_c()}, and, while
#' the rule is not met, adds the same number of observations to
#' \emph{every} group and calls it again. The function prints nothing and
#' keeps no state.
#'
#' \strong{Estimate.} The contrast is estimated by
#' \eqn{\hat{\psi}_n = \sum_k c_k \bar{X}_{kn}}, the weighted sum of the
#' group means (their Equation 16). The article writes the weights as
#' summing to one; the procedure does not depend on that, and any weights
#' that are finite and not all zero are accepted (weights that sum to zero
#' give a comparison among groups).
#'
#' \strong{Stopping rule.} With \eqn{n} observations in each group, sampling
#' stops at the first \eqn{n}, at or after the pilot size, with
#' \deqn{n \ge \left\lceil \frac{4 z_{\alpha/2}^2}{\omega^2}
#'   \left(\hat{\xi}_n^2 + \frac{1}{n}\right) \right\rceil}
#' (their Equation 22), where
#' \eqn{\hat{\xi}_n^2 = \sum_k c_k^2 s_{kn}^2} estimates the asymptotic
#' variance of \eqn{\sqrt{n}(\hat{\psi}_n - \psi)}. The reported interval is
#' \eqn{\hat{\psi}_n \pm z_{\alpha/2} \hat{\xi}_n / \sqrt{n}} (their
#' Equation 23). Whenever the rule is met, the width of that interval is at
#' most \code{width}.
#'
#' \strong{Balanced data.} The function implements the balanced procedure
#' of the article's Example 6, which adds an observation to every group at
#' each stage, so every group must have the same number of observations.
#' Groups of unequal size are an error rather than a fallback: a rule
#' applied group by group cannot guarantee the requested width. Bring every
#' group up to the common size and call the function again. The article
#' also outlines an extension to unbalanced designs (its Equations 39 to
#' 47), in which the optimal group sizes are proportional to
#' \eqn{c_k \sigma_k} and are to be estimated sequentially; the article
#' prints no stopping rule for it, and that extension is not implemented.
#'
#' \strong{Units.} Unlike a standardized effect size, the contrast and
#' \code{width} are in the units of the outcome, while the correction term
#' \eqn{1/n} in the stopping rule and the lower bound \eqn{n^{-3}} on the
#' variance estimate are pure numbers, as printed in the article. They play
#' their intended roles (a small correction, and a bound far below any
#' plausible estimate) when \eqn{\xi^2} is not many orders of magnitude from
#' one. If the outcome is recorded in units that make its variance very
#' small or very large, rescale it (and \code{width} with it) first. When
#' the variance estimate is positive but below \eqn{n^{-3}}, which on
#' ordinary scales a weighted sum of sample variances essentially never is,
#' the bound is used in its place, as the article prescribes, so the
#' interval is wider and the required sample size larger than the estimate
#' itself would give; a note on the returned object says that this happened
#' and advises rescaling.
#'
#' \strong{A group with no variation.} When every observation in a group
#' with a nonzero weight is the same (values that agree to ten significant
#' digits count as the same), the variance of that group has not been
#' estimated, and so neither has the variance of the contrast: the contrast
#' is still estimated, but the check reports the rule as not met and a note
#' names the group (see the section below), as \code{\link{seq_budget}}
#' does with the same data. Summing the variances of the other groups
#' instead would leave the constant group out of the variance of the
#' contrast. In simulated studies of three groups answering a 5 point item,
#' one of them with 70\% of its answers at the top of the scale (weights 1,
#' -1/2, -1/2, \code{width = 2}, 4,000 studies), 23\% of the studies would
#' have stopped with that group constant, after 4 to 7 observations per
#' group, and 74\% of their intervals covered the contrast; with the group
#' treated as not yet estimated, coverage over all studies was 0.934
#' against 0.899. The price is paid where a group stays constant: for binary
#' outcomes with proportions 0.30, 0.02, and 0.02 (\code{width = 0.5}, 2,000
#' studies) the mean final sample size was 75 per group against 18, and
#' coverage 0.939 against 0.891. If a group shows no variation well past the
#' pilot sample, reconsider the measure. A group whose weight is zero takes
#' no part in the contrast and is not looked at.
#'
#' \strong{Names.} \code{groups} and \code{contrast_weights} are matched by
#' position. If both are named, the names must be the same and in the same
#' order; otherwise the function stops, so that weights named in another
#' order (for example against a list made by \code{split()}, which is in
#' alphabetical order) are not applied to the wrong groups.
#'
#' \strong{What the guarantee is.} The width guarantee is exact. The
#' coverage guarantee is asymptotic: coverage approaches \code{conf_level}
#' as \code{width} shrinks (their Theorem 2). In simulated studies with a
#' 95\% interval and one observation per group added at each stage,
#' coverage for three normal groups with unit variances and weights 1,
#' -1/2, -1/2 was 0.948 at \code{width = 0.4} (6,000 studies, about 144 per
#' group at the end) and 0.945 at \code{width = 0.2} (2,000 studies, about
#' 577 per group); for three groups with skewed distributions and unequal
#' variances it was 0.949 and 0.951 (Monte Carlo standard errors 0.003 and
#' 0.005). At a wide target the final sample is small and coverage falls a
#' little short: with means 0.5, 0, and 0 it was 0.936 at \code{width = 1}
#' (10,000 studies, about 23 per group). The variance estimate, a weighted
#' sum of sample variances, cannot be negative, and it was never below its
#' lower bound in these studies.
#'
#' @inheritSection seq_aipe_smd Estimates That Are Not Yet Usable
#'
#' @param groups A list of K numeric vectors (K at least 2), one per group:
#'   the observations collected so far from each population. Every vector
#'   needs at least 4 observations and all must have the same length;
#'   missing or infinite values are an error. A matrix or data frame is not
#'   accepted (a data frame in "long" format would be misread); with one
#'   column per group use \code{as.list(as.data.frame(m))}, and with an
#'   outcome \code{y} and a grouping variable \code{g} use
#'   \code{split(y, g)}.
#' @param contrast_weights A numeric vector of K contrast weights
#'   \eqn{c_k}, in the order of \code{groups} (see "Names" in Details).
#'   They must be finite and not all zero.
#' @param width The target full width of the confidence interval for the
#'   contrast (\eqn{\omega} in the cited article): the upper limit minus the
#'   lower limit, in the units of the outcome. A single positive number.
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
#'     \item{\code{contrast}}{the estimated contrast;}
#'     \item{\code{upper_limit}}{the upper confidence limit (\code{NA}
#'       whenever \code{lower_limit} is);}
#'     \item{\code{xi_squared}}{the estimated asymptotic variance
#'       \eqn{\hat{\xi}_n^2}, at least \eqn{n^{-3}} (the bound itself
#'       when the estimate is not yet usable);}
#'     \item{\code{criterion_n}}{the per group sample size the stopping rule
#'       currently requires (at least \code{current_n + 1} when the
#'       variance estimate is not yet usable or the estimate cannot be
#'       computed);}
#'     \item{\code{current_n}}{the current per group sample size (the total
#'       sample size is K times this).}
#'   }
#'   The object also carries these attributes:
#'   \describe{
#'     \item{\code{is_satisfied}}{\code{TRUE} when the stopping rule is
#'       met, otherwise \code{FALSE}; never \code{NA}, and always
#'       identical to \code{current_n >= criterion_n};}
#'     \item{\code{method}}{the string \code{"contrast"};}
#'     \item{\code{n_per_group}}{an integer vector with the size of each
#'       group, named when \code{groups} is;}
#'     \item{\code{contrast_weights}}{the weights, as doubles;}
#'     \item{\code{current_n}, \code{criterion_n}}{the values in the
#'       rows of those names (both per group);}
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
#' # Does a new onboarding program change first year performance ratings
#' # (standardized) relative to the average of two existing programs? The
#' # aim is a 95% confidence interval for that contrast no wider than 0.80.
#' set.seed(113)
#' program_means <- c(new = 0.5, existing_a = 0, existing_b = 0)
#' weights <- c(1, -0.5, -0.5)
#' m <- seq_aipe_c_pilot(width = 0.80)
#' ratings <- lapply(program_means, function(mu) rnorm(m, mean = mu))
#' out <- seq_aipe_c(ratings, weights, width = 0.80)
#' out
#'
#' # The rule is not met at the pilot. Add one employee to every program and
#' # check again, until it is.
#' while (!attr(out, "is_satisfied")) {
#'   ratings <- Map(function(r, mu) c(r, rnorm(1, mean = mu)),
#'                  ratings, program_means)
#'   out <- seq_aipe_c(ratings, weights, width = 0.80)
#' }
#' out
#' attr(out, "n_per_group")
#'
#' @seealso \code{\link{seq_aipe_c_pilot}} for the pilot size, and
#'   \code{\link{seq_budget}} for a sequential interval for a contrast when
#'   the groups have different sampling costs and there is a budget. In
#'   \pkg{DMAR}: \code{\link[DMAR]{ss_aipe_c}}, the fixed sample size analog
#'   for when planning values are available, \code{\link[DMAR]{ss_aipe_sc}}
#'   for the standardized contrast, and \code{\link[DMAR]{ci_c}}.
#' @family seq_aipe_c
#' @export
seq_aipe_c <- function(groups, contrast_weights, width, conf_level = 0.95) {
  if (missing(groups)) {
    .smsd_stop("'groups' must be supplied: a list with one numeric vector ",
               "per group.")
  }
  if (!is.list(groups) || is.data.frame(groups)) {
    .smsd_stop("'groups' must be a list with one numeric vector per group; ",
               "got ", .smsd_describe(groups), ". With one column per ",
               "group use as.list(as.data.frame(m)); with an outcome y and ",
               "a grouping variable g use split(y, g).")
  }
  if (length(groups) < 2L) {
    .smsd_stop("'groups' must contain at least 2 groups; it has ",
               length(groups), ".")
  }
  group_names <- names(groups)
  groups <- lapply(seq_along(groups), function(k) {
    .check_data(groups[[k]], sprintf("groups[[%d]]", k), min_n = 4L,
                na.rm = NULL)
  })
  names(groups) <- group_names
  # The weights are matched to the groups by position. The names are read
  # before the check, which returns the weights without them.
  weights <- .aipe_check_weights(contrast_weights, length(groups))
  .budget_check_names(list(groups = group_names,
                           contrast_weights = names(contrast_weights)))
  contrast_weights <- weights
  width <- .check_width(width)
  conf_level <- .check_conf_level(conf_level)

  n_per_group <- lengths(groups)
  if (!all(n_per_group == n_per_group[1L])) {
    .smsd_stop("All groups must have the same per group sample size (the ",
               "sequential contrast procedure adds an observation to ",
               "every group at each stage); the groups currently have ",
               "sizes ", paste(n_per_group, collapse = ", "),
               ". Bring every group up to the common size and call again.")
  }
  # The common per group size, without the name of the first group (which
  # would otherwise travel into current_n, criterion_n, and the row names
  # of summary()).
  n <- unname(n_per_group[1L])

  contrast <- sum(contrast_weights * vapply(groups, mean, numeric(1)))
  # Balanced rule of Kelley, Darku, and Chattopadhyay (2018, Equation 22)
  # for the contrast of their Example 6: the criterion bounds the variance
  # of the whole contrast, so n >= criterion guarantees a full width of at
  # most `width`.
  raw <- .var_c(groups, contrast_weights, floor = FALSE)
  # A group with a nonzero weight and no variation: its variance has not
  # been estimated, so neither has the variance of the contrast (.var_c()
  # returns NaN), and the note says which group it is.
  flat <- which(contrast_weights != 0 &
                  vapply(groups, .is_constant, logical(1)))
  flat_note <- if (length(flat)) {
    .note_group_constant(.budget_group_label(flat), length(flat) > 1L)
  }

  # The contrast has units, the bound n^-3 does not. When the variance
  # estimate is positive but below the bound (which a weighted sum of
  # sample variances on ordinary scales essentially never is), the bound is
  # used, as for every effect size, and the rule is evaluated; the note
  # says that this happened and that it is a matter of units.
  bound_note <- .note_bound_applied(
    effect = "the contrast is in the units of the outcome",
    cause = paste0("this is what happens when the outcome is recorded in ",
                   "units that make its variance very small, or when the ",
                   "contrast weights are very small"),
    remedy = paste0("Rescale the outcome (to standard deviation units, for ",
                    "example) or the weights, and width with them,"),
    topic = "seq_aipe_c"
  )

  .aipe_result(estimate = contrast, raw = raw, n = n, width = width,
               conf_level = conf_level, term = "contrast",
               method = "contrast",
               extra = list(n_per_group = n_per_group,
                            contrast_weights = contrast_weights),
               bound_note = bound_note, unusable_note = flat_note)
}
