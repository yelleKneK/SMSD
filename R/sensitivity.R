#' Monte Carlo Sensitivity Analysis for Sequential AIPE Procedures
#'
#' Each function runs a sequential accuracy in parameter estimation (AIPE)
#' procedure from its pilot sample to its stopping point \code{G} times,
#' against a data generating mechanism whose true effect size is known, and
#' summarizes the final sample size, the width of the final confidence
#' interval, and how often that interval contains the true value. Use them
#' before data collection to see what a planned study should expect: how
#' large the sample is likely to become, how variable that is, and whether
#' the confidence interval can be trusted at the planned width.
#'
#' @details
#' One call does what the pilot, collect, and check loop does by hand: it
#' draws the pilot sample that the matching pilot function gives, evaluates
#' the stopping rule with the matching check function, adds \code{batch}
#' observations, and repeats until the rule is met or the sample reaches
#' \code{max_n}.
#'
#' \strong{Data generating mechanisms.} The true value is exact for each
#' mechanism.
#' \itemize{
#'   \item \code{seq_aipe_smd_sensitivity()}: two independent groups with
#'     unit variance, the first with mean \code{true_value} and the second
#'     with mean 0.
#'   \item \code{seq_aipe_r_sensitivity()}: pairs \eqn{(x, y)} with
#'     \eqn{y = \rho x + \sqrt{1 - \rho^2}\, e}{y = rho x + sqrt(1 - rho^2) e},
#'     where \eqn{x} and \eqn{e} are independent with unit variance and
#'     \eqn{\rho}{rho} is \code{true_value}. The procedure studied is the
#'     Pearson correlation.
#'   \item \code{seq_aipe_reg_coef_sensitivity()}:
#'     \eqn{y = \beta x + e}{y = beta x + e}, where \eqn{x} and \eqn{e} are
#'     independent with unit variance and \eqn{\beta}{beta} is
#'     \code{true_value}.
#'   \item \code{seq_aipe_sm_sensitivity()}: one sample with mean
#'     \code{true_value} and unit variance.
#'   \item \code{seq_aipe_c_sensitivity()}: independent groups with means
#'     \code{true_means} and standard deviations \code{true_sds}. The true
#'     value is the contrast \eqn{\sum_k c_k \mu_k}{sum of c_k mu_k}. The
#'     groups are balanced, as \code{\link{seq_aipe_c}} requires: each
#'     addition is \code{batch} observations from every group.
#'   \item \code{seq_aipe_cv_sensitivity()}: gamma data with shape
#'     \eqn{1/\kappa^2}{1 / kappa^2}, or lognormal data with log standard
#'     deviation \eqn{\sqrt{\log(1 + \kappa^2)}}{sqrt(log(1 + kappa^2))},
#'     where \eqn{\kappa}{kappa} is \code{true_value}.
#'   \item \code{seq_aipe_gini_sensitivity()}: lognormal data with log
#'     standard deviation
#'     \eqn{\sqrt{2}\,\Phi^{-1}((g + 1)/2)}{sqrt(2) qnorm((g + 1) / 2)},
#'     where \eqn{g} is \code{true_value} and \eqn{\Phi^{-1}}{qnorm} is the
#'     standard normal quantile function.
#' }
#' For the first five, \code{distribution = "normal"} draws the unit
#' variance variables from the standard normal distribution and
#' \code{distribution = "skewed"} draws them from a chi square distribution
#' with 3 degrees of freedom, standardized to mean 0 and variance 1
#' (skewness about 1.63).
#'
#' \strong{Batch size.} With \code{batch = 1}, the default, the stopping
#' rule is checked after every observation. That is the purely sequential
#' procedure, and it is what the simulations of Kelley, Darku, and
#' Chattopadhyay (2018, 2019) use; \code{batch} is the quantity they call
#' \eqn{m'} (2018, Figure 1). A study that can check only after every 5, 10,
#' or 20 observations stops later than the first sample size at which the
#' rule would have been met; set \code{batch} to that number to see the
#' overshoot to expect. For a standardized mean difference of 0.5 at
#' \code{width = 0.4}, for example, the mean final sample size per group was
#' 199.1 with \code{batch = 1}, 201.2 with 5, 203.6 with 10, and 209.5 with
#' 20 (4,700 to 6,000 replications each, standard errors about 0.1), and
#' coverage stayed between 0.949 and 0.952: checking less often cost about
#' half a batch per group.
#'
#' \strong{Coverage is an asymptotic property.} The procedures deliver the
#' requested width by construction: the final interval is never wider than
#' \code{width} unless \code{max_n} ends the replication first. Coverage is
#' a different matter. The articles prove that it approaches
#' \code{conf_level} as \code{width} shrinks (Kelley, Darku, &
#' Chattopadhyay, 2018, p. 232; 2019, p. 495), that is, as the sample that
#' the stopping rule asks for grows. When that sample is small, coverage can
#' fall well short of \code{conf_level}. How small it is depends on the true
#' value as well as on \code{width}: a strong correlation is estimated with
#' so little sampling error that the rule is met in a small sample even at a
#' width that would be narrow for a moderate correlation.
#'
#' The table gives what these functions measured at
#' \code{conf_level = 0.95}, with the defaults and normal data unless a
#' distribution is named: the mean final sample size (the row
#' \code{mean_n}), the coverage, and the percentage of replications in which
#' the variance estimate was unusable at one or more checks (the row
#' \code{pct_variance_unusable}, explained below). Each row rests on
#' 800 to 34,000 replications; the Monte Carlo standard error of the
#' coverage is below 0.008 in every row but three with a mean \emph{n} above
#' 800 (the correlation of 0.9 at 0.025, the slope of 2 at 0.1, and the
#' lognormal coefficient of variation at 0.05), where it is 0.009 to 0.012.
#' \tabular{lrrrr}{
#'   Effect size and true value \tab Width \tab Mean \emph{n} \tab
#'     Coverage \tab Percent unusable \cr
#'   Standardized mean difference, 0.5 \tab
#'     0.8 \tab 51 \tab 0.95 \tab <0.1 \cr
#'   Standardized mean difference, 0.5 \tab
#'     0.4 \tab 199 \tab 0.95 \tab 0.0 \cr
#'   Standardized mean, 0.4 \tab 0.5 \tab 67 \tab 0.95 \tab <0.1 \cr
#'   Standardized mean, 0.4 \tab 0.25 \tab 267 \tab 0.95 \tab 0.0 \cr
#'   Contrast, 0.5 \tab 1 \tab 23 \tab 0.94 \tab 0.0 \cr
#'   Contrast, 0.5 \tab 0.5 \tab 92 \tab 0.95 \tab 0.0 \cr
#'   Correlation, 0.4 \tab 0.4 \tab 55 \tab 0.78 \tab 16.5 \cr
#'   Correlation, 0.4 \tab 0.2 \tab 263 \tab 0.93 \tab 0.9 \cr
#'   Correlation, 0.4 \tab 0.1 \tab 1,079 \tab 0.95 \tab 0.0 \cr
#'   Correlation, 0.9 \tab 0.1 \tab 75 \tab 0.45 \tab 79.5 \cr
#'   Correlation, 0.9 \tab 0.05 \tab 153 \tab 0.63 \tab 40.4 \cr
#'   Correlation, 0.9 \tab 0.025 \tab 819 \tab 0.90 \tab 1.3 \cr
#'   Regression slope, 0.5 \tab 0.4 \tab 80 \tab 0.83 \tab 20.1 \cr
#'   Regression slope, 0.5 \tab 0.2 \tab 377 \tab 0.94 \tab 1.3 \cr
#'   Regression slope, 2 \tab 0.2 \tab 217 \tab 0.61 \tab 65.3 \cr
#'   Regression slope, 2 \tab 0.1 \tab 1,382 \tab 0.87 \tab 22.4 \cr
#'   Coefficient of variation, 0.5, gamma \tab
#'     0.2 \tab 52 \tab 0.88 \tab 0.0 \cr
#'   Coefficient of variation, 0.5, gamma \tab
#'     0.1 \tab 217 \tab 0.92 \tab 0.0 \cr
#'   Coefficient of variation, 0.5, lognormal \tab
#'     0.2 \tab 62 \tab 0.79 \tab 0.0 \cr
#'   Coefficient of variation, 0.5, lognormal \tab
#'     0.1 \tab 322 \tab 0.86 \tab 0.0 \cr
#'   Coefficient of variation, 0.5, lognormal \tab
#'     0.05 \tab 1,561 \tab 0.90 \tab 0.0 \cr
#'   Gini index, 0.3 \tab 0.2 \tab 30 \tab 0.89 \tab 0.0 \cr
#'   Gini index, 0.3 \tab 0.1 \tab 88 \tab 0.92 \tab 0.0}
#' (The contrast has weights 1, -0.5, and -0.5, means 0.5, 0, and 0, and
#' unit standard deviations. The coverage for the Gini index is that of the
#' interval \code{\link{seq_aipe_gini}} reports. The fixed width interval of
#' Chattopadhyay and De, 2016, the estimate plus and minus half the width,
#' contained the index in 0.986 and 0.960 of the same studies; see
#' \code{\link{seq_aipe_gini}}.) These values hold for the true values shown,
#' not for an effect size in general: the rows for a correlation of 0.9 and
#' for a slope of 2 are far from those for a correlation of 0.4 and a slope
#' of 0.5, at widths that are no wider. Skewed data can lower coverage
#' further: for a slope of 0.5 at \code{width = 0.4} it was 0.71 with the
#' skewed data.
#'
#' Read the \code{coverage} row before relying on an interval at a planned
#' width, and run the function at each true value that is plausible for the
#' study, the strongest included. When coverage is low, plan a narrower
#' \code{width} (a larger study), raise the minimum pilot size \code{m_0},
#' or both. With \code{m_0 = 40} at \code{width = 0.4}, coverage rose from
#' 0.78 to 0.90 for a correlation of 0.4 and from 0.83 to 0.92 for a slope
#' of 0.5, while the mean final sample size rose from 55 to 64 and from 80
#' to 90. A larger pilot sample is a weaker remedy when the true value is
#' strong: for a correlation of 0.9 at \code{width = 0.1}, \code{m_0 = 100}
#' raised coverage only from 0.45 to 0.69.
#'
#' \strong{When the variance estimate is unusable.} In a small sample the
#' raw estimate of the asymptotic variance can be zero, negative, or not
#' finite: no variance has been estimated. The check functions report such a
#' check as not met and set the attribute \code{variance_unusable} to
#' \code{TRUE}, so sampling continues (see \code{\link{smsd_check_methods}};
#' this is a deliberate departure from the stopping rule as printed in the
#' articles, under which the procedure would stop at once with an interval
#' far narrower than the data support). A positive estimate, however small,
#' is used as the articles prescribe and is not counted here.
#'
#' The row \code{pct_variance_unusable} is the percentage of replications in
#' which this happened at one or more checks; it is the last column of the
#' table. Like coverage, it depends on the true value as well as on the
#' width, and it is not confined to wide targets or to one effect size. It
#' was 16 percent for a correlation of 0.4 at \code{width = 0.4} but 79
#' percent for a correlation of 0.9 at \code{width = 0.1}; 20 percent for a
#' slope of 0.5 at \code{width = 0.4} but 65 percent for a slope of 2 at
#' \code{width = 0.2}; and 15 percent for a standardized mean of 1 at
#' \code{width = 1}, where coverage was 0.86. With the skewed data it was 37
#' percent for the slope of 0.5 at \code{width = 0.4}.
#'
#' Continuing past such a check does not by itself make the final interval
#' trustworthy. The estimate on which a replication finally stops can still
#' be far too small, which makes the final interval too narrow: for the
#' correlation of 0.9 at \code{width = 0.1}, the final intervals had a mean
#' width of 0.028, about a third of the 0.086 that the true variance gives
#' at the mean final sample size, and their coverage was 0.45. With a very
#' strong relationship the estimate can also stay unusable for a long time.
#' For a correlation of 0.99 at \code{width = 0.1} (a pilot sample of 40),
#' 20 replications ended at a mean of 7,500 pairs, each taking about 6.5
#' seconds on the computer used for these runs, and 1 of their 20 final
#' intervals contained the true value. Set \code{max_n} when exploring such
#' conditions, and read \code{pct_hit_max_n}.
#'
#' \strong{Monte Carlo error.} Every summary is an estimate from \code{G}
#' replications. The standard error of \code{coverage} is about
#' \eqn{\sqrt{p(1 - p)/G}}{sqrt(p (1 - p) / G)} for a coverage of \eqn{p}:
#' 0.007 at the default \code{G = 1000} and 0.002 at \code{G = 10000} when
#' \eqn{p} is near 0.95. The standard error of \code{mean_n} is \code{sd_n}
#' divided by \eqn{\sqrt{G}}{sqrt(G)}.
#'
#' @param true_value The population value of the effect size under the data
#'   generating mechanism: the standardized mean difference (any finite
#'   number), the correlation (strictly between -1 and 1), the coefficient
#'   of variation (positive), the regression slope (any finite number), the
#'   standardized mean (any finite number), or the Gini index (strictly
#'   between 0 and 1). Within those ranges, a value at which the generated
#'   data would have too little variation to be used is an error. The check
#'   functions take a sample whose values agree to ten significant digits to
#'   have no variation, and the sensitivity functions keep the generated
#'   data a factor of 100 clear of that: the standardized mean difference
#'   and the standardized mean must be less than 1e8 in absolute value (the
#'   data have unit variance), the coefficient of variation must be greater
#'   than 1e-8, and the Gini index greater than about 5.6e-9. The regression
#'   slope must be less than \code{1 / .Machine$double.eps} (about 4.5e15)
#'   in absolute value, beside which the unit variance of the errors cannot
#'   be represented. A coefficient of variation or Gini index so large that
#'   its generating parameter (see Details) is not a positive finite number
#'   is an error as well.
#' @param true_means For the contrast: the population means, one per group
#'   (at least two groups), each less than 1e8 times its group's standard
#'   deviation in absolute value (for the reason given under
#'   \code{true_value}).
#' @param contrast_weights For the contrast: the contrast coefficients
#'   \eqn{c_k}, one per group and not all zero, as in
#'   \code{\link{seq_aipe_c}}. \code{true_means}, \code{contrast_weights},
#'   and \code{true_sds} (when it has one entry per group) are matched to
#'   the groups by position; if more than one of them is named, the names
#'   must be the same and in the same order, or the function stops.
#' @param true_sds For the contrast: the population standard deviations,
#'   either one positive number shared by every group (default 1) or one
#'   per group. Their squares must be positive finite numbers.
#' @param width The target full width of the confidence interval (the
#'   quantity \eqn{\omega}{omega} of Kelley, Darku, and Chattopadhyay, 2018,
#'   2019), in the units of the effect size. For the contrast these are the
#'   units of the response. For the correlation it must be less than 2, and
#'   for the Gini index less than 1, the lengths of their ranges.
#' @param conf_level The confidence level, strictly between 0 and 1
#'   (default 0.95).
#' @param distribution The data generating mechanism (see Details).
#'   \code{"normal"} or \code{"skewed"} for the standardized mean
#'   difference, the correlation, the regression slope, the standardized
#'   mean, and the contrast; \code{"gamma"} or \code{"lognormal"} for the
#'   coefficient of variation. The Gini index function always generates
#'   lognormal data and has no such argument.
#' @param G The number of Monte Carlo replications, a whole number of at
#'   least 1 (default 1000). See Details for the Monte Carlo error this
#'   implies.
#' @param seed \code{NULL} (the default) or a single whole number. With a
#'   number the result is reproducible, and the caller's random number
#'   stream is left exactly as it was found. With \code{NULL} the draws come
#'   from the current stream.
#' @param m_0 The minimum pilot sample size, a whole number of at least 4
#'   (default 4), passed to the matching \code{seq_aipe_*_pilot()} function.
#'   The pilot sample size is the larger of \code{m_0} and
#'   \eqn{\lceil 2 z/\omega \rceil}{ceiling(2 z / omega)}, where \eqn{z} is
#'   the standard normal quantile for \code{conf_level} (1.96 at 0.95).
#'   Like every sample size here (\code{max_n}, \code{batch}, and the rows
#'   \code{mean_n}, \code{median_n}, and \code{sd_n}), it counts
#'   observations per group for the standardized mean difference and the
#'   contrast, pairs for the correlation and the regression slope, and
#'   observations otherwise.
#' @param max_n The largest sample size a replication may reach, in the
#'   unit described under \code{m_0} (per group, for example, for the
#'   standardized mean difference and the contrast, not in total): a whole
#'   number no smaller than the pilot sample size (default 100,000). A
#'   replication that has not met the stopping rule at \code{max_n} is
#'   stopped there, exactly, whatever the batch size (the last batch is
#'   shortened). It is an error for \code{max_n} to be below the pilot
#'   sample size. The row \code{pct_hit_max_n} reports how often the cap,
#'   not the rule, ended a replication; when it is above zero the other
#'   rows describe a truncated procedure and should be read with care.
#' @param batch The number of observations added between checks of the
#'   stopping rule after the pilot sample, a whole number of at least 1
#'   (default 1; see Details), in the unit described under \code{m_0}: for
#'   the standardized mean difference and the contrast it is the number
#'   added to every group; for the correlation and the regression slope it
#'   is the number of pairs.
#'
#' @return An object of class \code{c("smsd_sensitivity", "data.frame")}: a
#'   data frame with columns \code{term} (character) and \code{value}
#'   (numeric) that can be used as any data frame, with a \code{print()}
#'   method that shows every value in fixed notation (see
#'   \code{\link{smsd_sensitivity_methods}}). It has these 13 rows, in this
#'   order:
#'   \describe{
#'     \item{\code{true_value}}{The population value of the effect size. For
#'       the contrast, the contrast of \code{true_means}.}
#'     \item{\code{target_width}}{The requested full width, \code{width}.}
#'     \item{\code{conf_level}}{The confidence level, a proportion.}
#'     \item{\code{replications}}{The number of replications run, \code{G}.}
#'     \item{\code{mean_n}, \code{median_n}, \code{sd_n}}{The mean, median,
#'       and standard deviation of the final sample size across
#'       replications: observations for the coefficient of variation, the
#'       standardized mean, and the Gini index; pairs for the correlation
#'       and the regression slope; observations per group for the
#'       standardized mean difference and the contrast. \code{sd_n} is
#'       \code{NA} when \code{G} is 1.}
#'     \item{\code{mean_width}, \code{median_width}}{The mean and median
#'       full width of the final confidence interval, in the units of the
#'       effect size, over the replications that ended with an interval
#'       (see below); \code{NA} when none did.}
#'     \item{\code{pct_width_within}}{The percentage (0 to 100) of
#'       replications whose final interval is no wider than \code{width}.}
#'     \item{\code{coverage}}{The proportion (0 to 1) of replications whose
#'       final interval contains the true value.}
#'     \item{\code{pct_hit_max_n}}{The percentage (0 to 100) of
#'       replications stopped by \code{max_n} before the stopping rule was
#'       met.}
#'     \item{\code{pct_variance_unusable}}{The percentage (0 to 100) of
#'       replications in which, at one or more checks, the check function
#'       reported its variance estimate as unusable: the raw estimate was
#'       zero, negative, or not finite, so that no variance had been
#'       estimated and the stopping rule was reported as not met at that
#'       check. A replication is counted once however many of its checks
#'       were of this kind, and the check at which \code{max_n} ended it
#'       counts as well. A positive estimate, however small, is usable and
#'       is not counted.}
#'   }
#'   A replication that the stopping rule ends always has an interval. A
#'   replication that \code{max_n} ends can finish on a check that gives no
#'   interval, because the estimate is not finite or the variance estimate
#'   is unusable (the check function then reports its limits as \code{NA}).
#'   Such a replication counts as not covering and as not within the target
#'   width, so that \code{coverage} and \code{pct_width_within} remain
#'   proportions of all \code{G} replications; it is left out of
#'   \code{mean_width} and \code{median_width}; and it is counted like any
#'   other in the remaining rows, \code{pct_hit_max_n} among them. No row
#'   becomes \code{NA} on this account, except the two width rows when not
#'   one replication ended with an interval.
#'
#' @references
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#' parameter estimation for a general class of effect sizes: A sequential
#' approach. \emph{Psychological Methods}, \emph{23}(2), 226--243.
#' \doi{10.1037/met0000127}
#'
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential accuracy
#' in parameter estimation for population correlation coefficients.
#' \emph{Psychological Methods}, \emph{24}(4), 492--515.
#' \doi{10.1037/met0000203}
#'
#' Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
#' pre-specified error bound. \emph{Econometrics}, \emph{4}(3), 30.
#' \doi{10.3390/econometrics4030030}
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#' for statistical power and accuracy in parameter estimation. \emph{Annual
#' Review of Psychology}, \emph{59}, 537--563.
#' \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu},
#' Francis Bilson Darku \email{fbilsond@@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{seq_aipe_smd}}, \code{\link{seq_aipe_r}},
#'     \code{\link{seq_aipe_cv}}, and \code{\link{seq_aipe_reg_coef}}: four
#'     of the check functions these run;
#'   \item \code{\link{seq_aipe_sm}}, \code{\link{seq_aipe_gini}}, and
#'     \code{\link{seq_aipe_c}}: the other three;
#'   \item \code{\link{smsd_sensitivity_methods}}: how the summary prints;
#'   \item \code{vignette("sensitivity-analysis", package = "SMSD")}: a
#'     worked sensitivity analysis;
#'   \item \code{\link[DMAR]{ss_aipe_smd_sensitivity}}: in \pkg{DMAR}, the
#'     sensitivity function for fixed sample size planning for the
#'     standardized mean difference;
#'   \item \code{\link[DMAR]{ss_aipe_r_sensitivity}}: the same for the
#'     correlation;
#'   \item \code{\link[DMAR]{ss_aipe_cv_sensitivity}}: for the coefficient
#'     of variation;
#'   \item \code{\link[DMAR]{ss_aipe_reg_coef_sensitivity}}: for the
#'     regression slope;
#'   \item \code{\link[DMAR]{ss_aipe_sm_sensitivity}}: for the standardized
#'     mean;
#'   \item \code{\link[DMAR]{ss_aipe_c_sensitivity}}: for the contrast.
#' }
#'
#' @name seq_aipe_sensitivity
NULL

#' @rdname seq_aipe_sensitivity
#' @examples
#' # G is kept very small in these examples so that they run in a fraction
#' # of a second, and the targets are wide for the same reason. A planning
#' # result that will be reported deserves G = 10000, at the width the
#' # study will use.
#'
#' # A training program compared with a control group; the true
#' # standardized mean difference in engagement is 0.5.
#' seq_aipe_smd_sensitivity(true_value = 0.5, width = 1.2, G = 10, seed = 113)
#'
#' @family seq_aipe_smd
#' @export
seq_aipe_smd_sensitivity <- function(true_value, width, conf_level = 0.95,
                                     distribution = c("normal", "skewed"),
                                     G = 1000, seed = NULL, m_0 = 4,
                                     max_n = 1e5L, batch = 1L) {
  .sensitivity_need(true_value, "true_value",
                    "the population standardized mean difference to simulate")
  distribution <- .check_choice(distribution, c("normal", "skewed"),
                                "distribution")
  a <- .sensitivity_check_args(true_value, width, conf_level, G, m_0,
                               max_n, batch, seed)
  delta <- a$true_value
  .sensitivity_check_location(delta)
  rand <- .sensitivity_std(distribution)
  draw <- function(data, k) {
    g1 <- rand(k) + delta
    g2 <- rand(k)
    if (is.null(data)) list(g1 = g1, g2 = g2)
    else list(g1 = c(data$g1, g1), g2 = c(data$g2, g2))
  }
  pilot <- seq_aipe_smd_pilot(width = a$width, conf_level = a$conf_level,
                              m_0 = a$m_0)
  evaluate <- function(d) {
    seq_aipe_smd(d$g1, d$g2, width = a$width, conf_level = a$conf_level)
  }
  .sensitivity_run(a$G, a$seed, pilot, draw, evaluate, delta, a$width,
                   a$conf_level, a$max_n, a$batch)
}

#' @rdname seq_aipe_sensitivity
#' @examples
#' # Job satisfaction and performance ratings with a true correlation of
#' # 0.4. Ten replications say little about coverage (its standard error
#' # here is about 0.15); in 10,000 replications at this very wide target
#' # it was 0.65, far below 0.95.
#' seq_aipe_r_sensitivity(true_value = 0.4, width = 0.8, G = 10, seed = 113)
#'
#' @family seq_aipe_r
#' @export
seq_aipe_r_sensitivity <- function(true_value, width, conf_level = 0.95,
                                   distribution = c("normal", "skewed"),
                                   G = 1000, seed = NULL, m_0 = 4,
                                   max_n = 1e5L, batch = 1L) {
  .sensitivity_need(true_value, "true_value",
                    "the population correlation to simulate")
  distribution <- .check_choice(distribution, c("normal", "skewed"),
                                "distribution")
  .check_number(true_value, "true_value", lower = -1, upper = 1,
                lower_open = TRUE, upper_open = TRUE)
  a <- .sensitivity_check_args(true_value, width, conf_level, G, m_0,
                               max_n, batch, seed, width_range = 2)
  rho <- a$true_value
  rand <- .sensitivity_std(distribution)
  draw <- function(data, k) {
    z1 <- rand(k)
    z2 <- rand(k)
    x <- z1
    y <- rho * z1 + sqrt(1 - rho^2) * z2
    if (is.null(data)) list(x = x, y = y)
    else list(x = c(data$x, x), y = c(data$y, y))
  }
  pilot <- seq_aipe_r_pilot(width = a$width, conf_level = a$conf_level,
                            m_0 = a$m_0)
  evaluate <- function(d) {
    seq_aipe_r(d$x, d$y, width = a$width, conf_level = a$conf_level,
               method = "pearson")
  }
  .sensitivity_run(a$G, a$seed, pilot, draw, evaluate, rho, a$width,
                   a$conf_level, a$max_n, a$batch)
}

#' @rdname seq_aipe_sensitivity
#' @examples
#' # Response times with a true coefficient of variation of 0.5.
#' seq_aipe_cv_sensitivity(true_value = 0.5, width = 0.5, G = 10, seed = 113)
#'
#' @family seq_aipe_cv
#' @export
seq_aipe_cv_sensitivity <- function(true_value, width, conf_level = 0.95,
                                    distribution = c("gamma", "lognormal"),
                                    G = 1000, seed = NULL, m_0 = 4,
                                    max_n = 1e5L, batch = 1L) {
  .sensitivity_need(true_value, "true_value",
                    "the population coefficient of variation to simulate")
  distribution <- .check_choice(distribution, c("gamma", "lognormal"),
                                "distribution")
  .check_positive(true_value, "true_value")
  a <- .sensitivity_check_args(true_value, width, conf_level, G, m_0,
                               max_n, batch, seed)
  cv <- a$true_value
  # The coefficient of variation is the standard deviation of the generated
  # data over their mean.
  .sensitivity_check_rel_sd(cv, format(1 / .sensitivity_mean_limit()), cv)
  if (distribution == "gamma") {
    shape <- .sensitivity_check_parameter(
      .sensitivity_gamma_shape_cv(cv), "gamma",
      "shape, 1 / true_value^2", cv)
    rand <- function(k) stats::rgamma(k, shape = shape, scale = 1)
  } else {
    sdlog <- .sensitivity_check_parameter(
      .sensitivity_sdlog_cv(cv), "lognormal",
      "log standard deviation, sqrt(log(1 + true_value^2))", cv)
    rand <- function(k) stats::rlnorm(k, meanlog = 0, sdlog = sdlog)
  }
  draw <- function(data, k) {
    x <- rand(k)
    if (is.null(data)) x else c(data, x)
  }
  pilot <- seq_aipe_cv_pilot(width = a$width, conf_level = a$conf_level,
                             m_0 = a$m_0)
  evaluate <- function(d) {
    seq_aipe_cv(d, width = a$width, conf_level = a$conf_level)
  }
  .sensitivity_run(a$G, a$seed, pilot, draw, evaluate, cv, a$width,
                   a$conf_level, a$max_n, a$batch)
}

#' @rdname seq_aipe_sensitivity
#' @examples
#' # Hours of tutoring predicting a test score, with a true slope of 0.5.
#' seq_aipe_reg_coef_sensitivity(true_value = 0.5, width = 0.8, G = 10,
#'                               seed = 113)
#'
#' @family seq_aipe_reg_coef
#' @export
seq_aipe_reg_coef_sensitivity <- function(true_value, width,
                                          conf_level = 0.95,
                                          distribution = c("normal", "skewed"),
                                          G = 1000, seed = NULL, m_0 = 4,
                                          max_n = 1e5L, batch = 1L) {
  .sensitivity_need(true_value, "true_value",
                    "the population regression slope to simulate")
  distribution <- .check_choice(distribution, c("normal", "skewed"),
                                "distribution")
  a <- .sensitivity_check_args(true_value, width, conf_level, G, m_0,
                               max_n, batch, seed)
  beta <- a$true_value
  .sensitivity_check_slope(beta)
  rand <- .sensitivity_std(distribution)
  draw <- function(data, k) {
    x <- rand(k)
    e <- rand(k)
    y <- beta * x + e
    if (is.null(data)) list(x = x, y = y)
    else list(x = c(data$x, x), y = c(data$y, y))
  }
  pilot <- seq_aipe_reg_coef_pilot(width = a$width,
                                   conf_level = a$conf_level, m_0 = a$m_0)
  evaluate <- function(d) {
    seq_aipe_reg_coef(d$x, d$y, width = a$width, conf_level = a$conf_level)
  }
  .sensitivity_run(a$G, a$seed, pilot, draw, evaluate, beta, a$width,
                   a$conf_level, a$max_n, a$batch)
}

#' @rdname seq_aipe_sensitivity
#' @examples
#' # A change score with a true standardized mean of 0.4.
#' seq_aipe_sm_sensitivity(true_value = 0.4, width = 1, G = 10, seed = 113)
#'
#' @family seq_aipe_sm
#' @export
seq_aipe_sm_sensitivity <- function(true_value, width, conf_level = 0.95,
                                    distribution = c("normal", "skewed"),
                                    G = 1000, seed = NULL, m_0 = 4,
                                    max_n = 1e5L, batch = 1L) {
  .sensitivity_need(true_value, "true_value",
                    "the population standardized mean to simulate")
  distribution <- .check_choice(distribution, c("normal", "skewed"),
                                "distribution")
  a <- .sensitivity_check_args(true_value, width, conf_level, G, m_0,
                               max_n, batch, seed)
  theta <- a$true_value
  .sensitivity_check_location(theta)
  rand <- .sensitivity_std(distribution)
  draw <- function(data, k) {
    x <- rand(k) + theta
    if (is.null(data)) x else c(data, x)
  }
  pilot <- seq_aipe_sm_pilot(width = a$width, conf_level = a$conf_level,
                             m_0 = a$m_0)
  evaluate <- function(d) {
    seq_aipe_sm(d, width = a$width, conf_level = a$conf_level)
  }
  .sensitivity_run(a$G, a$seed, pilot, draw, evaluate, theta, a$width,
                   a$conf_level, a$max_n, a$batch)
}

#' @rdname seq_aipe_sensitivity
#' @examples
#' # Household incomes with a true Gini index of 0.3.
#' seq_aipe_gini_sensitivity(true_value = 0.3, width = 0.4, G = 10, seed = 113)
#'
#' @family seq_aipe_gini
#' @export
seq_aipe_gini_sensitivity <- function(true_value, width, conf_level = 0.95,
                                      G = 1000, seed = NULL, m_0 = 4,
                                      max_n = 1e5L, batch = 1L) {
  .sensitivity_need(true_value, "true_value",
                    "the population Gini index to simulate")
  .check_number(true_value, "true_value", lower = 0, upper = 1,
                lower_open = TRUE, upper_open = TRUE)
  a <- .sensitivity_check_args(true_value, width, conf_level, G, m_0,
                               max_n, batch, seed, width_range = 1)
  gini <- a$true_value
  sdlog <- .sensitivity_sdlog_gini(gini)
  # Lognormal data with log standard deviation s have a standard deviation
  # of sqrt(exp(s^2) - 1) times their mean. The smallest Gini index allowed
  # is that of the lognormal distribution at the limit of that ratio.
  least <- .sensitivity_gini_lognormal(
    .sensitivity_sdlog_cv(1 / .sensitivity_mean_limit()))
  .sensitivity_check_rel_sd(
    sqrt(expm1(sdlog^2)),
    paste0("about ", format(least, digits = 2), " (the Gini index of ",
           "lognormal data with a coefficient of variation of ",
           format(1 / .sensitivity_mean_limit()), ")"),
    gini)
  sdlog <- .sensitivity_check_parameter(
    sdlog, "lognormal",
    "log standard deviation, sqrt(2) * qnorm((true_value + 1) / 2)", gini)
  draw <- function(data, k) {
    x <- stats::rlnorm(k, meanlog = 0, sdlog = sdlog)
    if (is.null(data)) x else c(data, x)
  }
  pilot <- seq_aipe_gini_pilot(width = a$width, conf_level = a$conf_level,
                               m_0 = a$m_0)
  evaluate <- function(d) {
    seq_aipe_gini(d, width = a$width, conf_level = a$conf_level)
  }
  .sensitivity_run(a$G, a$seed, pilot, draw, evaluate, gini, a$width,
                   a$conf_level, a$max_n, a$batch)
}

#' @rdname seq_aipe_sensitivity
#' @examples
#' # A new teaching method against the average of two existing methods,
#' # with exam scores that have a standard deviation of 10. The class sizes
#' # can only be checked after every 5 students per method.
#' seq_aipe_c_sensitivity(true_means = c(75, 70, 70),
#'                        contrast_weights = c(1, -0.5, -0.5),
#'                        width = 12, true_sds = 10, G = 10, seed = 113,
#'                        batch = 5)
#'
#' @family seq_aipe_c
#' @export
seq_aipe_c_sensitivity <- function(true_means, contrast_weights, width,
                                   conf_level = 0.95, true_sds = 1,
                                   distribution = c("normal", "skewed"),
                                   G = 1000, seed = NULL, m_0 = 4,
                                   max_n = 1e5L, batch = 1L) {
  .sensitivity_need(true_means, "true_means",
                    "the population means to simulate, one per group")
  .sensitivity_need(contrast_weights, "contrast_weights",
                    "one contrast weight per group")
  distribution <- .check_choice(distribution, c("normal", "skewed"),
                                "distribution")
  # The entries are matched to the groups by position. The names are read
  # here, before the checks, which return the vectors without them; a single
  # shared standard deviation belongs to no one group, so its name is not
  # compared.
  labels <- list(true_means = names(true_means),
                 contrast_weights = names(contrast_weights),
                 true_sds = if (length(true_sds) > 1L) names(true_sds))
  true_means <- .check_vector(true_means, "true_means",
                              "one population mean per group")
  n_groups <- length(true_means)
  if (n_groups < 2L) {
    .smsd_stop("'true_means' must give the population means of at least 2 ",
               "groups; got 1.")
  }
  contrast_weights <- .check_vector(
    contrast_weights, "contrast_weights", "one contrast weight per group",
    len = n_groups, len_what = "one per group in 'true_means'")
  if (all(contrast_weights == 0)) {
    .smsd_stop("'contrast_weights' must not all be zero: a contrast needs ",
               "at least one group with a nonzero weight.")
  }
  # One standard deviation shared by every group, or one per group.
  true_sds <- .check_vector(
    true_sds, "true_sds",
    "one population standard deviation, or one per group",
    len = if (length(true_sds) == 1L) 1L else n_groups,
    len_what = paste0("one per group in 'true_means', or else a single ",
                      "entry shared by every group"),
    positive = TRUE)
  if (any(!is.finite(true_sds^2) | true_sds^2 == 0)) {
    .smsd_stop("'true_sds' must have squares (the population variances) ",
               "that are positive finite numbers; ",
               sum(!is.finite(true_sds^2) | true_sds^2 == 0), " of ",
               length(true_sds), " would overflow or underflow.")
  }
  true_sds <- rep_len(true_sds, n_groups)
  .budget_check_names(labels)
  .sensitivity_check_location(true_means, true_sds, "true_means")
  psi <- sum(contrast_weights * true_means)
  if (!is.finite(psi)) {
    .smsd_stop("The contrast implied by 'true_means' and ",
               "'contrast_weights' is not finite; check both arguments.")
  }
  a <- .sensitivity_check_args(psi, width, conf_level, G, m_0, max_n,
                               batch, seed)
  rand <- .sensitivity_std(distribution)
  draw <- function(data, k) {
    new <- lapply(seq_len(n_groups), function(j) {
      true_means[j] + true_sds[j] * rand(k)
    })
    if (is.null(data)) new else Map(c, data, new)
  }
  pilot <- seq_aipe_c_pilot(width = a$width, conf_level = a$conf_level,
                            m_0 = a$m_0)
  evaluate <- function(d) {
    seq_aipe_c(d, contrast_weights = contrast_weights, width = a$width,
               conf_level = a$conf_level)
  }
  .sensitivity_run(a$G, a$seed, pilot, draw, evaluate, psi, a$width,
                   a$conf_level, a$max_n, a$batch)
}
