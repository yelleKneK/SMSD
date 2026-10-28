#' SMSD: Sequential Methods for Study Design
#'
#' Sequential procedures for planning a study whose sample size cannot be
#' fixed in advance, because the population values that a sample size
#' formula needs (a standard deviation, a correlation, the effect size
#' itself) are unknown. Data are collected in stages, and after each stage a
#' stopping rule computed from the data in hand says whether sampling can
#' stop. The package covers accuracy in parameter estimation (AIPE;
#' Maxwell, Kelley, & Rausch, 2008), in which the goal is a confidence
#' interval as narrow as requested; minimum risk point estimation (MRPE),
#' in which the goal is the smallest combined cost of estimation error and
#' of sampling; a procedure that spends a fixed budget across several
#' groups; and the classical fixed width, bounded risk, and two stage
#' procedures for a mean, a median, or a location parameter.
#'
#' @section Getting Started:
#' Three vignettes introduce the package.
#' \itemize{
#'   \item \code{vignette("introduction", package = "SMSD")}: the overview.
#'   \item \code{vignette("step-by-step-sequential-workflow", package = "SMSD")}:
#'     one study, from its pilot sample to its stopping point.
#'   \item \code{vignette("complete-tour", package = "SMSD")}: the procedures
#'     one by one.
#' }
#' A reader who wants one function to begin with should begin with
#' \code{\link{seq_aipe_smd}}, the standardized mean difference between two
#' groups; its help page has a complete example.
#'
#' @section How the Functions Fit Together:
#' Every procedure has a function that returns its pilot (first stage)
#' sample size, named \code{*_pilot()}, and a check function that takes the
#' data collected so far and evaluates the stopping rule once. The user
#' collects the pilot sample, calls the check function, and, while
#' \code{attr(out, "is_satisfied")} is \code{FALSE}, collects more
#' observations and calls it again. A check function does not prompt,
#' print, or keep anything between calls: the user owns the sampling loop,
#' and each call is a calculation on the data it is given.
#'
#' A check function returns an object of class \code{smsd_check}: a data
#' frame with the columns \code{term} and \code{value}, holding the
#' estimate, its confidence limits where the procedure has an interval,
#' and the sample size in hand beside the one the rule calls for, with the
#' stopping decision and the targets as attributes. See
#' \code{\link{smsd_check_methods}}.
#'
#' @section A Narrow Confidence Interval for an Effect Size:
#' Sequential AIPE: sampling stops at the first sample size that meets
#' the stopping rule, and the confidence interval at that point is no
#' wider than the requested full width. Each effect size has a pilot
#' function, a check function, and a function that simulates the whole
#' procedure.
#' \tabular{lll}{
#'   \strong{Pilot sample size} \tab \strong{Check} \tab \strong{Simulation} \cr
#'   \code{\link{seq_aipe_smd_pilot}} \tab \code{\link{seq_aipe_smd}} \tab \code{\link{seq_aipe_smd_sensitivity}} \cr
#'   \code{\link{seq_aipe_cv_pilot}} \tab \code{\link{seq_aipe_cv}} \tab \code{\link{seq_aipe_cv_sensitivity}} \cr
#'   \code{\link{seq_aipe_reg_coef_pilot}} \tab \code{\link{seq_aipe_reg_coef}} \tab \code{\link{seq_aipe_reg_coef_sensitivity}} \cr
#'   \code{\link{seq_aipe_sm_pilot}} \tab \code{\link{seq_aipe_sm}} \tab \code{\link{seq_aipe_sm_sensitivity}} \cr
#'   \code{\link{seq_aipe_c_pilot}} \tab \code{\link{seq_aipe_c}} \tab \code{\link{seq_aipe_c_sensitivity}} \cr
#'   \code{\link{seq_aipe_r_pilot}} \tab \code{\link{seq_aipe_r}} \tab \code{\link{seq_aipe_r_sensitivity}} \cr
#'   \code{\link{seq_aipe_R2_pilot}} \tab \code{\link{seq_aipe_R2}} \tab \code{\link{seq_aipe_R2_sensitivity}} \cr
#'   \code{\link{seq_aipe_gini_pilot}} \tab \code{\link{seq_aipe_gini}} \tab \code{\link{seq_aipe_gini_sensitivity}} \cr
#' }
#' The effect sizes are the standardized mean difference (\code{smd}), the
#' coefficient of variation (\code{cv}), the simple linear regression slope
#' (\code{reg_coef}), the standardized mean (\code{sm}), and a contrast of
#' group means (\code{c}), from Kelley, Darku, and Chattopadhyay (2018); the
#' Pearson, Kendall, and Spearman correlations (\code{r}) and the squared
#' multiple correlation (\code{R2}), from Kelley, Darku, and Chattopadhyay
#' (2019); and the Gini index (\code{gini}), from Chattopadhyay and De
#' (2016). The squared multiple correlation procedure assumes multivariate
#' normality; the others assume no particular distribution. The point
#' estimates \code{\link{gini_index}}, \code{\link{gini_mean_diff}}, and
#' \code{\link{slr_slope}} are exported for use on their own.
#'
#' @section Minimum Risk Point Estimation:
#' Sampling stops when the estimated risk, the cost of estimation error
#' plus the cost of the observations, is as small as it can be made.
#' \tabular{lll}{
#'   \strong{Pilot sample size} \tab \strong{Check} \tab \strong{Simulation} \cr
#'   \code{\link{seq_mr_cv_pilot}} \tab \code{\link{seq_mr_cv}} \tab \code{\link{seq_mr_cv_sensitivity}} \cr
#'   \code{\link{seq_mr_smd_pilot}} \tab \code{\link{seq_mr_smd}} \tab \code{\link{seq_mr_smd_sensitivity}} \cr
#' }
#' The coefficient of variation follows Chattopadhyay and Kelley (2016). The
#' standardized mean difference follows Chattopadhyay and Kelley (2017),
#' with one sampling cost or a different cost in each group; that article
#' assumes normal populations with a common variance.
#'
#' @section A Budget Across Groups:
#' A contrast of several group means is estimated as precisely as a fixed
#' budget allows when the groups differ in what an observation costs and
#' in their unknown variances (Chattopadhyay, Bandyopadhyay, Kelley, &
#' Padalunkal, 2025).
#' \tabular{ll}{
#'   \code{\link{seq_budget_optimal}} \tab The allocation if the standard deviations were known \cr
#'   \code{\link{seq_budget_pilot}} \tab The pilot sample size \cr
#'   \code{\link{seq_budget}} \tab The check on the data collected so far \cr
#'   \code{\link{seq_budget_next}} \tab The group that should receive the next observation \cr
#'   \code{\link{seq_budget_sensitivity}} \tab A simulation of the whole procedure \cr
#'   \code{\link{decide_seq_budget_noninferiority}} \tab The noninferiority test \cr
#'   \code{\link{decide_seq_budget_superiority}} \tab The superiority test \cr
#'   \code{\link{decide_seq_budget_inferiority}} \tab The inferiority test \cr
#'   \code{\link{decide_seq_budget_equivalence}} \tab The equivalence test \cr
#' }
#' The four tests are one-sided and are applied once sampling has ended.
#'
#' @section Classical Procedures for a Mean, a Median, or a Location:
#' These procedures are from Mukhopadhyay and de Silva (2009): confidence
#' intervals of fixed width or bounded length, and point estimates with
#' bounded or minimum risk. The purely sequential ones:
#' \tabular{lll}{
#'   \strong{Procedure} \tab \strong{Pilot sample size} \tab \strong{Check} \cr
#'   Fixed width, normal mean \tab \code{\link{seq_fwd_norm_mean_pilot}} \tab \code{\link{seq_fwd_norm_mean}} \cr
#'   Fixed width, any distribution \tab \code{\link{seq_fwd_dfr_mean_pilot}} \tab \code{\link{seq_fwd_dfr_mean}} \cr
#'   Fixed width, negative exponential \tab \code{\link{seq_fwd_nexp_mean_pilot}} \tab \code{\link{seq_fwd_nexp_mean}} \cr
#'   Bounded length, median \tab \code{\link{seq_blci_dfr_median_pilot}} \tab \code{\link{seq_blci_dfr_median}} \cr
#'   Bounded risk, normal mean \tab \code{\link{seq_br_norm_mean_pilot}} \tab \code{\link{seq_br_norm_mean}} \cr
#'   Minimum risk, normal mean \tab \code{\link{seq_mr_norm_mean_pilot}} \tab \code{\link{seq_mr_norm_mean}} \cr
#'   Minimum risk, exponential mean \tab \code{\link{seq_mr_exp_mean_pilot}} \tab \code{\link{seq_mr_exp_mean}} \cr
#'   Minimum risk, negative exponential \tab \code{\link{seq_mr_nexp_mean_pilot}} \tab \code{\link{seq_mr_nexp_mean}} \cr
#'   Minimum risk, any distribution \tab \code{\link{seq_mr_dfr_mean_arb_pilot}} \tab \code{\link{seq_mr_dfr_mean_arb}} \cr
#' }
#' The two stage ones, in which one pilot sample fixes the total sample
#' size:
#' \tabular{lll}{
#'   \strong{Procedure} \tab \strong{Pilot sample size} \tab \strong{Check} \cr
#'   Fixed width, normal mean (Stein) \tab \code{\link{ts_fwd_norm_mean_pilot}} \tab \code{\link{ts_fwd_norm_mean}} \cr
#'   Fixed width, any distribution \tab \code{\link{ts_fwd_dfr_mean_pilot}} \tab \code{\link{ts_fwd_dfr_mean}} \cr
#'   Fixed width, negative exponential \tab \code{\link{ts_fwd_nexp_mean_pilot}} \tab \code{\link{ts_fwd_nexp_mean}} \cr
#'   Bounded risk, normal mean \tab \code{\link{ts_br_norm_mean_pilot}} \tab \code{\link{ts_br_norm_mean}} \cr
#'   Bounded risk, exponential mean \tab \code{\link{ts_br_exp_mean_pilot}} \tab \code{\link{ts_br_exp_mean}} \cr
#' }
#' \code{\link{seq_fwd_norm_mean_optimal}} gives the sample size that a
#' known variance would give for the fixed width interval for a normal
#' mean, and \code{\link{power_z}} gives the Type II error rate and the
#' sample size of a one sample, one-sided \emph{z} test. The negative
#' exponential procedures estimate the location of that distribution.
#'
#' @section Studying a Design Before Collecting Data:
#' The \code{*_sensitivity()} functions linked above run a whole procedure,
#' from its pilot sample to its stopping point, many times on simulated
#' data with a known population value. They report the final sample size
#' to expect and how much it varies, and, for the interval procedures, the
#' width of the final interval and how often it contains the population
#' value. The \code{batch} argument of the \code{seq_aipe_*} and
#' \code{seq_mr_*} sensitivity functions shows what it costs to check the
#' rule only after every several observations. Each help page lists the
#' rows of its summary, which prints through
#' \code{\link{smsd_sensitivity_methods}}.
#'
#' @section Reading and Reporting a Result:
#' Every check object has methods for \code{print()}, \code{summary()},
#' \code{format()}, and \code{as.data.frame()}, and for the \code{tidy()}
#' and \code{glance()} generics of the \pkg{generics} package; see
#' \code{\link{smsd_check_methods}}. \pkg{SMSD} does not export those
#' two generics, so write \code{generics::tidy(out)} and
#' \code{generics::glance(out)}, or attach \pkg{generics} or \pkg{broom}
#' first.
#' \code{\link{seq_path}} applies a check function to data already in hand
#' at each successive sample size and returns the whole stopping path;
#' \code{\link{plot.smsd_path}} draws it, and
#' \code{\link{smsd_path_methods}} documents its other methods.
#' \code{\link{seq_report}} writes the paragraph, with its reference, that
#' reports a check or a path in a methods or results section.
#'
#' @section What the Procedures Guarantee:
#' The \code{seq_aipe_*} procedures deliver the requested width by
#' construction: whenever the stopping rule is met, the interval is no
#' wider than \code{width}. Coverage is a different matter. It is an
#' asymptotic property, which approaches the confidence level as the
#' target narrows and the final sample grows, and it can fall well short
#' of the confidence level when the target is wide, when the relationship
#' studied is strong, or when a squared multiple correlation has many
#' predictors. In 2,000 simulated studies each at the 95\% level, with
#' normal data unless noted, the final interval contained the population
#' value in 0.94 of the studies for a standardized mean difference of 0.5
#' at a width of 0.8, in 0.78 for a Pearson correlation of 0.4 and in 0.70
#' for one of 0.8 at a width of 0.4, in 0.84 for a regression slope of 0.5
#' at a width of 0.4, in 0.79 for a coefficient of variation of 0.5 at a
#' width of 0.2 with lognormal data, and in 0.36 for a squared multiple
#' correlation of 0.5 with 15 predictors at a width of 0.3 (0.85 with 5
#' predictors). These figures have Monte Carlo standard errors between
#' 0.005 and 0.011, so the figure that another help page reports for the
#' same design, from another run, can differ in the second decimal. The
#' help page of each check function states what its guarantee is, exact or
#' asymptotic; the sensitivity functions measure it for the design at
#' hand.
#'
#' @section An Unusable Variance Estimate:
#' In a small sample the estimate of the variance that a stopping rule
#' needs can be zero, negative, or not finite. A check function then
#' reports the stopping rule as not met and attaches a note, and its
#' attribute \code{variance_unusable} is \code{TRUE}, so sampling
#' continues; see \code{\link{smsd_check_methods}}. This is a
#' deliberate departure from the stopping rules as printed in the source
#' articles, under which a lower bound would stand in for such an estimate
#' and the rule would be met at once with an interval far narrower than
#' requested. A positive estimate, however small, is used as the articles
#' prescribe. (In \code{\link{ts_fwd_norm_mean}},
#' \code{\link{ts_fwd_dfr_mean}}, \code{\link{ts_fwd_nexp_mean}}, and
#' \code{\link{ts_br_norm_mean}}, a first stage with no variation is an
#' error, because there the first stage variance or scale alone fixes the
#' total sample size. \code{\link{ts_br_exp_mean}} uses the first stage
#' mean and is not affected.)
#'
#' @section Relation to DMAR:
#' \pkg{SMSD} depends on the package \pkg{DMAR}, which is attached with
#' it and provides the fixed sample size counterparts of the confidence
#' interval procedures: \code{\link[DMAR]{ss_aipe_smd}},
#' \code{\link[DMAR]{ss_aipe_r}}, \code{\link[DMAR]{ss_aipe_R2}},
#' \code{\link[DMAR]{ss_aipe_cv}}, \code{\link[DMAR]{ss_aipe_reg_coef}},
#' \code{\link[DMAR]{ss_aipe_sm}}, and \code{\link[DMAR]{ss_aipe_c}}
#' (there is none for the Gini index). If the population values are known
#' well enough to fix the sample size in advance, use \pkg{DMAR}; if they
#' are not, use \pkg{SMSD}.
#'
#' @references
#' Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
#' (2025). A sequential approach for noninferiority or equivalence of a
#' linear contrast under cost constraints. \emph{Psychological Methods},
#' \emph{30}(2), 425--439. \doi{10.1037/met0000570}
#'
#' Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
#' pre-specified error bound. \emph{Econometrics}, \emph{4}(3), 30.
#' \doi{10.3390/econometrics4030030}
#'
#' Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient
#' of variation with minimum risk: A sequential method for minimizing
#' sampling error and study cost. \emph{Multivariate Behavioral
#' Research}, \emph{51}(5), 627--648. \doi{10.1080/00273171.2016.1203279}
#'
#' Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized
#' mean difference with minimum risk: Maximizing accuracy and minimizing
#' cost with sequential estimation. \emph{Psychological Methods},
#' \emph{22}(1), 94--113. \doi{10.1037/met0000089}
#'
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#' parameter estimation for a general class of effect sizes: A
#' sequential approach. \emph{Psychological Methods}, \emph{23}(2),
#' 226--243. \doi{10.1037/met0000127}
#'
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
#' accuracy in parameter estimation for population correlation
#' coefficients. \emph{Psychological Methods}, \emph{24}(4), 492--515.
#' \doi{10.1037/met0000203}
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#' for statistical power and accuracy in parameter estimation. \emph{Annual
#' Review of Psychology}, \emph{59}, 537--563.
#' \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @keywords package
"_PACKAGE"
