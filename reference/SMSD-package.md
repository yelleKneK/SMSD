# SMSD: Sequential Methods for Study Design

Sequential procedures for planning a study whose sample size cannot be
fixed in advance, because the population values that a sample size
formula needs (a standard deviation, a correlation, the effect size
itself) are unknown. Data are collected in stages, and after each stage
a stopping rule computed from the data in hand says whether sampling can
stop. The package covers accuracy in parameter estimation (AIPE;
Maxwell, Kelley, & Rausch, 2008), in which the goal is a confidence
interval as narrow as requested; minimum risk point estimation (MRPE),
in which the goal is the smallest combined cost of estimation error and
of sampling; a procedure that spends a fixed budget across several
groups; and the classical fixed width, bounded risk, and two stage
procedures for a mean, a median, or a location parameter.

## Getting Started

Three vignettes introduce the package.

- [`vignette("introduction", package = "SMSD")`](https://yelleknek.github.io/SMSD/articles/introduction.md):
  the overview.

- [`vignette("step-by-step-sequential-workflow", package = "SMSD")`](https://yelleknek.github.io/SMSD/articles/step-by-step-sequential-workflow.md):
  one study, from its pilot sample to its stopping point.

- [`vignette("complete-tour", package = "SMSD")`](https://yelleknek.github.io/SMSD/articles/complete-tour.md):
  the procedures one by one.

A reader who wants one function to begin with should begin with
[`seq_aipe_smd`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md),
the standardized mean difference between two groups; its help page has a
complete example.

## How the Functions Fit Together

Every procedure has a function that returns its pilot (first stage)
sample size, named `*_pilot()`, and a check function that takes the data
collected so far and evaluates the stopping rule once. The user collects
the pilot sample, calls the check function, and, while
`attr(out, "is_satisfied")` is `FALSE`, collects more observations and
calls it again. A check function does not prompt, print, or keep
anything between calls: the user owns the sampling loop, and each call
is a calculation on the data it is given.

A check function returns an object of class `smsd_check`: a data frame
with the columns `term` and `value`, holding the estimate, its
confidence limits where the procedure has an interval, and the sample
size in hand beside the one the rule calls for, with the stopping
decision and the targets as attributes. See
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md).

## A Narrow Confidence Interval for an Effect Size

Sequential AIPE: sampling stops at the first sample size that meets the
stopping rule, and the confidence interval at that point is no wider
than the requested full width. Each effect size has a pilot function, a
check function, and a function that simulates the whole procedure.

|  |  |  |
|----|----|----|
| **Pilot sample size** | **Check** | **Simulation** |
| [`seq_aipe_smd_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd_pilot.md) | [`seq_aipe_smd`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md) | [`seq_aipe_smd_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md) |
| [`seq_aipe_cv_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv_pilot.md) | [`seq_aipe_cv`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md) | [`seq_aipe_cv_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md) |
| [`seq_aipe_reg_coef_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef_pilot.md) | [`seq_aipe_reg_coef`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md) | [`seq_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md) |
| [`seq_aipe_sm_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm_pilot.md) | [`seq_aipe_sm`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm.md) | [`seq_aipe_sm_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md) |
| [`seq_aipe_c_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c_pilot.md) | [`seq_aipe_c`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md) | [`seq_aipe_c_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md) |
| [`seq_aipe_r_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r_pilot.md) | [`seq_aipe_r`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md) | [`seq_aipe_r_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md) |
| [`seq_aipe_R2_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_pilot.md) | [`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md) | [`seq_aipe_R2_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_sensitivity.md) |
| [`seq_aipe_gini_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini_pilot.md) | [`seq_aipe_gini`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md) | [`seq_aipe_gini_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md) |

The effect sizes are the standardized mean difference (`smd`), the
coefficient of variation (`cv`), the simple linear regression slope
(`reg_coef`), the standardized mean (`sm`), and a contrast of group
means (`c`), from Kelley, Darku, and Chattopadhyay (2018); the Pearson,
Kendall, and Spearman correlations (`r`) and the squared multiple
correlation (`R2`), from Kelley, Darku, and Chattopadhyay (2019); and
the Gini index (`gini`), from Chattopadhyay and De (2016). The squared
multiple correlation procedure assumes multivariate normality; the
others assume no particular distribution. The point estimates
[`gini_index`](https://yelleknek.github.io/SMSD/reference/gini_index.md),
[`gini_mean_diff`](https://yelleknek.github.io/SMSD/reference/gini_mean_diff.md),
and
[`slr_slope`](https://yelleknek.github.io/SMSD/reference/slr_slope.md)
are exported for use on their own.

## Minimum Risk Point Estimation

Sampling stops when the estimated risk, the cost of estimation error
plus the cost of the observations, is as small as it can be made.

|  |  |  |
|----|----|----|
| **Pilot sample size** | **Check** | **Simulation** |
| [`seq_mr_cv_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv_pilot.md) | [`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md) | [`seq_mr_cv_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md) |
| [`seq_mr_smd_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd_pilot.md) | [`seq_mr_smd`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md) | [`seq_mr_smd_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md) |

The coefficient of variation follows Chattopadhyay and Kelley (2016).
The standardized mean difference follows Chattopadhyay and Kelley
(2017), with one sampling cost or a different cost in each group; that
article assumes normal populations with a common variance.

## A Budget Across Groups

A contrast of several group means is estimated as precisely as a fixed
budget allows when the groups differ in what an observation costs and in
their unknown variances (Chattopadhyay, Bandyopadhyay, Kelley, &
Padalunkal, 2025).

|  |  |
|----|----|
| [`seq_budget_optimal`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md) | The allocation if the standard deviations were known |
| [`seq_budget_pilot`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md) | The pilot sample size |
| [`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md) | The check on the data collected so far |
| [`seq_budget_next`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md) | The group that should receive the next observation |
| [`seq_budget_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md) | A simulation of the whole procedure |
| [`decide_seq_budget_noninferiority`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md) | The noninferiority test |
| [`decide_seq_budget_superiority`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md) | The superiority test |
| [`decide_seq_budget_inferiority`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md) | The inferiority test |
| [`decide_seq_budget_equivalence`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md) | The equivalence test |

The four tests are one-sided and are applied once sampling has ended.

## Classical Procedures for a Mean, a Median, or a Location

These procedures are from Mukhopadhyay and de Silva (2009): confidence
intervals of fixed width or bounded length, and point estimates with
bounded or minimum risk. The purely sequential ones:

|  |  |  |
|----|----|----|
| **Procedure** | **Pilot sample size** | **Check** |
| Fixed width, normal mean | [`seq_fwd_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_pilot.md) | [`seq_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md) |
| Fixed width, any distribution | [`seq_fwd_dfr_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean_pilot.md) | [`seq_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md) |
| Fixed width, negative exponential | [`seq_fwd_nexp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean_pilot.md) | [`seq_fwd_nexp_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean.md) |
| Bounded length, median | [`seq_blci_dfr_median_pilot`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median_pilot.md) | [`seq_blci_dfr_median`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md) |
| Bounded risk, normal mean | [`seq_br_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean_pilot.md) | [`seq_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md) |
| Minimum risk, normal mean | [`seq_mr_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean_pilot.md) | [`seq_mr_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean.md) |
| Minimum risk, exponential mean | [`seq_mr_exp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean_pilot.md) | [`seq_mr_exp_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean.md) |
| Minimum risk, negative exponential | [`seq_mr_nexp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean_pilot.md) | [`seq_mr_nexp_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean.md) |
| Minimum risk, any distribution | [`seq_mr_dfr_mean_arb_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb_pilot.md) | [`seq_mr_dfr_mean_arb`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md) |

The two stage ones, in which one pilot sample fixes the total sample
size:

|  |  |  |
|----|----|----|
| **Procedure** | **Pilot sample size** | **Check** |
| Fixed width, normal mean (Stein) | [`ts_fwd_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean_pilot.md) | [`ts_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md) |
| Fixed width, any distribution | [`ts_fwd_dfr_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean_pilot.md) | [`ts_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md) |
| Fixed width, negative exponential | [`ts_fwd_nexp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean_pilot.md) | [`ts_fwd_nexp_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md) |
| Bounded risk, normal mean | [`ts_br_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean_pilot.md) | [`ts_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md) |
| Bounded risk, exponential mean | [`ts_br_exp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean_pilot.md) | [`ts_br_exp_mean`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md) |

[`seq_fwd_norm_mean_optimal`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_optimal.md)
gives the sample size that a known variance would give for the fixed
width interval for a normal mean, and
[`power_z`](https://yelleknek.github.io/SMSD/reference/power_z.md) gives
the Type II error rate and the sample size of a one sample, one-sided
*z* test. The negative exponential procedures estimate the location of
that distribution.

## Studying a Design Before Collecting Data

The `*_sensitivity()` functions linked above run a whole procedure, from
its pilot sample to its stopping point, many times on simulated data
with a known population value. They report the final sample size to
expect and how much it varies, and, for the interval procedures, the
width of the final interval and how often it contains the population
value. The `batch` argument of the `seq_aipe_*` and `seq_mr_*`
sensitivity functions shows what it costs to check the rule only after
every several observations. Each help page lists the rows of its
summary, which prints through
[`smsd_sensitivity_methods`](https://yelleknek.github.io/SMSD/reference/smsd_sensitivity_methods.md).

## Reading and Reporting a Result

Every check object has methods for
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`format()`](https://rdrr.io/r/base/format.html), and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), and for
the `tidy()` and `glance()` generics of the generics package; see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md).
SMSD does not export those two generics, so write `generics::tidy(out)`
and `generics::glance(out)`, or attach generics or broom first.
[`seq_path`](https://yelleknek.github.io/SMSD/reference/seq_path.md)
applies a check function to data already in hand at each successive
sample size and returns the whole stopping path;
[`plot.smsd_path`](https://yelleknek.github.io/SMSD/reference/plot.smsd_path.md)
draws it, and
[`smsd_path_methods`](https://yelleknek.github.io/SMSD/reference/smsd_path_methods.md)
documents its other methods.
[`seq_report`](https://yelleknek.github.io/SMSD/reference/seq_report.md)
writes the paragraph, with its reference, that reports a check or a path
in a methods or results section.

## What the Procedures Guarantee

The `seq_aipe_*` procedures deliver the requested width by construction:
whenever the stopping rule is met, the interval is no wider than
`width`. Coverage is a different matter. It is an asymptotic property,
which approaches the confidence level as the target narrows and the
final sample grows, and it can fall well short of the confidence level
when the target is wide, when the relationship studied is strong, or
when a squared multiple correlation has many predictors. In 2,000
simulated studies each at the 95% level, with normal data unless noted,
the final interval contained the population value in 0.94 of the studies
for a standardized mean difference of 0.5 at a width of 0.8, in 0.78 for
a Pearson correlation of 0.4 and in 0.70 for one of 0.8 at a width of
0.4, in 0.84 for a regression slope of 0.5 at a width of 0.4, in 0.79
for a coefficient of variation of 0.5 at a width of 0.2 with lognormal
data, and in 0.36 for a squared multiple correlation of 0.5 with 15
predictors at a width of 0.3 (0.85 with 5 predictors). These figures
have Monte Carlo standard errors between 0.005 and 0.011, so the figure
that another help page reports for the same design, from another run,
can differ in the second decimal. The help page of each check function
states what its guarantee is, exact or asymptotic; the sensitivity
functions measure it for the design at hand.

## An Unusable Variance Estimate

In a small sample the estimate of the variance that a stopping rule
needs can be zero, negative, or not finite. A check function then
reports the stopping rule as not met and attaches a note, and its
attribute `variance_unusable` is `TRUE`, so sampling continues; see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md).
This is a deliberate departure from the stopping rules as printed in the
source articles, under which a lower bound would stand in for such an
estimate and the rule would be met at once with an interval far narrower
than requested. A positive estimate, however small, is used as the
articles prescribe. (In
[`ts_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md),
[`ts_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md),
[`ts_fwd_nexp_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md),
and
[`ts_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md),
a first stage with no variation is an error, because there the first
stage variance or scale alone fixes the total sample size.
[`ts_br_exp_mean`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md)
uses the first stage mean and is not affected.)

## Relation to DMAR

SMSD depends on the package DMAR, which is attached with it and provides
the fixed sample size counterparts of the confidence interval
procedures:
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.html),
[`ss_aipe_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.html),
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.html),
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.html),
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.html),
[`ss_aipe_sm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sm.html),
and
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.html)
(there is none for the Gini index). If the population values are known
well enough to fix the sample size in advance, use DMAR; if they are
not, use SMSD.

## References

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods*,
*30*(2), 425–439.
[doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
pre-specified error bound. *Econometrics*, *4*(3), 30.
[doi:10.3390/econometrics4030030](https://doi.org/10.3390/econometrics4030030)

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
difference with minimum risk: Maximizing accuracy and minimizing cost
with sequential estimation. *Psychological Methods*, *22*(1), 94–113.
[doi:10.1037/met0000089](https://doi.org/10.1037/met0000089)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
accuracy in parameter estimation for population correlation
coefficients. *Psychological Methods*, *24*(4), 492–515.
[doi:10.1037/met0000203](https://doi.org/10.1037/met0000203)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology*, *59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

Useful links:

- <https://kenkelley.org>

- <https://github.com/yelleKneK/SMSD>

- Report bugs at <https://github.com/yelleKneK/SMSD/issues>

## Author

**Maintainer**: Ken Kelley <kkelley@nd.edu>
([ORCID](https://orcid.org/0000-0002-4756-8360))

Authors:

- Ken Kelley <kkelley@nd.edu>
  ([ORCID](https://orcid.org/0000-0002-4756-8360))

- Francis Bilson Darku <fbilsond@nd.edu>

- Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>
