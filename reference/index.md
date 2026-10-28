# Package index

## A Narrow Confidence Interval for an Effect Size

Sequential accuracy in parameter estimation (AIPE). Each effect size has
a function for the pilot sample size and a function that takes the data
collected so far and reports the interval, the stopping criterion, and
whether sampling can stop.

- [`seq_aipe_R2()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md)
  : Sequential AIPE Check for the Squared Multiple Correlation
  Coefficient
- [`seq_aipe_R2_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_pilot.md)
  : Pilot Sample Size for Sequential AIPE: Squared Multiple Correlation
  Coefficient
- [`seq_aipe_c()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md)
  : Sequential AIPE Check for a Contrast
- [`seq_aipe_c_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c_pilot.md)
  : Pilot Sample Size for Sequential AIPE: Contrast
- [`seq_aipe_cv()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md)
  : Sequential AIPE Check for the Coefficient of Variation
- [`seq_aipe_cv_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv_pilot.md)
  : Pilot Sample Size for Sequential AIPE: Coefficient of Variation
- [`seq_aipe_gini()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md)
  : Sequential AIPE Check for the Gini Index
- [`seq_aipe_gini_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini_pilot.md)
  : Pilot Sample Size for Sequential AIPE: Gini Index
- [`seq_aipe_r()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md)
  : Sequential AIPE Check for a Correlation Coefficient
- [`seq_aipe_r_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r_pilot.md)
  : Pilot Sample Size for Sequential AIPE: Correlation Coefficient
- [`seq_aipe_reg_coef()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md)
  : Sequential AIPE Check for the Simple Linear Regression Slope
- [`seq_aipe_reg_coef_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef_pilot.md)
  : Pilot Sample Size for Sequential AIPE: Simple Linear Regression
  Slope
- [`seq_aipe_sm()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm.md)
  : Sequential AIPE Check for the Standardized Mean
- [`seq_aipe_sm_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm_pilot.md)
  : Pilot Sample Size for Sequential AIPE: Standardized Mean
- [`seq_aipe_smd()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md)
  : Sequential AIPE Check for the Standardized Mean Difference
- [`seq_aipe_smd_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd_pilot.md)
  : Pilot Sample Size for Sequential AIPE: Standardized Mean Difference

### The Effect Size Estimates

Point estimates used by the procedures above.

- [`gini_index()`](https://yelleknek.github.io/SMSD/reference/gini_index.md)
  : Gini Index
- [`gini_mean_diff()`](https://yelleknek.github.io/SMSD/reference/gini_mean_diff.md)
  : Gini's Mean Difference
- [`slr_slope()`](https://yelleknek.github.io/SMSD/reference/slr_slope.md)
  : Simple Linear Regression Slope

## Minimum Risk Point Estimation

Sequential procedures that stop when the estimated risk, the cost of
estimation error plus the cost of sampling, is as small as it can be
made: the coefficient of variation and the standardized mean difference.
The minimum risk procedures for a mean are with the classical
procedures.

- [`seq_mr_cv()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
  : Minimum Risk Point Estimation of the Coefficient of Variation
- [`seq_mr_cv_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv_pilot.md)
  : Pilot Sample Size for Minimum Risk Point Estimation of the
  Coefficient of Variation
- [`seq_mr_smd()`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)
  : Minimum Risk Point Estimation of the Standardized Mean Difference
- [`seq_mr_smd_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd_pilot.md)
  : Pilot Sample Sizes for Minimum Risk Point Estimation of the
  Standardized Mean Difference

## A Budget Across Groups

A sequential procedure for a linear contrast of several group means
under a fixed budget, with unequal sampling costs and unequal variances,
and the decisions it supports.

- [`decide_seq_budget_noninferiority()`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md)
  [`decide_seq_budget_superiority()`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md)
  [`decide_seq_budget_inferiority()`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md)
  [`decide_seq_budget_equivalence()`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md)
  : Noninferiority and Equivalence Decisions for a Budget Constrained
  Contrast
- [`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
  : Sequential Check for a Budget Constrained K-Population Contrast
- [`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
  : Which Group to Sample Next Under a Budget Constraint
- [`seq_budget_optimal()`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md)
  : Optimal Sample Sizes for a Budget Constrained Contrast
- [`seq_budget_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md)
  : Pilot Sample Size for the Budget Constrained Sequential Procedure

## Classical Procedures for a Mean, a Median, or a Location

Fixed width, bounded length, bounded risk, minimum risk, and two stage
procedures from Mukhopadhyay and de Silva (2009).

### Fixed Width and Bounded Length Confidence Intervals

- [`seq_blci_dfr_median()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md)
  : Bounded Length Confidence Interval for a Median, Distribution Free
- [`seq_blci_dfr_median_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median_pilot.md)
  : Pilot Sample Size for a Bounded Length Confidence Interval for a
  Median, Distribution Free
- [`seq_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md)
  : Fixed Width Confidence Interval for a Mean, Distribution Free
- [`seq_fwd_dfr_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean_pilot.md)
  : Pilot Sample Size for a Fixed Width Confidence Interval for a Mean,
  Distribution Free
- [`seq_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean.md)
  : Fixed Width Confidence Interval for the Location of a Negative
  Exponential Distribution
- [`seq_fwd_nexp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean_pilot.md)
  : Pilot Sample Size for a Fixed Width Confidence Interval for the
  Location of a Negative Exponential Distribution
- [`seq_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md)
  : Fixed Width Confidence Interval for a Normal Mean
- [`seq_fwd_norm_mean_optimal()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_optimal.md)
  : Known Variance Sample Size for a Fixed Width Interval for a Normal
  Mean
- [`seq_fwd_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_pilot.md)
  : Pilot Sample Size for a Fixed Width Confidence Interval for a Normal
  Mean

### Bounded Risk Point Estimation

- [`seq_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md)
  : Bounded Risk Point Estimation of a Normal Mean
- [`seq_br_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean_pilot.md)
  : Pilot Sample Size for Bounded Risk Point Estimation of a Normal Mean

### Minimum Risk Point Estimation of a Mean or a Location

- [`seq_mr_dfr_mean_arb()`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md)
  : Minimum Risk Point Estimation of a Mean, Distribution Free
- [`seq_mr_dfr_mean_arb_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb_pilot.md)
  : Pilot Sample Size for Minimum Risk Point Estimation of a Mean,
  Distribution Free
- [`seq_mr_exp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean.md)
  : Minimum Risk Point Estimation of an Exponential Mean
- [`seq_mr_exp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean_pilot.md)
  : Pilot Sample Size for Minimum Risk Point Estimation of an
  Exponential Mean
- [`seq_mr_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean.md)
  : Minimum Risk Point Estimation of the Location of a Negative
  Exponential Distribution
- [`seq_mr_nexp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean_pilot.md)
  : Pilot Sample Size for Minimum Risk Point Estimation of the Location
  of a Negative Exponential Distribution
- [`seq_mr_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean.md)
  : Minimum Risk Point Estimation of a Normal Mean
- [`seq_mr_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean_pilot.md)
  : Pilot Sample Size for Minimum Risk Point Estimation of a Normal Mean

### Two Stage Procedures

- [`ts_br_exp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md)
  : Two Stage Bounded Risk Point Estimation of an Exponential Mean
- [`ts_br_exp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean_pilot.md)
  : Pilot Sample Size for Two Stage Bounded Risk Estimation of an
  Exponential Mean
- [`ts_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md)
  : Two Stage Bounded Risk Point Estimation of a Normal Mean
- [`ts_br_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean_pilot.md)
  : Pilot Sample Size for Two Stage Bounded Risk Estimation of a Normal
  Mean
- [`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md)
  : Two Stage Fixed Width Confidence Interval for a Distribution Free
  Mean
- [`ts_fwd_dfr_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean_pilot.md)
  : Pilot Sample Size for a Two Stage Fixed Width Interval for a
  Distribution Free Mean
- [`ts_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md)
  : Two Stage Fixed Width Interval for a Negative Exponential Location
- [`ts_fwd_nexp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean_pilot.md)
  : Pilot Sample Size for a Two Stage Fixed Width Interval for a
  Negative Exponential Location
- [`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md)
  : Two Stage Fixed Width Confidence Interval for a Normal Mean
- [`ts_fwd_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean_pilot.md)
  : Pilot Sample Size for a Two Stage Fixed Width Interval for a Normal
  Mean

### Supporting Calculations

- [`power_z()`](https://yelleknek.github.io/SMSD/reference/power_z.md) :
  Type II Error and Sample Size for a One Sample, One Sided z Test

## Studying a Design Before Collecting Data

Monte Carlo runs of a whole sequential procedure against a known data
generating mechanism: the final sample size to expect and, for the
interval procedures, the interval width and the coverage.

- [`seq_aipe_R2_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_sensitivity.md)
  : Monte Carlo Sensitivity Analysis for Sequential AIPE: Squared
  Multiple Correlation Coefficient
- [`seq_aipe_smd_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
  [`seq_aipe_r_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
  [`seq_aipe_cv_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
  [`seq_aipe_reg_coef_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
  [`seq_aipe_sm_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
  [`seq_aipe_gini_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
  [`seq_aipe_c_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
  : Monte Carlo Sensitivity Analysis for Sequential AIPE Procedures
- [`seq_budget_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)
  : Monte Carlo Sensitivity Analysis for the Budget Constrained
  Procedure
- [`seq_mr_cv_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)
  [`seq_mr_smd_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)
  : Monte Carlo Sensitivity Analysis for Minimum Risk Point Estimation
- [`format(`*`<smsd_sensitivity>`*`)`](https://yelleknek.github.io/SMSD/reference/smsd_sensitivity_methods.md)
  [`print(`*`<smsd_sensitivity>`*`)`](https://yelleknek.github.io/SMSD/reference/smsd_sensitivity_methods.md)
  : Print a Sensitivity Summary

## Reading and Reporting a Result

The object every check function returns, with its print, summary, tidy,
and glance methods; the path a procedure takes through the stages of the
data in hand, with its plot; and a report of a procedure in words.

- [`plot(`*`<smsd_path>`*`)`](https://yelleknek.github.io/SMSD/reference/plot.smsd_path.md)
  : Plot the Stopping Path of a Sequential Procedure
- [`seq_report()`](https://yelleknek.github.io/SMSD/reference/seq_report.md)
  [`format(`*`<smsd_report>`*`)`](https://yelleknek.github.io/SMSD/reference/seq_report.md)
  [`print(`*`<smsd_report>`*`)`](https://yelleknek.github.io/SMSD/reference/seq_report.md)
  : Report a Sequential Procedure in Words
- [`format(`*`<smsd_check>`*`)`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
  [`print(`*`<smsd_check>`*`)`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
  [`summary(`*`<smsd_check>`*`)`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
  [`as.data.frame(`*`<smsd_check>`*`)`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
  [`tidy(`*`<smsd_check>`*`)`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
  [`glance(`*`<smsd_check>`*`)`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
  : Methods for Sequential Check Objects
- [`format(`*`<smsd_path>`*`)`](https://yelleknek.github.io/SMSD/reference/smsd_path_methods.md)
  [`print(`*`<smsd_path>`*`)`](https://yelleknek.github.io/SMSD/reference/smsd_path_methods.md)
  [`summary(`*`<smsd_path>`*`)`](https://yelleknek.github.io/SMSD/reference/smsd_path_methods.md)
  [`as.data.frame(`*`<smsd_path>`*`)`](https://yelleknek.github.io/SMSD/reference/smsd_path_methods.md)
  : Methods for Sequential Path Objects
- [`seq_path()`](https://yelleknek.github.io/SMSD/reference/seq_path.md)
  : Stopping Path of a Sequential Procedure for Data in Hand

## The Package

- [`SMSD`](https://yelleknek.github.io/SMSD/reference/SMSD-package.md)
  [`SMSD-package`](https://yelleknek.github.io/SMSD/reference/SMSD-package.md)
  : SMSD: Sequential Methods for Study Design
