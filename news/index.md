# Changelog

## SMSD 1.0.0

### Release Preparation

- Added two installed vignettes: replications of worked examples and new
  simulation rows from published sequential methods papers; and a study
  of coverage and sample size in difficult Pearson and multiple
  R-squared designs, with individual-study data, uncertainty intervals,
  and figures. Simulation caches record seeds and implementation
  fingerprints.
- AIPE check displays and generated reports identify their confidence
  level as asymptotic, separately from satisfying the interval width
  criterion. The multiple R-squared help page and coverage vignette
  distinguish interior population values from the null boundary, where
  the normal approximation does not acquire nominal coverage with
  increasing n. A reproducible normal model reference calculation
  compares alternative intervals with fixed sample sizes that meet the
  coverage and width targets.
- Fixed overflow in the standardized mean difference used by
  [`power_z()`](https://yelleknek.github.io/SMSD/reference/power_z.md)
  when opposite signed means have a finite difference in standard
  deviation units. Two-sided normal and t quantiles retain finite values
  at extreme valid confidence levels.
- Rank correlation variance calculations use sorting and cumulative
  counts, with O(n log n) work and O(n) working memory, instead of
  quadratic comparisons or comparison matrices. Tied ranks retain the
  same formulas; the existing statistical warning about ties remains.
- Added a site builder that renders a committed checkout, excluding
  ignored working notes. Continuous integration can run on preparation
  branches and on demand, and includes Windows R-devel.

First CRAN release. SMSD (“Sequential Methods for Study Design”) is for
planning a study when the population values that a sample size formula
needs (a standard deviation, a correlation, the effect size itself) are
unknown. Data are collected in stages, and after each stage a stopping
rule computed from the data in hand says whether sampling can stop. SMSD
is the sequential companion of ‘DMAR’, which plans fixed sample sizes
and is installed with it.

Every procedure has a `*_pilot()` function, which returns the pilot
sample size, and a check function, which takes the data collected so far
and evaluates the stopping rule once. The accuracy in parameter
estimation procedures, the minimum risk procedures for the coefficient
of variation and the standardized mean difference, and the budget
procedure also have a `*_sensitivity()` function, which simulates the
whole procedure for a planned design. The user owns the sampling loop: a
check function does not prompt, print, or keep anything between calls.

### The Procedures

#### A Narrow Confidence Interval for an Effect Size

- **Sequential accuracy in parameter estimation (AIPE).** Sampling stops
  at the first sample size that meets the stopping rule, and the
  confidence interval at that point is no wider than the requested full
  width. The standardized mean difference
  ([`seq_aipe_smd()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md)),
  the coefficient of variation
  ([`seq_aipe_cv()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md)),
  the simple linear regression slope
  ([`seq_aipe_reg_coef()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md)),
  the standardized mean
  ([`seq_aipe_sm()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm.md)),
  and a contrast of group means
  ([`seq_aipe_c()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md))
  follow Kelley, Darku, and Chattopadhyay (2018, *Psychological Methods,
  23*(2), 226–243, <doi:10.1037/met0000127>). The Pearson, Kendall, and
  Spearman correlations
  ([`seq_aipe_r()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md))
  follow Kelley, Darku, and Chattopadhyay (2019, *Psychological Methods,
  24*(4), 492–515, <doi:10.1037/met0000203>). The Gini index
  ([`seq_aipe_gini()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md))
  follows Chattopadhyay and De (2016, *Econometrics, 4*(3), 30,
  <doi:10.3390/econometrics4030030>). None of these assumes a particular
  distribution for the data.

- **The squared multiple correlation coefficient.**
  [`seq_aipe_R2_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_pilot.md),
  [`seq_aipe_R2()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md),
  and
  [`seq_aipe_R2_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_sensitivity.md)
  implement the procedure of Kelley, Darku, and Chattopadhyay (2019,
  pp. 503–507), with the Wald interval (the default) or the log
  transformed interval of Bonett and Wright (2011). Unlike the other
  AIPE procedures, this one assumes multivariate normality. The
  procedure reproduces the article’s Tables 11 to 14. For the Wald
  interval the stopping rule uses `(R2 + 1/n) * (1 - R2 + 1/n)^2`, with
  1/n added to both factors of the variance. Equation 63 as printed has
  `(1 - (R2 + 1/n))^2` as its second factor, which does not reproduce
  the published tables and can leave the interval wider than the target.
  The help page documents the difference.

- **What is guaranteed, and what is not.** Whenever a `seq_aipe_*`
  stopping rule is met, the interval is no wider than the target.
  Coverage is an asymptotic property under the relevant regularity
  conditions: it approaches the confidence level as the target narrows
  and the final sample grows, and it can fall well short in small
  samples. Each help page says so and reports what was measured for that
  effect size. Coverage at wide targets is below the confidence level
  for the Pearson correlation, the regression slope, and the coefficient
  of variation, and for any correlation when the association is strong;
  for the squared multiple correlation it falls badly when the sample is
  small relative to the number of predictors, and
  [`seq_aipe_R2()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md)
  attaches a note whenever it reports an interval from fewer than
  `p * max(20, 2 * p)` observations. A narrower target or a larger
  minimum pilot sample (`m_0`) can help, and the sensitivity functions
  show what a planned design can expect. At the null multiple R-squared
  boundary, increasing n does not restore nominal Wald coverage; a
  different interval construction is needed.

- **The correction to Equation 26.** The asymptotic variance of the
  standardized mean difference restores a factor of 1/2 that Equation 26
  of Kelley, Darku, and Chattopadhyay (2018) omits. The two forms differ
  only when the means differ and the third central moments of the groups
  differ. The authors’ correction is installed with the package
  (`system.file("eq26_correction.pdf", package = "SMSD")`), and
  [`?seq_aipe_smd`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md)
  gives its web address.

- [`seq_aipe_r()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md)
  takes its paired data as `x` and `y`, like
  [`seq_aipe_reg_coef()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md).
  With `method = "kendall"` or `"spearman"` it warns when the data
  contain ties, because those variance estimators assume continuous
  data.

- [`seq_aipe_c()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md)
  requires balanced groups, as in Example 6 of Kelley, Darku, and
  Chattopadhyay (2018); `criterion_n` and `current_n` are per group. The
  unbalanced extension of the article (Equations 39 to 47) is not
  implemented.

- [`seq_aipe_reg_coef()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md)
  and
  [`seq_aipe_c()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md)
  estimate quantities that have units, and their help pages have a Units
  paragraph: when the variables are recorded on a scale that makes the
  variance estimate tiny, a note says so and advises rescaling.

- The point estimates
  [`gini_index()`](https://yelleknek.github.io/SMSD/reference/gini_index.md),
  [`gini_mean_diff()`](https://yelleknek.github.io/SMSD/reference/gini_mean_diff.md),
  and
  [`slr_slope()`](https://yelleknek.github.io/SMSD/reference/slr_slope.md)
  are exported for use on their own.

#### Minimum Risk Point Estimation

- **Minimum risk point estimation (MRPE).** Sampling stops when the
  estimated risk, the cost of estimation error plus the cost of the
  observations, is as small as it can be made:
  [`seq_mr_cv()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
  for the coefficient of variation (Chattopadhyay & Kelley, 2016,
  *Multivariate Behavioral Research, 51*(5), 627–648,
  <doi:10.1080/00273171.2016.1203279>) and
  [`seq_mr_smd()`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)
  for the standardized mean difference (Chattopadhyay & Kelley, 2017,
  *Psychological Methods, 22*(1), 94–113, <doi:10.1037/met0000089>),
  each with its pilot function.

- [`seq_mr_smd()`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)
  takes one sampling cost `c` for the equal cost procedure (Equation 15
  of the article) or `c = c(c1, c2)` for the article’s unequal cost
  procedure (“A General Scenario”, Equations 17 to 21). The rows
  `next_n1` and `next_n2` give the size each group should have at the
  next check, so the sampling loop is the same whichever group costs
  more. The article assumes normal populations with a common variance,
  and
  [`seq_mr_smd_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)
  shows what a skewed population costs.
  [`?seq_mr_smd`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)
  lists each point at which the function departs from the article as
  printed.

- [`seq_mr_cv()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
  has `variance = c("corrected", "published")`. The default is the
  estimator of Kelley, Darku, and Chattopadhyay (2018, Equations 30 to
  32), which is consistent for any population with a finite fourth
  moment. `"published"` is Equation 22 of Chattopadhyay and Kelley
  (2016), offered so that the article’s worked example and tables, and
  [`DMAR::mr_cv()`](https://yelleknek.github.io/DMAR/reference/mr_cv.html),
  can be reproduced. The two agree for normal data, and for any
  population with no excess kurtosis, and can differ otherwise. The
  derivation is installed with the package
  (`system.file("cv_variance_note.md", package = "SMSD")`).

- [`seq_mr_cv_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)
  and
  [`seq_mr_smd_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)
  run a minimum risk procedure to its stopping point on simulated data
  and report the best fixed sample size and its risk, the distribution
  of the final sample size, the risk efficiency with its Monte Carlo
  standard error, and the regret.

#### A Budget Across Groups

- **A contrast under a fixed budget** (Chattopadhyay, Bandyopadhyay,
  Kelley, & Padalunkal, 2025, *Psychological Methods, 30*(2), 425–439,
  <doi:10.1037/met0000570>).
  [`seq_budget_optimal()`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md)
  returns the optimal sample sizes of Theorem 1;
  [`seq_budget_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md)
  returns the pilot size;
  [`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
  evaluates the Equation 18 criterion on the data in hand;
  [`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
  names the group that should receive the next observation; and
  [`seq_budget_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)
  simulates the whole procedure for a planned design.

- [`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
  never names a group the remaining budget cannot pay for, so a study
  that follows it one observation at a time cannot overspend. It returns
  `NA` when sampling has ended, and a sampling loop should end on that
  `NA`.

- [`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
  reports how sampling ended in three attributes: `rule_met` (every
  group reached its criterion), `budget_exhausted` (the budget cannot
  pay for another needed observation), and `is_satisfied`. Data that
  already cost more than the budget are flagged (`budget_exceeded`). By
  default every group is judged on the current data at each check; the
  `closed` argument, with `order = "rounds"` in
  [`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md),
  gives the package’s reading of the article’s first passage rule, and
  [`?seq_budget_next`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
  lists the three details in which it differs from the article.

- **Four decisions, one family.**
  [`decide_seq_budget_noninferiority()`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md),
  [`decide_seq_budget_superiority()`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md),
  [`decide_seq_budget_inferiority()`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md),
  and
  [`decide_seq_budget_equivalence()`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md)
  are one-sided tests at `alpha_level` (default 0.05; two one-sided
  tests for equivalence), computed from the contrast estimate and its
  standard error as in the article and in
  [`DMAR::ss_seq_c()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.html).
  None reads the two-sided interval stored on the check. The help states
  each null and alternative hypothesis in symbols and in words and the
  sign convention of each margin.

- The help for the budget functions states how far coverage and the
  level of the tests fall short of nominal, by pilot size.

#### Classical Procedures for a Mean, a Median, or a Location

- Procedures from Mukhopadhyay and de Silva (2009, *Sequential methods
  and their applications*, <ISBN:9780367386535>): fixed width confidence
  intervals
  ([`seq_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md),
  [`seq_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md),
  [`seq_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean.md)),
  a bounded length interval for a median
  ([`seq_blci_dfr_median()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md)),
  bounded risk point estimation
  ([`seq_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md)),
  minimum risk point estimation
  ([`seq_mr_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean.md),
  [`seq_mr_exp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean.md),
  [`seq_mr_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean.md),
  [`seq_mr_dfr_mean_arb()`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md)),
  and the two stage procedures
  ([`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md),
  [`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md),
  [`ts_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md),
  [`ts_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md),
  [`ts_br_exp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md)),
  each with its pilot function.

- Each help page says what is estimated (the `*_nexp_mean` functions
  estimate the location of a negative exponential distribution, not its
  mean), what `d` is in that procedure (a half width, a full width, or
  half the largest length allowed), every returned row and attribute,
  and whether the guarantee is exact or asymptotic, with the coverage or
  risk measured by simulation where it is asymptotic. Stein’s two stage
  interval for a normal mean
  ([`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md))
  is exact; the distribution free two stage interval
  ([`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md))
  is asymptotic and falls short of its nominal level, and its help page
  gives the size of the shortfall.

- [`seq_fwd_norm_mean_optimal()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_optimal.md)
  gives the sample size for a fixed width interval for a normal mean
  when the variance is known, the benchmark the sequential procedures
  are compared with.
  [`power_z()`](https://yelleknek.github.io/SMSD/reference/power_z.md)
  gives the Type II error rate and the sample size of a one sample,
  one-sided *z* test.

#### Studying a Design Before Collecting Data

- Every `seq_aipe_*` procedure, the contrast and the squared multiple
  correlation included, both minimum risk procedures, and the budget
  procedure have a `*_sensitivity()` function. Each runs the procedure
  from its pilot sample to its stopping point `G` times against a known
  data generating mechanism and returns a `term`/`value` summary of
  class `smsd_sensitivity`, which prints in fixed notation.

- The `seq_aipe_*` and `seq_mr_*` sensitivity functions have
  `batch = 1L`, the number of observations added between checks of the
  stopping rule. The default is the purely sequential procedure of the
  articles; a larger value shows the overshoot to expect when a study
  can check only after every several observations. `max_n` is an exact
  cap for any batch size.

- The AIPE summary has 13 rows, among them `pct_hit_max_n`, the
  percentage of replications that the cap stopped before the rule was
  met, and `pct_variance_unusable`, the percentage in which the variance
  estimate was unusable at one or more checks.

- `seed = NULL` is the default everywhere. With a number the result is
  reproducible and the caller’s random number stream is left as it was
  found.

- A true value at which the simulated data could not be used (for
  example a standardized mean difference of 1e8 or more, or a
  coefficient of variation of 1e-8 or less) is an error that names the
  argument, raised before anything is simulated.

#### Reading and Reporting a Result

- **The whole stopping path in one call.**
  `seq_path(check, data, ..., start = NULL, by = 1L)` applies a
  sequential check function to data already in hand at each successive
  sample size and returns an `smsd_path` object: one row per stage and
  the first sample size at which the stopping rule is met (`stop_n`). It
  replaces the loop a user would otherwise write by hand, starts at the
  procedure’s pilot sample size by default, and has
  [`print()`](https://rdrr.io/r/base/print.html),
  [`format()`](https://rdrr.io/r/base/format.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods. The
  plot shows the estimate with its interval, and the current sample size
  against the criterion, at every stage.

- The path shows what the rule would have said at each sample size. A
  sequential design reports the estimate at the first sample size at
  which the rule is met; the later stages are shown for illustration
  only, and the help page and the print method say so. All of the data
  go in `data`, and
  [`seq_path()`](https://yelleknek.github.io/SMSD/reference/seq_path.md)
  stops when a second group or variable is passed outside it, where it
  would be seen in full from the first stage on. Data with missing
  values are refused. The two stage and budget procedures, which have no
  path indexed by one common sample size, are refused with advice.

- **A paragraph for the write-up.**
  [`seq_report()`](https://yelleknek.github.io/SMSD/reference/seq_report.md)
  turns a check or a path into a paragraph for a methods or results
  section: the procedure and its source, the target, the sample size at
  which the rule was met (or that it has not yet been met), and the
  estimate with its interval. The full reference is returned with it.
  For the budget procedure it says how sampling ended, and a spent
  budget is never reported as a stopping rule that was met.

### What Every Check Function Shares

- **The object.** A check function returns an object of class
  `smsd_check`: a data frame with the columns `term` and `value`. The
  stopping decision is `attr(out, "is_satisfied")`, which is always
  `TRUE` or `FALSE`, never `NA`; the method, the current and criterion
  sample sizes, the targets, and any notes are attributes as well.
  Stored values are never rounded.

- **Its methods.** [`print()`](https://rdrr.io/r/base/print.html),
  [`format()`](https://rdrr.io/r/base/format.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), and
  methods for the `tidy()` and `glance()` generics of the ‘generics’
  package. SMSD does not export those two generics: call
  `generics::tidy(out)` and `generics::glance(out)`, or attach
  ‘generics’ or ‘broom’ first. `tidy()` uses the column layout of
  ‘DMAR’. Subsetting with `[` returns a plain data frame whenever the
  result is no longer the complete check, so a partial object is never
  printed as a complete one.

- **An unusable variance estimate means the rule is not met.** In a
  small sample the variance estimate that a stopping rule needs can be
  zero, negative, or not finite. The check then reports the rule as not
  met, sets `attr(out, "variance_unusable")` to `TRUE`, and attaches a
  note that [`print()`](https://rdrr.io/r/base/print.html) shows, so
  sampling continues. This is a deliberate departure from the stopping
  rules as printed in the source articles. Under them a lower bound of
  n^-3 stands in for such an estimate and the rule is met at once, with
  an interval about a tenth of the requested width. At a target width of
  0.4 that happened at the pilot sample in 14 percent of simulated
  studies of a Pearson correlation of 0.4 and in 16 percent of studies
  of a regression slope of 0.5. A positive estimate, however small, is
  used as the articles prescribe.

- **No variation is judged the same way everywhere.** A sample whose
  values agree to ten significant digits is treated as having no
  variation, so scores that are equal on paper but differ by rounding
  error (gains of 0.3 computed as posttest minus pretest) are not taken
  for data with a tiny variance. In the purely sequential functions that
  estimate a variance or a scale, such a sample leaves the rule not met.
  (A predictor that does not vary is an error in
  [`seq_aipe_R2()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md),
  because the squared multiple correlation is then not defined.) In
  [`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md),
  [`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md),
  [`ts_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md),
  and
  [`ts_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md)
  it is an error, because there the first stage variance or scale alone
  fixes the total sample size. Three functions estimate neither and are
  not affected:
  [`seq_blci_dfr_median()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md),
  where constant data give an interval of length zero and so meet the
  rule, and
  [`seq_mr_exp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean.md)
  and
  [`ts_br_exp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md),
  whose rules use only the sample mean.

- **Results are right at any magnitude of the data.** The `seq_aipe_*`
  variance estimators,
  [`gini_index()`](https://yelleknek.github.io/SMSD/reference/gini_index.md),
  [`gini_mean_diff()`](https://yelleknek.github.io/SMSD/reference/gini_mean_diff.md),
  [`slr_slope()`](https://yelleknek.github.io/SMSD/reference/slr_slope.md),
  [`seq_mr_cv()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md),
  [`seq_mr_smd()`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md),
  and the two stage fixed width rules lose no precision on any scale a
  double can hold. Data recorded at 1e-150 or 1e150 times their size
  give the same unit free effect size and the same stopping decision; an
  estimate that has units is rescaled with the data; and a two stage
  rule gives the same total sample size when its target is stated in the
  same units. (The slope and the contrast have units, and their stopping
  rules depend on the units in which the variables are recorded, as in
  the article; see the Units paragraphs of
  [`?seq_aipe_reg_coef`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md)
  and
  [`?seq_aipe_c`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md).)
  Where a classical procedure cannot represent a variance at an extreme
  scale, it stops with a message that says to change the units. No
  result holds `Inf` or `NaN`: a sample size, a limit, or a risk too
  large to represent is an error whose message names the arguments to
  check.

- **Every argument is checked.** A missing, non-finite, wrongly sized,
  or wrongly typed argument stops with a message that names the argument
  and says what was received. Counts such as `m_0`, `G`, `max_n`, and
  `batch` must be whole numbers.

- **One set of names.** Arguments are in snake_case, the confidence
  level is always `conf_level = 0.95`, and the data arguments are `x`
  for one sample, `group_1` and `group_2` for two independent groups,
  `x` and `y` for paired variables, `groups` for several groups, and
  `pilot_data` and `additional_data` for the two stage procedures. Pilot
  functions take the method’s own arguments, then `m_0`, then `gamma`.

### Corrections to Earlier Development Versions

These changed results relative to the development versions distributed
on GitHub before this release.

- [`seq_mr_smd()`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)
  applies Equation 15 of Chattopadhyay and Kelley (2017) exactly.
  Earlier versions placed the correction term inside the square root and
  stopped early (at 71 per group instead of 75 in the article’s worked
  example); a still earlier one divided the mean difference by the
  pooled variance where the pooled standard deviation belongs. The
  article’s worked example and its Tables 1 to 4 now reproduce. With two
  costs the function implements the article’s own procedure; Equations
  18 and 42 are implemented with the squared term that the article’s
  first order conditions (Equation 41) give, where the article prints
  the first power.

- [`seq_mr_cv()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md):
  earlier versions used the corrected variance estimator while citing
  Equation 22 of Chattopadhyay and Kelley (2016), so results for data
  that are not normal differed from the article and from
  [`DMAR::mr_cv()`](https://yelleknek.github.io/DMAR/reference/mr_cv.html)
  without explanation. The `variance` argument now makes the choice
  explicit.

- The pilot sample sizes of the classical procedures follow one rule,
  `max(m_0, ceiling(K^(1 / (1 + gamma))))`, where `K` is the procedure’s
  criterion constant, with `m_0` and `gamma` as arguments. Several
  earlier pilots equaled the optimal sample size at a unit scale, which
  was several times too large for data with a small standard deviation
  and no protection for data with a large one. The two stage pilots have
  larger default minimums (10 for
  [`ts_fwd_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean_pilot.md)
  and
  [`ts_fwd_nexp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean_pilot.md),
  30 for
  [`ts_fwd_dfr_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean_pilot.md),
  20 for
  [`ts_br_exp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean_pilot.md)),
  because a very small first stage costs many times the optimal sample
  size or, for the distribution free interval, much of the coverage.

- [`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md)
  and
  [`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md)
  use the two-sided quantile for their two-sided intervals, as do the
  `seq_fwd_*` and `seq_blci_*` pilots. Earlier versions used the
  one-sided quantile, which at `conf_level = 0.95` made the total sample
  size about 30 percent too small, and more than that for a small first
  stage in Stein’s procedure.

- [`seq_mr_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean.md)
  uses the scale estimate to the power `k` in its stopping criterion and
  reports the risk at the current sample size; its `location` row is the
  sample minimum. The negative exponential pilots use the criterion’s
  own constant.

- The classical stopping criteria round the criterion sample size up.
  Earlier versions truncated it, which declared the rule met one
  observation early.

- [`seq_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md)
  and
  [`ts_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md)
  report the risk at the current sample size and compute their constants
  in closed form, so that the reported risk no longer sits a shade above
  the bound at the sample size the rule calls for, and
  [`ts_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md)
  no longer fails for a first stage above 344 observations.

- [`seq_blci_dfr_median()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md)
  requires the smallest sample size at which an order statistic interval
  can reach the confidence level (6 observations at 95%). Earlier
  versions could stop at two observations. The confidence level holds at
  a fixed sample size; a check whose rule is met carries a note that the
  coverage of the sequential procedure is lower, and the help page gives
  the measured figures.

- Constant data no longer satisfy the stopping rule of a procedure that
  estimates a variance or a scale. Earlier versions of
  [`seq_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md),
  [`seq_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md),
  [`seq_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean.md),
  [`seq_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md),
  [`seq_mr_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean.md),
  [`seq_mr_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean.md),
  and
  [`seq_mr_dfr_mean_arb()`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md)
  reported an optimal sample size of 0 or 1 for such data and the rule
  as met, and
  [`seq_aipe_c()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md)
  reported a met rule with an interval of zero width.
  [`seq_blci_dfr_median()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md)
  and
  [`seq_mr_exp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean.md)
  estimate no variance, and constant data can meet their rules, as
  before.
  [`seq_mr_exp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean.md)
  now requires positive data.

- [`seq_aipe_c()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md)
  requires balanced groups. The removed unbalanced branch could report a
  satisfied stopping rule with an interval wider than requested.

- The central moment estimators behind the variance estimates are
  computed on centered data. The earlier raw power sum forms lost
  precision when the mean was large relative to the standard deviation.

- [`seq_aipe_r()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md):
  the Pearson procedure is defined when the sample covariance is exactly
  zero, as happens with ratings, and the Kendall and Spearman variance
  estimators are computed in exact arithmetic on the ranks, so that an
  estimate that is zero on paper is exactly zero. Perfectly discordant
  pairs under Kendall’s tau, and two dichotomous variables in exact
  opposition under Spearman’s rho, are therefore reported as not yet
  usable and not as a met rule.

- [`seq_aipe_cv()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md)
  and
  [`seq_mr_cv()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
  report a sample whose mean is zero, negative, or positive only by
  rounding error (centered or standardized scores) as not yet estimable.
  Earlier versions could return a negative or an enormous coefficient of
  variation.

- [`gini_index()`](https://yelleknek.github.io/SMSD/reference/gini_index.md)
  and
  [`gini_mean_diff()`](https://yelleknek.github.io/SMSD/reference/gini_mean_diff.md):
  integer input no longer overflows, a missing value with
  `na.rm = FALSE` is an error and not a wrong number, negative data are
  an error for the Gini index, and the index can no longer exceed 1
  through rounding.

- [`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md):
  only the groups that still need observations count toward
  `budget_exhausted`; with unequal costs a check could previously sit at
  “not satisfied” with no affordable observation left. The superiority,
  inferiority, and equivalence decisions test at a one-sided level of
  0.05 by default, where they previously read the two-sided interval on
  the check (a one-sided level of 0.025 on a 95 percent interval), and
  their margin arguments are `margin`, `margin_lower`, and
  `margin_upper`.

- [`power_z()`](https://yelleknek.github.io/SMSD/reference/power_z.md)
  checks that `mu_1` exceeds `mu_0`, returns the sample size as the
  smallest whole number whose Type II error rate is at most `beta`, and
  handles a very small `beta`.

- The sensitivity functions grow each simulated sample one observation
  at a time unless `batch` says otherwise; `pct_hit_max_n` counts only
  replications that the cap stopped before the rule was met; and
  coverage and the percentage rows are formed from counts.

- Renamed for accuracy: `seq_blci_dfr_mean()` is
  [`seq_blci_dfr_median()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md)
  (it estimates a median); `calculate_c()` is
  `seq_fwd_norm_mean_optimal(variance, d, conf_level)` and returns a
  whole sample size of at least 1; the data argument of the one sample
  classical procedures is `x`;
  [`seq_aipe_r()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md)
  takes `x` and `y`; the pilot functions take `m_0` before `gamma`; and
  the estimate rows of the standardized mean and the slope are named
  `sm` and `reg_coef`.

- Functions of the earliest versions under their present names:
  `seq_aipe_cor()` is
  [`seq_aipe_r()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md),
  `seq_aipe_slr_slope()` is
  [`seq_aipe_reg_coef()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md),
  `seq_aipe_std_mean()` is
  [`seq_aipe_sm()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm.md),
  and `seq_aipe_contrast()` is
  [`seq_aipe_c()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md).
  [`power()`](https://rdrr.io/r/stats/power.html) and the three
  functions it called (`sample_size()`, `type_two_calc()`, and
  `power_seq()`) are
  [`power_z()`](https://yelleknek.github.io/SMSD/reference/power_z.md)
  with its `method` argument. The helper functions those versions
  exported (`V2_contrast()` and the other `V2_*` functions, `mu3()`,
  `mu4()`, `mu13()`, `mu22()`, `mu31()`, `s()`, `chi_prob()`,
  `chi_prob_2()`, and `chi_prob_3()`) are no longer exported.

- [`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md),
  [`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md),
  and
  [`ts_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md)
  no longer take `gamma`, which is an argument of their pilot functions.
  [`ts_br_exp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean_pilot.md)
  takes the loss constant `A` as its first argument:
  `ts_br_exp_mean_pilot(A, w, m_0 = 20, gamma = 1)`.

- The interactive interface of the earliest versions (an `alpha`
  argument, a `pilot = TRUE/FALSE` switch, console messages, and lists
  with title cased names) is gone, with no deprecation shims.
  `seq_ubdes_norm()` was removed;
  [`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
  covers the design it was meant for. The package no longer needs the
  `copula`, `matrixStats`, or `MASS` packages.

### Held Back for a Later Version

- `seq_mr_dfr_mean()`, `ts_mr_dfr_mean()`, and their pilot functions are
  not part of this release. Whether their known lower bound is a bound
  on the standard deviation or on the variance, and how it should set
  the first stage size, could not be confirmed against Mukhopadhyay and
  de Silva (2009).
  [`seq_mr_dfr_mean_arb()`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md),
  which places no restriction on the variance, remains.

- `prob_gn()`, `prob_hn()`, `prob_kn()`, and `seven_1()` are not part of
  this release either. They could not be checked against Mukhopadhyay
  and de Silva (2009), and they return once they have been.

### Documentation

- Vignettes: `introduction`, the overview;
  `step-by-step-sequential-workflow`, one study followed from its pilot
  sample to its stopping point; `complete-tour`, the procedures one by
  one, with a cached simulation study of coverage and sample size;
  `seq-aipe-smd` and `seq-aipe-correlations`; `mrpe-cv` and `mrpe-smd`,
  which reproduce the published worked examples; `seq-budget`;
  `sensitivity-analysis`; and `special-cases-and-simulation`, on the
  classical procedures.

- The package website adds a validation article
  (`vignettes/articles/validation.Rmd` in the source repository), which
  sets the package’s results beside the published tables.

- `citation("SMSD")` gives the package and the six articles whose
  methods it implements, and the same metadata is in `CITATION.cff`.

### Acknowledgment

- Douglas A. Fitts caught an error in the point estimate of an earlier
  version of
  [`seq_aipe_smd()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md).
