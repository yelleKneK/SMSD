# Monte Carlo Sensitivity Analysis for Sequential AIPE Procedures

Each function runs a sequential accuracy in parameter estimation (AIPE)
procedure from its pilot sample to its stopping point `G` times, against
a data generating mechanism whose true effect size is known, and
summarizes the final sample size, the width of the final confidence
interval, and how often that interval contains the true value. Use them
before data collection to see what a planned study should expect: how
large the sample is likely to become, how variable that is, and whether
the confidence interval can be trusted at the planned width.

## Usage

``` r
seq_aipe_smd_sensitivity(
  true_value,
  width,
  conf_level = 0.95,
  distribution = c("normal", "skewed"),
  G = 1000,
  seed = NULL,
  m_0 = 4,
  max_n = 100000L,
  batch = 1L
)

seq_aipe_r_sensitivity(
  true_value,
  width,
  conf_level = 0.95,
  distribution = c("normal", "skewed"),
  G = 1000,
  seed = NULL,
  m_0 = 4,
  max_n = 100000L,
  batch = 1L
)

seq_aipe_cv_sensitivity(
  true_value,
  width,
  conf_level = 0.95,
  distribution = c("gamma", "lognormal"),
  G = 1000,
  seed = NULL,
  m_0 = 4,
  max_n = 100000L,
  batch = 1L
)

seq_aipe_reg_coef_sensitivity(
  true_value,
  width,
  conf_level = 0.95,
  distribution = c("normal", "skewed"),
  G = 1000,
  seed = NULL,
  m_0 = 4,
  max_n = 100000L,
  batch = 1L
)

seq_aipe_sm_sensitivity(
  true_value,
  width,
  conf_level = 0.95,
  distribution = c("normal", "skewed"),
  G = 1000,
  seed = NULL,
  m_0 = 4,
  max_n = 100000L,
  batch = 1L
)

seq_aipe_gini_sensitivity(
  true_value,
  width,
  conf_level = 0.95,
  G = 1000,
  seed = NULL,
  m_0 = 4,
  max_n = 100000L,
  batch = 1L
)

seq_aipe_c_sensitivity(
  true_means,
  contrast_weights,
  width,
  conf_level = 0.95,
  true_sds = 1,
  distribution = c("normal", "skewed"),
  G = 1000,
  seed = NULL,
  m_0 = 4,
  max_n = 100000L,
  batch = 1L
)
```

## Arguments

- true_value:

  The population value of the effect size under the data generating
  mechanism: the standardized mean difference (any finite number), the
  correlation (strictly between -1 and 1), the coefficient of variation
  (positive), the regression slope (any finite number), the standardized
  mean (any finite number), or the Gini index (strictly between 0 and
  1). Within those ranges, a value at which the generated data would
  have too little variation to be used is an error. The check functions
  take a sample whose values agree to ten significant digits to have no
  variation, and the sensitivity functions keep the generated data a
  factor of 100 clear of that: the standardized mean difference and the
  standardized mean must be less than 1e8 in absolute value (the data
  have unit variance), the coefficient of variation must be greater than
  1e-8, and the Gini index greater than about 5.6e-9. The regression
  slope must be less than `1 / .Machine$double.eps` (about 4.5e15) in
  absolute value, beside which the unit variance of the errors cannot be
  represented. A coefficient of variation or Gini index so large that
  its generating parameter (see Details) is not a positive finite number
  is an error as well.

- width:

  The target full width of the confidence interval (the quantity
  \\\omega\\ of Kelley, Darku, and Chattopadhyay, 2018, 2019), in the
  units of the effect size. For the contrast these are the units of the
  response. For the correlation it must be less than 2, and for the Gini
  index less than 1, the lengths of their ranges.

- conf_level:

  The confidence level, strictly between 0 and 1 (default 0.95).

- distribution:

  The data generating mechanism (see Details). `"normal"` or `"skewed"`
  for the standardized mean difference, the correlation, the regression
  slope, the standardized mean, and the contrast; `"gamma"` or
  `"lognormal"` for the coefficient of variation. The Gini index
  function always generates lognormal data and has no such argument.

- G:

  The number of Monte Carlo replications, a whole number of at least 1
  (default 1000). See Details for the Monte Carlo error this implies.

- seed:

  `NULL` (the default) or a single whole number. With a number the
  result is reproducible, and the caller's random number stream is left
  exactly as it was found. With `NULL` the draws come from the current
  stream.

- m_0:

  The minimum pilot sample size, a whole number of at least 4 (default
  4), passed to the matching `seq_aipe_*_pilot()` function. The pilot
  sample size is the larger of `m_0` and \\\lceil 2 z/\omega \rceil\\,
  where \\z\\ is the standard normal quantile for `conf_level` (1.96 at
  0.95). Like every sample size here (`max_n`, `batch`, and the rows
  `mean_n`, `median_n`, and `sd_n`), it counts observations per group
  for the standardized mean difference and the contrast, pairs for the
  correlation and the regression slope, and observations otherwise.

- max_n:

  The largest sample size a replication may reach, in the unit described
  under `m_0` (per group, for example, for the standardized mean
  difference and the contrast, not in total): a whole number no smaller
  than the pilot sample size (default 100,000). A replication that has
  not met the stopping rule at `max_n` is stopped there, exactly,
  whatever the batch size (the last batch is shortened). It is an error
  for `max_n` to be below the pilot sample size. The row `pct_hit_max_n`
  reports how often the cap, not the rule, ended a replication; when it
  is above zero the other rows describe a truncated procedure and should
  be read with care.

- batch:

  The number of observations added between checks of the stopping rule
  after the pilot sample, a whole number of at least 1 (default 1; see
  Details), in the unit described under `m_0`: for the standardized mean
  difference and the contrast it is the number added to every group; for
  the correlation and the regression slope it is the number of pairs.

- true_means:

  For the contrast: the population means, one per group (at least two
  groups), each less than 1e8 times its group's standard deviation in
  absolute value (for the reason given under `true_value`).

- contrast_weights:

  For the contrast: the contrast coefficients \\c_k\\, one per group and
  not all zero, as in
  [`seq_aipe_c`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md).
  `true_means`, `contrast_weights`, and `true_sds` (when it has one
  entry per group) are matched to the groups by position; if more than
  one of them is named, the names must be the same and in the same
  order, or the function stops.

- true_sds:

  For the contrast: the population standard deviations, either one
  positive number shared by every group (default 1) or one per group.
  Their squares must be positive finite numbers.

## Value

An object of class `c("smsd_sensitivity", "data.frame")`: a data frame
with columns `term` (character) and `value` (numeric) that can be used
as any data frame, with a [`print()`](https://rdrr.io/r/base/print.html)
method that shows every value in fixed notation (see
[`smsd_sensitivity_methods`](https://yelleknek.github.io/SMSD/reference/smsd_sensitivity_methods.md)).
It has these 13 rows, in this order:

- `true_value`:

  The population value of the effect size. For the contrast, the
  contrast of `true_means`.

- `target_width`:

  The requested full width, `width`.

- `conf_level`:

  The confidence level, a proportion.

- `replications`:

  The number of replications run, `G`.

- `mean_n`, `median_n`, `sd_n`:

  The mean, median, and standard deviation of the final sample size
  across replications: observations for the coefficient of variation,
  the standardized mean, and the Gini index; pairs for the correlation
  and the regression slope; observations per group for the standardized
  mean difference and the contrast. `sd_n` is `NA` when `G` is 1.

- `mean_width`, `median_width`:

  The mean and median full width of the final confidence interval, in
  the units of the effect size, over the replications that ended with an
  interval (see below); `NA` when none did.

- `pct_width_within`:

  The percentage (0 to 100) of replications whose final interval is no
  wider than `width`.

- `coverage`:

  The proportion (0 to 1) of replications whose final interval contains
  the true value.

- `pct_hit_max_n`:

  The percentage (0 to 100) of replications stopped by `max_n` before
  the stopping rule was met.

- `pct_variance_unusable`:

  The percentage (0 to 100) of replications in which, at one or more
  checks, the check function reported its variance estimate as unusable:
  the raw estimate was zero, negative, or not finite, so that no
  variance had been estimated and the stopping rule was reported as not
  met at that check. A replication is counted once however many of its
  checks were of this kind, and the check at which `max_n` ended it
  counts as well. A positive estimate, however small, is usable and is
  not counted.

A replication that the stopping rule ends always has an interval. A
replication that `max_n` ends can finish on a check that gives no
interval, because the estimate is not finite or the variance estimate is
unusable (the check function then reports its limits as `NA`). Such a
replication counts as not covering and as not within the target width,
so that `coverage` and `pct_width_within` remain proportions of all `G`
replications; it is left out of `mean_width` and `median_width`; and it
is counted like any other in the remaining rows, `pct_hit_max_n` among
them. No row becomes `NA` on this account, except the two width rows
when not one replication ended with an interval.

## Details

One call does what the pilot, collect, and check loop does by hand: it
draws the pilot sample that the matching pilot function gives, evaluates
the stopping rule with the matching check function, adds `batch`
observations, and repeats until the rule is met or the sample reaches
`max_n`.

**Data generating mechanisms.** The true value is exact for each
mechanism.

- `seq_aipe_smd_sensitivity()`: two independent groups with unit
  variance, the first with mean `true_value` and the second with mean 0.

- `seq_aipe_r_sensitivity()`: pairs \\(x, y)\\ with \\y = \rho x +
  \sqrt{1 - \rho^2}\\ e\\, where \\x\\ and \\e\\ are independent with
  unit variance and \\\rho\\ is `true_value`. The procedure studied is
  the Pearson correlation.

- `seq_aipe_reg_coef_sensitivity()`: \\y = \beta x + e\\, where \\x\\
  and \\e\\ are independent with unit variance and \\\beta\\ is
  `true_value`.

- `seq_aipe_sm_sensitivity()`: one sample with mean `true_value` and
  unit variance.

- `seq_aipe_c_sensitivity()`: independent groups with means `true_means`
  and standard deviations `true_sds`. The true value is the contrast
  \\\sum_k c_k \mu_k\\. The groups are balanced, as
  [`seq_aipe_c`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md)
  requires: each addition is `batch` observations from every group.

- `seq_aipe_cv_sensitivity()`: gamma data with shape \\1/\kappa^2\\, or
  lognormal data with log standard deviation \\\sqrt{\log(1 +
  \kappa^2)}\\, where \\\kappa\\ is `true_value`.

- `seq_aipe_gini_sensitivity()`: lognormal data with log standard
  deviation \\\sqrt{2}\\\Phi^{-1}((g + 1)/2)\\, where \\g\\ is
  `true_value` and \\\Phi^{-1}\\ is the standard normal quantile
  function.

For the first five, `distribution = "normal"` draws the unit variance
variables from the standard normal distribution and
`distribution = "skewed"` draws them from a chi square distribution with
3 degrees of freedom, standardized to mean 0 and variance 1 (skewness
about 1.63).

**Batch size.** With `batch = 1`, the default, the stopping rule is
checked after every observation. That is the purely sequential
procedure, and it is what the simulations of Kelley, Darku, and
Chattopadhyay (2018, 2019) use; `batch` is the quantity they call \\m'\\
(2018, Figure 1). A study that can check only after every 5, 10, or 20
observations stops later than the first sample size at which the rule
would have been met; set `batch` to that number to see the overshoot to
expect. For a standardized mean difference of 0.5 at `width = 0.4`, for
example, the mean final sample size per group was 199.1 with
`batch = 1`, 201.2 with 5, 203.6 with 10, and 209.5 with 20 (4,700 to
6,000 replications each, standard errors about 0.1), and coverage stayed
between 0.949 and 0.952: checking less often cost about half a batch per
group.

**Coverage is an asymptotic property.** The procedures deliver the
requested width by construction: the final interval is never wider than
`width` unless `max_n` ends the replication first. Coverage is a
different matter. The articles prove that it approaches `conf_level` as
`width` shrinks (Kelley, Darku, & Chattopadhyay, 2018, p. 232; 2019, p.
495), that is, as the sample that the stopping rule asks for grows. When
that sample is small, coverage can fall well short of `conf_level`. How
small it is depends on the true value as well as on `width`: a strong
correlation is estimated with so little sampling error that the rule is
met in a small sample even at a width that would be narrow for a
moderate correlation.

The table gives what these functions measured at `conf_level = 0.95`,
with the defaults and normal data unless a distribution is named: the
mean final sample size (the row `mean_n`), the coverage, and the
percentage of replications in which the variance estimate was unusable
at one or more checks (the row `pct_variance_unusable`, explained
below). Each row rests on 800 to 34,000 replications; the Monte Carlo
standard error of the coverage is below 0.008 in every row but three
with a mean *n* above 800 (the correlation of 0.9 at 0.025, the slope of
2 at 0.1, and the lognormal coefficient of variation at 0.05), where it
is 0.009 to 0.012.

|  |  |  |  |  |
|----|----|----|----|----|
| Effect size and true value | Width | Mean *n* | Coverage | Percent unusable |
| Standardized mean difference, 0.5 | 0.8 | 51 | 0.95 | \<0.1 |
| Standardized mean difference, 0.5 | 0.4 | 199 | 0.95 | 0.0 |
| Standardized mean, 0.4 | 0.5 | 67 | 0.95 | \<0.1 |
| Standardized mean, 0.4 | 0.25 | 267 | 0.95 | 0.0 |
| Contrast, 0.5 | 1 | 23 | 0.94 | 0.0 |
| Contrast, 0.5 | 0.5 | 92 | 0.95 | 0.0 |
| Correlation, 0.4 | 0.4 | 55 | 0.78 | 16.5 |
| Correlation, 0.4 | 0.2 | 263 | 0.93 | 0.9 |
| Correlation, 0.4 | 0.1 | 1,079 | 0.95 | 0.0 |
| Correlation, 0.9 | 0.1 | 75 | 0.45 | 79.5 |
| Correlation, 0.9 | 0.05 | 153 | 0.63 | 40.4 |
| Correlation, 0.9 | 0.025 | 819 | 0.90 | 1.3 |
| Regression slope, 0.5 | 0.4 | 80 | 0.83 | 20.1 |
| Regression slope, 0.5 | 0.2 | 377 | 0.94 | 1.3 |
| Regression slope, 2 | 0.2 | 217 | 0.61 | 65.3 |
| Regression slope, 2 | 0.1 | 1,382 | 0.87 | 22.4 |
| Coefficient of variation, 0.5, gamma | 0.2 | 52 | 0.88 | 0.0 |
| Coefficient of variation, 0.5, gamma | 0.1 | 217 | 0.92 | 0.0 |
| Coefficient of variation, 0.5, lognormal | 0.2 | 62 | 0.79 | 0.0 |
| Coefficient of variation, 0.5, lognormal | 0.1 | 322 | 0.86 | 0.0 |
| Coefficient of variation, 0.5, lognormal | 0.05 | 1,561 | 0.90 | 0.0 |
| Gini index, 0.3 | 0.2 | 30 | 0.89 | 0.0 |
| Gini index, 0.3 | 0.1 | 88 | 0.92 | 0.0 |

(The contrast has weights 1, -0.5, and -0.5, means 0.5, 0, and 0, and
unit standard deviations. The coverage for the Gini index is that of the
interval
[`seq_aipe_gini`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md)
reports. The fixed width interval of Chattopadhyay and De, 2016, the
estimate plus and minus half the width, contained the index in 0.986 and
0.960 of the same studies; see
[`seq_aipe_gini`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md).)
These values hold for the true values shown, not for an effect size in
general: the rows for a correlation of 0.9 and for a slope of 2 are far
from those for a correlation of 0.4 and a slope of 0.5, at widths that
are no wider. Skewed data can lower coverage further: for a slope of 0.5
at `width = 0.4` it was 0.71 with the skewed data.

Read the `coverage` row before relying on an interval at a planned
width, and run the function at each true value that is plausible for the
study, the strongest included. When coverage is low, plan a narrower
`width` (a larger study), raise the minimum pilot size `m_0`, or both.
With `m_0 = 40` at `width = 0.4`, coverage rose from 0.78 to 0.90 for a
correlation of 0.4 and from 0.83 to 0.92 for a slope of 0.5, while the
mean final sample size rose from 55 to 64 and from 80 to 90. A larger
pilot sample is a weaker remedy when the true value is strong: for a
correlation of 0.9 at `width = 0.1`, `m_0 = 100` raised coverage only
from 0.45 to 0.69.

**When the variance estimate is unusable.** In a small sample the raw
estimate of the asymptotic variance can be zero, negative, or not
finite: no variance has been estimated. The check functions report such
a check as not met and set the attribute `variance_unusable` to `TRUE`,
so sampling continues (see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md);
this is a deliberate departure from the stopping rule as printed in the
articles, under which the procedure would stop at once with an interval
far narrower than the data support). A positive estimate, however small,
is used as the articles prescribe and is not counted here.

The row `pct_variance_unusable` is the percentage of replications in
which this happened at one or more checks; it is the last column of the
table. Like coverage, it depends on the true value as well as on the
width, and it is not confined to wide targets or to one effect size. It
was 16 percent for a correlation of 0.4 at `width = 0.4` but 79 percent
for a correlation of 0.9 at `width = 0.1`; 20 percent for a slope of 0.5
at `width = 0.4` but 65 percent for a slope of 2 at `width = 0.2`; and
15 percent for a standardized mean of 1 at `width = 1`, where coverage
was 0.86. With the skewed data it was 37 percent for the slope of 0.5 at
`width = 0.4`.

Continuing past such a check does not by itself make the final interval
trustworthy. The estimate on which a replication finally stops can still
be far too small, which makes the final interval too narrow: for the
correlation of 0.9 at `width = 0.1`, the final intervals had a mean
width of 0.028, about a third of the 0.086 that the true variance gives
at the mean final sample size, and their coverage was 0.45. With a very
strong relationship the estimate can also stay unusable for a long time.
For a correlation of 0.99 at `width = 0.1` (a pilot sample of 40), 20
replications ended at a mean of 7,500 pairs, each taking about 6.5
seconds on the computer used for these runs, and 1 of their 20 final
intervals contained the true value. Set `max_n` when exploring such
conditions, and read `pct_hit_max_n`.

**Monte Carlo error.** Every summary is an estimate from `G`
replications. The standard error of `coverage` is about \\\sqrt{p(1 -
p)/G}\\ for a coverage of \\p\\: 0.007 at the default `G = 1000` and
0.002 at `G = 10000` when \\p\\ is near 0.95. The standard error of
`mean_n` is `sd_n` divided by \\\sqrt{G}\\.

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
accuracy in parameter estimation for population correlation
coefficients. *Psychological Methods*, *24*(4), 492–515.
[doi:10.1037/met0000203](https://doi.org/10.1037/met0000203)

Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
pre-specified error bound. *Econometrics*, *4*(3), 30.
[doi:10.3390/econometrics4030030](https://doi.org/10.3390/econometrics4030030)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology*, *59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

## See also

- [`seq_aipe_smd`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md),
  [`seq_aipe_r`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md),
  [`seq_aipe_cv`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md),
  and
  [`seq_aipe_reg_coef`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md):
  four of the check functions these run;

- [`seq_aipe_sm`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm.md),
  [`seq_aipe_gini`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md),
  and
  [`seq_aipe_c`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md):
  the other three;

- [`smsd_sensitivity_methods`](https://yelleknek.github.io/SMSD/reference/smsd_sensitivity_methods.md):
  how the summary prints;

- [`vignette("sensitivity-analysis", package = "SMSD")`](https://yelleknek.github.io/SMSD/articles/sensitivity-analysis.md):
  a worked sensitivity analysis;

- [`ss_aipe_smd_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd_sensitivity.html):
  in DMAR, the sensitivity function for fixed sample size planning for
  the standardized mean difference;

- [`ss_aipe_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.html):
  the same for the correlation;

- [`ss_aipe_cv_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.html):
  for the coefficient of variation;

- [`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.html):
  for the regression slope;

- [`ss_aipe_sm_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sm_sensitivity.html):
  for the standardized mean;

- [`ss_aipe_c_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.html):
  for the contrast.

Other seq_aipe_smd:
[`seq_aipe_smd()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md),
[`seq_aipe_smd_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd_pilot.md)

Other seq_aipe_r:
[`seq_aipe_r()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md),
[`seq_aipe_r_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r_pilot.md)

Other seq_aipe_cv:
[`seq_aipe_cv()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md),
[`seq_aipe_cv_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv_pilot.md)

Other seq_aipe_reg_coef:
[`seq_aipe_reg_coef()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md),
[`seq_aipe_reg_coef_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef_pilot.md)

Other seq_aipe_sm:
[`seq_aipe_sm()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm.md),
[`seq_aipe_sm_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm_pilot.md)

Other seq_aipe_gini:
[`seq_aipe_gini()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md),
[`seq_aipe_gini_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini_pilot.md)

Other seq_aipe_c:
[`seq_aipe_c()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md),
[`seq_aipe_c_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# G is kept very small in these examples so that they run in a fraction
# of a second, and the targets are wide for the same reason. A planning
# result that will be reported deserves G = 10000, at the width the
# study will use.

# A training program compared with a control group; the true
# standardized mean difference in engagement is 0.5.
seq_aipe_smd_sensitivity(true_value = 0.5, width = 1.2, G = 10, seed = 113)
#> term                   value
#> true_value             0.5
#> target_width           1.2
#> conf_level             0.95
#> replications           10
#> mean_n                 23.8
#> median_n               23
#> sd_n                   1.932184
#> mean_width             1.170544
#> median_width           1.171495
#> pct_width_within       100
#> coverage               1
#> pct_hit_max_n          0
#> pct_variance_unusable  0

# Job satisfaction and performance ratings with a true correlation of
# 0.4. Ten replications say little about coverage (its standard error
# here is about 0.15); in 10,000 replications at this very wide target
# it was 0.65, far below 0.95.
seq_aipe_r_sensitivity(true_value = 0.4, width = 0.8, G = 10, seed = 113)
#> term                   value
#> true_value             0.4
#> target_width           0.8
#> conf_level             0.95
#> replications           10
#> mean_n                 21.7
#> median_n               21
#> sd_n                   8.459708
#> mean_width             0.7284184
#> median_width           0.737364
#> pct_width_within       100
#> coverage               1
#> pct_hit_max_n          0
#> pct_variance_unusable  50

# Response times with a true coefficient of variation of 0.5.
seq_aipe_cv_sensitivity(true_value = 0.5, width = 0.5, G = 10, seed = 113)
#> term                   value
#> true_value             0.5
#> target_width           0.5
#> conf_level             0.95
#> replications           10
#> mean_n                 11.1
#> median_n               11
#> sd_n                   2.078995
#> mean_width             0.2920956
#> median_width           0.3157522
#> pct_width_within       100
#> coverage               1
#> pct_hit_max_n          0
#> pct_variance_unusable  0

# Hours of tutoring predicting a test score, with a true slope of 0.5.
seq_aipe_reg_coef_sensitivity(true_value = 0.5, width = 0.8, G = 10,
                              seed = 113)
#> term                   value
#> true_value             0.5
#> target_width           0.8
#> conf_level             0.95
#> replications           10
#> mean_n                 14.8
#> median_n               14.5
#> sd_n                   8.978988
#> mean_width             0.5470225
#> median_width           0.6225386
#> pct_width_within       100
#> coverage               0.8
#> pct_hit_max_n          0
#> pct_variance_unusable  50

# A change score with a true standardized mean of 0.4.
seq_aipe_sm_sensitivity(true_value = 0.4, width = 1, G = 10, seed = 113)
#> term                   value
#> true_value             0.4
#> target_width           1
#> conf_level             0.95
#> replications           10
#> mean_n                 18.3
#> median_n               17
#> sd_n                   3.267687
#> mean_width             0.9510451
#> median_width           0.9529519
#> pct_width_within       100
#> coverage               0.9
#> pct_hit_max_n          0
#> pct_variance_unusable  10

# Household incomes with a true Gini index of 0.3.
seq_aipe_gini_sensitivity(true_value = 0.3, width = 0.4, G = 10, seed = 113)
#> term                   value
#> true_value             0.3
#> target_width           0.4
#> conf_level             0.95
#> replications           10
#> mean_n                 12.6
#> median_n               12.5
#> sd_n                   1.429841
#> mean_width             0.2101984
#> median_width           0.2146468
#> pct_width_within       100
#> coverage               0.9
#> pct_hit_max_n          0
#> pct_variance_unusable  0

# A new teaching method against the average of two existing methods,
# with exam scores that have a standard deviation of 10. The class sizes
# can only be checked after every 5 students per method.
seq_aipe_c_sensitivity(true_means = c(75, 70, 70),
                       contrast_weights = c(1, -0.5, -0.5),
                       width = 12, true_sds = 10, G = 10, seed = 113,
                       batch = 5)
#> term                   value
#> true_value             5
#> target_width           12
#> conf_level             0.95
#> replications           10
#> mean_n                 18
#> median_n               19
#> sd_n                   2.108185
#> mean_width             10.84834
#> median_width           10.91791
#> pct_width_within       100
#> coverage               1
#> pct_hit_max_n          0
#> pct_variance_unusable  0
```
