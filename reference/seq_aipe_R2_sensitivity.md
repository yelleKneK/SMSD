# Monte Carlo Sensitivity Analysis for Sequential AIPE: Squared Multiple Correlation Coefficient

Runs the sequential accuracy in parameter estimation (AIPE) procedure
for the squared multiple correlation coefficient
([`seq_aipe_R2_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_pilot.md),
then
[`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md)
until the stopping rule is met) `G` times on multivariate normal data
with a known population value, and summarizes the realized sample size,
the realized confidence interval width, and the empirical coverage. Use
it before a study to see what sample size to expect and how close the
coverage is to the confidence level at the planned width and number of
predictors: the procedure's guarantee is asymptotic, and its coverage
falls, in some conditions far below the confidence level, as the number
of predictors grows relative to the sample size the width calls for (see
[`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md),
which attaches a note to an interval reported from fewer than \\p
\max(20, 2p)\\ observations; a design whose `mean_n` and `median_n` here
are below that number is one in which a study will usually end with that
note).

## Usage

``` r
seq_aipe_R2_sensitivity(
  true_value,
  width,
  p,
  conf_level = 0.95,
  ci_method = c("wald", "log_transformed"),
  G = 1000,
  seed = NULL,
  m_0 = p + 2,
  max_n = 100000L,
  batch = 1L
)
```

## Arguments

- true_value:

  The population squared multiple correlation coefficient under the data
  generating mechanism. A single number of at least 0 and less than 1.

- width:

  The target full width of the confidence interval (\\\omega\\ in the
  article). A single number greater than 0 and less than 1.

- p:

  The number of predictors (\\k\\ in the article). A whole number of at
  least 1.

- conf_level:

  The confidence level (default 0.95).

- ci_method:

  `"wald"` (the default) or `"log_transformed"`, as in
  [`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md).

- G:

  The number of Monte Carlo replications (default 1000). A whole number
  of at least 1.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single whole number that
  [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts to make the
  result reproducible; the caller's random number stream is then left as
  it was found.

- m_0:

  The smallest pilot sample size allowed, passed to
  [`seq_aipe_R2_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_pilot.md).
  A whole number of at least `p + 2`, which is also the default.

- max_n:

  The largest sample size a replication may reach: a replication in
  which the rule is still not met with `max_n` observations stops there.
  It must be at least the pilot sample size, and at least `p + 3` with
  the log transformed interval (which is not defined with `p + 2`
  observations). The cap is exact for any `batch` (the last batch is
  shortened).

- batch:

  The number of observations added between checks of the stopping rule
  after the pilot sample (default 1).

## Value

An object of class `c("smsd_sensitivity", "data.frame")`: a data frame
with columns `term` and `value`, which can be used as any data frame and
which prints its values in fixed notation. See
[`smsd_sensitivity_methods`](https://yelleknek.github.io/SMSD/reference/smsd_sensitivity_methods.md)
for the print method. It has one row per summary, in this order:

- `true_value`:

  The population squared multiple correlation supplied.

- `target_width`:

  The target width supplied.

- `conf_level`:

  The confidence level supplied.

- `replications`:

  The number of replications run.

- `mean_n`, `median_n`, `sd_n`:

  The mean, median, and standard deviation of the realized sample size.
  `sd_n` is `NA` when `G` is 1.

- `mean_width`, `median_width`:

  The mean and median of the realized confidence interval width, over
  the replications that ended with an interval.

- `pct_width_within`:

  The percentage (0 to 100) of replications whose realized width is at
  most the target.

- `coverage`:

  The proportion (0 to 1) of replications whose interval contains
  `true_value`. A replication stopped by `max_n` at a check that has no
  interval (see
  [`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md))
  counts as not covering and as not within the target width; this cannot
  happen with the data generated here.

- `pct_hit_max_n`:

  The percentage (0 to 100) of replications stopped by `max_n` before
  the stopping rule was met (a replication in which the rule is met with
  exactly `max_n` observations is not counted). When it is above zero
  the other summaries are censored and should be read with care.

- `pct_variance_unusable`:

  The percentage (0 to 100) of replications in which, at one or more
  checks up to the stopping point, the variance estimate was unusable:
  zero (a sample value of 0 or 1 to within the precision of the data) or
  not finite (see
  [`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md)).
  For the data generated here this is 0 in practice.

## Details

Each replication draws \\p\\ independent standard normal predictors and
a response equal to \\\sqrt{\mathrm{P}^2 / p}\\ times their sum plus
independent normal error with variance \\1 - \mathrm{P}^2\\, where
\\\mathrm{P}^2\\ is `true_value`. This is the design of the simulations
in Kelley, Darku, and Chattopadhyay (2019, p. 506). Under multivariate
normality the distribution of the sample squared multiple correlation
depends only on \\\mathrm{P}^2\\, \\p\\, and the sample size, so the
results apply to any means, variances, and correlations among the
predictors with the same population value.

The stopping rule is checked at the pilot sample size and then after
every `batch` further observations. With `batch = 1` this is the purely
sequential procedure of the article; a larger value shows the overshoot
to expect from a study that checks the rule only after every `batch`
observations. The function gives the results that a loop of
[`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md)
calls on the same observations would give, but it evaluates the rule
from running sums, so a narrow target that needs thousands of
observations per replication remains practical.

The function generates multivariate normal data only, the assumption
under which the procedure is derived; see
[`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md)
for measured coverage when that assumption fails.

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
accuracy in parameter estimation for population correlation
coefficients. *Psychological Methods*, *24*(4), 492–515.
[doi:10.1037/met0000203](https://doi.org/10.1037/met0000203)

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research*, *43*(4),
524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

## See also

[`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md),
[`seq_aipe_R2_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_pilot.md);
[`smsd_sensitivity_methods`](https://yelleknek.github.io/SMSD/reference/smsd_sensitivity_methods.md)
for how the result prints; in DMAR,
[`ss_aipe_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2_sensitivity.html)
for the fixed sample size analog.

Other seq_aipe_R2:
[`seq_aipe_R2()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md),
[`seq_aipe_R2_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Three predictors, a population squared multiple correlation of .40, and
# a 95% confidence interval no wider than 0.25. G is small to keep the
# example quick; a reported analysis deserves G = 10000.
seq_aipe_R2_sensitivity(true_value = 0.40, width = 0.25, p = 3, G = 100,
                        seed = 113)
#> term                   value
#> true_value             0.4
#> target_width           0.25
#> conf_level             0.95
#> replications           100
#> mean_n                 143.38
#> median_n               148
#> sd_n                   13.39213
#> mean_width             0.2440061
#> median_width           0.2444849
#> pct_width_within       100
#> coverage               0.94
#> pct_hit_max_n          0
#> pct_variance_unusable  0

# The same study checking the stopping rule after every 20 participants.
seq_aipe_R2_sensitivity(true_value = 0.40, width = 0.25, p = 3, G = 100,
                        seed = 113, batch = 20)
#> term                   value
#> true_value             0.4
#> target_width           0.25
#> conf_level             0.95
#> replications           100
#> mean_n                 149.8
#> median_n               152
#> sd_n                   6.90191
#> mean_width             0.2401216
#> median_width           0.2419052
#> pct_width_within       100
#> coverage               0.93
#> pct_hit_max_n          0
#> pct_variance_unusable  0
```
