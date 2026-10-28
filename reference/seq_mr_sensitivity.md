# Monte Carlo Sensitivity Analysis for Minimum Risk Point Estimation

These functions run a minimum risk point estimation (MRPE) procedure
from its pilot sample to its stopping point `G` times, on data simulated
from a population whose true value is known, and summarize what
happened: how many observations the procedure took compared with the
best fixed sample size, and how the realized risk compares with the risk
of that best fixed design. Use them before collecting data, to see what
a choice of the loss constant `A` and the sampling cost `c` implies, and
to see how the procedure behaves when checking in batches or when the
population is not the one the method assumes.

## Usage

``` r
seq_mr_cv_sensitivity(
  true_value,
  A,
  c,
  distribution = c("gamma", "lognormal", "normal"),
  variance = c("corrected", "published"),
  G = 1000,
  seed = NULL,
  m_0 = 4,
  gamma = 0.49,
  max_n = 100000L,
  batch = 1L
)

seq_mr_smd_sensitivity(
  true_value,
  A,
  c,
  distribution = c("normal", "skewed"),
  G = 1000,
  seed = NULL,
  m_0 = 4,
  gamma = 0.49,
  max_n = 100000L,
  batch = 1L
)
```

## Arguments

- true_value:

  The population value. For the coefficient of variation it must be
  positive; for the standardized mean difference it must be less than
  100 million in absolute value. A value so extreme that the population
  cannot be generated, or that the benchmark is not a positive finite
  number, is an error.

- A:

  The loss constant, the price of one unit of squared estimation error.
  A single positive number.

- c:

  The cost of sampling one observation. A single positive number or, for
  the standardized mean difference, two positive numbers, `c(c1, c2)`,
  the cost in each group.

- distribution:

  The population to simulate from. For the coefficient of variation:
  `"gamma"` (the default), `"lognormal"`, or `"normal"` (mean 1), each
  with coefficient of variation `true_value`. For the standardized mean
  difference: `"normal"` (the default; unit variance in both groups,
  means `true_value` and 0) or `"skewed"` (the same means and variance,
  with a standardized chi square distribution on 3 degrees of freedom in
  place of the normal).

- variance:

  For `seq_mr_cv_sensitivity()`, the variance estimator passed to
  [`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md):
  `"corrected"` (the default) or `"published"`.

- G:

  The number of Monte Carlo replications: a whole number of at least 2.
  The default is 1000. The examples use far fewer so that they run
  quickly; a reported analysis deserves at least the default.

- seed:

  A single whole number that
  [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts, to make
  the result reproducible, or `NULL` (the default) to use the current
  state of the random number generator. When a seed is given, the
  generator's state is restored afterward.

- m_0:

  The smallest pilot size allowed (per group for the standardized mean
  difference): a whole number of at least 4, the default.

- gamma:

  The tuning constant of the correction term, strictly between 0 and
  1/2; default 0.49.

- max_n:

  A cap on the sample size (on each group's size for the standardized
  mean difference): a replication that reaches it is stopped there. A
  whole number no smaller than the pilot size. The cap is exact: the
  last batch is shortened if needed. The default is large; see "Run
  time" in Details before using it with constants whose optimal sample
  size is in the thousands.

- batch:

  The number of observations added between checks of the stopping rule:
  a whole number of at least 1. The default, 1, is the purely sequential
  procedure of the articles. For the standardized mean difference it is
  the number added per group with one cost, and the number added to the
  less costly group with two costs.

## Value

An object of class `c("smsd_sensitivity", "data.frame")` (see
[`smsd_sensitivity_methods`](https://yelleknek.github.io/SMSD/reference/smsd_sensitivity_methods.md)):
a data frame with columns `term` (character) and `value` (numeric), one
row per summary, that prints in fixed notation. The values are stored
unrounded. For `seq_mr_cv_sensitivity()` the rows are

- `true_value`, `loss_A`, `cost`, `replications`:

  the inputs;

- `optimal_n`, `minimum_risk`:

  the benchmark: the best fixed sample size and the risk there;

- `mean_n`, `median_n`, `sd_n`:

  the mean, median, and standard deviation of the final sample size;

- `ratio_n`:

  `mean_n` divided by `optimal_n`;

- `mean_estimate`:

  the mean of the final estimates;

- `mean_risk`:

  the mean realized loss;

- `risk_efficiency`:

  `mean_risk` divided by `minimum_risk`;

- `se_risk_efficiency`:

  the Monte Carlo standard error of `risk_efficiency`;

- `regret`:

  `mean_risk` minus `minimum_risk`, in the currency of `c`;

- `pct_hit_max_n`:

  the percentage (0 to 100) of replications stopped by `max_n`.

For `seq_mr_smd_sensitivity()` the rows are the same, except that `cost`
is `cost_1` and `cost_2`, and each sample size row is reported for both
groups: `optimal_n1`, `optimal_n2`, `mean_n1`, `mean_n2`, `median_n1`,
`median_n2`, `sd_n1`, `sd_n2`, `ratio_n1`, and `ratio_n2`. A replication
stopped by `max_n` whose final data give no estimate (a sample mean that
is not positive, or constant data) is left out of `mean_estimate` and
the four risk rows, which are `NA` when every replication is left out
(`pct_hit_max_n` is then 100).

## Details

The function `seq_mr_cv_sensitivity()` runs
[`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
(Chattopadhyay & Kelley, 2016), and the function
`seq_mr_smd_sensitivity()` runs
[`seq_mr_smd`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)
(Chattopadhyay & Kelley, 2017). Each replication applies the same
stopping rule as the check function, to the same kind of growing sample
a user would collect.

**What a replication does.** It draws the pilot sample
([`seq_mr_cv_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv_pilot.md)
or
[`seq_mr_smd_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd_pilot.md)),
evaluates the stopping rule, and, while the rule is not met, adds
`batch` observations and evaluates it again. For the standardized mean
difference with one cost, each stage adds `batch` observations to each
group. With two costs it follows the article's order: each stage adds
`batch` observations to the less costly group and then brings the
costlier group up to the size the article pairs with it. With
`batch = 1` these are the sizes
[`seq_mr_smd`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)
reports in its rows `next_n1` and `next_n2`, so a replication is the
loop that function's help describes. A stage at which the check function
would report the rule as not met for want of an estimate (a sample mean
that is not positive, constant data, a variance estimate that is zero,
negative, or not finite), or at which the criterion is too large to
represent, is a stage at which the rule is not met: sampling goes on.
The realized loss of a replication is the squared error of its final
estimate priced at `A`, plus what its observations cost.

**Run time.** The stopping rule is evaluated on all the data collected
so far after every batch, so the time for one replication grows faster
than its final sample size, and toward the square of it for large
samples. On the machine used to check the package, with `batch = 1`, one
replication of the coefficient of variation procedure took about 0.004
seconds at a final sample size near 100, 0.07 seconds near 1,100, and
0.5 seconds near 4,000 (about nine minutes for the default `G`);
`batch = 50` brought the last figure to 0.01 seconds. When the optimal
sample size is in the thousands, use a larger `batch`, a smaller `G`, or
a smaller `max_n`.

**The benchmark.** `optimal_n` (or `optimal_n1` and `optimal_n2`) is the
fixed sample size that minimizes the large sample risk for the simulated
population, and `minimum_risk` is the risk there. They are not rounded.
For the coefficient of variation they are \\n_c = (A\\\xi^2/c)^{1/2}\\
and \\2 c\\ n_c\\ (Equations 15 and 16 of Chattopadhyay & Kelley, 2016),
with \\\xi^2\\ the asymptotic variance of the simulated population
(Equation 30 of Kelley, Darku, & Chattopadhyay, 2018) whichever
`variance` is requested. For the standardized mean difference they
minimize the approximate risk of Chattopadhyay and Kelley (2017,
Equation 9, or Equation 17 with two costs). That risk assumes normal
populations. For `distribution = "skewed"` the benchmark instead uses
the large sample mean squared error of the estimate for the skewed
population, \\1/n_1 + 1/n_2 + 3\delta^2 / \\2(n_1 + n_2)\\\\, so that
the risk efficiency shows what the normality assumption costs.

**Reading the result.** A procedure that works as intended has `ratio_n`
and `risk_efficiency` near 1. Both guarantees are asymptotic (they hold
as `A / c` grows). At small sample sizes the correction term makes the
procedure take more than the optimal number of observations (for the
standardized mean difference with `A = 500` and `c = 2`, about 19 per
group against 16), and for a strongly skewed population the coefficient
of variation procedure stops early (for gamma data with a coefficient of
variation of 1.3, `A = 900`, and `c = 1`, about 26 observations against
45). The benchmark is a large sample approximation, so at small sample
sizes `risk_efficiency` can fall below 1; that reflects the
approximation, not a procedure better than the best fixed design.
`risk_efficiency` is also an average of very variable losses, so judge
it against `se_risk_efficiency`. If `pct_hit_max_n` is above zero, some
replications were stopped by `max_n` and not by the rule, and the other
rows are censored.

## References

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

## See also

[`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
and
[`seq_mr_smd`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md),
the procedures simulated. The accuracy in parameter estimation
counterparts:
[`seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md).
How the result prints:
[`smsd_sensitivity_methods`](https://yelleknek.github.io/SMSD/reference/smsd_sensitivity_methods.md).

Other seq_mr_cv:
[`seq_mr_cv()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md),
[`seq_mr_cv_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv_pilot.md)

Other seq_mr_smd:
[`seq_mr_smd()`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md),
[`seq_mr_smd_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A survey team prices squared error at A = 2,000 and pays 1 per
# respondent. If household spending has a coefficient of variation of
# 0.5, how many respondents will the procedure take?
out <- seq_mr_cv_sensitivity(true_value = 0.5, A = 2000, c = 1, G = 100,
                             seed = 113)
out
#> term                value
#> true_value          0.5
#> loss_A              2000
#> cost                1
#> replications        100
#> optimal_n           17.67767
#> minimum_risk        35.35534
#> mean_n              18.3
#> median_n            18
#> sd_n                3.026816
#> ratio_n             1.035204
#> mean_estimate       0.489766
#> mean_risk           36.83699
#> risk_efficiency     1.041908
#> se_risk_efficiency  0.07084226
#> regret              1.481655
#> pct_hit_max_n       0
out$value[out$term == "ratio_n"]
#> [1] 1.035204

# A training study prices squared error at A = 500 and pays 2 per
# participant in each group. For a true standardized mean difference of
# 0.5, checking after every pair of participants:
seq_mr_smd_sensitivity(true_value = 0.5, A = 500, c = 2, G = 100,
                       seed = 113)
#> term                value
#> true_value          0.5
#> loss_A              500
#> cost_1              2
#> cost_2              2
#> replications        100
#> optimal_n1          16.05654
#> optimal_n2          16.05654
#> minimum_risk        128.4523
#> mean_n1             19.24
#> mean_n2             19.24
#> median_n1           19
#> median_n2           19
#> sd_n1               0.4739688
#> sd_n2               0.4739688
#> ratio_n1            1.198266
#> ratio_n2            1.198266
#> mean_estimate       0.5115002
#> mean_risk           130.8696
#> risk_efficiency     1.018819
#> se_risk_efficiency  0.0610766
#> regret              2.417313
#> pct_hit_max_n       0

# The control group costs twice as much per participant, and the rule is
# checked after every fifth participant in the less costly group.
seq_mr_smd_sensitivity(true_value = 0.5, A = 500, c = c(2, 4), G = 100,
                       seed = 113, batch = 5)
#> term                value
#> true_value          0.5
#> loss_A              500
#> cost_1              2
#> cost_2              4
#> replications        100
#> optimal_n1          16.14988
#> optimal_n2          11.29811
#> minimum_risk        154.9844
#> mean_n1             22.15
#> mean_n2             13.03
#> median_n1           22
#> median_n2           13
#> sd_n1               0.857233
#> sd_n2               0.1714466
#> ratio_n1            1.371528
#> ratio_n2            1.15329
#> mean_estimate       0.5792
#> mean_risk           182.6061
#> risk_efficiency     1.178222
#> se_risk_efficiency  0.08753002
#> regret              27.62165
#> pct_hit_max_n       0
```
