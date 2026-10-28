# Monte Carlo Sensitivity Analysis for the Budget Constrained Procedure

Runs the budget constrained sequential procedure of Chattopadhyay,
Bandyopadhyay, Kelley, and Padalunkal (2025) from pilot to stopping `G`
times, for population means, standard deviations, costs, and a budget
that the user supplies, and summarizes what the procedure delivers: how
the budget is spent, how close the final sample sizes and the variance
of the contrast estimator come to the optimal allocation of
[`seq_budget_optimal`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md),
the coverage of the confidence interval, and, for a given margin, how
often the noninferiority test rejects. Use it before a study to see what
a budget buys.

## Usage

``` r
seq_budget_sensitivity(
  true_means,
  true_sds,
  contrast_weights,
  cost_per_group,
  budget,
  conf_level = 0.95,
  distribution = c("normal", "skewed"),
  margin = NULL,
  alpha_level = 0.05,
  G = 1000,
  seed = NULL,
  m_0 = 5,
  procedure = c("reallocate", "first_passage")
)
```

## Arguments

- true_means:

  Numeric vector of length *K* (at least 2): the population means. In a
  group with a nonzero weight the mean must be less than `1e6` standard
  deviations (`true_sds`) from zero. The package treats scores that
  agree to ten significant digits as showing no variation, and beside a
  larger mean the simulated studies would begin to depend on where the
  scores are located. Subtracting a common constant from the means
  leaves a contrast whose weights sum to zero unchanged.

- true_sds:

  Numeric vector of length *K*: the population standard deviations, each
  positive.

- contrast_weights:

  Numeric vector of length *K* giving the contrast coefficients, not all
  zero.

- cost_per_group:

  Numeric vector of length *K* giving the cost of one observation from
  each population, each positive. The four vector arguments are matched
  by position; if more than one is named, the names must be the same and
  in the same order.

- budget:

  The total budget, pilot included: a single positive number. The pilot,
  which costs `m_0 * sum(cost_per_group)`, must not cost more than the
  budget, and the budget must not pay for more than one million
  observations in the cheapest group (see Details).

- conf_level:

  The confidence level of the two-sided interval whose coverage is
  reported (default 0.95).

- distribution:

  The shape of the populations: `"normal"` (default) or `"skewed"`; an
  unambiguous abbreviation is accepted.

- margin:

  Optional noninferiority margin: the size of the margin, a single
  positive number, as in
  [`decide_seq_budget_noninferiority`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md).
  When supplied, the rejection rate of that test at this margin is
  reported. It is the power of the test when the true contrast is above
  `-margin` and its Type I error rate when the true contrast is at or
  below `-margin`.

- alpha_level:

  The one-sided level of the noninferiority test: a single number
  strictly between 0 and 0.5 (default 0.05, the default of the decision
  functions). It is used only when `margin` is supplied, but is checked
  either way.

- G:

  Number of Monte Carlo replications: a whole number of at least 2
  (default 1000).

- seed:

  Optional whole number that
  [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts. When
  supplied the replications are reproducible and the random number
  generator state of the session is restored on exit; `NULL` (default)
  uses the current state.

- m_0:

  Pilot sample size per population: a whole number of at least 2
  (default 5, as in
  [`seq_budget_pilot`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md)).

- procedure:

  Which reading of the procedure to run (an unambiguous abbreviation is
  accepted). `"reallocate"` (default) is the package's: every group is
  judged on the current data at each check and the next observation goes
  to the group furthest below its estimated optimal size. The other
  choice, `"first_passage"`, follows the article's rule: a group is
  closed for good the first time it meets its criterion, and the open
  groups are visited in rounds (`closed` in
  [`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
  with `order = "rounds"` in
  [`seq_budget_next`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)).

## Value

An object of class `c("smsd_sensitivity", "data.frame")`: a data frame
with columns `term` and `value`, one row per quantity, that prints in
fixed notation (see
[`smsd_sensitivity_methods`](https://yelleknek.github.io/SMSD/reference/smsd_sensitivity_methods.md))
and can otherwise be used as any data frame. With *K* populations the
rows are, in order:

- `true_contrast`:

  the population contrast, `sum(contrast_weights * true_means)`.

- `budget`, `conf_level`, `replications`:

  the inputs.

- `coverage`:

  the proportion (0 to 1) of replications whose interval contains the
  true contrast.

- `mean_total_cost`, `sd_total_cost`:

  the mean and standard deviation of the final total cost.

- `pct_within_budget`:

  the percentage (0 to 100) of replications whose final cost is no more
  than the budget.

- `pct_rule_met`:

  the percentage (0 to 100) of replications that ended with every group
  having met its criterion; the others ended because the budget could
  not pay for another needed observation.

- `mean_n_1`, ..., `mean_n_K`:

  the mean final sample size of each group.

- `optimal_n_1`, ..., `optimal_n_K`:

  for comparison, the optimal sample sizes for the true standard
  deviations, as
  [`seq_budget_optimal`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md)
  returns them.

- `mean_estimated_variance`:

  the mean over replications of the estimated variance of the contrast
  estimator at stopping (the `variance` row of the final check).

- `mean_true_variance`:

  the mean over replications of the true variance of the contrast
  estimator at the final sample sizes,
  `sum(contrast_weights^2 * true_sds^2 / n)`. It is never below
  `optimal_variance`; the excess is the price of not knowing the
  standard deviations in advance.

- `empirical_variance`:

  the mean squared deviation of the contrast estimate from the true
  contrast over replications: the variance the procedure actually
  delivered. It has more Monte Carlo error than the two rows above.

- `optimal_variance`:

  the smallest variance the budget allows with known standard
  deviations, \\V_o\\.

- `margin`, `alpha_level`:

  the inputs of the noninferiority test (`NA` when no margin was
  supplied).

- `noninferiority_rejection_rate`:

  the proportion (0 to 1) of replications in which inferiority was
  rejected (`NA` when no margin was supplied).

## Details

Each replication draws `m_0` observations per population, then adds one
observation at a time to the group that
[`seq_budget_next`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
names, re-evaluating
[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
after each, until no group that still needs an observation can be paid
for. The total cost can therefore never pass the budget;
`pct_within_budget` reports this and is 100 by construction.

Observations are generated as `true_means[i] + true_sds[i] * Z`, where
\\Z\\ is standard normal (`distribution = "normal"`) or a standardized
right skewed chi square variable with 3 degrees of freedom (`"skewed"`;
skewness about 1.63).

The interval and the test are asymptotic, and the size of the shortfall
for a design of interest is what this function reports. As a guide, in
runs of 10,000 replications (standard error about 0.002) over five
designs with a pilot of 10 and about 15 to 120 observations per group,
coverage of the 95 percent interval was 0.936 to 0.949 for normal
populations and 0.925 to 0.941 for skewed ones, and the rejection rate
of the noninferiority test at the margin (nominal 0.05) was 0.050 to
0.061 for normal populations and 0.040 to 0.085 for skewed ones. The
default pilot of 5 does a little worse than that when the populations
are skewed and the budget is small, and a pilot of 2 does much worse;
[`seq_budget_pilot`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md)
has the figures, and `m_0` lets the comparison be made for the design at
hand.

`procedure = "first_passage"` is the package's reading of the article's
sequence, and differs from it in three details that the help page of
[`seq_budget_next`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
sets out; the main one is that every open group is judged after each
observation, where the article judges a population only after an
observation drawn from it.

Each step of a replication recomputes the standard deviation of the
group that received the observation, so the time for one replication
grows faster than the number of observations the budget buys: on one
core of the computer used to check the package, the default `G = 1000`
took about 6 seconds when a study ended with about 200 observations in
all, about 2 minutes with about 3,200, and about 11 minutes with about
12,800. For a budget that buys thousands of observations, start with a
small `G`. A budget that pays for more than one million observations in
the cheapest group (`budget / min(cost_per_group)`) is refused before
anything is simulated: at that limit, with costs of 38, 15, and 10, a
study ended with about 415,000 observations and one replication took
more than two minutes.

## References

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods*,
*30*(2), 425–439.
[doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

## See also

[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md),
[`seq_budget_next`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md),
[`seq_budget_optimal`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md),
[`decide_seq_budget_noninferiority`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md).
In DMAR,
[`ss_seq_c_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_seq_c_sensitivity.html)
is the corresponding study for the sequential fixed-width interval for a
contrast.

Other seq_budget:
[`decide_seq_budget`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md),
[`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md),
[`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md),
[`seq_budget_optimal()`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md),
[`seq_budget_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# A new training format (group 1) against the average of two existing
# ones, with unequal costs and a budget of 1,500. G = 100 keeps the
# example quick; a reported analysis deserves the default G = 1000 or
# more.
seq_budget_sensitivity(
  true_means = c(50, 50, 50),
  true_sds = c(9, 6, 6),
  contrast_weights = c(1, -0.5, -0.5),
  cost_per_group = c(38, 15, 10),
  budget = 1500,
  margin = 5,
  G = 100,
  seed = 113,
  m_0 = 10
)
#> term                           value
#> true_contrast                  0
#> budget                         1500
#> conf_level                     0.95
#> replications                   100
#> coverage                       0.92
#> mean_total_cost                1490.74
#> sd_total_cost                  8.192902
#> pct_within_budget              100
#> pct_rule_met                   0
#> mean_n_1                       28.28
#> mean_n_2                       15.22
#> mean_n_3                       18.78
#> optimal_n_1                    28.59535
#> optimal_n_2                    15.17121
#> optimal_n_3                    18.58086
#> mean_estimated_variance        3.707542
#> mean_true_variance             3.989052
#> empirical_variance             4.883586
#> optimal_variance               3.910227
#> margin                         5
#> alpha_level                    0.05
#> noninferiority_rejection_rate  0.75
```
