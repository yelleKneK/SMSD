# Sequential Check for a Budget Constrained K-Population Contrast

Given the data collected so far from each of *K* populations, the cost
of an observation in each, and a total budget, reports whether sampling
should continue under the sequential procedure of Chattopadhyay,
Bandyopadhyay, Kelley, and Padalunkal (2025), and returns the current
estimate of the contrast with its confidence interval. Use it when a
contrast of means has to be estimated as precisely as a fixed budget
allows and the population variances that determine the best allocation
are unknown.
[`seq_budget_next`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
says which group the next observation should come from.

## Usage

``` r
seq_budget(
  groups,
  contrast_weights,
  cost_per_group,
  budget,
  conf_level = 0.95,
  closed = NULL
)
```

## Arguments

- groups:

  A list of *K* (at least 2) numeric vectors, one per population,
  holding the observations collected so far. Each needs at least 2
  observations, so that a sample variance exists, and no missing or
  infinite values. The groups may have different sizes.

- contrast_weights:

  Numeric vector of length *K* giving the contrast coefficients \\c_i\\,
  in the order of `groups`, not all zero. Code the contrast so that
  larger values are better if a noninferiority or superiority decision
  is planned. The article takes the weights to sum to zero; the sum is
  not checked.

- cost_per_group:

  Numeric vector of length *K* giving the cost \\a_i\\ of one
  observation from each population, each positive, in the same units as
  `budget` and in the order of `groups`.

- budget:

  The total budget \\A_0\\ for all observations, the pilot included: a
  single positive number.

- conf_level:

  The confidence level of the two-sided interval (default 0.95). It
  plays no part in the decision functions, which take a one-sided level
  of their own.

- closed:

  Optional; leave `NULL` (the default) to judge every group on the
  current data. For the first passage rule, a logical vector of length
  *K* that is `TRUE` for each group from which sampling has already
  stopped: the `met_per_group` attribute of the previous check. A closed
  group stays closed, unless its variance is unusable.

## Value

An object of class `c("smsd_check", "data.frame")` with columns `term`
and `value` and six rows, in this order: `lower_limit`, `contrast` (the
estimate \\\hat{\psi}\\), `upper_limit`, `variance` (the estimated
variance \\\hat{V}\\ of the contrast estimator), `current_total_n`, and
`current_total_cost`. A group whose variance is unusable adds nothing to
the variance row, so the two limits are then too close together and are
not an interval to report (see Details). Its attributes are:

- `method`:

  the label of the procedure,
  `"budget constrained K-population contrast"`. The next observation and
  decision functions recognize the object by it.

- `is_satisfied`, `rule_met`, `budget_exhausted`:

  the three stopping attributes described under Details.

- `budget_exceeded`:

  `TRUE` when the data already cost more than the budget.

- `variance_unusable`:

  `TRUE` when the sample variance of at least one group with a nonzero
  weight is unusable (zero, or computed from observations that agree to
  about ten significant digits); `FALSE` otherwise.

- `current_n`, `criterion_n`:

  the totals over groups of the current sample sizes and of the
  criteria.

- `conf_level`:

  the confidence level of the interval.

- `notes`:

  a character vector with the notes that
  [`print()`](https://rdrr.io/r/base/print.html) shows (empty when there
  are none).

- `n_per_group`, `mean_per_group`, `sd_per_group`:

  the sample size, mean, and standard deviation of each group, in the
  order of `groups`.

- `target_per_group`:

  the right hand side of Equation 18 for each group, not rounded.

- `criterion_per_group`:

  the smallest whole number at or above each target, a target within one
  part in 1e10 above a whole number counting as that number (an integer
  vector); for a group whose variance is unusable, its current size plus
  one.

- `met_per_group`:

  `TRUE` for each group that has met its criterion or, under the first
  passage rule, was closed earlier.

- `not_estimable_per_group`:

  `TRUE` for each group with a nonzero weight whose variance is
  unusable.

- `cost_per_group`, `contrast_weights`, `budget`:

  the inputs.

- `remaining_budget`:

  the budget minus the current total cost; negative when the budget is
  exceeded.

## Details

**The criterion.** The contrast \\\psi = \sum\_{i=1}^{K} c_i \mu_i\\ is
estimated by \\\hat{\psi} = \sum\_{i=1}^{K} c_i \bar{X}\_i\\, whose
variance is estimated by \\\hat{V} = \sum\_{i=1}^{K} c_i^2 s_i^2 / N_i\\
(Equation 5 of the article with sample variances). The variance is
smallest, for a budget \\A_0\\ and unit costs \\a_i\\, at the allocation
of
[`seq_budget_optimal`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md),
which depends on the unknown standard deviations. The sequential
procedure substitutes the current sample standard deviations: group
\\i\\ has met its criterion when (Equation 18) \$\$N_i \ge \frac{A_0 \\
\|c_i\| \\ s_i} {\sqrt{a_i} \\ \sum\_{l=1}^{K} \|c_l\| \\ s_l \\
\sqrt{a_l}}.\$\$ The right hand side is stored unrounded as
`target_per_group`; the smallest whole number that satisfies the
inequality is `criterion_per_group`. As for every criterion sample size
in the package, a target that exceeds a whole number by no more than one
part in 1e10 is taken to be that number, so that representation error in
a target that is a whole number does not add an observation. The
interval is \\\hat{\psi} \pm z \sqrt{\hat{V}}\\ with the standard normal
quantile \\z\\ for `conf_level`.

**How close the coverage is to nominal.** The coverage of the interval
is asymptotic. In samples of realistic size it falls short of nominal,
because the estimated variance runs below the true one, and the
shortfall depends on the pilot size. In simulations with a pilot of 10
per group (the article's simulation setting) and a budget that bought
about 40 to 100 observations per group, the 95 percent interval covered
the contrast in 0.942 to 0.949 of studies for normal populations and
0.936 to 0.941 for skewed ones; with 15 to 30 per group the figures were
0.936 and 0.925. A pilot of 5 gave about the same with the larger
budgets and somewhat less with the small one (0.931 and 0.907 in another
set of runs), and the minimum pilot of 2 gave markedly less;
[`seq_budget_pilot`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md)
has those figures.
[`seq_budget_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)
measures the coverage for a design of interest.

**When sampling ends.** Sampling ends when (a) every group has met its
criterion or (b) the budget cannot pay for another needed observation
(the article, p. 430). Three attributes state which:

- `rule_met`:

  `TRUE` when every group with a nonzero weight has met its criterion.

- `budget_exhausted`:

  `TRUE` when at least one group is still below its criterion and what
  is left of the budget cannot pay for one more observation from any
  such group. Only the groups that still need observations count: a
  cheaper group that has already met its criterion does not keep the
  procedure going.

- `is_satisfied`:

  `TRUE` when sampling has ended, by (a) or by (b), and the variance of
  every group could be used. It is `FALSE` while sampling should
  continue. It is also `FALSE` whenever the variance of a group is
  unusable, which the attribute `variance_unusable` reports (see below).

The function
[`seq_budget_next`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
has no group to name, and returns `NA`, exactly when one of the first
two attributes is `TRUE`; a loop that adds observations should end on
that `NA`. Unless a variance is unusable, that is also exactly when
`is_satisfied` is `TRUE`.

The targets of Equation 18 always cost exactly the budget in total
(multiply each by its \\a_i\\ and add). Judged on the current data,
therefore, every group can be at or above its criterion without passing
the budget only when each sample size equals its target exactly. So with
the default `closed = NULL` the procedure nearly always ends by (b),
with less of the budget left than the next needed observation costs and
one or more groups still below their criteria. The printed object then
carries a note that gives the amount left and names those groups, and
its stopping line reads
`Stopping rule met: FALSE (sampling has ended: the budget is exhausted)`,
because the rule itself (every group at its criterion) was not met.
`is_satisfied` is `TRUE` all the same: it says that sampling has ended.
A fourth attribute, `budget_exceeded`, is `TRUE` when the data supplied
already cost more than the budget (for example a pilot the budget cannot
pay for); the result is then flagged with a note, because the budget
constraint does not hold for it.

**Two readings of "has met its criterion".** In the article a population
stops being sampled the first time it meets its criterion and is never
sampled again (Stage 1 and Stage 2, pp. 429 to 430), and the open
populations are visited one at a time, in rounds. A check sees only the
current data, so by default it judges every group afresh: a group that
met its criterion earlier and has since fallen below it, because the
other standard deviations moved, is sampled again. This default is the
package's choice, not the article's procedure. In simulations of three
designs from Table 1 of the article (budget 5,000, pilot 10, normal and
skewed populations) it left about 10 to 12 of the budget unspent where
the first passage rule left 19 to 42, and the true variance of the
estimator at the final sample sizes was 0.2 to 1.0 percent smaller;
coverage differed by at most 0.006, in either direction. For the first
passage rule, pass the `met_per_group` attribute of the previous check
as `closed`, and ask
[`seq_budget_next`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
for `order = "rounds"`; the help page of that function lists the three
respects in which the result still differs from the article's procedure.
Under the first passage rule `rule_met` can become `TRUE` with part of
the budget unspent.

**A group whose variance is unusable.** When the sample variance of a
group with a nonzero weight is zero, no variance has been estimated for
that group. As everywhere in the package when a variance estimate is
unusable, the stopping rule is then reported as not met (`rule_met` and
`is_satisfied` are `FALSE`), `attr(out, "variance_unusable")` is `TRUE`,
and a note says so; see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md).
This happens when every observation in the group is the same, which a
small pilot and coarsely scored data can produce, or the same to about
ten significant digits, so that whatever variance is computed is
rounding error. (Scores so large that a sample variance overflows are
not this case: the function stops with a message that says to rescale
them.) The group's share of the budget cannot be estimated and the
variance of the contrast is understated; its target is reported as zero
and its criterion as one more than its current size. The limits are
still computed, from the understated variance, so they are too close
together (when every group is constant both equal the estimate) and are
not an interval to report.
[`summary()`](https://rdrr.io/r/base/summary.html), `glance()`, and
`tidy()` return them without the note, so test
`attr(out, "variance_unusable")` before using the limits. This is a
deliberate departure from the stopping rule as printed in the article,
under which such a group has met its criterion (the right hand side of
Equation 18 is zero) and is never sampled again, so that a pilot of
equal scores would fix that group's variance at zero for the rest of the
study.

[`seq_budget_next`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
names such a group whenever no group whose variance can be used has
fewer observations than it has, so it is sampled at the pace of the
least sampled of the other groups until its observations differ. When
the equal scores were an accident of a small pilot, one or two further
observations usually settle the matter. A group that really is constant,
or nearly so (a ceiling or a floor on the measure), keeps that pace for
as long as its observations stay equal, which can be to the end of the
study. It then receives about as many observations as the least sampled
of the other groups, far more than its share when it is expensive, and
the study ends with a variance that is still unusable. In a simulated
example in which the most expensive of three groups scored at its
ceiling 98 percent of the time, that group took 49 percent of the budget
where its optimal share was 14 percent, the variance of the estimator
was 1.6 times the optimum, and 61 percent of the studies ended with the
group's variance still unusable. If a group shows no variation after
more than a few observations beyond the pilot, stop and reconsider the
measure before spending the rest of the budget.

Sampling still ends when the budget cannot pay for another observation:
`budget_exhausted` is `TRUE` and
[`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
returns `NA`. `is_satisfied` stays `FALSE`, the printed notes say why,
and the decision functions refuse the check, because the interval rests
on a variance of zero for at least one group. (For large scores that
differ only far down their digits, subtracting the same constant from
every score restores a usable variance; it leaves a contrast whose
weights sum to zero unchanged.)

**Names.** `groups`, `contrast_weights`, and `cost_per_group` are
matched by position. If more than one of them is named, the names must
be the same and in the same order; otherwise the function stops, so that
weights or costs named in another order (for example against a list made
by [`split()`](https://rdrr.io/r/base/split.html), which is in
alphabetical order) are not applied to the wrong groups.

## References

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods*,
*30*(2), 425–439.
[doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

Luh, W. M., & Guo, J. H. (2016). Sample size planning for the
noninferiority or equivalence of a linear contrast with cost
considerations. *Psychological Methods*, *21*(1), 13–34.
[doi:10.1037/met0000039](https://doi.org/10.1037/met0000039)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`seq_budget_next`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
for the next observation,
[`seq_budget_pilot`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md)
and
[`seq_budget_optimal`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md)
for planning,
[`decide_seq_budget`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md)
for the tests, and
[`seq_budget_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)
for a simulation of a planned design. The vignette is
[`vignette("seq-budget", package = "SMSD")`](https://yelleknek.github.io/SMSD/articles/seq-budget.md).

In DMAR,
[`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.html)
is the sequential fixed-width interval for a contrast, and
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.html)
and
[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.html)
are the fixed sample size cousins.

Other seq_budget:
[`decide_seq_budget`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md),
[`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md),
[`seq_budget_optimal()`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md),
[`seq_budget_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md),
[`seq_budget_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Scores after a new training format (group 1) and two existing formats.
# An observation costs 38, 15, and 10, and the budget is 1,500.
set.seed(113)
weights <- c(1, -0.5, -0.5)
costs <- c(38, 15, 10)
scores <- list(
  rnorm(10, mean = 52, sd = 9),
  rnorm(10, mean = 50, sd = 6),
  rnorm(10, mean = 50, sd = 6)
)
out <- seq_budget(scores, weights, cost_per_group = costs, budget = 1500)
out
#> Sequential check (budget constrained K-population contrast)
#> 
#>   Current n       : 30
#>   Criterion n     : 65
#>   Confidence level: 95%
#> 
#>   lower_limit       : -9.311685
#>   contrast          : -3.364913
#>   upper_limit       : 2.581859
#>   variance          : 9.205903
#>   current_total_n   : 30
#>   current_total_cost: 630
#> 
#>   Stopping rule met: FALSE
attr(out, "criterion_per_group")
#> [1] 29 14 22
attr(out, "remaining_budget")
#> [1] 870

# Add one observation at a time, to the group seq_budget_next() names,
# until it has none to name. The cost never passes the budget.
true_means <- c(52, 50, 50)
true_sds <- c(9, 6, 6)
repeat {
  out <- seq_budget(scores, weights, cost_per_group = costs, budget = 1500)
  g <- seq_budget_next(out)
  if (is.na(g)) break
  scores[[g]] <- c(scores[[g]], rnorm(1, true_means[g], true_sds[g]))
}
out
#> Sequential check (budget constrained K-population contrast)
#> 
#>   Current n       : 65
#>   Criterion n     : 67
#>   Confidence level: 95%
#> 
#>   lower_limit       : -3.175372
#>   contrast          : 0.5710176
#>   upper_limit       : 4.317407
#>   variance          : 3.653673
#>   current_total_n   : 65
#>   current_total_cost: 1499
#> 
#>   Stopping rule met: FALSE (sampling has ended: the budget is exhausted)
#> 
#>   Note: Sampling has ended because the budget cannot pay for another
#>         observation: 1 of 1500 is left, and one more observation from a
#>         group still below its criterion (groups 2 and 3) costs at least
#>         10.
attr(out, "n_per_group")
#> [1] 28 13 24
```
