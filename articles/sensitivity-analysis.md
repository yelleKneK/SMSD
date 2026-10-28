# Sensitivity Analysis: Studying a Design Before Collecting Data

## Why Study a Design First

In a sequential procedure the final sample size is not chosen; it is
whatever the data make it. The same is true of the width of the final
confidence interval, which never exceeds the target but can fall below
it, and of whether that interval contains the population value. Theory
describes these only in the limit. As the target width narrows, the
final sample size approaches the one a planner who knew the population
values would have fixed in advance, and the coverage of the interval
approaches the confidence level (Kelley, Darku, & Chattopadhyay, 2018,
Lemma 3 and Theorem 2; 2019, Theorem 1). A particular study is not at
the limit. A sensitivity analysis shows what the procedure does for one
population value, one target, and one way of checking, by running the
whole procedure many times on simulated data whose population value is
known.

DMAR pairs most of its `ss_aipe_*()` fixed sample size planners with a
sensitivity function
([`ss_aipe_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd_sensitivity.html)
beside
[`ss_aipe_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.html),
for example), which shows what happens when a planning value is wrong. A
sequential procedure has no planning value. Its sensitivity analysis
shows what the stopping rule does in samples of the size the rule itself
produces.

## What Can Be Studied

``` r

grep("^seq_.*_sensitivity$", ls("package:SMSD"), value = TRUE)
#>  [1] "seq_aipe_c_sensitivity"        "seq_aipe_cv_sensitivity"      
#>  [3] "seq_aipe_gini_sensitivity"     "seq_aipe_r_sensitivity"       
#>  [5] "seq_aipe_R2_sensitivity"       "seq_aipe_reg_coef_sensitivity"
#>  [7] "seq_aipe_sm_sensitivity"       "seq_aipe_smd_sensitivity"     
#>  [9] "seq_budget_sensitivity"        "seq_mr_cv_sensitivity"        
#> [11] "seq_mr_smd_sensitivity"
```

Every sequential AIPE procedure has a sensitivity function, the contrast
and the squared multiple correlation included; the one for the
correlation studies the Pearson coefficient. The two minimum risk
procedures,
[`seq_mr_cv()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
and
[`seq_mr_smd()`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md),
have one each, and so does the budget constrained procedure. The
classical procedures for a mean have none; the section [A Population of
Your Own](#a-population-of-your-own) shows how to write the loop by
hand, which is also the way to study Kendall’s tau, Spearman’s rho, or a
population that the built in choices do not cover.

Each function does what a study would do: it draws the pilot sample,
checks the stopping rule, adds observations, and checks again, until the
rule is met. It does this `G` times and returns a summary.

## One Call

A training program is to be compared with a control group, and the
researcher wants a 95% confidence interval for the standardized mean
difference. What should a study with a target width of 0.80 expect if
the population value is 0.5?

``` r

smd <- seq_aipe_smd_sensitivity(true_value = 0.5, width = 0.80, G = 100,
                                seed = 113)
smd
#> term                   value
#> true_value             0.5
#> target_width           0.8
#> conf_level             0.95
#> replications           100
#> mean_n                 49.97
#> median_n               49
#> sd_n                   3.534691
#> mean_width             0.7912251
#> median_width           0.7916159
#> pct_width_within       100
#> coverage               0.96
#> pct_hit_max_n          0
#> pct_variance_unusable  0
```

The rows are the same for every `seq_aipe_*_sensitivity()` function.

- `true_value`, `target_width`, `conf_level`, `replications`: what was
  asked for.
- `mean_n`, `median_n`, `sd_n`: the final sample size, per group for
  this effect size. A study like this one ends with about 50 per group,
  give or take a few. With the population values known, the sample size
  for this width would be
  `4 * qnorm(0.975)^2 / 0.80^2 * (2 + 0.5^2 / 4)` = 49.5 per group for
  normal data, so the procedure spends almost nothing for not knowing
  them.
- `mean_width`, `median_width`, `pct_width_within`: the width of the
  final interval, and the percentage of replications in which it was no
  wider than the target. That percentage is 100 whenever the stopping
  rule ended every replication.
- `coverage`: the proportion of final intervals that contained the
  population value.
- `pct_hit_max_n`: the percentage of replications that the cap `max_n`
  (100,000 by default) stopped before the rule was met. When it is above
  zero the other rows describe a truncated procedure.
- `pct_variance_unusable`: the percentage of replications in which, at
  one or more checks, the variance estimate was zero, negative, or not
  finite, so that the rule was reported as not met there and sampling
  continued
  ([`vignette("seq-aipe-smd", package = "SMSD")`](https://yelleknek.github.io/SMSD/articles/seq-aipe-smd.md)
  shows such a check).

The result is an object of class `smsd_sensitivity`. It prints every
value in fixed notation, and it is an ordinary data frame with columns
`term` and `value`, so a number can be taken out by name. The one line
function `value_of()` defined here does that in the rest of this
vignette.

``` r

class(smd)
#> [1] "smsd_sensitivity" "data.frame"
smd$value[smd$term == "coverage"]
#> [1] 0.96

value_of <- function(out, term) out$value[out$term == term]
value_of(smd, "mean_n")
#> [1] 49.97
```

Every row is an estimate from `G` replications. The Monte Carlo standard
error of a coverage of *p* is the square root of *p*(1 − *p*) / *G*:

``` r

mc_se <- function(out) {
  p <- value_of(out, "coverage")
  G <- value_of(out, "replications")
  sqrt(p * (1 - p) / G)
}
mc_se(smd)
#> [1] 0.01959592
```

With `G = 100` the coverage above is known only to within about two
standard errors, 0.039, either way. The runs in this vignette use a
small `G` and wide targets so that it builds in seconds. Unless the text
gives a standard error and draws a conclusion from it, they show how to
read the output and are not evidence about the procedures. A result that
will be reported deserves `G = 10000`, at the width the study will use,
which puts the standard error of a coverage near 0.95 at about 0.002. At
a narrow target such a run can take from minutes to an hour or more, and
nothing is printed while it runs, so time a call with a small `G` first.
The time grows with the final sample size and falls as `batch` rises,
but a larger `batch` describes a study that checks less often ([Checking
in Batches](#checking-in-batches)), so use the batch size the study will
use. The `seed` argument makes a run reproducible and leaves the
caller’s random number stream as it found it.

Kelley, Darku, and Chattopadhyay (2018, pp. 232–233, Tables 1 and 2)
report simulations of this kind for the standardized mean difference,
with one observation per group added at each stage and 5,000
replications, under the heading “Characteristics of the final sample
size: An empirical demonstration.”

The other functions are called the same way, with two exceptions that
follow from what they study.
[`seq_aipe_R2_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_sensitivity.md)
also takes the number of predictors, `p`, and is used below.
[`seq_aipe_c_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
has no single population value to be given: it takes the population
means of the groups, the contrast weights, and the population standard
deviations (1 by default), and its `true_value` row is the contrast of
the means. Here a new training format is compared with the average of
two existing ones, the population contrast is half a standard deviation,
and the target width is 1:

``` r

contrast <- seq_aipe_c_sensitivity(true_means = c(0.5, 0, 0),
                                   contrast_weights = c(1, -0.5, -0.5),
                                   width = 1, G = 100, seed = 113)
contrast
#> term                   value
#> true_value             0.5
#> target_width           1
#> conf_level             0.95
#> replications           100
#> mean_n                 22.92
#> median_n               23
#> sd_n                   4.877655
#> mean_width             0.9654817
#> median_width           0.9682345
#> pct_width_within       100
#> coverage               0.94
#> pct_hit_max_n          0
#> pct_variance_unusable  0
```

The sample sizes in this summary are per group, and the width of a
contrast is in the units of the response. The coverage row again rests
on 100 replications; its standard error is 0.024.

## Checking in Batches

By default the stopping rule is checked after every observation (per
group, for a two group design), which is the purely sequential procedure
of the articles. Few studies can do that. The `batch` argument is the
number of observations added between checks, and a sensitivity analysis
shows what checking less often costs.

``` r

batches <- c(10, 20)
runs <- lapply(batches, function(b) {
  seq_aipe_smd_sensitivity(true_value = 0.5, width = 0.80, G = 100,
                           seed = 113, batch = b)
})
runs <- c(list(smd), runs)

overshoot <- sapply(runs, function(out) {
  c(mean_n = value_of(out, "mean_n"),
    mean_width = value_of(out, "mean_width"),
    coverage = value_of(out, "coverage"),
    se_coverage = mc_se(out))
})
colnames(overshoot) <- paste("batch", c(1, batches))
round(overshoot, 3)
#>             batch 1 batch 10 batch 20
#> mean_n       49.970   54.900   64.000
#> mean_width    0.791    0.758    0.701
#> coverage      0.960    0.930    0.940
#> se_coverage   0.020    0.026    0.024
```

A study that checks only after every batch cannot stop between checks,
so it ends with more observations and a somewhat narrower interval. Here
the mean final sample size per group rises from 50.0 with a check after
every observation to 64.0 with a check after every 20. The coverage row
cannot separate the three designs: its entries differ by at most 0.030,
and each has a Monte Carlo standard error of 0.020 to 0.026 (the last
row). In larger runs at a target width of 0.4 the mean final sample size
per group was 199.1, 201.2, 203.6, and 209.5 for batches of 1, 5, 10,
and 20, about half a batch more per group, and coverage stayed between
0.949 and 0.952 (4,700 to 6,000 replications each;
[`?seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)).

## Reading the Coverage Row

The width of the final interval is guaranteed. Its coverage is not: the
articles prove that coverage approaches the confidence level as the
target narrows, and say nothing about any particular width. Whether a
target is narrow enough depends on the effect size, on the population
value, and on the distribution of the data.

For a standardized mean difference of 0.5 and normal data, even the wide
target used above is narrow enough: with many more replications the
coverage at a width of 0.8 was 0.95
([`?seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md);
the row is in the table below). The same design with right skewed data
is a different matter. The argument `distribution = "skewed"` draws both
groups from a chi square distribution with 3 degrees of freedom, shifted
and scaled so that the means and the variance are those of the normal
case. This run checks after every five participants per group, which
keeps 1,000 replications quick:

``` r

skewed <- seq_aipe_smd_sensitivity(true_value = 0.5, width = 0.80,
                                   distribution = "skewed", G = 1000,
                                   seed = 113, batch = 5)
round(c(mean_n = value_of(skewed, "mean_n"),
        coverage = value_of(skewed, "coverage"),
        se_coverage = mc_se(skewed)), 3)
#>      mean_n    coverage se_coverage 
#>      52.285       0.925       0.008
```

The coverage is 0.925, which is 3.0 standard errors below the 0.95 that
was asked for, more than Monte Carlo error accounts for. The interval is
as narrow as requested, but it contains the population value somewhat
less often than its confidence level says.

At narrower targets the shortfall under skew is much smaller, though at
a width of 0.4 it has not vanished: coverage was 0.944 at a width of 0.4
and 0.949 at 0.2, against 0.952 and 0.946 for normal data
([`?seq_aipe_smd`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md)).
The verification study in
[`vignette("complete-tour", package = "SMSD")`](https://yelleknek.github.io/SMSD/articles/complete-tour.md),
with 10,000 replications for each of the population values 0.2, 0.5, and
0.8 at a width of 0.4, found 0.942 to 0.946 for skewed data and 0.949 to
0.951 for normal data: about half a percentage point short of 0.95 under
skew, which is more than Monte Carlo error.

In the other direction, coverage is lower again at a target wider still,
for normal data as well, and the more so the larger the population
value. A design at the wide end deserves a run of its own under each
population that is plausible.

For a correlation, a target half as wide as the one above may not be
narrow enough even with normal data. This study of a correlation of 0.4
checks after every five pairs:

``` r

r <- seq_aipe_r_sensitivity(true_value = 0.4, width = 0.4, G = 100,
                            seed = 113, batch = 5)
r
#> term                   value
#> true_value             0.4
#> target_width           0.4
#> conf_level             0.95
#> replications           100
#> mean_n                 62.25
#> median_n               60
#> sd_n                   22.64543
#> mean_width             0.3714632
#> median_width           0.3826196
#> pct_width_within       100
#> coverage               0.82
#> pct_hit_max_n          0
#> pct_variance_unusable  12
```

In this run 82 of the 100 final intervals contained the population
correlation, a coverage of 0.820 with a standard error of 0.038 where
0.95 was asked for, and in 12 percent of the replications the variance
estimate was unusable at some check. A run of this size cannot give the
second decimal place, but the shortfall is real: with one pair added at
a time and many more replications, the coverage at this population value
and width is 0.78 (the third row of the table below). The reason is the
estimate of the asymptotic variance in small samples. It is a
combination of ratios of sample moments that tends to run low when the
sample is small, so the rule is met too early and the interval is too
narrow (Kelley, Darku, & Chattopadhyay, 2019, Appendix D, p. 514).

The table collects what these functions measured for the standardized
mean difference and for some of the harder cases, at a 95% confidence
level with the default settings and normal data unless a distribution is
named. It is taken from
[`?seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md),
where the replication counts and standard errors are given (each row
rests on 800 to 34,000 replications).

| Effect size and population value | Target width | Mean final n | Coverage | Percent of replications with an unusable variance estimate |
|----|----|----|----|----|
| Standardized mean difference, 0.5 | 0.8 | 51 | 0.95 | below 0.1 |
| Standardized mean difference, 0.5 | 0.4 | 199 | 0.95 | 0.0 |
| Correlation, 0.4 | 0.4 | 55 | 0.78 | 16.5 |
| Correlation, 0.4 | 0.2 | 263 | 0.93 | 0.9 |
| Correlation, 0.4 | 0.1 | 1,079 | 0.95 | 0.0 |
| Correlation, 0.9 | 0.1 | 75 | 0.45 | 79.5 |
| Correlation, 0.9 | 0.025 | 819 | 0.90 | 1.3 |
| Regression slope, 0.5 | 0.4 | 80 | 0.83 | 20.1 |
| Regression slope, 0.5 | 0.2 | 377 | 0.94 | 1.3 |
| Regression slope, 2 | 0.2 | 217 | 0.61 | 65.3 |
| Coefficient of variation, 0.5, lognormal | 0.2 | 62 | 0.79 | 0.0 |
| Coefficient of variation, 0.5, lognormal | 0.05 | 1,561 | 0.90 | 0.0 |
| Gini index, 0.3 | 0.2 | 30 | 0.89 | 0.0 |

Three things in the table matter for planning.

First, coverage depends on the population value as well as on the width.
A width of 0.1 is ample for a correlation of 0.4 and far from enough for
a correlation of 0.9, because a strong correlation is estimated with so
little sampling error that the rule is met in a small sample. Run the
analysis at every population value that is plausible for the study, the
strongest included.

Second, the remedies are a narrower target, a larger minimum pilot
sample (`m_0`), or both. With `m_0 = 40` at a width of 0.4, coverage for
a correlation of 0.4 rose from 0.78 to 0.90 while the mean final sample
size rose from 55 to 64, and coverage for a slope of 0.5 rose from 0.83
to 0.92. A larger pilot is a weaker remedy when the relationship is
strong: for a correlation of 0.9 at a width of 0.1, `m_0 = 100` raised
coverage only from 0.45 to 0.69.

Third, an unusable variance estimate at an early check is not the
problem itself; the check functions report the rule as not met there and
sampling continues. It is a sign that the samples the rule is working
with are small for this effect size, and the estimate on which a
replication finally stops can still be too low. When the relationship is
very strong the estimate can also stay unusable for a long time, so set
`max_n` when exploring such conditions and read `pct_hit_max_n`.

## Several Conditions in One Table: The Squared Multiple Correlation

The summaries are data frames, so a table over conditions takes a few
lines. For the squared multiple correlation the condition that matters
most is the number of predictors. The function generates multivariate
normal data, the assumption under which that procedure is derived; `p`
is the number of predictors.

``` r

predictors <- c(2, 5, 10)
r2_runs <- lapply(predictors, function(p) {
  seq_aipe_R2_sensitivity(true_value = 0.50, width = 0.30, p = p, G = 200,
                          seed = 113)
})
data.frame(
  predictors = predictors,
  mean_n = sapply(r2_runs, value_of, "mean_n"),
  coverage = sapply(r2_runs, value_of, "coverage"),
  se_coverage = sapply(r2_runs, mc_se)
)
#>   predictors mean_n coverage se_coverage
#> 1          2  84.39    0.880  0.02297825
#> 2          5  81.97    0.870  0.02378024
#> 3         10  72.47    0.705  0.03224709
```

The stopping rule asks for no more observations when there are more
predictors, and the coverage falls: the sample squared multiple
correlation overestimates the population value by more the more
predictors there are, and the rule does not correct for it. The standard
errors in the last column are large because `G` is small, but the trend
is not a Monte Carlo accident. With 10,000 replications per condition
the coverage at this population value and width was 0.900 with 2
predictors, 0.845 with 5, 0.668 with 10, and 0.361 with 15
([`?seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md)).
[`vignette("seq-aipe-correlations", package = "SMSD")`](https://yelleknek.github.io/SMSD/articles/seq-aipe-correlations.md)
discusses what to do about it.

## Minimum Risk

The minimum risk procedures (Chattopadhyay & Kelley, 2016, 2017) have no
interval and no coverage. Their sensitivity functions compare the
procedure with the best fixed sample size, the one that would minimize
the risk if the population value were known. For a standardized mean
difference of 0.5, a price of 8,000 per unit of squared error, and a
cost of 4 per participant:

``` r

mr <- seq_mr_smd_sensitivity(true_value = 0.5, A = 8000, c = 4, G = 200,
                             seed = 113)
mr
#> term                value
#> true_value          0.5
#> loss_A              8000
#> cost_1              4
#> cost_2              4
#> replications        200
#> optimal_n1          45.41476
#> optimal_n2          45.41476
#> minimum_risk        726.6361
#> mean_n1             50.635
#> mean_n2             50.635
#> median_n1           51
#> median_n2           51
#> sd_n1               0.6958419
#> sd_n2               0.6958419
#> ratio_n1            1.114946
#> ratio_n2            1.114946
#> mean_estimate       0.5020693
#> mean_risk           726.8934
#> risk_efficiency     1.000354
#> se_risk_efficiency  0.04153699
#> regret              0.257265
#> pct_hit_max_n       0
```

The rows to read are `ratio_n1` and `ratio_n2`, the mean final sample
size of each group over the optimal one, and `risk_efficiency`, the mean
realized loss over the minimum risk, with its Monte Carlo standard error
in the next row. Each should be near 1. Here the procedure takes about
11 percent more than the optimal 45.4 per group, because the correction
term in its stopping rule guards against stopping too early, and its
risk efficiency is 1.000 (standard error 0.042): the risk is nearly flat
around its minimum, so the extra observations cost little.
[`seq_mr_cv_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)
reports the same quantities for the coefficient of variation, with one
row where the two group function has two: `cost`, `optimal_n`, `mean_n`,
`median_n`, `sd_n`, and `ratio_n`. Both functions take `batch`.

## A Budget

[`seq_budget_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)
runs the budget constrained procedure of Chattopadhyay, Bandyopadhyay,
Kelley, and Padalunkal (2025) for population means, standard deviations,
costs, and a budget that the user supplies. It adds one observation at a
time to the group that
[`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
names, so it has no `batch` argument. Here a new training format is
compared with the average of two existing ones; the population contrast
is 2, and a noninferiority margin of 5 is of interest.

``` r

budget <- seq_budget_sensitivity(
  true_means = c(52, 50, 50),
  true_sds = c(9, 6, 6),
  contrast_weights = c(1, -0.5, -0.5),
  cost_per_group = c(38, 15, 10),
  budget = 1500,
  margin = 5,
  G = 200,
  seed = 113,
  m_0 = 10
)
budget
#> term                           value
#> true_contrast                  2
#> budget                         1500
#> conf_level                     0.95
#> replications                   200
#> coverage                       0.93
#> mean_total_cost                1490.24
#> sd_total_cost                  8.400311
#> pct_within_budget              100
#> pct_rule_met                   0
#> mean_n_1                       28.23
#> mean_n_2                       15.22
#> mean_n_3                       18.92
#> optimal_n_1                    28.59535
#> optimal_n_2                    15.17121
#> optimal_n_3                    18.58086
#> mean_estimated_variance        3.721556
#> mean_true_variance             3.989493
#> empirical_variance             4.425726
#> optimal_variance               3.910227
#> margin                         5
#> alpha_level                    0.05
#> noninferiority_rejection_rate  0.975
```

The summary answers the questions a budget raises. How much of it is
spent: `mean_total_cost`, and `pct_within_budget`, which is 100 by
construction. How sampling ends: `pct_rule_met` is the percentage of
replications in which every group reached its criterion, and in the rest
the budget could not pay for another needed observation, which is the
usual ending. How close the allocation comes to the best one: the
`mean_n_` rows beside the `optimal_n_` rows, and `mean_true_variance`
beside `optimal_variance`. How well the interval and the test perform:
`coverage`, with a standard error here of 0.018, and
`noninferiority_rejection_rate`, which is the power of the
noninferiority test at this margin because the population contrast is
above the negative of the margin.

## A Population of Your Own

The sensitivity functions generate data from a few named populations.
When the population of interest is not among them, or the procedure has
no sensitivity function, the loop is short to write, because the check
functions do all of the work. The example is a standardized mean
difference between two groups rated on a 7 point scale: discrete,
bounded, and skewed in opposite directions in the two groups.

``` r

p_control <- c(0.10, 0.20, 0.25, 0.20, 0.12, 0.08, 0.05)
p_program <- rev(p_control)

rating_mean <- function(p) sum(1:7 * p)
rating_var <- function(p) sum((1:7 - rating_mean(p))^2 * p)
true_delta <- (rating_mean(p_program) - rating_mean(p_control)) /
  sqrt((rating_var(p_program) + rating_var(p_control)) / 2)
true_delta
#> [1] 0.6487847

one_study <- function(width = 0.80, batch = 10) {
  m <- seq_aipe_smd_pilot(width = width)
  program <- sample(7, m, replace = TRUE, prob = p_program)
  control <- sample(7, m, replace = TRUE, prob = p_control)
  repeat {
    out <- seq_aipe_smd(program, control, width = width)
    if (attr(out, "is_satisfied")) break
    program <- c(program, sample(7, batch, replace = TRUE, prob = p_program))
    control <- c(control, sample(7, batch, replace = TRUE, prob = p_control))
  }
  summary(out)
}

set.seed(113)
studies <- do.call(rbind, replicate(100, one_study(), simplify = FALSE))
covered <- with(studies, lower_limit <= true_delta & true_delta <= upper_limit)
c(mean_n = mean(studies$current_n), coverage = mean(covered))
#>   mean_n coverage 
#>    61.60     0.97
```

[`summary()`](https://rdrr.io/r/base/summary.html) turns each final
check into one row, so the replications stack into a data frame with one
row per simulated study. The coverage of 0.970 rests on 100 studies and
has a standard error of about 0.017. A loop for Kendall’s tau or
Spearman’s rho through
[`seq_aipe_r()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md)
has the same shape.

So does a loop for a classical procedure, none of which has a
sensitivity function. The purely sequential fixed width interval for a
normal mean,
[`seq_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md),
takes the half width `d`. With a half width of 0.5 and a population
standard deviation of 2, a planner who knew the standard deviation would
take (`qnorm(0.975) * 2 / 0.5`)² = 61.5 observations.

``` r

one_mean_study <- function(d = 0.5, mu = 10, sigma = 2) {
  x <- rnorm(seq_fwd_norm_mean_pilot(d = d), mean = mu, sd = sigma)
  repeat {
    out <- seq_fwd_norm_mean(x, d = d)
    if (attr(out, "is_satisfied")) break
    x <- c(x, rnorm(1, mean = mu, sd = sigma))
  }
  summary(out)
}

set.seed(113)
mean_studies <- do.call(rbind, replicate(100, one_mean_study(),
                                         simplify = FALSE))
mean_covered <- with(mean_studies, lower_limit <= 10 & 10 <= upper_limit)
c(mean_n = mean(mean_studies$current_n), coverage = mean(mean_covered))
#>   mean_n coverage 
#>    59.25     0.89
```

The 100 studies ended with 59.2 observations on average, and 0.890 of
their intervals contained the population mean (standard error 0.031).
The confidence level of this procedure is asymptotic, and a shortfall at
a sample size of this order is real: for the same design
[`?seq_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md)
reports a coverage of 0.938 from 40,000 replications where 0.95 was
asked for. A run of this size cannot measure a shortfall of one
percentage point; it shows the shape of the loop. The help pages of the
classical procedures report runs of this kind with tens of thousands of
replications, and
[`vignette("special-cases-and-simulation", package = "SMSD")`](https://yelleknek.github.io/SMSD/articles/special-cases-and-simulation.md)
collects them.

## Practical Guidance

- Run the analysis at the width, confidence level, pilot size, and batch
  size the study will use, at every population value that is plausible,
  and under a skewed population as well as a normal one where the
  function offers the choice.
- Use `G = 10000` for anything that will be reported, give the seed, and
  report the Monte Carlo standard error with the coverage. Time the call
  with a small `G` first: at a narrow target ten thousand replications
  can take an hour or more.
- Read `coverage` before relying on an interval at a planned width, and
  `pct_hit_max_n` before reading anything else when `max_n` was set.
- When coverage is short of the confidence level, plan a narrower
  target, raise `m_0`, or both, and run the analysis again.
- A preregistration or a proposal is the place for the resulting table:
  it tells a reviewer how large the study is likely to become and what
  the interval at the end can be trusted to do.

## References

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods, 30*(2),
425–439. [doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research, 51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
difference with minimum risk: Maximizing accuracy and minimizing cost
with sequential estimation. *Psychological Methods, 22*(1), 94–113.
[doi:10.1037/met0000089](https://doi.org/10.1037/met0000089)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods, 23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
accuracy in parameter estimation for population correlation
coefficients. *Psychological Methods, 24*(4), 492–515.
[doi:10.1037/met0000203](https://doi.org/10.1037/met0000203)
