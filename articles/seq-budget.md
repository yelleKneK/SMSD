# Budget Constrained Sequential Procedures for K Populations

## The Research Question

A large organization trains new employees face to face at a central
site. Two cheaper formats are proposed: one-on-one video calls, and a
hybrid of recorded material followed by a shorter call. Face to face
training costs far more per trainee than either alternative, because of
travel and lodging. The organization will spend a fixed amount on a
study and wants to know whether the average of the two alternatives is
no worse than the face to face format by more than a margin that
matters, or is equivalent to it within such a margin. This is the
application that motivates Chattopadhyay, Bandyopadhyay, Kelley, and
Padalunkal (2025), and the same structure arises whenever several groups
with unequal sampling costs are compared under a fixed budget: multisite
program evaluations in education, marketing experiments across regions,
surveys that compare strata.

With μ₁ the mean score after face to face training and μ₂ and μ₃ the
means after the two alternatives, the quantity of interest is the
contrast (Equation 1 of the article)

ψ = −μ₁ + ½ μ₂ + ½ μ₃

How should the budget be divided among the three groups so that ψ is
estimated as precisely as possible? The answer depends on the population
variances, which are unknown before the study. A fixed sample size plan
has to take them from a guess or from a separate pilot study (the
approach of Luh & Guo, 2016). The sequential procedure learns them from
the data as the data arrive, in the purely sequential tradition surveyed
by Mukhopadhyay and de Silva (2009), and it counts the pilot
observations against the same budget and uses them in the analysis.

## The Method

Take a contrast ψ of *K* population means with weights *c*_(*i*), the
sum over the populations of *c*_(*i*) μ_(*i*). Its estimate is the same
weighted sum of the sample means, and the variance of the estimate is
the sum of *c*_(*i*)² σ_(*i*)² / *n*_(*i*), where σ_(*i*) is the
standard deviation of population *i* and *n*_(*i*) its sample size. An
observation from population *i* costs *a*_(*i*), and the total budget is
*A*₀. Theorem 1 of the article gives the sample sizes that minimize the
variance when the total cost, the sum of *a*_(*i*)*n*_(*i*), equals the
budget (its Equation 12),

*n*_(*i*0) = *A*₀ \|*c*_(*i*)\| σ_(*i*) / (√*a*_(*i*) × Σ_(*l*)
\|*c*_(*l*)\| σ_(*l*) √*a*_(*l*))

and the minimized variance (Equation 13),

*V*_(*o*) = (Σ_(*i*) \|*c*_(*i*)\| σ_(*i*) √*a*_(*i*))² / *A*₀

A group receives more observations when its scores vary more, when its
weight in the contrast is larger, and when it is cheaper to sample.
Minimizing the variance is what maximizes the power of the tests
described below.

Because the σ_(*i*) are unknown, the sequential procedure replaces them
with the current sample standard deviations *s*_(*i*). Group *i* has met
its criterion when its sample size *N*_(*i*) satisfies (Equation 18)

*N*_(*i*) ≥ *A*₀ \|*c*_(*i*)\| *s*_(*i*) / (√*a*_(*i*) × Σ_(*l*)
\|*c*_(*l*)\| *s*_(*l*) √*a*_(*l*))

and sampling ends when every group has met its criterion or when the
budget cannot pay for another needed observation. The functions divide
the work as follows.

| Step | Function |
|----|----|
| Plan with assumed standard deviations | [`seq_budget_optimal()`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md) |
| Choose the pilot size | [`seq_budget_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md) |
| Check the data in hand | [`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md) |
| Choose the group for the next observation | [`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md) |
| Decide noninferiority, superiority, inferiority, or equivalence | `decide_seq_budget_*()` |
| Study the design by simulation before collecting data | [`seq_budget_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md) |

## A Priori Planning

The article’s own illustration (p. 430) has three groups, a cost of 30
per observation in each, and a budget of 10,000. When the products of
weight and standard deviation, \|*c*_(*i*)\| σ_(*i*), are equal, the
optimal size is the same in every group:

``` r

plan <- seq_budget_optimal(
  contrast_weights = c(-1, 0.5, 0.5),
  sigmas           = c(1, 2, 2),
  cost_per_group   = c(30, 30, 30),
  budget           = 10000
)
plan$n_per_group
#> [1] 111.1111 111.1111 111.1111
plan$minimum_variance
#> [1] 0.027
```

The sizes are not whole numbers, and they cannot simply be rounded up:
112 per group would cost 10,080. The article rounds to 111 per group,
333 observations in all, at a cost of 9,990. A few units of the budget
go unspent whenever the unit costs do not divide what is left, and this
is how the sequential procedure ends as well.

Costs are rarely equal. Suppose face to face training costs 500 per
trainee (the figure the article uses), a video call 150, and the hybrid
format 250, with a budget of 50,000 (again the article’s figure) and
scores on a 30 point assessment whose standard deviations are thought to
be about 3, 3.5, and 2.5:

``` r

weights <- c(-1, 0.5, 0.5)
costs   <- c(500, 150, 250)
budget  <- 50000

plan <- seq_budget_optimal(
  contrast_weights = weights,
  sigmas           = c(3, 3.5, 2.5),
  cost_per_group   = costs,
  budget           = budget
)
plan$n_per_group
#> [1] 61.95278 65.98070 36.50602
sum(costs * plan$n_per_group)
#> [1] 50000
sqrt(plan$minimum_variance)
#> [1] 0.4842398
```

The face to face group carries twice the weight of either alternative,
yet it does not receive twice the observations, because each of its
observations costs two to more than three times as much. The last line
is the smallest standard error of the contrast estimate that this budget
can buy if the assumed standard deviations are right.

## The Sequential Workflow

The standard deviations above were guesses. The sequential procedure
does not need them. It starts with a pilot from every group, which is
paid for from the same budget. The article requires at least two
observations per group and advises more, because variance estimates from
very small samples are unstable (p. 431); its simulations used 10.

``` r

m_0 <- seq_budget_pilot(m_0 = 10)
m_0 * sum(costs)
#> [1] 9000
```

The pilot costs 9,000 of the 50,000. To have data to work with, the
scores are simulated here from populations whose true contrast is zero
(the alternatives are exactly as effective as face to face training on
average); in a study each new score would come from the next trainee.

``` r

set.seed(113)
true_means <- c(15, 15, 15)
true_sds   <- c(3, 3.5, 2.5)
new_score  <- function(g, k = 1) rnorm(k, true_means[g], true_sds[g])

scores <- lapply(1:3, new_score, k = m_0)

out <- seq_budget(scores, weights, cost_per_group = costs, budget = budget)
out
#> Sequential check (budget constrained K-population contrast)
#> 
#>   Current n       : 30
#>   Criterion n     : 164
#>   Confidence level: 95%
#> 
#>   lower_limit       : -0.01524134
#>   contrast          : 2.176975
#>   upper_limit       : 4.369191
#>   variance          : 1.251038
#>   current_total_n   : 30
#>   current_total_cost: 9000
#> 
#>   Stopping rule met: FALSE
attr(out, "n_per_group")
#> [1] 10 10 10
attr(out, "criterion_per_group")
#> [1] 61 60 43
attr(out, "remaining_budget")
#> [1] 41000
seq_budget_next(out)
#> [1] 1
```

The weights and costs are matched to the groups by position. If they are
named and the list of groups is named too, the names must agree in
order;
[`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
stops otherwise, so that, for example, weights typed in the order of the
design are not applied to a list that
[`split()`](https://rdrr.io/r/base/split.html) has put in alphabetical
order.

After the pilot every group is below its criterion, so sampling
continues.
[`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
names the group that should receive the next observation: by default the
one furthest below its estimated optimal size, among the groups whose
next observation the remaining budget can pay for. Because it never
names a group that cannot be paid for, a study that follows it one
observation at a time cannot overspend. The whole study is then a short
loop (the helper `term()` pulls one row out of a check, here to record
the running cost):

``` r

term <- function(check, name) check$value[check$term == name]

cost_path <- term(out, "current_total_cost")
repeat {
  out <- seq_budget(scores, weights, cost_per_group = costs, budget = budget)
  g <- seq_budget_next(out)
  if (is.na(g)) break
  scores[[g]] <- c(scores[[g]], new_score(g))
  cost_path <- c(cost_path, term(out, "current_total_cost") + costs[g])
}
out
#> Sequential check (budget constrained K-population contrast)
#> 
#>   Current n       : 164
#>   Criterion n     : 166
#>   Confidence level: 95%
#> 
#>   lower_limit       : -0.1118478
#>   contrast          : 0.8450936
#>   upper_limit       : 1.802035
#>   variance          : 0.2383826
#>   current_total_n   : 164
#>   current_total_cost: 50000
#> 
#>   Stopping rule met: FALSE (sampling has ended: the budget is exhausted)
#> 
#>   Note: Sampling has ended because the budget cannot pay for another
#>         observation: 0 of 50000 is left, and one more observation from
#>         a group still below its criterion (groups 1 and 2) costs at
#>         least 150.
```

``` r

attr(out, "n_per_group")
#> [1] 62 65 37
attr(out, "criterion_per_group")
#> [1] 63 66 37
max(cost_path) <= budget
#> [1] TRUE
```

``` r

plot(seq_along(cost_path) - 1, cost_path, type = "s",
     xlab = "Observations added after the pilot",
     ylab = "Total cost so far", ylim = c(0, budget))
abline(h = budget, lty = 2)
```

![Running total cost of the study against the number of observations
added after the pilot, rising in steps to the budget line without
crossing it.](seq-budget_files/figure-html/cost-path-1.png)

### Why Sampling Ends With Sizes Below the Criteria

Sampling has ended (`is_satisfied` is `TRUE`) although the current total
sample size (164) is below the total of the criteria (166). The note
under the printed check gives the reason: sampling ended because the
budget could not pay for another observation that was still needed, not
because every group reached its criterion. The attributes separate the
two ways of ending:

``` r

attr(out, "rule_met")
#> [1] FALSE
attr(out, "budget_exhausted")
#> [1] TRUE
attr(out, "is_satisfied")
#> [1] TRUE
attr(out, "remaining_budget")
#> [1] 0
```

`rule_met` is `TRUE` when every group has met its criterion.
`budget_exhausted` is `TRUE` when at least one group is still below its
criterion and what is left of the budget cannot pay for one more
observation from any such group.
[`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
returns `NA` exactly when one of the two holds, and `is_satisfied` is
`TRUE` then too, with one exception that is taken up below: a study in
which the variance of some group is unusable is never reported as
satisfied.
[`seq_report()`](https://yelleknek.github.io/SMSD/reference/seq_report.md)
words the result for a methods section and names the ending for what it
is:

``` r

seq_report(out)
#> Sample size was determined with the sequential budget constrained
#> procedure for a linear contrast of population means (Chattopadhyay,
#> Bandyopadhyay, Kelley, & Padalunkal, 2025), with a 95% confidence
#> interval and a total budget of 50,000. Sampling ended at a total sample
#> size of N = 164 because the budget could not pay for another
#> observation from a group still below its criterion sample size: 0 of
#> the 50,000 was left (the criterion sample sizes at that point came to
#> 166 in total). At that sample size the estimate of the contrast was
#> 0.845, 95% CI [-0.112, 1.802].
#> 
#> Reference
#> Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
#>     (2025). A sequential approach for noninferiority or equivalence of
#>     a linear contrast under cost constraints. Psychological Methods,
#>     30(2), 425–439. doi:10.1037/met0000570
```

Ending by the budget is the normal case, not a failure. Multiply each
right hand side of Equation 18 by its unit cost and add them up: the sum
is *A*₀, whatever the standard deviations are. Judged on the current
data, therefore, every group can be at or above its criterion without
passing the budget only when each sample size equals its target exactly,
which requires every target to be a whole number. Otherwise the
procedure runs until what is left of the budget (here 0) is less than
the cost of any observation still needed, just as 333 observations at 30
each leave 10 of 10,000 unspent in the article’s illustration.

Only the groups that still need observations count. At the end of this
run groups 1 and 2 were still below the criterion, and a cheaper group
that had already met its criterion would not have kept the study going.
Conversely, when the most expensive group can no longer be paid for but
a cheaper one still needs observations, the study continues with the
cheaper one.

If the data handed to
[`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
already cost more than the budget (a pilot the budget cannot pay for, or
a loop that added observations without asking
[`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)),
the check says so in a note and sets the attribute `budget_exceeded`:

``` r

over <- seq_budget(scores, weights, cost_per_group = costs, budget = 40000)
attr(over, "budget_exceeded")
#> [1] TRUE
writeLines(strwrap(attr(over, "notes"), width = 72))
#> The data already cost 50000, which is more than the budget of 40000, so
#> the budget constraint does not hold for this result. Raise the budget
#> or use a smaller pilot; once sampling is under way, choosing each new
#> observation with seq_budget_next() keeps the cost within the budget.
```

### When Observations Arrive in Batches

[`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
names one group for one observation. A study that collects several
observations at a time takes them from the group it names, but no more
than the remaining budget pays for and no more than the group is still
short of its criterion, and then checks again. The loop below does that
for batches of up to 10, starting from the same pilot scores as the loop
above.

``` r

batch <- 10

set.seed(113)
batched <- lapply(1:3, new_score, k = m_0)
repeat {
  out_b <- seq_budget(batched, weights, cost_per_group = costs,
                      budget = budget)
  g <- seq_budget_next(out_b)
  if (is.na(g)) break
  affordable <- floor(attr(out_b, "remaining_budget") / costs[g])
  short <- ceiling(attr(out_b, "criterion_per_group")[g] -
                     attr(out_b, "n_per_group")[g])
  batched[[g]] <- c(batched[[g]],
                    new_score(g, min(batch, affordable, max(1, short))))
}
attr(out_b, "n_per_group")
#> [1] 66 60 32
term(out_b, "current_total_cost")
#> [1] 50000
attr(out_b, "budget_exceeded")
#> [1] FALSE
```

The study ends within the budget, with 66, 60, 32 observations in the
three groups where the loop that checked after every observation ended
with 62, 65, 37. A loop that adds a whole batch whatever the remaining
budget is can overspend; the check then sets `budget_exceeded`, as in
the example above, and the budget constraint no longer holds for its
result.

### When a Group Shows No Variation

With a small pilot and coarsely scored data, every observation in a
group can be the same. Its sample variance is then zero, Equation 18
gives it a target of zero, and the stopping rule as printed in the
article would end the sampling from that group with a variance estimate
of zero. A variance of zero is not an estimate, and
[`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
treats it as every check in the package treats an unusable variance,
which is a deliberate departure from the article: the rule is reported
as not met, the attribute `variance_unusable` is `TRUE`, and the printed
notes say so. (The same holds when the scores of a group agree to about
ten significant digits, so that whatever variance is computed is
rounding error.)

``` r

ratings <- list(c(4, 4, 4), c(3, 5, 4), c(2, 4, 5))
flat <- seq_budget(ratings, weights, cost_per_group = costs, budget = budget)
flat
#> Sequential check (budget constrained K-population contrast)
#> 
#>   Current n       : 9
#>   Criterion n     : 250
#>   Confidence level: 95%
#> 
#>   lower_limit       : -1.199658
#>   contrast          : -0.1666667
#>   upper_limit       : 0.8663251
#>   variance          : 0.2777778
#>   current_total_n   : 9
#>   current_total_cost: 2700
#> 
#>   Stopping rule met: FALSE
#> 
#>   Note: The variance estimate is not yet usable at this sample size (it
#>         is zero, negative, or not finite), so the stopping rule is
#>         reported as not met. Collect more data and check again.
#> 
#>   Note: The observations in group 1 are all equal, or differ too little
#>         relative to their size for a variance to be computed, so its
#>         share of the budget cannot be estimated yet. seq_budget_next()
#>         names such a group when no group whose variance can be
#>         estimated has fewer observations.
attr(flat, "variance_unusable")
#> [1] TRUE
seq_budget_next(flat)
#> [1] 1
```

[`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
names such a group whenever no group whose variance can be used has
fewer observations than it has, so the group is sampled at the pace of
the least sampled of the others until its scores differ; one or two more
observations usually settle it. A group that really is constant (a
ceiling or a floor on the measure) keeps that pace to the end, which can
take a large part of the budget when the group is expensive, and the
study ends with `is_satisfied` still `FALSE` and a result that the
decision functions refuse. If a group shows no variation after more than
a few observations beyond the pilot, stop and reconsider the measure
before spending the rest of the budget.

## Hypothesis Decisions

When sampling has ended, the estimate and its standard error are used
for the test the study was designed for. The four decision functions are
one family. Each is a one-sided test at level `alpha_level`, 0.05 by
default as in the article’s tables, computed from the contrast estimate
and its standard error, the square root of the row `variance`;
equivalence is two such tests. None of them reads the interval printed
on the check, whose confidence level is two-sided, so the confidence
level of the check does not enter a decision. Below, SE stands for that
standard error and *z*_(α) for the upper 100α percent point of the
standard normal distribution (1.645 for α = 0.05).

The article formulates noninferiority with inferiority as the null
hypothesis (its Equations 6 and 7): for a margin δ \> 0 fixed in advance
(the article’s \|ψ₀\|), the null hypothesis is ψ ≤ −δ and the
alternative is ψ \> −δ. The null hypothesis is rejected at one-sided
level α when the estimate is above −δ + *z*_(α) SE (Equation C3). The
formulation assumes that a larger response is better, which is why the
contrast was written with the face to face format entering negatively.
The `margin` argument is the size of the margin, a positive number. The
article illustrates equivalence with a margin of 1.25 points; the same
margin is used here.

``` r

decide_seq_budget_noninferiority(out, margin = 1.25)
#> [1] TRUE
```

The estimate is 0.85 with a standard error of 0.49, so the estimate
minus 1.645 standard errors, the lower limit of a 90 percent interval,
is 0.04, which is above −1.25: inferiority is rejected, and
noninferiority is declared.

Equivalence is two one-sided tests, each at level α (Equations 11 and
C15): the contrast is declared to lie between the margins ψ₁ and ψ₂ when
the estimate is above ψ₁ + *z*_(α) SE and below ψ₂ − *z*_(α) SE. The two
margins are given as they are, with their signs.

``` r

decide_seq_budget_equivalence(out, margin_lower = -1.25, margin_upper = 1.25)
#> [1] FALSE
```

The 90 percent interval runs from 0.04 to 1.65. Its upper limit is above
1.25, so equivalence within 1.25 points is not declared. The interval is
not too wide for the margins: it is 1.61 wide, and the band between them
is 2.5 wide. It crosses a margin because the estimate in this run lies
1.7 standard errors above the true contrast of zero. The scores were
simulated from formats that are equally effective, so the failure to
declare equivalence is a Type II error of this run and not a verdict on
the budget.

How often a study on this budget declares equivalence when the formats
are equally effective can be worked out before any data are collected.
The test declares it when the estimate is within 1.25 − *z*_(α) SE of
zero. Take for SE the smallest standard error this budget allows, from
the plan in the section on a priori planning, and take the estimate to
be normally distributed around a true contrast of zero:

``` r

se_best <- sqrt(plan$minimum_variance)    # smallest standard error
reach <- 1.25 - qnorm(0.95) * se_best     # declared when |estimate| < reach
c(reach = reach,
  noninferiority = pnorm(reach / se_best),
  equivalence = 2 * pnorm(reach / se_best) - 1)
#>          reach noninferiority    equivalence 
#>      0.4534964      0.8254952      0.6509904
```

By this approximation the budget declares noninferiority in about 83
percent of studies and equivalence in about 65 percent when the true
contrast is zero. Equivalence asks more of a study than noninferiority
does, because the estimate has to stay clear of both margins.

The other two functions are the one-sided tests against a single margin,
which is given with its sign. For a margin *m*,
[`decide_seq_budget_superiority()`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md)
tests the null hypothesis ψ ≤ *m* against the alternative ψ \> *m*. With
a positive margin it is the article’s superiority test (Equations 8 and
9), with a margin of zero it asks whether the alternatives are better
than face to face training at all, and noninferiority is its special
case `margin = -abs(delta)`:

``` r

decide_seq_budget_superiority(out, margin = -abs(1.25))
#> [1] TRUE
decide_seq_budget_superiority(out, margin = 0)
#> [1] TRUE
```

The first line repeats the noninferiority decision, because a contrast
above −δ and noninferiority with margin δ are the same hypothesis. The
second declares the alternatives better than face to face training,
because the lower limit of the 90 percent interval (0.04) is just above
zero. These scores were simulated with a true contrast of zero, so this
is a Type I error, which a test at level 0.05 makes in about one study
in twenty when the true contrast equals the margin. It has the same
cause as the failure to declare equivalence: one estimate that came out
well above the true value.
[`decide_seq_budget_inferiority()`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md)
is the mirror image, the null hypothesis ψ ≥ *m* against the alternative
ψ \< *m*: it concludes that the contrast is below the margin, with that
conclusion as the alternative hypothesis. It is not the article’s
noninferiority test, in which inferiority is the null hypothesis.

``` r

decide_seq_budget_inferiority(out, margin = 0)
#> [1] FALSE
```

These two tests are the halves of the equivalence test. In the article
they are its Equations C11 to C14: equivalence between ψ₁ and ψ₂ is
declared exactly when the superiority test at ψ₁ and the inferiority
test at ψ₂ both reject.

``` r

decide_seq_budget_superiority(out, margin = -1.25) &&
  decide_seq_budget_inferiority(out, margin = 1.25)
#> [1] FALSE
```

A stricter test is a smaller `alpha_level`. The limits that a test at
level α compares with its margins are those of a check computed with
`conf_level = 1 - 2 * alpha_level`; for the default level that is the 90
percent interval quoted above.

``` r

decide_seq_budget_superiority(out, margin = 0, alpha_level = 0.025)
#> [1] FALSE

out_90 <- seq_budget(scores, weights, cost_per_group = costs,
                     budget = budget, conf_level = 1 - 2 * 0.05)
c(term(out_90, "lower_limit"), term(out_90, "upper_limit"))
#> [1] 0.04200302 1.64818417
```

At one-sided level 0.025 the alternatives are no longer declared better
than face to face training.

The level of each test is stated for one decision, made when sampling
has ended. Asked for a decision earlier, the functions return it with a
warning that it is provisional.

## What to Expect Before the Study

[`seq_budget_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)
runs the whole procedure many times for assumed means, standard
deviations, costs, and budget, and reports what it delivers. Here the
alternatives are taken to be exactly as effective as face to face
training, so the rejection rate of the noninferiority test is its power
at a true contrast of zero.

``` r

sens <- seq_budget_sensitivity(
  true_means       = true_means,
  true_sds         = true_sds,
  contrast_weights = weights,
  cost_per_group   = costs,
  budget           = budget,
  margin           = 1.25,
  G                = 400,
  seed             = 113,
  m_0              = 10
)
print(sens, digits = 4)
#> term                           value
#> true_contrast                  0
#> budget                         50000
#> conf_level                     0.95
#> replications                   400
#> coverage                       0.93
#> mean_total_cost                49910
#> sd_total_cost                  92.09
#> pct_within_budget              100
#> pct_rule_met                   0
#> mean_n_1                       61.59
#> mean_n_2                       66.38
#> mean_n_3                       36.63
#> optimal_n_1                    61.95
#> optimal_n_2                    65.98
#> optimal_n_3                    36.51
#> mean_estimated_variance        0.2319
#> mean_true_variance             0.2363
#> empirical_variance             0.2336
#> optimal_variance               0.2345
#> margin                         1.25
#> alpha_level                    0.05
#> noninferiority_rejection_rate  0.845
```

No replication went over the budget (`pct_within_budget` is 100), which
holds by construction. The mean final sizes are close to the optimal
ones. Three rows describe the variance of the estimator and are to be
read against `optimal_variance`, the smallest the budget allows (0.234):
the mean of the estimated variance at stopping (0.232), the mean true
variance at the final sample sizes (0.236, which no allocation within
the budget can push below the optimum), and the variance actually
observed over the replications (0.234, the noisiest of the three). The
interval covered the true contrast in 93.0 percent of the 400
replications (Monte Carlo standard error 1.3 percentage points), and
with a margin of 1.25 this budget rejects inferiority in about 84.5
percent of studies when the formats are in fact equally effective (Monte
Carlo standard error 1.8 percentage points), which agrees with the 83
percent that the normal approximation gave in the section on hypothesis
decisions.

Two cautions. First, the interval and the tests are asymptotic. The
estimated variance tends to run slightly below the true one (sample
standard deviations underestimate population ones on average, and the
allocation follows them). The Details of
[`?seq_budget_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)
report runs of 10,000 replications of five designs with a pilot of 10
and about 15 to 120 observations per group: coverage of the 95 percent
interval was 0.936 to 0.949 for normal populations and 0.925 to 0.941
for skewed ones, and the noninferiority test, with the true contrast
exactly at the margin and a nominal level of 0.05, rejected in 0.050 to
0.061 of studies for normal populations and in 0.040 to 0.085 for skewed
ones. A smaller pilot makes the shortfall larger, and the minimum pilot
of two makes it much larger; the help page of
[`seq_budget_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md)
gives measured figures, and the `m_0` argument of
[`seq_budget_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)
shows the effect for a design of interest. Second, the procedure makes
the most of a budget; it does not promise a given power. A study that
needs more must raise the budget or accept a wider margin. (The values
above are shown to four significant digits; the stored values are not
rounded. A reported analysis deserves more replications than the 400
used to keep this vignette quick.)

## The Article’s First Passage Rule

In the article a population stops being sampled the first time it meets
its criterion and is never sampled again, and the populations still open
receive one observation each per stage, in decreasing order of their
standard deviations (pp. 429 to 430). A check sees only the data in
hand, so by default
[`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
judges every group afresh at each check: a group that met its criterion
earlier and has fallen below it again, because the other standard
deviations moved, receives further observations, and the next
observation always goes to the group furthest below its estimated
optimal size. That default is this package’s choice, not the article’s
procedure. For the article’s first passage rule, hand each check the
groups that the previous check closed and ask for the article’s order:

``` r

set.seed(113)
scores_fp <- lapply(1:3, new_score, k = m_0)
closed <- NULL
repeat {
  out_fp <- seq_budget(scores_fp, weights, cost_per_group = costs,
                       budget = budget, closed = closed)
  closed <- attr(out_fp, "met_per_group")
  g <- seq_budget_next(out_fp, order = "rounds")
  if (is.na(g)) break
  scores_fp[[g]] <- c(scores_fp[[g]], new_score(g))
}
attr(out_fp, "n_per_group")
#> [1] 61 71 33
attr(out_fp, "rule_met")
#> [1] TRUE
attr(out_fp, "remaining_budget")
#> [1] 600
```

Under this rule every group can have met its criterion, each at its own
time, with part of the budget still unspent. (The loop is the package’s
reading of the article’s sequence rather than a step for step copy:
[`?seq_budget_next`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
lists three details that differ, the main one being that a check judges
every open group each time, where the article judges a population only
after an observation drawn from it.) The two readings can be compared
for a design with the `procedure` argument of
[`seq_budget_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md):

``` r

sens_fp <- seq_budget_sensitivity(
  true_means, true_sds, weights, costs, budget,
  margin = 1.25, G = 400, seed = 113, m_0 = 10,
  procedure = "first_passage"
)
rows <- c("mean_total_cost", "pct_rule_met", "mean_n_1", "mean_n_2",
          "mean_n_3", "mean_true_variance", "coverage")
data.frame(
  term = rows,
  reallocate = round(sens$value[match(rows, sens$term)], 3),
  first_passage = round(sens_fp$value[match(rows, sens_fp$term)], 3)
)
#>                 term reallocate first_passage
#> 1    mean_total_cost  49909.875     49849.625
#> 2       pct_rule_met      0.000        24.750
#> 3           mean_n_1     61.587        61.748
#> 4           mean_n_2     66.382        65.598
#> 5           mean_n_3     36.635        36.545
#> 6 mean_true_variance      0.236         0.237
#> 7           coverage      0.930         0.938
```

Under the article’s rule about 25 percent of the studies ended with
every group having met its criterion and money left over; under the
default none did, and on average 60 more of the budget was spent. The
two readings give intervals of practically the same quality. In the
simulations reported in the Details of
[`?seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
(three designs from Table 1 of the article, each with a budget of 5,000
and a pilot of 10, for normal and for skewed populations), the default
left about 10 to 12 of the budget unspent where the first passage rule
left 19 to 42, the true variance of the estimator at the final sample
sizes was 0.2 to 1.0 percent smaller under the default, and coverage
differed by at most 0.006, in either direction. The difference in
coverage in the table above (0.930 against 0.938) is within Monte Carlo
error at 400 replications, where each figure has a standard error of
about 0.013.

## References

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods, 30*(2),
425–439. [doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

Luh, W. M., & Guo, J. H. (2016). Sample size planning for the
noninferiority or equivalence of a linear contrast with cost
considerations. *Psychological Methods, 21*(1), 13–34.
[doi:10.1037/met0000039](https://doi.org/10.1037/met0000039)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. <ISBN:9780367386535>
