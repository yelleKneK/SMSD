# Which Group to Sample Next Under a Budget Constraint

Given a check from
[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md),
returns the group that should receive the next observation, or `NA` when
sampling has ended. Only a group whose next observation the remaining
budget can pay for is ever returned, so a study that adds one
observation at a time to the group named here cannot spend more than its
budget.

## Usage

``` r
seq_budget_next(check, order = c("shortfall", "rounds"))
```

## Arguments

- check:

  An object returned by
  [`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md).

- order:

  How the next group is chosen among those that still need an
  observation and can be paid for: `"shortfall"` (default) or
  `"rounds"`; an unambiguous abbreviation is accepted. See Details.

## Value

A length one integer: the position in `groups` (as passed to
[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md))
of the group that should receive the next observation. It is `NA` when
sampling has ended, which is exactly when the check's `rule_met` or
`budget_exhausted` attribute is `TRUE`: either every group has met its
criterion, or the remaining budget cannot pay for an observation from
any group that has not. (That is also when `is_satisfied` is `TRUE`,
unless the variance of a group is unusable.)

## Details

A group is a candidate when it has not met its criterion (it still needs
observations) and one more observation from it costs no more than what
is left of the budget. If there is no candidate the function returns
`NA`, and the check's `rule_met` and `budget_exhausted` attributes say
why. A group that needs observations but cannot be paid for is passed
over in favor of one that can, so the budget is used as fully as the
unit costs allow.

Among the candidates:

- `order = "shortfall"`:

  (the default, and the package's choice) the next observation goes to
  the group furthest below its estimated optimal size, that is, the one
  with the smallest ratio of its current size to its target in Equation
  18 of Chattopadhyay, Bandyopadhyay, Kelley, and Padalunkal (2025).
  This keeps the allocation as close as possible to the estimated
  optimal one at every step.

- `order = "rounds"`:

  the next observation goes to the candidate with the fewest
  observations. With equal pilot sizes this visits the open groups one
  at a time, in rounds, in decreasing order of their standard
  deviations, which is the order the article describes (Stage 2, pp. 429
  to 430; one observation per group per stage). Combine it with the
  `closed` argument of
  [`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
  for the first passage rule.

Ties go to the group with the larger sample standard deviation and then
to the one listed first.

**A group whose variance is unusable**, because every observation in it
is equal, has no target to be measured against (the check marks it in
`not_estimable_per_group`). Under either order it takes its turn by
rounds: it is named when no group whose variance can be used has fewer
observations than it has, and otherwise only when it is the one group
that can still be sampled. It is never sampled ahead of the others, so a
group that stays constant cannot draw the whole budget to itself; it
does keep pace with the least sampled of the other groups, as the help
for
[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
explains. This rule is the package's choice: the article does not treat
the case.

**How the first passage option differs from the article.** With `closed`
and `order = "rounds"` the package follows the article's rule that a
population is not sampled again once it has met its criterion, but three
details differ. (1) The article orders the populations by their first
stage standard deviations; a check does not know those, so the current
ones are used. (2) The article's flowchart (Figure 1) ends the procedure
at the first observation that cannot be paid for, and its text does not
say whether a cheaper open group may still be sampled when a more
expensive one can no longer be afforded; the package continues with the
cheaper group, which uses the budget more fully and still cannot exceed
it. (3) The article judges a population against its criterion only after
an observation has been drawn from that population (Equations 16 and
17); a check judges every open group each time, so a group can be closed
after an observation drawn from another group, sooner than the article
would close it. The differences are small: in simulations of the normal
row of the article's Table 1, the option ended with every group closed
in 31 percent of studies where the article's own sequence, programmed
separately, did in 29 percent; the mean final sample sizes differed by
0.3 or less and the mean cost by about 1 in a budget of 5,000. Even so,
the option is the package's reading of the article's procedure, not a
step for step copy of it.

## References

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods*,
*30*(2), 425–439.
[doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

## See also

[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)

Other seq_budget:
[`decide_seq_budget`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md),
[`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md),
[`seq_budget_optimal()`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md),
[`seq_budget_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md),
[`seq_budget_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Scores after a new training format (group 1) and two existing formats;
# an observation costs 38, 15, and 10.
set.seed(113)
scores <- list(
  rnorm(10, mean = 52, sd = 9),
  rnorm(10, mean = 50, sd = 6),
  rnorm(10, mean = 50, sd = 6)
)
out <- seq_budget(scores, contrast_weights = c(1, -0.5, -0.5),
                  cost_per_group = c(38, 15, 10), budget = 1500)
seq_budget_next(out)
#> [1] 1

# Near the end of a smaller study, 24 is left of a budget of 680. Group 1
# is still below its criterion but costs 38 per observation, so the next
# observation goes to group 3, which also needs one and costs 10.
set.seed(113)
scores <- list(
  rnorm(12, mean = 52, sd = 9),
  rnorm(8, mean = 50, sd = 6),
  rnorm(8, mean = 50, sd = 6)
)
out <- seq_budget(scores, contrast_weights = c(1, -0.5, -0.5),
                  cost_per_group = c(38, 15, 10), budget = 680)
attr(out, "remaining_budget")
#> [1] 24
attr(out, "met_per_group")
#> [1] FALSE  TRUE FALSE
seq_budget_next(out)
#> [1] 3
```
