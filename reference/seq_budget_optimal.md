# Optimal Sample Sizes for a Budget Constrained Contrast

For planning before any data are collected: given assumed population
standard deviations, returns the sample sizes that estimate a linear
contrast of *K* population means as precisely as a fixed budget allows
when an observation costs a different amount in each population. This is
Theorem 1 of Chattopadhyay, Bandyopadhyay, Kelley, and Padalunkal
(2025).

## Usage

``` r
seq_budget_optimal(contrast_weights, sigmas, cost_per_group, budget)
```

## Arguments

- contrast_weights:

  Numeric vector of length *K* (at least 2) giving the contrast
  coefficients \\c_i\\, one per population, not all zero. The article
  takes them to sum to zero (a contrast); the allocation holds for any
  linear combination, so the sum is not checked. A population with a
  weight of zero is allocated no observations.

- sigmas:

  Numeric vector of length *K* giving the assumed population standard
  deviations \\\sigma_i\\, each positive, in the order of the weights.

- cost_per_group:

  Numeric vector of length *K* giving the cost \\a_i\\ of one
  observation from each population, each positive, in the units of
  `budget`. The order is that of the weights.

- budget:

  The total budget \\A_0\\: a single positive number.

## Value

A list with two components: `n_per_group`, a numeric vector of length
*K* with the optimal sample sizes \\n\_{i0}\\ (not rounded; their total
cost equals `budget`), and `minimum_variance`, the minimized variance
\\V_o\\ of the contrast estimator.

## Details

Write the contrast as \\\psi = \sum\_{i=1}^{K} c_i \mu_i\\, estimated by
\\\hat{\psi} = \sum\_{i=1}^{K} c_i \bar{X}\_i\\ with variance
\\\sum\_{i=1}^{K} c_i^2 \sigma_i^2 / n_i\\ (Equations 3 to 5 of the
article). With a cost of \\a_i\\ per observation from population \\i\\
and a budget \\A_0\\, the sample sizes that minimize that variance
subject to \\\sum\_{i=1}^{K} a_i n_i = A_0\\ are (Equation 12)
\$\$n\_{i0} = \frac{A_0 \\ \|c_i\| \\ \sigma_i} {\sqrt{a_i} \\
\sum\_{l=1}^{K} \|c_l\| \\ \sigma_l \\ \sqrt{a_l}}, \quad i = 1, \ldots,
K,\$\$ and the minimized variance is (Equation 13) \$\$V_o =
\frac{1}{A_0} \left( \sum\_{i=1}^{K} \|c_i\| \\ \sigma_i \\ \sqrt{a_i}
\right)^2.\$\$ A population receives more observations when its scores
vary more, when its weight in the contrast is larger, and when it is
cheaper to sample. Minimizing the variance is what maximizes the power
of the article's noninferiority, superiority, and equivalence tests for
a given budget.

The sizes are returned as they come out of the formula, not rounded.
Rounding every size up would cost more than the budget: in the article's
example (three populations, 30 per observation, a budget of 10,000; p.
430) each \\n\_{i0}\\ is 111.11, which holds when the products \\\|c_i\|
\sigma_i\\ are equal, and 111 per group costs 9,990 where 112 per group
would cost 10,080.

In practice the standard deviations are not known.
[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
replaces them with the sample standard deviations and updates the
allocation as observations arrive, and
[`seq_budget_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)
shows by simulation how close that procedure comes to the allocation
returned here.

The three vector arguments are matched by position. If more than one of
them is named, the names must be the same and in the same order;
otherwise the function stops, so that a vector named in another order is
not applied to the wrong populations.

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

## See also

[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
for the sequential procedure when the standard deviations are unknown.
In DMAR,
[`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.html)
reports the same cost optimal allocation for a sequential fixed-width
interval.

Other seq_budget:
[`decide_seq_budget`](https://yelleknek.github.io/SMSD/reference/decide_seq_budget.md),
[`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md),
[`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md),
[`seq_budget_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md),
[`seq_budget_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Two new training formats (groups 1 and 3) against the current one
# (group 2): the weights 1, -2, 1 give twice the difference between the
# average of the new formats and the current format. An observation costs
# 38, 15, and 10 in the three groups, and the budget is 5,000 (the normal
# row of Table 1 in the article).
seq_budget_optimal(
  contrast_weights = c(1, -2, 1),
  sigmas = c(4.5, 0.9, 1),
  cost_per_group = c(38, 15, 10),
  budget = 5000
)
#> $n_per_group
#> [1] 96.37295 61.35660 41.74788
#> 
#> $minimum_variance
#> [1] 0.2868806
#> 

# The article's example on p. 430: equal costs, and equal products of
# weight and standard deviation.
seq_budget_optimal(
  contrast_weights = c(-1, 0.5, 0.5),
  sigmas = c(1, 2, 2),
  cost_per_group = c(30, 30, 30),
  budget = 10000
)$n_per_group
#> [1] 111.1111 111.1111 111.1111
```
