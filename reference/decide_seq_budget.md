# Noninferiority and Equivalence Decisions for a Budget Constrained Contrast

Four tests of the contrast estimated by
[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md),
for use once sampling has ended. Each is a one-sided test at level
`alpha_level` (equivalence is two one-sided tests, each at that level),
computed from the contrast estimate and its standard error. The
noninferiority and equivalence functions are the tests of Chattopadhyay,
Bandyopadhyay, Kelley, and Padalunkal (2025). The superiority and
inferiority functions are the one-sided tests against a single margin:
the article's superiority test and its noninferiority test are special
cases of the first, and the two together are the two halves of its
equivalence test.

## Usage

``` r
decide_seq_budget_noninferiority(check, margin, alpha_level = 0.05)

decide_seq_budget_superiority(check, margin, alpha_level = 0.05)

decide_seq_budget_inferiority(check, margin, alpha_level = 0.05)

decide_seq_budget_equivalence(
  check,
  margin_lower,
  margin_upper,
  alpha_level = 0.05
)
```

## Arguments

- check:

  An object returned by
  [`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md).
  A check from any other procedure is refused.

- margin:

  A single finite number on the scale of the contrast. For
  noninferiority it is the size of the noninferiority margin and must be
  positive; inferiority, the null hypothesis, is that the contrast is at
  or below `-margin`. For superiority and inferiority it is the value
  the contrast is compared with, given with its sign (negative, zero, or
  positive).

- alpha_level:

  The one-sided level of the test: a single number strictly between 0
  and 0.5 (default 0.05). For equivalence it is the level of each of the
  two one-sided tests.

- margin_lower, margin_upper:

  For equivalence: the lower and upper equivalence margins, single
  finite numbers on the scale of the contrast given with their signs,
  with `margin_lower < margin_upper` (for example `-3` and `3`).

## Value

A length one logical: `TRUE` if the null hypothesis is rejected, so that
the alternative (noninferiority, superiority, inferiority, or
equivalence) is supported, and `FALSE` otherwise. `FALSE` is not
evidence for the null hypothesis.

## Details

Write \\\psi\\ for the population contrast, \\\hat{\psi}\\ for its
estimate (the `contrast` row of the check), \\\hat{V}\\ for the
estimated variance of the estimate (the `variance` row), and
\\z\_{\alpha}\\ for the upper \\100\alpha\\ percent point of the
standard normal distribution, with \\\alpha\\ the one-sided level
`alpha_level`.

**The level.** Every decision is computed from those two rows of the
check and `alpha_level`. None is read from the confidence limits printed
on the check, which are two-sided at the check's `conf_level`, so a
decision does not depend on that confidence level. The limits a decision
compares with its margin, \\\hat{\psi} \pm z\_{\alpha} \sqrt{\hat{V}}\\,
are those of a \\1 - 2\alpha\\ interval: a 90 percent interval for the
default `alpha_level = 0.05`, the level of the article's tables. This is
how the sibling functions
[`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.html)
and
[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.html)
in DMAR state the level.

**The direction.** The article assumes that a larger response is better
(p. 428), so code the contrast so that larger values of \\\psi\\ favor
the treatment of interest. If a smaller response is better, reverse the
sign of every contrast weight first.

**How close the level is to nominal.** All four tests are asymptotic
(the article relies on the central limit theorem, Appendix A), and their
level is not exact in samples of realistic size. In simulations of the
noninferiority test at nominal level 0.05 (10,000 studies per condition,
standard error about 0.002; pilot of 10 per group; the true contrast
exactly at the margin), the rejection rate was 0.050 to 0.057 for normal
populations when the budget bought about 30 to 120 observations per
group. In a smaller design with 15 to 28 per group it was about 0.058
(0.052 to 0.061 over seven such runs). A one-sided test is more
sensitive to skewness than a two-sided interval: for skewed populations
the rate ranged from 0.040 to 0.085, the largest where the groups with
negative weights were small (about 29 and 48 observations). The
inferiority test rejects in the other tail. In two of these designs (the
normal row of the article's Table 1 with a budget of 5,000, and the
smaller design) its rate at the margin was 0.054 and 0.059 for normal
populations and 0.066 and 0.076 for right skewed ones; the
noninferiority test in the same runs had 0.056 and 0.059, and 0.043 and
0.055. The equivalence test, with the true contrast on one of its
margins, rejects at the rate of the corresponding one-sided test when
the other margin is far away, and less often when the two margins are
close. A smaller pilot raises the rates further;
[`seq_budget_pilot`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md)
has the figures for the noninferiority test. For normal populations the
simulated power agreed with the article's Equation C9 (0.80 to 0.81 in
five designs against 0.804).
[`seq_budget_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)
reports the rejection rate for a design of interest.

**When the decision is made.** The level of each test is stated for one
decision, made when sampling has ended. If the check says that sampling
has not ended (`is_satisfied` is `FALSE`), the decision is still
returned, with a warning that it is provisional. If the variance of a
group is unusable (the check's `variance_unusable` attribute is `TRUE`:
the observations in some group are all equal), the functions stop,
whether or not the budget has ended the sampling, because the decision
would rest on a variance of zero for that group.

## Noninferiority

`decide_seq_budget_noninferiority()` is the article's test for
noninferiority (Equations 6 and 7; in Appendix C, Equations C1 and C3).
With \\\delta\\ the margin, it tests \\H_0: \psi \le -\delta\\ against
\\H_A: \psi \> -\delta\\. In words, the null hypothesis is inferiority:
the treatment of interest falls short of the comparison by the margin or
more. The alternative is noninferiority: any shortfall is smaller than
the margin. `margin` is the size of the margin, a positive number (the
article's \\\|\psi_0\|\\, "the least relevant value of the contrast");
the hypotheses are about `-margin`. The null hypothesis is rejected, and
noninferiority is declared, when (Equation C3) \$\$\hat{\psi} \>
-\delta + z\_{\alpha} \sqrt{\hat{V}}.\$\$ Noninferiority is the special
case `margin = -abs(delta)` of the superiority test below: with `delta`
the noninferiority margin, these two calls return the same decision at
any `alpha_level`.

    decide_seq_budget_noninferiority(check, margin = delta)
    decide_seq_budget_superiority(check, margin = -abs(delta))

## Superiority

`decide_seq_budget_superiority()` is the one-sided test that the
contrast is above a margin \\m\\: \\H_0: \psi \le m\\ against \\H_A:
\psi \> m\\. In words, the null hypothesis is that the contrast does not
exceed the margin, and the alternative is that it does. `margin` is the
value \\m\\ itself, with its sign: it may be positive, zero, or
negative. The null hypothesis is rejected when \$\$\hat{\psi} \> m +
z\_{\alpha} \sqrt{\hat{V}}.\$\$ With a positive margin, \\m =
\|\psi_0\|\\, this is the article's test for superiority (Equations 8
and 9; Equations C2 and C4): the treatment of interest is better than
the comparison by more than the margin. With `margin = 0` it is the
ordinary one-sided test that the contrast is positive. With a negative
margin it is the noninferiority test above. With \\m\\ the lower
equivalence margin \\\psi_1\\, it is the first of the two one-sided
tests that make up the article's test for equivalence (Equations C11 and
C13).

## Inferiority

`decide_seq_budget_inferiority()` is the one-sided test that the
contrast is below a margin \\m\\: \\H_0: \psi \ge m\\ against \\H_A:
\psi \< m\\. In words, the null hypothesis is that the contrast is at or
above the margin, and the alternative, the conclusion the data must
establish, is that it is below the margin. `margin` is the value \\m\\
itself, with its sign. The null hypothesis is rejected when
\$\$\hat{\psi} \< m - z\_{\alpha} \sqrt{\hat{V}}.\$\$ This is the mirror
image of the superiority test: reversing the sign of every contrast
weight and of the margin turns one into the other. In the article these
hypotheses and this rejection rule appear, with \\m\\ the upper
equivalence margin \\\psi_2\\, as the second of the two one-sided tests
that make up the test for equivalence (Equations C12 and C14). As a test
on its own the article states the hypotheses once, in its footnote 1,
with \\m = -\|\psi_0\|\\, as the form its superiority test takes when a
smaller response is better. When a larger response is better,
`margin = -delta` declares the treatment of interest inferior by more
than `delta`. This is not the article's noninferiority test, in which
inferiority is the null hypothesis.

## Equivalence

`decide_seq_budget_equivalence()` is the article's test for equivalence
(Equation 11; in Appendix C, Equations C10 to C15). With \\\psi_1\\ the
lower and \\\psi_2\\ the upper equivalence margin, it tests \\H_0: \psi
\le \psi_1\\ or \\\psi \ge \psi_2\\ against \\H_A: \psi_1 \< \psi \<
\psi_2\\. In words, the null hypothesis is nonequivalence (the contrast
is at or beyond one of the margins) and the alternative is equivalence
(the contrast lies between them). `margin_lower` and `margin_upper` are
the margins \\\psi_1\\ and \\\psi_2\\ themselves, with their signs: for
equivalence within 3 points either way they are `-3` and `3`. The
article takes \\\psi_1 \< 0 \< \psi_2\\; only
`margin_lower < margin_upper` is required here. The null hypothesis is
rejected when (Equation C15) \$\$\psi_1 + z\_{\alpha} \sqrt{\hat{V}} \<
\hat{\psi} \< \psi_2 - z\_{\alpha} \sqrt{\hat{V}}.\$\$ These are two
one-sided tests, each at level \\\alpha\\ (Equations C11 to C14):
equivalence is declared exactly when the superiority test at
`margin_lower` and the inferiority test at `margin_upper`, each at the
same `alpha_level`, both reject.

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

[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md).
In DMAR,
[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.html)
gives the fixed sample size equivalence and noninferiority tests for a
contrast with the same one-sided `alpha_level`, and
[`ss_seq_c`](https://yelleknek.github.io/DMAR/reference/ss_seq_c.html)
plans a sequential fixed-width interval for those decisions.

Other seq_budget:
[`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md),
[`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md),
[`seq_budget_optimal()`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md),
[`seq_budget_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_budget_pilot.md),
[`seq_budget_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Is a cheaper online training format (group 1) no worse, by more than
# 3 points, than the average of two established formats?
set.seed(113)
weights <- c(1, -0.5, -0.5)
costs <- c(10, 38, 15)
true_means <- c(50, 50, 50)
scores <- lapply(true_means, function(m) rnorm(10, mean = m, sd = 6))
repeat {
  out <- seq_budget(scores, weights, cost_per_group = costs, budget = 3000)
  g <- seq_budget_next(out)
  if (is.na(g)) break
  scores[[g]] <- c(scores[[g]], rnorm(1, mean = true_means[g], sd = 6))
}
out
#> Sequential check (budget constrained K-population contrast)
#> 
#>   Current n       : 192
#>   Criterion n     : 194
#>   Confidence level: 95%
#> 
#>   lower_limit       : -0.8960654
#>   contrast          : 0.8640861
#>   upper_limit       : 2.624238
#>   variance          : 0.8064992
#>   current_total_n   : 192
#>   current_total_cost: 2995
#> 
#>   Stopping rule met: FALSE (sampling has ended: the budget is exhausted)
#> 
#>   Note: Sampling has ended because the budget cannot pay for another
#>         observation: 5 of 3000 is left, and one more observation from a
#>         group still below its criterion (groups 1 and 3) costs at least
#>         10.

# Noninferiority with a margin of 3 points, at one-sided level 0.05.
decide_seq_budget_noninferiority(out, margin = 3)
#> [1] TRUE

# The same test, as the special case margin = -abs(3) of superiority.
decide_seq_budget_superiority(out, margin = -abs(3))
#> [1] TRUE

# Is the online format better than the average of the other two at all,
# or worse than it?
decide_seq_budget_superiority(out, margin = 0)
#> [1] FALSE
decide_seq_budget_inferiority(out, margin = 0)
#> [1] FALSE

# Equivalence within 3 points either way: two one-sided tests at 0.05,
# and then at a much stricter level.
decide_seq_budget_equivalence(out, margin_lower = -3, margin_upper = 3)
#> [1] TRUE
decide_seq_budget_equivalence(out, margin_lower = -3, margin_upper = 3,
                              alpha_level = 0.005)
#> [1] FALSE
```
