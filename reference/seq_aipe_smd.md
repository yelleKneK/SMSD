# Sequential AIPE Check for the Standardized Mean Difference

Checks, with the data collected so far from two independent groups,
whether enough observations are in hand for a confidence interval for
the population standardized mean difference \\\delta = (\mu_1 -
\mu_2)/\sigma\\ that is no wider than a chosen width, without assuming
normality and without a planning value for \\\delta\\. This is the
sequential accuracy in parameter estimation (AIPE) procedure of Kelley,
Darku, and Chattopadhyay (2018).

## Usage

``` r
seq_aipe_smd(group_1, group_2, width, conf_level = 0.95)
```

## Arguments

- group_1, group_2:

  Numeric vectors: the observations collected so far from the first and
  the second group. Each needs at least 4 observations, the two must
  have the same length (the procedure adds the same number of
  observations to each group at every stage), and missing or infinite
  values are an error.

- width:

  The target full width of the confidence interval for \\\delta\\
  (\\\omega\\ in the cited article): the upper limit minus the lower
  limit. A single positive number.

- conf_level:

  The confidence level of the interval, a single number strictly between
  0 and 1. Default 0.95.

## Value

An object of class `c("smsd_check", "data.frame")` (see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
for its [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and `tidy()`
methods) with one row per quantity and columns `term` and `value`. The
rows are, in order:

- `lower_limit`:

  the lower confidence limit (`NA` when the variance estimate is not yet
  usable or the estimate cannot be computed; see the section below);

- `smd`:

  the estimated standardized mean difference (`group_1` minus `group_2`,
  in pooled standard deviation units; `NA` when neither group has any
  variation);

- `upper_limit`:

  the upper confidence limit (`NA` whenever `lower_limit` is);

- `xi_squared`:

  the estimated asymptotic variance \\\hat{\xi}\_n^2\\, at least
  \\n^{-3}\\ (the bound itself when the estimate is not yet usable);

- `criterion_n`:

  the per group sample size the stopping rule currently requires (at
  least `current_n + 1` when the variance estimate is not yet usable or
  the estimate cannot be computed);

- `current_n`:

  the current per group sample size.

The object also carries these attributes:

- `is_satisfied`:

  `TRUE` when the stopping rule is met, otherwise `FALSE`; never `NA`,
  and always identical to `current_n >= criterion_n`;

- `method`:

  the string `"standardized mean difference"`;

- `current_n`, `criterion_n`:

  the values in the rows of those names;

- `conf_level`, `width`:

  the values supplied;

- `variance_unusable`:

  `TRUE` when the raw variance estimate is zero, negative, or not finite
  (see the section below), otherwise `FALSE`;

- `notes`:

  a character vector, empty unless something needs explaining.

## Details

The function evaluates the stopping rule once and returns. The user
collects the pilot sample (see
[`seq_aipe_smd_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd_pilot.md)),
calls `seq_aipe_smd()`, and, while the rule is not met, adds one or more
observations to *each* group and calls it again. The function prints
nothing and keeps no state.

**Estimate.** The standardized mean difference is \\d_n =
(\bar{X}\_{1n} - \bar{X}\_{2n}) / s\_{pn}\\, with \\s\_{pn}^2 =
(s\_{1n}^2 + s\_{2n}^2)/2\\ the pooled variance of the two groups of
common size \\n\\ (their Equation 24), computed by
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.html).

**Stopping rule.** Sampling stops at the first \\n\\ (per group), at or
after the pilot size, with \$\$n \ge \left\lceil \frac{4
z\_{\alpha/2}^2}{\omega^2} \left(\hat{\xi}\_n^2 + \frac{1}{n}\right)
\right\rceil\$\$ (their Equation 22), where \\\hat{\xi}\_n^2\\ estimates
the asymptotic variance of \\\sqrt{n}(d_n - \delta)\\ from the sample
means, the pooled variance, and unbiased estimates of the third and
fourth central moments of each group (their Equation 27). The reported
interval is \\d_n \pm z\_{\alpha/2} \hat{\xi}\_n / \sqrt{n}\\ (their
Equation 28). Whenever the rule is met, the width of that interval is at
most `width`.

**Correction to Equation 26.** Equation 26 of the article, and the
sample version printed under its Equation 27, omit a factor of \\1/2\\
in the term that carries the third central moments. The package
implements the corrected variance, \$\$\xi^2 = 2 - \frac{(\mu_1 -
\mu_2)(\mu\_{13} - \mu\_{23})}{2 \sigma^4} + \frac{(\mu_1 - \mu_2)^2}{4
\sigma^6} \left(\frac{\mu\_{14} + \mu\_{24}}{4} -
\frac{\sigma^4}{2}\right),\$\$ where \\\mu\_{kj}\\ is the \\j\\th
central moment of group \\k\\. The two forms differ only when the means
differ and the third central moments differ; they agree when the groups
have the same third central moment (in particular for symmetric
distributions) and when the two means are equal. The authors' correction
is posted at
<https://kenkelley.org/wp-content/uploads/articles/Kelley_Darku_Chattopadhyay_Equation_26.pdf>
and is installed with the package, at the path returned by

    system.file("eq26_correction.pdf", package = "SMSD")

**What the guarantee is.** The width guarantee is exact. The coverage
guarantee is asymptotic: the coverage of the interval approaches
`conf_level` as `width` shrinks, and so as the final sample size grows
(their Theorem 4). In simulated studies with a 95% interval, a true
\\\delta\\ of 0.5, and one observation per group added at each stage,
coverage was 0.952 at `width = 0.4` (6,000 studies, about 200 per group
at the end) and 0.946 at `width = 0.2` (2,000 studies, about 790 per
group) for normal data, and 0.944 and 0.949 for right skewed data (a
standardized chi square distribution with 3 degrees of freedom); the
Monte Carlo standard errors are 0.003 and 0.005. The variance estimate
was never zero, negative, or below its lower bound in these studies. At
a wide target the final sample is small, and coverage falls short for
skewed data and, to a lesser degree, for a large \\\delta\\: at
`width = 0.8` (10,000 studies each, about 50 to 105 per group at the
end) coverage was 0.950, 0.946, and 0.943 for normal data and 0.926,
0.926, and 0.919 for the skewed data at \\\delta\\ of 0.5, 1, and 2. Use
[`seq_aipe_smd_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
to study other conditions before collecting data.

**A group with no variation.** When every observation in one group is
the same (values that agree to ten significant digits count as the same)
while the other group varies, the estimate is still reported, but the
variance estimate is not yet usable: the variance of that group has not
been estimated, and the pooled standard deviation rests on the other
group alone, so that it is too small and the estimate too large in size.
The check reports the rule as not met and a note names the group (see
the section below), as
[`seq_aipe_c`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md)
and
[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
do with a constant group. In simulated studies of a 5 point item with
70% of one group's answers at the top of the scale (\\\delta = 1.08\\,
`width = 1.5`, 4,000 studies), 2.4% of the studies would have stopped
with that group constant, and 41% of their intervals covered \\\delta\\;
with the group treated as not yet estimated, coverage over all studies
was 0.927 against 0.917. The price is paid where a group stays constant:
for a binary outcome with proportions 0.30 and 0.02 (`width = 0.8`,
2,000 studies) the mean final sample size was 64 per group against 33,
and coverage 0.967 against 0.943. If a group shows no variation well
past the pilot sample, reconsider the measure.

## Estimates That Are Not Yet Usable

When the raw variance estimate is zero (exactly, or to within rounding
error), negative, or not finite, which is common at small sample sizes
for some effect sizes, or when the point estimate cannot be computed (a
sample with no variation, where values that agree to ten significant
digits count as equal), the check reports the stopping rule as not met,
sets `attr(out, "variance_unusable")` to `TRUE`, and attaches a note
that [`print()`](https://rdrr.io/r/base/print.html) shows; see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md).
This is a deliberate departure from the stopping rule as printed in the
source articles, under which the lower bound \\n^{-3}\\ would stand in
for such an estimate and the rule would be met at once with an interval
far narrower than requested; a positive estimate, however small, is used
as the articles prescribe (the larger of it and \\n^{-3}\\).

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods*, *11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology*, *59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

## See also

[`seq_aipe_smd_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd_pilot.md)
for the pilot size and
[`seq_aipe_smd_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
to simulate the procedure before collecting data. In DMAR:
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.html),
the fixed sample size analog for when a planning value of \\\delta\\ is
available, and
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.html), the point
estimator used here.

Other seq_aipe_smd:
[`seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md),
[`seq_aipe_smd_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Does a training program raise employee engagement (rated 1 to 7)? The
# aim is a 95% confidence interval for the standardized mean difference
# that is no wider than 0.60.
set.seed(113)
m <- seq_aipe_smd_pilot(width = 0.60)
trained <- rnorm(m, mean = 4.9, sd = 1)
control <- rnorm(m, mean = 4.5, sd = 1)
out <- seq_aipe_smd(trained, control, width = 0.60)
out
#> Sequential check (standardized mean difference)
#> 
#>   Current n       : 7
#>   Criterion n     : 91
#>   Target width    : 0.6
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : -1.143786
#>   smd               : -0.1011596
#>   upper_limit       : 0.9414664
#>   xi_squared        : 1.980884
#>   criterion_n       : 91
#>   current_n         : 7
#> 
#>   Stopping rule met: FALSE

# The rule is not met at the pilot. Add one employee to each group and
# check again, until it is.
while (!attr(out, "is_satisfied")) {
  trained <- c(trained, rnorm(1, mean = 4.9, sd = 1))
  control <- c(control, rnorm(1, mean = 4.5, sd = 1))
  out <- seq_aipe_smd(trained, control, width = 0.60)
}
out
#> Sequential check (standardized mean difference)
#> 
#>   Current n       : 85
#>   Criterion n     : 85
#>   Target width    : 0.6
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : 0.05840176
#>   smd               : 0.3564726
#>   upper_limit       : 0.6545434
#>   xi_squared        : 1.965901
#>   criterion_n       : 85
#>   current_n         : 85
#> 
#>   Stopping rule met: TRUE
```
