# Sequential AIPE Check for a Correlation Coefficient

Checks, with the paired data collected so far, whether enough pairs are
in hand for a confidence interval for a population correlation that is
no wider than a chosen width, without assuming bivariate normality and
without a planning value for the correlation. The Pearson product moment
correlation, Kendall's tau, and Spearman's rho are available. This is
the sequential accuracy in parameter estimation (AIPE) procedure of
Kelley, Darku, and Chattopadhyay (2019).

## Usage

``` r
seq_aipe_r(
  x,
  y,
  width,
  conf_level = 0.95,
  method = c("pearson", "kendall", "spearman")
)
```

## Arguments

- x, y:

  Numeric vectors: the two variables observed so far on the same units,
  paired by position. They must have the same length, with at least 4
  pairs; missing or infinite values are an error.

- width:

  The target full width of the confidence interval for the population
  correlation (\\\omega\\ in the cited article): the upper limit minus
  the lower limit. A single positive number less than 2: a correlation
  lies between -1 and 1, so an interval 2 or more wide would say
  nothing, and such a target is an error.

- conf_level:

  The confidence level of the interval, a single number strictly between
  0 and 1. Default 0.95.

- method:

  Which correlation coefficient to estimate: `"pearson"` (the default),
  `"kendall"`, or `"spearman"`.

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

- `r`:

  the estimated correlation (Pearson's *r*, Kendall's tau, or Spearman's
  rho, according to `method`, as
  [`cor`](https://rdrr.io/r/stats/cor.html) computes it, so that with
  ties Kendall's tau is tau-b; `NA` when either variable has no
  variation);

- `upper_limit`:

  the upper confidence limit (`NA` whenever `lower_limit` is);

- `xi_squared`:

  the estimated asymptotic variance \\\hat{\xi}\_n^2\\, at least
  \\n^{-3}\\ (the bound itself when the estimate is not yet usable);

- `criterion_n`:

  the number of pairs the stopping rule currently requires (at least
  `current_n + 1` when the variance estimate is not yet usable or the
  estimate cannot be computed);

- `current_n`:

  the current number of pairs.

The object also carries these attributes:

- `is_satisfied`:

  `TRUE` when the stopping rule is met, otherwise `FALSE`; never `NA`,
  and always identical to `current_n >= criterion_n`;

- `method`:

  a string such as `"Pearson correlation"`, with the name of the
  coefficient that `method` selects;

- `cor_method`:

  the value of the argument `method`;

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
[`seq_aipe_r_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r_pilot.md)),
calls `seq_aipe_r()`, and, while the rule is not met, adds one or more
pairs and calls it again. The function prints nothing and keeps no
state.

**Stopping rule.** Sampling stops at the first \\n\\, at or after the
pilot size, with \$\$n \ge \left\lceil \frac{4
z\_{\alpha/2}^2}{\omega^2} \left(\hat{\xi}\_n^2 + \frac{1}{n}\right)
\right\rceil\$\$ (their Equation 10 for the Pearson correlation,
Equation 30 for Kendall's tau, and Equation 54 for Spearman's rho),
where \\\hat{\xi}\_n^2\\ estimates the asymptotic variance of
\\\sqrt{n}\\ times the sample correlation. The reported interval is
\\r_n \pm z\_{\alpha/2} \hat{\xi}\_n / \sqrt{n}\\ (their Equation 11).
It is symmetric on the correlation scale: the article does not transform
to Fisher's *z*, and neither does the package, so a limit can fall
outside \\\[-1, 1\]\\ at small sample sizes. Whenever the rule is met,
the width of the interval is at most `width`.

**Variance estimators.** For the Pearson correlation, \\\hat{\xi}\_n^2\\
is the moment based estimator of the article's Equation 70, a sample
version of the distribution free variance in its Equation 3 (which
equals \\(1 - \rho^2)^2\\ under bivariate normality). Equation 70 is
printed with divisions by the sample covariance; the factor \\r_n^2\\ in
front cancels them, and the estimate is computed with that cancellation
carried out, so that it is defined when the sample covariance is exactly
zero (which happens with discrete data such as ratings). For Kendall's
tau it is the estimator of the article's Equations 27 to 29. For
Spearman's rho it is the estimator of Kojadinovic and Yan (2010) given
as the article's Equations 51 to 53. (The article states its Spearman
stopping rule with the estimator of Genest and Favre, 2007, its Equation
45, and offers Equation 51 as a simpler consistent alternative; Equation
51 is what the package computes, and it reproduces the article's Table
10.)

**Tied values.** The Kendall and Spearman variance estimators assume
continuous data. With tied values (ratings on a 5 point scale, for
example) they are too small, so the rule stops too early and the
interval is too narrow: with bivariate normal data (\\\rho = 0.5\\) cut
into five categories, the estimates averaged 78% (Kendall) and 85%
(Spearman) of the actual sampling variance at \\n = 200\\, and 95%
intervals at that fixed sample size covered 91% and 92% of the time
(6,000 simulated samples). The function warns, once per call, when
`method` is `"kendall"` or `"spearman"` and either variable has tied
values. With ties the Kendall estimate reported is also a different
coefficient from the one the variance estimator belongs to: it is tau-b,
as [`cor`](https://rdrr.io/r/stats/cor.html) computes it, whereas the
article's estimator (its Equation 23) is the average over pairs of
observations of the product of the signs of their differences (tau-a).
The two are the same number when there are no ties and can differ
materially when there are: 0.334 against 0.265 for 200 bivariate normal
pairs with \\\rho = 0.5\\, both variables cut into five categories.

**A perfect sample correlation.** A sample correlation of 1 or \\-1\\
(to within rounding error) says nothing about its own sampling variance.
For the Pearson correlation the variance estimate is then zero, and so
is the Kendall estimate for pairs ordered in opposite ways on the two
variables; the Kendall estimate for pairs ordered the same way, and the
Spearman estimate for pairs ordered either way, are positive, and with
them the rule as printed would be met at once with a limit outside
\\\[-1, 1\]\\ (\\-1.04\\ and \\-0.96\\ for Spearman's rho with 12 pairs
in exactly reversed order at `width = 0.4`). The function treats the
three alike: the variance estimate of a perfect sample correlation is
not yet usable (see the section below), so the rule is reported as not
met until a pair breaks the perfect order.

**What the guarantee is.** The width guarantee is exact. The coverage
guarantee is asymptotic: coverage approaches `conf_level` as `width`
shrinks (their Theorem 1). At wide targets, and so small final sample
sizes, coverage falls short, and the article says so: its Table 2
reports coverage of 0.90 for a 95% interval of width 0.2 when \\\rho =
0.5\\, and its Appendix D attributes the shortfall to the bias of the
ratios of moment estimates in the variance estimator at small \\n\\,
which makes the procedure stop early. In simulated studies for this
package (bivariate normal data, a 95% Pearson interval, the default
`m_0`, one pair added at each stage), coverage was 0.93 for \\\rho =
0.4\\ and 0.69 for \\\rho = 0.8\\ at `width = 0.2` (5,000 and 10,000
studies), and 0.78 and 0.69 at `width = 0.4` (10,000 studies each). A
large correlation is the hard case: its interval is short, so the rule
can be met with few pairs (about 43 on average for \\\rho = 0.8\\ at
`width = 0.2`), and in samples that small the variance estimate is often
unusable (it was zero or negative in 42% of the pilot samples in that
condition, and in 14% for \\\rho = 0.4\\ at `width = 0.4`) and otherwise
tends to run low. When it is positive but below \\n^{-3}\\ the bound is
used, as the article prescribes, and the rule is then met with an
interval far narrower than requested: in other runs of 6,000 studies at
`width = 0.4`, 58 with \\\rho = 0.8\\ and 10 with \\\rho = 0.4\\ ended
that way, and in runs of 2,000 at `width = 0.2`, 10 and none. A larger
minimum pilot helps: with `m_0 = 50` the coverages were 0.93 and 0.83 at
`width = 0.2` and 0.92 and 0.83 at `width = 0.4`. For Kendall's tau and
Spearman's rho with \\\rho = 0.4\\, coverage was 0.94 and 0.94 at
`width = 0.2` and 0.91 and 0.87 at `width = 0.4`; their variance
estimates were never zero or below the lower bound. A stronger
association is the hard case for them too: with \\\rho = 0.8\\ or
\\-0.8\\, coverage was 0.91 to 0.92 for Kendall's tau and 0.87 to 0.88
for Spearman's rho at `width = 0.2`, and 0.89 to 0.91 and 0.84 to 0.85
at `width = 0.4` (0.91 and 0.85 for \\\rho = 0.8\\ in 10,000 studies).
For a negative association the Kendall interval covers a little less
often at a wide target than for a positive one of the same size: at
`width = 0.4` coverage was 0.90 for \\\rho = -0.4\\, and 0.89 for \\\rho
= -0.95\\ against 0.93 for \\\rho = 0.95\\ (6,000 studies each); at
`width = 0.2` it was 0.94 for \\\rho = -0.4\\ and for \\\rho = -0.95\\
(2,000 studies each). A strongly negative association can also make a
small sample perfectly discordant (any two pairs ordered in opposite
ways on the two variables). Kendall's tau is then \\-1\\ and its
variance estimate is zero, so the check is not yet usable until a
concordant pair has been observed; this happened in 91 of the 6,000
pilot samples of 10 pairs with \\\rho = -0.95\\, in 5 with \\\rho =
-0.8\\, and in none of the pilot samples of 20 pairs. Use
[`seq_aipe_r_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
to study a planned Pearson correlation study before collecting data, and
prefer a narrow target or a larger `m_0` when the correlation may be
large.

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

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
accuracy in parameter estimation for population correlation
coefficients. *Psychological Methods*, *24*(4), 492–515.
[doi:10.1037/met0000203](https://doi.org/10.1037/met0000203)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Genest, C., & Favre, A.-C. (2007). Everything you always wanted to know
about copula modeling but were afraid to ask. *Journal of Hydrologic
Engineering*, *12*, 347–368.
[doi:10.1061/(ASCE)1084-0699(2007)12:4(347)](https://doi.org/10.1061/%28ASCE%291084-0699%282007%2912%3A4%28347%29)

Kojadinovic, I., & Yan, J. (2010). Comparison of three semiparametric
methods for estimating dependence parameters in copula models.
*Insurance: Mathematics and Economics*, *47*, 52–63.
[doi:10.1016/j.insmatheco.2010.03.008](https://doi.org/10.1016/j.insmatheco.2010.03.008)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology*, *59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

## See also

[`seq_aipe_r_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r_pilot.md)
for the pilot size and
[`seq_aipe_r_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
to simulate the Pearson procedure before collecting data. In DMAR:
[`ss_aipe_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.html),
the fixed sample size analog for when a planning value of the
correlation is available, and
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.html),
a fixed sample size confidence interval.

Other seq_aipe_r:
[`seq_aipe_r_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r_pilot.md),
[`seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# How strongly are citizenship behavior at work and neuroticism related?
# The aim is a 95% confidence interval for the correlation that is no
# wider than 0.30.
set.seed(113)
new_pairs <- function(k, rho = -0.3) {
  citizenship <- rnorm(k)
  cbind(citizenship, rho * citizenship + sqrt(1 - rho^2) * rnorm(k))
}
scores <- new_pairs(seq_aipe_r_pilot(width = 0.30))
out <- seq_aipe_r(scores[, 1], scores[, 2], width = 0.30)
out
#> Sequential check (Pearson correlation)
#> 
#>   Current n       : 14
#>   Criterion n     : 102
#>   Target width    : 0.3
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : -0.7589524
#>   r                 : -0.3790672
#>   upper_limit       : 0.0008179968
#>   xi_squared        : 0.5259405
#>   criterion_n       : 102
#>   current_n         : 14
#> 
#>   Stopping rule met: FALSE

# The rule is not met at the pilot. Add one employee at a time and check
# again, until it is.
while (!attr(out, "is_satisfied")) {
  scores <- rbind(scores, new_pairs(1))
  out <- seq_aipe_r(scores[, 1], scores[, 2], width = 0.30)
}
out
#> Sequential check (Pearson correlation)
#> 
#>   Current n       : 125
#>   Criterion n     : 124
#>   Target width    : 0.3
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : -0.3967731
#>   r                 : -0.2487647
#>   upper_limit       : -0.1007562
#>   xi_squared        : 0.7128316
#>   criterion_n       : 124
#>   current_n         : 125
#> 
#>   Stopping rule met: TRUE

# Kendall's tau from the same data and target
seq_aipe_r(scores[, 1], scores[, 2], width = 0.30, method = "kendall")
#> Sequential check (Kendall correlation)
#> 
#>   Current n       : 125
#>   Criterion n     : 70
#>   Target width    : 0.3
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : -0.2759581
#>   r                 : -0.1654194
#>   upper_limit       : -0.05488058
#>   xi_squared        : 0.397597
#>   criterion_n       : 70
#>   current_n         : 125
#> 
#>   Stopping rule met: TRUE
```
