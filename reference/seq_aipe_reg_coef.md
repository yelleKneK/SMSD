# Sequential AIPE Check for the Simple Linear Regression Slope

Checks, with the paired data collected so far, whether enough pairs are
in hand for a confidence interval for the slope of the simple linear
regression of `y` on `x`, \\\beta_1 = \sigma\_{XY} / \sigma_X^2\\, that
is no wider than a chosen width, without assuming normal or constant
variance errors and without planning values. This is the sequential
accuracy in parameter estimation (AIPE) procedure of Kelley, Darku, and
Chattopadhyay (2018).

## Usage

``` r
seq_aipe_reg_coef(x, y, width, conf_level = 0.95)
```

## Arguments

- x:

  A numeric vector: the predictor values observed so far.

- y:

  A numeric vector: the outcome values observed so far, paired by
  position with `x`. The two must have the same length, with at least 4
  pairs; missing or infinite values are an error.

- width:

  The target full width of the confidence interval for \\\beta_1\\
  (\\\omega\\ in the cited article): the upper limit minus the lower
  limit, in units of `y` per unit of `x`. A single positive number.

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

- `reg_coef`:

  the estimated slope (`NA` when `x` has no variation);

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

  the string `"simple linear regression slope"`;

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
[`seq_aipe_reg_coef_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef_pilot.md)),
calls `seq_aipe_reg_coef()`, and, while the rule is not met, adds one or
more pairs and calls it again. The function prints nothing and keeps no
state.

**Estimate.** The slope is the ordinary least squares estimate \\b\_{1n}
= s\_{XYn} / s\_{Xn}^2\\ (their Equation 35), as returned by
[`slr_slope`](https://yelleknek.github.io/SMSD/reference/slr_slope.md).

**Stopping rule.** Sampling stops at the first \\n\\, at or after the
pilot size, with \$\$n \ge \left\lceil \frac{4
z\_{\alpha/2}^2}{\omega^2} \left(\hat{\xi}\_n^2 + \frac{1}{n}\right)
\right\rceil\$\$ (their Equation 22), where \\\hat{\xi}\_n^2\\ is the
estimator in their Equations 84 and 85 of the asymptotic variance of
\\\sqrt{n}(b\_{1n} - \beta_1)\\ given in their Equation 37, built from
estimates of the joint central moments of `x` and `y`. The reported
interval is \\b\_{1n} \pm z\_{\alpha/2} \hat{\xi}\_n / \sqrt{n}\\ (their
Equation 38). Whenever the rule is met, the width of that interval is at
most `width`.

**Units.** Unlike a standardized effect size, the slope and `width` are
in units of `y` per unit of `x`, while the correction term \\1/n\\ in
the stopping rule and the lower bound \\n^{-3}\\ on the variance
estimate are pure numbers, as printed in the article. They play their
intended roles (a small correction, and a bound far below any plausible
estimate) when \\\xi^2\\ is not many orders of magnitude from one, that
is, when the residual standard deviation of `y` is of the same order as
the standard deviation of `x`. If the two variables are recorded in very
different units, rescale one of them (and `width` with it) first. When
the variance estimate is positive but below \\n^{-3}\\, the bound is
used in its place, as the article prescribes, so the interval is wider
and the required sample size larger than the estimate itself would give;
a note on the returned object says that this happened and advises
rescaling. On ordinary scales that is a rare chance event of a small
sample: with \\\xi^2 = 1\\ it happened in 0.05% of samples of 10 pairs
and 0.002% of samples of 20 pairs (200,000 simulated samples of normal
data each), against 16% and 1% in which the estimate was zero or
negative.

**What the guarantee is.** The width guarantee is exact. The coverage
guarantee is asymptotic: coverage approaches `conf_level` as `width`
shrinks (their Theorem 7). At wide targets, and so small final sample
sizes, it falls short. In simulated studies (a 95% interval, the default
`m_0`, one pair added at each stage, a slope of 0.5, and `x` and the
errors with unit variance, so that \\\xi^2 = 1\\), coverage for normal
`x` and errors was 0.94 at `width = 0.2` (5,000 studies) and 0.83 at
`width = 0.4` (10,000 studies); for right skewed `x` and errors (a
standardized chi square distribution with 3 degrees of freedom) it was
0.86 (2,000 studies) and 0.71 (10,000 studies). The reason is the
variance estimate at small sample sizes. It was zero or negative in
16.1% (normal) and 26.5% (skewed) of pilot samples of 10 pairs, the
pilot at `width = 0.4`, and in 1.3% and 9.9% of pilot samples of 20, the
pilot at `width = 0.2` (200,000 simulated pilot samples each); and when
it is usable it tends to run low. When it is positive but below
\\n^{-3}\\ the bound is used and the rule is then met with an interval
far narrower than requested: in other runs, 7 (normal) and 15 (skewed)
of 6,000 studies ended that way at `width = 0.4`, and none of 2,000 at
`width = 0.2`. A larger minimum pilot helps: with `m_0 = 50` the
coverages were 0.95 and 0.91 at `width = 0.2` and 0.92 and 0.89 at
`width = 0.4`. Use
[`seq_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
to study a planned study before collecting data, and prefer a narrow
target or a larger `m_0`.

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

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods*, *8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology*, *59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

## See also

[`seq_aipe_reg_coef_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef_pilot.md)
for the pilot size,
[`seq_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
to simulate the procedure before collecting data, and
[`slr_slope`](https://yelleknek.github.io/SMSD/reference/slr_slope.md)
for the point estimate.

[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.html)
in DMAR is the fixed sample size analog for when planning values are
available, and
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.html)
gives a fixed sample size interval.

Other seq_aipe_reg_coef:
[`seq_aipe_reg_coef_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef_pilot.md),
[`seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# How much does job satisfaction change per standard deviation of
# perceived autonomy (both standardized)? The aim is a 95% confidence
# interval for the slope that is no wider than 0.40.
set.seed(113)
m <- seq_aipe_reg_coef_pilot(width = 0.40, m_0 = 20)
autonomy <- rnorm(m)
satisfaction <- 0.5 * autonomy + rnorm(m, sd = sqrt(0.75))
out <- seq_aipe_reg_coef(autonomy, satisfaction, width = 0.40)
out
#> Sequential check (simple linear regression slope)
#> 
#>   Current n       : 20
#>   Criterion n     : 103
#>   Target width    : 0.4
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : -0.2773156
#>   reg_coef          : 0.1640938
#>   upper_limit       : 0.6055032
#>   xi_squared        : 1.014418
#>   criterion_n       : 103
#>   current_n         : 20
#> 
#>   Stopping rule met: FALSE

# The rule is not met at the pilot. Add one employee at a time and check
# again, until it is.
while (!attr(out, "is_satisfied")) {
  new_autonomy <- rnorm(1)
  autonomy <- c(autonomy, new_autonomy)
  satisfaction <- c(satisfaction,
                    0.5 * new_autonomy + rnorm(1, sd = sqrt(0.75)))
  out <- seq_aipe_reg_coef(autonomy, satisfaction, width = 0.40)
}
out
#> Sequential check (simple linear regression slope)
#> 
#>   Current n       : 89
#>   Criterion n     : 89
#>   Target width    : 0.4
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : 0.1485017
#>   reg_coef          : 0.347154
#>   upper_limit       : 0.5458064
#>   xi_squared        : 0.9142845
#>   criterion_n       : 89
#>   current_n         : 89
#> 
#>   Stopping rule met: TRUE
```
