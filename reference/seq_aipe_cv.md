# Sequential AIPE Check for the Coefficient of Variation

Checks, with the data collected so far, whether enough observations are
in hand for a confidence interval for the population coefficient of
variation \\\kappa = \sigma / \mu\\ (the standard deviation relative to
the mean) that is no wider than a chosen width, without assuming a
distribution and without a planning value for \\\kappa\\. This is the
sequential accuracy in parameter estimation (AIPE) procedure of Kelley,
Darku, and Chattopadhyay (2018).

## Usage

``` r
seq_aipe_cv(x, width, conf_level = 0.95)
```

## Arguments

- x:

  A numeric vector: the observations collected so far, at least 4 of
  them. Missing or infinite values are an error.

- width:

  The target full width of the confidence interval for \\\kappa\\
  (\\\omega\\ in the cited article): the upper limit minus the lower
  limit, on the scale of the coefficient of variation as a proportion
  (not a percentage). A single positive number.

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

- `cv`:

  the estimated coefficient of variation, as a proportion (`NA` when the
  sample mean is zero, negative, or positive only by rounding error);

- `upper_limit`:

  the upper confidence limit (`NA` whenever `lower_limit` is);

- `xi_squared`:

  the estimated asymptotic variance \\\hat{\xi}\_n^2\\, at least
  \\n^{-3}\\ (the bound itself when the estimate is not yet usable);

- `criterion_n`:

  the sample size the stopping rule currently requires (at least
  `current_n + 1` when the variance estimate is not yet usable or the
  estimate cannot be computed);

- `current_n`:

  the current sample size.

The object also carries these attributes:

- `is_satisfied`:

  `TRUE` when the stopping rule is met, otherwise `FALSE`; never `NA`,
  and always identical to `current_n >= criterion_n`;

- `method`:

  the string `"coefficient of variation"`;

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
[`seq_aipe_cv_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv_pilot.md)),
calls `seq_aipe_cv()`, and, while the rule is not met, adds one or more
observations and calls it again. The function prints nothing and keeps
no state.

**Estimate.** The coefficient of variation is estimated by \\k_n = s_n /
\bar{X}\_n\\ (their Equation 11). It is defined for a population with a
positive mean, and is most meaningful for ratio scale measurements
(times, amounts, counts). Individual observations may be negative, but a
sample mean that is zero or negative cannot estimate \\\kappa\\; because
that can happen by chance in a small early sample, it is reported as
"not yet estimable" rather than as an error (see the section below). So
is a mean that differs from zero only by rounding error, as the mean of
centered or standardized scores does (the coefficient of variation is
not meaningful for such scores); the same test is applied by
[`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md).

**Stopping rule.** Sampling stops at the first \\n\\, at or after the
pilot size, with \$\$n \ge \left\lceil \frac{4
z\_{\alpha/2}^2}{\omega^2} \left(\hat{\xi}\_n^2 + \frac{1}{n}\right)
\right\rceil\$\$ (their Equation 22), where \\\hat{\xi}\_n^2\\ is the
distribution free estimator of the asymptotic variance of
\\\sqrt{n}(k_n - \kappa)\\ in their Equations 31 and 32, built from the
sample mean and variance and unbiased estimates of the third and fourth
central moments. The reported interval is \\k_n \pm z\_{\alpha/2}
\hat{\xi}\_n / \sqrt{n}\\ (their Equation 33). Whenever the rule is met,
the width of that interval is at most `width`.

Chattopadhyay and Kelley (2016) treat the same estimator under a minimum
risk criterion (see
[`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)).
The variance in their Equations 9 and 22 agrees with the one used here
for normal data and differs otherwise; the 2018 form used here is the
distribution free variance.

**What the guarantee is.** The width guarantee is exact. The coverage
guarantee is asymptotic: coverage approaches `conf_level` as `width`
shrinks (their Theorem 5). At the settings of the article's Table 4
(\\\kappa = 0.2\\, `width = 0.04`, 95%, one observation added at each
stage) coverage is close to nominal: 0.936 for normal data and 0.933 for
gamma data in 5,000 simulated studies each (the article reports 0.942
and 0.934), with about 240 observations at the end. At targets that are
wide relative to \\\kappa\\ the final sample is small and coverage falls
short. With \\\kappa = 0.5\\, coverage was 0.88 for gamma data and 0.79
for lognormal data at `width = 0.2` (10,000 studies each, about 52 and
62 observations at the end), and 0.79 and 0.69 at `width = 0.4` (6,000
studies each, about 16 observations). The variance estimate was zero or
negative in at most 2% of the pilot samples in these conditions, and
positive but below \\n^{-3}\\, so that the bound was used and the rule
met with an interval far narrower than requested, in 12 and 20 of the
6,000 studies at `width = 0.4` and in none of 2,000 at `width = 0.2`;
the shortfall comes from the small final sample rather than from the
bound. Use
[`seq_aipe_cv_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
to study a planned study before collecting data.

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

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Kelley, K. (2007). Sample size planning for the coefficient of variation
from the accuracy in parameter estimation approach. *Behavior Research
Methods*, *39*(4), 755–766.
[doi:10.3758/BF03192966](https://doi.org/10.3758/BF03192966)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology*, *59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

## See also

[`seq_aipe_cv_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv_pilot.md)
for the pilot size,
[`seq_aipe_cv_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
to simulate the procedure before collecting data, and
[`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
for the minimum risk point estimation counterpart. In DMAR:
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.html),
the fixed sample size analog for when a planning value of \\\kappa\\ is
available,
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.html), and
[`cv`](https://yelleknek.github.io/DMAR/reference/cv.html).

Other seq_aipe_cv:
[`seq_aipe_cv_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv_pilot.md),
[`seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# How variable are response times (in milliseconds) on a decision task,
# relative to their mean? The aim is a 95% confidence interval for the
# coefficient of variation that is no wider than 0.15.
set.seed(113)
m <- seq_aipe_cv_pilot(width = 0.15)
response_time <- rgamma(m, shape = 6.25, scale = 80)
out <- seq_aipe_cv(response_time, width = 0.15)
out
#> Sequential check (coefficient of variation)
#> 
#>   Current n       : 27
#>   Criterion n     : 59
#>   Target width    : 0.15
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : 0.3024828
#>   cv                : 0.3862728
#>   upper_limit       : 0.4700628
#>   xi_squared        : 0.04934597
#>   criterion_n       : 59
#>   current_n         : 27
#> 
#>   Stopping rule met: FALSE

# The rule is not met at the pilot. Add one participant at a time and
# check again, until it is.
while (!attr(out, "is_satisfied")) {
  response_time <- c(response_time, rgamma(1, shape = 6.25, scale = 80))
  out <- seq_aipe_cv(response_time, width = 0.15)
}
out
#> Sequential check (coefficient of variation)
#> 
#>   Current n       : 53
#>   Criterion n     : 53
#>   Target width    : 0.15
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : 0.3167059
#>   cv                : 0.3814155
#>   upper_limit       : 0.4461252
#>   xi_squared        : 0.05777204
#>   criterion_n       : 53
#>   current_n         : 53
#> 
#>   Stopping rule met: TRUE
```
