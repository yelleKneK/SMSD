# Sequential AIPE Check for the Standardized Mean

Checks, with the data collected so far, whether enough observations are
in hand for a confidence interval for the population standardized mean
\\\mu/\sigma\\ (the mean in standard deviation units; for difference
scores, the standardized mean change) that is no wider than a chosen
width, without assuming normality and without a planning value. The
standardized mean is one of the effect sizes covered by the sequential
accuracy in parameter estimation (AIPE) procedure of Kelley, Darku, and
Chattopadhyay (2018, Example 3).

## Usage

``` r
seq_aipe_sm(x, width, conf_level = 0.95)
```

## Arguments

- x:

  A numeric vector: the observations collected so far, at least 4 of
  them. Missing or infinite values are an error.

- width:

  The target full width of the confidence interval for \\\mu/\sigma\\
  (\\\omega\\ in the cited article): the upper limit minus the lower
  limit, in standard deviation units. A single positive number.

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

- `sm`:

  the estimated standardized mean (`NA` when the sample has no
  variation);

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

  the string `"standardized mean"`;

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
[`seq_aipe_sm_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm_pilot.md)),
calls `seq_aipe_sm()`, and, while the rule is not met, adds one or more
observations and calls it again. The function prints nothing and keeps
no state.

**Estimate.** The standardized mean is estimated by \\\bar{X}\_n /
s_n\\. Any real mean is allowed (the estimate is negative when the
sample mean is); only a sample with no variation at all cannot be used.

**Stopping rule.** Sampling stops at the first \\n\\, at or after the
pilot size, with \$\$n \ge \left\lceil \frac{4
z\_{\alpha/2}^2}{\omega^2} \left(\hat{\xi}\_n^2 + \frac{1}{n}\right)
\right\rceil\$\$ (their Equation 22), and the reported interval is
\\\bar{X}\_n / s_n \pm z\_{\alpha/2} \hat{\xi}\_n / \sqrt{n}\\ (their
Equation 23). Whenever the rule is met, the width of that interval is at
most `width`. The article names the standardized mean as a member of its
class of effect sizes but does not print its asymptotic variance. It
follows from the article's Theorem 1 by the delta method: \$\$\xi^2 =
1 - \frac{\mu \mu_3}{\sigma^4} + \frac{\mu^2 (\mu_4 - \sigma^4)}{4
\sigma^6},\$\$ where \\\mu_3\\ and \\\mu_4\\ are the third and fourth
central moments (for normal data this is \\1 + (\mu/\sigma)^2/2\\). The
estimate \\\hat{\xi}\_n^2\\ substitutes the sample mean and variance and
the unbiased central moment estimators of the article's Equations 78 and
79.

**What the guarantee is.** The width guarantee is exact. The coverage
guarantee is asymptotic: coverage approaches `conf_level` as `width`
shrinks (their Theorem 2). In simulated studies with a 95% interval and
one observation added at each stage, coverage for normal data with
\\\mu/\sigma = 0.4\\ was 0.949 at `width = 0.4` (6,000 studies, about
105 observations at the end) and 0.956 at `width = 0.2` (2,000 studies,
about 415 observations). For right skewed data (a standardized chi
square distribution with 3 degrees of freedom, shifted so that
\\\mu/\sigma = 1\\) it was 0.893 at `width = 0.4` and 0.940 at
`width = 0.2`: with skewed data the variance estimate runs low in small
samples, so the procedure stops early on average (at 62 observations
where the known variance would require 84, and at 289 where it would
require 334). The variance estimate was zero or negative in fewer than
1% of the pilot samples, and positive but below \\n^{-3}\\, so that the
bound was used and the rule met with an interval far narrower than
requested, in 2 of the 6,000 studies with skewed data at `width = 0.4`
and in no other. Use
[`seq_aipe_sm_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
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

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology*, *59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

## See also

[`seq_aipe_sm_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm_pilot.md)
for the pilot size and
[`seq_aipe_sm_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
to simulate the procedure before collecting data. In DMAR:
[`ss_aipe_sm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sm.html),
the fixed sample size analog for when a planning value is available, and
[`ci_sm`](https://yelleknek.github.io/DMAR/reference/ci_sm.html).

Other seq_aipe_sm:
[`seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md),
[`seq_aipe_sm_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# By how many standard deviations does a reading score change from fall
# to spring? The data are the students' change scores, and the aim is a
# 95% confidence interval for the standardized mean change that is no
# wider than 0.50.
set.seed(113)
m <- seq_aipe_sm_pilot(width = 0.50)
change <- rnorm(m, mean = 4, sd = 10)
out <- seq_aipe_sm(change, width = 0.50)
out
#> Sequential check (standardized mean)
#> 
#>   Current n       : 8
#>   Criterion n     : 69
#>   Target width    : 0.5
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : -0.6051389
#>   sm                : 0.0815806
#>   upper_limit       : 0.7683001
#>   xi_squared        : 0.9820929
#>   criterion_n       : 69
#>   current_n         : 8
#> 
#>   Stopping rule met: FALSE

# The rule is not met at the pilot. Add one student at a time and check
# again, until it is.
while (!attr(out, "is_satisfied")) {
  change <- c(change, rnorm(1, mean = 4, sd = 10))
  out <- seq_aipe_sm(change, width = 0.50)
}
out
#> Sequential check (standardized mean)
#> 
#>   Current n       : 75
#>   Criterion n     : 75
#>   Target width    : 0.5
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : 0.200412
#>   sm                : 0.4477031
#>   upper_limit       : 0.6949941
#>   xi_squared        : 1.193938
#>   criterion_n       : 75
#>   current_n         : 75
#> 
#>   Stopping rule met: TRUE
```
