# Two Stage Fixed Width Confidence Interval for a Normal Mean

Evaluates Stein's two stage procedure for a confidence interval of fixed
width for the mean \\\mu\\ of a normal distribution whose variance is
unknown. The first stage (`pilot_data`) fixes the total sample size; the
second stage (`additional_data`) completes it. The interval is
\\\bar{X}\_N \pm d\\, where `d` is the half width (the margin of error)
and \\\bar{X}\_N\\ is the mean of all observations from both stages.

## Usage

``` r
ts_fwd_norm_mean(
  pilot_data,
  additional_data = NULL,
  d,
  conf_level = 0.95,
  na.rm = TRUE
)
```

## Arguments

- pilot_data:

  Numeric vector of the first stage observations, at least two. The
  whole first stage must be supplied at every call, because it alone
  determines the total sample size.

- additional_data:

  Numeric vector of the second stage observations collected so far, or
  `NULL` (the default) when none have been collected.

- d:

  The half width of the confidence interval (the margin of error), in
  the units of the data. A single positive number. The full width is
  \\2d\\.

- conf_level:

  The confidence level, strictly between 0 and 1. Default 0.95.

- na.rm:

  Logical. If `TRUE` (the default), missing values are removed from both
  stages before anything is computed; if `FALSE`, a missing value is an
  error.

## Value

An object of class `c("smsd_check", "data.frame")` (see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
for its [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and `tidy()`
methods) with columns `term` and `value` and five rows:

- `lower_limit`:

  \\\bar{X} - d\\, with \\\bar{X}\\ the mean of all observations
  supplied.

- `mean`:

  \\\bar{X}\\, the estimate of \\\mu\\.

- `upper_limit`:

  \\\bar{X} + d\\.

- `optimal_n`:

  The total sample size \\N\\ fixed by the first stage.

- `current_n`:

  The number of observations supplied, both stages together.

Attributes: `is_satisfied` (`TRUE` once `current_n` is at least
`optimal_n`), `current_n`, `criterion_n` (equal to `optimal_n`),
`conf_level`, `width` (the full width \\2d\\), `method`, and `notes`
(always empty for this function). There is no `variance_unusable`
attribute: a first stage whose variance cannot be estimated is an error
(see Details).

## Details

With \\m\\ first stage observations and their sample variance \\S_m^2\\,
the total sample size is \$\$N = \max\\m,\\ \lceil t\_{m-1,\alpha/2}^2
\\ S_m^2 / d^2 \rceil\\,\$\$ where \\t\_{m-1,\alpha/2}\\ is the upper
\\\alpha/2\\ point of the *t* distribution with \\m - 1\\ degrees of
freedom and \\\alpha = 1 - \\ `conf_level`. If \\N = m\\ the first stage
is enough; otherwise \\N - m\\ further observations are collected, once.
\\N\\ depends on the first stage alone and does not change as second
stage observations arrive. The row `optimal_n` is this total \\N\\, not
the number still to collect.

**The guarantee is exact.** For normal data the interval \\\bar{X}\_N
\pm d\\ covers \\\mu\\ with probability at least `conf_level` for every
\\\mu\\, every \\\sigma\\, every `d`, and every first stage size \\m \ge
2\\ (Stein, 1945; Mukhopadhyay and de Silva, 2009, discuss the procedure
in their Section 6.2.1). The price is sample size: the expected \\N\\ is
about \\t\_{m-1,\alpha/2}^2 / z\_{\alpha/2}^2\\ times the known variance
sample size \\C = z\_{\alpha/2}^2 \sigma^2 / d^2\\, which is why
[`ts_fwd_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean_pilot.md)
does not return a very small first stage. The guarantee rests on
normality; for data that are not normal see
[`ts_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md),
whose level is only asymptotic.

The interval rows are reported at every call so that progress can be
followed, but the interval carries the stated confidence level only when
the stopping rule is met, that is, when the current sample size has
reached `optimal_n`.

**A first stage with no variation is an error.** If every first stage
observation is equal, or equal up to rounding error (a range of no more
than \\10^{-10}\\ times the largest absolute value, the rule used
throughout the package), the first stage variance cannot be estimated
and the function stops with an error that says what to do. The purely
sequential functions report such data as a stopping rule that is not yet
met, with the attribute `variance_unusable`; here that would mislead,
because the total sample size is fixed by the first stage alone, so no
later observation could turn a rule reported as not met into one that is
met. Data whose real variation is that small beside their size (times
recorded in seconds since 1970 that differ by hundredths of a second,
for example) should have a common value subtracted from both stages
first; the total sample size does not depend on it.

## References

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

Stein, C. (1945). A two-sample test for a linear hypothesis whose power
is independent of the variance. *The Annals of Mathematical Statistics*,
*16*(3), 243–258.
[doi:10.1214/aoms/1177731088](https://doi.org/10.1214/aoms/1177731088)

## See also

[`ts_fwd_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean_pilot.md)
for the first stage size;
[`seq_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md)
for the purely sequential procedure;
[`seq_fwd_norm_mean_optimal`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_optimal.md)
for the sample size when the variance is known;
[`ts_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md)
when normality is not assumed.

Other ts_fwd_norm_mean:
[`ts_fwd_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# A training director wants the mean score on a 100 point skills
# assessment to within 2 points, with 95% confidence.
set.seed(113)
m <- ts_fwd_norm_mean_pilot(d = 2, conf_level = 0.95)
first_stage <- rnorm(m, mean = 70, sd = 8)
check <- ts_fwd_norm_mean(first_stage, d = 2)
check
#> Sequential check (two stage fixed width interval (normal mean))
#> 
#>   Current n       : 10
#>   Criterion n     : 75
#>   Target width    : 4
#>   Confidence level: 95%
#> 
#>   lower_limit       : 65.49654
#>   mean              : 67.49654
#>   upper_limit       : 69.49654
#>   optimal_n         : 75
#>   current_n         : 10
#> 
#>   Stopping rule met: FALSE

# The first stage fixed the total sample size; collect the remainder.
n_more <- attr(check, "criterion_n") - attr(check, "current_n")
second_stage <- rnorm(n_more, mean = 70, sd = 8)
ts_fwd_norm_mean(first_stage, second_stage, d = 2)
#> Sequential check (two stage fixed width interval (normal mean))
#> 
#>   Current n       : 75
#>   Criterion n     : 75
#>   Target width    : 4
#>   Confidence level: 95%
#> 
#>   lower_limit       : 68.66762
#>   mean              : 70.66762
#>   upper_limit       : 72.66762
#>   optimal_n         : 75
#>   current_n         : 75
#> 
#>   Stopping rule met: TRUE
```
