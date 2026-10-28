# Two Stage Fixed Width Confidence Interval for a Distribution Free Mean

Evaluates the two stage procedure for a confidence interval of fixed
width for a population mean \\\mu\\ when no distribution is assumed
(only a finite variance). The first stage (`pilot_data`) fixes the total
sample size; the second stage (`additional_data`) completes it. The
interval is \\\bar{X}\_N \pm d\\, where `d` is the half width (the
margin of error) and \\\bar{X}\_N\\ is the mean of all observations from
both stages.

## Usage

``` r
ts_fwd_dfr_mean(
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

  The nominal confidence level, strictly between 0 and 1. Default 0.95.
  The level is asymptotic; see Details for the size of the shortfall.

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
the total sample size is \$\$N = \max\\m,\\ \lceil z\_{\alpha/2}^2 \\
S_m^2 / d^2 \rceil\\,\$\$ where \\z\_{\alpha/2}\\ is the upper
\\\alpha/2\\ point of the standard normal distribution and \\\alpha =
1 - \\ `conf_level`. If \\N = m\\ the first stage is enough; otherwise
\\N - m\\ further observations are collected, once. \\N\\ depends on the
first stage alone. The row `optimal_n` is this total \\N\\, not the
number still to collect.

**The confidence level is asymptotic, and the shortfall can be large.**
The rule plugs the first stage variance into the known variance sample
size \\C = z\_{\alpha/2}^2 \sigma^2 / d^2\\ with no allowance for having
estimated it. Once \\C\\ is more than about 1.2 times the first stage,
the realized coverage is therefore below the nominal level, and it
approaches that level only as the first stage grows. (The figure 1.2 is
for normal data, a first stage of 30, and a nominal 0.95; it is between
about 1.1 and 1.3 for first stages of 5 to 200 and nominal levels of
0.90 to 0.99. When the first stage is nearly as large as \\C\\, or
larger, more observations are taken than the target calls for and
coverage is at or above nominal.) For normal data and a narrow target
the coverage is \\P(\|t\_{m-1}\| \le z\_{\alpha/2})\\: at a nominal
0.95, 0.918 for a first stage of 10, 0.940 for 30, and 0.947 for 100.
For skewed data it is lower. With the first stage that
[`ts_fwd_dfr_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean_pilot.md)
returns by default (30 observations unless the target is very narrow),
the coverage at a nominal 0.95 for populations with standard deviation 2
is:

|                                  |       |       |       |       |
|----------------------------------|-------|-------|-------|-------|
| Half width `d`                   | 1     | 0.5   | 0.2   | 0.1   |
| Known variance sample size \\C\\ | 15.4  | 61.5  | 384   | 1537  |
| Normal data                      | 0.994 | 0.942 | 0.941 | 0.940 |
| Exponential data (skewness 2)    | 0.996 | 0.919 | 0.915 | 0.917 |
| Lognormal data (skewness 4)      | 0.998 | 0.928 | 0.874 | 0.876 |

The normal row is exact; the other rows are from 20,000 simulated
studies per cell (standard error at most 0.0023). At `d = 1` the first
stage of 30 already exceeds the sample size needed, which is why
coverage there is above nominal. Larger first stages help slowly: with
100 first stage observations the narrow target coverages are 0.947
(normal), 0.936 (exponential), and 0.918 (lognormal). A first stage of
only two observations gives 0.78 for normal data at `d = 1`.

If the data can be taken as normal, use
[`ts_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md)
instead: Stein's procedure replaces \\z\_{\alpha/2}\\ by a *t* quantile
and its coverage is exact for any first stage size. If they cannot,
consider the purely sequential
[`seq_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md),
which keeps updating the variance estimate as observations arrive (its
level is asymptotic too). This two stage rule is for studies in which
only one further round of data collection is possible; give it as large
a first stage as the study allows.

The interval rows are reported at every call so that progress can be
followed, but the interval carries even its asymptotic confidence level
only when the stopping rule is met, that is, when the current sample
size has reached `optimal_n`.

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

## See also

[`ts_fwd_dfr_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean_pilot.md)
for the first stage size;
[`ts_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md)
for the exact procedure under normality;
[`seq_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md)
for the purely sequential procedure.

Other ts_fwd_dfr_mean:
[`ts_fwd_dfr_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# An HR analyst wants mean weekly hours of unpaid overtime to within
# half an hour, with a nominal 95% confidence level; the hours are
# right skewed, so normality is not assumed.
set.seed(113)
m <- ts_fwd_dfr_mean_pilot(d = 0.5, conf_level = 0.95)
first_stage <- rgamma(m, shape = 2, rate = 0.5)
check <- ts_fwd_dfr_mean(first_stage, d = 0.5)
check
#> Sequential check (two stage fixed width interval (distribution free mean))
#> 
#>   Current n       : 30
#>   Criterion n     : 91
#>   Target width    : 1
#>   Confidence level: 95%
#> 
#>   lower_limit       : 2.652349
#>   mean              : 3.152349
#>   upper_limit       : 3.652349
#>   optimal_n         : 91
#>   current_n         : 30
#> 
#>   Stopping rule met: FALSE

# The first stage fixed the total sample size; collect the remainder.
n_more <- attr(check, "criterion_n") - attr(check, "current_n")
second_stage <- rgamma(n_more, shape = 2, rate = 0.5)
ts_fwd_dfr_mean(first_stage, second_stage, d = 0.5)
#> Sequential check (two stage fixed width interval (distribution free mean))
#> 
#>   Current n       : 91
#>   Criterion n     : 91
#>   Target width    : 1
#>   Confidence level: 95%
#> 
#>   lower_limit       : 2.986092
#>   mean              : 3.486092
#>   upper_limit       : 3.986092
#>   optimal_n         : 91
#>   current_n         : 91
#> 
#>   Stopping rule met: TRUE
```
