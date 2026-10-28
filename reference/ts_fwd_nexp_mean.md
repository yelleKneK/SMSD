# Two Stage Fixed Width Interval for a Negative Exponential Location

Evaluates the two stage procedure for a confidence interval of fixed
width for the location \\\theta\\ of a negative exponential
distribution: the model \\X = \theta + \sigma E\\ with \\E\\ a standard
exponential variable, so that \\\theta\\ is the smallest value the
outcome can take (a threshold, a guarantee time, a minimum response
time) and \\\sigma\\ is the scale. Despite the function's name, the
quantity estimated is the location \\\theta\\, not the mean \\\theta +
\sigma\\. The first stage (`pilot_data`) fixes the total sample size;
the second stage (`additional_data`) completes it. The interval is
\\\[X\_{(1)} - d,\\ X\_{(1)}\]\\, where \\X\_{(1)}\\ is the smallest of
all observations from both stages and `d` is the full width.

## Usage

``` r
ts_fwd_nexp_mean(
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

  The full width of the confidence interval, in the units of the data. A
  single positive number.

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

  \\X\_{(1)} - d\\, with \\X\_{(1)}\\ the smallest of all observations
  supplied.

- `location`:

  \\X\_{(1)}\\, the estimate of \\\theta\\.

- `upper_limit`:

  \\X\_{(1)}\\; the interval is one-sided, so its upper limit is the
  estimate itself.

- `optimal_n`:

  The total sample size \\N\\ fixed by the first stage.

- `current_n`:

  The number of observations supplied, both stages together.

Attributes: `is_satisfied` (`TRUE` once `current_n` is at least
`optimal_n`), `current_n`, `criterion_n` (equal to `optimal_n`),
`conf_level`, `width` (the full width `d`), `method`, and `notes`
(always empty for this function). There is no `variance_unusable`
attribute: a first stage whose scale cannot be estimated is an error
(see Details).

## Details

With \\m\\ first stage observations, their minimum \\X\_{(1:m)}\\, and
the scale estimate \\U_m = \sum\_{i=1}^{m} (X_i - X\_{(1:m)}) / (m -
1)\\, the total sample size is \$\$N = \max\\m,\\ \lceil
F\_{2,2m-2,\alpha} \\ U_m / d \rceil\\,\$\$ where \\F\_{2,2m-2,\alpha}\\
is the upper \\\alpha\\ point of the *F* distribution with 2 and \\2m -
2\\ degrees of freedom and \\\alpha = 1 - \\ `conf_level`. If \\N = m\\
the first stage is enough; otherwise \\N - m\\ further observations are
collected, once. \\N\\ depends on the first stage alone. The row
`optimal_n` is this total \\N\\, not the number still to collect.

The quantile is one-sided because the interval is: the smallest
observation can never fall below \\\theta\\, so the interval misses only
when \\X\_{(1)} - \theta \> d\\.

**The guarantee is exact.** For negative exponential data the interval
covers \\\theta\\ with probability at least `conf_level` for every
\\\theta\\, every \\\sigma\\, every `d`, and every first stage size \\m
\ge 2\\. (Given \\N\\, the interval misses with probability \\\exp(-N d
/ \sigma)\\; \\N\\ depends only on \\U_m\\, with \\2(m - 1) U_m /
\sigma\\ distributed as chi square on \\2m - 2\\ degrees of freedom; and
the expectation of \\\exp(-F U_m / \sigma)\\ equals \\\alpha\\ exactly
at the *F* quantile used.) The price is sample size: the expected \\N\\
is about \\F\_{2,2m-2,\alpha} / \log(1/\alpha)\\ times the known scale
sample size \\C = \log(1/\alpha) \sigma / d\\, which is why
[`ts_fwd_nexp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean_pilot.md)
does not return a very small first stage. The guarantee rests on the
negative exponential model.

The interval rows are reported at every call so that progress can be
followed, but the interval carries the stated confidence level only when
the stopping rule is met, that is, when the current sample size has
reached `optimal_n`.

**A first stage with no variation is an error.** If every first stage
observation is equal, or equal up to rounding error (a range of no more
than \\10^{-10}\\ times the largest absolute value, the rule used
throughout the package), the scale cannot be estimated and the function
stops with an error that says what to do. The purely sequential
functions report such data as a stopping rule that is not yet met, with
the attribute `variance_unusable`; here that would mislead, because the
total sample size is fixed by the first stage alone, so no later
observation could turn a rule reported as not met into one that is met.
Data whose real variation is that small beside their size (times
recorded in seconds since 1970 that differ by hundredths of a second,
for example) should have a common value subtracted from both stages
first; the total sample size does not depend on it.

## References

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`ts_fwd_nexp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean_pilot.md)
for the first stage size;
[`seq_fwd_nexp_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean.md)
for the purely sequential procedure.

Other ts_fwd_nexp_mean:
[`ts_fwd_nexp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# A usability researcher models the time to complete a task as a
# shortest possible time plus an exponential delay, and wants the
# shortest possible time to within half a second, with 95% confidence.
set.seed(113)
m <- ts_fwd_nexp_mean_pilot(d = 0.5, conf_level = 0.95)
first_stage <- 3 + rexp(m, rate = 1 / 4)
check <- ts_fwd_nexp_mean(first_stage, d = 0.5)
check
#> Sequential check (two stage fixed width interval (negative exponential location))
#> 
#>   Current n       : 10
#>   Criterion n     : 29
#>   Target width    : 0.5
#>   Confidence level: 95%
#> 
#>   lower_limit       : 2.924348
#>   location          : 3.424348
#>   upper_limit       : 3.424348
#>   optimal_n         : 29
#>   current_n         : 10
#> 
#>   Stopping rule met: FALSE

# The first stage fixed the total sample size; collect the remainder.
n_more <- attr(check, "criterion_n") - attr(check, "current_n")
second_stage <- 3 + rexp(n_more, rate = 1 / 4)
ts_fwd_nexp_mean(first_stage, second_stage, d = 0.5)
#> Sequential check (two stage fixed width interval (negative exponential location))
#> 
#>   Current n       : 29
#>   Criterion n     : 29
#>   Target width    : 0.5
#>   Confidence level: 95%
#> 
#>   lower_limit       : 2.669607
#>   location          : 3.169607
#>   upper_limit       : 3.169607
#>   optimal_n         : 29
#>   current_n         : 29
#> 
#>   Stopping rule met: TRUE
```
