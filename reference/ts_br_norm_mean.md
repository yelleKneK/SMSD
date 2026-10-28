# Two Stage Bounded Risk Point Estimation of a Normal Mean

Evaluates the two stage bounded risk point estimation procedure for the
mean \\\mu\\ of a normal distribution whose variance is unknown. The
loss in estimating \\\mu\\ by the sample mean is \\A\\\|\bar{X} -
\mu\|^k\\ (absolute error for \\k = 1\\, squared error for \\k = 2\\),
and the goal is a sample size at which the expected loss, the risk, is
at most `w`. The first stage (`pilot_data`) fixes the total sample size;
the second stage (`additional_data`) completes it. The estimate is the
mean of all observations from both stages. The procedure is risk based,
so there is no confidence level.

## Usage

``` r
ts_br_norm_mean(pilot_data, additional_data = NULL, A, k, w, na.rm = TRUE)
```

## Arguments

- pilot_data:

  Numeric vector of the first stage observations, at least \\k + 2\\.
  The whole first stage must be supplied at every call, because it alone
  determines the total sample size.

- additional_data:

  Numeric vector of the second stage observations collected so far, or
  `NULL` (the default) when none have been collected.

- A:

  The loss constant: the loss is \\A\\\|\bar{X} - \mu\|^k\\. A single
  positive number, in the units of `w` per unit of error to the power
  `k`.

- k:

  The error index: 1 for absolute error loss, 2 for squared error loss.

- w:

  The risk bound: the largest acceptable expected loss. A single
  positive number.

- na.rm:

  Logical. If `TRUE` (the default), missing values are removed from both
  stages before anything is computed; if `FALSE`, a missing value is an
  error.

## Value

An object of class `c("smsd_check", "data.frame")` (see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
for its [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and `tidy()`
methods) with columns `term` and `value` and four rows:

- `mean`:

  The mean of all observations supplied, the estimate of \\\mu\\.

- `risk`:

  The estimated risk of that mean at the current sample size, \\B\\(S^2
  / n)^{k/2}\\, with \\S^2\\ the sample variance of all observations
  supplied. A descriptive estimate; the guarantee is the one stated in
  Details.

- `optimal_n`:

  The total sample size \\N\\ fixed by the first stage.

- `current_n`:

  The number of observations supplied, both stages together.

Attributes: `is_satisfied` (`TRUE` once `current_n` is at least
`optimal_n`), `current_n`, `criterion_n` (equal to `optimal_n`),
`loss_A`, `k`, `risk_bound`, `method`, and `notes` (always empty for
this function). There is no `conf_level` or `width`, and no
`variance_unusable` attribute: a first stage whose variance cannot be
estimated is an error (see Details).

## Details

If the standard deviation \\\sigma\\ were known, the risk of the mean of
\\n\\ observations would be \\B\\(\sigma^2 / n)^{k/2}\\ with \\B =
A\\2^{k/2}\\\Gamma((k + 1)/2)/\sqrt{\pi}\\ (which is \\A\\ for \\k = 2\\
and \\A\sqrt{2/\pi}\\ for \\k = 1\\), and the smallest sample size
meeting the bound would be \\n^\* = (B / w)^{2/k} \sigma^2\\. With \\m\\
first stage observations and their sample variance \\S_m^2\\, the total
sample size is \$\$N = \max\\m,\\ \lceil (b_m B / w)^{2/k} \\ S_m^2
\rceil\\, \qquad b_m = \left(\frac{m - 1}{2}\right)^{k/2}
\frac{\Gamma((m - k - 1)/2)}{\Gamma((m - 1)/2)}.\$\$ If \\N = m\\ the
first stage is enough; otherwise \\N - m\\ further observations are
collected, once. \\N\\ depends on the first stage alone. The row
`optimal_n` is this total \\N\\, not the number still to collect.

The constant \\b_m\\ is the expectation of \\(S_m^2 /
\sigma^2)^{-k/2}\\, which exists only when \\m \ge k + 2\\; for \\k =
2\\ it is \\(m - 1)/(m - 3)\\. A smaller first stage is an error.

**The guarantee is exact.** For normal data the risk of the final mean
is at most `w` for every \\\mu\\, every \\\sigma\\, and every first
stage size \\m \ge k + 2\\: \\N\\ depends only on \\S_m^2\\, which is
independent of the sample means, so the risk is \\B \sigma^k
E\[N^{-k/2}\]\\, and \\N \ge (b_m B / w)^{2/k} S_m^2\\ bounds that by
`w`. The price is sample size: the expected \\N\\ is about \\b_m^{2/k}\\
times \\n^\*\\ (for \\k = 2\\, \\(m - 1)/(m - 3)\\), which is why
[`ts_br_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean_pilot.md)
does not return a very small first stage. The guarantee rests on
normality.

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

Data on a scale so large or so small that a variance cannot be
represented in double precision (observations beyond about \\10^{154}\\,
or differing by less than about \\10^{-154}\\) are an error that says to
express them in other units; no row of the result is ever infinite.

## References

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`ts_br_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean_pilot.md)
for the first stage size;
[`seq_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md)
for the purely sequential procedure.

Other ts_br_norm_mean:
[`ts_br_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# A clinic manager wants the mean waiting time (minutes) estimated with
# a mean squared error of at most 4, that is, squared error loss with
# A = 1 and risk bound w = 4.
set.seed(113)
m <- ts_br_norm_mean_pilot(A = 1, k = 2, w = 4)
first_stage <- rnorm(m, mean = 45, sd = 12)
check <- ts_br_norm_mean(first_stage, A = 1, k = 2, w = 4)
check
#> Sequential check (two stage bounded risk estimation (normal mean))
#> 
#>   Current n       : 10
#>   Criterion n     : 43
#> 
#>   mean              : 41.24481
#>   risk              : 13.11412
#>   optimal_n         : 43
#>   current_n         : 10
#> 
#>   Stopping rule met: FALSE

# The first stage fixed the total sample size; collect the remainder.
n_more <- attr(check, "criterion_n") - attr(check, "current_n")
second_stage <- rnorm(n_more, mean = 45, sd = 12)
ts_br_norm_mean(first_stage, second_stage, A = 1, k = 2, w = 4)
#> Sequential check (two stage bounded risk estimation (normal mean))
#> 
#>   Current n       : 43
#>   Criterion n     : 43
#> 
#>   mean              : 47.41062
#>   risk              : 3.74368
#>   optimal_n         : 43
#>   current_n         : 43
#> 
#>   Stopping rule met: TRUE
```
