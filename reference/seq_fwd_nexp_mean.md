# Fixed Width Confidence Interval for the Location of a Negative Exponential Distribution

Given the data collected so far, checks the stopping rule of the purely
sequential fixed width confidence interval for the location parameter of
a negative exponential distribution (Mukhopadhyay and de Silva, 2009).
The negative exponential (two parameter exponential) distribution has
density \\\sigma^{-1}\exp\\-(x - \theta)/\sigma\\\\ for \\x \> \theta\\:
\\\theta\\ is the location, the smallest value the variable can take (a
minimum response time, a guarantee period), and \\\sigma\\ is the scale.
Despite the function's name, the quantity estimated is the location
\\\theta\\, not the mean (which is \\\theta + \sigma\\).

## Usage

``` r
seq_fwd_nexp_mean(x, d, conf_level = 0.95, na.rm = TRUE)
```

## Arguments

- x:

  Numeric vector of the observations collected so far (at least 2).

- d:

  The full width of the confidence interval, in the units of the data:
  the interval reported runs from the sample minimum minus `d` up to the
  sample minimum. A positive number. (In the fixed width procedures for
  a mean, `d` is a half width; here the interval is one-sided, so `d` is
  its whole width.)

- conf_level:

  The confidence level (default 0.95).

- na.rm:

  Logical; should `NA` values be removed first? Default `TRUE`. With
  `FALSE`, missing values are an error.

## Value

An object of class `c("smsd_check", "data.frame")` with columns `term`
and `value` and these rows:

- `lower_limit`:

  The sample minimum minus `d`.

- `location`:

  The sample minimum, the point estimate of the location \\\theta\\.

- `upper_limit`:

  The sample minimum.

- `optimal_n`:

  The sample size the current scale estimate calls for, \\\lceil
  \log(1/\alpha)\\ U_n / d \rceil\\ (`NA` when the data are constant).

- `current_n`:

  The number of observations used.

Attributes: `is_satisfied` (`TRUE` once `current_n` is at least
`optimal_n`), `current_n`, `criterion_n` (equal to `optimal_n`),
`conf_level`, `width` (the full width `d`), `method`,
`variance_unusable` (`TRUE` when the data are constant, so that no scale
has been estimated; the attribute has the same name in every procedure
of the package), and `notes`.

## Details

The function does not collect data and keeps no state. Call it after
each new observation (or batch), starting once
[`seq_fwd_nexp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean_pilot.md)
observations are in hand, and stop sampling the first time
`attr(out, "is_satisfied")` is `TRUE`.

The location is estimated by the sample minimum \\X\_{(1)}\\, which can
only overestimate it, so the interval is \\\[X\_{(1)} - d,\\
X\_{(1)}\]\\, of width `d`. Because \\X\_{(1)} - \theta\\ is exponential
with scale \\\sigma/n\\, the interval covers \\\theta\\ with probability
\\1 - \exp(-n d/\sigma)\\, which is at least \\1 - \alpha\\ when \\n \ge
n^\* = \log(1/\alpha)\\\sigma/d\\. The sequential rule replaces
\\\sigma\\ by \\U_n = \sum (X_i - X\_{(1)})/(n - 1)\\: sampling stops at
the first \\n\\ (at or beyond the pilot size) with \$\$n \ge
\log(1/\alpha)\\ U_n / d .\$\$ The `optimal_n` row is the right hand
side rounded up, evaluated at the current data.

**The confidence level is asymptotic.** It is attained in the limit as
`d` shrinks, not at every `d`; for wide targets sampling tends to stop
early on an underestimate of the scale. In simulation (scale 2, nominal
0.95, the pilot of
[`seq_fwd_nexp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean_pilot.md)
with its defaults, 40,000 replications, Monte Carlo standard error
0.001), coverage was 0.938 at `d = 0.4` (where the required sample size
\\n^\*\\ is about 15). It was lowest, about 0.930, when \\n^\*\\ was 25,
and rose from there: 0.941 at \\n^\*\\ of 40, 0.944 at `d = 0.1`
(\\n^\*\\ about 60), and 0.948 at `d = 0.025` (\\n^\*\\ about 240;
400,000 replications).
[`ts_fwd_nexp_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md)
is the two stage alternative, whose confidence level is guaranteed for
negative exponential data.

**Data.** Any finite values are accepted, because the location may be
negative. The model is right skewed with a sharp lower endpoint; the
function cannot check that.

**Constant data.** If every observation is equal, or equal up to
rounding error (a range of no more than \\10^{-10}\\ times the largest
absolute value), \\U_n\\ is zero, or nothing but rounding error, and the
rule as printed would be met at once. The function departs from the
printed rule deliberately: no scale has been estimated, so the rule is
reported as not met, `optimal_n` is `NA`,
`attr(out, "variance_unusable")` is `TRUE`, and a note says so (see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md));
data whose real spread is that small beside their size should have a
common value subtracted first. The limits are still reported (the sample
minimum minus `d`, and the sample minimum, which involve no estimate of
the scale), but they are then not an interval to report:
[`summary()`](https://rdrr.io/r/base/summary.html), `glance()`, and
`tidy()` return them without the note, so test
`attr(out, "variance_unusable")` before using them.

**Extreme units.** Data that do vary, but whose spread is at the very
edge of double precision, so that the scale estimate cannot be
represented (differences among the observations near \\10^{308}\\ or
near \\10^{-323}\\), are an error that says to change the units. That
case is not reported as a rule not yet met, because collecting more data
cannot help.

## References

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_fwd_nexp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean_pilot.md);
[`ts_fwd_nexp_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md)
for the two stage procedure;
[`seq_mr_nexp_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean.md)
for minimum risk point estimation of the same location.

Other seq_fwd_nexp_mean:
[`seq_fwd_nexp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Shortest possible response time on a reaction task (seconds), to
# within 0.05 seconds with 95% confidence. Times are simulated as a
# minimum of 0.20 seconds plus an exponential delay averaging 0.30.
set.seed(113)
m <- seq_fwd_nexp_mean_pilot(d = 0.05, conf_level = 0.95)
times <- 0.20 + rexp(m, rate = 1 / 0.30)
out <- seq_fwd_nexp_mean(times, d = 0.05, conf_level = 0.95)
out
#> Sequential check (fixed width interval (negative exponential location))
#> 
#>   Current n       : 10
#>   Criterion n     : 19
#>   Target width    : 0.05
#>   Confidence level: 95%
#> 
#>   lower_limit       : 0.1818261
#>   location          : 0.2318261
#>   upper_limit       : 0.2318261
#>   optimal_n         : 19
#>   current_n         : 10
#> 
#>   Stopping rule met: FALSE

# Add one trial at a time until the stopping rule is met.
while (!attr(out, "is_satisfied")) {
  times <- c(times, 0.20 + rexp(1, rate = 1 / 0.30))
  out <- seq_fwd_nexp_mean(times, d = 0.05, conf_level = 0.95)
}
out
#> Sequential check (fixed width interval (negative exponential location))
#> 
#>   Current n       : 17
#>   Criterion n     : 16
#>   Target width    : 0.05
#>   Confidence level: 95%
#> 
#>   lower_limit       : 0.1818261
#>   location          : 0.2318261
#>   upper_limit       : 0.2318261
#>   optimal_n         : 16
#>   current_n         : 17
#> 
#>   Stopping rule met: TRUE
```
