# Fixed Width Confidence Interval for a Normal Mean

Given the data collected so far, checks the stopping rule of the purely
sequential fixed width confidence interval for the mean of a normal
population with unknown variance (Mukhopadhyay and de Silva, 2009). The
aim is an interval of half width `d` (full width \\2d\\) around the
sample mean with the stated confidence level. No fixed sample size can
promise that, because the size needed depends on the unknown variance,
so observations are collected until the sample size reaches the size the
current variance estimate calls for.

## Usage

``` r
seq_fwd_norm_mean(x, d, conf_level = 0.95, na.rm = TRUE)
```

## Arguments

- x:

  Numeric vector of the observations collected so far (at least 2).

- d:

  The half width of the confidence interval, in the units of the data:
  the interval reported is the sample mean plus and minus `d`, so its
  full width is \\2d\\. A positive number.

- conf_level:

  The confidence level (default 0.95).

- na.rm:

  Logical; should `NA` values be removed first? Default `TRUE`. With
  `FALSE`, missing values are an error.

## Value

An object of class `c("smsd_check", "data.frame")` with columns `term`
and `value` and these rows:

- `lower_limit`:

  The sample mean minus `d`.

- `mean`:

  The sample mean, the point estimate.

- `upper_limit`:

  The sample mean plus `d`.

- `optimal_n`:

  The sample size the current variance estimate calls for, \\\lceil
  z\_{\alpha/2}^2 S_n^2 / d^2 \rceil\\ (`NA` when the data are
  constant).

- `current_n`:

  The number of observations used.

Attributes: `is_satisfied` (`TRUE` once `current_n` is at least
`optimal_n`), `current_n`, `criterion_n` (equal to `optimal_n`),
`conf_level`, `width` (the full width \\2d\\), `method`,
`variance_unusable` (`TRUE` when the data are constant, so that no
variance has been estimated), and `notes`.

## Details

The function does not collect data and keeps no state. Call it after
each new observation (or batch), starting once
[`seq_fwd_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_pilot.md)
observations are in hand, and stop sampling the first time
`attr(out, "is_satisfied")` is `TRUE`.

With the population standard deviation \\\sigma\\ known, the interval
\\\bar{X}\_n \pm d\\ has confidence level \\1 - \alpha\\ when \\n \ge
n^\* = z\_{\alpha/2}^2 \sigma^2 / d^2\\
([`seq_fwd_norm_mean_optimal`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_optimal.md)
returns that sample size, rounded up to a whole number, for a given
variance). The sequential rule replaces \\\sigma^2\\ by the sample
variance \\S_n^2\\: sampling stops at the first \\n\\ (at or beyond the
pilot size) with \$\$n \ge z\_{\alpha/2}^2 S_n^2 / d^2 .\$\$ The
`optimal_n` row is the right hand side rounded up, evaluated at the
current data. The reported interval is always \\\bar{X}\_n \pm d\\; it
has the stated confidence level only once the rule is met.

**The confidence level is asymptotic.** It is attained in the limit as
`d` shrinks (so that the required sample size grows), not at every `d`.
When the target is wide relative to the standard deviation, sampling
tends to stop early on an underestimate of the variance, and coverage
falls below the stated level. In simulation (normal data with standard
deviation 2, nominal 0.95, the pilot of
[`seq_fwd_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_pilot.md)
with its defaults, which was 10 at every target used, 40,000
replications, Monte Carlo standard error 0.001), coverage was 0.933 at
`d = 1` (where the required sample size \\n^\*\\ is about 15). It was
lowest, about 0.922, when \\n^\*\\ was between 24 and 28, and rose from
there: 0.930 at \\n^\*\\ of 40, 0.938 at `d = 0.5` (\\n^\*\\ about 61),
0.943 at \\n^\*\\ of 100, and 0.949 at `d = 0.25` (\\n^\*\\ about 246;
400,000 replications). The pilot matters: with `m_0 = 4` in place of the
default 10, coverage at `d = 1` was 0.880.
[`ts_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md)
is the two stage alternative, whose confidence level is guaranteed for
normal data at the price of a larger expected sample size.

**Constant data.** If every observation is equal, or equal up to
rounding error (a range of no more than \\10^{-10}\\ times the largest
absolute value, as with gain scores of 0.3 computed as posttest minus
pretest), the sample variance is zero, or nothing but rounding error,
and the rule as printed would be met at once. The function departs from
the printed rule deliberately: no variance has been estimated, so the
rule is reported as not met, `optimal_n` is `NA`,
`attr(out, "variance_unusable")` is `TRUE`, and a note says so (see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md));
data whose real spread is that small beside their size should have a
common value subtracted first. The limits are still reported (the sample
mean minus and plus `d`, which involve no variance), but they are then
not an interval to report:
[`summary()`](https://rdrr.io/r/base/summary.html), `glance()`, and
`tidy()` return them without the note, so test
`attr(out, "variance_unusable")` before using them.

**Extreme units.** Data that do vary, but whose spread is so large or so
small that the sample variance cannot be represented in double precision
(differences among the observations beyond about \\10^{154}\\ or below
about \\10^{-162}\\), are an error that says to change the units. That
case is not reported as a rule not yet met, because collecting more data
cannot help.

## References

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_fwd_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_pilot.md);
[`seq_fwd_norm_mean_optimal`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_optimal.md)
for the sample size when the variance is known, the benchmark this
procedure aims for;
[`seq_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md),
which makes no normality assumption;
[`ts_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md)
for the two stage procedure.

Other seq_fwd_norm_mean:
[`seq_fwd_norm_mean_optimal()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_optimal.md),
[`seq_fwd_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Mean score on a 100 point exam, to within 2 points with 95% confidence.
# Scores are simulated here with a standard deviation of 12.
set.seed(113)
m <- seq_fwd_norm_mean_pilot(d = 2, conf_level = 0.95)
scores <- rnorm(m, mean = 70, sd = 12)
out <- seq_fwd_norm_mean(scores, d = 2, conf_level = 0.95)
out
#> Sequential check (fixed width interval (normal mean))
#> 
#>   Current n       : 10
#>   Criterion n     : 126
#>   Target width    : 4
#>   Confidence level: 95%
#> 
#>   lower_limit       : 64.24481
#>   mean              : 66.24481
#>   upper_limit       : 68.24481
#>   optimal_n         : 126
#>   current_n         : 10
#> 
#>   Stopping rule met: FALSE

# Add one student at a time until the stopping rule is met.
while (!attr(out, "is_satisfied")) {
  scores <- c(scores, rnorm(1, mean = 70, sd = 12))
  out <- seq_fwd_norm_mean(scores, d = 2, conf_level = 0.95)
}
out
#> Sequential check (fixed width interval (normal mean))
#> 
#>   Current n       : 146
#>   Criterion n     : 145
#>   Target width    : 4
#>   Confidence level: 95%
#> 
#>   lower_limit       : 67.16453
#>   mean              : 69.16453
#>   upper_limit       : 71.16453
#>   optimal_n         : 145
#>   current_n         : 146
#> 
#>   Stopping rule met: TRUE
```
