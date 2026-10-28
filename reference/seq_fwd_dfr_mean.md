# Fixed Width Confidence Interval for a Mean, Distribution Free

Given the data collected so far, checks the stopping rule of the purely
sequential fixed width confidence interval for a population mean when
the distribution is not specified (Mukhopadhyay and de Silva, 2009). The
aim is an interval of half width `d` (full width \\2d\\) around the
sample mean with approximately the stated confidence level. The sample
size needed depends on the unknown variance, so observations are
collected until the sample size reaches the size the current variance
estimate calls for.

## Usage

``` r
seq_fwd_dfr_mean(x, d, conf_level = 0.95, na.rm = TRUE)
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

  The sample size the current data call for, \\\lceil z\_{\alpha/2}^2
  (S_n^2 + 1/n) / d^2 \rceil\\ (`NA` when the data are constant).

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
[`seq_fwd_dfr_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean_pilot.md)
observations are in hand, and stop sampling the first time
`attr(out, "is_satisfied")` is `TRUE`.

Sampling stops at the first \\n\\ (at or beyond the pilot size) with
\$\$n \ge \frac{z\_{\alpha/2}^2}{d^2}\left(S_n^2 +
\frac{1}{n}\right),\$\$ where \\S_n^2\\ is the sample variance. The
\\1/n\\ term keeps the rule from being met on a variance estimate that
happens to be very small early on; it vanishes as the sample grows. The
`optimal_n` row is the right hand side rounded up, evaluated at the
current data. The reported interval is always \\\bar{X}\_n \pm d\\; it
has the stated confidence level, approximately, only once the rule is
met. The rule has the form of the stopping rule of Kelley, Darku, and
Chattopadhyay (2018, Equation 22), with the sample variance as the
estimate of the asymptotic variance and their full width \\\omega\\
equal to \\2d\\.

**The confidence level is asymptotic.** It rests on the central limit
theorem and on the consistency of \\S_n^2\\, and is attained in the
limit as `d` shrinks. For wide targets, and for skewed data even at
fairly narrow ones, coverage is below the stated level. In simulation
(standard deviation 2, nominal 0.95, the pilot of
[`seq_fwd_dfr_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean_pilot.md)
with its defaults, 40,000 replications, Monte Carlo standard error at
most 0.002), coverage for normal data was 0.935 at `d = 1` (where the
required sample size \\n^\* = z\_{\alpha/2}^2\sigma^2/d^2\\ is about
15), 0.942 at `d = 0.5` (\\n^\*\\ about 61), and 0.948 at `d = 0.25`
(\\n^\*\\ about 246). For exponential data with the same standard
deviation it was 0.960 at `d = 1` (where the pilot of 10 is already most
of the 15 required) and then far lower: 0.879 at \\n^\*\\ of 30, and
lowest, about 0.855, when \\n^\*\\ was between 40 and 61 (`d = 0.5`). It
recovered slowly: 0.886 at \\n^\*\\ of 100, 0.934 at `d = 0.25`
(\\n^\*\\ about 246), 0.946 at `d = 0.125` (\\n^\*\\ about 983; 400,000
replications), and still 0.946 at `d = 0.0625` (\\n^\*\\ about 3,934;
10,000 replications, standard error 0.002). Skewness slows the approach
to the stated level considerably.

**The units of the data.** The \\1/n\\ term is added to a variance, so
it is in the squared units of the data, and the procedure depends on the
units in which the data are recorded. It adds up to about \\1/\sigma^2\\
observations to the sample size that would otherwise be needed: a
negligible number when the standard deviation is 1 or more in the units
of the data, and a substantial one when it is much smaller. In
simulation with normal data, \\\sigma = 0.25\\, and `d = 0.1` (\\n^\*\\
about 24), the average final sample size was 35 and coverage was 0.98.
The interval is then conservative, not wrong.

**Constant data.** If every observation is equal, or equal up to
rounding error (a range of no more than \\10^{-10}\\ times the largest
absolute value, as with gain scores of 0.3 computed as posttest minus
pretest), the sample variance is zero, or nothing but rounding error,
and under the rule as printed the \\1/n\\ term alone would eventually
stop the sampling. The function departs from the printed rule
deliberately: no variance has been estimated, so the rule is reported as
not met, `optimal_n` is `NA`, `attr(out, "variance_unusable")` is
`TRUE`, and a note says so (see
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

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_fwd_dfr_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean_pilot.md);
[`seq_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md)
for the normal theory rule without the \\1/n\\ term;
[`ts_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md)
for the two stage procedure;
[`seq_aipe_sm`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm.md)
for a mean standardized by its standard deviation.

Other seq_fwd_dfr_mean:
[`seq_fwd_dfr_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Mean weekly hours of remote work in a firm, to within half an hour
# with 90% confidence. Hours are right skewed; simulated here.
set.seed(113)
m <- seq_fwd_dfr_mean_pilot(d = 0.5, conf_level = 0.90)
hours <- rgamma(m, shape = 4, scale = 2)
out <- seq_fwd_dfr_mean(hours, d = 0.5, conf_level = 0.90)
out
#> Sequential check (fixed width interval (distribution free mean))
#> 
#>   Current n       : 10
#>   Criterion n     : 112
#>   Target width    : 1
#>   Confidence level: 90%
#> 
#>   lower_limit       : 6.478901
#>   mean              : 6.978901
#>   upper_limit       : 7.478901
#>   optimal_n         : 112
#>   current_n         : 10
#> 
#>   Stopping rule met: FALSE

# Add one employee at a time until the stopping rule is met.
while (!attr(out, "is_satisfied")) {
  hours <- c(hours, rgamma(1, shape = 4, scale = 2))
  out <- seq_fwd_dfr_mean(hours, d = 0.5, conf_level = 0.90)
}
out
#> Sequential check (fixed width interval (distribution free mean))
#> 
#>   Current n       : 150
#>   Criterion n     : 150
#>   Target width    : 1
#>   Confidence level: 90%
#> 
#>   lower_limit       : 7.177634
#>   mean              : 7.677634
#>   upper_limit       : 8.177634
#>   optimal_n         : 150
#>   current_n         : 150
#> 
#>   Stopping rule met: TRUE
```
