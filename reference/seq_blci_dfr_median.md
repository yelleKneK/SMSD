# Bounded Length Confidence Interval for a Median, Distribution Free

Given the data collected so far, checks the stopping rule of the purely
sequential bounded length confidence interval for a population median
when the distribution is not specified (Mukhopadhyay and de Silva,
2009). The interval runs between two order statistics of the sample, and
observations are collected until its length is at most \\2d\\. The
quantity estimated is the population median, not the mean.

## Usage

``` r
seq_blci_dfr_median(x, d, conf_level = 0.95, na.rm = TRUE)
```

## Arguments

- x:

  Numeric vector of the observations collected so far (at least 1; an
  interval is reported once there are enough for the confidence level,
  see Details).

- d:

  Half of the largest interval length allowed, in the units of the data:
  sampling continues until the length of the interval is at most \\2d\\.
  A positive number.

- conf_level:

  The confidence level (default 0.95).

- na.rm:

  Logical; should `NA` values be removed first? Default `TRUE`. With
  `FALSE`, missing values are an error.

## Value

An object of class `c("smsd_check", "data.frame")` with columns `term`
and `value` and these rows:

- `lower_limit`:

  The order statistic \\X\_{(b)}\\.

- `median`:

  The sample median, the point estimate.

- `upper_limit`:

  The order statistic \\X\_{(n-b+1)}\\.

- `interval_length`:

  The upper limit minus the lower limit.

- `current_n`:

  The number of observations used.

The two limits and the length are `NA` when there are too few
observations for the confidence level. Attributes: `is_satisfied`
(`TRUE` when the interval length is at most \\2d\\); `current_n`;
`criterion_n`, which is `NA` because the rule is stated in terms of the
interval length and not a target sample size; `conf_level`; `width` (the
bound \\2d\\); `method`; `min_n` (the smallest sample size at which the
confidence level can be reached); `ranks` (the ranks \\b\\ and \\n - b +
1\\ of the two order statistics); and `notes` (when the rule is met, the
reminder that a sequential procedure covers the median less often than
stated; when there are too few observations, how many are needed;
otherwise empty).

## Details

The function does not collect data and keeps no state. Call it after
each new observation (or batch), starting once
[`seq_blci_dfr_median_pilot`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median_pilot.md)
observations are in hand, and stop sampling the first time
`attr(out, "is_satisfied")` is `TRUE`.

With the \\n\\ observations sorted as \\X\_{(1)} \le \dots \le
X\_{(n)}\\, the interval is \\\[X\_{(b)},\\ X\_{(n-b+1)}\]\\. For a
continuous distribution it covers the median exactly when the number of
observations below the median, a binomial count with \\n\\ trials and
probability 1/2, lies between \\b\\ and \\n - b\\. The rank \\b\\ starts
from the continuity corrected normal approximation to that binomial, \\b
= \max\\1, \lceil (n - z\_{\alpha/2}\sqrt{n} - 1)/2 \rceil\\\\, and is
lowered if necessary until the exact binomial coverage is at least
`conf_level`. (At confidence levels of 0.80 and above the approximation
needed no lowering for any \\n\\ up to 3,000; at lower levels it is
occasionally short by a small amount.) At the 90%, 95%, and 99% levels
this is the largest rank with that property, that is, the narrowest such
interval, at more than 97% of the sample sizes up to 2,000, and one rank
lower at the rest. Sampling stops at the first \\n\\ (at or beyond the
pilot size) at which \\X\_{(n-b+1)} - X\_{(b)} \le 2d\\. Unlike the
fixed width intervals for a mean, this interval is not symmetric about
the point estimate and its length is not fixed in advance; \\2d\\ is an
upper bound on it.

**Minimum sample size.** Below a certain sample size no pair of order
statistics can reach the confidence level: even the interval from the
smallest to the largest observation has coverage \\1 - 2^{1-n}\\ (0.50,
0.75, 0.875, and 0.9375 for \\n\\ = 2 to 5). With fewer than that many
observations (6 at 95%), the rule is reported as not met, the limits and
the interval length are `NA`, and a note says how many observations are
needed. `attr(out, "min_n")` holds that number.

**The confidence level is asymptotic, and the shortfall is not small.**
At any fixed \\n\\ the interval has at least the stated coverage. The
sequential procedure, however, stops when the interval happens to be
short, and short intervals miss the median more often, so its coverage
is below the stated level and approaches it only slowly as `d` shrinks.
In simulation (the pilot of
[`seq_blci_dfr_median_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median_pilot.md)
with its defaults, 40,000 replications in each condition, Monte Carlo
standard error at most 0.002), with `d` chosen so that the sample size
needed (\\n^\*\\ in the Details of
[`seq_blci_dfr_median_pilot`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median_pilot.md))
was 25, 50, 100, 200, 400, and 800, a nominal 0.95 interval covered the
median with probability 0.936, 0.928, 0.926, 0.929, 0.935, and 0.938 for
standard normal data, 0.947, 0.925, 0.918, 0.925, 0.929, and 0.937 for
unit exponential data, and 0.947, 0.937, 0.936, 0.936, 0.938, and 0.937
for standard Cauchy data. The average final sample size was between 0.94
and 1.05 times \\n^\*\\ once \\n^\*\\ was 100 or more. Expect coverage
of roughly 0.92 to 0.94 from a nominal 0.95 interval, unless the pilot
is already a large part of the sample size needed. The shortfall is
larger at a lower level and smaller at a higher one: for the same three
distributions with \\n^\*\\ of 50, 100, 200, and 800, a nominal 0.90
gave between 0.842 and 0.879 and a nominal 0.99 between 0.983 and 0.989.
Because the shortfall is not negligible, a check whose rule is met
carries a note saying so, which
[`print()`](https://rdrr.io/r/base/print.html) shows.

**If the stated level must be met,** ask for more than is needed. With
`d` held at the value for which the 0.95 procedure needs about 100 and
about 400 observations (normal data; 20,000 replications at the first
target and 6,000 at the second, standard error at most 0.003), a nominal
0.96 was not enough (coverage 0.941 and 0.948). A nominal 0.97 gave
0.955 and 0.962, with an average final sample size of 1.20 and 1.17
times those 100 and 400 observations, and a nominal 0.975 gave 0.961 and
0.970, at 1.29 and 1.25 times as many. For exponential data at the first
of the two targets, a nominal 0.97 gave 0.952.

**Ties and constant data.** Ties are permitted. With heavily tied data
the two order statistics can coincide, giving an interval of length
zero; that is a valid (if uninformative about precision) interval for
the median of a discrete population, and the rule is then met. (This
procedure estimates no variance, so the package's rule for a variance
estimate that cannot be used does not arise here.)

## References

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_blci_dfr_median_pilot`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median_pilot.md);
[`seq_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md)
for a fixed width interval for a mean.

Other seq_blci_dfr_median:
[`seq_blci_dfr_median_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Median household income in a district (thousands of dollars), with an
# interval no longer than 10 thousand, 95% confidence. Incomes are right
# skewed; simulated here from a lognormal distribution.
set.seed(113)
m <- seq_blci_dfr_median_pilot(d = 5, conf_level = 0.95)
income <- rlnorm(m, meanlog = log(60), sdlog = 0.5)
out <- seq_blci_dfr_median(income, d = 5, conf_level = 0.95)
out
#> Sequential check (bounded length confidence interval (distribution free))
#> 
#>   Current n       : 10
#>   Criterion n     : NA
#>   Target width    : 10
#>   Confidence level: 95%
#> 
#>   lower_limit       : 31.4192
#>   median            : 51.19733
#>   upper_limit       : 87.24341
#>   interval_length   : 55.82421
#>   current_n         : 10
#> 
#>   Stopping rule met: FALSE

# Add one household at a time until the stopping rule is met. The note
# printed with the result is a reminder that a completed sequential
# procedure covers the median somewhat less often than the stated 95%
# (see Details).
while (!attr(out, "is_satisfied")) {
  income <- c(income, rlnorm(1, meanlog = log(60), sdlog = 0.5))
  out <- seq_blci_dfr_median(income, d = 5, conf_level = 0.95)
}
out
#> Sequential check (bounded length confidence interval (distribution free))
#> 
#>   Current n       : 199
#>   Criterion n     : NA
#>   Target width    : 10
#>   Confidence level: 95%
#> 
#>   lower_limit       : 56.09514
#>   median            : 60.46713
#>   upper_limit       : 65.82721
#>   interval_length   : 9.732074
#>   current_n         : 199
#> 
#>   Stopping rule met: TRUE
#> 
#>   Note: The confidence level holds for a sample size fixed in advance.
#>         A sequential procedure that stops at this check covers the
#>         median less often than stated: in simulation, about 0.84 to
#>         0.88 for a nominal 0.90, 0.92 to 0.94 for 0.95, and 0.98 to
#>         0.99 for 0.99. See the Details of ?seq_blci_dfr_median.

# Too few observations for 95% confidence: not met, with a note.
seq_blci_dfr_median(income[1:4], d = 5, conf_level = 0.95)
#> Sequential check (bounded length confidence interval (distribution free))
#> 
#>   Current n       : 4
#>   Criterion n     : NA
#>   Target width    : 10
#>   Confidence level: 95%
#> 
#>   lower_limit       : NA
#>   median            : 75.69022
#>   upper_limit       : NA
#>   interval_length   : NA
#>   current_n         : 4
#> 
#>   Stopping rule met: FALSE
#> 
#>   Note: With 4 observations no order statistic interval for the median
#>         can reach 95% confidence (at least 6 are needed), so the
#>         stopping rule is reported as not met. Collect more data and
#>         check again.
```
