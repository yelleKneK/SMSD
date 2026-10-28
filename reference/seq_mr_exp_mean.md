# Minimum Risk Point Estimation of an Exponential Mean

Given the data collected so far, checks the stopping rule of the purely
sequential minimum risk point estimation (MRPE) procedure for the mean
of an exponential population, such as a mean waiting time or time to
failure (Mukhopadhyay and de Silva, 2009). Estimation error and sampling
both cost something; the procedure looks for the sample size at which
their expected total, the risk, is smallest. That sample size depends on
the unknown mean, so observations are collected until the sample size
reaches the size the current estimate calls for. The procedure is risk
based, so there is no confidence level.

## Usage

``` r
seq_mr_exp_mean(x, A, c, na.rm = TRUE)
```

## Arguments

- x:

  Numeric vector of the observations collected so far (at least 1). All
  values must be positive.

- A:

  The loss constant: the price of one unit of squared estimation error,
  in the same currency as `c`. A positive number.

- c:

  The cost of sampling one observation. A positive number.

- na.rm:

  Logical; should `NA` values be removed first? Default `TRUE`. With
  `FALSE`, missing values are an error.

## Value

An object of class `c("smsd_check", "data.frame")` with columns `term`
and `value` and these rows:

- `mean`:

  The sample mean, the point estimate.

- `risk`:

  The estimated risk at the current sample size, \\A\bar{X}\_n^2/n +
  c\\n\\.

- `optimal_n`:

  The sample size the current estimate of the mean calls for, \\\lceil
  (A/c)^{1/2}\\\bar{X}\_n \rceil\\.

- `current_n`:

  The number of observations used.

Attributes: `is_satisfied` (`TRUE` once `current_n` is at least
`optimal_n`), `current_n`, `criterion_n` (equal to `optimal_n`),
`method`, `loss_A`, `cost_c`, and `notes`. `conf_level` and `width` are
`NULL`.

## Details

The function does not collect data and keeps no state. Call it after
each new observation (or batch), starting once
[`seq_mr_exp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean_pilot.md)
observations are in hand, and stop sampling the first time
`attr(out, "is_satisfied")` is `TRUE`.

The loss from estimating the mean \\\lambda\\ by the sample mean of
\\n\\ observations is \\A(\bar{X}\_n - \lambda)^2 + c\\n\\: \\A\\ prices
a unit of squared estimation error and \\c\\ is the cost of one
observation. For exponential data the variance of the sample mean is
\\\lambda^2/n\\, so the risk for a fixed \\n\\ is \\A\lambda^2/n +
c\\n\\, which is smallest at \\n^\* = (A/c)^{1/2}\lambda\\, where it
equals \\2 c\\n^\*\\. The sequential rule replaces \\\lambda\\ by the
sample mean: sampling stops at the first \\n\\ (at or beyond the pilot
size) with \$\$n \ge (A/c)^{1/2}\\\bar{X}\_n .\$\$ The `optimal_n` row
is the right hand side rounded up, and the `risk` row is the estimated
risk at the current sample size, \\A\bar{X}\_n^2/n + c\\n\\; both are
evaluated at the current data.

**The optimality is asymptotic.** As the cost per observation falls (so
that \\n^\*\\ grows), the expected final sample size approaches \\n^\*\\
and the risk of the procedure approaches the minimum \\2 c\\n^\*\\. When
\\n^\*\\ is small the risk is noticeably above the minimum, more so than
for a normal mean, because the sample mean of skewed data is itself the
estimate of spread and a low early mean stops the sampling. In
simulation (exponential data, the pilot of
[`seq_mr_exp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean_pilot.md)
with its defaults, 40,000 replications, Monte Carlo standard error at
most 0.006), with a mean of 1 the ratio of the risk to the minimum was
1.08, 1.16, 1.15, 1.05, and 1.02 when \\n^\*\\ was 10, 15, 20, 50, and
100 (that is, when \\A/c\\ was 100, 225, 400, 2,500, and 10,000, for
which the pilot is 4, 4, 5, 8, and 10). The ratio is not a function of
\\n^\*\\ alone: the pilot depends only on \\A/c\\, so the same \\n^\*\\
reached with another mean starts from a different pilot. With a mean of
0.25 (pilots of 7 to 20) the five ratios were 0.95, 1.04, 1.07, 1.04,
and 1.01, and with a mean of 4 (pilots of 4 or 5) they were 1.08, 1.16,
1.16, 1.06, and 1.02. (A ratio below 1 is possible because the final
sample size depends on the data.)

**Data.** The exponential model is for positive data. Zero or negative
values are an error: they are not possible under the model, and a
negative sample mean would give a negative sample size.

## References

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_mr_exp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean_pilot.md);
[`ts_br_exp_mean`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md)
for two stage bounded risk estimation of the same mean;
[`seq_mr_dfr_mean_arb`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md),
which makes no distributional assumption.

Other seq_mr_exp_mean:
[`seq_mr_exp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean_pilot.md)

## Author

Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>, Ken Kelley
<kkelley@nd.edu>

## Examples

``` r
# Mean time to resolve a support ticket, in minutes. An error of one
# minute (squared) is valued at 100 dollars and timing one more ticket
# costs 0.04 dollars. Times are simulated here with a mean of 6
# minutes, so the best sample size is sqrt(100 / 0.04) * 6 = 300.
set.seed(113)
m <- seq_mr_exp_mean_pilot(A = 100, c = 0.04)
minutes <- rexp(m, rate = 1 / 6)
out <- seq_mr_exp_mean(minutes, A = 100, c = 0.04)
out
#> Sequential check (minimum risk point estimation (exponential mean))
#> 
#>   Current n       : 8
#>   Criterion n     : 233
#> 
#>   mean              : 4.642232
#>   risk              : 269.699
#>   optimal_n         : 233
#>   current_n         : 8
#> 
#>   Stopping rule met: FALSE

# Add one ticket at a time until the stopping rule is met.
while (!attr(out, "is_satisfied")) {
  minutes <- c(minutes, rexp(1, rate = 1 / 6))
  out <- seq_mr_exp_mean(minutes, A = 100, c = 0.04)
}
out
#> Sequential check (minimum risk point estimation (exponential mean))
#> 
#>   Current n       : 295
#>   Criterion n     : 295
#> 
#>   mean              : 5.897602
#>   risk              : 23.59041
#>   optimal_n         : 295
#>   current_n         : 295
#> 
#>   Stopping rule met: TRUE
```
