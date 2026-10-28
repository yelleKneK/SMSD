# Minimum Risk Point Estimation of the Location of a Negative Exponential Distribution

Given the data collected so far, checks the stopping rule of the purely
sequential minimum risk point estimation (MRPE) procedure for the
location parameter of a negative exponential distribution (Mukhopadhyay
and de Silva, 2009). The negative exponential (two parameter
exponential) distribution has density \\\sigma^{-1}\exp\\-(x -
\theta)/\sigma\\\\ for \\x \> \theta\\: \\\theta\\ is the location, the
smallest value the variable can take (a minimum response time, a
guarantee period), and \\\sigma\\ is the scale. Despite the function's
name, the quantity estimated is the location \\\theta\\, not the mean
(which is \\\theta + \sigma\\). The procedure is risk based, so there is
no confidence level.

## Usage

``` r
seq_mr_nexp_mean(x, A, c, k, t, na.rm = TRUE)
```

## Arguments

- x:

  Numeric vector of the observations collected so far (at least 2).

- A:

  The loss constant: the price of one unit of estimation error raised to
  the power `k`, in the same currency as `c`. A positive number.

- c:

  The cost constant: sampling \\n\\ observations costs \\c\\n^t\\. A
  positive number.

- k:

  The power applied to the estimation error in the loss: 1 for absolute
  error, 2 for squared error. Any positive number.

- t:

  The power applied to the sample size in the cost: 1 when each
  observation costs the same, more than 1 when later observations cost
  more. Any positive number.

- na.rm:

  Logical; should `NA` values be removed first? Default `TRUE`. With
  `FALSE`, missing values are an error.

## Value

An object of class `c("smsd_check", "data.frame")` with columns `term`
and `value` and these rows:

- `location`:

  The sample minimum, the point estimate of the location \\\theta\\.

- `risk`:

  The estimated risk at the current sample size,
  \\A\\\Gamma(k+1)\\(U_n/n)^k + c\\n^t\\.

- `optimal_n`:

  The sample size the current scale estimate calls for, \\\lceil
  \\A\\k\\\Gamma(k+1)\\U_n^k/(c\\t)\\^{1/(t+k)} \rceil\\.

- `current_n`:

  The number of observations used.

`risk` and `optimal_n` are `NA` when the data are constant. Attributes:
`is_satisfied` (`TRUE` when `current_n` has reached `optimal_n`),
`current_n`, `criterion_n` (equal to `optimal_n`), `method`, `loss_A`,
`cost_c`, `k`, `t`, `variance_unusable` (`TRUE` when the data are
constant, so that no scale has been estimated; the attribute has the
same name in every procedure of the package), and `notes`. `conf_level`
and `width` are `NULL`.

## Details

The function does not collect data and keeps no state. Call it after
each new observation (or batch), starting once
[`seq_mr_nexp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean_pilot.md)
observations are in hand, and stop sampling the first time
`attr(out, "is_satisfied")` is `TRUE`.

The location is estimated by the sample minimum \\X\_{(1)}\\, and the
loss from \\n\\ observations is \$\$A\\(X\_{(1)} - \theta)^k + c\\n^t
.\$\$ Here \\A\\ prices a unit of estimation error raised to the power
\\k\\, and \\c\\n^t\\ is the cost of sampling. Because \\X\_{(1)} -
\theta\\ is exponential with scale \\\sigma/n\\, the risk for a fixed
\\n\\ is \$\$A\\\Gamma(k+1)\\(\sigma/n)^k + c\\n^t ,\$\$ which is
smallest at \\n^\* = \\A\\k\\\Gamma(k+1)\\\sigma^k/(c\\t)\\^{1/(t+k)}\\.
The sequential rule replaces \\\sigma\\ by \\U_n = \sum (X_i -
X\_{(1)})/(n - 1)\\: sampling stops at the first \\n\\ (at or beyond the
pilot size) with \$\$n \ge \\A\\k\\\Gamma(k+1)\\U_n^k/(c\\t)\\^{1/(t+k)}
.\$\$ The `optimal_n` row is the right hand side rounded up, and the
`risk` row is the estimated risk at the current sample size,
\\A\\\Gamma(k+1)\\(U_n/n)^k + c\\n^t\\; both are evaluated at the
current data.

**The optimality is asymptotic.** As the cost constant falls (so that
\\n^\*\\ grows), the expected final sample size approaches \\n^\*\\ and
the risk of the procedure approaches the minimum. In simulation
(`A = 500`, `c = 0.01`, each combination of `k` and `t` equal to 1 or 2,
scales 0.25, 1, and 3, the pilot of
[`seq_mr_nexp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean_pilot.md)
with its defaults, 40,000 replications, Monte Carlo standard error at
most 0.007), the ratio of the risk to the minimum was between 1.00 and
1.03 in ten of the twelve conditions, in which \\n^\*\\ ran from 18 to
387. The two exceptions both had `k = 2` and `t = 2`: the ratio was 1.04
at scale 1, where \\n^\*\\ is 18, and 1.09 at scale 0.25, where \\n^\*\\
is 9 (400,000 replications each, standard error 0.002).

**Data.** Any finite values are accepted, because the location may be
negative. The model is right skewed with a sharp lower endpoint; the
function cannot check that.

**Constant data.** If every observation is equal, or equal up to
rounding error (a range of no more than \\10^{-10}\\ times the largest
absolute value), \\U_n\\ is zero, or nothing but rounding error, and the
rule as printed would be met at once. The function departs from the
printed rule deliberately: no scale has been estimated, so the rule is
reported as not met, `risk` and `optimal_n` are `NA`,
`attr(out, "variance_unusable")` is `TRUE`, and a note says so (see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md));
data whose real spread is that small beside their size should have a
common value subtracted first.

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

[`seq_mr_nexp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean_pilot.md);
[`seq_fwd_nexp_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean.md)
for a fixed width confidence interval for the same location.

Other seq_mr_nexp_mean:
[`seq_mr_nexp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean_pilot.md)

## Author

Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>, Ken Kelley
<kkelley@nd.edu>

## Examples

``` r
# Shortest possible delivery time for an online retailer, in days. An
# error of one day is valued at 500 dollars and tracking one more order
# costs 0.25 dollars. Times are simulated here as a minimum of 2 days
# plus an exponential delay averaging 1.5 days, so the best sample size
# is sqrt(500 * 1.5 / 0.25), about 55.
set.seed(113)
m <- seq_mr_nexp_mean_pilot(A = 500, c = 0.25, k = 1, t = 1)
days <- 2 + rexp(m, rate = 1 / 1.5)
out <- seq_mr_nexp_mean(days, A = 500, c = 0.25, k = 1, t = 1)
out
#> Sequential check (minimum risk point estimation (negative exponential location))
#> 
#>   Current n       : 7
#>   Criterion n     : 50
#> 
#>   location          : 2.15913
#>   risk              : 90.03639
#>   optimal_n         : 50
#>   current_n         : 7
#> 
#>   Stopping rule met: FALSE

# Add one order at a time until the stopping rule is met.
while (!attr(out, "is_satisfied")) {
  days <- c(days, 2 + rexp(1, rate = 1 / 1.5))
  out <- seq_mr_nexp_mean(days, A = 500, c = 0.25, k = 1, t = 1)
}
out
#> Sequential check (minimum risk point estimation (negative exponential location))
#> 
#>   Current n       : 53
#>   Criterion n     : 53
#> 
#>   location          : 2.001436
#>   risk              : 26.17606
#>   optimal_n         : 53
#>   current_n         : 53
#> 
#>   Stopping rule met: TRUE
```
