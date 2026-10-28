# Minimum Risk Point Estimation of a Mean, Distribution Free

Given the data collected so far, checks the stopping rule of the purely
sequential minimum risk point estimation (MRPE) procedure for a
population mean when the distribution is not specified and the variance
is unknown and unrestricted (Mukhopadhyay and de Silva, 2009).
Estimation error and sampling both cost something; the procedure looks
for the sample size at which their expected total, the risk, is
smallest. That sample size depends on the unknown variance, so
observations are collected until the sample size reaches the size the
current estimate calls for. The procedure is risk based, so there is no
confidence level.

## Usage

``` r
seq_mr_dfr_mean_arb(x, A, c, gamma = 1, na.rm = TRUE)
```

## Arguments

- x:

  Numeric vector of the observations collected so far (at least 2).

- A:

  The loss constant: the price of one unit of squared estimation error,
  in the same currency as `c`. A positive number.

- c:

  The cost of sampling one observation. A positive number.

- gamma:

  The positive constant of the correction term \\n^{-\gamma}\\ in the
  stopping rule; default 1. See Details for the tradeoff. Use the same
  value in
  [`seq_mr_dfr_mean_arb_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb_pilot.md).

- na.rm:

  Logical; should `NA` values be removed first? Default `TRUE`. With
  `FALSE`, missing values are an error.

## Value

An object of class `c("smsd_check", "data.frame")` with columns `term`
and `value` and these rows:

- `mean`:

  The sample mean, the point estimate.

- `risk`:

  The estimated risk at the current sample size, \\A S_n^2/n + c\\n\\.

- `optimal_n`:

  The sample size the current data call for, \\\lceil (A/c)^{1/2}(S_n +
  n^{-\gamma}) \rceil\\.

- `current_n`:

  The number of observations used.

`risk` and `optimal_n` are `NA` when the data are constant. Attributes:
`is_satisfied` (`TRUE` when `current_n` has reached `optimal_n`),
`current_n`, `criterion_n` (equal to `optimal_n`), `method`, `loss_A`,
`cost_c`, `gamma`, `variance_unusable` (`TRUE` when the data are
constant, so that no variance has been estimated), and `notes`.
`conf_level` and `width` are `NULL`.

## Details

The function does not collect data and keeps no state. Call it after
each new observation (or batch), starting once
[`seq_mr_dfr_mean_arb_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb_pilot.md)
observations are in hand, and stop sampling the first time
`attr(out, "is_satisfied")` is `TRUE`.

The loss from estimating the mean \\\mu\\ by the sample mean of \\n\\
observations is \\A(\bar{X}\_n - \mu)^2 + c\\n\\: \\A\\ prices a unit of
squared estimation error and \\c\\ is the cost of one observation. Its
expectation for a fixed \\n\\, the risk, is \\A\sigma^2/n + c\\n\\ for
any distribution with variance \\\sigma^2\\; it is smallest at \\n^\* =
(A/c)^{1/2}\sigma\\, where it equals \\2 c\\n^\*\\. Sampling stops at
the first \\n\\ (at or beyond the pilot size) with \$\$n \ge
(A/c)^{1/2}\\(S_n + n^{-\gamma}),\$\$ where \\S_n\\ is the sample
standard deviation. The correction term \\n^{-\gamma}\\ keeps the rule
from being met on a standard deviation estimate that happens to be very
small early on; it vanishes as the sample grows. The `optimal_n` row is
the right hand side rounded up, and the `risk` row is the estimated risk
at the current sample size, \\A S_n^2/n + c\\n\\; both are evaluated at
the current data.

**The choice of `gamma`, and the units of the data.** The correction
term is added to a standard deviation, so it is in the units of the
data, and the procedure depends on the units in which the data are
recorded. It inflates the final sample size \\n\\ by a factor of roughly
\\1 + n^{-\gamma}/\sigma\\: negligible when the standard deviation is 1
or more in the units of the data, and not negligible when it is much
smaller unless `gamma` is large. A smaller `gamma` gives more protection
against stopping too early at the price of more observations.
Chattopadhyay and Kelley (2016, p. 636) use 0.49 for a correction of the
same kind in their rule for the coefficient of variation, and note that
a value above 1/2 falls outside part of their theory. The default here
is 1 because it behaved well across scales in simulation, where 0.49
oversampled at small scales (see below).

**The optimality is asymptotic.** As the cost per observation falls (so
that \\n^\*\\ grows), the expected final sample size approaches \\n^\*\\
and the risk approaches the minimum \\2 c\\n^\*\\. In simulation
(exponential data, `A = 100`, `c = 0.01`, so that \\n^\*\\ is 100 times
the standard deviation, the pilot of
[`seq_mr_dfr_mean_arb_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb_pilot.md)
with matching `gamma`, 40,000 replications, Monte Carlo standard errors
at most 0.004), with `gamma = 1` the expected final sample size was
1.10, 0.99, and 0.99 times \\n^\*\\ at standard deviations of 0.25, 1,
and 3, and the risk was 1.02, 1.02, and 1.01 times the minimum. With
`gamma = 0.49` the sample size ratios were 1.64, 1.08, and 1.01 and the
risk ratios 1.13, 1.01, and 1.01.

**Constant data.** If every observation is equal, or equal up to
rounding error (a range of no more than \\10^{-10}\\ times the largest
absolute value, as with gain scores of 0.3 computed as posttest minus
pretest), the sample variance is zero, or nothing but rounding error,
and under the rule as printed the correction term alone would eventually
stop the sampling. The function departs from the printed rule
deliberately: no variance has been estimated, so the rule is reported as
not met, `risk` and `optimal_n` are `NA`,
`attr(out, "variance_unusable")` is `TRUE`, and a note says so (see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md));
data whose real spread is that small beside their size should have a
common value subtracted first.

**Extreme units.** Data that do vary, but whose spread is so large or so
small that the sample variance cannot be represented in double precision
(differences among the observations beyond about \\10^{154}\\ or below
about \\10^{-162}\\), are an error that says to change the units. That
case is not reported as a rule not yet met, because collecting more data
cannot help.

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_mr_dfr_mean_arb_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb_pilot.md);
[`seq_mr_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean.md)
for the normal theory rule without the correction term;
[`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
and
[`seq_mr_smd`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)
for minimum risk estimation of effect sizes.

Other seq_mr_dfr_mean_arb:
[`seq_mr_dfr_mean_arb_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb_pilot.md)

## Author

Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>, Ken Kelley
<kkelley@nd.edu>

## Examples

``` r
# Mean order value for an online retailer, in dollars. An error of one
# dollar (squared) is worth 400 dollars and auditing one more order
# costs 0.25 dollars. Order values are right skewed; they are simulated
# here with a standard deviation of 10, so the best sample size is 400
# (the square root of 400 / 0.25, times 10).
set.seed(113)
m <- seq_mr_dfr_mean_arb_pilot(A = 400, c = 0.25)
orders <- rgamma(m, shape = 4, scale = 5)
out <- seq_mr_dfr_mean_arb(orders, A = 400, c = 0.25)
out
#> Sequential check (minimum risk point estimation (distribution free mean, arbitrary variance))
#> 
#>   Current n       : 7
#>   Criterion n     : 345
#> 
#>   mean              : 19.85849
#>   risk              : 4092.808
#>   optimal_n         : 345
#>   current_n         : 7
#> 
#>   Stopping rule met: FALSE

# Add one order at a time until the stopping rule is met.
while (!attr(out, "is_satisfied")) {
  orders <- c(orders, rgamma(1, shape = 4, scale = 5))
  out <- seq_mr_dfr_mean_arb(orders, A = 400, c = 0.25)
}
out
#> Sequential check (minimum risk point estimation (distribution free mean, arbitrary variance))
#> 
#>   Current n       : 401
#>   Criterion n     : 401
#> 
#>   mean              : 20.27692
#>   risk              : 200.2885
#>   optimal_n         : 401
#>   current_n         : 401
#> 
#>   Stopping rule met: TRUE
```
