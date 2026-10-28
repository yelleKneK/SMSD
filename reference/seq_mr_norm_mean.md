# Minimum Risk Point Estimation of a Normal Mean

Given the data collected so far, checks the stopping rule of the purely
sequential minimum risk point estimation (MRPE) procedure for the mean
of a normal population with unknown variance (Mukhopadhyay and de Silva,
2009, Section 6.4, who credit the original formulation to Robbins,
1959). Estimation error and sampling both cost something; the procedure
looks for the sample size at which their expected total, the risk, is
smallest. That sample size depends on the unknown variance, so
observations are collected until the sample size reaches the size the
current variance estimate calls for. The procedure is risk based, so
there is no confidence level.

## Usage

``` r
seq_mr_norm_mean(x, A = NULL, c, price = NULL, epsilon = NULL, na.rm = TRUE)
```

## Arguments

- x:

  Numeric vector of the observations collected so far (at least 2).

- A:

  The loss constant: the price of one unit of squared estimation error,
  in the same currency as `c`. A positive number. Either supply `A`
  directly, or supply `price` and `epsilon`, in which case \\A =
  \mathrm{price}/\mathrm{epsilon}^2\\.

- c:

  The cost of sampling one observation. A positive number.

- price, epsilon:

  An alternative way to state `A` (Chattopadhyay and Kelley, 2016):
  `price` is the amount one is willing to pay for an estimate within
  `epsilon` of the mean, and \\A = \mathrm{price}/\mathrm{epsilon}^2\\.
  Both positive numbers; supply both or neither.

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

  The sample size the current estimate of the standard deviation calls
  for, \\\lceil (A/c)^{1/2} S_n \rceil\\.

- `current_n`:

  The number of observations used.

`risk` and `optimal_n` are `NA` when the data are constant. Attributes:
`is_satisfied` (`TRUE` when `current_n` has reached `optimal_n`),
`current_n`, `criterion_n` (equal to `optimal_n`), `method`, `loss_A`,
`cost_c`, `variance_unusable` (`TRUE` when the data are constant, so
that no variance has been estimated), and `notes`. `conf_level` and
`width` are `NULL`.

## Details

The function does not collect data and keeps no state. Call it after
each new observation (or batch), starting once
[`seq_mr_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean_pilot.md)
observations are in hand, and stop sampling the first time
`attr(out, "is_satisfied")` is `TRUE`.

The loss from estimating the mean \\\mu\\ by the sample mean of \\n\\
observations is \$\$L_n = A(\bar{X}\_n - \mu)^2 + c\\n\$\$ (Mukhopadhyay
and de Silva, 2009, Equation 6.4.1): \\A\\ prices a unit of squared
estimation error and \\c\\ is the cost of one observation. Its
expectation for a fixed \\n\\, the risk, is \\A\sigma^2/n + c\\n\\
(their Equation 6.4.2), which is smallest at \\n^\* =
(A/c)^{1/2}\sigma\\ (their Equation 6.4.3), where it equals \\2
c\\n^\*\\. The sequential rule replaces \\\sigma\\ by the sample
standard deviation \\S_n\\: sampling stops at the first \\n\\ (at or
beyond the pilot size) with \$\$n \ge (A/c)^{1/2} S_n .\$\$ The
`optimal_n` row is the right hand side rounded up, and the `risk` row is
the estimated risk at the current sample size, \\A S_n^2/n + c\\n\\;
both are evaluated at the current data.

**The optimality is asymptotic.** As the cost per observation falls (so
that \\n^\*\\ grows), the expected final sample size approaches \\n^\*\\
and the risk of the procedure approaches the minimum \\2 c\\n^\*\\. When
\\n^\*\\ is small the procedure sometimes stops early on an
underestimate of the standard deviation, and its risk is somewhat above
the minimum. In simulation (normal data, the pilot of
[`seq_mr_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean_pilot.md)
with its defaults, 40,000 replications, Monte Carlo standard error at
most 0.005), with a standard deviation of 1 the ratio of the risk to the
minimum was 1.07, 1.06, 1.04, 1.01, and 1.00 when \\n^\*\\ was 10, 15,
20, 50, and 100 (that is, when \\A/c\\ was 100, 225, 400, 2,500, and
10,000, for which the pilot is 4, 4, 5, 8, and 10). The ratio is not a
function of \\n^\*\\ alone: the pilot depends only on \\A/c\\, so the
same \\n^\*\\ reached with another standard deviation starts from a
different pilot. With a standard deviation of 0.25 (pilots of 7 to 20)
the five ratios were 1.02, 1.03, 1.02, 1.01, and 1.00, and with a
standard deviation of 4 (pilots of 4 or 5) they were 1.07, 1.07, 1.04,
1.00, and 1.00.

**Constant data.** If every observation is equal, or equal up to
rounding error (a range of no more than \\10^{-10}\\ times the largest
absolute value, as with gain scores of 0.3 computed as posttest minus
pretest), the standard deviation estimate is zero, or nothing but
rounding error, and the rule as printed would be met at once. The
function departs from the printed rule deliberately: no variance has
been estimated, so the rule is reported as not met, `risk` and
`optimal_n` are `NA`, `attr(out, "variance_unusable")` is `TRUE`, and a
note says so (see
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

Robbins, H. (1959). Sequential estimation of the mean of a normal
population. In U. Grenander (Ed.), *Probability and statistics: The
Harald Cramer volume* (pp. 235–245). Almqvist & Wiksell.

## See also

[`seq_mr_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean_pilot.md)
for the pilot. With no normality assumption:
[`seq_mr_dfr_mean_arb`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md).
To bound the risk of the same kind of loss instead of trading it against
cost:
[`seq_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md).

Other seq_mr_norm_mean:
[`seq_mr_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean_pilot.md)

## Author

Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>, Ken Kelley
<kkelley@nd.edu>

## Examples

``` r
# A school district wants its mean test score. An error of one point
# (squared) is valued at 2,500 dollars and testing one more student
# costs 4 dollars. Scores are simulated here with a standard deviation
# of 12, so the best sample size is sqrt(2500 / 4) * 12 = 300.
set.seed(113)
m <- seq_mr_norm_mean_pilot(A = 2500, c = 4)
scores <- rnorm(m, mean = 70, sd = 12)
out <- seq_mr_norm_mean(scores, A = 2500, c = 4)
out
#> Sequential check (minimum risk point estimation (normal mean))
#> 
#>   Current n       : 5
#>   Criterion n     : 316
#> 
#>   mean              : 70.97121
#>   risk              : 79636.06
#>   optimal_n         : 316
#>   current_n         : 5
#> 
#>   Stopping rule met: FALSE

# Add one student at a time until the stopping rule is met.
while (!attr(out, "is_satisfied")) {
  scores <- c(scores, rnorm(1, mean = 70, sd = 12))
  out <- seq_mr_norm_mean(scores, A = 2500, c = 4)
}
out
#> Sequential check (minimum risk point estimation (normal mean))
#> 
#>   Current n       : 299
#>   Criterion n     : 299
#> 
#>   mean              : 69.45986
#>   risk              : 2386.488
#>   optimal_n         : 299
#>   current_n         : 299
#> 
#>   Stopping rule met: TRUE
```
