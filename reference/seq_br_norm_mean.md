# Bounded Risk Point Estimation of a Normal Mean

Given the data collected so far, checks the stopping rule of the purely
sequential bounded risk point estimator of the mean of a normal
population with unknown variance (Mukhopadhyay and de Silva, 2009). The
aim is a sample mean whose expected loss (its risk) is no larger than a
bound `w` chosen in advance. No fixed sample size can promise that,
because the risk depends on the unknown variance, so observations are
collected until the sample size reaches the size the current variance
estimate calls for. The procedure is risk based, so there is no
confidence level.

## Usage

``` r
seq_br_norm_mean(x, A, k, w, na.rm = TRUE)
```

## Arguments

- x:

  Numeric vector of the observations collected so far (at least 2).

- A:

  The loss constant: the loss from estimating the mean \\\mu\\ by the
  sample mean is \\A\\\|\bar{X}\_n - \mu\|^k\\. A positive number; it
  sets the units in which the risk bound `w` is expressed (use `A = 1`
  to bound the mean absolute or mean squared error itself).

- k:

  The power in the loss: 1 for absolute error, 2 for squared error.

- w:

  The risk bound: the largest expected loss that is acceptable. A
  positive number, in the units of `A` times the data's units to the
  power `k`.

- na.rm:

  Logical; should `NA` values be removed first? Default `TRUE`. With
  `FALSE`, missing values are an error.

## Value

An object of class `c("smsd_check", "data.frame")` with columns `term`
and `value` and these rows:

- `mean`:

  The sample mean, the point estimate.

- `sd`:

  The sample standard deviation.

- `risk`:

  The estimated risk at the current sample size, \\B\\(S_n^2/n)^{k/2}\\.

- `optimal_n`:

  The sample size the current variance estimate calls for, \\\lceil
  (B/w)^{2/k} S_n^2 \rceil\\.

- `current_n`:

  The number of observations used.

`risk` and `optimal_n` are `NA` when the data are constant. Attributes:
`is_satisfied` (`TRUE` when `current_n` has reached `optimal_n`, which
is when the estimated risk is at most `w`; the rounding up of
`optimal_n` allows for rounding error of one part in \\10^{10}\\, so
when the estimated risk and `w` agree to about ten significant digits
the two statements can differ, and `is_satisfied` is the authority),
`current_n`, `criterion_n` (equal to `optimal_n`), `method`, `loss_A`,
`k`, `risk_bound`, `variance_unusable` (`TRUE` when the data are
constant, so that no variance has been estimated), and `notes`.
`conf_level` and `width` are `NULL`.

## Details

The function does not collect data and keeps no state. Call it after
each new observation (or batch), starting once
[`seq_br_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean_pilot.md)
observations are in hand, and stop sampling the first time
`attr(out, "is_satisfied")` is `TRUE`.

The loss from estimating the mean \\\mu\\ by the sample mean is
\\A\\\|\bar{X}\_n - \mu\|^k\\ with \\k\\ equal to 1 (absolute error) or
2 (squared error). For normal data and a fixed \\n\\ its expectation is
\$\$B\\(\sigma^2/n)^{k/2}, \qquad B =
A\\2^{k/2}\\\Gamma\\(k+1)/2\\/\sqrt{\pi},\$\$ that is, \\B =
A\sqrt{2/\pi}\\ for \\k = 1\\ and \\B = A\\ for \\k = 2\\. The risk is
at most `w` when \\n \ge n^\* = (B/w)^{2/k}\sigma^2\\. The sequential
rule replaces \\\sigma^2\\ by the sample variance \\S_n^2\\: sampling
stops at the first \\n\\ (at or beyond the pilot size) with \$\$n \ge
(B/w)^{2/k} S_n^2 .\$\$ The `optimal_n` row is the right hand side
rounded up, and the `risk` row is the estimated risk at the current
sample size, \\B\\(S_n^2/n)^{k/2}\\; both are evaluated at the current
data.

**The risk bound is asymptotic.** The risk of the sequential procedure
approaches `w` as `w` shrinks (so that the required sample size grows),
and it approaches it from above: for squared error and normal data the
risk is \\A\sigma^2 E(1/N)\\, where \\N\\ is the final sample size,
which exceeds `w` whenever the expected final sample size is at or below
\\n^\*\\. With a loose bound the excess is not small. In simulation
(normal data with standard deviation 3, `A = 2`, the pilot of
[`seq_br_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean_pilot.md)
with its defaults, which was 10 at every bound used, 40,000
replications), the ratio of the risk to `w` for squared error (`k = 2`)
was 1.11 when the bound called for 15 observations, highest, about 1.22,
when it called for 25 to 30, and then 1.20, 1.05, and 1.01 when it
called for 36, 100, and 400 (`w` of 0.5, 0.18, and 0.045). These were
computed from the identity above, with a Monte Carlo standard error of
at most 0.003. For absolute error (`k = 1`) the ratio was 1.04 at 15
observations, highest, between 1.08 and 1.09, at 25 to 30, and 1.07,
1.02, and 1.00 at 36, 100, and 400 (standard error 0.004). A larger
`m_0` in the pilot reduces the excess for loose bounds.
[`ts_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md)
is the two stage alternative, whose risk stays at or under the bound for
normal data at the price of a larger expected sample size.

**Constant data.** If every observation is equal, or equal up to
rounding error (a range of no more than \\10^{-10}\\ times the largest
absolute value, as with gain scores of 0.3 computed as posttest minus
pretest), the sample variance is zero, or nothing but rounding error,
and the rule as printed would be met at once. The function departs from
the printed rule deliberately: no variance has been estimated, so the
rule is reported as not met, `risk` and `optimal_n` are `NA`,
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

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_br_norm_mean_pilot`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean_pilot.md);
[`ts_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md)
for the two stage procedure;
[`seq_mr_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean.md),
which trades the same loss against the cost of sampling instead of
bounding it.

Other seq_br_norm_mean:
[`seq_br_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Mean rating on a 7 point satisfaction scale, with a mean absolute
# error of at most 0.25 points. Ratings are simulated here as normal
# with a standard deviation of 1.2.
set.seed(113)
m <- seq_br_norm_mean_pilot(A = 1, k = 1, w = 0.25)
ratings <- rnorm(m, mean = 4.8, sd = 1.2)
out <- seq_br_norm_mean(ratings, A = 1, k = 1, w = 0.25)
out
#> Sequential check (bounded risk estimation (normal mean))
#> 
#>   Current n       : 10
#>   Criterion n     : 14
#> 
#>   mean              : 4.424481
#>   sd                : 1.145169
#>   risk              : 0.2889413
#>   optimal_n         : 14
#>   current_n         : 10
#> 
#>   Stopping rule met: FALSE

# Add one respondent at a time until the stopping rule is met.
while (!attr(out, "is_satisfied")) {
  ratings <- c(ratings, rnorm(1, mean = 4.8, sd = 1.2))
  out <- seq_br_norm_mean(ratings, A = 1, k = 1, w = 0.25)
}
out
#> Sequential check (bounded risk estimation (normal mean))
#> 
#>   Current n       : 16
#>   Criterion n     : 15
#> 
#>   mean              : 4.654325
#>   sd                : 1.193281
#>   risk              : 0.2380251
#>   optimal_n         : 15
#>   current_n         : 16
#> 
#>   Stopping rule met: TRUE
```
