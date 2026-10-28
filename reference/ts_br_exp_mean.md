# Two Stage Bounded Risk Point Estimation of an Exponential Mean

Evaluates the two stage bounded risk point estimation procedure for the
mean \\\theta\\ of an exponential distribution (a mean time to failure,
a mean waiting time, a mean time between events). The loss in estimating
\\\theta\\ by the sample mean is \\A\\(\bar{X} - \theta)^2\\, and the
goal is a sample size at which the expected loss, the risk, is at most
`w` whatever the value of \\\theta\\. The first stage (`pilot_data`)
fixes the total sample size; the second stage (`additional_data`)
completes it. The estimate is the mean of all observations from both
stages. The procedure is risk based, so there is no confidence level.

## Usage

``` r
ts_br_exp_mean(pilot_data, additional_data = NULL, A, w, na.rm = TRUE)
```

## Arguments

- pilot_data:

  Numeric vector of the first stage observations, at least three, all
  positive. The whole first stage must be supplied at every call,
  because it alone determines the total sample size.

- additional_data:

  Numeric vector of the second stage observations collected so far, all
  positive, or `NULL` (the default) when none have been collected.

- A:

  The loss constant: the loss is \\A\\(\bar{X} - \theta)^2\\. A single
  positive number, in the units of `w` per squared unit of the data.

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

  The mean of all observations supplied, the estimate of \\\theta\\.

- `risk`:

  \\A \bar{X}^2 / n\\: the risk \\A \theta^2 / n\\ of the mean of a
  fixed sample of the current size, with \\\theta\\ replaced by its
  estimate. A descriptive estimate; the guarantee is the one stated in
  Details.

- `optimal_n`:

  The total sample size \\N\\ fixed by the first stage.

- `current_n`:

  The number of observations supplied, both stages together.

Attributes: `is_satisfied` (`TRUE` once `current_n` is at least
`optimal_n`), `current_n`, `criterion_n` (equal to `optimal_n`),
`loss_A`, `risk_bound`, `method`, and `notes` (always empty for this
function). There is no `conf_level` or `width`, and no
`variance_unusable` attribute: the total sample size depends on the
first stage mean, not on a variance, so a first stage of equal
observations is accepted.

## Details

The procedure is the two stage methodology that Mukhopadhyay and de
Silva (2009, Section 8.3.1) attribute to Mukhopadhyay and Pepe (2006),
with the improved constant that their Section 8.3.2 attributes to Zacks
and Mukhopadhyay (2006). With \\m \ge 3\\ first stage observations and
their mean \\\bar{X}\_m\\, the total sample size is (the book's Equation
8.3.9) \$\$N = \max\\m,\\ \lceil B \\ \bar{X}\_m^2 / w \rceil\\.\$\$ If
\\N = m\\ the first stage is enough; otherwise \\N - m\\ further
observations are collected, once, and \\\theta\\ is estimated by the
mean \\\bar{X}\_N\\ of all \\N\\ observations. \\N\\ depends on the
first stage alone. The row `optimal_n` is this total \\N\\, not the
number still to collect. The constant is (Equations 8.3.11 and 8.3.13)
\$\$B = 0.565 \times \frac{2 m (m + 1) A}{(m - 1)(m - 2)}.\$\$ (The book
writes the rounding as one more than the largest integer below \\B
\bar{X}\_m^2 / w\\, which is the ceiling.)

**What is guaranteed.** For the constant without the factor 0.565,
Theorem 8.3.1 of Mukhopadhyay and de Silva (2009) states that the risk
is at most `w` for every \\\theta\\, every \\m \ge 3\\, and every `A`;
but that procedure uses about twice the sample size \\n^\* = A \theta^2
/ w\\ at which a fixed sample has risk `w`. The factor 0.565 (Equation
8.3.13) is the improvement of Zacks and Mukhopadhyay (2006), found by
exact computation of the risk; the book reports that it keeps the risk
"only a shade under" the bound. That is a computed result, not a
theorem, and it was recomputed for this package: the risk relative to
`w` depends only on \\m\\ and \\n^\*\\, and exact calculation over first
stage sizes from 3 to 100 and values of \\n^\*\\ from 1 to 2,000 gives a
largest value of 0.990 (for a first stage of 14 or 15 observations with
\\n^\*\\ near 44). The bound holds throughout that range, with little to
spare by design. The guarantee rests on the exponential model.

**What it costs.** The expected total sample size exceeds \\n^\*\\ by a
factor that depends on the first stage size: about \\0.565 \times 2 (m +
1)^2 / ((m - 1)(m - 2))\\, which is 9.04 for \\m = 3\\, 1.90 for \\m =
10\\, 1.46 for \\m = 20\\, and tends to 1.13. See
[`ts_br_exp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean_pilot.md)
for the table and for the choice of first stage size, which the source
leaves open.

The data must be positive, as exponential observations are. Data on a
scale so large or so small that the square of a mean cannot be
represented in double precision (a mean beyond about \\10^{154}\\ or
below about \\10^{-154}\\) are an error that says to express them in
other units; no row of the result is ever infinite.

## References

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

Mukhopadhyay, N., & Pepe, W. (2006). Exact bounded risk estimation when
the terminal sample size and estimator are dependent: The exponential
case. *Sequential Analysis*, *25*(1), 85–101.
[doi:10.1080/07474940500452254](https://doi.org/10.1080/07474940500452254)
(As described in Mukhopadhyay and de Silva, 2009, Section 8.3.1,
Equations 8.3.9 and 8.3.11.)

Zacks, S., & Mukhopadhyay, N. (2006). Bounded risk estimation of the
exponential parameter in a two-stage sampling. *Sequential Analysis*,
*25*(4), 437–452.
[doi:10.1080/07474940600934896](https://doi.org/10.1080/07474940600934896)
(As described in Mukhopadhyay and de Silva, 2009, Section 8.3.2,
Equation 8.3.13.)

## See also

[`ts_br_exp_mean_pilot`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean_pilot.md)
for the first stage size;
[`seq_mr_exp_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean.md)
for purely sequential minimum risk estimation of an exponential mean.

Other ts_br_exp_mean:
[`ts_br_exp_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# An operations analyst wants the mean time between service calls
# (hours) estimated under squared error loss with A = 100, keeping the
# risk at or below 4.
set.seed(113)
m <- ts_br_exp_mean_pilot(A = 100, w = 4)
first_stage <- rexp(m, rate = 1 / 2)
check <- ts_br_exp_mean(first_stage, A = 100, w = 4)
check
#> Sequential check (two stage bounded risk estimation (exponential mean))
#> 
#>   Current n       : 20
#>   Criterion n     : 116
#> 
#>   mean              : 1.825926
#>   risk              : 16.67003
#>   optimal_n         : 116
#>   current_n         : 20
#> 
#>   Stopping rule met: FALSE

# The first stage fixed the total sample size; collect the remainder.
n_more <- attr(check, "criterion_n") - attr(check, "current_n")
second_stage <- rexp(n_more, rate = 1 / 2)
ts_br_exp_mean(first_stage, second_stage, A = 100, w = 4)
#> Sequential check (two stage bounded risk estimation (exponential mean))
#> 
#>   Current n       : 116
#>   Criterion n     : 116
#> 
#>   mean              : 1.977535
#>   risk              : 3.371244
#>   optimal_n         : 116
#>   current_n         : 116
#> 
#>   Stopping rule met: TRUE
```
