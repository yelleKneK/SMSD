# Known Variance Sample Size for a Fixed Width Interval for a Normal Mean

Returns the sample size at which a confidence interval for the mean of a
normal distribution has half width `d` when the population variance is
known. This is the benchmark that the purely sequential
([`seq_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md))
and two stage
([`ts_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md))
fixed width procedures aim for without knowing the variance. It is
useful for planning with an assumed variance and for judging how much a
sequential procedure costs relative to knowing the variance in advance.

## Usage

``` r
seq_fwd_norm_mean_optimal(variance, d, conf_level = 0.95)
```

## Arguments

- variance:

  The population variance \\\sigma^2\\ (the square of the standard
  deviation, not the standard deviation itself). A single positive
  number, in the squared units of the data.

- d:

  The half width of the confidence interval (the margin of error), in
  the units of the data. A single positive number. The full width is
  \\2d\\.

- conf_level:

  The confidence level, strictly between 0 and 1. Default 0.95.

## Value

A length one integer of at least 1: the smallest sample size \\n\\ for
which the known variance interval \\\bar{X}\_n \pm d\\ has confidence
level at least `conf_level`, that is, \\\lceil C \rceil\\.

## Details

With known variance \\\sigma^2\\ the interval \\\bar{X}\_n \pm d\\
covers the mean with probability \\2\Phi(d \sqrt{n} / \sigma) - 1\\,
which reaches the confidence level \\1 - \alpha\\ when \\n\\ is at least
\$\$C = \frac{z\_{\alpha/2}^2 \\ \sigma^2}{d^2},\$\$ where
\\z\_{\alpha/2}\\ is the upper \\\alpha/2\\ point of the standard normal
distribution. \\C\\ is rarely a whole number; the function returns the
smallest whole sample size at which the coverage is at least the
confidence level, which is \\C\\ rounded up and never less than 1. (When
\\C\\ is within rounding error of a whole number, the coverage formula
decides between that number and the next.) The unrounded value is
`qnorm((1 - conf_level) / 2, lower.tail = FALSE)^2 * variance / d^2`.

Because the result is proportional to the variance, a planning value
that is too small gives a sample that is too small by the same factor.
When the variance is not known, use one of the procedures that estimate
it as the data arrive.

## References

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md)
and
[`ts_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md),
which reach the same target without knowing the variance;
[`seq_budget_optimal`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md)
for the corresponding known parameter allocation in the budget
constrained procedure.

Other seq_fwd_norm_mean:
[`seq_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md),
[`seq_fwd_norm_mean_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_pilot.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Mean score on a 100 point skills assessment to within 2 points, with
# 95% confidence, if the standard deviation were known to be 8.
seq_fwd_norm_mean_optimal(variance = 8^2, d = 2)
#> [1] 62

# The same target at 99% confidence.
seq_fwd_norm_mean_optimal(variance = 8^2, d = 2, conf_level = 0.99)
#> [1] 107
```
