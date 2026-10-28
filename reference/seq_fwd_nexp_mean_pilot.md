# Pilot Sample Size for a Fixed Width Confidence Interval for the Location of a Negative Exponential Distribution

Returns the number of observations to collect before the stopping rule
of
[`seq_fwd_nexp_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean.md)
is checked for the first time. That procedure is the purely sequential
fixed width confidence interval for the location parameter (the lower
endpoint, or minimum guaranteed value) of a negative exponential
distribution (Mukhopadhyay and de Silva, 2009).

## Usage

``` r
seq_fwd_nexp_mean_pilot(d, conf_level = 0.95, m_0 = 10, gamma = 1)
```

## Arguments

- d:

  The full width of the confidence interval, in the units of the data:
  the interval reported runs from the sample minimum minus `d` up to the
  sample minimum. A positive number. (In the fixed width procedures for
  a mean, `d` is a half width; here the interval is one-sided, so `d` is
  its whole width.)

- conf_level:

  The confidence level (default 0.95).

- m_0:

  The smallest pilot sample size allowed. A whole number of at least 2;
  default 10. A larger value protects against stopping too early when
  the target is wide (see the Details of
  [`seq_fwd_nexp_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean.md));
  it costs at most `m_0` observations.

- gamma:

  Positive constant that sets how fast the pilot grows as the target
  narrows; see Details. Default 1.

## Value

A length one integer giving the pilot sample size.

## Details

If the scale \\\sigma\\ of the distribution were known, the sample size
that gives the interval \\\[X\_{(1)} - d,\\ X\_{(1)}\]\\ the stated
confidence level would be \\n^\* = \log(1/\alpha)\\\sigma / d\\. The
pilot must be chosen before anything is known about \\\sigma\\, so it
uses only the known constant \\\log(1/\alpha)/d\\: \$\$m = \max\\m_0,\\
\lceil (\log(1/\alpha)/d)^{1/(1 + \gamma)} \rceil\\.\$\$ With the
default \\\gamma = 1\\ this is the square root of that constant, which
grows more slowly than \\n^\*\\ as `d` shrinks and so does not presume a
value of \\\sigma\\.

This is the package's pilot rule, of the same form as the pilot sample
sizes of the package's other sequential procedures (Chattopadhyay and
Kelley, 2016, Equation 24). It should not be read as the pilot size
prescribed by Mukhopadhyay and de Silva (2009).

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_fwd_nexp_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean.md)

Other seq_fwd_nexp_mean:
[`seq_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Shortest possible response time (seconds), to within 0.05 seconds
seq_fwd_nexp_mean_pilot(d = 0.05, conf_level = 0.95)
#> [1] 10
seq_fwd_nexp_mean_pilot(d = 0.005, conf_level = 0.95)
#> [1] 25
```
