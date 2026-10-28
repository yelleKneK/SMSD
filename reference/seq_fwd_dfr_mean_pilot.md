# Pilot Sample Size for a Fixed Width Confidence Interval for a Mean, Distribution Free

Returns the number of observations to collect before the stopping rule
of
[`seq_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md)
is checked for the first time. That procedure is the purely sequential
fixed width confidence interval for the mean of a population whose
distribution is not specified (Mukhopadhyay and de Silva, 2009).

## Usage

``` r
seq_fwd_dfr_mean_pilot(d, conf_level = 0.95, m_0 = 10, gamma = 1)
```

## Arguments

- d:

  The half width of the confidence interval, in the units of the data:
  the interval reported is the sample mean plus and minus `d`, so its
  full width is \\2d\\. A positive number.

- conf_level:

  The confidence level (default 0.95).

- m_0:

  The smallest pilot sample size allowed. A whole number of at least 2;
  default 10. A larger value protects against stopping too early when
  the target is wide (see the Details of
  [`seq_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md));
  it costs at most `m_0` observations.

- gamma:

  Positive constant that sets how fast the pilot grows as the target
  narrows; see Details. Default 1.

## Value

A length one integer giving the pilot sample size.

## Details

If the population standard deviation \\\sigma\\ were known, the sample
size that gives the interval \\\bar{X}\_n \pm d\\ approximately the
stated confidence level would be \\n^\* = z\_{\alpha/2}^2 \sigma^2 /
d^2\\. The pilot must be chosen before anything is known about
\\\sigma\\, so it uses only the known constant \\(z\_{\alpha/2}/d)^2\\:
\$\$m = \max\\m_0,\\ \lceil (z\_{\alpha/2}/d)^{2/(1 + \gamma)}
\rceil\\.\$\$ With the default \\\gamma = 1\\ this is \\\max\\m_0,
\lceil z\_{\alpha/2}/d \rceil\\\\. Because of the \\1/n\\ term in the
stopping rule of
[`seq_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md),
that rule can never be met with fewer than \\z\_{\alpha/2}/d\\
observations, so nothing is lost by collecting that many at the outset.

With \\\gamma = 1\\ this is the pilot sample size of Kelley, Darku, and
Chattopadhyay (2018, Equation 21), whose full width \\\omega\\ is
\\2d\\; those authors present it as following Mukhopadhyay and de Silva
(2009). The argument `gamma` and the default value of `m_0` are the
package's.

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md)

Other seq_fwd_dfr_mean:
[`seq_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Mean weekly hours of remote work, to within half an hour, 90% confidence
seq_fwd_dfr_mean_pilot(d = 0.5, conf_level = 0.90)
#> [1] 10
seq_fwd_dfr_mean_pilot(d = 0.1, conf_level = 0.90)
#> [1] 17
```
