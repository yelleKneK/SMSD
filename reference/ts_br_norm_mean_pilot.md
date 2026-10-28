# Pilot Sample Size for Two Stage Bounded Risk Estimation of a Normal Mean

Returns the first stage (pilot) sample size for the two stage bounded
risk point estimation procedure for the mean of a normal distribution,
evaluated by
[`ts_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md).
The loss in estimating the mean \\\mu\\ by the sample mean is
\\A\\\|\bar{X} - \mu\|^k\\, and the procedure keeps its expectation, the
risk, at or below `w`.

## Usage

``` r
ts_br_norm_mean_pilot(A, k, w, m_0 = 10, gamma = 1)
```

## Arguments

- A:

  The loss constant: the loss is \\A\\\|\bar{X} - \mu\|^k\\. A single
  positive number, in the units of `w` per unit of error to the power
  `k`.

- k:

  The error index: 1 for absolute error loss, 2 for squared error loss.

- w:

  The risk bound: the largest acceptable expected loss. A single
  positive number.

- m_0:

  The smallest first stage sample size to return. A whole number of at
  least \\k + 2\\. Default 10; see Details.

- gamma:

  Positive tuning constant of the pilot rule. Larger values give a
  smaller first stage. Default 1.

## Value

A length one integer: the first stage sample size.

## Details

The pilot size is \$\$m = \max\\m_0,\\ \lceil ((B / w)^{2/k})^{1/(1 +
\gamma)} \rceil\\, \qquad B = A \\ \frac{2^{k/2}}{\sqrt{\pi}} \\
\Gamma\\\left(\frac{k + 1}{2}\right),\$\$ where \\B = A\\E\|Z\|^k\\ for
a standard normal \\Z\\, which is \\A\\ for \\k = 2\\ and
\\A\sqrt{2/\pi}\\ for \\k = 1\\. (For \\k = 2\\ and the default \\\gamma
= 1\\ the rule is therefore \\\lceil \sqrt{A / w} \rceil\\.) The sample
size that would be used if the variance \\\sigma^2\\ were known is
\\n^\* = (B / w)^{2/k} \sigma^2\\; the rule raises the constant \\(B /
w)^{2/k}\\ to the power \\1/(1 + \gamma)\\, so the first stage grows
with the size of the problem but more slowly than \\n^\*\\ does. This is
the form of pilot that Chattopadhyay and Kelley (2016, Equation 24)
recommend, following Mukhopadhyay and de Silva (2009), for a purely
sequential procedure. Using it here, with \\n^\*\\ for \\\sigma = 1\\ as
the constant, is this package's choice and not part of the two stage
procedure. A pilot equal to \\(B / w)^{2/k}\\ itself, the optimal sample
size for \\\sigma = 1\\, would exceed \\n^\*\\ whenever \\\sigma\\ is
below 1. The rule's first stage exceeds \\n^\*\\ only when \\\sigma^2\\
is below \\((B / w)^{2/k})^{-\gamma/(1 + \gamma)}\\, a threshold that
shrinks as the problem grows.

The minimum `m_0` matters more than the rule. The procedure is defined
only for a first stage of at least \\k + 2\\ observations, and it pays
for not knowing the variance through a constant that is large when the
first stage is small: the expected total sample size is about
\\b_m^{2/k}\\ times \\n^\*\\ (see
[`ts_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md)
for \\b_m\\). For squared error (\\k = 2\\) that ratio is \\(m - 1)/(m -
3)\\: 3 for \\m = 4\\, 2 for \\m = 5\\, 1.29 for \\m = 10\\, and 1.12
for \\m = 20\\. For absolute error (\\k = 1\\) it is 3.14 for \\m = 3\\,
1.57 for \\m = 5\\, 1.20 for \\m = 10\\, and 1.09 for \\m = 20\\. The
default `m_0 = 10` holds the expected total to about 1.3 times \\n^\*\\.
A first stage larger than \\n^\*\\ is wasted, so lower `m_0` only when
even ten observations may be more than the study needs.

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`ts_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md)

Other ts_br_norm_mean:
[`ts_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Mean waiting time, squared error loss, mean squared error at most 4.
ts_br_norm_mean_pilot(A = 1, k = 2, w = 4)
#> [1] 10

# A demanding bound: the rule exceeds the minimum.
ts_br_norm_mean_pilot(A = 1, k = 2, w = 0.001)
#> [1] 32
```
