# Pilot Sample Size for Minimum Risk Point Estimation of the Location of a Negative Exponential Distribution

Returns the number of observations to collect before the stopping rule
of
[`seq_mr_nexp_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean.md)
is checked for the first time. That procedure is the purely sequential
minimum risk point estimation (MRPE) procedure for the location
parameter (the lower endpoint) of a negative exponential distribution
(Mukhopadhyay and de Silva, 2009).

## Usage

``` r
seq_mr_nexp_mean_pilot(A, c, k, t, m_0 = 4, gamma = 1)
```

## Arguments

- A:

  The loss constant: the price of one unit of estimation error raised to
  the power `k`, in the same currency as `c`, so that the loss is
  \\A\\(X\_{(1)} - \theta)^k + c\\n^t\\. A positive number.

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

- m_0:

  The smallest pilot sample size allowed. A whole number of at least 2;
  default 4.

- gamma:

  Positive constant that sets how fast the pilot grows as the cost per
  observation falls; see Details. Default 1.

## Value

A length one integer giving the pilot sample size.

## Details

If the scale \\\sigma\\ of the distribution were known, the sample size
that minimizes the risk would be \\n^\* = K\\\sigma^{k/(t+k)}\\ with \\K
= \\A\\k\\\Gamma(k+1)/(c\\t)\\^{1/(t+k)}\\ (see
[`seq_mr_nexp_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean.md)).
The pilot must be chosen before anything is known about \\\sigma\\, so
it uses only the known constant \\K\\: \$\$m = \max\\m_0,\\ \lceil
K^{1/(1 + \gamma)} \rceil\\.\$\$ With the default \\\gamma = 1\\ the
second term is the square root of \\K\\, which grows more slowly than
\\n^\*\\ as the cost per observation falls. It therefore does not
presume a value of \\\sigma\\: a pilot equal to \\n^\*\\ evaluated at
\\\sigma = 1\\ already exceeds the best sample size whenever the scale
is below 1 in the units of the data.

This is the package's pilot rule. It has the form of the pilot sample
size of Chattopadhyay and Kelley (2016, Equation 24). It should not be
read as the pilot size prescribed by Mukhopadhyay and de Silva (2009).

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_mr_nexp_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean.md)

Other seq_mr_nexp_mean:
[`seq_mr_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean.md)

## Author

Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>, Ken Kelley
<kkelley@nd.edu>

## Examples

``` r
# Shortest possible delivery time (days): an error of one day is valued
# at 500 dollars and tracking one more order costs 0.25 dollars.
seq_mr_nexp_mean_pilot(A = 500, c = 0.25, k = 1, t = 1)
#> [1] 7
```
