# Pilot Sample Size for Minimum Risk Point Estimation of a Normal Mean

Returns the number of observations to collect before the stopping rule
of
[`seq_mr_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean.md)
is checked for the first time. That procedure is the purely sequential
minimum risk point estimation (MRPE) procedure for the mean of a normal
population with unknown variance (Robbins, 1959; Mukhopadhyay and de
Silva, 2009, Section 6.4).

## Usage

``` r
seq_mr_norm_mean_pilot(
  A = NULL,
  c,
  price = NULL,
  epsilon = NULL,
  m_0 = 4,
  gamma = 1
)
```

## Arguments

- A:

  The loss constant: the price of one unit of squared estimation error,
  in the same currency as `c`, so that the loss is \\A(\bar{X}\_n -
  \mu)^2 + c\\n\\. A positive number. Either supply `A` directly, or
  supply `price` and `epsilon`, in which case \\A =
  \mathrm{price}/\mathrm{epsilon}^2\\.

- c:

  The cost of sampling one observation. A positive number.

- price, epsilon:

  An alternative way to state `A` (Chattopadhyay and Kelley, 2016):
  `price` is the amount one is willing to pay for an estimate within
  `epsilon` of the mean, and \\A = \mathrm{price}/\mathrm{epsilon}^2\\.
  Both positive numbers; supply both or neither.

- m_0:

  The smallest pilot sample size allowed. A whole number of at least 2;
  default 4.

- gamma:

  Positive constant that sets how fast the pilot grows as the cost per
  observation falls; see Details. Default 1.

## Value

A length one integer giving the pilot sample size.

## Details

If the population standard deviation \\\sigma\\ were known, the sample
size that minimizes the risk would be \\n^\* = (A/c)^{1/2}\sigma\\
(Mukhopadhyay and de Silva, 2009, Equation 6.4.3). The pilot must be
chosen before anything is known about \\\sigma\\, so it uses only the
known constant \\(A/c)^{1/2}\\: \$\$m = \max\\m_0,\\ \lceil
(A/c)^{1/(2 + 2\gamma)} \rceil\\.\$\$ With the default \\\gamma = 1\\
the second term is \\(A/c)^{1/4}\\, which grows more slowly than
\\n^\*\\ as the cost per observation falls. It therefore does not
presume a value of \\\sigma\\: a pilot equal to \\n^\*\\ evaluated at
\\\sigma = 1\\ already exceeds the best sample size whenever the
standard deviation is below 1 in the units of the data, and every
observation beyond the best sample size adds cost without a matching
gain in accuracy.

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

Robbins, H. (1959). Sequential estimation of the mean of a normal
population. In U. Grenander (Ed.), *Probability and statistics: The
Harald Cramer volume* (pp. 235–245). Almqvist & Wiksell.

## See also

[`seq_mr_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean.md)

Other seq_mr_norm_mean:
[`seq_mr_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean.md)

## Author

Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>, Ken Kelley
<kkelley@nd.edu>

## Examples

``` r
# An error of one point (squared) in a mean test score is valued at
# 2,500 dollars; testing one more student costs 4 dollars.
seq_mr_norm_mean_pilot(A = 2500, c = 4)
#> [1] 5
# The same loss constant stated as 10,000 dollars for accuracy within 2
seq_mr_norm_mean_pilot(c = 4, price = 10000, epsilon = 2)
#> [1] 5
```
