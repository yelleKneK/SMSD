# Pilot Sample Size for Bounded Risk Point Estimation of a Normal Mean

Returns the number of observations to collect before the stopping rule
of
[`seq_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md)
is checked for the first time. That procedure is the purely sequential
bounded risk point estimator of the mean of a normal population with
unknown variance (Mukhopadhyay and de Silva, 2009).

## Usage

``` r
seq_br_norm_mean_pilot(A, k, w, m_0 = 10, gamma = 1)
```

## Arguments

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

- m_0:

  The smallest pilot sample size allowed. A whole number of at least 2;
  default 10. A larger value reduces the amount by which the risk can
  exceed the bound when the bound is loose (see the Details of
  [`seq_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md));
  it costs at most `m_0` observations.

- gamma:

  Positive constant that sets how fast the pilot grows as the bound
  tightens; see Details. Default 1.

## Value

A length one integer giving the pilot sample size.

## Details

If the population standard deviation \\\sigma\\ were known, the smallest
sample size with risk at most `w` would be \\n^\* =
(B/w)^{2/k}\sigma^2\\, where \\B =
A\\2^{k/2}\\\Gamma\\(k+1)/2\\/\sqrt{\pi}\\ (see
[`seq_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md)).
The pilot must be chosen before anything is known about \\\sigma\\, so
it uses only the known constant \\(B/w)^{2/k}\\: \$\$m = \max\\m_0,\\
\lceil \\(B/w)^{2/k}\\^{1/(1 + \gamma)} \rceil\\.\$\$ With the default
\\\gamma = 1\\ the second term is the square root of that constant,
which grows more slowly than \\n^\*\\ as `w` shrinks. It therefore does
not presume a value of \\\sigma\\: a pilot equal to \\n^\*\\ evaluated
at \\\sigma = 1\\ is far too large when the standard deviation is well
below 1 in the units of the data, and far too small to protect against
early stopping when it is well above 1.

This is the package's pilot rule, of the same form as the pilot sample
sizes of the package's other sequential procedures (Chattopadhyay and
Kelley, 2016, Equation 24; Kelley, Darku, and Chattopadhyay, 2018,
Equation 21). It should not be read as the pilot size prescribed by
Mukhopadhyay and de Silva (2009).

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_br_norm_mean`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md)

Other seq_br_norm_mean:
[`seq_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Mean absolute error of at most 0.25 points on a satisfaction scale
seq_br_norm_mean_pilot(A = 1, k = 1, w = 0.25)
#> [1] 10
# Mean squared error of at most 0.01
seq_br_norm_mean_pilot(A = 1, k = 2, w = 0.01)
#> [1] 10
```
