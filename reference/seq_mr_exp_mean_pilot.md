# Pilot Sample Size for Minimum Risk Point Estimation of an Exponential Mean

Returns the number of observations to collect before the stopping rule
of
[`seq_mr_exp_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean.md)
is checked for the first time. That procedure is the purely sequential
minimum risk point estimation (MRPE) procedure for the mean of an
exponential population (Mukhopadhyay and de Silva, 2009).

## Usage

``` r
seq_mr_exp_mean_pilot(A, c, m_0 = 4, gamma = 1)
```

## Arguments

- A:

  The loss constant: the price of one unit of squared estimation error,
  in the same currency as `c`, so that the loss is \\A(\bar{X}\_n -
  \lambda)^2 + c\\n\\. A positive number.

- c:

  The cost of sampling one observation. A positive number.

- m_0:

  The smallest pilot sample size allowed. A whole number of at least 2;
  default 4.

- gamma:

  Positive constant that sets how fast the pilot grows as the cost per
  observation falls; see Details. Default 1.

## Value

A length one integer giving the pilot sample size.

## Details

If the population mean \\\lambda\\ were known, the sample size that
minimizes the risk would be \\n^\* = (A/c)^{1/2}\lambda\\ (the standard
deviation of an exponential distribution equals its mean). The pilot
must be chosen before anything is known about \\\lambda\\, so it uses
only the known constant \\(A/c)^{1/2}\\: \$\$m = \max\\m_0,\\ \lceil
(A/c)^{1/(2 + 2\gamma)} \rceil\\.\$\$ With the default \\\gamma = 1\\
the second term is \\(A/c)^{1/4}\\, which grows more slowly than
\\n^\*\\ as the cost per observation falls. It therefore does not
presume a value of \\\lambda\\: a pilot equal to \\n^\*\\ evaluated at
\\\lambda = 1\\ already exceeds the best sample size whenever the mean
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

[`seq_mr_exp_mean`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean.md)

Other seq_mr_exp_mean:
[`seq_mr_exp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean.md)

## Author

Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>, Ken Kelley
<kkelley@nd.edu>

## Examples

``` r
# An error of one minute (squared) in the mean time to resolve a support
# ticket is valued at 100 dollars; timing one more ticket costs 0.04.
seq_mr_exp_mean_pilot(A = 100, c = 0.04)
#> [1] 8
```
