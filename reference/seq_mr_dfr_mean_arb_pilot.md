# Pilot Sample Size for Minimum Risk Point Estimation of a Mean, Distribution Free

Returns the number of observations to collect before the stopping rule
of
[`seq_mr_dfr_mean_arb`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md)
is checked for the first time. That procedure is the purely sequential
minimum risk point estimation (MRPE) procedure for the mean of a
population whose distribution is not specified and whose variance is
unknown and unrestricted (Mukhopadhyay and de Silva, 2009).

## Usage

``` r
seq_mr_dfr_mean_arb_pilot(A, c, m_0 = 4, gamma = 1)
```

## Arguments

- A:

  The loss constant: the price of one unit of squared estimation error,
  in the same currency as `c`, so that the loss is \\A(\bar{X}\_n -
  \mu)^2 + c\\n\\. A positive number.

- c:

  The cost of sampling one observation. A positive number.

- m_0:

  The smallest pilot sample size allowed. A whole number of at least 2;
  default 4.

- gamma:

  The positive constant of the correction term \\n^{-\gamma}\\ in the
  stopping rule; default 1. Use the value passed to
  [`seq_mr_dfr_mean_arb`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md).

## Value

A length one integer giving the pilot sample size.

## Details

The pilot is \$\$m = \max\\m_0,\\ \lceil (A/c)^{1/(2 + 2\gamma)}
\rceil\\,\$\$ the pilot sample size of Chattopadhyay and Kelley (2016,
Equation 24), who give it for a stopping rule with a correction term of
the same kind and present it as following Mukhopadhyay and de Silva
(2009). The stopping rule of
[`seq_mr_dfr_mean_arb()`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md),
\\n \ge (A/c)^{1/2}(S_n + n^{-\gamma})\\, can never be met with fewer
than \\(A/c)^{1/(2 + 2\gamma)}\\ observations, whatever the data, so
nothing is lost by collecting that many at the outset. Use the same
`gamma` here and in
[`seq_mr_dfr_mean_arb()`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md).
The default values of `m_0` and `gamma` are the package's.

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_mr_dfr_mean_arb`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md)

Other seq_mr_dfr_mean_arb:
[`seq_mr_dfr_mean_arb()`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md)

## Author

Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>, Ken Kelley
<kkelley@nd.edu>

## Examples

``` r
# An error of one dollar (squared) in mean order value is worth 400
# dollars to a retailer; auditing one more order costs 0.25 dollars.
seq_mr_dfr_mean_arb_pilot(A = 400, c = 0.25)
#> [1] 7
```
