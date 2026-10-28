# Pilot Sample Size for Minimum Risk Point Estimation of the Coefficient of Variation

Returns the pilot (first stage) sample size for the purely sequential
minimum risk point estimation (MRPE) procedure for the population
coefficient of variation \\\kappa = \sigma / \mu\\ of Chattopadhyay and
Kelley (2016). It needs only the two cost constants, so it can be
computed before any data are collected.

## Usage

``` r
seq_mr_cv_pilot(A, c, m_0 = 4, gamma = 0.49)
```

## Arguments

- A:

  The loss constant: the price of one unit of squared estimation error,
  in the same currency as `c`. The article builds it as the amount one
  is willing to pay to estimate the coefficient of variation to within
  \\\epsilon\\, divided by \\\epsilon^2\\ (for example \\1000 / 0.05^2 =
  400000\\). A single positive number.

- c:

  The cost of sampling one observation. A single positive number.

- m_0:

  The smallest pilot size allowed: a whole number of at least 4 (the
  default), because the variance estimator uses the fourth central
  moment, which needs four observations.

- gamma:

  The tuning constant of the stopping rule's correction term, strictly
  between 0 and 1/2. The default, 0.49, is the value Chattopadhyay and
  Kelley (2016) suggest. Use the same value in
  [`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md).

## Value

A length one integer, the pilot sample size.

## Details

The pilot size is Equation 24 of Chattopadhyay and Kelley (2016), \$\$m
= \max\left\\m_0,\\ \left\lceil (A/c)^{1/(2 + 2\gamma)}
\right\rceil\right\\,\$\$ where \\\lceil\cdot\rceil\\ rounds up to the
next whole number. The stopping rule of
[`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
can never be met below \\(A/c)^{1/(2 + 2\gamma)}\\ observations,
whatever the data, which is why the pilot starts there. In the article's
worked example (\\A\\ = 400,000, \\c\\ = 75, \\\gamma\\ = 0.49; p. 637)
the pilot size is 18.

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

## See also

[`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
for the stopping rule and
[`seq_mr_cv_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)
to study the procedure before collecting data; the accuracy in parameter
estimation companion
[`seq_aipe_cv`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md);
in DMAR,
[`mr_cv`](https://yelleknek.github.io/DMAR/reference/mr_cv.html).

Other seq_mr_cv:
[`seq_mr_cv()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md),
[`seq_mr_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)

## Author

Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>, Ken Kelley
<kkelley@nd.edu>

## Examples

``` r
# The worked example of Chattopadhyay and Kelley (2016, p. 637): a home
# visit costs 75 dollars, and the team would pay 1,000 dollars to
# estimate the coefficient of variation to within 0.05.
seq_mr_cv_pilot(A = 1000 / 0.05^2, c = 75)
#> [1] 18
```
