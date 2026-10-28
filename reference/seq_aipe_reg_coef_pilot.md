# Pilot Sample Size for Sequential AIPE: Simple Linear Regression Slope

Returns the number of paired observations to collect before the
sequential accuracy in parameter estimation (AIPE) stopping rule for the
slope of a simple linear regression, \\\beta_1 =
\sigma\_{XY}/\sigma_X^2\\, can be evaluated for the first time with
[`seq_aipe_reg_coef`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md).

## Usage

``` r
seq_aipe_reg_coef_pilot(width, conf_level = 0.95, m_0 = 4)
```

## Arguments

- width:

  The target full width of the confidence interval for the slope
  (\\\omega\\ in the cited article): the upper limit minus the lower
  limit, in units of the outcome per unit of the predictor. A single
  positive number.

- conf_level:

  The confidence level of the interval, a single number strictly between
  0 and 1. Default 0.95.

- m_0:

  The smallest pilot size allowed, a whole number of at least 4 (the
  variance estimator needs four pairs). Default 4. Raising it makes an
  unusable early variance estimate less likely; see
  [`seq_aipe_reg_coef`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md).

## Value

A single integer: the pilot number of pairs.

## Details

The pilot size is \$\$m = \max\\m_0, \lceil 2 z\_{\alpha/2} / \omega
\rceil\\,\$\$ where \\\omega\\ is `width`, \\z\_{\alpha/2}\\ is the
upper \\\alpha/2\\ quantile of the standard normal distribution with
\\\alpha\\ equal to one minus `conf_level`, and \\m_0\\ is `m_0`
(Kelley, Darku, and Chattopadhyay, 2018, Equation 21). Whatever the
data, the stopping rule cannot be met with fewer than \\\lceil 2
z\_{\alpha/2} / \omega \rceil\\ pairs, which is why the pilot is at
least that large. Because the slope is not scale free, the pilot size
depends on the units of `width`: see
[`seq_aipe_reg_coef`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md).

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

## See also

[`seq_aipe_reg_coef`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md)

Other seq_aipe_reg_coef:
[`seq_aipe_reg_coef()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md),
[`seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# A 95% interval for a slope no wider than 0.20
seq_aipe_reg_coef_pilot(width = 0.20)
#> [1] 20

# The same target with at least 30 pairs in the pilot
seq_aipe_reg_coef_pilot(width = 0.20, m_0 = 30)
#> [1] 30
```
