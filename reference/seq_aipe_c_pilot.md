# Pilot Sample Size for Sequential AIPE: Contrast

Returns the number of observations to collect from each group before the
sequential accuracy in parameter estimation (AIPE) stopping rule for a
contrast of group means, \\\psi = \sum\_{k=1}^{K} c_k \mu_k\\, can be
evaluated for the first time with
[`seq_aipe_c`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md).

## Usage

``` r
seq_aipe_c_pilot(width, conf_level = 0.95, m_0 = 4)
```

## Arguments

- width:

  The target full width of the confidence interval for the contrast
  (\\\omega\\ in the cited article): the upper limit minus the lower
  limit, in the units of the outcome. A single positive number.

- conf_level:

  The confidence level of the interval, a single number strictly between
  0 and 1. Default 0.95.

- m_0:

  The smallest pilot size allowed per group, a whole number of at
  least 4. Default 4.

## Value

A single integer: the pilot sample size per group.

## Details

The pilot size is \$\$m = \max\\m_0, \lceil 2 z\_{\alpha/2} / \omega
\rceil\\,\$\$ where \\\omega\\ is `width`, \\z\_{\alpha/2}\\ is the
upper \\\alpha/2\\ quantile of the standard normal distribution with
\\\alpha\\ equal to one minus `conf_level`, and \\m_0\\ is `m_0`
(Kelley, Darku, and Chattopadhyay, 2018, Equation 21). Whatever the
data, the stopping rule cannot be met with fewer than \\\lceil 2
z\_{\alpha/2} / \omega \rceil\\ observations per group, which is why the
pilot is at least that large. Because a contrast is not scale free, the
pilot size depends on the units of `width`: see
[`seq_aipe_c`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md).

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

## See also

[`seq_aipe_c`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md)

Other seq_aipe_c:
[`seq_aipe_c()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md),
[`seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# A 95% interval for a contrast no wider than 0.40 outcome units
seq_aipe_c_pilot(width = 0.40)
#> [1] 10

# The same target with at least 20 per group in the pilot
seq_aipe_c_pilot(width = 0.40, m_0 = 20)
#> [1] 20
```
