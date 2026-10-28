# Pilot Sample Size for Sequential AIPE: Correlation Coefficient

Returns the number of paired observations to collect before the
sequential accuracy in parameter estimation (AIPE) stopping rule for a
population correlation coefficient can be evaluated for the first time
with
[`seq_aipe_r`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md).
The same pilot size applies to the Pearson, Kendall, and Spearman
coefficients, because it depends only on the target width and the
confidence level.

## Usage

``` r
seq_aipe_r_pilot(width, conf_level = 0.95, m_0 = 4)
```

## Arguments

- width:

  The target full width of the confidence interval for the population
  correlation (\\\omega\\ in the cited article): the upper limit minus
  the lower limit. A single positive number less than 2, the length of
  the range of a correlation.

- conf_level:

  The confidence level of the interval, a single number strictly between
  0 and 1. Default 0.95.

- m_0:

  The smallest pilot size allowed, a whole number of at least 4 (the
  variance estimator needs four pairs). Default 4. Raising it makes an
  unusable early variance estimate less likely; see
  [`seq_aipe_r`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md).

## Value

A single integer: the pilot number of pairs.

## Details

The pilot size is \$\$m = \max\\m_0, \lceil 2 z\_{\alpha/2} / \omega
\rceil\\,\$\$ where \\\omega\\ is `width`, \\z\_{\alpha/2}\\ is the
upper \\\alpha/2\\ quantile of the standard normal distribution with
\\\alpha\\ equal to one minus `conf_level`, and \\m_0\\ is `m_0`
(Kelley, Darku, and Chattopadhyay, 2019, Equation 9, where \\m_0 = 4\\).
The article's worked example, a 95% interval of width 0.10, has a pilot
of 40 pairs. Whatever the data, the stopping rule cannot be met with
fewer than \\\lceil 2 z\_{\alpha/2} / \omega \rceil\\ pairs, which is
why the pilot is at least that large.

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
accuracy in parameter estimation for population correlation
coefficients. *Psychological Methods*, *24*(4), 492–515.
[doi:10.1037/met0000203](https://doi.org/10.1037/met0000203)

## See also

[`seq_aipe_r`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md)

Other seq_aipe_r:
[`seq_aipe_r()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md),
[`seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# The worked example of Kelley, Darku, and Chattopadhyay (2019): a 95%
# interval no wider than 0.10
seq_aipe_r_pilot(width = 0.10)
#> [1] 40

# A wider target, with at least 30 pairs in the pilot
seq_aipe_r_pilot(width = 0.40, m_0 = 30)
#> [1] 30
```
