# Pilot Sample Size for Sequential AIPE: Standardized Mean

Returns the number of observations to collect before the sequential
accuracy in parameter estimation (AIPE) stopping rule for the population
standardized mean \\\mu/\sigma\\ can be evaluated for the first time
with
[`seq_aipe_sm`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm.md).

## Usage

``` r
seq_aipe_sm_pilot(width, conf_level = 0.95, m_0 = 4)
```

## Arguments

- width:

  The target full width of the confidence interval for the population
  standardized mean (\\\omega\\ in the cited article): the upper limit
  minus the lower limit, in standard deviation units. A single positive
  number.

- conf_level:

  The confidence level of the interval, a single number strictly between
  0 and 1. Default 0.95.

- m_0:

  The smallest pilot size allowed, a whole number of at least 4 (the
  variance estimator needs four observations). Default 4. Raising it
  makes an unusable early variance estimate less likely; see
  [`seq_aipe_sm`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm.md).

## Value

A single integer: the pilot sample size.

## Details

The pilot size is \$\$m = \max\\m_0, \lceil 2 z\_{\alpha/2} / \omega
\rceil\\,\$\$ where \\\omega\\ is `width`, \\z\_{\alpha/2}\\ is the
upper \\\alpha/2\\ quantile of the standard normal distribution with
\\\alpha\\ equal to one minus `conf_level`, and \\m_0\\ is `m_0`
(Kelley, Darku, and Chattopadhyay, 2018, Equation 21). Whatever the
data, the stopping rule cannot be met with fewer than \\\lceil 2
z\_{\alpha/2} / \omega \rceil\\ observations, which is why the pilot is
at least that large.

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

## See also

[`seq_aipe_sm`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm.md)

Other seq_aipe_sm:
[`seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md),
[`seq_aipe_sm()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# A 95% interval for a standardized mean no wider than 0.20
seq_aipe_sm_pilot(width = 0.20)
#> [1] 20

# A wider target at 90% confidence, with at least 15 in the pilot
seq_aipe_sm_pilot(width = 0.50, conf_level = 0.90, m_0 = 15)
#> [1] 15
```
