# Pilot Sample Size for Sequential AIPE: Gini Index

Returns the number of observations to collect before the sequential
accuracy in parameter estimation (AIPE) stopping rule for the population
Gini index can be evaluated for the first time with
[`seq_aipe_gini`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md).

## Usage

``` r
seq_aipe_gini_pilot(width, conf_level = 0.95, m_0 = 4)
```

## Arguments

- width:

  The target full width of the confidence interval for the population
  Gini index: the upper limit minus the lower limit (twice the half
  width \\d\\ of Chattopadhyay and De, 2016). A single positive number
  less than 1, the length of the range of the Gini index.

- conf_level:

  The confidence level of the interval, a single number strictly between
  0 and 1. Default 0.95.

- m_0:

  The smallest pilot size allowed, a whole number of at least 4. Default
  4.

## Value

A single integer: the pilot sample size.

## Details

The pilot size is \$\$m = \max\\m_0, \lceil 2 z\_{\alpha/2} / \omega
\rceil\\,\$\$ where \\\omega\\ is `width`, \\z\_{\alpha/2}\\ is the
upper \\\alpha/2\\ quantile of the standard normal distribution with
\\\alpha\\ equal to one minus `conf_level`, and \\m_0\\ is `m_0`. This
is the pilot size of Chattopadhyay and De (2016), \\\max\\4, \lceil
z\_{\alpha/2} / d \rceil\\\\, written with the full width \\\omega =
2d\\ in place of their half width \\d\\. Whatever the data, the stopping
rule cannot be met with fewer than \\\lceil 2 z\_{\alpha/2} / \omega
\rceil\\ observations, which is why the pilot is at least that large.

## References

Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
pre-specified error bound. *Econometrics*, *4*(3), 30.
[doi:10.3390/econometrics4030030](https://doi.org/10.3390/econometrics4030030)

## See also

[`seq_aipe_gini`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md)

Other seq_aipe_gini:
[`seq_aipe_gini()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md),
[`seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# A 95% interval for a Gini index no wider than 0.10
seq_aipe_gini_pilot(width = 0.10)
#> [1] 40

# A 90% interval no wider than 0.04, a half width of 0.02
seq_aipe_gini_pilot(width = 0.04, conf_level = 0.90)
#> [1] 83
```
