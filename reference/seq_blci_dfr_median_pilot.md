# Pilot Sample Size for a Bounded Length Confidence Interval for a Median, Distribution Free

Returns the number of observations to collect before the stopping rule
of
[`seq_blci_dfr_median`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md)
is checked for the first time. That procedure is the purely sequential
bounded length confidence interval for the median of a population whose
distribution is not specified (Mukhopadhyay and de Silva, 2009).

## Usage

``` r
seq_blci_dfr_median_pilot(d, conf_level = 0.95, m_0 = 10, gamma = 1)
```

## Arguments

- d:

  Half of the largest interval length allowed, in the units of the data:
  sampling continues until the length of the interval is at most \\2d\\.
  A positive number.

- conf_level:

  The confidence level (default 0.95).

- m_0:

  The smallest pilot sample size allowed, in addition to the minimum the
  confidence level itself requires. A whole number of at least 2;
  default 10.

- gamma:

  Positive constant that sets how fast the pilot grows as the target
  narrows; see Details. Default 1.

## Value

A length one integer giving the pilot sample size.

## Details

The pilot is the largest of three numbers: \$\$m = \max\\n\_{\min},\\
m_0,\\ \lceil (z\_{\alpha/2}/d)^{2/(1 + \gamma)} \rceil\\.\$\$ The
first, \\n\_{\min}\\, is the smallest sample size at which an interval
formed from order statistics can reach the confidence level at all. The
widest such interval, from the smallest to the largest observation,
misses the median only when every observation falls on the same side of
it, so its coverage is \\1 - 2^{1-n}\\; \\n\_{\min}\\ is the smallest
\\n\\ for which this is at least `conf_level` (5, 6, and 8 observations
at 90%, 95%, and 99%).

The third is of lower order than the sample size the procedure will
need, which is approximately \\n^\* = (z\_{\alpha/2}/d)^2 / \\4
f(\theta)^2\\\\ for a population with density \\f\\ at its median
\\\theta\\. The density is unknown before sampling, so the pilot uses
only the known constant \\(z\_{\alpha/2}/d)^2\\; with the default
\\\gamma = 1\\ the third term is \\\lceil z\_{\alpha/2}/d \rceil\\. A
pilot equal to \\(z\_{\alpha/2}/d)^2\\ itself would be several times too
large for a population that is concentrated near its median in the units
of the data.

This is the package's pilot rule, of the same form as the pilot sample
size of Kelley, Darku, and Chattopadhyay (2018, Equation 21). It should
not be read as the pilot size prescribed by Mukhopadhyay and de Silva
(2009).

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_blci_dfr_median`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md)

Other seq_blci_dfr_median:
[`seq_blci_dfr_median()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Median household income (thousands of dollars), interval no longer
# than 10 thousand
seq_blci_dfr_median_pilot(d = 5, conf_level = 0.95)
#> [1] 10
# A tight target and a high confidence level
seq_blci_dfr_median_pilot(d = 0.1, conf_level = 0.99)
#> [1] 26
```
