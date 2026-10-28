# Pilot Sample Size for a Two Stage Fixed Width Interval for a Negative Exponential Location

Returns the first stage (pilot) sample size for the two stage fixed
width confidence interval for the location (the threshold, or minimum
possible value) of a negative exponential distribution, evaluated by
[`ts_fwd_nexp_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md).
The interval is \\\[X\_{(1)} - d,\\ X\_{(1)}\]\\, with \\X\_{(1)}\\ the
smallest observation: `d` is the full width of the interval.

## Usage

``` r
ts_fwd_nexp_mean_pilot(d, conf_level = 0.95, m_0 = 10, gamma = 1)
```

## Arguments

- d:

  The full width of the confidence interval, in the units of the data. A
  single positive number.

- conf_level:

  The confidence level, strictly between 0 and 1. Default 0.95.

- m_0:

  The smallest first stage sample size to return. A whole number of at
  least 2 (the scale estimate needs two observations). Default 10; see
  Details.

- gamma:

  Positive tuning constant of the pilot rule. Larger values give a
  smaller first stage. Default 1, for which the rule is \\\lceil \sqrt{a
  / d} \rceil\\.

## Value

A length one integer: the first stage sample size.

## Details

The pilot size is \$\$m = \max\\m_0,\\ \lceil (a / d)^{1/(1 + \gamma)}
\rceil\\, \qquad a = \log(1/\alpha),\$\$ with \\\alpha = 1 - \\
`conf_level`. The sample size that would be used if the scale \\\sigma\\
were known is \\C = a \sigma / d\\; the rule raises the constant \\a /
d\\ to the power \\1/(1 + \gamma)\\, so the first stage grows as the
target narrows but more slowly than \\C\\ does. This is the form of
pilot that Chattopadhyay and Kelley (2016, Equation 24) recommend,
following Mukhopadhyay and de Silva (2009), for a purely sequential
procedure. Using it here, with the known scale sample size for \\\sigma
= 1\\ as the constant, is this package's choice and not part of the two
stage procedure, which is valid for any first stage of at least two
observations.

The minimum `m_0` matters more than the rule. The procedure pays for not
knowing the scale through an *F* quantile on \\2\\ and \\2m - 2\\
degrees of freedom, so the expected total sample size is about
\\F\_{2,2m-2,\alpha} / a\\ times \\C\\. At a confidence level of 0.95
that ratio is 6.34 for \\m = 2\\, 2.32 for \\m = 3\\, 1.49 for \\m =
5\\, 1.19 for \\m = 10\\, and 1.08 for \\m = 20\\. The default
`m_0 = 10` holds the expected total to about 1.2 times \\C\\ at 0.95
(about 1.3 times \\C\\ at 0.99). A first stage larger than \\C\\ is
wasted, so lower `m_0` only when even ten observations may be more than
the study needs.

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`ts_fwd_nexp_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md)

Other ts_fwd_nexp_mean:
[`ts_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Shortest possible task completion time, to within half a second.
ts_fwd_nexp_mean_pilot(d = 0.5, conf_level = 0.95)
#> [1] 10

# A narrow target: the rule exceeds the minimum.
ts_fwd_nexp_mean_pilot(d = 0.01, conf_level = 0.99)
#> [1] 22
```
