# Pilot Sample Size for a Two Stage Fixed Width Interval for a Distribution Free Mean

Returns the first stage (pilot) sample size for the two stage fixed
width confidence interval for a mean without a distributional
assumption, evaluated by
[`ts_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md).
The interval is \\\bar{X} \pm d\\: `d` is the half width (the margin of
error), so the full width is \\2d\\.

## Usage

``` r
ts_fwd_dfr_mean_pilot(d, conf_level = 0.95, m_0 = 30, gamma = 1)
```

## Arguments

- d:

  The half width of the confidence interval (the margin of error), in
  the units of the data. A single positive number.

- conf_level:

  The nominal confidence level, strictly between 0 and 1. Default 0.95.

- m_0:

  The smallest first stage sample size to return. A whole number of at
  least 2. Default 30; see Details before lowering it.

- gamma:

  Positive tuning constant of the pilot rule. Larger values give a
  smaller first stage. Default 1, for which the rule is \\\lceil
  z\_{\alpha/2} / d \rceil\\.

## Value

A length one integer: the first stage sample size.

## Details

The pilot size is \$\$m = \max\\m_0,\\ \lceil (z\_{\alpha/2} /
d)^{2/(1 + \gamma)} \rceil\\,\$\$ where \\z\_{\alpha/2}\\ is the
standard normal quantile for the confidence level. The sample size that
would be used if the variance \\\sigma^2\\ were known is \\C =
z\_{\alpha/2}^2 \sigma^2 / d^2\\; the rule raises the constant
\\(z\_{\alpha/2} / d)^2\\ to the power \\1/(1 + \gamma)\\, so the first
stage grows as the target narrows but more slowly than \\C\\ does. This
is the form of pilot that Chattopadhyay and Kelley (2016, Equation 24)
recommend, following Mukhopadhyay and de Silva (2009), for a purely
sequential procedure. Using it here, with the known variance sample size
for \\\sigma = 1\\ as the constant, is this package's choice.

**The first stage size sets the coverage of this procedure**, so the
minimum `m_0` is not a detail. The procedure's confidence level is
asymptotic: it treats the first stage variance as if it were the
population variance. For normal data and a narrow target the coverage of
the resulting interval is \\P(\|t\_{m-1}\| \le z\_{\alpha/2})\\, with
\\t\_{m-1}\\ a *t* variable on \\m - 1\\ degrees of freedom. At a
nominal level of 0.95 that is 0.700 for \\m = 2\\, 0.878 for \\m = 5\\,
0.918 for \\m = 10\\, 0.935 for \\m = 20\\, 0.940 for \\m = 30\\, 0.944
for \\m = 50\\, and 0.947 for \\m = 100\\; for skewed data it is lower
(see
[`ts_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md)).
The default `m_0 = 30` keeps the shortfall for normal data to about one
percentage point at 0.95. Lowering `m_0` lowers the coverage; raise it
when the data may be skewed or heavy tailed.

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`ts_fwd_dfr_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md)

Other ts_fwd_dfr_mean:
[`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Mean weekly hours of unpaid overtime, to within half an hour.
ts_fwd_dfr_mean_pilot(d = 0.5, conf_level = 0.95)
#> [1] 30

# A narrow target on a standardized scale: the rule exceeds the minimum.
ts_fwd_dfr_mean_pilot(d = 0.05, conf_level = 0.95)
#> [1] 40
```
