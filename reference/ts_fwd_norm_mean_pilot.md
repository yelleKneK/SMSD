# Pilot Sample Size for a Two Stage Fixed Width Interval for a Normal Mean

Returns the first stage (pilot) sample size for Stein's two stage fixed
width confidence interval for the mean of a normal distribution,
evaluated by
[`ts_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md).
The interval is \\\bar{X} \pm d\\: `d` is the half width (the margin of
error), so the full width is \\2d\\.

## Usage

``` r
ts_fwd_norm_mean_pilot(d, conf_level = 0.95, m_0 = 10, gamma = 1)
```

## Arguments

- d:

  The half width of the confidence interval (the margin of error), in
  the units of the data. A single positive number.

- conf_level:

  The confidence level, strictly between 0 and 1. Default 0.95.

- m_0:

  The smallest first stage sample size to return. A whole number of at
  least 2 (a variance needs two observations). Default 10; see Details.

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
z\_{\alpha/2}^2 \sigma^2 / d^2\\ (see
[`seq_fwd_norm_mean_optimal`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_optimal.md));
the rule raises the constant \\(z\_{\alpha/2} / d)^2\\ to the power
\\1/(1 + \gamma)\\, so the first stage grows as the target narrows but
more slowly than \\C\\ does. This is the form of pilot that
Chattopadhyay and Kelley (2016, Equation 24) recommend, following
Mukhopadhyay and de Silva (2009), for a purely sequential procedure.
Using it here, with the known variance sample size for \\\sigma = 1\\ as
the constant, is this package's choice and not part of Stein's
procedure, which is valid for any first stage of at least two
observations.

The minimum `m_0` matters more than the rule. Stein's procedure pays for
not knowing the variance through the *t* quantile on \\m - 1\\ degrees
of freedom, so the expected total sample size is about
\\t\_{m-1,\alpha/2}^2 / z\_{\alpha/2}^2\\ times \\C\\. At a confidence
level of 0.95 that ratio is 42.0 for \\m = 2\\, 4.82 for \\m = 3\\, 2.01
for \\m = 5\\, 1.33 for \\m = 10\\, 1.14 for \\m = 20\\, and 1.09 for
\\m = 30\\. The default `m_0 = 10` holds the expected total to about
1.33 times \\C\\ at 0.95 (about 1.6 times \\C\\ at 0.99). A first stage
larger than \\C\\ is wasted, so lower `m_0` only when even ten
observations may be more than the study needs.

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

Stein, C. (1945). A two-sample test for a linear hypothesis whose power
is independent of the variance. *The Annals of Mathematical Statistics*,
*16*(3), 243–258.
[doi:10.1214/aoms/1177731088](https://doi.org/10.1214/aoms/1177731088)

## See also

[`ts_fwd_norm_mean`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md);
[`seq_fwd_norm_mean_optimal`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean_optimal.md)
for the known variance sample size.

Other ts_fwd_norm_mean:
[`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Mean score on a 100 point skills assessment, to within 2 points.
ts_fwd_norm_mean_pilot(d = 2, conf_level = 0.95)
#> [1] 10

# A narrow target on a standardized scale: the rule exceeds the minimum.
ts_fwd_norm_mean_pilot(d = 0.1, conf_level = 0.95)
#> [1] 20
```
