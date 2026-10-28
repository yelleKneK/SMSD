# Pilot Sample Size for Two Stage Bounded Risk Estimation of an Exponential Mean

Returns the first stage (pilot) sample size for the two stage bounded
risk point estimation procedure for the mean of an exponential
distribution, evaluated by
[`ts_br_exp_mean`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md).
The loss in estimating the mean \\\theta\\ by the sample mean is
\\A\\(\bar{X} - \theta)^2\\, and the procedure keeps its expectation,
the risk, at or below `w`.

## Usage

``` r
ts_br_exp_mean_pilot(A, w, m_0 = 20, gamma = 1)
```

## Arguments

- A:

  The loss constant: the loss is \\A\\(\bar{X} - \theta)^2\\. A single
  positive number, in the units of `w` per squared unit of the data.

- w:

  The risk bound: the largest acceptable expected loss. A single
  positive number.

- m_0:

  The smallest first stage sample size to return. A whole number of at
  least 3. Default 20; see Details.

- gamma:

  Positive tuning constant of the pilot rule. Larger values give a
  smaller first stage. Default 1, for which the rule is \\\lceil \sqrt{A
  / w} \rceil\\.

## Value

A length one integer: the first stage sample size.

## Details

The source of the procedure (Mukhopadhyay and de Silva, 2009, Section
8.3.1) requires only that the first stage have at least three
observations and gives no rule for choosing its size. The choice does
not affect the risk bound, but it decides how many observations the
procedure uses. Writing \\n^\* = A \theta^2 / w\\ for the sample size at
which the risk \\A \theta^2 / n\\ of a fixed sample equals `w`, the
expected total sample size for a first stage of \\m\\ observations is,
when \\n^\*\\ is well above \\m\\, \$\$E\[N\] / n^\* \approx 0.565
\times \frac{2 (m + 1)^2}{(m - 1)(m - 2)},\$\$ where 0.565 is the
constant of Zacks and Mukhopadhyay (2006; the book's Equation 8.3.13)
and the rest is the expectation of the stopping variable (Equations
8.3.9 and 8.3.11).

|                   |      |      |      |      |      |      |      |
|-------------------|------|------|------|------|------|------|------|
| \\m\\             | 3    | 5    | 10   | 15   | 20   | 30   | 50   |
| \\E\[N\] / n^\*\\ | 9.04 | 3.39 | 1.90 | 1.59 | 1.46 | 1.34 | 1.25 |

The ratio falls toward 1.13 as \\m\\ grows. (Exact calculation agrees
with the approximation to within 1 percent once \\n^\*\\ is at least
five times \\m\\, and to within 6 percent when \\n^\*\\ is as small as
\\m\\.) A first stage of three observations therefore uses about nine
times the optimal sample size. On the other side, a first stage larger
than \\n^\*\\ is itself the waste: the total can never be smaller than
\\m\\. The first stage should be as large as one can be confident does
not exceed \\n^\*\\.

The pilot size returned is \$\$m = \max\\m_0,\\ \lceil (A / w)^{1/(1 +
\gamma)} \rceil\\.\$\$ With the default `m_0 = 20`, a first stage of 20
keeps the expected total within about 1.5 times \\n^\*\\ whenever
\\n^\*\\ is at least 20 (1.54 at \\n^\* = 20\\, 1.47 at 50). The second
term raises \\A / w\\, which is \\n^\*\\ for \\\theta = 1\\, to the
power \\1/(1 + \gamma)\\, so that the first stage grows with the size of
the problem but more slowly than \\n^\*\\ does. This is the form of
pilot that Chattopadhyay and Kelley (2016, Equation 24) recommend,
following Mukhopadhyay and de Silva (2009), for a purely sequential
procedure. Using it here is this package's choice; as said above, the
source of this two stage procedure gives no rule. If the mean is
expected to be so small that \\n^\* = A \theta^2 / w\\ may be below 20,
lower `m_0`.

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

Zacks, S., & Mukhopadhyay, N. (2006). Bounded risk estimation of the
exponential parameter in a two-stage sampling. *Sequential Analysis*,
*25*(4), 437–452.
[doi:10.1080/07474940600934896](https://doi.org/10.1080/07474940600934896)
(As described in Mukhopadhyay and de Silva, 2009, Section 8.3.2,
Equation 8.3.13.)

## See also

[`ts_br_exp_mean`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md)

Other ts_br_exp_mean:
[`ts_br_exp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md)

## Author

Ken Kelley <kkelley@nd.edu>, Bhargab Chattopadhyay
<bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Mean time between service calls, squared error loss with A = 100,
# risk at most 4.
ts_br_exp_mean_pilot(A = 100, w = 4)
#> [1] 20

# A demanding bound: the rule exceeds the minimum.
ts_br_exp_mean_pilot(A = 100, w = 0.01)
#> [1] 100

# A smaller first stage, accepted knowingly (see the table in Details).
ts_br_exp_mean_pilot(A = 100, w = 4, m_0 = 10)
#> [1] 10
```
