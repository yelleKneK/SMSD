# Sequential AIPE Check for the Gini Index

Checks, with the data collected so far, whether enough observations are
in hand for a confidence interval for the population Gini index of
inequality (of incomes, wealth, sales, and the like) that is no wider
than a chosen width, without assuming a distribution and without a
planning value. This is the sequential fixed width interval procedure of
Chattopadhyay and De (2016), expressed in the package's accuracy in
parameter estimation (AIPE) form with the full width of the interval as
the target.

## Usage

``` r
seq_aipe_gini(x, width, conf_level = 0.95)
```

## Arguments

- x:

  A numeric vector of nonnegative values: the observations collected so
  far, at least 4 of them. Missing, infinite, or negative values are an
  error.

- width:

  The target full width of the confidence interval for the Gini index:
  the upper limit minus the lower limit (twice the half width \\d\\ of
  Chattopadhyay and De, 2016). A single positive number less than 1: the
  Gini index lies between 0 and 1, so an interval 1 or more wide would
  say nothing, and such a target is an error.

- conf_level:

  The confidence level of the interval, a single number strictly between
  0 and 1. Default 0.95.

## Value

An object of class `c("smsd_check", "data.frame")` (see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
for its [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and `tidy()`
methods) with one row per quantity and columns `term` and `value`. The
rows are, in order:

- `lower_limit`:

  the lower confidence limit (`NA` when the variance estimate is not yet
  usable or the estimate cannot be computed; see the section below);

- `gini`:

  the estimated Gini index (`NA` when every value is zero);

- `upper_limit`:

  the upper confidence limit (`NA` whenever `lower_limit` is);

- `xi_squared`:

  the estimated asymptotic variance \\\hat{\xi}\_n^2\\, at least
  \\n^{-3}\\ (the bound itself when the estimate is not yet usable);

- `criterion_n`:

  the sample size the stopping rule currently requires (at least
  `current_n + 1` when the variance estimate is not yet usable or the
  estimate cannot be computed);

- `current_n`:

  the current sample size.

The object also carries these attributes:

- `is_satisfied`:

  `TRUE` when the stopping rule is met, otherwise `FALSE`; never `NA`,
  and always identical to `current_n >= criterion_n`;

- `method`:

  the string `"Gini index"`;

- `current_n`, `criterion_n`:

  the values in the rows of those names;

- `conf_level`, `width`:

  the values supplied;

- `variance_unusable`:

  `TRUE` when the raw variance estimate is zero, negative, or not finite
  (see the section below), otherwise `FALSE`;

- `notes`:

  a character vector, empty unless something needs explaining.

## Details

The function evaluates the stopping rule once and returns. The user
collects the pilot sample (see
[`seq_aipe_gini_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini_pilot.md)),
calls `seq_aipe_gini()`, and, while the rule is not met, adds one or
more observations and calls it again. The function prints nothing and
keeps no state.

**Estimate.** The Gini index \\G = \Delta / (2\mu)\\, with \\\Delta =
E\|X_1 - X_2\|\\ Gini's mean difference, is estimated by \\G_n =
\hat{\Delta}\_n / (2 \bar{X}\_n)\\, as returned by
[`gini_index`](https://yelleknek.github.io/SMSD/reference/gini_index.md).
It is defined for nonnegative data with a positive mean, so a negative
value in `x` is an error. A sample in which every value is zero cannot
estimate the index; because that can happen by chance in a small early
sample, it is reported as "not yet estimable" rather than as an error
(see the section below).

**Stopping rule.** Sampling stops at the first \\n\\, at or after the
pilot size, with \$\$n \ge \left\lceil \frac{4
z\_{\alpha/2}^2}{\omega^2} \left(\hat{\xi}\_n^2 + \frac{1}{n}\right)
\right\rceil,\$\$ which is the stopping rule of Chattopadhyay and De
(2016, Equation 12), \\n \ge (z\_{\alpha/2}/d)^2 (V_n^2 + 1/n)\\, with
\\\omega = 2d\\. Here \\\hat{\xi}\_n^2\\ is their estimator \\V_n^2\\
(Equation 11) of the asymptotic variance of \\\sqrt{n}(G_n - G)\\
(Equation 6), built from the sample mean and variance, the Gini mean
difference, a second U-statistic, and a jackknife variance. The reported
interval is \\G_n \pm z\_{\alpha/2} \hat{\xi}\_n / \sqrt{n}\\. Whenever
the rule is met, the width of that interval is at most `width`. The
article bounds \\V_n^2\\ below by zero; the package bounds a positive
estimate below by \\n^{-3}\\, as for the other effect sizes, and treats
an estimate that is zero or negative as not yet usable (see the section
below).

**Which interval.** The interval reported is built as for every effect
size in the family, from the estimated standard error. It is not the
interval of Chattopadhyay and De (2016), which is the fixed width
interval \\(G_n - d, G_n + d)\\, always exactly \\\omega\\ wide, whose
coverage their stopping rule is designed to secure. The interval
reported is never wider than that one once the rule is met, and it
covers less often: in simulated studies with lognormal data, a Gini
index of 0.3, and a 95% level (10,000 studies each, the default `m_0`,
one observation added at each stage), the interval reported covered the
Gini index in 0.890 of the studies at `width = 0.2` and 0.916 at
`width = 0.1`, and the article's interval in 0.986 and 0.960. For the
article's interval, take the estimate (the row `gini`) plus and minus
half of `width`.

Kelley, Darku, and Chattopadhyay (2018) cite this procedure and give the
general sequential AIPE framework that the rest of the `seq_aipe_*`
family implements; their article does not itself contain a Gini index
procedure.

**What the guarantee is.** The width guarantee is exact. The coverage
guarantee is asymptotic: coverage approaches `conf_level` as `width`
shrinks. In simulated studies with a 95% interval and one observation
added at each stage, coverage for exponential data (a Gini index of 0.5)
was 0.951 at `width = 0.1` (3,000 studies, about 142 observations at the
end), 0.937 at `width = 0.2` (2,000 studies, about 44 observations), and
0.945 at `width = 0.4` (6,000 studies, about 17 observations). For
lognormal data with a Gini index of 0.3 it was 0.890 at `width = 0.2`
(10,000 studies, about 30 observations) and 0.870 at `width = 0.4`
(about 13 observations). The variance estimate was never zero, negative,
or below its lower bound in these studies; the shortfall comes from the
small final samples that wide targets allow. Use
[`seq_aipe_gini_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
to study a planned study before collecting data.

## Estimates That Are Not Yet Usable

When the raw variance estimate is zero (exactly, or to within rounding
error), negative, or not finite, which is common at small sample sizes
for some effect sizes, or when the point estimate cannot be computed (a
sample with no variation, where values that agree to ten significant
digits count as equal), the check reports the stopping rule as not met,
sets `attr(out, "variance_unusable")` to `TRUE`, and attaches a note
that [`print()`](https://rdrr.io/r/base/print.html) shows; see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md).
This is a deliberate departure from the stopping rule as printed in the
source articles, under which the lower bound \\n^{-3}\\ would stand in
for such an estimate and the rule would be met at once with an interval
far narrower than requested; a positive estimate, however small, is used
as the articles prescribe (the larger of it and \\n^{-3}\\).

## References

Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
pre-specified error bound. *Econometrics*, *4*(3), 30.
[doi:10.3390/econometrics4030030](https://doi.org/10.3390/econometrics4030030)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

## See also

[`seq_aipe_gini_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini_pilot.md)
for the pilot size,
[`seq_aipe_gini_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
to simulate the procedure before collecting data, and
[`gini_index`](https://yelleknek.github.io/SMSD/reference/gini_index.md)
and
[`gini_mean_diff`](https://yelleknek.github.io/SMSD/reference/gini_mean_diff.md)
for the point estimates. (DMAR has no fixed sample size analog for the
Gini index.)

Other seq_aipe_gini:
[`seq_aipe_gini_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini_pilot.md),
[`seq_aipe_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# How unequal are household incomes (in thousands of dollars) in a
# region? The aim is a 95% confidence interval for the Gini index that
# is no wider than 0.15.
set.seed(113)
m <- seq_aipe_gini_pilot(width = 0.15)
income <- rlnorm(m, meanlog = 4, sdlog = 0.6)
out <- seq_aipe_gini(income, width = 0.15)
out
#> Sequential check (Gini index)
#> 
#>   Current n       : 27
#>   Criterion n     : 47
#>   Target width    : 0.15
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : 0.2540293
#>   gini              : 0.3200881
#>   upper_limit       : 0.386147
#>   xi_squared        : 0.03067111
#>   criterion_n       : 47
#>   current_n         : 27
#> 
#>   Stopping rule met: FALSE

# The rule is not met at the pilot. Add one household at a time and check
# again, until it is.
while (!attr(out, "is_satisfied")) {
  income <- c(income, rlnorm(1, meanlog = 4, sdlog = 0.6))
  out <- seq_aipe_gini(income, width = 0.15)
}
out
#> Sequential check (Gini index)
#> 
#>   Current n       : 47
#>   Criterion n     : 47
#>   Target width    : 0.15
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : 0.2835468
#>   gini              : 0.3452086
#>   upper_limit       : 0.4068704
#>   xi_squared        : 0.04651941
#>   criterion_n       : 47
#>   current_n         : 47
#> 
#>   Stopping rule met: TRUE
```
