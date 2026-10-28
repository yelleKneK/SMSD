# Pilot Sample Size for Sequential AIPE: Squared Multiple Correlation Coefficient

Returns the number of observations to collect before the sequential
accuracy in parameter estimation (AIPE) stopping rule for the population
squared multiple correlation coefficient is evaluated for the first time
with
[`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md).

## Usage

``` r
seq_aipe_R2_pilot(
  width,
  p,
  conf_level = 0.95,
  ci_method = c("wald", "log_transformed"),
  m_0 = p + 2
)
```

## Arguments

- width:

  The target full width of the confidence interval for the population
  squared multiple correlation coefficient (\\\omega\\ in the article).
  A single number greater than 0 and less than 1.

- p:

  The number of predictors in the regression (\\k\\ in the article). A
  whole number of at least 1.

- conf_level:

  The confidence level (default 0.95). A single number strictly between
  0 and 1.

- ci_method:

  Which confidence interval the procedure is built on: `"wald"` (the
  default) or `"log_transformed"`. See
  [`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md).

- m_0:

  The smallest pilot sample size allowed. A whole number of at least
  `p + 2`, which is also the default. (The other pilot sample size
  functions in this package default to 4; here the smallest usable
  sample depends on the number of predictors.)

## Value

A length one integer: the pilot sample size.

## Details

Kelley, Darku, and Chattopadhyay (2019, p. 505, and Appendix C,
Equations 95 and 97) give a pilot sample size for each of the two
confidence intervals they study. With \\p\\ predictors (\\k\\ in the
article), target full width \\\omega\\, and standard normal quantile
\\z\_{1-\alpha/2}\\:

- `ci_method = "wald"`: \\\max\\p, 4 z\_{1-\alpha/2}/\omega\\\\, rounded
  up to a whole number;

- `ci_method = "log_transformed"`: \\p + 2\\.

The value returned is the larger of that size and `m_0`. The article
prints \\4 z\_{1-\alpha/2}/\omega\\ without saying how it is made a
whole number; rounding up is what reproduces the mean final sample sizes
of its Tables 13 and 14 in the conditions where most replications stop
at the pilot sample size.

With \\p + 1\\ or fewer observations the sample squared multiple
correlation is identically 1, so `m_0` is at least \\p + 2\\ and the
returned pilot sample size is never smaller than that. This departs from
the article's expression for the Wald interval only when \\\lceil 4
z\_{1-\alpha/2}/\omega \rceil \< p + 2\\, that is, for a wide target
with many predictors; it never happens in the conditions of the
article's tables.

For the log transformed interval the article's pilot sample size of
\\p + 2\\ lets the procedure stop within a few observations of the pilot
when the sample value happens to be near 0 or near 1, and the article's
Tables 11 and 12 show the loss of coverage that follows. With either
interval, a pilot sample that is not much larger than the number of
predictors gives a sample value that is far too high, and the procedure
can then stop at once with an interval that misses the population value.
Raising `m_0` guards against both.
[`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md)
attaches a note to an interval reported from fewer than \\p \max(20,
2p)\\ observations (20 per predictor, or \\2p\\ per predictor with more
than 10 predictors); an `m_0` of that size keeps every check out of that
range, at the price of a larger sample whenever the stopping rule would
have been met earlier. In simulations with 5, 10, and 20 predictors,
population values of .30 to .90, widths of .20 and .30, and a 95% level,
coverage with the Wald interval was then between .896 and .932, where
the default pilot sample size gave between .019 and .932. See the
measured coverage in the help for
[`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md),
and use
[`seq_aipe_R2_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_sensitivity.md)
to check a planned design.

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
accuracy in parameter estimation for population correlation
coefficients. *Psychological Methods*, *24*(4), 492–515.
[doi:10.1037/met0000203](https://doi.org/10.1037/met0000203)

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research*, *43*(4),
524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

## See also

[`seq_aipe_R2`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md),
[`seq_aipe_R2_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_sensitivity.md);
in DMAR, the fixed sample size planner
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.html).

Other seq_aipe_R2:
[`seq_aipe_R2()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md),
[`seq_aipe_R2_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Three predictors of job satisfaction; a 95% confidence interval for the
# population squared multiple correlation no wider than 0.20.
seq_aipe_R2_pilot(width = 0.20, p = 3)
#> [1] 40

# The log transformed interval starts from p + 2 unless m_0 raises it.
seq_aipe_R2_pilot(width = 0.20, p = 3, ci_method = "log_transformed")
#> [1] 5
seq_aipe_R2_pilot(width = 0.20, p = 3, ci_method = "log_transformed",
                  m_0 = 30)
#> [1] 30
```
