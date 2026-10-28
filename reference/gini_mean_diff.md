# Gini's Mean Difference

Estimates Gini's mean difference, the average absolute difference
between two randomly chosen members of a population. It is a measure of
spread on the scale of the data (like the standard deviation, but based
on absolute rather than squared differences) and is the numerator of the
Gini index.

## Usage

``` r
gini_mean_diff(x, na.rm = TRUE)
```

## Arguments

- x:

  A numeric vector with at least two observations. Any real values are
  allowed.

- na.rm:

  Should missing values be removed before the calculation? Default
  `TRUE`. With `na.rm = FALSE` a missing value is an error.

## Value

A single nonnegative number, the estimated Gini mean difference, in the
units of `x`.

## Details

The population quantity is \\\Delta = E\|X_1 - X_2\|\\ for two
independent draws. The estimate returned is the U-statistic
\$\$\hat{\Delta}\_n = \binom{n}{2}^{-1} \sum\_{i \< j} \|x_i -
x_j\|,\$\$ the average of \\\|x_i - x_j\|\\ over all pairs of distinct
observations, which is unbiased for \\\Delta\\. It does not change when
a constant is added to every value. Kelley, Darku, and Chattopadhyay
(2018, Example 7) name Gini's mean difference among the parameters
covered by their general sequential procedure; Chattopadhyay and De
(2016) use it to estimate the Gini index, \\\hat{\Delta}\_n / (2
\bar{X}\_n)\\ (see
[`gini_index`](https://yelleknek.github.io/SMSD/reference/gini_index.md)).

## References

Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
pre-specified error bound. *Econometrics*, *4*(3), 30.
[doi:10.3390/econometrics4030030](https://doi.org/10.3390/econometrics4030030)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

## See also

[`gini_index`](https://yelleknek.github.io/SMSD/reference/gini_index.md),
[`seq_aipe_gini`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Annual salaries (in thousands of dollars) of ten employees of a firm:
# two employees chosen at random differ by this much on average.
salary <- c(28, 31, 35, 38, 42, 47, 55, 64, 90, 170)
gini_mean_diff(salary)
#> [1] 42.04444

# The Gini index is the mean difference over twice the mean.
gini_mean_diff(salary) / (2 * mean(salary))
#> [1] 0.3503704
gini_index(salary)
#> [1] 0.3503704
```
