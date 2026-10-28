# Gini Index

Estimates the Gini index of inequality from a sample of nonnegative
values such as incomes, wealth, sales, or citation counts. The Gini
index is 0 when every unit has the same amount and approaches 1 when a
single unit holds everything.

## Usage

``` r
gini_index(x, na.rm = TRUE)
```

## Arguments

- x:

  A numeric vector of nonnegative values with at least two observations,
  not all zero.

- na.rm:

  Should missing values be removed before the calculation? Default
  `TRUE`. With `na.rm = FALSE` a missing value is an error.

## Value

A single number between 0 and 1, the estimated Gini index.

## Details

The population Gini index is \\G = \Delta / (2\mu)\\, where \\\Delta =
E\|X_1 - X_2\|\\ is Gini's mean difference (see
[`gini_mean_diff`](https://yelleknek.github.io/SMSD/reference/gini_mean_diff.md))
and \\\mu\\ is the mean. It is defined for nonnegative data with a
positive mean. The estimate returned is \$\$G_n =
\frac{\hat{\Delta}\_n}{2 \bar{X}\_n},\$\$ the sample Gini mean
difference (the average of \\\|x_i - x_j\|\\ over all pairs) divided by
twice the sample mean, which is the estimator used in the sequential
procedure of Chattopadhyay and De (2016) that
[`seq_aipe_gini`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md)
implements. Because the average is over pairs of distinct observations,
\\G_n\\ equals \\n / (n - 1)\\ times the "plug in" Gini index that some
software reports; the two agree as the sample grows. The estimate does
not change when every value is multiplied by the same positive constant.

## References

Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
pre-specified error bound. *Econometrics*, *4*(3), 30.
[doi:10.3390/econometrics4030030](https://doi.org/10.3390/econometrics4030030)

## See also

[`seq_aipe_gini`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md)
for the sequential confidence interval;
[`gini_mean_diff`](https://yelleknek.github.io/SMSD/reference/gini_mean_diff.md).

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Annual salaries (in thousands of dollars) of ten employees of a firm.
salary <- c(28, 31, 35, 38, 42, 47, 55, 64, 90, 170)
gini_index(salary)
#> [1] 0.3503704

# Equal salaries give 0; the index does not depend on the currency unit.
gini_index(rep(50, 10))
#> [1] 0
gini_index(1000 * salary)
#> [1] 0.3503704
```
