# Simple Linear Regression Slope

Calculates the ordinary least squares slope of the simple linear
regression of `y` on `x`: the expected change in `y` for a one unit
increase in `x`. This is the point estimate used by the sequential
procedure
[`seq_aipe_reg_coef`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md).

## Usage

``` r
slr_slope(x, y, na.rm = TRUE)
```

## Arguments

- x:

  A numeric vector, the predictor.

- y:

  A numeric vector, the outcome. It must have the same length as `x`;
  the two are paired by position.

- na.rm:

  Should a pair be removed when either of its values is missing? Default
  `TRUE`. With `na.rm = FALSE` a missing value is an error.

## Value

A single number, the estimated slope, in units of `y` per unit of `x`.
At least two complete pairs are needed, and `x` must vary: a predictor
with no variation is an error, because the slope is then undefined.
Values of `x` that agree to ten significant digits count as no variation
(such differences are usually rounding error); if they are real,
subtract a constant from `x` first. A slope too large for a double
(which takes units of `y` and `x` more than about 300 orders of
magnitude apart) is an error that says so.

## Details

The slope is \\b_1 = s\_{xy} / s_x^2\\, the sample covariance of `x` and
`y` divided by the sample variance of `x` (Kelley, Darku, and
Chattopadhyay, 2018, Equation 35). It is the same number as
`coef(lm(y ~ x))[2]`.

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

## See also

[`seq_aipe_reg_coef`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md);
in DMAR,
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.html)
for a fixed sample size confidence interval.

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# Hours of training and a performance rating for eight employees.
hours  <- c(2, 4, 5, 7, 8, 10, 12, 15)
rating <- c(48, 55, 53, 61, 66, 64, 73, 80)
slr_slope(hours, rating)
#> [1] 2.403056
```
