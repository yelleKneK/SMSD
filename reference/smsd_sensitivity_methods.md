# Print a Sensitivity Summary

The `*_sensitivity()` functions,
[`seq_aipe_smd_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md)
among them, return a summary that these methods print. The object is a
data frame with columns `term` and `value` and can be used as one; these
methods only change how it is displayed. Each value is formatted on its
own, in fixed notation, so that a table mixing a replication count, a
coverage proportion, and a mean sample size does not fall into
scientific notation as the default data frame print would have it.

## Usage

``` r
# S3 method for class 'smsd_sensitivity'
format(x, digits = getOption("digits"), ...)

# S3 method for class 'smsd_sensitivity'
print(x, digits = getOption("digits"), ...)
```

## Arguments

- x:

  An object of class `"smsd_sensitivity"`.

- digits:

  Number of significant digits to show. Defaults to
  `getOption("digits")`. Stored values are never rounded; this affects
  display only.

- ...:

  Not used.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector, one element per printed line.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# G is kept small so the example runs quickly; a study to be reported
# deserves G = 10000 or more.
out <- seq_aipe_smd_sensitivity(true_value = 0.5, width = 0.8, G = 20,
                                seed = 113)
out
#> term                   value
#> true_value             0.5
#> target_width           0.8
#> conf_level             0.95
#> replications           20
#> mean_n                 50.6
#> median_n               49.5
#> sd_n                   4.429566
#> mean_width             0.7906421
#> median_width           0.7906787
#> pct_width_within       100
#> coverage               1
#> pct_hit_max_n          0
#> pct_variance_unusable  0
out$value[out$term == "coverage"]
#> [1] 1
```
