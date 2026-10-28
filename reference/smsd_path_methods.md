# Methods for Sequential Path Objects

Methods for the `smsd_path` objects returned by
[`seq_path`](https://yelleknek.github.io/SMSD/reference/seq_path.md):
the stage by stage record of a sequential procedure applied to data in
hand.

## Usage

``` r
# S3 method for class 'smsd_path'
format(x, digits = getOption("digits"), rows = 3L, ...)

# S3 method for class 'smsd_path'
print(x, digits = getOption("digits"), rows = 3L, ...)

# S3 method for class 'smsd_path'
summary(object, ...)

# S3 method for class 'smsd_path'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x, object:

  An object of class `"smsd_path"`.

- digits:

  Number of significant digits to show: a whole number from 1 to 22.
  Defaults to `getOption("digits")`. Stored values are never rounded;
  this affects display only.

- rows:

  The number of stages shown at the start and at the end of the printed
  table. A whole number of at least 1; default 3. A path with no more
  than `2 * rows + 1` stages is shown in full.

- ...:

  Not used.

- row.names, optional:

  Ignored; present for compatibility with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector, one element per printed line.
[`summary()`](https://rdrr.io/r/base/summary.html) returns a one row
`data.frame` with columns `method`, `estimate_term`, `stop_n` (`NA` when
the rule was never met), `n` (the stage described), `estimate`,
`lower_limit`, `upper_limit`, `current_n`, `criterion_n`,
`is_satisfied`, `variance_unusable`, `conf_level`, `width` (`NA` where
the procedure has no such quantity), and `stages` (the number of stages
checked). [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
returns a data frame (class `"data.frame"` only) with one row per stage
and the columns `n`, `estimate`, `lower_limit`, `upper_limit`,
`current_n`, `criterion_n`, `is_satisfied`, and `variance_unusable`.

## Details

[`print()`](https://rdrr.io/r/base/print.html) shows a compact header
(the method, the stages checked, the target, the first sample size at
which the stopping rule was met, the estimate with its interval there,
and any note the check attached at that stage), then the first and last
few stages, and closes with a reminder of which estimate is the one to
report. With `options(OutDec = ",")` the two limits are separated by a
semicolon. To keep the table narrow it leaves out columns that carry no
information for the path at hand: the confidence limits of a procedure
that has none, `current_n` when it equals `n` at every stage, and
`variance_unusable` when it is never `TRUE`. When it is, the header
counts the stages and a line under the table names their sample sizes,
because they need not be among the rows shown (stages too scattered to
list are counted, with the range they fall in), and says what the column
means: at those stages the variance estimate was zero, negative, or not
finite, so no variance had yet been estimated and the stopping rule was
reported as not met (see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)).
The header gives the estimate at such a stage without its interval. For
a procedure whose interval is bounded in length and not aimed at a width
(the bounded length interval for a median) the header gives the bound as
the largest length allowed.
[`format()`](https://rdrr.io/r/base/format.html) returns the same lines
as a character vector.
[`summary()`](https://rdrr.io/r/base/summary.html) returns the stopping
stage as a one row data frame: the first stage at which the rule was met
or, when it was never met, the last stage checked.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the plain table of stages without the class or the path attributes.
Subsetting with `[` returns a plain data frame (or vector) whenever the
result is no longer the complete path, so part of a path is never
printed as if it were the whole of one; a subset that keeps every stage
and every column is still the path.

Other operations can leave the class on something that is no longer a
path: [`rbind()`](https://rdrr.io/r/base/cbind.html) of two paths, a
data frame verb of another package that keeps the class on some of the
rows, or an attribute changed by hand.
[`print()`](https://rdrr.io/r/base/print.html),
[`format()`](https://rdrr.io/r/base/format.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html), and
[`seq_report`](https://yelleknek.github.io/SMSD/reference/seq_report.md)
stop with a message that says so, because anything they reported would
describe a path that was never computed.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) always
works and gives the plain table.

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

## See also

[`seq_path`](https://yelleknek.github.io/SMSD/reference/seq_path.md),
[`plot.smsd_path`](https://yelleknek.github.io/SMSD/reference/plot.smsd_path.md),
[`seq_report`](https://yelleknek.github.io/SMSD/reference/seq_report.md)

Other seq_path:
[`plot.smsd_path()`](https://yelleknek.github.io/SMSD/reference/plot.smsd_path.md),
[`seq_path()`](https://yelleknek.github.io/SMSD/reference/seq_path.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
# Weekly hours of remote work reported by 90 employees, in survey order.
hours <- rgamma(90, shape = 9, rate = 0.6)
path <- seq_path(seq_aipe_cv, hours, width = 0.15)
print(path, digits = 4, rows = 2)
#> Sequential path: coefficient of variation
#> 
#>   Stages checked  : 64 (n = 27 to 90)
#>   First check at  : the pilot sample size
#>   Target width    : 0.15
#>   Confidence level: 95%
#>   Rule first met  : n = 43
#>   Estimate there  : cv = 0.3161 [0.2579, 0.3743]
#> 
#>    n  estimate  lower_limit  upper_limit  criterion_n  is_satisfied
#>   27    0.3175       0.2485       0.3865           49         FALSE
#>   28    0.3162       0.2487       0.3837           48         FALSE
#>   ...
#>   89    0.3277       0.2862       0.3692           35          TRUE
#>   90    0.3269       0.2861       0.3678           35          TRUE
#> 
#>   The estimate to report is the one at n = 43. Later rows show what the
#>   rule would have said had sampling continued; they are not grounds for
#>   choosing a different stopping point.
summary(path)
#>                     method estimate_term stop_n  n  estimate lower_limit
#> 1 coefficient of variation            cv     43 43 0.3160575   0.2578572
#>   upper_limit current_n criterion_n is_satisfied variance_unusable conf_level
#> 1   0.3742578        43          42         TRUE             FALSE       0.95
#>   width stages
#> 1  0.15     64
head(as.data.frame(path))
#>    n  estimate lower_limit upper_limit current_n criterion_n is_satisfied
#> 1 27 0.3174939   0.2484683   0.3865195        27          49        FALSE
#> 2 28 0.3162329   0.2487346   0.3837312        28          48        FALSE
#> 3 29 0.3114342   0.2456160   0.3772523        29          46        FALSE
#> 4 30 0.3363939   0.2643931   0.4083947        30          51        FALSE
#> 5 31 0.3336129   0.2623672   0.4048586        31          51        FALSE
#> 6 32 0.3279930   0.2569223   0.3990638        32          51        FALSE
#>   variance_unusable
#> 1             FALSE
#> 2             FALSE
#> 3             FALSE
#> 4             FALSE
#> 5             FALSE
#> 6             FALSE
path[path$is_satisfied, c("n", "estimate")][1, ]
#>     n  estimate
#> 17 43 0.3160575
```
