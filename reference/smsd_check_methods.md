# Methods for Sequential Check Objects

Methods for the `smsd_check` objects returned by the sequential check
functions (the `seq_aipe_*` accuracy in parameter estimation family, the
`seq_mr_*` and `seq_br_*` minimum risk and bounded risk families, the
fixed width and two stage procedures, and the budget constrained
method).

## Usage

``` r
# S3 method for class 'smsd_check'
format(x, digits = getOption("digits"), ...)

# S3 method for class 'smsd_check'
print(x, digits = getOption("digits"), ...)

# S3 method for class 'smsd_check'
summary(object, ...)

# S3 method for class 'smsd_check'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'smsd_check'
tidy(x, ...)

# S3 method for class 'smsd_check'
glance(x, ...)
```

## Arguments

- x, object:

  An object of class `"smsd_check"`.

- digits:

  Number of significant digits to show for numeric values; by default
  the session's `getOption("digits")`. Stored values are never rounded;
  this affects display only.

- ...:

  Not used.

- row.names, optional:

  Ignored; present for compatibility with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector, one element per printed line.
[`summary()`](https://rdrr.io/r/base/summary.html) and `glance()` return
a one row `data.frame` with columns `method`, `estimate_term`,
`estimate`, `lower_limit`, `upper_limit`, `current_n`, `criterion_n`,
`is_satisfied`, `conf_level`, and `width` (`NA` where the procedure has
no such quantity). `tidy()` returns a one row `data.frame` with columns
`term`, `estimate`, `ci_lower`, `ci_upper`, and `conf_level`.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
`data.frame` with columns `term` and `value`.

## Details

[`print()`](https://rdrr.io/r/base/print.html) shows the method, the
current and criterion sample sizes, the target width and confidence
level when the procedure has them, every term, whether the stopping rule
is met, and any notes (for example that the variance estimate is not yet
usable, in which case the rule is reported as not met). For a check from
[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
at which sampling has ended because the budget cannot pay for another
needed observation, with a group still below its criterion, the stopping
line reads `FALSE (sampling has ended: the budget is exhausted)`,
because the rule itself was not met; the attribute `is_satisfied`, which
says whether sampling has ended, is `TRUE`.
[`format()`](https://rdrr.io/r/base/format.html) returns the same lines
as a character vector.
[`summary()`](https://rdrr.io/r/base/summary.html) and `glance()` return
a one row data frame describing the check; `tidy()` returns the point
estimate with its confidence limits in the column layout used by the
sibling package DMAR.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the plain `term`/`value` data frame without the class or the sequential
attributes. Subsetting with `[` returns a plain data frame (or vector)
whenever the result is no longer the complete check, so a partial object
is never printed as if it were a complete one.

**Width and coverage.** An AIPE check carries
`attr(x, "coverage_guarantee") = "asymptotic"`. Its display and
[`seq_report`](https://yelleknek.github.io/SMSD/reference/seq_report.md)
distinguish the asymptotic confidence level from satisfying the width
criterion. The latter does not establish nominal coverage at the
observed sample size. See
[`vignette("coverage-and-design", package = "SMSD")`](https://yelleknek.github.io/SMSD/articles/coverage-and-design.md)
for measured operating characteristics and simulation uncertainty.

**A check that cannot be evaluated.** In a small sample the estimate of
the variance that a stopping rule needs can be zero, negative, or not
finite. No variance has then been estimated: the stopping flag is
`FALSE`, `attr(x, "variance_unusable")` is `TRUE`, and a note says so.
This is a deliberate departure from the stopping rules as printed in the
source articles, under which such an estimate, or the lower bound
\\n^{-3}\\ that some of them put in its place, would make the rule be
met at once. A positive estimate, however small, is used as the articles
prescribe (with that bound where an article has one). When the point
estimate itself is not finite, the flag is `FALSE` as well and a note
gives the reason. (The two stage `ts_*` functions whose total sample
size the first stage variance or scale fixes stop with an error on a
first stage with no variation instead; see
[`SMSD`](https://yelleknek.github.io/SMSD/reference/SMSD-package.md).)
Test `attr(x, "is_satisfied")`, which is always `TRUE` or `FALSE`, to
drive a sampling loop. What the limit rows then hold depends on the
procedure. The `seq_aipe_*` functions report `NA` for both limits. The
`seq_fwd_*` functions still report their limits, which are computed from
the estimate and `d` alone and involve no variance.
[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
still reports an interval, computed with a variance of zero for each
group whose variance is unusable, which is too narrow. Neither of the
last two is an interval to report.
[`summary()`](https://rdrr.io/r/base/summary.html), `glance()`, and
`tidy()` return the limits as stored, without the notes or
`variance_unusable`, so test `attr(x, "variance_unusable")` before using
them.

**An object that is no longer one check.** Some operations leave the
class on something that is not a check as the function returned it:
[`rbind()`](https://rdrr.io/r/base/cbind.html) of two or more checks
(which keeps the attributes of the first), a data frame verb of another
package that binds checks, or an attribute removed by hand.
[`print()`](https://rdrr.io/r/base/print.html),
[`format()`](https://rdrr.io/r/base/format.html),
[`summary()`](https://rdrr.io/r/base/summary.html), `tidy()`,
`glance()`, and
[`seq_report`](https://yelleknek.github.io/SMSD/reference/seq_report.md)
stop with a message that says so, rather than describe the object as if
it were one complete check.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) always
works and gives the plain table. To keep a record of the stages of a
study, bind the one row
[`summary()`](https://rdrr.io/r/base/summary.html) of each check, or use
[`seq_path`](https://yelleknek.github.io/SMSD/reference/seq_path.md),
which records every stage.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
group_1 <- rnorm(40, mean = 0.5)
group_2 <- rnorm(40, mean = 0)
out <- seq_aipe_smd(group_1, group_2, width = 0.50)
out
#> Sequential check (standardized mean difference)
#> 
#>   Current n       : 40
#>   Criterion n     : 125
#>   Target width    : 0.5
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : 0.08665111
#>   smd               : 0.5257976
#>   upper_limit       : 0.964944
#>   xi_squared        : 2.008087
#>   criterion_n       : 125
#>   current_n         : 40
#> 
#>   Stopping rule met: FALSE
summary(out)
#>                         method estimate_term  estimate lower_limit upper_limit
#> 1 standardized mean difference           smd 0.5257976  0.08665111    0.964944
#>   current_n criterion_n is_satisfied conf_level width
#> 1        40         125        FALSE       0.95   0.5
generics::tidy(out)
#>   term  estimate   ci_lower ci_upper conf_level
#> 1  smd 0.5257976 0.08665111 0.964944       0.95
as.data.frame(out)
#>          term        value
#> 1 lower_limit   0.08665111
#> 2         smd   0.52579756
#> 3 upper_limit   0.96494402
#> 4  xi_squared   2.00808721
#> 5 criterion_n 125.00000000
#> 6   current_n  40.00000000
out[out$term == "smd", "value"]
#> [1] 0.5257976
```
