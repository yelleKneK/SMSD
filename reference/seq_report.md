# Report a Sequential Procedure in Words

Turns the result of a sequential procedure into a short paragraph that
can be pasted into a methods or results section: the procedure and its
source, its target (the confidence interval width and confidence level,
or the loss constants), the sample size at which the stopping rule was
met (or that it has not yet been met), and the estimate with its
confidence interval. The full reference for the procedure is returned
with the paragraph.

## Usage

``` r
seq_report(x, digits = 3L)

# S3 method for class 'smsd_report'
format(x, ...)

# S3 method for class 'smsd_report'
print(x, width = 0.9 * getOption("width"), ...)
```

## Arguments

- x:

  For `seq_report()`: an object of class `"smsd_check"`, as returned by
  a sequential check function (for example
  [`seq_aipe_smd`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md)
  or
  [`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)),
  or of class `"smsd_path"`, as returned by
  [`seq_path`](https://yelleknek.github.io/SMSD/reference/seq_path.md).
  For [`format()`](https://rdrr.io/r/base/format.html) and
  [`print()`](https://rdrr.io/r/base/print.html): the report that
  `seq_report()` returned.

- digits:

  The number of decimal places shown for the estimate and its confidence
  limits. A whole number from 0 to 15; default 3. This affects the text
  only: nothing stored in `x` is rounded. Constants the user supplied
  (the target width, the confidence level, the loss constants) are shown
  as supplied, to at most 7 significant digits.

- ...:

  Not used.

- width:

  The line width, in characters, at which
  [`print()`](https://rdrr.io/r/base/print.html) wraps the text (not a
  confidence interval width). Default: nine tenths of
  `getOption("width")`.

## Value

`seq_report()` returns an object of class
`c("smsd_report", "character")`: a single character string holding the
paragraph, with the attributes `reference` (a character vector with the
full reference for the procedure; of length zero when the method label
is not known) and `method` (the method label of `x`). Use
[`as.character()`](https://rdrr.io/r/base/character.html) for the bare
string.

[`print()`](https://rdrr.io/r/base/print.html) shows the paragraph,
wrapped to `width` characters, followed by the reference, and returns
its input invisibly. [`format()`](https://rdrr.io/r/base/format.html)
returns the bare paragraph as a plain character vector of the same
length as its input (one string), so a report can be held in a column of
a data frame.

## Details

Given a check object, the paragraph describes that one check: the rule
was met at the current sample size, or it has not yet been met. Whether
the current sample size is the *first* at which the rule was met is
something only the user's sampling record can say. Given a path from
[`seq_path`](https://yelleknek.github.io/SMSD/reference/seq_path.md),
the paragraph describes the stopping stage (the first sample size at
which the rule was met or, when it was never met, the last stage
checked) and adds where checking began and how often the rule was
checked.

For an AIPE result whose confidence level is asymptotic, the paragraph
states that satisfying the width criterion does not establish nominal
coverage at the observed sample size. The simulation vignette
[`vignette("coverage-and-design", package = "SMSD")`](https://yelleknek.github.io/SMSD/articles/coverage-and-design.md)
gives operating characteristics for difficult designs.

The paragraph is a starting point. It reports what the object holds and
cannot know how the data were collected, so read it against the study
before using it. Notes carried by the object are appended, each
introduced by "Note:". One such note says that the variance estimate was
zero, negative, or not finite at that sample size, so that no variance
had yet been estimated and the stopping rule was reported as not met
(see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)).
For an object flagged in that way (its attribute `variance_unusable` is
`TRUE`) the paragraph gives the estimate without the criterion sample
size and without a confidence interval, because neither rests on an
estimated variance; the object still holds what the check function
returned. For a path, the paragraph also says at how many of the stages
before the one it reports the variance estimate was unusable, because
the published stopping rules do not set such stages aside.

The procedure and its reference are looked up from the `method` label of
the object. The sources credited are Kelley, Darku, and Chattopadhyay
(2018) for the sequential accuracy in parameter estimation procedures
for the coefficient of variation, the standardized mean, the
standardized mean difference, the regression slope, and the contrast;
Chattopadhyay and De (2016) for the Gini index; Kelley, Darku, and
Chattopadhyay (2019) for the correlation coefficients and the squared
multiple correlation coefficient; Chattopadhyay and Kelley (2016, 2017)
for the minimum risk procedures for the coefficient of variation and the
standardized mean difference; Chattopadhyay, Bandyopadhyay, Kelley, and
Padalunkal (2025) for the budget constrained procedure; and Mukhopadhyay
and de Silva (2009) for the classical procedures. An object whose label
is not known to the package is still reported, without a citation and
without a reference. The constants of a risk based procedure are
reported when the object carries them, under the attributes `loss_A`
(the loss constant), `cost_c` or `cost_1` and `cost_2` (the sampling
cost), `risk_bound`, and, where the procedure has them, `k` (the power
of the estimation error in the loss) and `t` (the power of the sample
size in the cost). A width is reported as the target full width of the
confidence interval, except for a procedure that bounds the length of
its interval and does not aim at a width (the bounded length interval
for a median), where it is reported as the largest length allowed.

Sampling under the budget constrained procedure ends in one of three
ways, and the paragraph says which: every group reached its criterion
sample size, the budget could not pay for another observation from a
group still below its criterion, or the data in hand already cost more
than the budget. The second is reported as what it is, with the amount
of the budget left and the total of the criterion sample sizes that was
not reached (left out when the variance of a group could not be
estimated), and never as a stopping rule that was met; the third with
what the data cost and the budget. The paragraph reads these endings
from the attributes `rule_met`, `budget_exhausted`, and
`budget_exceeded` that
[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
sets, and the amounts from `budget` and `remaining_budget`; for a budget
constrained check without them it says only that sampling ended, which
is true of every ending. The note of
[`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
that says the budget is spent, or already exceeded, is not appended when
the paragraph has given the amounts, since it would say the same again.

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
accuracy in parameter estimation for population correlation
coefficients. *Psychological Methods*, *24*(4), 492–515.
[doi:10.1037/met0000203](https://doi.org/10.1037/met0000203)

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
difference with minimum risk: Maximizing accuracy and minimizing cost
with sequential estimation. *Psychological Methods*, *22*(1), 94–113.
[doi:10.1037/met0000089](https://doi.org/10.1037/met0000089)

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods*,
*30*(2), 425–439.
[doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
pre-specified error bound. *Econometrics*, *4*(3), 30.
[doi:10.3390/econometrics4030030](https://doi.org/10.3390/econometrics4030030)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`seq_path`](https://yelleknek.github.io/SMSD/reference/seq_path.md);
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
for the other methods of a check object.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
# Engagement after a training program against a control group.
trained <- rnorm(100, mean = 4.8, sd = 0.9)
control <- rnorm(100, mean = 4.5, sd = 0.9)

# One check, on the data collected so far.
out <- seq_aipe_smd(trained, control, width = 0.60)
seq_report(out)
#> Sample size was determined with the sequential accuracy in parameter
#> estimation procedure for the standardized mean difference (Kelley,
#> Darku, & Chattopadhyay, 2018), with a target full width of 0.6 for the
#> 95% confidence interval. The stopping rule was met at n = 100 per
#> group. At that sample size the estimate of the standardized mean
#> difference was 0.483, 95% CI [0.198, 0.768]. The confidence level is
#> asymptotic; meeting the width target does not establish nominal
#> coverage at this sample size.
#> 
#> Reference
#> Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#>     parameter estimation for a general class of effect sizes: A
#>     sequential approach. Psychological Methods, 23(2), 226–243.
#>     doi:10.1037/met0000127

# A path: the report describes the first sample size at which the rule
# was met, and says how often the rule was checked.
path <- seq_path(seq_aipe_smd, list(trained, control), width = 0.60,
                 by = 5)
report <- seq_report(path, digits = 2)
report
#> Sample size was determined with the sequential accuracy in parameter
#> estimation procedure for the standardized mean difference (Kelley,
#> Darku, & Chattopadhyay, 2018), with a target full width of 0.6 for the
#> 95% confidence interval. The stopping rule was evaluated from the pilot
#> sample size (n = 7 per group) onward, after every 5 additional
#> observations per group, and was first met at n = 92 per group. At that
#> sample size the estimate of the standardized mean difference was 0.52,
#> 95% CI [0.22, 0.82]. The confidence level is asymptotic; meeting the
#> width target does not establish nominal coverage at this sample size.
#> 
#> Reference
#> Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#>     parameter estimation for a general class of effect sizes: A
#>     sequential approach. Psychological Methods, 23(2), 226–243.
#>     doi:10.1037/met0000127
attr(report, "reference")
#> [1] "Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in parameter estimation for a general class of effect sizes: A sequential approach. Psychological Methods, 23(2), 226–243. doi:10.1037/met0000127"
```
