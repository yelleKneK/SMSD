# Stopping Path of a Sequential Procedure for Data in Hand

Applies a sequential check function to data already collected, at each
successive sample size, and returns the whole path in one object: the
estimate, its confidence limits, the criterion sample size, and the
stopping decision at every stage, together with the first sample size at
which the stopping rule is met. One call replaces the collect and check
loop that would otherwise be written by hand, and the
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method draws
the picture a preregistration or a methods section needs.

## Usage

``` r
seq_path(check, data, ..., start = NULL, by = 1L)
```

## Arguments

- check:

  A sequential check function: any function that returns a sequential
  check object (class `"smsd_check"`), such as
  [`seq_aipe_cv`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md),
  [`seq_aipe_smd`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md),
  [`seq_aipe_r`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md),
  or
  [`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md).

- data:

  The data collected so far, in the order collected: a vector, a list
  with one element for each data argument of `check`, or a data frame
  whose columns are named after those arguments. See *Supplying the
  Data*.

- ...:

  The other arguments of `check`, passed unchanged and whole at every
  stage (for example `width` and `conf_level`, or `A` and `c`). Data do
  not belong here: an argument that looks like data is refused unless it
  is wrapped in [`I()`](https://rdrr.io/r/base/AsIs.html) (see
  *Supplying the Data*). `seq_path()` matches its own arguments `check`
  and `data` by exact name or by position only, so an argument such as
  `c = 75` or `d = 0.5` reaches the check function and is not taken, by
  partial matching, for `check` or `data`.

- start:

  The first sample size checked. `NULL` (the default) is the procedure's
  pilot sample size when `check` is one of this package's check
  functions, and the smallest sample size `check` accepts otherwise.
  Other choices: a single whole number; a pilot function (called with
  those arguments in `...` that it shares with `check`); `"pilot"` (the
  package's pilot function for `check`; an error if there is none); or
  `"earliest"` (the smallest sample size `check` accepts). See *Where
  the Path Starts*.

- by:

  The number of observations (per data element) added between checks.
  The default, 1, is the purely sequential procedure of the cited
  articles. A larger value shows what a study that checks after every
  batch of `by` observations would have seen; observations beyond the
  last complete batch are not used.

## Value

An object of class `c("smsd_path", "data.frame")` with one row per stage
and the columns

- `n`:

  the number of observations taken from each data element at this stage
  (integer).

- `estimate`:

  the point estimate.

- `lower_limit`, `upper_limit`:

  the confidence limits (`NA` for the risk based procedures, which have
  no interval).

- `current_n`:

  the sample size as `check` reports it.

- `criterion_n`:

  the criterion sample size at this stage (`NA` for a procedure that has
  none).

- `is_satisfied`:

  whether the stopping rule is met at this stage.

- `variance_unusable`:

  whether the variance estimate was zero, negative, or not finite at
  this stage. No variance has then been estimated, and the check
  functions report the stopping rule as not met, whatever the other
  columns show (see
  [`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)):
  the criterion sample size and any confidence limits in that row rest
  on no estimated variance. `NA` when `check` does not report it.

The attributes are `stop_n` (the first `n` at which the rule is met;
`NA` if it is never met), `method`, `estimate_term` (the name of the
quantity estimated), `conf_level` and `width` (absent for procedures
without them), `start` (the first sample size checked), `start_rule`
(how it was chosen: `"pilot"`, `"earliest"`, or `"supplied"`), `by`,
`n_observations` (the number of observations in hand per data element),
`arguments` (the list passed through `...`, as `check` received it), and
`final_check`: the `smsd_check` object at `stop_n` or, when the rule is
never met, at the last stage. Stored values are never rounded. A warning
that `check` gives at several stages is reported once, with the number
of stages.

## What the Path Shows, and What It Does Not

The path shows what the stopping rule *would have said* at each sample
size for these data in this order. A sequential design decides at each
stage from the observations collected so far, without looking ahead: it
stops at the first sample size at which the rule is met, and the
estimate and interval it reports are the ones at that sample size
(`stop_n`), not at a later, more favorable one. The rows after `stop_n`
are shown for illustration only. Choosing a later stage because its
estimate or interval looks better is not the sequential procedure, and
the properties the procedure was derived to have (the width and coverage
of the interval, or the risk) do not hold for an estimate chosen that
way. For the same reason the observations must be in the order in which
they were, or will be, collected.

## Supplying the Data

At each stage `data` is cut to its first \\n\\ observations and passed
to `check`.

- A vector is passed as the first argument of `check` (the one sample
  procedures).

- A list holds one element for each data argument of `check`: either
  named after those arguments, as in `list(group_1 = a, group_2 = b)`
  for two independent groups or `list(x = x, y = y)` for a predictor and
  a response, or unnamed and in the order `check` expects them (the help
  page of each check function gives the names of its data arguments).
  Each element is a vector, a matrix or data frame whose rows are
  observations, or a list of these (for example the `groups` of
  [`seq_aipe_c`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md)),
  every member of which is cut to its first \\n\\ observations.

- A data frame is read as a named list of its columns, so the column
  names must be argument names of `check`.

Every element must hold the same number of observations, because stage
\\n\\ uses the first \\n\\ of each: \\n\\ is the sample size for one
sample procedures, the per group sample size for group designs, and the
number of pairs for paired data. The column `current_n` gives the sample
size on the scale `check` itself reports (for example the total over
both groups for
[`seq_mr_smd`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)),
which is the scale of `criterion_n`.

*All* of the data go in `data`. Only `data` is cut to its first \\n\\
observations; everything in `...` is passed whole at every stage. A
second group or a second variable written after `data` (for example the
second group in `seq_path(seq_mr_smd, group_1, group_2, ...)`) would
therefore be seen in full from the first stage on, which is the looking
ahead a sequential design rules out. `seq_path()` stops with a message
when an argument in `...` looks like data: it holds two or more
observations and `check` still runs when it is cut along with `data`. In
the rare case that such an argument really is not data, wrap it in
[`I()`](https://rdrr.io/r/base/AsIs.html) to pass it whole.

Missing values are refused. The path is indexed by the number of
observations collected, and several check functions drop missing values
by default (the minimum risk, bounded risk, fixed width, and bounded
length functions, whose `na.rm` is `TRUE` unless set otherwise; the
`seq_aipe_*` functions stop on a missing value), so with missing values
the sample size a stage reports would be smaller than the stage and the
first check would come before the pilot sample size had been reached.
Remove the incomplete observations first (both members of a pair, the
whole row of a matrix or data frame), and trim the elements to a common
length if they then differ.

## Where the Path Starts

A sequential procedure begins checking at its pilot sample size, not
before, and the path does the same. By default (`start = NULL`) the
first stage is the procedure's pilot sample size whenever `check` is one
of this package's check functions: its `_pilot` companion (for example
[`seq_aipe_cv_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv_pilot.md)
for
[`seq_aipe_cv`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md))
is called with those arguments in `...` that it shares with `check`, and
with its own defaults for the rest. For any other function (a wrapper,
or a check function of the user's own) the default is the smallest
sample size the function accepts, found by trying. `start` can also be
given as a number, as a pilot function, or as `"earliest"` to see the
rule from the smallest sample size the check accepts. A pilot function
that needs an argument `check` does not take (the minimum pilot sample
size `m_0`, for example, or a number of predictors) cannot be given that
argument through `...`: call the pilot function yourself and pass the
number it returns as `start`.

Where the path starts can change where it stops. The stopping rules of
[`seq_aipe_cv`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md),
[`seq_aipe_sm`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm.md),
[`seq_aipe_smd`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md),
[`seq_aipe_r`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md),
[`seq_aipe_reg_coef`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md),
[`seq_aipe_c`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md),
[`seq_aipe_gini`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md),
[`seq_mr_cv`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md),
and
[`seq_mr_smd`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)
(with one sampling cost) contain a correction term that keeps the rule
from being met below the sample size the pilot formula gives. The pilot
sample size is the larger of that and a minimum (4 by default), so for
these procedures an earlier start can change `stop_n` only when the
target is so loose that the formula gives less than the minimum, and
then only if the check function accepts fewer observations than the
minimum. Several of the classical procedures (for example the fixed
width, minimum risk, and bounded risk procedures for a normal mean) have
no correction term at all: started below the pilot, their rule can be
met by chance at a very small sample size, which is not a sample size at
which the procedure itself would have stopped. If a study used a pilot
sample size other than the default (a larger minimum, for example),
supply that number.

## Procedures That Do Not Fit

The two stage procedures (`ts_*`) take one pilot sample and one further
sample whose size the pilot determines, so they have no stage by stage
path: call the `ts_*` function once on the pilot sample to learn the
total sample size required, and once more when the additional
observations are in. The budget constrained procedure
([`seq_budget`](https://yelleknek.github.io/SMSD/reference/seq_budget.md))
allocates observations unequally across groups according to their costs
and variances, so its path is not indexed by one common sample size:
call
[`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
after each new observation. The minimum risk procedure for the
standardized mean difference with two costs
([`seq_mr_smd`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)
with `c = c(c1, c2)`) is refused for the same reason: its pilot returns
one size per group. `seq_path()` stops with a message saying so when
given one of these.

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

## See also

[`plot.smsd_path`](https://yelleknek.github.io/SMSD/reference/plot.smsd_path.md)
for the plot,
[`smsd_path_methods`](https://yelleknek.github.io/SMSD/reference/smsd_path_methods.md)
for the other methods, and
[`seq_report`](https://yelleknek.github.io/SMSD/reference/seq_report.md)
for a sentence that reports the result.

Other seq_path:
[`plot.smsd_path()`](https://yelleknek.github.io/SMSD/reference/plot.smsd_path.md),
[`smsd_path_methods`](https://yelleknek.github.io/SMSD/reference/smsd_path_methods.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
# Response times (seconds) of 120 participants, in the order tested. At
# what sample size does a 95% confidence interval of full width 0.10 for
# the coefficient of variation first meet the stopping rule?
response_time <- rgamma(120, shape = 16, rate = 20)
path <- seq_path(seq_aipe_cv, response_time, width = 0.10)
path
#> Sequential path: coefficient of variation
#> 
#>   Stages checked  : 81 (n = 40 to 120)
#>   First check at  : the pilot sample size
#>   Target width    : 0.1
#>   Confidence level: 95%
#>   Rule first met  : n = 61
#>   Estimate there  : cv = 0.2530166 [0.2152415, 0.2907917]
#> 
#>     n   estimate  lower_limit  upper_limit  criterion_n  is_satisfied
#>    40  0.2322022    0.1843518    0.2800525           76         FALSE
#>    41  0.2332873    0.1859133    0.2806612           75         FALSE
#>    42  0.2326324    0.1873538    0.2779110           72         FALSE
#>   ...
#>   118  0.2503885    0.2223552    0.2784219           51          TRUE
#>   119  0.2504532    0.2228470    0.2780593           50          TRUE
#>   120  0.2493359    0.2217831    0.2768887           50          TRUE
#> 
#>   The estimate to report is the one at n = 61. Later rows show what the
#>   rule would have said had sampling continued; they are not grounds for
#>   choosing a different stopping point.
attr(path, "stop_n")
#> [1] 61

# Engagement after a training program against a control group, checked
# after every 10 participants per group. This study's pilot was 20 per
# group, larger than the default pilot, so it is supplied.
trained <- rnorm(120, mean = 4.8, sd = 0.9)
control <- rnorm(120, mean = 4.5, sd = 0.9)
seq_path(seq_aipe_smd, list(group_1 = trained, group_2 = control),
         width = 0.60, start = 20, by = 10)
#> Sequential path: standardized mean difference
#> 
#>   Stages checked  : 11 (n = 20 to 120 by 10)
#>   First check at  : the sample size supplied as 'start'
#>   Target width    : 0.6
#>   Confidence level: 95%
#>   Rule first met  : n = 100
#>   Estimate there  : smd = 0.4527979 [0.1591645, 0.7464314]
#> 
#>     n   estimate  lower_limit  upper_limit  criterion_n  is_satisfied
#>    20  0.2513720  -0.37484398    0.8775880           90         FALSE
#>    30  0.5121502  -0.02351928    1.0478197           98         FALSE
#>    40  0.4407556  -0.01942790    0.9009391           96         FALSE
#>   ...
#>   100  0.4527979   0.15916448    0.7464314           97          TRUE
#>   110  0.4843968   0.20437793    0.7644156           97          TRUE
#>   120  0.4383967   0.17313872    0.7036547           95          TRUE
#> 
#>   The estimate to report is the one at n = 100. Later rows show what
#>   the rule would have said had sampling continued; they are not grounds
#>   for choosing a different stopping point.

# A minimum risk procedure: the loss constant A and the unit cost c
# replace the width and the confidence level.
seq_path(seq_mr_cv, response_time, A = 400000, c = 75)
#> Sequential path: minimum risk point estimation (coefficient of variation)
#> 
#>   Stages checked  : 103 (n = 18 to 120)
#>   First check at  : the pilot sample size
#>   Arguments       : A = 400000, c = 75
#>   Rule first met  : n = 20
#>   Estimate there  : cv = 0.2195011
#> 
#>     n   estimate  criterion_n  is_satisfied
#>    18  0.2274672           20         FALSE
#>    19  0.2230411           20         FALSE
#>    20  0.2195011           20          TRUE
#>   ...
#>   118  0.2503885           14          TRUE
#>   119  0.2504532           14          TRUE
#>   120  0.2493359           14          TRUE
#> 
#>   The estimate to report is the one at n = 20. Later rows show what the
#>   rule would have said had sampling continued; they are not grounds for
#>   choosing a different stopping point.
```
