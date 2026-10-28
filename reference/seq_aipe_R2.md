# Sequential AIPE Check for the Squared Multiple Correlation Coefficient

Given the predictors and the response collected so far, evaluates
whether the sequential accuracy in parameter estimation (AIPE) stopping
rule of Kelley, Darku, and Chattopadhyay (2019) has been satisfied for
the population squared multiple correlation coefficient, and returns the
current estimate with its confidence interval. Use it when the goal is a
confidence interval for the population squared multiple correlation that
is no wider than a chosen width and the population value needed to plan
a fixed sample size is unknown.

## Usage

``` r
seq_aipe_R2(
  x,
  y,
  width,
  conf_level = 0.95,
  ci_method = c("wald", "log_transformed")
)
```

## Arguments

- x:

  The predictors collected so far: a numeric matrix or a data frame of
  numeric columns with one row per observation and one column per
  predictor, or a numeric vector for a single predictor. Do not include
  a column of ones; the intercept is part of the model.

- y:

  Numeric vector. The response collected so far, one value per row of
  `x`.

- width:

  The target full width of the confidence interval for the population
  squared multiple correlation coefficient (\\\omega\\ in the article).
  A single number greater than 0 and less than 1.

- conf_level:

  The confidence level (default 0.95). A single number strictly between
  0 and 1.

- ci_method:

  Which confidence interval to report and to build the stopping rule on:
  `"wald"` (the default; the symmetric interval with the Olkin and Finn,
  1995, variance) or `"log_transformed"` (the Bonett and Wright, 2011,
  interval).

## Value

An object of class `c("smsd_check", "data.frame")` with columns `term`
and `value` and six rows:

- `lower_limit`:

  The lower confidence limit.

- `R2`:

  The sample squared multiple correlation coefficient.

- `upper_limit`:

  The upper confidence limit.

- `xi_squared`:

  The estimated asymptotic variance of \\R^2\\ under multivariate
  normality, \\4 R^2 (1 - R^2)^2\\ (Olkin and Finn, 1995). The Wald
  standard error is its square root over \\\sqrt{n}\\. It is reported as
  computed; no lower bound is applied to it.

- `criterion_n`:

  The criterion sample size of the stopping rule, given the current
  data. The rule is met when `current_n` is at least this number. When
  no decision can be made it is at least `current_n + 1`.

- `current_n`:

  The number of observations supplied.

A quantity that is not defined for the current data is `NA`: the two
limits whenever no decision can be made (the variance estimate is
unusable, the estimate cannot be computed, or
`ci_method = "log_transformed"` and \\n = p + 2\\), and `R2` and
`xi_squared` when the response does not vary. The attributes are:

- `method`:

  The label of the method, which names the interval:
  `"squared multiple correlation, Wald interval"` or
  `"squared multiple correlation, log transformed interval"`.

- `is_satisfied`:

  `TRUE` when the stopping rule is met, otherwise `FALSE`; never `NA`.
  It is `TRUE` exactly when `current_n` is at least `criterion_n`.

- `current_n`, `criterion_n`:

  The two sample sizes, as integers.

- `conf_level`, `width`:

  As supplied.

- `notes`:

  A character vector, of length zero when there is nothing to say, that
  [`print()`](https://rdrr.io/r/base/print.html) shows under the
  stopping line. The notes that say why no decision can be made (the
  estimate cannot be computed; the variance estimate is unusable; the
  log transformed interval is not yet defined) accompany a rule that is
  not met. The note that the sample is small relative to the number of
  predictors can accompany either decision and does not change it.

- `ci_method`:

  The interval used, `"wald"` or `"log_transformed"`.

- `p`:

  The number of predictors.

- `variance_unusable`:

  `TRUE` when the variance estimate is unusable (zero, or not finite
  because the estimate cannot be computed), otherwise `FALSE`.

## Details

**This procedure assumes multivariate normality.** Unlike the
correlation procedures of the same article
([`seq_aipe_r`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md))
and the other `seq_aipe_*` functions, it is not distribution free: the
article develops it for multivariate normal data only (p. 503), because
no distribution free interval for the squared multiple correlation was
available. See the section on assumptions below.

**The coverage guarantee is asymptotic and can fail badly in small
samples, above all with many predictors.** Coverage approaches the
confidence level as the target width narrows and the sample grows for
interior population values, \\0 \< \mathrm{P}^2 \< 1\\, under the stated
assumptions. At the null boundary \\\mathrm{P}^2 = 0\\, the normal
approximation is degenerate; increasing the sample size does not restore
the usual Wald calibration. See the section on the null boundary. At a
wide target it falls short even with few predictors: with five
predictors, a 95% level, and a target width of .30 it was .851 for a
population value of .50 and .797 for a population value of .90. It falls
much further when the number of predictors is not small relative to the
sample size that the width calls for. The sample \\R^2\\ overestimates
the population value by roughly \\(1 - \mathrm{P}^2) p / n\\, both
intervals are built around it, and the stopping rule does not correct
for it, so the rule is met with an interval that tends to lie above the
population value. In simulations of the full procedure with multivariate
normal data and a 95% level (10,000 replications per condition),
coverage with the default interval for a population value of .50 was
.865 with 10 predictors and .790 with 15 at a target width of .20, and
.668 with 10 and .361 with 15 at a target width of .30. With 20
predictors unrelated to the response (a population value of 0) and a
target width of .30, the final interval contained 0 in 4% of
replications; in 36% of data sets the rule was met at the very first
check, with 27 observations, a sample value near .88, and a lower limit
near .80. The article offers no remedy, and the function does not alter
the stopping rule to correct for it. What it does is say so: whenever an
interval is reported from a sample that is small relative to the number
of predictors (fewer than \\p \max(20, 2p)\\ observations), the object
returned carries a note (see the section on small samples below). The
article's own simulations cover 2, 5, and 10 predictors at widths of .05
and .10. Before relying on the procedure, run
[`seq_aipe_R2_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_sensitivity.md)
at the planned width and number of predictors; a larger pilot sample
(`m_0` in
[`seq_aipe_R2_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_pilot.md))
helps. See the sections on the asymptotic guarantee and on small samples
below.

The function is stateless. After collecting the pilot sample, whose size
[`seq_aipe_R2_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_pilot.md)
gives, the user calls `seq_aipe_R2()` and, while the rule is not met,
adds observations and calls it again on all the data collected so far.

[`seq_path`](https://yelleknek.github.io/SMSD/reference/seq_path.md)
applies the function to data already in hand at every sample size from
the pilot sample size on. For most procedures it works out the pilot
sample size itself. For this one it cannot, because the pilot sample
size depends on the number of predictors, which the check function reads
from the data and
[`seq_path()`](https://yelleknek.github.io/SMSD/reference/seq_path.md)
has no way to hand to the pilot function; called without `start` it
stops with a message. Supply the start, as in the last example:

    seq_path(seq_aipe_R2, list(x = x, y = y), width = width,
             start = seq_aipe_R2_pilot(width, p = ncol(x)))

The printed path gives the stopping point and the interval there but not
the notes of the check at that point, so a path that stops while the
sample is small relative to the number of predictors (see the section on
small samples below) prints without that note. The check at the stopping
point, with its notes, is `attr(path, "final_check")`; look at it before
the interval is reported.

## Estimator and Confidence Intervals

With \\n\\ observations on \\p\\ predictors (\\k\\ in the article) and a
response, the estimate is the sample squared multiple correlation
coefficient \\R^2\\ (Equation 57 of the article), the value
`summary(lm(y ~ x))$r.squared` reports. Write \\R\\ for its positive
square root, \\\omega\\ for `width`, and \\z = z\_{1-\alpha/2}\\ for the
standard normal quantile.

The article studies two intervals, both with \\R\\ in place of the
population value in the standard error:

- `ci_method = "wald"` (Equation 59), built on the asymptotic variance
  \\4 \mathrm{P}^2 (1 - \mathrm{P}^2)^2 / n\\ of Olkin and Finn
  (1995, p. 161, Equation 10): \$\$R^2 \pm z \\ \frac{2 R (1 -
  R^2)}{\sqrt{n}}.\$\$

- `ci_method = "log_transformed"` (Equation 58), the interval of Bonett
  and Wright (2011), formed for \\\ln(1 - R^2)\\ and transformed back:
  \$\$1 - \exp\left\\\ln(1 - R^2) \pm z \\ \frac{2 R}{\sqrt{n - p -
  2}}\right\\.\$\$

Neither interval is truncated to \\\[0, 1\]\\ in the article, and the
limits are returned as computed: with either interval the lower limit
can be negative, and with the Wald interval the upper limit can exceed 1
when the sample is very small (\\n \< 4 z^2 R^2\\). Once the stopping
rule is met the interval is no wider than \\\omega\\. The log
transformed interval is not defined at \\n = p + 2\\, where its limits
are returned as `NA`.

The Wald interval is the default because the article concludes from its
simulations (p. 506) that the procedure built on it did better than the
one built on the log transformed interval in terms of the sample size
required; its coverage is also closer to the confidence level in 59 of
the 60 conditions of the article's Tables 11 to 14 (see the section on
the asymptotic guarantee).

## Stopping Rule

The rule is met when \\n\\ is at least `criterion_n`, the smallest whole
number not below

- for `"wald"`: \$\$\frac{16 z^2}{\omega^2} \left(R^2 +
  \frac{1}{n}\right) \left(1 - R^2 + \frac{1}{n}\right)^2;\$\$

- for `"log_transformed"` (Equation 62): \$\$p + 2 + \frac{4 R^2
  z^2}{L^2}, \quad L = \ln\left\[\frac{1}{2}\left(\frac{\omega}{1 -
  R^2} + \sqrt{\frac{\omega^2}{(1 - R^2)^2} + 4}\right)\right\].\$\$

These estimate the sample sizes that would be needed if the population
value \\\mathrm{P}^2\\ were known (Equations 61 and 60): \\16 z^2
\mathrm{P}^2 (1 - \mathrm{P}^2)^2 / \omega^2\\ for the Wald interval,
and the log transformed expression with \\\mathrm{P}^2\\ in place of
\\R^2\\. The terms \\1/n\\ in the Wald rule are the correction that
keeps sampling from stopping too early; the log transformed rule has no
such correction in the article.

**The Wald rule differs from Equation 63 as printed.** The printed
equation has \\(1 - (R^2 + 1/n))^2\\ as its last factor, which *lowers*
the criterion whenever \\R^2\\ is above about one third, allows the
interval at the stopping point to be wider than \\\omega\\, and does not
reproduce the article's Tables 13 and 14. In simulations of all 60
conditions of those tables the printed form gave a mean final sample
size below the tabled value in every condition, by more than three
standard errors in 49 of them (and by as much as 39.7), and a mean
interval width above \\\omega\\ in 27. The form above, with both factors
increased by \\1/n\\, reproduced all 60 conditions within Monte Carlo
error; it never lowers the criterion, and it guarantees a width of at
most \\\omega\\ at the stopping point. The package implements the form
that reproduces the published results.

**When no decision can be made the rule is "not met".** As in the other
`seq_aipe_*` check functions, when the variance estimate is unusable or
the estimate cannot be computed, the rule is reported as not met, the
limits are `NA`, `criterion_n` is at least one more than `current_n`,
the attribute `variance_unusable` is `TRUE`, and a note says why (see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)).
This is a deliberate departure from the rules as printed, which would be
met at once with an interval of zero width.

For this procedure the variance estimate, \\4 R^2 (1 - R^2)^2\\, is
unusable when the sample \\R^2\\ is 0 or 1 to within the precision of
the data (for example because the response is an exact linear function
of the predictors), which makes the estimate zero, and it cannot be
computed when the response does not vary. A positive estimate, however
small, is used as it stands: the lower bound of \\n^{-3}\\ that the
distribution free procedures in this package place on their variance
estimates is not part of this procedure. The rule is also reported as
not met, with limits of `NA` and a note, while the log transformed
interval is not yet defined (\\n = p + 2\\).

A variable whose values agree to ten significant digits counts as not
varying, as everywhere in this package: such a predictor is an error,
and such a response is the case in which the estimate cannot be
computed. Whether \\R^2\\ is 0 or 1 is likewise judged to the precision
of the data. For ordinary data that is the precision of the arithmetic.
When a variable is a large constant plus a small spread (a mean more
than about 700,000 times the standard deviation, as with long
identifiers or time stamps), the last digits of its stored values are
rounding error, and a share of the sum of squares small enough to be
that error alone is treated as zero. Subtracting a constant from such a
variable before calling avoids both and does not change \\R^2\\.

## Assumptions

The predictors and the response are assumed to be independent draws from
a multivariate normal distribution (so the predictors are random, not
fixed by design). Both standard errors come from the normal theory
variance of \\R^2\\, which is too small for heavy tailed data, so the
intervals are then too narrow however large the sample. In simulations
of the full procedure with three predictors, a population value of .50,
a target width of .10, and a 95% level (40,000 replications per
condition), coverage with the Wald interval was .946 for multivariate
normal data, .902 for a multivariate *t* distribution with 10 degrees of
freedom, .788 for a multivariate *t* distribution with 5 degrees of
freedom, and .860 when the predictors and the error were independent
right skewed variables (standardized chi square with 3 degrees of
freedom). A narrower target does not help: at a width of .05 coverage
was .948 for normal data and .782 for the *t* distribution with 5
degrees of freedom. The log transformed interval behaved the same way
(.938, .896, .773, and .851 at a width of .10). The article offers no
remedy, and this function should not be used for data that are far from
multivariate normal.

## The Null Boundary

When all \\p\\ predictors and the response are independent normal
variables and the population squared multiple correlation is zero, the
fixed sample distribution is \\R^2 \sim \mathrm{Beta}(p/2, (n-p-1)/2)\\.
In particular, \\n R^2\\ approaches a chi squared distribution with
\\p\\ degrees of freedom, rather than the normal approximation with
positive variance used by the Wald interval. Its lower limit includes
zero only when \\n R^2 \leq 4 z^2 (1-R^2)^2\\. Its limiting fixed sample
coverage is therefore \\\Pr(\chi^2_p \leq 4 z^2)\\, which is about
24.49% for \\p=20\\ and a nominal 95% level. This is a fixed sample
limit, not the coverage of the sequential stopping rule. It shows why a
larger pilot alone cannot establish nominal coverage at the null
boundary. Moreover, when \\p=20\\, no legal fixed sample size gives 95%
coverage with this Wald interval. Put \\c=4z^2\\. The inclusion event
implies \\R^2 \leq c/n\\. The beta distribution's mean and variance,
together with a one sided variance inequality, bound its probability by
\\40/(40+(20-c)^2)\\, or about 65.07%, for every \\n\>21\\. This loose
bound rules out 95% coverage at any fixed sample size; the limit above
gives the sharper large sample result.

## Asymptotic Guarantee

For fixed \\p\\ and an interior population value \\0 \< \mathrm{P}^2 \<
1\\, the coverage of the interval at the stopping point approaches the
confidence level as the target width narrows; it is not guaranteed at
any particular width, and the article notes (p. 508) that the procedure
works well as the sample size increases but not for small sample sizes.
With five predictors, multivariate normal data, and a 95% level (10,000
replications per condition), coverage with the Wald interval for a
population value of .50 was .851 at a target width of .30, .912 at .20,
.936 at .10, and .949 at .05; with the log transformed interval it was
.831, .896, .930, and .939. Among these interior conditions, Wald
coverage is lowest when the population value is large: for .90 it was
.797 at a width of .30, .848 at .20, .901 at .10, and .913 at .05. The
log transformed interval, which can stop within a few observations of
its pilot sample size, does not improve steadily as the target narrows:
for a population value of .90 it gave .925, .906, .774, and .847 at the
same four widths. The article's Tables 11 to 14 show the same pattern at
widths of .05 and .10: at the 95% level its coverage ranges from .830 to
.954 for the Wald interval and from .646 to .949 for the log transformed
interval.

Coverage also falls, much further, when the sample is small relative to
the number of predictors; see the next section. (The Monte Carlo
standard error of each coverage from a simulation of the full procedure
reported on this page for a single condition is at most .011, and at
most .005 where 10,000 replications are stated; a condition that appears
in two places comes from two separate runs.)

## Small Samples Relative to the Number of Predictors

The sample \\R^2\\ overestimates the population value by roughly \\(1 -
\mathrm{P}^2) p / n\\. In standard errors that is about \\p / (2
\mathrm{P} \sqrt{n})\\, and both intervals are built around the sample
value. Their coverage is therefore below the confidence level until
\\n\\ is large relative to \\p\\ and, with many predictors, relative to
\\p^2\\. Olkin and Finn (1995, p. 161) give their variance as an
approximation for samples in which, in the notation used here, \\2p +
5\\ is small relative to \\n\\. The stopping rule takes no account of
this: with a wide target or many predictors it can be met at the pilot
sample size.

**The note.** Whenever an interval is reported from fewer than \\p
\max(20, 2p)\\ observations, the object returned carries a note that
says so. That is fewer than 20 observations per predictor or, with more
than 10 predictors, fewer than \\2p\\ per predictor: 40 observations
with 2 predictors, 100 with 5, 200 with 10, 800 with 20. The note
depends on \\n\\ and \\p\\ alone. It does not change `is_satisfied`: the
article has no remedy for the shortfall, and none is invented here. It
marks the sample sizes at which the interval should not be relied on
without checking the design with
[`seq_aipe_R2_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_sensitivity.md).

The threshold is where the coverage of the Wald interval at a *fixed*
sample size and a 95% level reaches about .90 for population values of
.30 or more. With 2 to 40 predictors and population values of .30, .50,
.70, and .90 (400,000 samples per condition from the exact distribution
of \\R^2\\ under multivariate normality, which Gatsonis and Sampson,
1989, p. 517, give), that coverage was between .902 and .928 at exactly
\\p \max(20, 2p)\\ observations, between .847 and .904 at half as many,
and between .714 and .855 at a quarter as many; at twice as many it was
between .927 and .938. (With one predictor it was between .864 and .891
at 20 observations. At the 90% and 99% levels it was between .834 and
.870 and between .953 and .981 at the threshold.) A smaller population
value needs a larger sample: for a population value of .10 the coverage
at the threshold was between .834 and .934. The log transformed interval
does better at a fixed sample size (between .912 and .950 at the
threshold and between .799 and .957 at a quarter of it), but the
sequential procedure built on it stops earlier.

In simulations of the full sequential procedure at a 95% level
(population values of .10, .30, .50, .70, and .90; 500 to 2,000
replications per condition; every condition given the same weight),
coverage among the replications that stopped with the note depended
above all on the number of predictors. With the Wald interval it was
.793 with 2 predictors, .817 with 5, .663 with 10, .428 with 20, and
.180 with 30. Among the replications that stopped without the note it
was .915 with 1 predictor (with which the note was never attached), .933
with 2, .942 with 5, .917 with 10, .927 with 20, and .892 with 30. With
the log transformed interval coverage with the note was .466, .732,
.807, .755, .615, and .422 with 1, 2, 5, 10, 20, and 30 predictors, and
without it between .923 and .951. The conditions are not fully crossed:
the target widths were .05, .07, .10, .15, .20, .30, and .40 with 1, 2,
and 5 predictors, .10 to .40 with 10 predictors, .07 to .40 with 20, and
.05, .20, and .30 with 30 (175 conditions for each interval). A single
figure for all replications that stopped with the note therefore depends
on which conditions are pooled, and none is offered: for the Wald
interval it is .530 over all 175 conditions, .811 over those with 1, 2,
and 5 predictors, and .439 over those with widths of .20 and .30, which
every number of predictors has. The figure without the note depends on
that far less (.926, .928, and .914).

By the sample size at the stopping point relative to the threshold, for
population values of .30 or more and over all of those conditions,
coverage with the Wald interval was .780 from a quarter to a half of the
threshold, .877 from a half to the threshold, .901 from the threshold to
twice it, .928 from twice to four times it, and .949 beyond; with the
log transformed interval it was .815, .911, .938, .946, and .951. Below
a quarter of the threshold it again depended on the number of
predictors: .601 with 5, .395 with 10, .243 with 20, and .115 with 30
for the Wald interval (with 1 and 2 predictors no replication stopped
there), and .239, .587, .643, .607, .506, and .399 with 1, 2, 5, 10, 20,
and 30 for the log transformed interval.

The note is a statement about \\n\\ and \\p\\ and nothing else, so it is
neither a verdict nor a guarantee. Coverage can be at the confidence
level with it: with five predictors, a population value of .10, and a
width of .30, 90% of replications stopped with the note and coverage was
.941. And coverage can be below the confidence level without it, for the
reasons given in the previous section: at a wide target and for a small
population value. Among the replications with 1, 2, and 5 predictors
that stopped without the note, Wald coverage for a population value of
.50 was .949 at a width of .10 and .918 at .30, and for a population
value of .10 it was between .864 and .942 across the seven widths. With
a single predictor, for which the pilot sample size at these widths is
20 or more, so that the note is never attached, coverage for a
population value of .10 was between .80 and .89 at widths of .10 to .40.

Some single conditions, each from a run of its own. For a population
value of .50 at the 95% level (10,000 replications per condition; 4,000
for 40 predictors at a width of .10), coverage with the Wald interval at
a target width of .30 was .900 with 2 predictors, .845 with 5, .668 with
10, .361 with 15, .071 with 20, and 0 with 40 (where every replication
stopped at the pilot sample size of 42); at a width of .20 it was .865
with 10 predictors and .790 with 15; at a width of .10 it was .948,
.942, .934, .903, and .778 with 2, 5, 10, 20, and 40. The log
transformed interval gave .885, .837, .709, .563, .427, and .076 at a
width of .30, .829 and .748 at .20, and .936, .926, .906, .851, and .671
at .10. When the predictors are unrelated to the response (a population
value of 0), coverage at a width of .30 was .971 with 5 predictors, .711
with 10, and .040 with 20 for the Wald interval, and .967, .847, and
.273 for the log transformed interval. This case does not improve with
the sample size, because the variance \\4 \mathrm{P}^2 (1 -
\mathrm{P}^2)^2\\ on which both intervals rest is zero there: the Wald
interval contains 0 when \\n R^2 \le 4 z^2 (1 - R^2)^2\\, so as \\n\\
grows its coverage tends to the probability that a chi square variable
with \\p\\ degrees of freedom is at most \\4 z^2\\, which at the 95%
level is .991, .881, .575, and .245 for 5, 10, 15, and 20 predictors (at
a fixed sample size of 400 observations per predictor, 100,000 samples
each: .990, .877, .569, and .240).

**A larger pilot sample.** The one thing in the user's hands is `m_0` in
[`seq_aipe_R2_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_pilot.md).
With 20 predictors, a population value of .50, and a width of .30, pilot
samples of 100, 200, and 400 gave coverage of .601, .789, and .877 with
the Wald interval. With the pilot sample size at the threshold itself,
`m_0 = p * max(20, 2 * p)`, so that no interval is ever reported with
the note, coverage with the Wald interval was between .896 and .932 for
5, 10, and 20 predictors, population values of .30 to .90, and widths of
.20 and .30 (2,000 replications per condition), where the default pilot
sample size gave between .019 and .932; for a population value of .10 it
was between .853 and .953. The log transformed interval gave between
.920 and .955 (default pilot sample size: .373 to .930), and between
.875 and .966 for a population value of .10. The price is the sample
size: with 20 predictors, a population value of .50, and a width of .30
that pilot sample is 800 observations, where the rule by itself stops
near 33. This is not part of the article's procedure and is not a
default; the article's own simulations cover 2, 5, and 10 predictors at
widths of .05 and .10. Use
[`seq_aipe_R2_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_sensitivity.md)
to see what to expect at the width, number of predictors, and pilot
sample size of a planned study before relying on the procedure.

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

Bonett, D. G., & Wright, T. A. (2011). Sample size requirements for
multiple regression interval estimation. *Journal of Organizational
Behavior*, *32*(6), 822–830.
[doi:10.1002/job.717](https://doi.org/10.1002/job.717)

Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact power
and sample size calculations. *Psychological Bulletin*, *106*(3),
516–524.

Olkin, I., & Finn, J. D. (1995). Correlations redux. *Psychological
Bulletin*, *118*(1), 155–164.
[doi:10.1037/0033-2909.118.1.155](https://doi.org/10.1037/0033-2909.118.1.155)

## See also

[`seq_aipe_R2_pilot`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_pilot.md),
[`seq_aipe_R2_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_sensitivity.md);
[`seq_path`](https://yelleknek.github.io/SMSD/reference/seq_path.md) for
the whole stopping path of data in hand;
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
for the methods of the returned object;
[`seq_aipe_r`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md)
for the distribution free procedures for a single correlation; in DMAR,
the fixed sample size planner
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.html)
and the confidence interval
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.html)
(Kelley, 2008), which invert the exact sampling distribution of \\R^2\\
rather than use an approximate standard error.

Other seq_aipe_R2:
[`seq_aipe_R2_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_pilot.md),
[`seq_aipe_R2_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# How well do autonomy, supervisor support, and workload together predict
# job satisfaction? The goal is a 95% confidence interval for the
# population squared multiple correlation no wider than 0.25.
set.seed(113)
survey <- function(n) {
  predictors <- cbind(autonomy = rnorm(n), support = rnorm(n),
                      workload = rnorm(n))
  satisfaction <- drop(predictors %*% c(0.45, 0.40, -0.30)) + rnorm(n)
  list(x = predictors, y = satisfaction)
}
m <- seq_aipe_R2_pilot(width = 0.25, p = 3)
employees <- survey(m)
out <- seq_aipe_R2(employees$x, employees$y, width = 0.25)
# Not met, with a note: 32 observations are few for three predictors.
out
#> Sequential check (squared multiple correlation, Wald interval)
#> 
#>   Current n       : 32
#>   Criterion n     : 104
#>   Target width    : 0.25
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : 0.42985
#>   R2                : 0.6324142
#>   upper_limit       : 0.8349784
#>   xi_squared        : 0.3418055
#>   criterion_n       : 104
#>   current_n         : 32
#> 
#>   Stopping rule met: FALSE
#> 
#>   Note: With 32 observations and 3 predictors the sample is small
#>         relative to the number of predictors (fewer than 60
#>         observations, 20 per predictor). The sample squared multiple
#>         correlation then overestimates the population value, and in
#>         simulations the interval of this procedure contained the
#>         population value less often than the confidence level, far less
#>         often the smaller the sample. This note does not change the
#>         stopping rule. See 'Small Samples Relative to the Number of
#>         Predictors' in the help of seq_aipe_R2().

# Survey ten more employees at a time until the stopping rule is met.
while (!attr(out, "is_satisfied")) {
  more <- survey(10)
  employees <- list(x = rbind(employees$x, more$x),
                    y = c(employees$y, more$y))
  out <- seq_aipe_R2(employees$x, employees$y, width = 0.25)
}
out
#> Sequential check (squared multiple correlation, Wald interval)
#> 
#>   Current n       : 152
#>   Criterion n     : 149
#>   Target width    : 0.25
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : 0.2709712
#>   R2                : 0.3920038
#>   upper_limit       : 0.5130364
#>   xi_squared        : 0.5796315
#>   criterion_n       : 149
#>   current_n         : 152
#> 
#>   Stopping rule met: TRUE

# The same data one employee at a time, as a path. seq_path() cannot work
# out the pilot sample size for this procedure (it depends on the number
# of predictors), so the start is supplied.
path <- seq_path(seq_aipe_R2, employees, width = 0.25,
                 start = seq_aipe_R2_pilot(width = 0.25, p = 3))
attr(path, "stop_n")
#> [1] 149
# The check at the stopping point. Its notes, if it has any, are shown
# here and not when the path itself is printed.
attr(path, "final_check")
#> Sequential check (squared multiple correlation, Wald interval)
#> 
#>   Current n       : 149
#>   Criterion n     : 149
#>   Target width    : 0.25
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : 0.2688058
#>   R2                : 0.3910916
#>   upper_limit       : 0.5133775
#>   xi_squared        : 0.5800192
#>   criterion_n       : 149
#>   current_n         : 149
#> 
#>   Stopping rule met: TRUE
```
