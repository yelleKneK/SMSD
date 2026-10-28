# Minimum Risk Point Estimation of the Coefficient of Variation

Given the data collected so far, evaluates the purely sequential minimum
risk point estimation (MRPE) stopping rule of Chattopadhyay and Kelley
(2016) for the population coefficient of variation \\\kappa = \sigma /
\mu\\. The procedure is for a study in which each observation costs
money and a less accurate estimate also has a cost: it stops sampling
when the estimated total of the two is as small as it can be. It targets
a risk, not a confidence interval, so there is no confidence level.

## Usage

``` r
seq_mr_cv(
  x,
  A,
  c,
  gamma = 0.49,
  variance = c("corrected", "published"),
  na.rm = TRUE
)
```

## Arguments

- x:

  A numeric vector: the observations collected so far, at least four.
  The coefficient of variation is meaningful for a quantity measured on
  a ratio scale with a positive mean (the article assumes nonnegative
  data); it is not meaningful for centered or standardized scores or for
  differences, whose mean can be zero.

- A:

  The loss constant: the price of one unit of squared estimation error,
  in the same currency as `c`. The article builds it as the amount one
  is willing to pay to estimate the coefficient of variation to within
  \\\epsilon\\, divided by \\\epsilon^2\\ (for example \\1000 / 0.05^2 =
  400000\\). A single positive number.

- c:

  The cost of sampling one observation. A single positive number.

- gamma:

  The tuning constant of the correction term, strictly between 0 and
  1/2; default 0.49. Use the same value as in
  [`seq_mr_cv_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv_pilot.md).

- variance:

  Which estimator of the asymptotic variance \\\xi^2\\ to use:
  `"corrected"` (the default), the consistent estimator of Kelley,
  Darku, and Chattopadhyay (2018); or `"published"`, Equation 22 of
  Chattopadhyay and Kelley (2016), for reproducing that article. See
  Details.

- na.rm:

  Should missing values be removed before the computation? `TRUE` (the
  default) removes them; `FALSE` stops with an error if any are present.

## Value

An object of class `c("smsd_check", "data.frame")` with columns `term`
and `value` and one row for each of

- `cv`:

  the sample coefficient of variation, \\s_n / \bar{X}\_n\\;

- `risk`:

  the estimated risk at the current sample size, \\A\\\hat{\xi}\_n^2 /
  n + c\\n\\;

- `xi_squared`:

  the variance estimate \\\hat{\xi}\_n^2\\ from the estimator chosen by
  `variance`, bounded below by \\n^{-3}\\ (it equals \\n^{-3}\\ when the
  variance cannot be estimated);

- `optimal_n`:

  the criterion sample size, rounded up: the stopping rule is met
  exactly when `current_n` is at least `optimal_n`. When the variance
  cannot be estimated it is at least `current_n + 1`;

- `current_n`:

  the number of observations supplied.

Attributes: `is_satisfied` (whether to stop sampling: `TRUE` or `FALSE`,
never `NA`), `current_n`, `criterion_n` (equal to `optimal_n`),
`method`, `loss_A`, `cost_c`, `gamma`, `variance` (the estimator used),
`variance_unusable`, and `notes` (a character vector, printed under the
result). The attribute `variance_unusable` is `TRUE` when the variance
estimate is zero, negative, or not finite (see Details). `conf_level`
and `width` are `NULL`. When the sample mean is not positive, `cv`,
`risk`, `xi_squared`, `optimal_n`, and `criterion_n` are `NA`. See
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
for the print, summary, and conversion methods.

## Details

The function does not collect data or keep state. Call
[`seq_mr_cv_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv_pilot.md)
for the starting sample size, collect that many observations, then call
`seq_mr_cv()` after each new observation (or batch) until
`attr(out, "is_satisfied")` is `TRUE`.

**The risk and the stopping rule.** The loss of estimating \\\kappa\\ by
the sample coefficient of variation \\k_n = s_n / \bar{X}\_n\\ from
\\n\\ observations is the squared error priced at \\A\\ plus the
sampling cost, and its expectation is the risk \\R_n(\kappa) =
A\\E\[(k_n - \kappa)^2\] + c\\n\\ (Equation 11 of Chattopadhyay and
Kelley, 2016). For samples that are not too small this is approximately
\\A\\\xi^2 / n + c\\n\\ (their Equation 14), where \\\xi^2\\ is the
asymptotic variance of \\\sqrt{n}\\(k_n - \kappa)\\. The approximate
risk is smallest at \\n_c = (A/c)^{1/2}\\\xi\\ (their Equation 15),
where it equals \\2 c\\ n_c\\ (their Equation 16). Because \\\xi\\ is
unknown, the purely sequential rule (their Equation 23) stops at the
first \\n \ge m\\ with \$\$n^2 \ge \frac{A}{c}\left(\hat{\xi}\_n^2 +
n^{-2\gamma}\right),\$\$ where \\m\\ is the pilot size,
\\\hat{\xi}\_n^2\\ estimates \\\xi^2\\ from the data, and
\\n^{-2\gamma}\\ with \\0 \< \gamma \< 1/2\\ keeps the rule from
stopping too early while the estimate is still noisy. The reported
`optimal_n` is the same criterion on the scale of \\n\\, \\\lceil
(A/c)^{1/2} (\hat{\xi}\_n^2 + n^{-2\gamma})^{1/2} \rceil\\, and the
reported `risk` is \\A\\\hat{\xi}\_n^2 / n + c\\n\\. The coefficient of
variation and the estimate \\\hat{\xi}\_n^2\\ are ratios in which the
units of the data cancel, so the result is the same whatever the units
(centimeters or meters, dollars or thousands of dollars), up to rounding
error.

**Two variance estimators.** The argument `variance` chooses
\\\hat{\xi}\_n^2\\.

- `"corrected"` (the default; use this one) is the estimator of Kelley,
  Darku, and Chattopadhyay (2018, Equations 31 and 32), the sample
  version of the asymptotic variance in their Equation 30, \$\$\xi^2 =
  \frac{\mu_4}{4\mu^2\sigma^2} - \frac{\sigma^2}{4\mu^2} -
  \frac{\mu_3}{\mu^3} + \frac{\sigma^4}{\mu^4},\$\$ where \\\mu_3\\ and
  \\\mu_4\\ are the third and fourth central moments. It is consistent
  for any population with a positive mean and a finite fourth moment,
  and it is the estimator
  [`seq_aipe_cv`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md)
  uses.

- `"published"` is Equation 22 of Chattopadhyay and Kelley (2016), the
  sample version of their Equation 9, \$\$\frac{\mu_4}{4\mu^4} +
  \frac{\sigma^4}{4\mu^4} - \frac{\mu_3}{\mu^3} +
  \frac{\sigma^2}{2\mu^2}.\$\$ It is offered so that the article's
  worked example and tables, and the results of
  [`mr_cv`](https://yelleknek.github.io/DMAR/reference/mr_cv.html) in
  DMAR and `mr.cv()` in MBESS, can be reproduced.

The two expressions differ by \\\frac{1}{4}\kappa^2 (1 -
\kappa^2)(\mu_4/\sigma^4 - 3)\\. They are therefore equal for a normal
population, and for any population with no excess kurtosis, and they
differ otherwise. For gamma data with \\\kappa = 0.8\\, for example,
\\n\\ times the mean squared error of \\k_n\\ at \\n = 400\\ was 0.516
in a simulation of 200,000 samples (standard error 0.002); the corrected
expression gives 0.525 and the published one 0.304. In the full
procedure for that population (\\A\\ = 5,000,000, \\c\\ = 10; 5,000
replications) the risk is smallest at 512 observations; the mean final
sample size was 479 with the corrected estimator and about 374 with the
published one (the article's Table 1 reports 374.0). The derivation,
with code that a reader can run, is in the file `cv_variance_note.md`
installed with the package
(`system.file("cv_variance_note.md", package = "SMSD")`). Both
estimators use the unbiased estimators of \\\mu_3\\ and \\\mu_4\\
(Equations 20 and 21 of the 2016 article), and a positive estimate from
either is bounded below by \\n^{-3}\\, as in Equation 31 of the 2018
article (Equation 22 of the 2016 article has no such bound; the package
applies it to the published estimator as well).

**When the variance cannot be estimated.** If the variance estimate is
zero, negative, or not finite (every observation so far is the same, or
the same to ten significant digits, or the sample gives a negative
estimate), the function reports the stopping rule as not met, sets
`attr(out, "variance_unusable")` to `TRUE`, and adds a note; see
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md).
This is a deliberate departure from the article's rule, which read
literally would stop at once on such an estimate. A positive estimate,
however small, is used as the articles prescribe, so a small coefficient
of variation does not by itself make an estimate unusable.

What the departure costs depends on how often the estimate is negative.
For normal data that is seldom. With a coefficient of variation of 0.10,
\\A\\ = 400, and \\c\\ = 1 (a pilot size of 8), the mean final sample
size in 2,000 simulated runs was 8.05, where the rule as printed, with
the bound in place of a negative estimate, gave 8.00. It is frequent for
a quantity that takes only two or three distinct values and has a small
coefficient of variation: the asymptotic variance of such a population
is close to zero (it is \\\kappa^4\\ for two equally likely values), the
default estimate is a difference of terms that nearly cancel, and it is
negative at many stages in a row. For observations equal to 9 or 11 with
equal probability (\\\kappa\\ = 0.10) and the same constants, the rule
as printed stopped at 8 in every one of 2,000 runs. This function
reports the rule as met at 8 with probability exactly one half (when 1,
2, 5, 6, or 7 of the 8 observations are 11) and keeps sampling
otherwise: in the same 2,000 runs it stopped at 8 in 48.5 percent, the
mean final sample size was 25, the 90th percentile 51, the 99th
percentile 244, and the largest 693. (With `variance = "published"` the
mean was 8.01, because Equation 22 seldom turns negative for such data;
for this population, however, its large sample value is about fifty
times the asymptotic variance.) If the quantity you measure takes very
few distinct values, decide on a largest sample size before sampling
begins.

**Agreement with
[`mr_cv()`](https://yelleknek.github.io/DMAR/reference/mr_cv.html) in
DMAR.** With `variance = "published"` this function and
[`mr_cv`](https://yelleknek.github.io/DMAR/reference/mr_cv.html) apply
the same rule to the same number whenever Equation 22 is at least
\\n^{-3}\\, as it is in the article's example.
[`mr_cv()`](https://yelleknek.github.io/DMAR/reference/mr_cv.html)
applies Equation 22 as printed, with no lower bound, so for a smaller
value, which occurs in small samples, the two can differ, in two ways.

- Equation 22 is zero or negative. This function reports the rule as not
  met;
  [`mr_cv()`](https://yelleknek.github.io/DMAR/reference/mr_cv.html)
  uses the value and can stop. In 4,000 samples of 8 observations from a
  lognormal population (log scale standard deviation 1, \\A\\ = 400,
  \\c\\ = 1), Equation 22 was zero or negative in 223 (5.6 percent), and
  [`mr_cv()`](https://yelleknek.github.io/DMAR/reference/mr_cv.html)
  stopped on each of them.

- Equation 22 is positive but below \\n^{-3}\\. This function uses
  \\n^{-3}\\ in its place, so its `risk` row is the larger of the two,
  and so is its criterion: the decisions differ when \\n^2\\ lies
  between \\(A/c)(V_n^2 + n^{-2\gamma})\\ and \\(A/c)(n^{-3} +
  n^{-2\gamma})\\, with \\V_n^2\\ for Equation 22. The band is narrow
  (with the default \\\gamma\\, at most 1.5 percent of \\n^2\\ at \\n\\
  = 8 and 0.6 percent at \\n\\ = 13), but when \\A/c\\ puts \\n^2\\
  inside it the two functions disagree for every sample with so small an
  estimate. For the eight heights 168.2, 171.9, 166.4, 173.5, 170.8,
  169.1, 174.6, and 167.3 (Equation 22 is 0.00015, and \\8^{-3}\\ is
  0.00195) with \\c\\ = 1, both functions report the rule as met at
  \\A\\ = 400 and as not met at \\A\\ = 500. At \\A\\ = 485,
  [`mr_cv()`](https://yelleknek.github.io/DMAR/reference/mr_cv.html)
  reports it as met (63.3 against \\n^2\\ = 64) and this function as not
  met (64.1 against 64). In 4,000 samples of 8 normal observations with
  a coefficient of variation of 0.02, the two decisions were the same in
  every sample at \\A\\ = 400 and different in every sample at \\A\\ =
  485.

**When the mean is not positive.** If the sample mean is zero or
negative, the coefficient of variation is not estimable: the estimate,
risk, variance, and criterion are `NA`, the rule is reported as not met,
and a note says why. A mean that differs from zero only by rounding
error (no larger than about 2e-15 times the largest absolute
observation, as for centered or standardized scores) is treated as zero,
whatever its sign. A positive mean that is very close to zero gives a
very large coefficient of variation, and the function stops with an
error if the criterion that follows is too large to represent.

**What to expect.** The article proves that the expected final sample
size and the risk approach their optimal values as \\A/c\\ grows; at
practical sizes the agreement is approximate. In simulations of the
default procedure with \\A\\ = 5,000,000 and \\c\\ = 10 (5,000
replications each), the ratio of the mean final sample size to \\n_c\\
was 1.18 for normal and for gamma data with \\\kappa = 0.2\\ (there
\\n_c\\ is near the pilot size, so the correction term matters most),
1.01 for normal data with \\\kappa = 0.4\\, and about 0.99, 0.95, and
0.94 for gamma data with \\\kappa\\ = 0.4, 0.6, and 0.8; the ratio of
the realized risk to the minimum risk was within about 0.03 of 1 in each
case (its Monte Carlo standard error is about 0.01). When \\A/c\\ is
small and the population is strongly skewed the procedure stops early:
the variance estimate rests on the third and fourth sample moments,
which are very variable in small samples, and sampling ends the first
time the criterion is met. For gamma data with \\\kappa = 1.3\\ and \\c
= 1\\ the ratio of sample sizes was about 0.58, 0.74, 0.85, and 0.93 at
\\A\\ = 900, 9,000, 90,000, and 900,000 (0.86 and 0.94 for the last two
in a second set of runs). Use
[`seq_mr_cv_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)
to see the behavior for your own constants before collecting data.

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research*, *51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

## See also

[`seq_mr_cv_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv_pilot.md)
for the pilot size and
[`seq_mr_cv_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)
to study the procedure before collecting data; the accuracy in parameter
estimation companion
[`seq_aipe_cv`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md);
in DMAR,
[`mr_cv`](https://yelleknek.github.io/DMAR/reference/mr_cv.html) (the
same procedure with the published variance), the estimator
[`cv`](https://yelleknek.github.io/DMAR/reference/cv.html), and the
fixed sample size planner
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.html).

Other seq_mr_cv:
[`seq_mr_cv_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv_pilot.md),
[`seq_mr_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)

## Author

Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>, Ken Kelley
<kkelley@nd.edu>

## Examples

``` r
# The worked example of Chattopadhyay and Kelley (2016, pp. 637-638): the
# number of age appropriate books in the homes of third graders. A home
# visit costs 75 dollars and A = 1000 / 0.05^2.
books <- c(36, 53, 19, 11, 10, 24, 14, 65, 18, 48, 25, 35, 13, 18, 3, 41,
           5, 3, 44, 26, 13, 39, 2, 3, 26, 22, 8, 15, 12, 22, 5, 21, 23,
           40, 18)
m <- seq_mr_cv_pilot(A = 400000, c = 75)
m
#> [1] 18

# After the pilot sample the rule is not met, so sampling continues.
seq_mr_cv(books[1:m], A = 400000, c = 75)
#> Sequential check (minimum risk point estimation (coefficient of variation))
#> 
#>   Current n       : 18
#>   Criterion n     : 36
#> 
#>   cv                : 0.7391157
#>   risk              : 5354.513
#>   xi_squared        : 0.1802031
#>   optimal_n         : 36
#>   current_n         : 18
#> 
#>   Stopping rule met: FALSE

# With all 35 observations the default (corrected) variance asks for one
# more; the article's own variance reproduces its printed result (risk
# 4891.284, stop at 35).
seq_mr_cv(books, A = 400000, c = 75)
#> Sequential check (minimum risk point estimation (coefficient of variation))
#> 
#>   Current n       : 35
#>   Criterion n     : 36
#> 
#>   cv                : 0.7013904
#>   risk              : 4977.259
#>   xi_squared        : 0.2058226
#>   optimal_n         : 36
#>   current_n         : 35
#> 
#>   Stopping rule met: FALSE
seq_mr_cv(books, A = 400000, c = 75, variance = "published")
#> Sequential check (minimum risk point estimation (coefficient of variation))
#> 
#>   Current n       : 35
#>   Criterion n     : 35
#> 
#>   cv                : 0.7013904
#>   risk              : 4891.284
#>   xi_squared        : 0.1982998
#>   optimal_n         : 35
#>   current_n         : 35
#> 
#>   Stopping rule met: TRUE
#> 
#>   Note: The variance estimate is Equation 22 of Chattopadhyay and
#>         Kelley (2016), as requested with variance = "published". It
#>         agrees with the default estimator for normal data and can
#>         differ otherwise; see ?seq_mr_cv.
```
