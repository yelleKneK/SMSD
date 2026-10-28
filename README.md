
<!-- README.md is generated from README.Rmd. Edit README.Rmd, then run
     devtools::build_readme(). Every number below is computed when the file
     is knit; the simulation chunk takes a few minutes. -->

# SMSD <img src="man/figures/logo.png" align="right" height="139" alt="SMSD hex logo" />

SMSD (Sequential Methods for Study Design) is for researchers who have
to plan a study without knowing the population values that a sample size
formula asks for: a standard deviation, a correlation, or the effect
size itself. A sequential procedure does not guess them. Data are
collected in stages, and after each stage a stopping rule computed from
the data in hand says whether the goal has been reached. The goal can be
a confidence interval as narrow as requested (accuracy in parameter
estimation, AIPE), an estimate whose combined cost of estimation error
and of sampling is as small as it can be made (minimum risk point
estimation, MRPE), or the most precise estimate of a contrast that a
fixed budget can buy. The methods are general; the articles that
developed them have psychology, education, management, and related
fields in view. The confidence interval procedures for effect sizes
assume no particular distribution for the data, with the exception of
the one for the squared multiple correlation.

## Installation

``` r
# The released version, from CRAN
install.packages("SMSD")

# Or the development version, from GitHub (one of the two, not both)
# install.packages("remotes")
# remotes::install_github("yelleKneK/SMSD")
```

Install the source tarball identified by the release record when you
need this exact version. A CRAN mirror or the public GitHub branch can
serve a different version, so verify `packageVersion("SMSD")` after
installation. The tarball can be installed with
`install.packages("SMSD_1.0.0.tar.gz", repos = NULL, type = "source")`;
install its required dependencies first. The package needs DMAR 1.0.0 or
later and R 4.0.0 or later.

For the replication of published examples, see
`vignette("published-examples", package = "SMSD")`. For measured
coverage, sample sizes, and figures from difficult designs, see
`vignette("coverage-and-design", package = "SMSD")`.

## Where to Start

Start with `seq_aipe_smd()`, the standardized mean difference between
two groups. Suppose a training program is compared with a control
condition on employee engagement, and the aim is a 95% confidence
interval for the standardized mean difference that is no wider than
0.60. The pilot function gives the first sample size, and the check
function says whether sampling can stop. (The data are simulated, with a
true difference of 0.4 standard deviations, so that the example runs on
its own.)

``` r
library(SMSD)
set.seed(113)
trained <- rnorm(seq_aipe_smd_pilot(width = 0.60), mean = 4.9)
control <- rnorm(seq_aipe_smd_pilot(width = 0.60), mean = 4.5)
seq_aipe_smd(trained, control, width = 0.60)
#> Sequential check (standardized mean difference)
#> 
#>   Current n       : 7
#>   Criterion n     : 91
#>   Target width    : 0.6
#>   Confidence level: 95%
#>   Coverage        : asymptotic; assess with sensitivity analysis
#> 
#>   lower_limit       : -1.143786
#>   smd               : -0.1011596
#>   upper_limit       : 0.9414664
#>   xi_squared        : 1.980884
#>   criterion_n       : 91
#>   current_n         : 7
#> 
#>   Stopping rule met: FALSE
```

With the pilot sample of 7 people per group the rule is not met: the
criterion asks for 91 per group as things stand, and that number is
updated at every check. A study therefore adds participants and checks
again, until the rule is met. `seq_report()` then writes the result up.

``` r
out <- seq_aipe_smd(trained, control, width = 0.60)
while (!attr(out, "is_satisfied")) {
  trained <- c(trained, rnorm(1, mean = 4.9))
  control <- c(control, rnorm(1, mean = 4.5))
  out <- seq_aipe_smd(trained, control, width = 0.60)
}
seq_report(out)
#> Sample size was determined with the sequential accuracy in parameter
#> estimation procedure for the standardized mean difference (Kelley,
#> Darku, & Chattopadhyay, 2018), with a target full width of 0.6 for the
#> 95% confidence interval. The stopping rule was met at n = 85 per group.
#> At that sample size the estimate of the standardized mean difference
#> was 0.356, 95% CI [0.058, 0.655]. The confidence level is asymptotic;
#> meeting the width target does not establish nominal coverage at this
#> sample size.
#> 
#> Reference
#> Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#>     parameter estimation for a general class of effect sizes: A
#>     sequential approach. Psychological Methods, 23(2), 226–243.
#>     doi:10.1037/met0000127
```

The final interval is 0.596 wide (computed from the limits before they
are rounded for the report), within the 0.60 requested, and it contains
the 0.4 used to generate the data. The check function does not prompt
and keeps nothing between calls: the loop belongs to the study, and each
call is a calculation on the data collected so far. For data already in
hand, `seq_path()` runs the same loop in one call and `plot()` draws it.
`vignette("introduction", package = "SMSD")` is the overview, and
`vignette("step-by-step-sequential-workflow", package = "SMSD")` follows
one study from its pilot sample to its stopping point.

## Which Function Do I Need?

| What you want | Function | Source |
|----|----|----|
| **A Narrow Confidence Interval for an Effect Size** |  |  |
| Standardized mean difference | `seq_aipe_smd()` | Kelley, Darku, & Chattopadhyay (2018) |
| Pearson, Kendall, or Spearman correlation | `seq_aipe_r()` | Kelley, Darku, & Chattopadhyay (2019) |
| Squared multiple correlation (multivariate normal data) | `seq_aipe_R2()` | Kelley, Darku, & Chattopadhyay (2019) |
| Coefficient of variation | `seq_aipe_cv()` | Kelley, Darku, & Chattopadhyay (2018) |
| Simple linear regression slope | `seq_aipe_reg_coef()` | Kelley, Darku, & Chattopadhyay (2018) |
| Standardized mean | `seq_aipe_sm()` | Kelley, Darku, & Chattopadhyay (2018) |
| Contrast of group means | `seq_aipe_c()` | Kelley, Darku, & Chattopadhyay (2018) |
| Gini index | `seq_aipe_gini()` | Chattopadhyay & De (2016) |
| **Minimum Risk Point Estimation** (the smallest combined cost of error and of sampling) |  |  |
| Coefficient of variation | `seq_mr_cv()` | Chattopadhyay & Kelley (2016) |
| Standardized mean difference, with one sampling cost or one per group | `seq_mr_smd()` | Chattopadhyay & Kelley (2017) |
| **A Budget Across Groups** |  |  |
| The check, and the group to sample next | `seq_budget()`, `seq_budget_next()` | Chattopadhyay, Bandyopadhyay, Kelley, & Padalunkal (2025) |
| The allocation if the standard deviations were known | `seq_budget_optimal()` | Chattopadhyay, Bandyopadhyay, Kelley, & Padalunkal (2025) |
| A noninferiority, superiority, inferiority, or equivalence decision | `decide_seq_budget_noninferiority()`, `decide_seq_budget_superiority()`, `decide_seq_budget_inferiority()`, `decide_seq_budget_equivalence()` | Chattopadhyay, Bandyopadhyay, Kelley, & Padalunkal (2025) |
| **Classical Procedures for a Mean, a Median, or a Location** |  |  |
| A fixed width interval for a mean | `seq_fwd_norm_mean()`, `seq_fwd_dfr_mean()`; in two stages, `ts_fwd_norm_mean()`, `ts_fwd_dfr_mean()` | Mukhopadhyay & de Silva (2009) |
| The location of a negative exponential distribution | `seq_fwd_nexp_mean()`, `ts_fwd_nexp_mean()`, `seq_mr_nexp_mean()` | Mukhopadhyay & de Silva (2009) |
| An interval of bounded length for a median | `seq_blci_dfr_median()` | Mukhopadhyay & de Silva (2009) |
| A mean estimated with bounded risk | `seq_br_norm_mean()`; in two stages, `ts_br_norm_mean()`, `ts_br_exp_mean()` | Mukhopadhyay & de Silva (2009) |
| A mean estimated with minimum risk | `seq_mr_norm_mean()`, `seq_mr_exp_mean()`, `seq_mr_dfr_mean_arb()` | Mukhopadhyay & de Silva (2009) |
| **Studying a Design Before Collecting Data** |  |  |
| The final sample size, interval width, and coverage to expect | `seq_aipe_smd_sensitivity()` and its siblings for every `seq_aipe_*` function | The article of the procedure simulated |
| The final sample size and risk to expect | `seq_mr_cv_sensitivity()`, `seq_mr_smd_sensitivity()` | The article of the procedure simulated |
| The cost, coverage, and decisions to expect under a budget | `seq_budget_sensitivity()` | Chattopadhyay, Bandyopadhyay, Kelley, & Padalunkal (2025) |
| **Reading and Reporting a Result** |  |  |
| The result of one check | `print()`, `summary()`, `generics::tidy()`, `generics::glance()` | This package |
| The whole stopping path for data in hand, and its picture | `seq_path()`, `plot()` | This package |
| A paragraph for the methods or results section | `seq_report()` | This package; the paragraph cites the article of the procedure |

Each of the 25 check functions has a companion that returns its pilot
sample size and is named by adding `_pilot`, as in
`seq_aipe_smd_pilot()`. In all, the package exports 74 functions, 11 of
them `_sensitivity()` functions. `help("SMSD-package")` lists every one
of them with a link.

The `seq_aipe_*` functions take the full width of the interval as
`width`. The classical interval procedures take `d`, and what `d`
measures depends on the procedure: the half width of the interval for a
mean (`seq_fwd_norm_mean()`, `seq_fwd_dfr_mean()`, `ts_fwd_norm_mean()`,
`ts_fwd_dfr_mean()`), the full width for the location of a negative
exponential distribution (`seq_fwd_nexp_mean()`, `ts_fwd_nexp_mean()`),
and half of the largest length allowed for the interval for a median
(`seq_blci_dfr_median()`). Likewise `gamma` has two meanings. In
`seq_mr_cv()`, `seq_mr_smd()`, `seq_mr_dfr_mean_arb()`, and their pilot
functions it is the exponent of the correction term *n*<sup>−γ</sup> in
the stopping rule; in the other pilot functions of the classical
procedures it is a constant of the pilot rule that sets the size of the
pilot sample. Each help page says which.

## What the Procedures Guarantee

The `seq_aipe_*` procedures deliver the requested width by construction:
when the stopping rule is met, the interval is no wider than the target.
Coverage is a different matter. Under the relevant regularity conditions,
it is an asymptotic property, reached as the target narrows and the final
sample grows, and it can fall well short of the confidence level in the
small samples that a wide target, a
strong relationship, or many predictors produce. The sensitivity
functions show what to expect before any data are collected. Each row
below is one call.

``` r
G <- 2000
runs <- list(
  "Standardized mean difference of 0.5, width 0.8" =
    seq_aipe_smd_sensitivity(0.5, width = 0.8, G = G, seed = 113),
  "Pearson correlation of 0.4, width 0.4" =
    seq_aipe_r_sensitivity(0.4, width = 0.4, G = G, seed = 113),
  "Pearson correlation of 0.8, width 0.4" =
    seq_aipe_r_sensitivity(0.8, width = 0.4, G = G, seed = 113),
  "Regression slope of 0.5, width 0.4" =
    seq_aipe_reg_coef_sensitivity(0.5, width = 0.4, G = G, seed = 113),
  "Coefficient of variation of 0.5, lognormal data, width 0.2" =
    seq_aipe_cv_sensitivity(0.5, width = 0.2, distribution = "lognormal",
                            G = G, seed = 113),
  "Squared multiple correlation of 0.5, 5 predictors, width 0.3" =
    seq_aipe_R2_sensitivity(0.5, width = 0.3, p = 5, G = G, seed = 113),
  "Squared multiple correlation of 0.5, 15 predictors, width 0.3" =
    seq_aipe_R2_sensitivity(0.5, width = 0.3, p = 15, G = G, seed = 113)
)
row_of <- function(term) vapply(runs, function(run) run$value[run$term == term], 0)
knitr::kable(
  data.frame(names(runs), row_of("mean_n"), row_of("pct_width_within"),
             row_of("coverage")),
  digits = c(0, 0, 0, 2), row.names = FALSE,
  col.names = c("Design", "Mean final sample size",
                "Intervals within the width (%)", "Coverage of the 95% interval")
)
```

| Design | Mean final sample size | Intervals within the width (%) | Coverage of the 95% interval |
|:---|---:|---:|---:|
| Standardized mean difference of 0.5, width 0.8 | 50 | 100 | 0.94 |
| Pearson correlation of 0.4, width 0.4 | 56 | 100 | 0.78 |
| Pearson correlation of 0.8, width 0.4 | 23 | 100 | 0.70 |
| Regression slope of 0.5, width 0.4 | 82 | 100 | 0.84 |
| Coefficient of variation of 0.5, lognormal data, width 0.2 | 63 | 100 | 0.79 |
| Squared multiple correlation of 0.5, 5 predictors, width 0.3 | 82 | 100 | 0.85 |
| Squared multiple correlation of 0.5, 15 predictors, width 0.3 | 52 | 100 | 0.36 |

Every interval in every row was within its target width. Coverage was
close to the nominal level for the standardized mean difference (0.94,
with a Monte Carlo standard error of 0.005) and well below it in the
other rows, from 0.85 down to 0.36 for the squared multiple correlation
with 15 predictors. For the Pearson correlation it was lower for the
strong relationship (0.70) than for the moderate one (0.78). Each figure
comes from 2,000 simulated studies and has a Monte Carlo standard error
of at most 0.011, so another run of the same design, such as the one a
help page reports, can differ from it in the second decimal. (The sample
sizes count observations per group for the standardized mean difference,
pairs for the correlation and the slope, and observations otherwise.) A
narrower target, which means a larger study, brings coverage toward the
nominal level under those regularity conditions; a larger minimum pilot
sample (`m_0`) can help as well. The null multiple R-squared boundary is an
exception: increasing the sample size does not restore the Wald interval's
nominal coverage. The installed vignette `vignette("coverage-and-design")`
explains this boundary and compares normal model reference intervals. The
help page of each `seq_aipe_*` function states what its guarantee is and
reports what was measured for it, and `?seq_aipe_sensitivity` has a
table over more conditions.

In a small sample the variance estimate that a stopping rule needs can
be zero, negative, or not finite. The check functions then report the
rule as not met, with a note, and sampling continues. This is a
deliberate departure from the rules as printed in the articles, under
which the study would stop at once with an interval far narrower than
requested; the section “An Unusable Variance Estimate” of
`help("SMSD-package")` describes it.

The classical procedures state their guarantees on their help pages as
well. Some are exact: Stein’s two stage interval, `ts_fwd_norm_mean()`,
covers a normal mean with at least the stated probability for any first
stage. Others are asymptotic: the distribution free two stage interval,
`ts_fwd_dfr_mean()`, covers a normal mean with probability 0.940 at a
nominal 0.95 when the first stage has 30 observations and the target is
narrow, and with less for skewed data.

## Validation

- **Tests that pin published values.** The test suite (37 test files)
  compares the package with numbers printed in the source articles: rows
  of Tables 2, 3, and 4 of Kelley, Darku, and Chattopadhyay (2018); rows
  of Tables 2, 10, 11, 13, and 14 of Kelley, Darku, and Chattopadhyay
  (2019); the worked examples and rows of Table 1 of Chattopadhyay and
  Kelley (2016) and of Table 3 of Chattopadhyay and Kelley (2017); and
  the allocation example and Table 1 of Chattopadhyay, Bandyopadhyay,
  Kelley, and Padalunkal (2025). One published row is not reproduced in
  full, and its test records the difference: for Kendall’s tau (Table 8
  of the 2019 article) coverage and interval width are close to the
  article’s, but the mean final sample size runs about 2 percent above
  the published one and varies nearly twice as much from study to study.
  The cause is not resolved. Where no table exists, a procedure is run
  to its stopping point and compared with its theoretical sample size or
  its exact guarantee.
- **A cached simulation study.** A study of the standardized mean
  difference and Pearson correlation procedures under normal and skewed
  data (18 conditions, 10,000 replications each) is installed with the
  package, at
  `system.file("extdata", "aipe_simulation_results.rds", package = "SMSD")`,
  and is displayed and discussed in
  `vignette("complete-tour", package = "SMSD")`. The script that
  produces it is in the `data-raw` folder of the source repository.
- **Two notes on departures from the articles as printed.** The variance
  of the standardized mean difference corrects an omitted factor of 1/2
  in Equation 26 of Kelley, Darku, and Chattopadhyay (2018). The
  correction is installed with the package
  (`system.file("eq26_correction.pdf", package = "SMSD")`) and is posted
  by the authors at
  <https://kenkelley.org/wp-content/uploads/articles/Kelley_Darku_Chattopadhyay_Equation_26.pdf>.
  In `seq_mr_cv()` the variance of the coefficient of variation is by
  default the estimator of Kelley, Darku, and Chattopadhyay (2018),
  which is consistent for any distribution with a finite fourth moment.
  Equation 22 of Chattopadhyay and Kelley (2016) agrees with it for
  normal data, and for any population with no excess kurtosis, and can
  differ otherwise; `variance = "published"` uses that equation, so that
  the article can be reproduced. The derivation is in a note installed
  with the package
  (`system.file("cv_variance_note.md", package = "SMSD")`), and
  `?seq_mr_cv` summarizes it.
- **A comparison with the published tables.** The validation article of
  the package website (`vignettes/articles/validation.Rmd` in the source
  repository) sets the package’s results beside the published ones,
  table by table.

## How to Cite

`citation("SMSD")` gives the package entry and the six articles listed
here, and `toBibtex(citation("SMSD"))` gives them as BibTeX entries.
Please cite the package and the article that developed the method you
used. The book is the source of the classical procedures.

- Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
  (2025). A sequential approach for noninferiority or equivalence of a
  linear contrast under cost constraints. *Psychological Methods,
  30*(2), 425–439.
  [doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)
- Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
  pre-specified error bound. *Econometrics, 4*(3), 30.
  [doi:10.3390/econometrics4030030](https://doi.org/10.3390/econometrics4030030)
- Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient
  of variation with minimum risk: A sequential method for minimizing
  sampling error and study cost. *Multivariate Behavioral Research,
  51*(5), 627–648.
  [doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)
- Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized
  mean difference with minimum risk: Maximizing accuracy and minimizing
  cost with sequential estimation. *Psychological Methods, 22*(1),
  94–113. [doi:10.1037/met0000089](https://doi.org/10.1037/met0000089)
- Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
  parameter estimation for a general class of effect sizes: A sequential
  approach. *Psychological Methods, 23*(2), 226–243.
  [doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)
- Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
  accuracy in parameter estimation for population correlation
  coefficients. *Psychological Methods, 24*(4), 492–515.
  [doi:10.1037/met0000203](https://doi.org/10.1037/met0000203)
- Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
  their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## License

GPL (\>= 3).
