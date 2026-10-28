# Validation Against Published Results

Every sequential procedure in SMSD comes from a published source, and
most of those sources report what the procedure does: a simulation
table, a worked example, or a result that can be computed exactly. This
article sets the package beside those reports. It shows where the
package reproduces them, where it does not, where it differs on purpose
and why, and, at the end, what the procedures do not deliver when the
target is wide.

Nothing is simulated when this page is built. The numbers come from one
script, `data-raw/validation_study.R` in the package’s repository, which
runs each procedure through the package’s exported functions and saves
the results. The page reads that file:

``` r

results <- readRDS(file.path("..", "..", "data-raw", "validation_results.rds"))
dim(results)
#> [1] 304  15
names(results)
#>  [1] "family"           "procedure"        "source"           "condition"       
#>  [5] "quantity"         "published"        "published_se"     "published_digits"
#>  [9] "value"            "se"               "replications"     "kind"            
#> [13] "mean_n"           "note"             "z"
c(date = attr(results, "date"), SMSD = attr(results, "package_version"),
  seed = attr(results, "seed"))
#>         date         SMSD         seed 
#> "2026-10-10"      "1.0.0"   "20261010"
```

The study took 94 minutes on 6 cores (493 minutes of processor time).
Each simulated condition was run in 20 parts with their own seeds, so
the results do not depend on the number of cores, and the script can be
run again to check them (see the last section).

## How to Read the Tables

Each table sets a value from a source beside the value the package gives
for the same condition.

- **The package’s value** comes from the exported functions only. Where
  a `*_sensitivity()` function generates the article’s data, that
  function ran the procedure. Otherwise the pilot function and the check
  function were called in the loop a study would run: start at the pilot
  sample size, add one observation, check, and stop the first time the
  check reports the stopping rule as met.
- **SE** is the Monte Carlo standard error of the package’s value.
- **z** is the package’s value minus the published value, divided by the
  standard error of that difference. The standard error combines three
  things: the Monte Carlo error of the package’s value, the Monte Carlo
  error of the published value where that value is itself a simulation
  mean, and the rounding of the published value to the decimals printed.
  For the Monte Carlo error of a published mean the study uses the
  standard error that the article prints. Where an article prints none
  (the mean widths of the 2018 tables and of the squared multiple
  correlation tables, for example), it uses the standard error that the
  spread found in this study gives a mean of the article’s 5,000
  replications. A value of z between -2 and 2 is what chance alone
  produces 95 times in 100. With several thousand replications on each
  side the comparison is sharp: a z of 5 can belong to a difference in
  the fourth decimal.
- **Replications** is the number of simulated studies behind the
  package’s value. It is empty where the value is not simulated.
- **N** is the final sample size: the number of observations at which
  the procedure stopped.

The reproduction sections below hold 175 comparisons of a simulated
value from the package with a published simulation value, a population
value, or an exact value, so a few values of z beyond 2 are expected by
chance alone even where the package and the source agree exactly. (The
cached file has 304 rows in all. The others are worked examples and
closed forms, which have no sampling error; rows in which the package
departs from its source on purpose, or in which a rule as printed is
simulated in place of the package’s; and rows set beside a nominal level
or a bound and not beside a published value, most of them in the last
section.) Where the reference value is exact, a difference of more than
three standard errors can only be chance or a defect, and the study
settles which by running that condition a second time, with new seeds
and ten times the replications; both runs are shown. A row that does not
reproduce is shown as it is.

## A Narrow Confidence Interval: Kelley, Darku, and Chattopadhyay (2018)

Tables 2 to 4 of Kelley, Darku, and Chattopadhyay (2018, pp. 233–234)
report the mean final sample size, the coverage, and the mean width of
the interval for the standardized mean difference (95% intervals no
wider than 0.2) and for the coefficient of variation (90% and 95%
intervals no wider than 0.04). The rows below are the two normal
conditions of Table 2 with standardized mean differences of 0.3 and 0.5,
and the normal, gamma, and lognormal populations with a coefficient of
variation of 0.2 from Tables 3 and 4. For the standardized mean
difference the sample size is per group.

| Condition | Replications | N, article | N, package (SE) | z | Coverage, article | Coverage, package (SE) | z | Width, article | Width, package | z |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| normal, means 10 and 9.7, variance 1; width 0.2, 95% | 3,000 | 778.3 | 778.0 (0.1) | -1.4 | 0.9538 | 0.9477 (0.0041) | -1.2 | 0.1999 | 0.1999 | -1.1 |
| normal, means 10 and 9.5, variance 1; width 0.2, 95% | 3,000 | 793.6 | 793.7 (0.2) | 0.3 | 0.9456 | 0.9450 (0.0042) | -0.1 | 0.1999 | 0.1999 | -1.3 |
| normal, mean 10, variance 4; width 0.04, 90% | 10,000 | 178.8 | 179.4 (0.2) | 1.5 | 0.8844 | 0.8819 (0.0032) | -0.5 | 0.0351 | 0.0351 | 1.3 |
| gamma, shape 25, scale 0.6; width 0.04, 90% | 10,000 | 173.0 | 173.2 (0.3) | 0.4 | 0.8794 | 0.8772 (0.0033) | -0.4 | 0.0347 | 0.0347 | 0.1 |
| lognormal, log mean 1, log SD 0.1980; width 0.04, 90% | 10,000 | 182.3 | 182.9 (0.4) | 0.9 | 0.8650 | 0.8749 (0.0033) | 1.7 | 0.0350 | 0.0350 | 0.9 |
| normal, mean 10, variance 4; width 0.04, 95% | 10,000 | 241.2 | 241.8 (0.3) | 1.1 | 0.9422 | 0.9382 (0.0024) | -1.0 | 0.0363 | 0.0363 | -0.2 |
| gamma, shape 25, scale 0.6; width 0.04, 95% | 10,000 | 233.0 | 233.0 (0.4) | 0.1 | 0.9342 | 0.9367 (0.0024) | 0.6 | 0.0359 | 0.0359 | 1.0 |
| lognormal, log mean 1, log SD 0.1980; width 0.04, 95% | 10,000 | 247.4 | 248.3 (0.5) | 1.0 | 0.9234 | 0.9304 (0.0025) | 1.5 | 0.0362 | 0.0363 | 1.5 |

Of the 24 comparisons, 24 are within two standard errors and the largest
difference is 1.7 standard errors. The package reproduces these tables.

One detail of the comparison needs a note. The article’s tables print a
column s(N) described as the standard error of the mean final sample
size, and its entries are far too small to be that. For the standardized
mean difference of 0.5 the entry is 0.0026, while the standard deviation
of the final sample size in this study is 13.0, which gives a mean of
5,000 replications a standard error of 0.18. (The printed entry is close
to that standard deviation divided by 5,000.) The standard error of each
published mean final sample size in this section is therefore computed
from the standard deviation found here and the article’s 5,000
replications. The tables print no standard error for the mean width, so
its standard error is computed in the same way.

## Correlation Coefficients: Kelley, Darku, and Chattopadhyay (2019)

The rows are from Table 2 (Pearson, p. 498), Table 8 (Kendall, p. 502),
and Table 10 (Spearman, p. 503) of Kelley, Darku, and Chattopadhyay
(2019): bivariate normal data and 95% intervals. The Pearson rows were
run by
[`seq_aipe_r_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sensitivity.md);
that function studies the Pearson correlation only, so the Kendall and
Spearman rows were run through
[`seq_aipe_r_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r_pilot.md)
and
[`seq_aipe_r()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md).
The coverage of a rank correlation is counted against its own population
value (0.3333 for Kendall and 0.4826 for Spearman when ρ = 0.5, as the
article gives them).

| Condition | Replications | N, article | N, package (SE) | z | Coverage, article | Coverage, package (SE) | z | Width, article | Width, package | z |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Pearson, rho = 0.1; width 0.2 | 5,000 | 372.2 | 372.2 (0.6) | 0.0 | 0.9396 | 0.9392 (0.0034) | -0.1 | 0.1992 | 0.1994 | 1.9 |
| Pearson, rho = 0.5; width 0.2 | 10,000 | 204.0 | 205.5 (0.5) | 1.7 | 0.9012 | 0.9119 (0.0028) | 2.1 | 0.1938 | 0.1959 | 4.9 |
| Pearson, rho = 0.5; width 0.1 | 4,000 | 857.0 | 856.8 (1.3) | -0.1 | 0.9464 | 0.9487 (0.0035) | 0.5 | 0.0998 | 0.0998 | -0.3 |
| Kendall, rho = 0.1; width 0.2 | 5,000 | 170.9 | 174.0 (0.2) | 13.5 | 0.9440 | 0.9420 (0.0033) | -0.4 | 0.1984 | 0.1981 | -9.8 |
| Kendall, rho = 0.5; width 0.2 | 10,000 | 133.2 | 135.9 (0.2) | 10.6 | 0.9318 | 0.9372 (0.0024) | 1.3 | 0.1972 | 0.1968 | -12.6 |
| Spearman, rho = 0.5; width 0.2 | 4,000 | 238.8 | 239.7 (0.6) | 1.1 | 0.9300 | 0.9350 (0.0039) | 0.9 | 0.1986 | 0.1987 | 1.9 |

**Pearson.** The mean final sample sizes are within 1.7 standard errors
of the article’s. At ρ = 0.5 and a width of 0.2 the package’s studies
run slightly longer, cover slightly more often, and end with slightly
wider intervals than the article’s (205.5 against 204.0 pairs, coverage
0.912 against 0.901, mean width 0.1959 against 0.1938, the last a
difference of 4.9 standard errors). All three are the effect of a
deliberate departure described below: the package does not stop on a
variance estimate that is zero or negative. Coverage of about 0.90 at a
nominal 0.95 in that row is the article’s own finding, and the package’s
too. It is taken up in the last section.

**Kendall.** The mean final sample sizes and mean widths of these rows
do not reproduce. The package’s mean final sample size is 135.9 against
the article’s 133.2 at ρ = 0.5 and 174.0 against 170.9 at ρ = 0.1,
differences of 11 and 14 standard errors. The final sample size also
varies far more from study to study than the article’s standard errors
imply: its standard deviation is 20.4 here against 10.7 at ρ = 0.5,
about 2 times as large, and 16.0 against 2.6 at ρ = 0.1, about 6 times
as large. (With the spread found here in place of the article’s standard
errors, the two means are still 8 and 10 standard errors from the
article’s.) The differences in the means are small in practical terms,
about 2 percent of the sample size and at most 0.0004 in mean width, but
they are far outside Monte Carlo error. Coverage does agree with the
article’s: 0.937 against 0.932 at ρ = 0.5 and 0.942 against 0.944 at ρ =
0.1, differences of 1.3 and -0.4 standard errors. The package computes
the variance estimator of the article’s Equations 27 to 29 (p. 500); its
tests compare the code with a transcription of those equations and with
the large sample value for bivariate normal data. The cause is therefore
most likely a difference between the printed equations and the code
behind the article’s table, and it cannot be settled without that code.

**Spearman.** The row is reproduced within 1.9 standard errors.

## The Squared Multiple Correlation: Kelley, Darku, and Chattopadhyay (2019)

Tables 11 to 14 of the same article (pp. 504–507) cover the squared
multiple correlation with multivariate normal data: the log transformed
interval at 90% (Table 11) and 95% (Table 12), and the Wald interval at
90% (Table 13) and 95% (Table 14). The rows below are those with a
target width of 0.10, 2 and 5 predictors, and population values of 0.1,
0.5, and 0.9, run by
[`seq_aipe_R2_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2_sensitivity.md).

| Condition | Replications | N, article | N, package (SE) | z | Coverage, article | Coverage, package (SE) | z | Width, article | Width, package | z |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Table 11: 2 predictors, population value 0.1, 90% | 20,000 | 297.69 | 299.38 (0.99) | 0.8 | 0.7404 | 0.7402 (0.0031) | 0.0 | 0.0959 | 0.0958 | -0.4 |
| Table 11: 2 predictors, population value 0.5, 90% | 20,000 | 533.22 | 532.06 (0.59) | -0.9 | 0.8858 | 0.8771 (0.0023) | -1.7 | 0.0992 | 0.0991 | -1.1 |
| Table 11: 2 predictors, population value 0.9, 90% | 20,000 | 34.84 | 34.67 (0.16) | -0.5 | 0.7688 | 0.7605 (0.0030) | -1.2 | 0.0908 | 0.0904 | -1.9 |
| Table 11: 5 predictors, population value 0.1, 90% | 20,000 | 364.23 | 359.86 (0.68) | -2.9 | 0.8802 | 0.8721 (0.0024) | -1.6 | 0.0989 | 0.0987 | -1.9 |
| Table 11: 5 predictors, population value 0.5, 90% | 20,000 | 524.81 | 524.78 (0.72) | 0.0 | 0.8716 | 0.8600 (0.0025) | -2.2 | 0.0987 | 0.0985 | -1.3 |
| Table 11: 5 predictors, population value 0.9, 90% | 20,000 | 30.90 | 31.02 (0.15) | 0.3 | 0.6668 | 0.6679 (0.0033) | 0.1 | 0.0871 | 0.0873 | 0.4 |
| Table 12: 2 predictors, population value 0.1, 95% | 20,000 | 435.52 | 433.30 (1.31) | -0.8 | 0.8132 | 0.8079 (0.0028) | -0.9 | 0.0968 | 0.0965 | -2.1 |
| Table 12: 2 predictors, population value 0.5, 95% | 20,000 | 760.73 | 761.32 (0.66) | 0.4 | 0.9386 | 0.9356 (0.0017) | -0.8 | 0.0995 | 0.0994 | -1.0 |
| Table 12: 2 predictors, population value 0.9, 95% | 20,000 | 50.55 | 51.55 (0.21) | 2.1 | 0.8432 | 0.8485 (0.0025) | 0.9 | 0.0936 | 0.0938 | 1.1 |
| Table 12: 5 predictors, population value 0.1, 95% | 20,000 | 513.15 | 507.65 (0.81) | -3.1 | 0.9334 | 0.9303 (0.0018) | -0.8 | 0.0993 | 0.0992 | -1.1 |
| Table 12: 5 predictors, population value 0.5, 95% | 20,000 | 754.62 | 756.17 (0.79) | 0.9 | 0.9284 | 0.9274 (0.0018) | -0.3 | 0.0991 | 0.0991 | 0.1 |
| Table 12: 5 predictors, population value 0.9, 95% | 20,000 | 45.27 | 46.13 (0.21) | 1.9 | 0.7616 | 0.7671 (0.0030) | 0.8 | 0.0910 | 0.0912 | 0.8 |
| Table 13: 2 predictors, population value 0.1, 90% | 20,000 | 347.98 | 348.84 (0.65) | 0.6 | 0.8460 | 0.8515 (0.0025) | 1.0 | 0.0962 | 0.0962 | 0.2 |
| Table 13: 2 predictors, population value 0.5, 90% | 20,000 | 543.16 | 542.93 (0.24) | -0.4 | 0.8962 | 0.8934 (0.0022) | -0.6 | 0.0993 | 0.0993 | 1.5 |
| Table 13: 2 predictors, population value 0.9, 90% | 20,000 | 68.66 | 68.58 (0.05) | -0.7 | 0.8678 | 0.8785 (0.0023) | 2.0 | 0.0733 | 0.0735 | 1.0 |
| Table 13: 5 predictors, population value 0.1, 90% | 20,000 | 375.77 | 373.94 (0.55) | -1.5 | 0.8888 | 0.8860 (0.0022) | -0.6 | 0.0976 | 0.0975 | -1.6 |
| Table 13: 5 predictors, population value 0.5, 90% | 20,000 | 539.22 | 540.17 (0.24) | 1.8 | 0.8854 | 0.8852 (0.0023) | 0.0 | 0.0993 | 0.0993 | 1.3 |
| Table 13: 5 predictors, population value 0.9, 90% | 20,000 | 67.94 | 67.95 (0.04) | 0.1 | 0.8406 | 0.8374 (0.0026) | -0.6 | 0.0710 | 0.0709 | -0.6 |
| Table 14: 2 predictors, population value 0.1, 95% | 20,000 | 495.01 | 492.29 (0.80) | -1.5 | 0.9152 | 0.9111 (0.0020) | -0.9 | 0.0976 | 0.0975 | -0.9 |
| Table 14: 2 predictors, population value 0.5, 95% | 20,000 | 769.43 | 770.22 (0.28) | 1.3 | 0.9522 | 0.9449 (0.0016) | -2.1 | 0.0995 | 0.0995 | 1.3 |
| Table 14: 2 predictors, population value 0.9, 95% | 20,000 | 84.17 | 84.04 (0.07) | -0.8 | 0.9222 | 0.9306 (0.0018) | 2.0 | 0.0790 | 0.0790 | -0.3 |
| Table 14: 5 predictors, population value 0.1, 95% | 20,000 | 522.35 | 520.56 (0.68) | -1.2 | 0.9442 | 0.9429 (0.0016) | -0.4 | 0.0983 | 0.0983 | -0.5 |
| Table 14: 5 predictors, population value 0.5, 95% | 20,000 | 766.17 | 767.26 (0.28) | 1.7 | 0.9450 | 0.9440 (0.0016) | -0.3 | 0.0995 | 0.0995 | 1.2 |
| Table 14: 5 predictors, population value 0.9, 95% | 20,000 | 82.94 | 82.90 (0.06) | -0.3 | 0.8988 | 0.9061 (0.0021) | 1.5 | 0.0769 | 0.0771 | 1.1 |

Of the 24 mean final sample sizes, 21 are within two standard errors of
the article’s and 23 within three (the largest difference is 3.1); of
the 24 coverages, 21 are within two and 24 within three (largest 2.2).
The one mean that differs by more than three standard errors is within
5.5 observations (1.1 percent) of the article’s. The mean widths agree
with the article’s within 0.0004 in every row: of the 24 differences, 23
are within two standard errors and none is more than three (the largest
is 2.1). (The tables print no standard error for the mean width; the one
used is described under How to Read the Tables.) The low coverages in
these tables are the article’s as well as the package’s: 0.667 in the
article and 0.668 here for the 5 predictors, population value 0.9 row of
the log transformed interval at 90%. The Wald rows are reproduced by a
stopping rule that differs in one sign from the article’s Equation 63 as
printed; see the departures below.

## Minimum Risk for the Coefficient of Variation: Chattopadhyay and Kelley (2016)

Table 1 of Chattopadhyay and Kelley (2016, p. 642) gives the mean final
sample size of the minimum risk procedure for populations with
coefficients of variation from 0.2 to 0.8 (A = 5,000,000 and c = 10).
The article’s procedure estimates the variance of the sample coefficient
of variation with its Equation 22, which
[`seq_mr_cv()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
offers as `variance = "published"`. With that option the package
reproduces the eight rows below:

| Population | Article | Package (SE) | Replications | z |
|:---|---:|---:|---:|---:|
| gamma, coefficient of variation 0.2 | 118.31 | 118.25 (0.05) | 10,000 | -0.7 |
| gamma, coefficient of variation 0.4 | 182.37 | 182.36 (0.11) | 10,000 | 0.0 |
| gamma, coefficient of variation 0.6 | 248.11 | 248.42 (0.18) | 10,000 | 1.0 |
| gamma, coefficient of variation 0.8 | 374.04 | 373.16 (0.60) | 10,000 | -0.8 |
| lognormal, coefficient of variation 0.4 | 158.27 | 158.21 (0.14) | 10,000 | -0.2 |
| normal, coefficient of variation 0.2 | 123.85 | 123.73 (0.06) | 10,000 | -1.2 |
| normal, coefficient of variation 0.4 | 233.41 | 233.98 (0.18) | 10,000 | 1.8 |
| normal, coefficient of variation 0.6 | 393.92 | 393.02 (0.31) | 10,000 | -1.6 |

All 8 are within 1.8 standard errors. The default of
[`seq_mr_cv()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
is a different variance, and gives different sample sizes for
populations that are not normal; that departure has its own section
below.

The article’s worked example (p. 638) prints the output of its function
for 35 observations with A = 400,000 and c = 75:

| Data | Quantity | Article | Package | Same to the digits printed |
|:---|:---|---:|---:|:---|
| the article’s data, first 18 observations; A = 400,000, c = 75 | estimated risk | 5,964.345 | 5,964.345 | yes |
| the article’s data, first 18 observations; A = 400,000, c = 75 | coefficient of variation | 0.7391157 | 0.7391157 | yes |
| the article’s data, first 19 observations; A = 400,000, c = 75 | estimated risk | 6,224.861 | 6,224.861 | yes |
| the article’s data, first 19 observations; A = 400,000, c = 75 | coefficient of variation | 0.7113385 | 0.7113385 | yes |
| the article’s data, first 35 observations; A = 400,000, c = 75 | estimated risk | 4,891.284 | 4,891.284 | yes |
| the article’s data, first 35 observations; A = 400,000, c = 75 | coefficient of variation | 0.7013904 | 0.7013904 | yes |
| the article’s 35 observations; A = 400,000, c = 75 | first sample size at which the rule is met | 35 | 35 | yes |

Every printed value is reproduced, and
[`seq_path()`](https://yelleknek.github.io/SMSD/reference/seq_path.md)
finds the rule first met at 35 observations, as the article reports.

## Minimum Risk for the Standardized Mean Difference: Chattopadhyay and Kelley (2017)

Tables 1 and 3 of Chattopadhyay and Kelley (2017, pp. 103–104) give the
mean final sample size per group of the equal cost procedure for two
normal populations, at two settings of the loss constant and the cost.
They were run by
[`seq_mr_smd_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md).

| Condition | Article | Package (SE) | Replications | z |
|:---|---:|---:|---:|---:|
| Table 1: normal, delta = 0; A = 50,000, c = 1 | 235.0418 | 235.0391 (0.0020) | 10,000 | -0.8 |
| Table 1: normal, delta = 0.1; A = 50,000, c = 1 | 235.1774 | 235.1767 (0.0040) | 10,000 | -0.1 |
| Table 1: normal, delta = 0.2; A = 50,000, c = 1 | 235.6276 | 235.6384 (0.0065) | 10,000 | 1.0 |
| Table 1: normal, delta = 0.4; A = 50,000, c = 1 | 237.3000 | 237.2947 (0.0105) | 10,000 | -0.3 |
| Table 1: normal, delta = 0.8; A = 50,000, c = 1 | 243.7522 | 243.7192 (0.0203) | 10,000 | -0.9 |
| Table 1: normal, delta = 1.6; A = 50,000, c = 1 | 267.8542 | 267.7745 (0.0383) | 10,000 | -1.2 |
| Table 3: normal, delta = 0; A = 1,000,000, c = 500 | 50.0210 | 50.0194 (0.0014) | 10,000 | -0.7 |
| Table 3: normal, delta = 0.1; A = 1,000,000, c = 500 | 50.0392 | 50.0359 (0.0019) | 10,000 | -1.0 |
| Table 3: normal, delta = 0.2; A = 1,000,000, c = 500 | 50.0964 | 50.0912 (0.0029) | 10,000 | -1.0 |
| Table 3: normal, delta = 0.4; A = 1,000,000, c = 500 | 50.3882 | 50.3915 (0.0056) | 10,000 | 0.3 |
| Table 3: normal, delta = 0.8; A = 1,000,000, c = 500 | 51.6882 | 51.7035 (0.0092) | 10,000 | 1.0 |
| Table 3: normal, delta = 1.6; A = 1,000,000, c = 500 | 56.4876 | 56.4215 (0.0167) | 10,000 | -2.3 |

Of the 12 means, 11 are within two standard errors of the article’s
(largest difference 2.3). In the worked example of p. 106 (a sample
standardized mean difference of 1.00, A = 10,000, c = 2.40) the article
reports that the rule is met at 75 per group and at no smaller sample
size;
[`seq_mr_smd()`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)
first reports it met at 75. The article has no simulation table for its
procedure with two different costs; the script
`data-raw/smd_mrpe_simulation.R` in the repository checks that procedure
against the sample sizes that minimize the article’s risk function.

## A Contrast Under a Budget: Chattopadhyay, Bandyopadhyay, Kelley, and Padalunkal (2025)

**Optimal allocation.**
[`seq_budget_optimal()`](https://yelleknek.github.io/SMSD/reference/seq_budget_optimal.md)
gives the sample sizes and the variance that the article derives for
known standard deviations. The example of p. 430 (three groups, 30 per
observation, a budget of 10,000) and the optimal sample sizes and
minimum variances of all four blocks of Table 1 (p. 431; a budget of
5,000 and weights 1, -2, 1) make 17 printed values, and the package
agrees with 16 of them to the last digit printed:

| Design | Quantity | Article | Package | Same to the digits printed |
|:---|:---|---:|---:|:---|
| three groups, 30 per observation, budget 10,000, equal standard deviations, weights of equal size | optimal sample size per group | 111.11 | 111.11 | yes |
| gamma populations, weights 1, -2, 1 | optimal sample size, group 1 | 41.2035 | 41.2035 | yes |
| gamma populations, weights 1, -2, 1 | optimal sample size, group 2 | 83.1149 | 83.1149 | yes |
| gamma populations, weights 1, -2, 1 | optimal sample size, group 3 | 99.7385 | 99.7385 | yes |
| gamma populations, weights 1, -2, 1 | minimum variance of the contrast | 0.1885 | 0.1885 | yes |
| logistic populations, weights 1, -2, 1 | optimal sample size, group 1 | 42.2309 | 42.2309 | yes |
| logistic populations, weights 1, -2, 1 | optimal sample size, group 2 | 45.1468 | 45.1468 | yes |
| logistic populations, weights 1, -2, 1 | optimal sample size, group 3 | 94.6731 | 94.6731 | yes |
| logistic populations, weights 1, -2, 1 | minimum variance of the contrast | 2.3058 | 2.3058 | yes |
| normal populations, weights 1, -2, 1 | optimal sample size, group 1 | 96.373 | 96.373 | yes |
| normal populations, weights 1, -2, 1 | optimal sample size, group 2 | 61.3566 | 61.3566 | yes |
| normal populations, weights 1, -2, 1 | optimal sample size, group 3 | 41.7479 | 41.7479 | yes |
| normal populations, weights 1, -2, 1 | minimum variance of the contrast | 0.2869 | 0.2869 | yes |
| normal mixture populations, weights 1, -2, 1 | optimal sample size, group 1 | 75.913 | 75.913 | yes |
| normal mixture populations, weights 1, -2, 1 | optimal sample size, group 2 | 81.5976 | 81.5976 | yes |
| normal mixture populations, weights 1, -2, 1 | optimal sample size, group 3 | 89.1344 | 89.1344 | yes |
| normal mixture populations, weights 1, -2, 1 | minimum variance of the contrast | 0.8486 | 8.4788 | no |

The one exception is the minimum variance of the mixture block, printed
as 0.8486. The article’s own formula gives 8.4788 for the standard
deviations that reproduce that block’s three sample sizes, so the
printed value appears to be a misprint.

**The sequential procedure.** For the normal block of Table 1 the
article reports the mean final sample size of each group and the mean of
the estimated variance over 5,000 studies with a pilot sample of 10 per
group (p. 432).
[`seq_budget_sensitivity()`](https://yelleknek.github.io/SMSD/reference/seq_budget_sensitivity.md)
ran the same design under both of the package’s procedures.
`procedure = "first_passage"` is the package’s reading of the article’s
rule, in which a group is closed for good the first time it meets its
criterion; `procedure = "reallocate"`, the default, judges every group
on the current data at every step. The first passage option is not a
step for step copy of the article’s sequence (the help page of
[`seq_budget_next()`](https://yelleknek.github.io/SMSD/reference/seq_budget_next.md)
lists three differences), so the study also ran that sequence itself,
written out step for step from pp. 429–430 of the article in base R with
no package code. It is the third block of the table.

| Procedure | Quantity | Article | Simulated (SE) | Replications | z |
|:---|:---|---:|---:|---:|---:|
| `first_passage` | mean final sample size, group 1 | 95.868 | 95.693 (0.030) | 10,000 | -1.8 |
| `first_passage` | mean final sample size, group 2 | 61.389 | 61.429 (0.058) | 10,000 | 0.2 |
| `first_passage` | mean final sample size, group 3 | 41.304 | 41.772 (0.061) | 10,000 | 2.8 |
| `first_passage` | mean estimated variance of the contrast | 0.2819 | 0.2857 (0.0003) | 10,000 | 6.5 |
| `first_passage` | mean total cost | 4,976.9 | 4,975.5 (0.2) | 10,000 | -3.4 |
| `first_passage` | percent of studies within the budget | 100 | 100 | 10,000 |  |
| `reallocate` | mean final sample size, group 1 | 95.868 | 96.035 (0.027) | 10,000 | 1.7 |
| `reallocate` | mean final sample size, group 2 | 61.389 | 61.473 (0.073) | 10,000 | 0.4 |
| `reallocate` | mean final sample size, group 3 | 41.304 | 41.757 (0.049) | 10,000 | 2.8 |
| `reallocate` | mean estimated variance of the contrast | 0.2819 | 0.2847 (0.0003) | 10,000 | 4.7 |
| `reallocate` | mean total cost | 4,976.9 | 4,989.0 (0.1) | 10,000 | 84.3 |
| `reallocate` | percent of studies within the budget | 100 | 100 | 10,000 |  |
| Article’s sequence, step for step | mean final sample size, group 1 | 95.868 | 95.568 (0.021) | 20,000 | -3.1 |
| Article’s sequence, step for step | mean final sample size, group 2 | 61.389 | 61.690 (0.044) | 20,000 | 1.5 |
| Article’s sequence, step for step | mean final sample size, group 3 | 41.304 | 41.951 (0.040) | 20,000 | 4.1 |
| Article’s sequence, step for step | mean estimated variance of the contrast | 0.2819 | 0.2856 (0.0002) | 20,000 | 6.7 |
| Article’s sequence, step for step | mean total cost | 4,976.9 | 4,976.5 (0.1) | 20,000 | -1.3 |

The article’s table prints no mean cost. The value in the “Article”
column of those rows is the cost its mean final sample sizes imply: the
sum over the groups of the cost of an observation (38, 15, and 10) times
the mean final sample size. For the estimated variance the table prints
a mean absolute deviation (0.0285) and no standard error; z uses a
standard error for the article’s mean taken from it by the relation that
holds for a normal variable.

The table is reproduced only in part, and the article’s own sequence
does not reproduce it either. With the first passage rule the mean final
sample sizes of the first two groups agree with the article’s (-1.8 and
0.2 standard errors) and that of the third is higher (41.77 against
41.30, 2.8 standard errors); the mean of the estimated variance is
0.2857 where the article reports 0.2819 (6.5 standard errors); the
smallest attainable variance for this design is 0.2869. The article’s
sequence, step for step, is no closer: its mean final sample sizes are
-3.1, 1.5, and 4.1 standard errors from the table’s, and its mean
estimated variance is 0.2856 (6.7). The differences are small, at most
1.3 percent for the first passage option and 1.6 percent for the
article’s sequence, but the largest are outside Monte Carlo error under
both, so the package’s reading of the sequence is not what separates it
from the table. What does could not be determined from the article. The
cost that the table implies is reproduced by the article’s sequence,
which spent 4,976.5 of the 5,000 on average against the implied 4,976.9
(-1.3 standard errors). The first passage option spent 4,975.5 (-3.4),
and the default 4,989.0.

What the article claims for the procedure does hold: no study exceeded
the budget (100 percent within it under the first passage rule and 100
under the default), and the coverage of the 95% interval was 0.948 and
0.948. The default is the package’s choice and not the article’s
procedure; the help page of
[`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)
gives the reasons.

## The Gini Index: Chattopadhyay and De (2016)

The interval for the Gini index is the procedure of Chattopadhyay and De
(2016). Their simulation study draws from three populations: a gamma
distribution with shape 2.649 and rate 0.84, a lognormal distribution
with log mean 2.185 and log standard deviation 0.562, and a Pareto
distribution with scale 20,000 and shape 5. For each, the study computed
the Gini index and the asymptotic variance of its estimator (the
article’s Equation 6) by numerical integration, with no package code.
The table compares those population values with the mean of
[`gini_index()`](https://yelleknek.github.io/SMSD/reference/gini_index.md)
and the mean of the variance estimate that
[`seq_aipe_gini()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md)
reports, over 2,000 samples of 2,000 observations.

| Condition | Quantity | Population value | Package (SE) | Replications | z |
|:---|:---|---:|---:|---:|---:|
| gamma, shape 2.649, rate 0.84; samples of 2,000 | Gini index, mean estimate | 0.3308 | 0.3309 (0.0001) | 2,000 | 1.5 |
| gamma, shape 2.649, rate 0.84; samples of 2,000 | asymptotic variance, mean estimate | 0.0468 | 0.0468 (0.0001) | 2,000 | -0.1 |
| lognormal, log mean 2.185, log SD 0.562; samples of 2,000 | Gini index, mean estimate | 0.3089 | 0.3090 (0.0001) | 2,000 | 1.1 |
| lognormal, log mean 2.185, log SD 0.562; samples of 2,000 | asymptotic variance, mean estimate | 0.0532 | 0.0531 (0.0001) | 2,000 | -0.7 |
| Pareto, scale 20,000, shape 5; samples of 2,000 | Gini index, mean estimate | 0.1111 | 0.1111 (0.0001) | 2,000 | -0.6 |
| Pareto, scale 20,000, shape 5; samples of 2,000 | asymptotic variance, mean estimate | 0.0253 | 0.0252 (0.0002) | 2,000 | -0.7 |

The largest difference is 1.5 standard errors: the package estimates the
quantities the article defines. The population variances can themselves
be checked against the posted version of the article discussed next,
whose table prints the optimal fixed sample size for a half width of
0.01 at 90%: these variances give 1,266, 1,441, and 686 for the three
populations, and the table prints 1,267, 1,424, and 686.

The sequential procedure itself is compared with the simulation table of
the version of the article that its authors posted in 2015
(Chattopadhyay & De, 2015, Table 1, p. 19): 2,000 studies for each
population of a 90% interval with a half width of 0.01, which is a full
width of 0.02 in the package’s terms, starting from a pilot sample of
165. The rows were run through
[`seq_aipe_gini_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini_pilot.md)
and
[`seq_aipe_gini()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md).

| Population | Quantity | Posted version | Package (SE) | Replications | z |
|:---|:---|---:|---:|---:|---:|
| gamma, shape 2.649, rate 0.84 | mean final sample size | 1,259.5 | 1,282.7 (1.9) | 2,000 | 4.9 |
| gamma, shape 2.649, rate 0.84 | coverage of the article’s interval | 0.8780 | 0.8870 (0.0071) | 2,000 | 0.9 |
| gamma, shape 2.649, rate 0.84 | coverage of the interval the package reports | 0.8780 | 0.8830 (0.0072) | 2,000 | 0.5 |
| lognormal, log mean 2.185, log SD 0.562 | mean final sample size | 1,429.3 | 1,431.8 (4.0) | 2,000 | 0.4 |
| lognormal, log mean 2.185, log SD 0.562 | coverage of the article’s interval | 0.9015 | 0.9065 (0.0065) | 2,000 | 0.5 |
| lognormal, log mean 2.185, log SD 0.562 | coverage of the interval the package reports | 0.9015 | 0.9025 (0.0066) | 2,000 | 0.1 |
| Pareto, scale 20,000, shape 5 | mean final sample size | 654.5 | 650.2 (4.6) | 2,000 | -0.7 |
| Pareto, scale 20,000, shape 5 | coverage of the article’s interval | 0.9018 | 0.8810 (0.0072) | 2,000 | -2.2 |
| Pareto, scale 20,000, shape 5 | coverage of the interval the package reports | 0.9018 | 0.8605 (0.0077) | 2,000 | -4.1 |

The article’s interval is the estimate plus and minus the half width,
0.01.
[`seq_aipe_gini()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md)
reports a different interval, the estimate plus and minus z times its
estimated standard error, as the package’s other procedures do. When the
rule is met that interval is never wider than the target, and the
smaller the sample in which the rule is met, the narrower it is. The
table therefore counts coverage twice in the same studies: for the
article’s interval (the `gini` row plus and minus 0.01), which is what
compares with the posted table, and for the interval the package
reports.

*Final sample size.* The lognormal and Pareto rows are reproduced (0.4
and -0.7 standard errors). The gamma row is not. The package’s mean
final sample size is 1,282.7 against 1,259.5, 1.8 percent higher and 4.9
standard errors away, and its standard deviation is 83 where the posted
standard error implies 195 (for the lognormal population the two are 180
and 185, and for the Pareto population 206 and 189). The largest final
sample sizes are similar (1,634 here and 1,594 in the posted version),
which suggests that the posted gamma studies include final sample sizes
well below the smallest found here (1,064). The cause cannot be settled
without the code behind the posted table.

*Coverage.* The article’s interval covers 0.8870, 0.9065, and 0.8810 of
the time for the gamma, lognormal, and Pareto populations against the
posted 0.878, 0.9015, and 0.9018 (0.9, 0.5, and -2.2 standard errors).
The interval the package reports covers 0.8830, 0.9025, and 0.8605. For
the gamma and lognormal populations the two intervals hardly differ: the
package’s has a mean width of 0.0198 and 0.0199 where the article’s is
0.02. For the Pareto population, where the rule is met in smaller
samples (650 observations on average), its mean width is 0.0191, and it
covers 2.0 percentage points less often than the article’s interval in
the same studies. The coverage of the procedure at wide targets is in
the last section.

## Two Stage Procedures With Exact Properties

Four of the two stage procedures from Mukhopadhyay and de Silva (2009)
have properties that hold exactly, not only in large samples: Stein’s
(1945) interval for a normal mean
([`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md)),
the interval for the location of a negative exponential distribution
([`ts_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md)),
and the bounded risk procedures for a normal mean
([`ts_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md))
and an exponential mean
([`ts_br_exp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md);
the book’s Section 8.3.2, Equation 8.3.13). In each the total sample
size depends on the first stage through one statistic whose distribution
is known, so the expected sample size and the coverage or the risk can
be computed exactly. The study does that calculation in base R, with no
package code, and compares it with what the package’s functions did in
simulated studies: a first stage, the total sample size the function
asks for, a second stage of that size, and the final interval or
estimate.

| Condition | Quantity | Exact | Package (SE) | Replications | z |
|:---|:---|---:|---:|---:|---:|
| [`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md) normal, standard deviation 2; half width 1, 95%, first stage of 5 | mean total sample size | 31.39 | 31.31 (0.07) | 100,000 | -1.1 |
| [`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md) normal, standard deviation 2; half width 1, 95%, first stage of 5 | coverage | 0.9578 | 0.9585 (0.0006) | 100,000 | 1.1 |
| [`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md) normal, standard deviation 2; half width 0.5, 95%, first stage of 10 | mean total sample size | 82.38 | 82.39 (0.12) | 100,000 | 0.1 |
| [`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md) normal, standard deviation 2; half width 0.5, 95%, first stage of 10 | coverage | 0.9510 | 0.9516 (0.0007) | 100,000 | 0.9 |
| [`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md) normal, standard deviation 2; half width 0.2, 95%, first stage of 30 | mean total sample size | 418.80 | 418.63 (0.35) | 100,000 | -0.5 |
| [`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md) normal, standard deviation 2; half width 0.2, 95%, first stage of 30 | coverage | 0.9502 | 0.9478 (0.0007) | 100,000 | -3.4 |
| [`ts_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md) negative exponential, scale 2; width 0.2, 95%, first stage of 10 | mean total sample size | 36.05 | 35.97 (0.04) | 100,000 | -2.0 |
| [`ts_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md) negative exponential, scale 2; width 0.2, 95%, first stage of 10 | coverage | 0.9525 | 0.9516 (0.0007) | 100,000 | -1.3 |
| [`ts_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md) negative exponential, scale 2; width 0.05, 95%, first stage of 20 | mean total sample size | 130.29 | 130.21 (0.09) | 100,000 | -0.9 |
| [`ts_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md) negative exponential, scale 2; width 0.05, 95%, first stage of 20 | coverage | 0.9506 | 0.9509 (0.0007) | 100,000 | 0.4 |
| [`ts_br_exp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md) exponential, mean 2; A = 100, n\* = 25, first stage of 10 | mean total sample size | 48.03 | 48.02 (0.07) | 200,000 | 0.0 |
| [`ts_br_exp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md) exponential, mean 2; A = 100, n\* = 25, first stage of 10 | risk relative to the bound | 0.9513 | 0.9536 (0.0039) | 200,000 | 0.6 |
| [`ts_br_exp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md) exponential, mean 2; A = 100, n\* = 44, first stage of 14 | mean total sample size | 72.22 | 72.27 (0.09) | 200,000 | 0.5 |
| [`ts_br_exp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md) exponential, mean 2; A = 100, n\* = 44, first stage of 14 | risk relative to the bound | 0.9901 | 0.9926 (0.0044) | 200,000 | 0.6 |
| [`ts_br_exp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md) exponential, mean 2; A = 100, n\* = 100, first stage of 20 | mean total sample size | 146.21 | 146.12 (0.15) | 200,000 | -0.6 |
| [`ts_br_exp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md) exponential, mean 2; A = 100, n\* = 100, first stage of 20 | risk relative to the bound | 0.9401 | 0.9398 (0.0041) | 200,000 | -0.1 |
| [`ts_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md) normal, standard deviation 2; squared error loss, A = 100, n\* = 50, first stage of 10 | mean total sample size | 64.79 | 64.79 (0.10) | 100,000 | 0.0 |
| [`ts_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md) normal, standard deviation 2; squared error loss, A = 100, n\* = 50, first stage of 10 | risk relative to the bound | 0.9838 | 0.9822 (0.0054) | 100,000 | -0.3 |
| [`ts_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md) normal, standard deviation 2; absolute error loss, A = 100, n\* = 50, first stage of 10 | mean total sample size | 60.37 | 60.39 (0.09) | 100,000 | 0.1 |
| [`ts_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md) normal, standard deviation 2; absolute error loss, A = 100, n\* = 50, first stage of 10 | risk relative to the bound | 0.9932 | 0.9935 (0.0026) | 100,000 | 0.1 |
| [`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md) normal, standard deviation 2; half width 0.2, 95%, first stage of 30 (second run) | mean total sample size | 418.80 | 418.90 (0.11) | 1,000,000 | 1.0 |
| [`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md) normal, standard deviation 2; half width 0.2, 95%, first stage of 30 (second run) | coverage | 0.9502 | 0.9502 (0.0002) | 1,000,000 | 0.4 |

Of the 20 simulated values of the first runs, 19 are within two standard
errors of the exact ones and 19 within three. Under the rule described
at the top of the page, a condition with a value more than three
standard errors from its exact value was run a second time. That
happened once: the coverage of
[`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md)
(half width 0.2, 95%, first stage of 30) was 0.9478 in the first 100,000
studies, 3.4 standard errors below the exact 0.9502, and 0.9502 in the
second 1,000,000 (z = 0.4). The second run agrees with the exact value,
so the first difference was chance. The exact values themselves show the
guarantees: coverage is at least 0.95 in every interval row (from 0.9502
to 0.9578), and the risk is at most the bound in every bounded risk row
(from 0.940 to 0.993 of it). The row with n\* = 44 and a first stage of
14 was chosen because it is where
[`ts_br_exp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_exp_mean.md)
comes closest to its bound: its help page reports an exact calculation
over first stages of 3 to 100 and values of n\* up to 2,000 in which the
risk never exceeded the bound and was largest, 0.990 of it, for a first
stage of 14 or 15 with n\* near 44. The exact risk in that row is 0.9901
of the bound. For this procedure the bound rests on computation of that
kind and not on a theorem, as the help page explains.

The other classical procedures (the purely sequential procedures for a
mean, the purely sequential interval for a median, and the distribution
free two stage interval) have asymptotic guarantees only. No published
table is compared for them here; the last section measures them at wide
targets.

## Where the Package Departs From an Article

In several places the package does something other than what a source
article prints. This section gives the evidence for four of them, each
deliberate and documented in the help page of the function concerned.
Three more are described elsewhere: the interval that
[`seq_aipe_gini()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md)
reports, which is not the interval of its source article (the section on
the Gini index above and the last section), the default allocation rule
of the budget procedure (the section on the budget above and the help
page of
[`seq_budget()`](https://yelleknek.github.io/SMSD/reference/seq_budget.md)),
and the power to which one factor is raised in the two cost sample sizes
of the 2017 article (the help page of
[`seq_mr_smd()`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)).

### A Variance Estimate That Is Zero or Negative

The estimate of the asymptotic variance on which a stopping rule rests
can be zero or negative in a small sample. The articles replace it with
a small positive number, n⁻³, and evaluate the rule, under which a study
can stop at once with an interval far narrower than its data support.
The package reports such a check as not met, so sampling continues (a
positive estimate, however small, is used as the articles prescribe).
The study followed the same 10,000 simulated studies under both rules
for the Pearson correlation at ρ = 0.5 and a width of 0.2, taking the
rule as printed from the flag that
[`seq_aipe_r()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_r.md)
sets when the variance estimate is unusable.

| Rule | Replications | N, article | N, simulated (SE) | z | Coverage, article | Coverage, simulated (SE) | z | Width, article | Width, simulated | z |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Rule as printed | 10,000 | 204.0 | 205.2 (0.5) | 1.3 | 0.9012 | 0.9021 (0.0030) | 0.2 | 0.1938 | 0.1940 | 0.4 |
| Package | 10,000 | 204.0 | 206.7 (0.5) | 3.0 | 0.9012 | 0.9089 (0.0029) | 1.5 | 0.1938 | 0.1961 | 5.5 |

Under the rule as printed, 1.73 percent of the studies stopped on an
unusable variance estimate. With those stops the rule as printed gives
the article’s row (mean final sample size at 1.3 standard errors,
coverage at 0.2, mean width at 0.4); without them the package’s studies
are a little longer, cover a little more often (0.909 against 0.902),
and end with intervals that are wider on average (0.1961 against
0.1940), because the intervals of the early stops were far narrower than
the target.

### The Variance of the Coefficient of Variation

Chattopadhyay and Kelley (2016, Equation 9, p. 631) and Kelley, Darku,
and Chattopadhyay (2018, Equation 30, p. 233) give different expressions
for the asymptotic variance of the sample coefficient of variation. The
2018 expression minus the 2016 one is (kurtosis − 3) κ² (1 − κ²) / 4,
where κ is the coefficient of variation: the two agree when the kurtosis
of the population is 3, as it is for a normal population (and when κ is
exactly 1), and differ otherwise. For a gamma population with a
coefficient of variation of 0.8 the two expressions give 0.3036 and
0.5248. Simulation settles which describes the estimator: n times its
mean squared error in samples of 2,000 was 0.5242 (standard error
0.0024).

| Compared with | Estimate | Formula | Estimate (SE) | Replications | z |
|:---|:---|---:|---:|---:|---:|
| Equation 30, p. 233 | sample coefficient of variation (base R); samples of 2,000 | 0.5248 | 0.5242 (0.0024) | 100,000 | -0.2 |
| Equation 9, p. 631 | sample coefficient of variation (base R); samples of 2,000 | 0.3036 | 0.5242 (0.0024) | 100,000 | 93.5 |
| Equation 30, p. 233 | seq_mr_cv(variance = “corrected”), row xi_squared; samples of 20,000 | 0.5248 | 0.5240 (0.0022) | 400 | -0.4 |
| Equation 9, p. 631 | seq_mr_cv(variance = “published”), row xi_squared; samples of 20,000 | 0.3036 | 0.3032 (0.0011) | 400 | -0.3 |

The simulated value agrees with Equation 30 of the 2018 article and is
94 standard errors from Equation 9 of the 2016 article. The default of
[`seq_mr_cv()`](https://yelleknek.github.io/SMSD/reference/seq_mr_cv.md)
therefore estimates the 2018 expression, and `variance = "published"`
estimates the 2016 one, as the last two rows of the table confirm. The
choice changes the sample size at which the procedure stops:

| Population | Mean N, article | Optimal n, article | Optimal n, package | Mean N, published variance | Mean N, default variance | Risk efficiency, published | Risk efficiency, default |
|:---|---:|---:|---:|---:|---:|---:|---:|
| gamma, coefficient of variation 0.2 | 118.3 | 97 | 102 | 118.3 | 120.4 | 1.021 | 1.015 |
| gamma, coefficient of variation 0.4 | 182.4 | 174 | 216 | 182.4 | 212.8 | 1.033 | 0.999 |
| gamma, coefficient of variation 0.6 | 248.1 | 246 | 350 | 248.4 | 333.1 | 1.042 | 0.997 |
| gamma, coefficient of variation 0.8 | 374.0 | 390 | 513 | 373.2 | 478.0 | 1.017 | 0.996 |
| lognormal, coefficient of variation 0.4 | 158.3 | 148 | 268 | 158.2 | 240.9 | 1.137 | 0.996 |
| normal, coefficient of variation 0.2 | 123.9 | 104 | 104 | 123.7 | 122.2 | 1.013 | 1.017 |
| normal, coefficient of variation 0.4 | 233.4 | 230 | 230 | 234.0 | 232.1 | 1.004 | 1.000 |
| normal, coefficient of variation 0.6 | 393.9 | 394 | 394 | 393.0 | 392.0 | 0.998 | 1.012 |

“Optimal n” is the sample size that minimizes the risk when the
population is known: the article computes it from its Equation 9 and the
package from Equation 30, rounded up. For the normal populations the two
agree and so do the two procedures, within 1.2 percent in mean final
sample size. For the gamma and lognormal populations the article’s
optimal sample size is too small, and its procedure stops early: at a
coefficient of variation of 0.8 the published variance stops at a mean
of 373 observations and the default at 478, where the risk is smallest
near 513. Risk efficiency is the realized risk divided by the smallest
attainable risk (standard errors of about 0.007). The risk is flat near
its minimum, so stopping early costs less than the sample sizes suggest:
with the published variance the risk efficiency is 1.02 to 1.04 for the
gamma populations and 1.14 for the lognormal one, against 1.00 to 1.01
with the default. The option exists so that the article’s results, and
those of
[`DMAR::mr_cv()`](https://yelleknek.github.io/DMAR/reference/mr_cv.html),
which implements the article’s rule, can be reproduced.

### Equation 26 of the 2018 Article

Equation 26 of Kelley, Darku, and Chattopadhyay (2018, p. 232), the
asymptotic variance of the standardized mean difference, omits a factor
of one half in its second term. The authors’ correction is posted at
<https://kenkelley.org/wp-content/uploads/articles/Kelley_Darku_Chattopadhyay_Equation_26.pdf>
and ships with the package
(`system.file("eq26_correction.pdf", package = "SMSD")`). The two forms
differ only when the groups differ in mean and in skewness. For two
populations with unit variance, a standardized mean difference of 0.5,
and skewness in opposite directions (a standardized chi square with 3
degrees of freedom and its mirror image), the corrected form gives
1.3710 and the printed form 0.5545.

| Compared with | Estimate | Formula | Estimate (SE) | Replications | z |
|:---|:---|---:|---:|---:|---:|
| Equation 26 as corrected | sample standardized mean difference (base R); 1,000 per group | 1.3710 | 1.3780 (0.0097) | 40,000 | 0.7 |
| Equation 26 as printed, p. 232 | sample standardized mean difference (base R); 1,000 per group | 0.5545 | 1.3780 (0.0097) | 40,000 | 85.1 |
| Equation 26 as corrected | seq_aipe_smd(), row xi_squared; 5,000 per group | 1.3710 | 1.3713 (0.0010) | 400 | 0.3 |
| Equation 26 as printed, p. 232 | seq_aipe_smd(), row xi_squared; 5,000 per group | 0.5545 | 1.3713 (0.0010) | 400 | 795.0 |

In samples of 1,000 per group, n times the variance of the sample
standardized mean difference was 1.378 (standard error 0.010), which is
0.7 standard errors from the corrected value and 85 from the printed
one. The variance estimate of
[`seq_aipe_smd()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md)
averaged 1.3713: it estimates the corrected form.

### The Stopping Rule for the Squared Multiple Correlation

For the Wald interval the article’s Equation 63 (p. 505), as printed,
stops at the first n that is at least (16 z² / w²) (R² + 1/n) (1 − (R² +
1/n))², where w is the target width and z the normal quantile for the
confidence level.
[`seq_aipe_R2()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md)
uses (R² + 1/n) (1 − R² + 1/n)² in its place: the correction 1/n is
added inside the second factor, where the printed form subtracts it. The
two forms differ most when the population value is large. The study
evaluated the printed form on the estimates that
[`seq_aipe_R2()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_R2.md)
reports at each stage, for the two rows with 2 predictors and a
population value of 0.9:

| Equation 63 as printed | Replications | N, article | N, simulated (SE) | z | Coverage, article | Coverage, simulated (SE) | z | Width, article | Width, simulated | z |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Table 13, 90% | 5,000 | 68.66 | 66.34 (0.03) | -22.8 | 0.8678 | 0.8766 (0.0047) | 1.3 | 0.0733 | 0.0759 | 7.9 |
| Table 14, 95% | 5,000 | 84.17 | 79.91 (0.06) | -26.4 | 0.9222 | 0.9366 (0.0034) | 2.8 | 0.0790 | 0.0823 | 10.9 |

The printed form gives mean final sample sizes of 66.34 and 79.91, which
are 23 and 26 standard errors below the article’s 68.66 and 84.17 and
barely above the pilot sample sizes of 66 and 79. The form the package
uses gives 68.58 and 84.04 for the same rows (-0.7 and -0.8 standard
errors from the article), and reproduces the other Wald rows shown
earlier. The package therefore implements the form that the article’s
tables were evidently computed with, and its help page says so.

## What the Procedures Do Not Guarantee at Wide Targets

The sequential intervals for effect sizes deliver the requested width by
construction: the final interval is never wider than the target. Their
coverage is a different kind of property. The articles prove that it
approaches the confidence level as the target width shrinks, that is, as
the sample the rule asks for grows. They do not prove, and it is not
true, that coverage is close to the confidence level at every width.
When the rule is met in a small sample, coverage can be far below
nominal.

The table gives the coverage of 95% intervals measured in this study,
with the default pilot sample size unless the row says otherwise, and
the mean final sample size (per group for the standardized mean
difference and the contrast, pairs for the correlations and the slope).
“Skewed” data are standardized chi square variables with 3 degrees of
freedom. The help pages of the functions report separate runs of some of
these designs, often rounded to two decimals. A figure there can differ
from the one here by up to about 0.02, several times the standard errors
shown here, because each run has its own Monte Carlo error and some of
those runs are smaller than these.

| Procedure and condition | Mean N | Coverage (SE) | Replications | Note |
|:---|---:|---:|---:|:---|
| [`seq_aipe_smd()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md) standardized mean difference 0.5, normal; width 0.8 | 50.6 | 0.950 (0.002) | 10,000 | variance estimate unusable at one or more checks in 0.06% of studies |
| [`seq_aipe_smd()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md) standardized mean difference 0.5, skewed; width 0.4 | 207.3 | 0.944 (0.002) | 10,000 |  |
| [`seq_aipe_sm()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_sm.md) standardized mean 0.4, normal; width 0.5 | 67.6 | 0.943 (0.002) | 10,000 | variance estimate unusable at one or more checks in 0.03% of studies |
| [`seq_aipe_c()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_c.md) contrast 0.5 (weights 1, -0.5, -0.5), normal; width 1 | 23.4 | 0.935 (0.002) | 10,000 |  |
| `seq_aipe_r(), Pearson` correlation 0.4, normal; width 0.4 | 56.1 | 0.782 (0.004) | 10,000 | variance estimate unusable at one or more checks in 15.80% of studies |
| `seq_aipe_r(), Pearson` correlation 0.4, normal; width 0.4, m_0 = 40 | 64.1 | 0.898 (0.003) | 10,000 |  |
| `seq_aipe_r(), Pearson` correlation 0.4, normal; width 0.2 | 263.5 | 0.932 (0.004) | 5,000 | variance estimate unusable at one or more checks in 0.88% of studies |
| `seq_aipe_r(), Pearson` correlation 0.4, skewed; width 0.4 | 56.7 | 0.706 (0.005) | 10,000 | variance estimate unusable at one or more checks in 25.87% of studies |
| `seq_aipe_r(), Pearson` correlation 0.8, normal; width 0.4 | 22.9 | 0.690 (0.005) | 10,000 | variance estimate unusable at one or more checks in 71.40% of studies |
| `seq_aipe_r(), Pearson` correlation 0.8, normal; width 0.2 | 42.7 | 0.691 (0.005) | 10,000 | variance estimate unusable at one or more checks in 45.30% of studies |
| `seq_aipe_r(), Pearson` correlation 0.9, normal; width 0.1 | 75.4 | 0.435 (0.005) | 10,000 | variance estimate unusable at one or more checks in 80.69% of studies |
| `seq_aipe_r(), Kendall` bivariate normal, rho = 0.8; width 0.4 | 25.0 | 0.907 (0.003) | 10,000 |  |
| `seq_aipe_r(), Spearman` bivariate normal, rho = 0.8; width 0.4 | 21.8 | 0.853 (0.004) | 10,000 |  |
| [`seq_aipe_reg_coef()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md) slope 0.5, normal; width 0.4 | 81.1 | 0.833 (0.004) | 10,000 | variance estimate unusable at one or more checks in 20.28% of studies |
| [`seq_aipe_reg_coef()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md) slope 0.5, normal; width 0.4, m_0 = 40 | 90.0 | 0.915 (0.003) | 10,000 |  |
| [`seq_aipe_reg_coef()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md) slope 0.5, normal; width 0.2 | 376.9 | 0.944 (0.003) | 5,000 | variance estimate unusable at one or more checks in 1.28% of studies |
| [`seq_aipe_reg_coef()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md) slope 0.5, skewed; width 0.4 | 64.8 | 0.712 (0.005) | 10,000 | variance estimate unusable at one or more checks in 36.90% of studies |
| [`seq_aipe_reg_coef()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_reg_coef.md) slope 2, normal; width 0.2 | 217.5 | 0.610 (0.007) | 5,000 | variance estimate unusable at one or more checks in 65.28% of studies |
| [`seq_aipe_cv()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md) coefficient of variation 0.5, gamma; width 0.2 | 52.2 | 0.880 (0.003) | 10,000 |  |
| [`seq_aipe_cv()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md) coefficient of variation 0.5, gamma; width 0.1 | 216.6 | 0.916 (0.003) | 10,000 |  |
| [`seq_aipe_cv()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md) coefficient of variation 0.5, lognormal; width 0.2 | 61.6 | 0.786 (0.004) | 10,000 |  |
| [`seq_aipe_cv()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_cv.md) coefficient of variation 0.5, lognormal; width 0.1 | 322.2 | 0.858 (0.005) | 5,000 |  |
| [`seq_aipe_gini()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md) Gini index 0.3, lognormal; width 0.2 | 30.5 | 0.896 (0.003) | 10,000 |  |
| [`seq_aipe_gini()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md) Gini index 0.3, lognormal; width 0.1 | 87.8 | 0.916 (0.003) | 10,000 |  |
| `seq_aipe_R2(), Wald interval` squared multiple correlation 0.5, 5 predictors; width 0.3 | 82.3 | 0.845 (0.004) | 10,000 |  |
| `seq_aipe_R2(), Wald interval` squared multiple correlation 0.5, 10 predictors; width 0.2 | 183.9 | 0.867 (0.003) | 10,000 |  |
| `seq_aipe_R2(), Wald interval` squared multiple correlation 0.5, 10 predictors; width 0.3 | 70.8 | 0.665 (0.005) | 10,000 |  |
| `seq_aipe_R2(), Wald interval` squared multiple correlation 0.5, 15 predictors; width 0.3 | 51.6 | 0.354 (0.005) | 10,000 |  |
| [`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md) normal, standard deviation 2; half width 0.5, first stage of 30 | 62.0 | 0.943 (0.002) | 20,000 |  |
| [`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md) normal, standard deviation 2; half width 0.2, first stage of 30 | 384.3 | 0.940 (0.002) | 20,000 |  |
| [`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md) exponential (skewness 2), standard deviation 2; half width 0.5, first stage of 30 | 62.7 | 0.917 (0.002) | 20,000 |  |
| [`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md) exponential (skewness 2), standard deviation 2; half width 0.2, first stage of 30 | 382.9 | 0.911 (0.002) | 20,000 |  |
| [`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md) lognormal (skewness 4), standard deviation 2; half width 0.5, first stage of 30 | 63.7 | 0.926 (0.002) | 20,000 |  |
| [`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md) lognormal (skewness 4), standard deviation 2; half width 0.2, first stage of 30 | 388.6 | 0.877 (0.002) | 20,000 |  |

What the table shows, stated plainly:

- **Where the shortfall is small.** For the standardized mean
  difference, the standardized mean, and the contrast, coverage is
  between 0.935 and 0.950 in the 4 conditions above. The lowest is the
  contrast at a width of 1, which stops with about 23 observations per
  group.
- **The Pearson correlation.** At a correlation of 0.4 coverage is 0.78
  at a width of 0.4 and 0.93 at 0.2. A stronger correlation is estimated
  with less sampling error, so the rule is met sooner and coverage is
  worse at the same width: 0.69 and 0.69 at 0.8, and 0.43 at a
  correlation of 0.9 with a width of 0.1. Skewed data lower it further
  (0.71 at 0.4 and a width of 0.4).
- **The regression slope** behaves like the correlation: 0.83 for a
  slope of 0.5 at a width of 0.4, 0.71 with skewed data, and 0.61 for a
  slope of 2 at a width of 0.2.
- **The rank correlations** fall short at a strong association too: 0.91
  (Kendall) and 0.85 (Spearman) at ρ = 0.8 and a width of 0.4.
- **The coefficient of variation** is sensitive to skewness. At a value
  of 0.5 coverage is 0.88 for gamma data and 0.79 for lognormal data at
  a width of 0.2, and still 0.86 for lognormal data at 0.1, with a mean
  of 322 observations.
- **The Gini index** covers 0.90 and 0.92 of the time at widths of 0.2
  and 0.1. That is the interval
  [`seq_aipe_gini()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md)
  reports, which is not the interval of Chattopadhyay and De (2016); the
  two are set side by side below.
- **The squared multiple correlation** fails badly when the number of
  predictors is not small beside the sample size the width calls for:
  for a population value of 0.5 at a width of 0.3 coverage is 0.85 with
  5 predictors, 0.66 with 10, and 0.35 with 15, and at a width of 0.2 it
  is 0.87 with 10. The procedure also assumes multivariate normal data.
- **The distribution free two stage interval for a mean**
  ([`ts_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_dfr_mean.md))
  plugs the first stage variance into the known variance sample size
  with no allowance for having estimated it. With a first stage of 30
  its coverage for a narrow target is 0.940 for normal data (the exact
  value is 0.941), 0.911 for exponential data, and 0.877 for lognormal
  data with a skewness of 4. Narrowing the target does not remove the
  shortfall while the first stage stays at 30; a larger first stage
  reduces it.

### The Gini Index: Two Intervals

The interval of Chattopadhyay and De (2016) is the estimate plus and
minus half the target width.
[`seq_aipe_gini()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md)
reports the estimate plus and minus z times its estimated standard
error, which the stopping rule keeps below the target, and well below it
when the rule is met in a small sample. The study ran the two wide
conditions above again through
[`seq_aipe_gini_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini_pilot.md)
and
[`seq_aipe_gini()`](https://yelleknek.github.io/SMSD/reference/seq_aipe_gini.md),
with new data, and counted the coverage of both intervals in the same
studies.

| Condition | Quantity | Mean N | Coverage (SE) | Replications | Note |
|:---|:---|---:|---:|---:|:---|
| lognormal; width 0.2 | coverage of the interval the package reports | 30.4 | 0.890 (0.003) | 10,000 | mean width 0.1405 |
| lognormal; width 0.2 | coverage of the article’s interval | 30.4 | 0.986 (0.001) | 10,000 |  |
| lognormal; width 0.1 | coverage of the interval the package reports | 88.0 | 0.916 (0.003) | 10,000 | mean width 0.0866 |
| lognormal; width 0.1 | coverage of the article’s interval | 88.0 | 0.960 (0.002) | 10,000 |  |

At a target width of 0.2 the interval the package reports has a mean
width of 0.1405 and covers 0.890 of the time; the article’s interval
covers 0.986 in the same studies. At 0.1 the two figures are 0.916 (mean
width 0.0866) and 0.960. At these targets the shortfall belongs to the
interval the package reports and not to the article’s procedure, whose
interval covers more often than its nominal level here. The article’s
interval can be formed from the package’s output: it is the `gini` row
plus and minus half of `width`, once the rule is met.

### The Purely Sequential Classical Procedures

The purely sequential procedures from Mukhopadhyay and de Silva (2009)
check their rule after every observation. Their confidence level, their
risk bound, or their optimality holds in the limit as the target
narrows, and none of them has a `*_sensitivity()` function, so the study
ran each through its pilot function and its check function, with the
default pilot, at a target near the least favorable one that its help
page reports. In the table n\* is the sample size that would be needed
if the population were known. For an interval the target is the nominal
0.95. For the bounded risk procedure it is a risk of at most 1 times the
bound, and for the minimum risk procedures it is a risk close to 1 times
the smallest attainable.

| Procedure and condition | Quantity | Mean N | Value (SE) | Replications |
|:---|:---|---:|---:|---:|
| [`seq_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md) normal, standard deviation 2; half width 0.77 (n\* = 26) | coverage | 24.2 | 0.923 (0.001) | 40,000 |
| [`seq_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md) normal, standard deviation 2; half width 0.25 (n\* = 246) | coverage | 245.0 | 0.948 (0.002) | 10,000 |
| [`seq_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md) normal, standard deviation 2; half width 0.5 (n\* = 61) | coverage | 59.9 | 0.939 (0.002) | 20,000 |
| [`seq_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md) exponential, standard deviation 2; half width 0.5 (n\* = 61) | coverage | 52.4 | 0.854 (0.002) | 20,000 |
| [`seq_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md) exponential, standard deviation 2; half width 0.25 (n\* = 246) | coverage | 238.8 | 0.934 (0.002) | 10,000 |
| [`seq_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean.md) negative exponential, scale 2; width 0.24 (n\* = 25) | coverage | 24.4 | 0.932 (0.001) | 40,000 |
| [`seq_blci_dfr_median()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md) standard normal; d = 0.2456 (n\* = 100) | coverage | 97.5 | 0.924 (0.003) | 10,000 |
| [`seq_blci_dfr_median()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md) unit exponential; d = 0.196 (n\* = 100) | coverage | 95.6 | 0.917 (0.003) | 10,000 |
| [`seq_blci_dfr_median()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md) standard normal; d = 0.08685 (n\* = 800) | coverage | 761.8 | 0.932 (0.004) | 5,000 |
| [`seq_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md) normal, standard deviation 3; squared error loss, A = 2, w = 0.65 (n\* = 28) | risk relative to the bound | 25.8 | 1.217 (0.009) | 50,000 |
| [`seq_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md) normal, standard deviation 3; absolute error loss, A = 2, w = 0.92 (n\* = 27) | risk relative to the bound | 25.2 | 1.086 (0.004) | 50,000 |
| [`seq_mr_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean.md) normal, standard deviation 1; A = 100, c = 1 (n\* = 10) | risk relative to the smallest attainable | 9.4 | 1.070 (0.004) | 50,000 |
| [`seq_mr_exp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean.md) exponential, mean 1; A = 225, c = 1 (n\* = 15) | risk relative to the smallest attainable | 14.5 | 1.150 (0.004) | 50,000 |

- **The fixed width interval for a normal mean**
  ([`seq_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_norm_mean.md))
  covers 0.923 of the time when about 26 observations are needed and
  0.948 when about 246 are. Stein’s two stage interval,
  [`ts_fwd_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_norm_mean.md),
  is the alternative whose level is guaranteed.
- **The distribution free interval for a mean**
  ([`seq_fwd_dfr_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_dfr_mean.md))
  covers 0.939 for normal data when about 61 observations are needed.
  For exponential data with the same standard deviation it covers 0.854
  at that target and 0.934 when about 246 are needed: skewness slows the
  approach to the stated level considerably.
- **The interval for the location of a negative exponential
  distribution**
  ([`seq_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_fwd_nexp_mean.md))
  covers 0.932 when about 25 observations are needed. Here too the two
  stage interval,
  [`ts_fwd_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/ts_fwd_nexp_mean.md),
  has a guaranteed level.
- **The interval for a median**
  ([`seq_blci_dfr_median()`](https://yelleknek.github.io/SMSD/reference/seq_blci_dfr_median.md))
  has at least its stated coverage at any fixed sample size, but the
  sequential procedure stops when the interval happens to be short. Its
  coverage is 0.917 to 0.932 in the three conditions above, and the
  shortfall does not go away quickly: it is 0.932 in the condition that
  ends with about 762 observations.
- **The bounded risk procedure for a normal mean**
  ([`seq_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_br_norm_mean.md))
  approaches its bound from above. When the bound calls for about 28
  observations under squared error loss the risk is 1.22 times the
  bound, and when it calls for about 27 under absolute error loss it is
  1.09 times the bound. The two stage procedure
  [`ts_br_norm_mean()`](https://yelleknek.github.io/SMSD/reference/ts_br_norm_mean.md)
  keeps its risk at or below the bound for normal data, as the exact
  results above show.
- **The minimum risk procedures** attain the smallest risk only in the
  limit: the risk of
  [`seq_mr_norm_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_norm_mean.md)
  is 1.07 times the smallest attainable when the best sample size is 10,
  and that of
  [`seq_mr_exp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_exp_mean.md)
  is 1.15 times when it is
  15. ([`seq_mr_nexp_mean()`](https://yelleknek.github.io/SMSD/reference/seq_mr_nexp_mean.md)
      and
      [`seq_mr_dfr_mean_arb()`](https://yelleknek.github.io/SMSD/reference/seq_mr_dfr_mean_arb.md)
      are not in the study; their help pages report the same kind of
      figures.)

### What Helps

A larger pilot sample helps where the shortfall comes from stopping
early. A minimum pilot sample size of 40 raised coverage from 0.78 to
0.90 for the correlation of 0.4 at a width of 0.4, and from 0.83 to 0.92
for the slope of 0.5, at a cost of 8 and 9 more pairs on average. It is
not a cure for a strong relationship or for skewed data.

The practical rule follows from the tables. Before relying on an
interval at a planned width, run the matching `*_sensitivity()` function
at each true value that is plausible for the study, the strongest
included, with data as skewed as the study’s are likely to be, and read
its `coverage` row. If coverage is low, plan a narrower width, a larger
minimum pilot sample size, or both. For a classical procedure, which has
no such function, the help page gives its measured coverage or risk over
a range of targets, and where a two stage procedure with an exact
guarantee exists it is named there.

## Reproducing This Article

The study is one script. From the root of the package’s repository, with
the package’s source loaded or the package installed:

``` r
Rscript -e 'devtools::load_all(quiet = TRUE); source("data-raw/validation_study.R")'
```

It writes `data-raw/validation_results.rds`, a data frame with one row
per comparison (the procedure, the source with its table and page, the
condition, the published value and its standard error, the package’s
value and its standard error, the number of replications, and z), with
the seed, the date, and the package version as attributes. The
replication counts are in the tables above; the cached file was made on
2026-10-10 with SMSD 1.0.0 under R version 4.6.1 (2026-06-24). Every
part of every condition sets its own seed, so a second run of the script
gives the same numbers. The script uses 6 cores unless the environment
variable `SMSD_VALIDATION_CORES` says otherwise. On Windows, where R
cannot run the parts as forked processes, it uses one core and takes
several hours. The script’s header describes settings for a shorter
trial run or for some conditions only; this article refuses a cache made
that way.

## References

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods, 30*(2),
425–439. [doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

Chattopadhyay, B., & De, S. K. (2015). *Estimation of Gini index within
pre-specified error bound*. arXiv.
[arXiv:1503.08156v1](https://arxiv.org/abs/1503.08156v1)

Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
pre-specified error bound. *Econometrics, 4*(3), 30.
[doi:10.3390/econometrics4030030](https://doi.org/10.3390/econometrics4030030)

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research, 51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
difference with minimum risk: Maximizing accuracy and minimizing cost
with sequential estimation. *Psychological Methods, 22*(1), 94–113.
[doi:10.1037/met0000089](https://doi.org/10.1037/met0000089)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods, 23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
accuracy in parameter estimation for population correlation
coefficients. *Psychological Methods, 24*(4), 492–515.
[doi:10.1037/met0000203](https://doi.org/10.1037/met0000203)

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. <ISBN:9780367386535>

Stein, C. (1945). A two-sample test for a linear hypothesis whose power
is independent of the variance. *The Annals of Mathematical Statistics,
16*(3), 243–258.
[doi:10.1214/aoms/1177731088](https://doi.org/10.1214/aoms/1177731088)
