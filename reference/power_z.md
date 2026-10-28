# Type II Error and Sample Size for a One Sample, One Sided z Test

For the one sample, upper tailed *z* test of \\H_0\\: \mu = \mu_0\\
against \\H_1\\: \mu = \mu_1 \> \mu_0\\ with a known population standard
deviation, returns either the sample size needed for a chosen Type II
error rate, the Type II error rate at a chosen sample size, or a table
of Type II error rates over a range of sample sizes. Power is one minus
the Type II error rate. The function is named `power_z` so as not to
mask [`stats::power`](https://rdrr.io/r/stats/power.html).

## Usage

``` r
power_z(
  mu_0,
  mu_1,
  sigma,
  conf_level = 0.95,
  n = NULL,
  n_min = NULL,
  n_max = NULL,
  n_values = NULL,
  beta = NULL,
  method = c("sample_size", "type_2", "type_2_seq")
)
```

## Arguments

- mu_0:

  The population mean under the null hypothesis. A single number.

- mu_1:

  The population mean under the alternative hypothesis. A single number
  greater than `mu_0`.

- sigma:

  The known population standard deviation. A single positive number.

- conf_level:

  One minus the one-sided Type I error rate, strictly between 0 and 1.
  Default 0.95, a test at the 0.05 level.

- n:

  For `method = "type_2"`: the sample size at which to compute the Type
  II error rate. A whole number of at least 1.

- n_min, n_max:

  For `method = "type_2_seq"`: the smallest and largest sample sizes of
  the table. Whole numbers with `n_max >= n_min >= 1`.

- n_values:

  For `method = "type_2_seq"`: the number of sample sizes in the table.
  They are equally spaced from `n_min` to `n_max`, and each is taken to
  the nearest whole number (a value halfway between two whole numbers
  goes to the larger); duplicates are dropped, so the table can have
  fewer rows than `n_values`. When `n_values` is at least the number of
  whole sample sizes from `n_min` to `n_max`, the table has every one of
  them.

- beta:

  For `method = "sample_size"`: the target Type II error rate (one minus
  the desired power), strictly between 0 and 1.

- method:

  What to compute: `"sample_size"`, which needs `beta`; `"type_2"`,
  which needs `n`; or `"type_2_seq"`, which needs `n_min`, `n_max`, and
  `n_values`. When `method` is not given, it is the method whose
  arguments are supplied (`beta` alone, `n` alone, or the three grid
  arguments); supplying arguments for more than one result without
  naming `method` is an error. When `method` is given, arguments it does
  not use are ignored. A unique abbreviation of a method name is
  accepted; `NULL` is the same as not giving `method`.

## Value

Depends on `method`:

- `"sample_size"`:

  A length one integer: the smallest sample size whose Type II error
  rate is at most `beta` (\\n\_\beta\\ of Details, rounded up).

- `"type_2"`:

  A length one number: the Type II error rate \\\beta(n)\\.

- `"type_2_seq"`:

  A `data.frame` with one row per sample size and two columns: `n`
  (integer, the sample size) and `type_2_error` (the Type II error rate
  \\\beta(n)\\).

## Details

The test rejects \\H_0\\ when \\\sqrt{n}\\(\bar{X} - \mu_0)/\sigma \>
z\_{1-\alpha}\\, with one-sided Type I error rate \\\alpha = 1 - \\
`conf_level`. Its Type II error rate at \\\mu = \mu_1\\ with \\n\\
observations is \$\$\beta(n) = \Phi\\\left(z\_{1-\alpha} -
\sqrt{n}\\\frac{\mu_1 - \mu_0}{\sigma}\right),\$\$ and the sample size
at which it equals a target \\\beta\\ is \$\$n\_\beta =
\left(\frac{(z\_{1-\alpha} + z\_{1-\beta})\\\sigma} {\mu_1 -
\mu_0}\right)^2.\$\$ \\n\_\beta\\ is rarely a whole number.
`method = "sample_size"` returns the smallest whole sample size whose
Type II error rate \\\beta(n)\\ is at most `beta`, which is \\n\_\beta\\
rounded up. (When \\n\_\beta\\ is within rounding error of a whole
number, \\\beta(n)\\ decides between that number and the next, so asking
for the Type II error rate that this function reports at \\n\\ returns
\\n\\.) When `beta` is at least `conf_level` (a power no greater than
the Type I error rate) every sample size meets the target and 1 is
returned.

Only the upper tailed alternative is computed, so `mu_1` must exceed
`mu_0`. For a lower tailed alternative the results are the same with
both means reflected: pass `-mu_0` and `-mu_1`. For an unknown standard
deviation, two groups, or a two-sided test, use the planning functions
in DMAR, for example
[`ss_power_smd`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.html).

## References

Mukhopadhyay, N., & de Silva, B. M. (2009). *Sequential methods and
their applications*. Chapman & Hall/CRC. ISBN:9780367386535

## See also

[`ss_power_smd`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.html)
in DMAR for two groups with an unknown standard deviation.

## Author

Ken Kelley <kkelley@nd.edu>, Francis Bilson Darku <fbilsond@nd.edu>,
Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>

## Examples

``` r
# A training program is expected to raise the mean score on a
# standardized assessment from 100 to 104; the standard deviation of
# the assessment is known to be 15.

# Sample size for a Type II error rate of 0.20 (power 0.80)
power_z(mu_0 = 100, mu_1 = 104, sigma = 15, beta = 0.20)
#> [1] 87

# Type II error rate with 60 trainees
power_z(mu_0 = 100, mu_1 = 104, sigma = 15, n = 60)
#> [1] 0.3369734

# Type II error rates for 20 to 140 trainees
power_z(mu_0 = 100, mu_1 = 104, sigma = 15, n_min = 20, n_max = 140,
        n_values = 7)
#>     n type_2_error
#> 1  20   0.67446781
#> 2  40   0.48337114
#> 3  60   0.33697339
#> 4  80   0.22956337
#> 5 100   0.15343470
#> 6 120   0.10091885
#> 7 140   0.06547211
```
