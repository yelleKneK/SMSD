# Verification simulation for the standardized mean difference minimum risk
# point estimation procedure, seq_mr_smd(), which implements Chattopadhyay
# and Kelley (2017): the equal cost rule (their Equation 15) and the two cost
# procedure of their section "A General Scenario" (pp. 106 to 107).
#
# It runs the SHIPPED procedure to stopping, through
# seq_mr_smd_sensitivity(), 5,000 times per condition, and
#
#   1. reproduces the mean final sample sizes of the article's Tables 1 and 3
#      (pp. 103 to 104; equal costs, one observation per group per stage);
#   2. checks the two cost procedure, and the equal cost procedure at large
#      effects, against the fixed sample sizes that minimize the article's
#      risk (Equation 17) and against the minimum risk.
#
# This script is not run when the package is built (data-raw is in
# .Rbuildignore). It documents the numbers quoted in the seq_mr_smd() help
# page and the mrpe-smd vignette. It takes about ten minutes on one core.
#
# Run with: Rscript data-raw/smd_mrpe_simulation.R

library(SMSD)

G <- 5000

# ---------------------------------------------------------------------------
# 1. Tables 1 and 3 of Chattopadhyay and Kelley (2017)
# ---------------------------------------------------------------------------
# article_N and article_se are the published mean final sample size per group
# and its standard error.
article <- rbind(
  data.frame(table = 1, A = 50000, c = 1,
             delta = c(0, 0.1, 0.2, 0.4, 0.8, 1.6),
             article_N = c(235.0418, 235.1774, 235.6276, 237.3000, 243.7522,
                           267.8542),
             article_se = c(.0029, .0056, .0093, .0150, .0287, .0535)),
  data.frame(table = 3, A = 1e6, c = 500,
             delta = c(0, 0.1, 0.2, 0.4, 0.8, 1.6),
             article_N = c(50.0210, 50.0392, 50.0964, 50.3882, 51.6882,
                           56.4876),
             article_se = c(.0020, .0028, .0043, .0078, .0129, .0239))
)

tables <- do.call(rbind, lapply(seq_len(nrow(article)), function(i) {
  a <- article[i, ]
  s <- seq_mr_smd_sensitivity(true_value = a$delta, A = a$A, c = a$c,
                              G = G, seed = 2017 + i)
  v <- stats::setNames(s$value, s$term)
  se <- v[["sd_n1"]] / sqrt(G)
  data.frame(table = a$table, delta = a$delta,
             article_N = a$article_N, article_se = a$article_se,
             mean_N = v[["mean_n1"]], se = se,
             z = (v[["mean_n1"]] - a$article_N) / sqrt(se^2 + a$article_se^2),
             optimal_n = v[["optimal_n1"]],
             risk_efficiency = v[["risk_efficiency"]],
             se_risk_efficiency = v[["se_risk_efficiency"]])
}))
print(tables, row.names = FALSE, digits = 7)

# Result of the run recorded on 2026-10-06 (R 4.6.1):
#
#  table delta article_N article_se   mean_N      se      z optimal_n risk_eff
#      1   0.0  235.0418     0.0029 235.0376 0.00273  -1.05   223.607   0.9951
#      1   0.1  235.1774     0.0056 235.1848 0.00574   0.92   223.747   1.0112
#      1   0.2  235.6276     0.0093 235.6274 0.00915  -0.02   224.165   0.9931
#      1   0.4  237.3000     0.0150 237.2942 0.01491  -0.27   225.832   1.0062
#      1   0.8  243.7522     0.0287 243.7242 0.02791  -0.70   232.379   0.9847
#      1   1.6  267.8542     0.0535 267.7942 0.05438  -0.79   256.905   1.0110
#      3   0.0   50.0210     0.0020  50.0218 0.00207   0.28    44.721   1.0204
#      3   0.1   50.0392     0.0028  50.0318 0.00255  -1.96    44.749   1.0095
#      3   0.2   50.0964     0.0043  50.0920 0.00423  -0.73    44.833   1.0183
#      3   0.4   50.3882     0.0078  50.3980 0.00802   0.88    45.166   1.0119
#      3   0.8   51.6882     0.0129  51.6812 0.01274  -0.39    46.476   0.9999
#      3   1.6   56.4876     0.0239  56.4946 0.02354   0.21    51.381   0.9991
#
# Every published mean is reproduced within Monte Carlo error. The standard
# error of each risk efficiency is about 0.009.

# ---------------------------------------------------------------------------
# 2. Two costs, and equal costs at a large effect
# ---------------------------------------------------------------------------
settings <- list(
  list(A = 50000, c = c(2, 8), delta = 0),
  list(A = 50000, c = c(2, 8), delta = 0.5),
  list(A = 50000, c = c(2, 8), delta = 1),
  list(A = 50000, c = c(2, 8), delta = 2),
  list(A = 50000, c = c(2, 8), delta = 3),
  list(A = 50000, c = c(8, 2), delta = 1),
  list(A = 20000, c = c(4, 5), delta = 0.8),
  list(A = 5000,  c = c(1, 4), delta = 0.5),
  list(A = 1e6,   c = c(500, 1000), delta = 0.5),
  list(A = 5000,  c = 1, delta = 3),
  list(A = 5000,  c = 1, delta = 0.5),
  list(A = 50000, c = c(2, 2.0001), delta = 1),
  list(A = 50000, c = 2, delta = 1)
)

runs <- do.call(rbind, lapply(seq_along(settings), function(i) {
  s <- settings[[i]]
  r <- seq_mr_smd_sensitivity(true_value = s$delta, A = s$A, c = s$c,
                              G = G, seed = 500 + i)
  v <- stats::setNames(r$value, r$term)
  data.frame(A = s$A, c1 = v[["cost_1"]], c2 = v[["cost_2"]],
             delta = s$delta,
             optimal_n1 = v[["optimal_n1"]], optimal_n2 = v[["optimal_n2"]],
             mean_n1 = v[["mean_n1"]], mean_n2 = v[["mean_n2"]],
             ratio_n1 = v[["ratio_n1"]], ratio_n2 = v[["ratio_n2"]],
             risk_efficiency = v[["risk_efficiency"]],
             se_risk_efficiency = v[["se_risk_efficiency"]])
}))
print(runs, row.names = FALSE, digits = 5)

# Result of the run recorded on 2026-10-06 (R 4.6.1):
#
#      A  c1     c2 delta optimal_n1 optimal_n2 mean_n1 mean_n2 ratio_n1 ratio_n2 risk_eff
#  50000   2      8   0.0     158.11     79.057  171.56  81.000   1.0850   1.0246   0.9962
#  50000   2      8   0.5     162.51     79.589  176.12  81.548   1.0838   1.0246   1.0097
#  50000   2      8   1.0     175.66     81.002  189.18  82.719   1.0770   1.0212   0.9880
#  50000   2      8   2.0     226.98     84.694  240.08  85.797   1.0577   1.0130   1.0142
#  50000   2      8   3.0     304.85     87.450  316.95  88.152   1.0397   1.0080   0.9986
#  50000   8      2   1.0      81.00    175.656   82.74 189.355   1.0214   1.0780   1.0004
#  20000   4      5   0.8      73.82     65.442   82.69  71.779   1.1202   1.0968   1.0180
#   5000   1      4   0.5      72.67     35.593   81.88  37.007   1.1267   1.0397   1.0168
#  1e+06 500   1000   0.5      45.68     31.956   52.87  34.597   1.1575   1.0826   0.9997
#   5000   1      1   3.0     103.08    103.078  108.77 108.768   1.0552   1.0552   1.0017
#   5000   1      1   0.5      71.81     71.807   78.36  78.360   1.0913   1.0913   1.0164
#  50000   2 2.0001   1.0     167.71    167.701  180.74 180.737   1.0777   1.0777   0.9894
#  50000   2      2   1.0     167.71    167.705  177.25 177.250   1.0569   1.0569   1.0134
#
# The realized risk is within two Monte Carlo standard errors (0.009 to
# 0.010) of the minimum risk in every row. The less costly group is sampled
# more heavily, and both groups end somewhat above their optimal sizes; the
# excess is the article's correction term, sqrt(A / c1) * n^(-gamma) with two
# costs and sqrt(A / (2c)) * n^(-gamma) with one. The last two rows show the
# step between the article's two cost procedure and its one cost procedure
# as the costs approach equality (180.7 against 177.3 per group).
