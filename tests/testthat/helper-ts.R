# Helpers shared by the tests of the two stage (ts_*) procedures.

# P(N = j) for N = max(m, ceiling(mult * V)) with V chi square on df degrees
# of freedom. Used to compute exact coverage and expected sample size of a
# two stage procedure whose N depends on the first stage only through a
# chi square statistic.
ts_n_dist <- function(m, df, mult) {
  j_max <- max(m, ceiling(mult * stats::qchisq(1e-14, df, lower.tail = FALSE)))
  j <- m:j_max
  cdf <- stats::pchisq(j / mult, df)
  cdf[length(cdf)] <- 1
  list(n = j, p = diff(c(0, cdf)))
}

# One replication of a two stage fixed width procedure run through the
# package: first stage, the total N it fixes, the second stage, and whether
# the final interval covers the true value. (The generator argument is
# named `rgen` so that a `d = ` passed through `...` cannot partially match
# it.)
ts_run_two_stage <- function(check, rgen, m, true_value, ...) {
  first <- rgen(m)
  total <- attr(check(first, ...), "criterion_n")
  out <- check(first, rgen(total - m), ...)
  stopifnot(isTRUE(attr(out, "is_satisfied")),
            attr(out, "current_n") == total)
  v <- stats::setNames(out$value, out$term)
  c(n = total,
    covered = v[["lower_limit"]] <= true_value &&
      true_value <= v[["upper_limit"]])
}

# Exact expected sample size and risk of the two stage bounded risk procedure
# for an exponential mean, N = max(m, ceiling(B Xbar_m^2 / w)) with
# B = factor * 2 m (m + 1) A / ((m - 1)(m - 2)) (Mukhopadhyay and de Silva,
# 2009, Equations 8.3.9, 8.3.11, and 8.3.13 for factor = 0.565). Written
# from the derivation, not from the package code. With theta = 1 (the
# procedure is scale equivariant) the first stage sum S is Gamma(m, 1); given
# S the remaining N - m observations are independent of it, so
#   E[(Xbar_N - 1)^2 | S] = ((S - m)^2 + N - m) / N^2,
# and with n* = A theta^2 / w the risk relative to the bound is
#   n* E[((S - m)^2 + N - m) / N^2],  N = max(m, ceiling(c n* (S / m)^2)),
# c = B / A. N is a step function of S, so each expectation is a sum of
# incomplete gamma integrals.
ts_br_exp_exact <- function(m, nstar, factor = 0.565) {
  cm <- factor * 2 * m * (m + 1) / ((m - 1) * (m - 2))
  s_max <- stats::qgamma(1e-14, m + 2, lower.tail = FALSE)
  j <- m:max(m, ceiling(cm * nstar * (s_max / m)^2))
  upper <- m * sqrt(j / (cm * nstar))      # N <= j exactly when S <= upper
  upper[length(upper)] <- Inf
  part <- function(r) {                    # E[S^r; N = j]
    exp(lgamma(m + r) - lgamma(m)) * diff(c(0, stats::pgamma(upper, m + r)))
  }
  p0 <- part(0); p1 <- part(1); p2 <- part(2)
  c(n_ratio = sum(j * p0) / nstar,
    risk_ratio = nstar * sum((p2 - 2 * m * p1 + m^2 * p0 + (j - m) * p0) / j^2))
}
