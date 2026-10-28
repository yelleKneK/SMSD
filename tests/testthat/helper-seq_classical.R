# Helpers for the tests of the purely sequential classical procedures
# (test-seq_fwd.R, test-seq_blci.R, test-seq_br.R, test-seq_mr_classical.R).
#
# The operating characteristic tests need many complete sequential runs.
# Calling a check function after every single observation of every run would
# take far too long for a test suite, so each run is located with running
# sums (the stopping rule coded independently, from its derivation), and the
# package's check function is then asked to confirm the result on the same
# data: met at the stopping point, not met one observation earlier, and, for
# a few runs, not met at earlier sample sizes from the pilot on (every one
# of them for a short run, 100 evenly spaced ones for a long run).

# Rescale x to have exactly the sample variance v (mean unchanged).
seqc_scale_var <- function(x, v) {
  (x - mean(x)) / stats::sd(x) * sqrt(v) + mean(x)
}

# Running mean and variance of x at each n in `n`.
seqc_run_mv <- function(x, n) {
  cs <- cumsum(x)
  cq <- cumsum(x^2)
  xb <- cs[n] / n
  list(mean = xb, var = (cq[n] - n * xb^2) / (n - 1))
}

# Running negative exponential scale estimate U_n = sum(x - min x) / (n - 1).
seqc_run_u <- function(x, n) {
  (cumsum(x)[n] - n * cummin(x)[n]) / (n - 1)
}

# One sequential run: draw observations with rgen(k) and return the first
# n >= m at which ok(x, n) holds, together with the data. `ok` must return a
# logical vector over the integer vector n, using only x[1:n] for each n.
# The first draw is a little longer than the sample size expected (n_guess);
# if the rule is not met within it, the path is extended until it is.
seqc_first_stop <- function(m, rgen, ok, n_guess) {
  x <- rgen(max(m, as.integer(ceiling(n_guess + 8 * sqrt(n_guess)))) + 20L)
  repeat {
    n <- m:length(x)
    j <- which(ok(x, n))[1L]
    if (!is.na(j)) return(list(N = n[j], x = x))
    x <- c(x, rgen(length(x)))
  }
}

# Does the package's check function agree that the run stops at N? It must
# report the rule as met at N and as not met at N - 1 (when N is beyond the
# pilot m). With full = TRUE it must also report not met at the sample sizes
# from m to N - 1: all of them when there are at most 100, otherwise 100
# evenly spaced ones (always including m and N - 1).
seqc_agrees <- function(run, m, check, full = FALSE) {
  x <- run$x
  N <- run$N
  ok <- isTRUE(attr(check(x[seq_len(N)]), "is_satisfied"))
  before <- if (N <= m) {
    integer(0)
  } else if (!full) {
    N - 1L
  } else if (N - m <= 100L) {
    m:(N - 1L)
  } else {
    unique(as.integer(round(seq(m, N - 1L, length.out = 100L))))
  }
  for (n in before) {
    ok <- ok && identical(attr(check(x[seq_len(n)]), "is_satisfied"), FALSE)
  }
  ok
}

# Run `reps` sequential runs. `measure(x, N)` returns the quantity for each run
# whose mean is the operating characteristic (a 0/1 coverage indicator or a
# loss). The package check confirms the first `n_confirm` runs at N and
# N - 1 and the first `n_full` runs at the earlier sample sizes as well. The
# final sample sizes are returned too (element N), for quantities such as
# E(1 / N).
seqc_simulate <- function(reps, m, rgen, ok, n_guess, measure, check,
                          n_confirm = 60L, n_full = 4L) {
  N <- integer(reps)
  y <- numeric(reps)
  agree <- TRUE
  for (i in seq_len(reps)) {
    run <- seqc_first_stop(m, rgen, ok, n_guess)
    N[i] <- run$N
    y[i] <- measure(run$x, run$N)
    if (i <= n_confirm) {
      agree <- agree && seqc_agrees(run, m, check, full = i <= n_full)
    }
  }
  list(mean_n = mean(N), mean_y = mean(y), se_y = stats::sd(y) / sqrt(reps),
       se_n = stats::sd(N) / sqrt(reps), N = N, agree = agree)
}
