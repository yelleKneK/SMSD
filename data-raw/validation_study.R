# Validation study for SMSD: every procedure that has a published row, a
# worked example, or an exact value is run through the package's exported
# functions and compared with it. The result is the cached table that the
# site article vignettes/articles/validation.Rmd displays.
#
# The study has three parts.
#
#   1. Reproductions. Simulation rows and worked examples of the source
#      articles, the optimal allocations of Chattopadhyay, Bandyopadhyay,
#      Kelley, and Padalunkal (2025), the population values and the simulation
#      rows of Chattopadhyay and De (2016; the rows are those of the version
#      its authors posted in 2015), and the exact operating characteristics of
#      the two stage procedures of Mukhopadhyay and de Silva (2009).
#   2. Departures. The places where the package deliberately differs from a
#      source article, each with the evidence for the difference: the
#      unusable variance rule, the variance of the coefficient of variation,
#      Equation 26 of Kelley, Darku, and Chattopadhyay (2018), and the
#      stopping rule for the squared multiple correlation. (Two more, the
#      default allocation rule of the budget procedure and the interval
#      reported for the Gini index, are measured with the reproductions of
#      those procedures.)
#   3. Wide targets. Coverage, or risk, of the procedures whose guarantee is
#      asymptotic, at targets wide enough that the guarantee has not yet
#      taken hold: the intervals for effect sizes, the distribution free two
#      stage interval, and the purely sequential classical procedures.
#
# Every simulated procedure is run with the exported functions only: the
# *_sensitivity() function where it covers the article's condition, and
# otherwise the pilot function and the check function in the collect and check
# loop a study would run (one observation per stage, as in the articles). The
# few values not produced by the package (a rule or a sequence as printed that
# the package deliberately does not implement, a sampling variance, an exact
# or population value) are computed here in base R and labeled as such in the
# column `kind`.
#
# Run from the package root, against the package source or an installed copy:
#
#   Rscript -e 'devtools::load_all(quiet = TRUE); source("data-raw/validation_study.R")'
#   Rscript data-raw/validation_study.R
#
# A simulated value that differs from an exact value by more than three
# standard errors is run a second time with new seeds and ten times the
# replications, and both runs are kept (see "Second runs" near the end).
#
# Optional settings, as environment variables:
#
#   SMSD_VALIDATION_CORES   number of cores (default 6; always 1 on Windows,
#                           where R cannot fork the session)
#   SMSD_VALIDATION_OUT     the file written (default
#                           data-raw/validation_results.rds)
#   SMSD_VALIDATION_RESUME  set to 1 to keep the first runs already in the
#                           output file and repeat only the second runs
#   SMSD_VALIDATION_SCALE   multiplies every replication count (default 1)
#   SMSD_VALIDATION_ONLY    a regular expression: run only the conditions
#                           whose identifier matches it
#
# The last two are for trial runs, and the article refuses a cache made with
# either. Every part of every condition sets its own seed, which depends only
# on the condition and the part, so the results do not depend on the number of
# cores, on which conditions are run, or on whether a run is resumed.
#
# A full run uses about 490 minutes of processor time. The run that made the
# cached file took 94 minutes on 6 cores of a 10 core machine that was busy
# with other work for the first half of it (493 minutes of processor time;
# R 4.6.1, 2026-10-10). The result records the time on the clock and the
# processor time of the run that made it (the attributes minutes and
# processor_minutes, and the same two for the second runs). The script is
# not run when the package is built (data-raw is in .Rbuildignore).

if (!file.exists("DESCRIPTION")) stop("Run from the package root.", call. = FALSE)
if (!"SMSD" %in% loadedNamespaces()) library(SMSD)

scale <- as.numeric(Sys.getenv("SMSD_VALIDATION_SCALE", "1"))
cores <- as.integer(Sys.getenv("SMSD_VALIDATION_CORES", "6"))
# parallel::mclapply() forks the R session, which Windows cannot do: there it
# stops with an error unless it is given one core. The results are the same
# on one core; the study takes several hours.
if (.Platform$OS.type == "windows") cores <- 1L
out_file <- Sys.getenv("SMSD_VALIDATION_OUT", "data-raw/validation_results.rds")
resume <- nzchar(Sys.getenv("SMSD_VALIDATION_RESUME", ""))
only <- Sys.getenv("SMSD_VALIDATION_ONLY", "")
master_seed <- 20261010L
chunks <- 20L        # each condition is run as 20 independent parts
started <- Sys.time()
# Processor time used so far by this session and by the parts it has run (the
# parts are child processes), in minutes. Unlike the time on the clock, it
# does not depend on what else the machine is doing.
processor_minutes <- function() {
  sum(proc.time()[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE) / 60
}
processor_started <- processor_minutes()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
z_of <- function(conf_level) stats::qnorm(1 - (1 - conf_level) / 2)
value_of <- function(check, term) check$value[check$term == term]
binom_se <- function(p, g) sqrt(p * (1 - p) / g)
batch_se <- function(x) stats::sd(x) / sqrt(length(x))
article_reps <- 5000L   # every simulation table reproduced here used 5,000

# One row of the result. `published_digits` is the number of decimals the
# source prints; it is used for the rounding term of the standardized
# difference (see the end of the script).
result_row <- function(family, procedure, source, condition, quantity,
                       published, value, se = NA_real_, replications = NA_integer_,
                       published_se = NA_real_, published_digits = NA_integer_,
                       kind = "published simulation", mean_n = NA_real_,
                       note = "") {
  data.frame(family = family, procedure = procedure, source = source,
             condition = condition, quantity = quantity,
             published = published, published_se = published_se,
             published_digits = as.integer(published_digits),
             value = value, se = se, replications = as.integer(replications),
             kind = kind, mean_n = mean_n, note = note,
             stringsAsFactors = FALSE)
}

# A condition is a list: an identifier, the number of replications, a function
# run(g, seed) that performs g replications, and a function rows(parts) that
# turns the list of its parts into result rows. A condition that is compared
# with an exact value also carries `key`, the procedure and condition labels
# of its rows, by which a second run finds it (see "Second runs" below).
condition <- function(id, G, run, rows, key = NULL) {
  list(id = id, G = G, run = run, rows = rows, key = key)
}

# run() for a *_sensitivity() function: its summary as a named vector.
sensitivity_run <- function(fun, args) {
  force(fun); force(args)
  function(g, seed) {
    s <- do.call(fun, c(args, list(G = g, seed = seed)))
    stats::setNames(s$value, s$term)
  }
}

# Pool the parts of a sensitivity run. Means are means of the part means (the
# parts are of equal size); the standard deviation of the final sample size is
# pooled exactly from the within and between part sums of squares; a summary
# with no standard deviation of its own gets a standard error from the
# variation of the part means.
pool_sd <- function(means, sds, g) {
  G <- sum(g)
  m <- sum(g * means) / G
  sqrt((sum((g - 1) * sds^2) + sum(g * (means - m)^2)) / (G - 1))
}
pool_sensitivity <- function(parts) {
  m <- do.call(rbind, parts)
  g <- m[, "replications"]
  sd_terms <- grep("^sd_n", colnames(m), value = TRUE)
  list(G = sum(g), mean = colMeans(m), batch_se = apply(m, 2, batch_se),
       sd = vapply(sd_terms, function(v) pool_sd(m[, sub("^sd_", "mean_", v)], m[, v], g),
                   numeric(1)),
       se_risk_efficiency = if ("se_risk_efficiency" %in% colnames(m))
         sqrt(mean(m[, "se_risk_efficiency"]^2) / nrow(m)) else NA_real_)
}

# Pool the parts of a hand written loop: one row per replication.
pool_loop <- function(parts) do.call(rbind, parts)

# One replication of a purely sequential procedure through an exported check
# function: start at the pilot sample size, add one observation (or pair) per
# stage, and stop at the first stage at which the check reports the stopping
# rule as met. `check(n)` evaluates the check on the first n observations.
run_to_stopping <- function(check, start, cap) {
  n <- as.integer(start)
  repeat {
    chk <- check(n)
    if (isTRUE(attr(chk, "is_satisfied")) || n >= cap) break
    n <- n + 1L
  }
  c(n = n, lower = value_of(chk, "lower_limit"), upper = value_of(chk, "upper_limit"),
    capped = as.numeric(!isTRUE(attr(chk, "is_satisfied"))))
}

# Summaries of a confidence interval procedure, from a sensitivity function
# or from a loop, and the three rows they give against a published row `pub`:
# a list with N, N_digits, N_se, p, w, and w_se. A published mean is itself
# the mean of 5,000 simulated studies, so it has a Monte Carlo standard error
# of its own. Where the article prints that standard error it is used as
# printed. Where it prints none (N_se or w_se is NA), the standard error is
# taken from the standard deviation found here and the article's 5,000
# replications: for the mean width, the standard error found here times the
# square root of (replications here / 5,000).
interval_summary_sensitivity <- function(parts) {
  s <- pool_sensitivity(parts)
  list(G = s$G, mean_n = s$mean[["mean_n"]], sd_n = s$sd[["sd_n"]],
       coverage = s$mean[["coverage"]], mean_width = s$mean[["mean_width"]],
       se_width = s$batch_se[["mean_width"]],
       pct_unusable = s$mean[["pct_variance_unusable"]],
       pct_capped = s$mean[["pct_hit_max_n"]])
}
interval_summary_loop <- function(parts, truth) {
  m <- pool_loop(parts)
  w <- m[, "upper"] - m[, "lower"]
  list(G = nrow(m), mean_n = mean(m[, "n"]), sd_n = stats::sd(m[, "n"]),
       coverage = mean(m[, "lower"] <= truth & truth <= m[, "upper"]),
       mean_width = mean(w), se_width = batch_se(w),
       pct_unusable = NA_real_, pct_capped = 100 * mean(m[, "capped"]))
}
interval_rows <- function(s, family, procedure, source, condition, pub,
                          kind = "published simulation", note = "") {
  if (s$pct_capped > 0) stop("A replication reached the cap in: ", condition, call. = FALSE)
  n_se <- if (is.na(pub$N_se)) s$sd_n / sqrt(article_reps) else pub$N_se
  w_se <- if (is.na(pub$w_se)) s$se_width * sqrt(s$G / article_reps) else pub$w_se
  one <- function(quantity, published, value, se, published_se, digits) {
    result_row(family, procedure, source, condition, quantity, published, value, se = se,
               replications = s$G, published_se = published_se, published_digits = digits,
               kind = kind, mean_n = s$mean_n, note = note)
  }
  rbind(
    one("mean final sample size", pub$N, s$mean_n, s$sd_n / sqrt(s$G), n_se, pub$N_digits),
    one("coverage", pub$p, s$coverage, binom_se(s$coverage, s$G),
        binom_se(pub$p, article_reps), 4L),
    one("mean interval width", pub$w, s$mean_width, s$se_width, w_se, 4L))
}

# ---------------------------------------------------------------------------
# 1. Kelley, Darku, and Chattopadhyay (2018), Tables 2 to 4 (pp. 233 to 234)
# ---------------------------------------------------------------------------
# Table 2: 95% interval for the standardized mean difference, target width
# 0.2, two normal populations with variance 1. Tables 3 and 4: 90% and 95%
# intervals for the coefficient of variation, target width 0.04, populations
# with a coefficient of variation of 0.2. The article's column s(N) is not
# used as the standard error of its mean final sample size: its entries are
# about the standard deviation of N divided by 5,000 (0.0026 where the
# standard deviation found here is 13.0, for example), not by the square
# root of 5,000. The standard error of each published mean is therefore
# computed from the standard deviation of N found here and the article's
# 5,000 replications. The tables print no standard error for the mean width,
# so its standard error is computed in the same way (see interval_rows()).
family_2018 <- "aipe_2018"
cite_2018 <- "Kelley, Darku, and Chattopadhyay (2018)"

smd_2018 <- list(
  list(delta = 0.3, label = "normal, means 10 and 9.7, variance 1",
       pub = list(N = 778.3, N_digits = 1L, N_se = NA, p = 0.9538, w = 0.1999, w_se = NA)),
  list(delta = 0.5, label = "normal, means 10 and 9.5, variance 1",
       pub = list(N = 793.6, N_digits = 1L, N_se = NA, p = 0.9456, w = 0.1999, w_se = NA)))
conditions_smd_2018 <- lapply(smd_2018, function(k) condition(
  id = paste0("2018_smd_", k$delta), G = 3000L,
  run = sensitivity_run(seq_aipe_smd_sensitivity,
                        list(true_value = k$delta, width = 0.2, conf_level = 0.95,
                             distribution = "normal")),
  rows = function(parts) interval_rows(
    interval_summary_sensitivity(parts), family_2018, "seq_aipe_smd()",
    paste0(cite_2018, ", Table 2, p. 233"), paste0(k$label, "; width 0.2, 95%"), k$pub)))

cv_2018 <- list(
  list(conf = 0.90, table = "Table 3", dist = "normal",
       pub = list(N = 178.8, N_digits = 1L, N_se = NA, p = 0.8844, w = 0.0351, w_se = NA)),
  list(conf = 0.90, table = "Table 3", dist = "gamma",
       pub = list(N = 173.0, N_digits = 1L, N_se = NA, p = 0.8794, w = 0.0347, w_se = NA)),
  list(conf = 0.90, table = "Table 3", dist = "lognormal",
       pub = list(N = 182.3, N_digits = 1L, N_se = NA, p = 0.8650, w = 0.0350, w_se = NA)),
  list(conf = 0.95, table = "Table 4", dist = "normal",
       pub = list(N = 241.2, N_digits = 1L, N_se = NA, p = 0.9422, w = 0.0363, w_se = NA)),
  list(conf = 0.95, table = "Table 4", dist = "gamma",
       pub = list(N = 233, N_digits = 0L, N_se = NA, p = 0.9342, w = 0.0359, w_se = NA)),
  list(conf = 0.95, table = "Table 4", dist = "lognormal",
       pub = list(N = 247.4, N_digits = 1L, N_se = NA, p = 0.9234, w = 0.0362, w_se = NA)))
cv_labels <- c(normal = "normal, mean 10, variance 4",
               gamma = "gamma, shape 25, scale 0.6",
               lognormal = "lognormal, log mean 1, log SD 0.1980")
# seq_aipe_cv_sensitivity() generates gamma or lognormal data, so the normal
# rows are run through the pilot and check functions.
cv_normal_run <- function(conf_level) {
  force(conf_level)
  function(g, seed) {
    set.seed(seed)
    start <- seq_aipe_cv_pilot(width = 0.04, conf_level = conf_level)
    t(replicate(g, {
      x <- stats::rnorm(3000, mean = 10, sd = 2)
      run_to_stopping(function(n) seq_aipe_cv(x[seq_len(n)], width = 0.04,
                                              conf_level = conf_level),
                      start, cap = 3000L)
    }))
  }
}
conditions_cv_2018 <- lapply(cv_2018, function(k) {
  src <- paste0(cite_2018, ", ", k$table, ", p. 234")
  cond <- sprintf("%s; width 0.04, %g%%", cv_labels[[k$dist]], 100 * k$conf)
  if (k$dist == "normal") {
    condition(
      id = paste0("2018_cv_normal_", k$conf), G = 10000L, run = cv_normal_run(k$conf),
      rows = function(parts) interval_rows(
        interval_summary_loop(parts, truth = 0.2), family_2018,
        "seq_aipe_cv_pilot(), seq_aipe_cv()", src, cond, k$pub))
  } else {
    condition(
      id = paste0("2018_cv_", k$dist, "_", k$conf), G = 10000L,
      run = sensitivity_run(seq_aipe_cv_sensitivity,
                            list(true_value = 0.2, width = 0.04, conf_level = k$conf,
                                 distribution = k$dist)),
      rows = function(parts) interval_rows(
        interval_summary_sensitivity(parts), family_2018, "seq_aipe_cv()", src, cond, k$pub))
  }
})

# ---------------------------------------------------------------------------
# 2. Kelley, Darku, and Chattopadhyay (2019): Pearson (Table 2, p. 498),
#    Kendall (Table 8, p. 502), Spearman (Table 10, p. 503)
# ---------------------------------------------------------------------------
# Bivariate normal data, 95% intervals. The article prints the standard error
# of each mean final sample size and of each mean width; they are used as
# printed.
family_2019 <- "correlations_2019"
cite_2019 <- "Kelley, Darku, and Chattopadhyay (2019)"
bivariate_normal <- function(n, rho) {
  x <- stats::rnorm(n)
  list(x = x, y = rho * x + sqrt(1 - rho^2) * stats::rnorm(n))
}

pearson_2019 <- list(
  list(rho = 0.1, width = 0.2, G = 5000L,
       pub = list(N = 372.2, N_digits = 1L, N_se = 0.5746, p = 0.9396, w = 0.1992, w_se = 8.16e-5)),
  list(rho = 0.5, width = 0.2, G = 10000L,
       pub = list(N = 204.0, N_digits = 1L, N_se = 0.7435, p = 0.9012, w = 0.1938, w_se = 3.83e-4)),
  list(rho = 0.5, width = 0.1, G = 4000L,
       pub = list(N = 857.0, N_digits = 1L, N_se = 1.1920, p = 0.9464, w = 0.0998, w_se = 7.91e-6)))
conditions_pearson_2019 <- lapply(pearson_2019, function(k) condition(
  id = sprintf("2019_pearson_%g_%g", k$rho, k$width), G = k$G,
  run = sensitivity_run(seq_aipe_r_sensitivity,
                        list(true_value = k$rho, width = k$width, conf_level = 0.95,
                             distribution = "normal")),
  rows = function(parts) {
    s <- interval_summary_sensitivity(parts)
    interval_rows(s, family_2019, "seq_aipe_r(), Pearson",
                  paste0(cite_2019, ", Table 2, p. 498"),
                  sprintf("bivariate normal, rho = %g; width %g, 95%%", k$rho, k$width), k$pub,
                  note = sprintf("variance estimate unusable at one or more checks in %.2f%% of studies",
                                 s$pct_unusable))
  }))

rank_2019 <- list(
  list(method = "kendall", rho = 0.1, table = "Table 8, p. 502", G = 5000L, cap = 1500L,
       pub = list(N = 170.9, N_digits = 1L, N_se = 0.0371, p = 0.9440, w = 0.1984, w_se = 3.028e-6)),
  list(method = "kendall", rho = 0.5, table = "Table 8, p. 502", G = 10000L, cap = 1500L,
       pub = list(N = 133.2, N_digits = 1L, N_se = 0.1514, p = 0.9318, w = 0.1972, w_se = 1.045e-5)),
  list(method = "spearman", rho = 0.5, table = "Table 10, p. 503", G = 4000L, cap = 2500L,
       pub = list(N = 238.8, N_digits = 1L, N_se = 0.5289, p = 0.9300, w = 0.1986, w_se = 3.026e-5)))
rank_names <- c(kendall = "Kendall", spearman = "Spearman")
rank_truth <- function(method, rho) {
  if (method == "kendall") 2 * asin(rho) / pi else 6 * asin(rho / 2) / pi
}
rank_run <- function(method, rho, width, cap) {
  force(method); force(rho); force(width); force(cap)
  function(g, seed) {
    set.seed(seed)
    start <- seq_aipe_r_pilot(width = width, conf_level = 0.95)
    t(replicate(g, {
      d <- bivariate_normal(cap, rho)
      run_to_stopping(function(n) seq_aipe_r(d$x[seq_len(n)], d$y[seq_len(n)],
                                             width = width, conf_level = 0.95,
                                             method = method),
                      start, cap)
    }))
  }
}
conditions_rank_2019 <- lapply(rank_2019, function(k) condition(
  id = sprintf("2019_%s_%g", k$method, k$rho), G = k$G,
  run = rank_run(k$method, k$rho, 0.2, k$cap),
  rows = function(parts) interval_rows(
    interval_summary_loop(parts, truth = rank_truth(k$method, k$rho)), family_2019,
    sprintf("seq_aipe_r_pilot(), seq_aipe_r(), %s", rank_names[[k$method]]),
    paste0(cite_2019, ", ", k$table),
    sprintf("bivariate normal, rho = %g; width 0.2, 95%%", k$rho), k$pub)))

# ---------------------------------------------------------------------------
# 3. Kelley, Darku, and Chattopadhyay (2019): the squared multiple correlation
#    (Tables 11 to 14, pp. 504 to 507)
# ---------------------------------------------------------------------------
# Target width 0.10, 2 and 5 predictors, population values 0.1, 0.5, and 0.9,
# multivariate normal data. Tables 11 and 12: the log transformed interval at
# 90% and 95%. Tables 13 and 14: the Wald interval at 90% and 95%. The tables
# print the standard error of each mean final sample size, used as printed,
# and none for the mean width (see interval_rows()).
family_r2 <- "r2_2019"
r2_2019 <- data.frame(
  table = rep(c(11L, 12L, 13L, 14L), each = 6L),
  page = rep(c(504L, 505L, 506L, 507L), each = 6L),
  conf = rep(c(0.90, 0.95, 0.90, 0.95), each = 6L),
  ci_method = rep(c("log_transformed", "log_transformed", "wald", "wald"), each = 6L),
  k = rep(rep(c(2L, 5L), each = 3L), 4L),
  P2 = rep(c(0.1, 0.5, 0.9), 8L),
  N = c(297.69, 533.22, 34.84, 364.23, 524.81, 30.90,
        435.52, 760.73, 50.55, 513.15, 754.62, 45.27,
        347.98, 543.16, 68.66, 375.77, 539.22, 67.94,
        495.01, 769.43, 84.17, 522.35, 766.17, 82.94),
  N_se = c(2.0105, 1.1170, 0.3222, 1.3204, 1.4068, 0.3098,
           2.5998, 1.3081, 0.4177, 1.5693, 1.6021, 0.4144,
           1.3031, 0.4694, 0.0963, 1.0788, 0.4790, 0.0812,
           1.5851, 0.5548, 0.1493, 1.3550, 0.5611, 0.1290),
  p = c(0.7404, 0.8858, 0.7688, 0.8802, 0.8716, 0.6668,
        0.8132, 0.9386, 0.8432, 0.9334, 0.9284, 0.7616,
        0.8460, 0.8962, 0.8678, 0.8888, 0.8854, 0.8406,
        0.9152, 0.9522, 0.9222, 0.9442, 0.9450, 0.8988),
  w = c(0.0959, 0.0992, 0.0908, 0.0989, 0.0987, 0.0871,
        0.0968, 0.0995, 0.0936, 0.0993, 0.0991, 0.0910,
        0.0962, 0.0993, 0.0733, 0.0976, 0.0993, 0.0710,
        0.0976, 0.0995, 0.0790, 0.0983, 0.0995, 0.0769),
  stringsAsFactors = FALSE)
r2_names <- c(wald = "Wald", log_transformed = "log transformed")
conditions_r2_2019 <- lapply(seq_len(nrow(r2_2019)), function(i) {
  a <- r2_2019[i, ]
  condition(
    id = sprintf("2019_r2_t%d_k%d_%g", a$table, a$k, a$P2), G = 20000L,
    run = sensitivity_run(seq_aipe_R2_sensitivity,
                          list(true_value = a$P2, width = 0.10, p = a$k,
                               conf_level = a$conf, ci_method = a$ci_method)),
    rows = function(parts) interval_rows(
      interval_summary_sensitivity(parts), family_r2,
      sprintf("seq_aipe_R2(), %s interval", r2_names[[a$ci_method]]),
      sprintf("%s, Table %d, p. %d", cite_2019, a$table, a$page),
      sprintf("%d predictors, population value %g; width 0.10, %g%%", a$k, a$P2, 100 * a$conf),
      list(N = a$N, N_digits = 2L, N_se = a$N_se, p = a$p, w = a$w, w_se = NA)))
})

# ---------------------------------------------------------------------------
# 4. Chattopadhyay and Kelley (2016): Table 1 (p. 642) and the worked example
#    (p. 638)
# ---------------------------------------------------------------------------
# A = 5,000,000 and c = 10. The article's procedure uses its Equation 22 as
# the variance estimate, which seq_mr_cv() offers as variance = "published";
# the default is the asymptotic variance of Kelley, Darku, and Chattopadhyay
# (2018, Equation 30). Both are run. The column n_c of the article is its
# optimal sample size with known parameters.
family_2016 <- "mrpe_cv_2016"
cite_2016 <- "Chattopadhyay and Kelley (2016)"
ck2016 <- data.frame(
  distribution = c("gamma", "gamma", "gamma", "gamma", "lognormal", "normal", "normal", "normal"),
  kappa = c(0.2, 0.4, 0.6, 0.8, 0.4, 0.2, 0.4, 0.6),
  N = c(118.3114, 182.3666, 248.1058, 374.0366, 158.2654, 123.8528, 233.4136, 393.9186),
  N_se = c(0.0660, 0.1598, 0.2592, 0.9041, 0.2056, 0.0870, 0.2601, 0.4454),
  n_c = c(97, 174, 246, 390, 148, 104, 230, 394),
  stringsAsFactors = FALSE)
ck2016_grid <- expand.grid(variance = c("published", "corrected"),
                           row = seq_len(nrow(ck2016)), stringsAsFactors = FALSE)
conditions_2016 <- lapply(seq_len(nrow(ck2016_grid)), function(i) {
  a <- ck2016[ck2016_grid$row[i], ]
  variance <- ck2016_grid$variance[i]
  condition(
    id = sprintf("2016_%s_%g_%s", a$distribution, a$kappa, variance), G = 10000L,
    run = sensitivity_run(seq_mr_cv_sensitivity,
                          list(true_value = a$kappa, A = 5e6, c = 10,
                               distribution = a$distribution, variance = variance)),
    rows = function(parts) {
      s <- pool_sensitivity(parts)
      if (s$mean[["pct_hit_max_n"]] > 0) stop("cap reached", call. = FALSE)
      src <- paste0(cite_2016, ", Table 1, p. 642")
      cond <- sprintf("%s, coefficient of variation %g; A = 5,000,000, c = 10",
                      a$distribution, a$kappa)
      kind <- if (variance == "published") "published simulation" else "deliberate departure"
      out <- result_row(
        family_2016, sprintf("seq_mr_cv(variance = \"%s\")", variance), src, cond,
        "mean final sample size", published = a$N, value = s$mean[["mean_n"]],
        se = s$sd[["sd_n"]] / sqrt(s$G), replications = s$G, published_se = a$N_se,
        published_digits = 4L, kind = kind, mean_n = s$mean[["mean_n"]],
        note = sprintf("risk efficiency %.3f (standard error %.3f)",
                       s$mean[["risk_efficiency"]], s$se_risk_efficiency))
      if (variance == "corrected") {
        # The article prints this sample size rounded up to a whole number and
        # the package's row is not rounded, so the difference has no scale of
        # the kind z assumes: no decimals are recorded and z is NA.
        out <- rbind(out, result_row(
          family_2016, "seq_mr_cv_sensitivity(), row optimal_n", src, cond,
          "optimal sample size with known parameters", published = a$n_c,
          value = s$mean[["optimal_n"]], kind = "deliberate departure",
          note = paste("the package uses the asymptotic variance of the 2018 article,",
                       "Equation 30; the article prints the value rounded up")))
      }
      out
    })
})

# The worked example: 35 observations (number of books read), A = 400,000,
# c = 75. The article prints the risk and the estimate at n = 18, 19, and 35,
# and reports n = 35 as the first sample size at which the rule is met.
ck2016_books <- c(36, 53, 19, 11, 10, 24, 14, 65, 18, 48, 25, 35, 13, 18, 3, 41,
                  5, 3, 44, 26, 13, 39, 2, 3, 26, 22, 8, 15, 12, 22, 5, 21, 23,
                  40, 18)
ck2016_example <- function() {
  src <- paste0(cite_2016, ", worked example, p. 638")
  proc <- "seq_mr_cv(variance = \"published\")"
  printed <- data.frame(n = c(18, 19, 35),
                        risk = c(5964.345, 6224.861, 4891.284),
                        cv = c(0.7391157, 0.7113385, 0.7013904))
  at_n <- do.call(rbind, lapply(seq_len(nrow(printed)), function(i) {
    chk <- seq_mr_cv(ck2016_books[seq_len(printed$n[i])], A = 400000, c = 75,
                     variance = "published")
    cond <- sprintf("the article's data, first %d observations; A = 400,000, c = 75",
                    printed$n[i])
    rbind(
      result_row(family_2016, proc, src, cond, "estimated risk", printed$risk[i],
                 value_of(chk, "risk"), published_digits = 3L, kind = "worked example"),
      result_row(family_2016, proc, src, cond, "coefficient of variation", printed$cv[i],
                 value_of(chk, "cv"), published_digits = 7L, kind = "worked example"))
  }))
  path <- seq_path(seq_mr_cv, ck2016_books, A = 400000, c = 75, variance = "published")
  rbind(at_n, result_row(
    family_2016, "seq_path(seq_mr_cv, variance = \"published\")", src,
    "the article's 35 observations; A = 400,000, c = 75",
    "first sample size at which the rule is met", 35, attr(path, "stop_n"),
    published_digits = 0L, kind = "worked example"))
}

# ---------------------------------------------------------------------------
# 5. Chattopadhyay and Kelley (2017): Tables 1 and 3 (pp. 103 to 104) and the
#    worked example (p. 106)
# ---------------------------------------------------------------------------
family_2017 <- "mrpe_smd_2017"
cite_2017 <- "Chattopadhyay and Kelley (2017)"
ck2017 <- rbind(
  data.frame(table = "Table 1, p. 103", A = 50000, c = 1,
             delta = c(0, 0.1, 0.2, 0.4, 0.8, 1.6),
             N = c(235.0418, 235.1774, 235.6276, 237.3000, 243.7522, 267.8542),
             N_se = c(0.0029, 0.0056, 0.0093, 0.0150, 0.0287, 0.0535),
             stringsAsFactors = FALSE),
  data.frame(table = "Table 3, p. 104", A = 1e6, c = 500,
             delta = c(0, 0.1, 0.2, 0.4, 0.8, 1.6),
             N = c(50.0210, 50.0392, 50.0964, 50.3882, 51.6882, 56.4876),
             N_se = c(0.0020, 0.0028, 0.0043, 0.0078, 0.0129, 0.0239),
             stringsAsFactors = FALSE))
conditions_2017 <- lapply(seq_len(nrow(ck2017)), function(i) {
  a <- ck2017[i, ]
  condition(
    id = sprintf("2017_%g_%g", a$A, a$delta), G = 10000L,
    run = sensitivity_run(seq_mr_smd_sensitivity,
                          list(true_value = a$delta, A = a$A, c = a$c,
                               distribution = "normal")),
    rows = function(parts) {
      s <- pool_sensitivity(parts)
      if (s$mean[["pct_hit_max_n"]] > 0) stop("cap reached", call. = FALSE)
      result_row(family_2017, "seq_mr_smd()", paste0(cite_2017, ", ", a$table),
                 sprintf("normal, delta = %g; A = %s, c = %g", a$delta,
                         format(a$A, big.mark = ",", scientific = FALSE), a$c),
                 "mean final sample size per group", published = a$N,
                 value = s$mean[["mean_n1"]], se = s$sd[["sd_n1"]] / sqrt(s$G),
                 replications = s$G, published_se = a$N_se, published_digits = 4L,
                 mean_n = s$mean[["mean_n1"]],
                 note = sprintf("risk efficiency %.3f (standard error %.3f)",
                                s$mean[["risk_efficiency"]], s$se_risk_efficiency))
    })
})

# The worked example: with d = 1.00, A = 10,000, and c = 2.40 the rule is met
# at n = 75 per group and at no smaller sample size. Two groups whose sample
# standardized mean difference is exactly 1 are built at each n.
ck2017_example <- function() {
  sizes <- 13:80
  met <- vapply(sizes, function(n) {
    q <- stats::qnorm(stats::ppoints(n))
    q <- (q - mean(q)) / stats::sd(q)
    isTRUE(attr(seq_mr_smd(q + 1, q, A = 10000, c = 2.40), "is_satisfied"))
  }, logical(1))
  result_row(family_2017, "seq_mr_smd()", paste0(cite_2017, ", worked example, p. 106"),
             "sample standardized mean difference 1.00 at every n; A = 10,000, c = 2.40",
             "first sample size per group at which the rule is met", 75,
             sizes[which(met)[1L]], published_digits = 0L, kind = "worked example")
}

# ---------------------------------------------------------------------------
# 6. Chattopadhyay, Bandyopadhyay, Kelley, and Padalunkal (2025): the example
#    of p. 430, Table 1 (p. 431), and the article's sequence step for step
# ---------------------------------------------------------------------------
family_2025 <- "budget_2025"
cite_2025 <- "Chattopadhyay, Bandyopadhyay, Kelley, and Padalunkal (2025)"
mixture_sd <- function(m1, s1, m2, s2, p = 0.95) {
  mu <- p * m1 + (1 - p) * m2
  sqrt(p * (s1^2 + m1^2) + (1 - p) * (s2^2 + m2^2) - mu^2)
}
table1_2025 <- list(
  list(label = "gamma populations", cost = c(10, 12, 36),
       sd = c(0.2 * sqrt(16), 0.125 * sqrt(50), 1.50001 * sqrt(6)),
       n = c(41.2035, 83.1149, 99.7385), n_digits = c(4L, 4L, 4L), V = 0.1885, note = ""),
  list(label = "logistic populations", cost = c(16, 14, 39), sd = c(2, 1, 7) * pi / sqrt(3),
       n = c(42.2309, 45.1468, 94.6731), n_digits = c(4L, 4L, 4L), V = 2.3058, note = ""),
  list(label = "normal populations", cost = c(38, 15, 10), sd = c(4.5, 0.9, 1),
       n = c(96.373, 61.3566, 41.7479), n_digits = c(3L, 4L, 4L), V = 0.2869, note = ""),
  list(label = "normal mixture populations", cost = c(38, 15, 10),
       sd = c(mixture_sd(1, 19, 20, 15), mixture_sd(2, 5, 20, 8), mixture_sd(3, 11, 20, 11)),
       n = c(75.913, 81.5976, 89.1344), n_digits = c(3L, 4L, 4L), V = 0.8486,
       note = "the printed 0.8486 appears to be a misprint: the article's formula gives 8.4788"))
budget_closed_forms <- function() {
  ex <- seq_budget_optimal(contrast_weights = c(1, 1, 1), sigmas = c(1, 1, 1),
                           cost_per_group = c(30, 30, 30), budget = 10000)
  example <- result_row(
    family_2025, "seq_budget_optimal()", paste0(cite_2025, ", example, p. 430"),
    "three groups, 30 per observation, budget 10,000, equal standard deviations, weights of equal size",
    "optimal sample size per group", 111.11, ex$n_per_group[1L],
    published_digits = 2L, kind = "closed form")
  src <- paste0(cite_2025, ", Table 1, p. 431")
  table1 <- do.call(rbind, lapply(table1_2025, function(b) {
    o <- seq_budget_optimal(contrast_weights = c(1, -2, 1), sigmas = b$sd,
                            cost_per_group = b$cost, budget = 5000)
    cond <- paste0(b$label, ", weights 1, -2, 1; budget 5,000")
    rbind(
      do.call(rbind, lapply(1:3, function(j) result_row(
        family_2025, "seq_budget_optimal()", src, cond,
        sprintf("optimal sample size, group %d", j), b$n[j], o$n_per_group[j],
        published_digits = b$n_digits[j], kind = "closed form"))),
      result_row(family_2025, "seq_budget_optimal()", src, cond,
                 "minimum variance of the contrast", b$V, o$minimum_variance,
                 published_digits = 4L, kind = "closed form", note = b$note))
  }))
  rbind(example, table1)
}
# The normal row of Table 1 by simulation, with the article's pilot sample of
# 10 per group (p. 432). procedure = "first_passage" is the package's reading
# of the article's rule; "reallocate" is the package's default. The standard
# errors of the article's means come from what its table prints: a standard
# deviation for each final sample size; and, for the estimated variance, a
# mean absolute deviation (0.0285), which is turned into a standard deviation
# by the relation that holds for a normal variable (the mean absolute
# deviation times the square root of pi / 2). The table prints no mean cost.
# The cost its mean final sample sizes imply, the sum over groups of cost
# times mean final sample size, is compared with the mean cost found here;
# its standard error is taken from the standard deviation of the cost found
# here and the article's 5,000 replications.
budget_design <- list(means = c(1, 2, 3), sds = c(4.5, 0.9, 1), weights = c(1, -2, 1),
                      cost = c(38, 15, 10), budget = 5000, pilot = 10L)
budget_pub <- list(N = c(95.868, 61.3886, 41.3042), N_digits = c(3L, 4L, 4L),
                   sd_N = c(6.7458, 13.4775, 10.8187), V = 0.2819, V_mad = 0.0285)
budget_condition <- "normal populations, weights 1, -2, 1; budget 5,000, pilot 10 per group"
# The rows of one procedure against the table. `n`, `v`, and `cost` are lists
# with the mean and its standard error (`n` for the three groups).
budget_rows <- function(proc, kind, G, n, v, cost, note = "") {
  src <- paste0(cite_2025, ", Table 1, p. 431")
  rbind(
    do.call(rbind, lapply(1:3, function(j) result_row(
      family_2025, proc, src, budget_condition, sprintf("mean final sample size, group %d", j),
      budget_pub$N[j], n$mean[j], se = n$se[j], replications = G,
      published_se = budget_pub$sd_N[j] / sqrt(article_reps),
      published_digits = budget_pub$N_digits[j], kind = kind, note = note))),
    result_row(family_2025, proc, src, budget_condition,
               "mean estimated variance of the contrast", budget_pub$V, v$mean, se = v$se,
               replications = G,
               published_se = budget_pub$V_mad * sqrt(pi / 2) / sqrt(article_reps),
               published_digits = 4L, kind = kind, note = note),
    result_row(family_2025, proc,
               paste0(src, " (implied by its mean final sample sizes)"), budget_condition,
               "mean total cost", sum(budget_design$cost * budget_pub$N), cost$mean,
               se = cost$se, replications = G,
               published_se = cost$se * sqrt(G / article_reps), kind = kind, note = note))
}
conditions_2025 <- lapply(c("first_passage", "reallocate"), function(procedure) condition(
  id = paste0("2025_normal_", procedure), G = 10000L,
  run = sensitivity_run(seq_budget_sensitivity,
                        list(true_means = budget_design$means, true_sds = budget_design$sds,
                             contrast_weights = budget_design$weights,
                             cost_per_group = budget_design$cost,
                             budget = budget_design$budget, conf_level = 0.95,
                             distribution = "normal", m_0 = budget_design$pilot,
                             procedure = procedure)),
  rows = function(parts) {
    s <- pool_sensitivity(parts)
    proc <- sprintf("seq_budget(), seq_budget_next(); procedure = \"%s\"", procedure)
    kind <- if (procedure == "first_passage") "published simulation" else "deliberate departure"
    rbind(
      budget_rows(proc, kind, s$G,
                  n = list(mean = unname(s$mean[paste0("mean_n_", 1:3)]),
                           se = unname(s$batch_se[paste0("mean_n_", 1:3)])),
                  v = list(mean = s$mean[["mean_estimated_variance"]],
                           se = s$batch_se[["mean_estimated_variance"]]),
                  cost = list(mean = s$mean[["mean_total_cost"]],
                              se = s$batch_se[["mean_total_cost"]])),
      # Not a mean with sampling error: the article reports that the cost
      # never exceeded the budget, and the package cannot exceed it. z is NA.
      result_row(family_2025, proc, paste0(cite_2025, ", p. 432"), budget_condition,
                 "percent of studies within the budget", 100, s$mean[["pct_within_budget"]],
                 replications = s$G, kind = kind,
                 note = sprintf("the rule was met in every group in %.1f%% of studies",
                                s$mean[["pct_rule_met"]])),
      result_row(family_2025, proc, "nominal level", budget_condition, "coverage", 0.95,
                 s$mean[["coverage"]], se = binom_se(s$mean[["coverage"]], s$G),
                 replications = s$G, kind = "nominal level"))
  }))

# The article's sequence itself, written out step for step in base R with no
# package code, for the same design (pp. 429 to 430 and Figure 1). After the
# pilot stage a population whose sample size meets its criterion is closed
# (Equation 14). The open populations are then visited in rounds, in the
# order of their first stage standard deviations, the largest first. One
# observation is drawn from the population visited, and that population
# alone is judged against its criterion, with the standard deviations current
# at that moment (Equations 16 to 18); once it meets the criterion it is
# closed for good. Sampling ends when every population is closed or at the
# first observation that the budget cannot pay for. The package's first
# passage option differs from this sequence in three details, which the help
# page of seq_budget_next() lists.
budget_literal_run <- function(g, seed) {
  set.seed(seed)
  b <- budget_design
  target <- function(i, s) {
    b$budget * abs(b$weights[i]) * s[i] / sum(abs(b$weights) * s * sqrt(b$cost * b$cost[i]))
  }
  t(replicate(g, {
    x <- lapply(1:3, function(i) stats::rnorm(b$pilot, b$means[i], b$sds[i]))
    n <- rep(b$pilot, 3L)
    s <- vapply(x, stats::sd, numeric(1))
    open <- vapply(1:3, function(i) n[i] < target(i, s), logical(1))
    visit <- order(-s)
    affordable <- TRUE
    while (any(open) && affordable) {
      for (i in visit[open[visit]]) {
        if (sum(b$cost * n) + b$cost[i] > b$budget) {
          affordable <- FALSE
          break
        }
        x[[i]] <- c(x[[i]], stats::rnorm(1L, b$means[i], b$sds[i]))
        n[i] <- n[i] + 1L
        s[i] <- stats::sd(x[[i]])
        if (n[i] >= target(i, s)) open[i] <- FALSE
      }
    }
    c(n1 = n[1L], n2 = n[2L], n3 = n[3L], cost = sum(b$cost * n),
      variance = sum(b$weights^2 * s^2 / n), closed = as.numeric(!any(open)))
  }))
}
conditions_2025_literal <- list(condition(
  id = "2025_normal_article_sequence", G = 20000L, run = budget_literal_run,
  rows = function(parts) {
    m <- pool_loop(parts)
    stopifnot(all(m[, "cost"] <= budget_design$budget))
    budget_rows("the article's sequence, step for step (base R, not the package's code)",
                "rule as printed", nrow(m),
                n = list(mean = unname(colMeans(m[, c("n1", "n2", "n3")])),
                         se = unname(apply(m[, c("n1", "n2", "n3")], 2, batch_se))),
                v = list(mean = mean(m[, "variance"]), se = batch_se(m[, "variance"])),
                cost = list(mean = mean(m[, "cost"]), se = batch_se(m[, "cost"])),
                note = sprintf("every population was closed in %.1f%% of studies",
                               100 * mean(m[, "closed"])))
  }))

# ---------------------------------------------------------------------------
# 7. Chattopadhyay and De (2016): population values of the Gini index and of
#    its asymptotic variance, and the sequential procedure against the
#    simulation table of the version posted in 2015
# ---------------------------------------------------------------------------
# Two populations of the article's simulation study: gamma with shape 2.649
# and rate 0.84, and lognormal with log mean 2.185 and log standard deviation
# 0.562. Their Gini index, G = Delta / (2 mu) with Delta = E|X_1 - X_2|, and
# the asymptotic variance of its estimator (the article's Equation 6),
#
#   xi^2 = Delta^2 sigma^2 / (4 mu^4) - Delta tau / mu^3 + Delta^2 / mu^2
#          + sigma_1^2 / mu^2,
#
# with tau = E(X_1 |X_1 - X_2|) and sigma_1^2 the variance of
# h(X_1) = E(|X_1 - X_2| given X_1), are computed here by numerical
# integration with no package code, from
#
#   h(x) = x (2 F(x) - 1) + mu - 2 P(x),  P(x) = integral from 0 to x of t f(t) dt,
#
# so that Delta = E h(X), tau = E[X h(X)], and sigma_1^2 = E[h(X)^2] - Delta^2.
# gini_index() and the variance estimate of seq_aipe_gini() (the row
# xi_squared) are averaged over samples of 2,000 and compared with them.
family_gini <- "gini_2016"
cite_gini <- "Chattopadhyay and De (2016)"
gini_population_values <- function(density, cdf, partial_mean, mu, sigma2, lower = 0) {
  h <- function(x) x * (2 * cdf(x) - 1) + mu - 2 * partial_mean(x)
  expect <- function(f) {
    stats::integrate(function(x) f(x) * density(x), lower, Inf, rel.tol = 1e-10,
                     subdivisions = 2000L)$value
  }
  delta <- expect(h)
  tau <- expect(function(x) x * h(x))
  sigma1_2 <- expect(function(x) h(x)^2) - delta^2
  c(gini = delta / (2 * mu),
    xi2 = delta^2 * sigma2 / (4 * mu^4) - delta * tau / mu^3 + delta^2 / mu^2 +
      sigma1_2 / mu^2)
}
gini_populations <- local({
  a <- 2.649; r <- 0.84              # gamma: shape and rate
  ml <- 2.185; sl <- 0.562           # lognormal: log mean and log standard deviation
  mu_l <- exp(ml + sl^2 / 2)
  list(
    list(id = "gamma", label = "gamma, shape 2.649, rate 0.84",
         draw = function(n) stats::rgamma(n, shape = a, rate = r),
         population = gini_population_values(
           function(x) stats::dgamma(x, a, r), function(x) stats::pgamma(x, a, r),
           function(x) (a / r) * stats::pgamma(x, a + 1, r), a / r, a / r^2)),
    list(id = "lognormal", label = "lognormal, log mean 2.185, log SD 0.562",
         draw = function(n) stats::rlnorm(n, ml, sl),
         population = gini_population_values(
           function(x) stats::dlnorm(x, ml, sl), function(x) stats::plnorm(x, ml, sl),
           function(x) mu_l * stats::pnorm((log(x) - ml - sl^2) / sl), mu_l,
           (exp(sl^2) - 1) * mu_l^2)))
})
gini_population_condition <- function(k) condition(
  id = paste0("gini_", k$id), G = 2000L,
  run = function(g, seed) {
    set.seed(seed)
    t(replicate(g, {
      x <- k$draw(2000)
      c(gini = gini_index(x), xi2 = value_of(seq_aipe_gini(x, width = 0.04), "xi_squared"))
    }))
  },
  rows = function(parts) {
    m <- pool_loop(parts)
    src <- paste0(cite_gini, ", Equation 6, by numerical integration")
    cond <- paste0(k$label, "; samples of 2,000")
    rbind(
      result_row(family_gini, "gini_index()", src, cond, "Gini index, mean estimate",
                 k$population[["gini"]], mean(m[, "gini"]), se = batch_se(m[, "gini"]),
                 replications = nrow(m), kind = "population value"),
      result_row(family_gini, "seq_aipe_gini(), row xi_squared", src, cond,
                 "asymptotic variance, mean estimate", k$population[["xi2"]],
                 mean(m[, "xi2"]), se = batch_se(m[, "xi2"]), replications = nrow(m),
                 kind = "population value"))
  })
conditions_gini <- lapply(gini_populations, gini_population_condition)

# The sequential procedure itself. The version of the article that its
# authors posted in 2015 (arXiv:1503.08156v1) reports, in its Table 1
# (p. 19), 2,000 simulated studies of a 90% interval with a half width of
# 0.01 for each of these populations and for a Pareto population (scale
# 20,000, shape 5), with the mean final sample size, its standard error, the
# largest final sample size, and the coverage with its standard error. The
# package's target is the full width, 0.02, and its pilot sample size for
# that target is the 165 of the posted version. Coverage is counted twice in
# the same studies. The article's interval is the estimate plus and minus
# 0.01, and its coverage at the package's stopping point is what compares
# with the posted table. seq_aipe_gini() reports another interval, the
# estimate plus and minus z times its estimated standard error, which is
# never wider than the target and so never covers more often; its coverage is
# recorded as a deliberate departure.
cite_gini_posted <- "Chattopadhyay and De (2015, arXiv:1503.08156v1)"
gini_posted <- list(
  gamma = list(N = 1259.492, N_digits = 3L, N_se = 4.3639, p = 0.878, p_digits = 3L,
               p_se = 0.0073, max_n = 1594),
  lognormal = list(N = 1429.349, N_digits = 3L, N_se = 4.1393, p = 0.9015, p_digits = 4L,
                   p_se = 0.0067, max_n = 2391),
  pareto = list(N = 654.5364, N_digits = 4L, N_se = 4.2151, p = 0.9018, p_digits = 4L,
                p_se = 0.0063, max_n = 1666))
gini_procedure_condition <- function(k) condition(
  id = paste0("gini_procedure_", k$id), G = 2000L,
  run = function(g, seed) {
    set.seed(seed)
    cap <- 8000L
    start <- seq_aipe_gini_pilot(width = 0.02, conf_level = 0.90)
    stopifnot(start == 165L)
    t(replicate(g, {
      x <- k$draw(cap)
      run_to_stopping(function(n) seq_aipe_gini(x[seq_len(n)], width = 0.02,
                                                conf_level = 0.90),
                      start, cap)
    }))
  },
  rows = function(parts) {
    m <- pool_loop(parts)
    if (any(m[, "capped"] > 0)) stop("A replication reached the cap in: ", k$id, call. = FALSE)
    pub <- gini_posted[[k$id]]
    truth <- k$population[["gini"]]
    middle <- (m[, "lower"] + m[, "upper"]) / 2
    width <- m[, "upper"] - m[, "lower"]
    stopifnot(all(width <= 0.02 * (1 + 1e-12)))
    in_article <- mean(abs(middle - truth) <= 0.01)
    in_package <- mean(m[, "lower"] <= truth & truth <= m[, "upper"])
    src <- paste0(cite_gini_posted, ", Table 1, p. 19")
    proc <- "seq_aipe_gini_pilot(), seq_aipe_gini()"
    cond <- paste0(k$label, "; width 0.02 (half width 0.01), 90%")
    rbind(
      result_row(family_gini, proc, src, cond, "mean final sample size", pub$N,
                 mean(m[, "n"]), se = batch_se(m[, "n"]), replications = nrow(m),
                 published_se = pub$N_se, published_digits = pub$N_digits,
                 mean_n = mean(m[, "n"]),
                 note = sprintf("final sample sizes from %d to %d (largest in the posted version: %d)",
                                min(m[, "n"]), max(m[, "n"]), pub$max_n)),
      result_row(family_gini, proc, src, cond,
                 "coverage of the article's interval (the estimate plus and minus 0.01)",
                 pub$p, in_article, se = binom_se(in_article, nrow(m)),
                 replications = nrow(m), published_se = pub$p_se,
                 published_digits = pub$p_digits, mean_n = mean(m[, "n"])),
      result_row(family_gini, proc, src, cond, "coverage of the interval the package reports",
                 pub$p, in_package, se = binom_se(in_package, nrow(m)),
                 replications = nrow(m), published_se = pub$p_se,
                 published_digits = pub$p_digits, kind = "deliberate departure",
                 mean_n = mean(m[, "n"]),
                 note = sprintf("mean width %.4f, never more than 0.02", mean(width))))
  })
conditions_gini_procedure <- lapply(gini_populations, gini_procedure_condition)

# The Pareto population of the posted table: scale 20,000 and shape 5. The
# Gini index and its asymptotic variance do not depend on the scale, so the
# population values are computed for scale 1, where the density is
# 5 x^-6 for x > 1, the mean is 5 / 4, and the variance is 5 / 48. The Gini
# index is 1 / 9, and the asymptotic variance gives the optimal fixed sample
# size that the posted table prints for this population, 686.
gini_pareto <- local({
  a <- 5; mu <- a / (a - 1)
  list(id = "pareto", label = "Pareto, scale 20,000, shape 5",
       draw = function(n) 20000 * stats::runif(n)^(-1 / a),
       population = gini_population_values(
         function(x) a * x^(-a - 1), function(x) 1 - x^(-a),
         function(x) mu * (1 - x^(1 - a)), mu, a / ((a - 1)^2 * (a - 2)), lower = 1))
})
stopifnot(floor(z_of(0.90)^2 * gini_pareto$population[["xi2"]] / 0.01^2) + 1 == 686)
conditions_gini_pareto <- list(gini_population_condition(gini_pareto),
                               gini_procedure_condition(gini_pareto))

# ---------------------------------------------------------------------------
# 8. The two stage procedures of Mukhopadhyay and de Silva (2009), against
#    their exact operating characteristics
# ---------------------------------------------------------------------------
# For each procedure the total sample size N depends on the first stage
# through one statistic with a known distribution, so the expected N and the
# coverage or risk can be computed exactly by summing over the distribution
# of N. These functions call no package code.
#
# Stein (normal mean): N = max(m, ceiling(t^2 S^2 / d^2)), S^2 independent of
#   every sample mean; given N the mean of all N observations is normal, so
#   the coverage is E[2 Phi(d sqrt(N) / sigma) - 1]. The same with a normal
#   quantile in place of t is the distribution free rule applied to normal data.
# Negative exponential location: N = max(m, ceiling(F U / d)), with
#   2 (m - 1) U / sigma chi square on 2m - 2 degrees of freedom; given N the
#   interval misses with probability exp(-N d / sigma).
# Bounded risk, normal mean: N = max(m, ceiling((b B / w)^(2/k) S^2)); the
#   risk relative to w is E[(n* / N)^(k/2)].
# Bounded risk, exponential mean: N = max(m, ceiling(B Xbar_m^2 / w)); given
#   the first stage sum S (gamma with shape m when theta = 1) the squared
#   error of the final mean has expectation ((S - m)^2 + N - m) / N^2.
exact_n_chisq <- function(m, df, mult) {
  top <- max(m, ceiling(mult * stats::qchisq(1e-15, df, lower.tail = FALSE)))
  j <- m:top
  cdf <- stats::pchisq(j / mult, df)
  cdf[length(cdf)] <- 1
  list(j = j, p = diff(c(0, cdf)))
}
exact_stein <- function(m, C, conf_level = 0.95, quantile = c("t", "z")) {
  quantile <- match.arg(quantile)
  z <- z_of(conf_level)
  q <- if (quantile == "t") stats::qt(1 - (1 - conf_level) / 2, m - 1) else z
  pn <- exact_n_chisq(m, m - 1, q^2 * C / (z^2 * (m - 1)))
  c(mean_n = sum(pn$j * pn$p),
    second = sum(pn$p * (2 * stats::pnorm(z * sqrt(pn$j / C)) - 1)))
}
exact_nexp <- function(m, C, conf_level = 0.95) {
  alpha <- 1 - conf_level
  a <- -log(alpha)
  f <- stats::qf(1 - alpha, 2, 2 * m - 2)
  pn <- exact_n_chisq(m, 2 * m - 2, f * C / (a * 2 * (m - 1)))
  c(mean_n = sum(pn$j * pn$p), second = sum(pn$p * (1 - exp(-pn$j * a / C))))
}
exact_br_norm <- function(m, n_star, k) {
  b <- exp((k / 2) * log((m - 1) / 2) + lgamma((m - k - 1) / 2) - lgamma((m - 1) / 2))
  pn <- exact_n_chisq(m, m - 1, b^(2 / k) * n_star / (m - 1))
  c(mean_n = sum(pn$j * pn$p), second = sum(pn$p * (n_star / pn$j)^(k / 2)))
}
exact_br_exp <- function(m, n_star, factor = 0.565) {
  cm <- factor * 2 * m * (m + 1) / ((m - 1) * (m - 2))
  top <- max(m, ceiling(cm * n_star *
                          (stats::qgamma(1e-15, m + 2, lower.tail = FALSE) / m)^2))
  j <- m:top
  upper <- m * sqrt(j / (cm * n_star))
  upper[length(upper)] <- Inf
  part <- function(r) exp(lgamma(m + r) - lgamma(m)) * diff(c(0, stats::pgamma(upper, m + r)))
  i0 <- part(0); i1 <- part(1); i2 <- part(2)
  c(mean_n = sum(j * i0),
    second = n_star * sum((i2 - 2 * m * i1 + m^2 * i0 + (j - m) * i0) / j^2))
}

family_classical <- "classical"
cite_book <- "Mukhopadhyay and de Silva (2009)"
# `m` has one row per study: the total sample size n and either ok (the
# interval covered) or loss (the realized loss divided by the risk bound).
two_stage_rows <- function(m, procedure, source, condition, exact, guarantee) {
  coverage <- "ok" %in% colnames(m)
  second <- if (coverage) mean(m[, "ok"]) else mean(m[, "loss"])
  rbind(
    result_row(family_classical, procedure, source, condition, "mean total sample size",
               exact[["mean_n"]], mean(m[, "n"]), se = batch_se(m[, "n"]),
               replications = nrow(m), kind = "exact calculation", mean_n = mean(m[, "n"])),
    result_row(family_classical, procedure, source, condition,
               if (coverage) "coverage" else "risk relative to the bound",
               exact[["second"]], second,
               se = if (coverage) binom_se(second, nrow(m)) else batch_se(m[, "loss"]),
               replications = nrow(m), kind = "exact calculation", mean_n = mean(m[, "n"]),
               note = guarantee))
}
second_stage <- function(n, draw) if (n > 0) draw(n) else NULL
covers <- function(check, truth) {
  as.numeric(value_of(check, "lower_limit") <= truth && truth <= value_of(check, "upper_limit"))
}

# Stein's two stage interval for a normal mean (mean 50, standard deviation 2).
conditions_stein <- lapply(
  list(list(m = 5L, d = 1), list(m = 10L, d = 0.5), list(m = 30L, d = 0.2)),
  function(k) {
    label <- sprintf("normal, standard deviation 2; half width %g, 95%%, first stage of %d",
                     k$d, k$m)
    condition(
      id = sprintf("ts_stein_%d", k$m), G = 100000L,
      run = function(g, seed) {
        set.seed(seed)
        draw <- function(r) stats::rnorm(r, 50, 2)
        t(replicate(g, {
          pilot <- draw(k$m)
          n <- attr(ts_fwd_norm_mean(pilot, d = k$d), "criterion_n")
          c(n = n,
            ok = covers(ts_fwd_norm_mean(pilot, second_stage(n - k$m, draw), d = k$d), 50))
        }))
      },
      rows = function(parts) two_stage_rows(
        pool_loop(parts), "ts_fwd_norm_mean()",
        paste0("exact calculation; Stein's procedure, in ", cite_book), label,
        exact_stein(k$m, z_of(0.95)^2 * 4 / k$d^2), "guaranteed: at least 0.95"),
      key = c(procedure = "ts_fwd_norm_mean()", condition = label))
  })

# The two stage interval for the location of a negative exponential
# distribution (location 5, scale 2).
conditions_nexp <- lapply(
  list(list(m = 10L, d = 0.2), list(m = 20L, d = 0.05)),
  function(k) {
    label <- sprintf("negative exponential, scale 2; width %g, 95%%, first stage of %d",
                     k$d, k$m)
    condition(
      id = sprintf("ts_nexp_%d", k$m), G = 100000L,
      run = function(g, seed) {
        set.seed(seed)
        draw <- function(r) 5 + stats::rexp(r, rate = 1 / 2)
        t(replicate(g, {
          pilot <- draw(k$m)
          n <- attr(ts_fwd_nexp_mean(pilot, d = k$d), "criterion_n")
          c(n = n,
            ok = covers(ts_fwd_nexp_mean(pilot, second_stage(n - k$m, draw), d = k$d), 5))
        }))
      },
      rows = function(parts) two_stage_rows(
        pool_loop(parts), "ts_fwd_nexp_mean()",
        paste0("exact calculation; procedure in ", cite_book), label,
        exact_nexp(k$m, -log(0.05) * 2 / k$d), "guaranteed: at least 0.95"),
      key = c(procedure = "ts_fwd_nexp_mean()", condition = label))
  })

# The two stage bounded risk procedure for an exponential mean (mean 2,
# A = 100): Section 8.3.2 of the book, Equation 8.3.13. The risk bound w is
# chosen to give the stated n* = A theta^2 / w.
conditions_br_exp <- lapply(
  list(list(m = 10L, n_star = 25), list(m = 14L, n_star = 44), list(m = 20L, n_star = 100)),
  function(k) {
    label <- sprintf("exponential, mean 2; A = 100, n* = %g, first stage of %d", k$n_star, k$m)
    condition(
      id = sprintf("ts_br_exp_%d", k$m), G = 200000L,
      run = function(g, seed) {
        set.seed(seed)
        w <- 100 * 2^2 / k$n_star
        draw <- function(r) stats::rexp(r, rate = 1 / 2)
        t(replicate(g, {
          pilot <- draw(k$m)
          n <- attr(ts_br_exp_mean(pilot, A = 100, w = w), "criterion_n")
          chk <- ts_br_exp_mean(pilot, second_stage(n - k$m, draw), A = 100, w = w)
          c(n = n, loss = 100 * (value_of(chk, "mean") - 2)^2 / w)
        }))
      },
      rows = function(parts) two_stage_rows(
        pool_loop(parts), "ts_br_exp_mean()",
        paste0("exact calculation; procedure in ", cite_book,
               ", Section 8.3.2, Equation 8.3.13"),
        label, exact_br_exp(k$m, k$n_star), "bound: at most 1"),
      key = c(procedure = "ts_br_exp_mean()", condition = label))
  })

# The two stage bounded risk procedure for a normal mean (mean 50, standard
# deviation 2, A = 100, n* = 50), with squared error and absolute error loss.
conditions_br_norm <- lapply(c(2, 1), function(k) {
  label <- sprintf(
    "normal, standard deviation 2; %s error loss, A = 100, n* = 50, first stage of 10",
    if (k == 2) "squared" else "absolute")
  condition(
    id = sprintf("ts_br_norm_%g", k), G = 100000L,
    run = function(g, seed) {
      set.seed(seed)
      B <- 100 * 2^(k / 2) * gamma((k + 1) / 2) / sqrt(pi)
      w <- B * (4 / 50)^(k / 2)
      draw <- function(r) stats::rnorm(r, 50, 2)
      t(replicate(g, {
        pilot <- draw(10)
        n <- attr(ts_br_norm_mean(pilot, A = 100, k = k, w = w), "criterion_n")
        chk <- ts_br_norm_mean(pilot, second_stage(n - 10, draw), A = 100, k = k, w = w)
        c(n = n, loss = 100 * abs(value_of(chk, "mean") - 50)^k / w)
      }))
    },
    rows = function(parts) two_stage_rows(
      pool_loop(parts), "ts_br_norm_mean()",
      paste0("exact calculation; procedure in ", cite_book), label,
      exact_br_norm(10, 50, k), "guaranteed: at most 1"),
    key = c(procedure = "ts_br_norm_mean()", condition = label))
})

# ---------------------------------------------------------------------------
# 9. Departures: where the package deliberately differs from an article
# ---------------------------------------------------------------------------
family_departures <- "departures"

# 9a. The unusable variance rule. Under the stopping rule as printed, a
# variance estimate that is zero or negative is replaced by n^-3 and the rule
# is evaluated, so a study can stop on it with a very narrow interval. The
# package reports such a check as not met. The same studies are followed
# under both rules: the package's check decides when the package stops; the
# rule as printed stops at the first stage at which either the check is met
# or the variance is unusable and n >= (2z / width)^2 (n^-3 + 1 / n).
unusable_run <- function(g, seed) {
  set.seed(seed)
  z <- z_of(0.95); a <- (2 * z / 0.2)^2; cap <- 3000L
  start <- seq_aipe_r_pilot(width = 0.2, conf_level = 0.95)
  t(replicate(g, {
    d <- bivariate_normal(cap, 0.5)
    n <- start; printed <- c(NA_real_, NA_real_, NA_real_)
    repeat {
      chk <- seq_aipe_r(d$x[seq_len(n)], d$y[seq_len(n)], width = 0.2, conf_level = 0.95)
      met <- isTRUE(attr(chk, "is_satisfied"))
      if (is.na(printed[1L])) {
        if (met) {
          printed <- c(n, value_of(chk, "lower_limit"), value_of(chk, "upper_limit"))
        } else if (isTRUE(attr(chk, "variance_unusable")) && n >= a * (n^-3 + 1 / n)) {
          half <- z * sqrt(n^-3 / n)
          printed <- c(n, value_of(chk, "r") - half, value_of(chk, "r") + half)
        }
      }
      if (met || n >= cap) break
      n <- n + 1L
    }
    c(n = n, lower = value_of(chk, "lower_limit"), upper = value_of(chk, "upper_limit"),
      capped = as.numeric(!met), printed_n = printed[1L], printed_lower = printed[2L],
      printed_upper = printed[3L])
  }))
}
conditions_unusable <- list(condition(
  id = "departure_unusable_pearson", G = 10000L, run = unusable_run,
  rows = function(parts) {
    m <- pool_loop(parts)
    pub <- list(N = 204.0, N_digits = 1L, N_se = 0.7435, p = 0.9012, w = 0.1938, w_se = 3.83e-4)
    src <- paste0(cite_2019, ", Table 2, p. 498")
    cond <- "bivariate normal, rho = 0.5; width 0.2, 95%"
    printed <- cbind(n = m[, "printed_n"], lower = m[, "printed_lower"],
                     upper = m[, "printed_upper"], capped = m[, "capped"])
    rbind(
      interval_rows(interval_summary_loop(list(printed), 0.5), family_departures,
                    "the stopping rule as printed (not the package's rule)", src, cond, pub,
                    kind = "rule as printed",
                    note = sprintf("stopped on an unusable variance estimate, before the package's rule was met, in %.2f%% of studies",
                                   100 * mean(m[, "printed_n"] < m[, "n"]))),
      interval_rows(interval_summary_loop(list(m), 0.5), family_departures,
                    "seq_aipe_r_pilot(), seq_aipe_r(), Pearson (the same studies)", src, cond,
                    pub, kind = "deliberate departure"))
  }))

# 9b. The variance of the coefficient of variation. For a gamma population
# with a coefficient of variation of 0.8, n times the mean squared error of
# the sample coefficient of variation (base R, samples of 2,000) is compared
# with the two published forms of its asymptotic variance, and so are the two
# variance estimates of seq_mr_cv() (samples of 20,000).
cv_forms <- function(kappa) {
  a <- 1 / kappa^2   # gamma shape; scale 1: mean a, variance a, mu3 2a, mu4 3a^2 + 6a
  c(eq30 = (3 * a^2 + 6 * a) / (4 * a^3) - 1 / (4 * a) - 2 / a^2 + 1 / a^2,
    eq9 = (3 * a^2 + 6 * a) / (4 * a^4) + 1 / (4 * a^2) - 2 / a^2 + 1 / (2 * a))
}
cite_eq30 <- paste0(cite_2018, ", Equation 30, p. 233")
cite_eq9 <- paste0(cite_2016, ", Equation 9, p. 631")
cv_departure_condition <- "gamma, coefficient of variation 0.8 (shape 1.5625)"
conditions_cv_departure <- list(
  condition(
    id = "departure_cv_mse", G = 100000L,
    run = function(g, seed) {
      set.seed(seed)
      replicate(g, {
        x <- stats::rgamma(2000, shape = 1.5625)
        2000 * (stats::sd(x) / mean(x) - 0.8)^2
      })
    },
    rows = function(parts) {
      v <- unlist(parts); f <- cv_forms(0.8)
      one <- function(src, published) result_row(
        family_departures, "sample coefficient of variation (base R)", src,
        paste0(cv_departure_condition, "; samples of 2,000"),
        "n times the mean squared error of the estimate", published, mean(v),
        se = batch_se(v), replications = length(v), kind = "sampling variance")
      rbind(one(cite_eq30, f[["eq30"]]), one(cite_eq9, f[["eq9"]]))
    }),
  condition(
    id = "departure_cv_estimators", G = 400L,
    run = function(g, seed) {
      set.seed(seed)
      t(replicate(g, {
        x <- stats::rgamma(20000, shape = 1.5625)
        c(corrected = value_of(seq_mr_cv(x, A = 5e6, c = 10), "xi_squared"),
          published = value_of(seq_mr_cv(x, A = 5e6, c = 10, variance = "published"),
                               "xi_squared"))
      }))
    },
    rows = function(parts) {
      m <- pool_loop(parts); f <- cv_forms(0.8)
      one <- function(variance, src, published) result_row(
        family_departures, sprintf("seq_mr_cv(variance = \"%s\"), row xi_squared", variance),
        src, paste0(cv_departure_condition, "; samples of 20,000"), "variance estimate, mean",
        published, mean(m[, variance]), se = batch_se(m[, variance]),
        replications = nrow(m), kind = "estimator")
      rbind(one("corrected", cite_eq30, f[["eq30"]]), one("published", cite_eq9, f[["eq9"]]))
    }))

# 9c. Equation 26 of the 2018 article. Two populations with unit variance and
# a standardized mean difference of 0.5: a standardized chi square with 3
# degrees of freedom and its mirror image (third central moments sqrt(8/3)
# and -sqrt(8/3), fourth central moments 7). n times the variance of the
# sample standardized mean difference (base R, 1,000 per group) and the mean
# variance estimate of seq_aipe_smd() (5,000 per group) are compared with
# Equation 26 as printed and as corrected.
eq26_forms <- c(
  printed = 2 - 0.5 * (2 * sqrt(8 / 3)) + (0.5^2 / 4) * ((7 + 7) / 4 - 1 / 2),
  corrected = 2 - 0.5 * (2 * sqrt(8 / 3)) / 2 + (0.5^2 / 4) * ((7 + 7) / 4 - 1 / 2))
eq26_draw <- function(n) {
  s <- function(r) (stats::rchisq(r, 3) - 3) / sqrt(6)
  list(g1 = 0.5 + s(n), g2 = -s(n))
}
cite_eq26_printed <- paste0(cite_2018, ", Equation 26 as printed, p. 232")
cite_eq26_corrected <- paste0(cite_2018, ", Equation 26 as corrected")
eq26_condition <- "skewed in opposite directions, unit variances, delta = 0.5"
conditions_eq26 <- list(
  condition(
    id = "departure_eq26_variance", G = 40000L,
    run = function(g, seed) {
      set.seed(seed)
      replicate(g, {
        d <- eq26_draw(1000)
        (mean(d$g1) - mean(d$g2)) / sqrt((stats::var(d$g1) + stats::var(d$g2)) / 2)
      })
    },
    rows = function(parts) {
      d <- unlist(parts)
      v <- 1000 * (d - mean(d))^2 * length(d) / (length(d) - 1)
      one <- function(src, published) result_row(
        family_departures, "sample standardized mean difference (base R)", src,
        paste0(eq26_condition, "; 1,000 per group"), "n times the variance of the estimate",
        published, mean(v), se = batch_se(v), replications = length(d),
        kind = "sampling variance")
      rbind(one(cite_eq26_corrected, eq26_forms[["corrected"]]),
            one(cite_eq26_printed, eq26_forms[["printed"]]))
    }),
  condition(
    id = "departure_eq26_estimator", G = 400L,
    run = function(g, seed) {
      set.seed(seed)
      replicate(g, {
        d <- eq26_draw(5000)
        value_of(seq_aipe_smd(d$g1, d$g2, width = 0.2), "xi_squared")
      })
    },
    rows = function(parts) {
      v <- unlist(parts)
      one <- function(src, published) result_row(
        family_departures, "seq_aipe_smd(), row xi_squared", src,
        paste0(eq26_condition, "; 5,000 per group"), "variance estimate, mean",
        published, mean(v), se = batch_se(v), replications = length(v), kind = "estimator")
      rbind(one(cite_eq26_corrected, eq26_forms[["corrected"]]),
            one(cite_eq26_printed, eq26_forms[["printed"]]))
    }))

# 9d. The stopping rule for the squared multiple correlation with the Wald
# interval. Equation 63 of the 2019 article (p. 505) as printed is
# n >= (16 z^2 / width^2) (R2 + 1/n) (1 - (R2 + 1/n))^2; the package
# implements (R2 + 1/n) (1 - R2 + 1/n)^2, which reproduces Tables 13 and 14
# (section 3 above). The rule as printed is evaluated here from the estimate
# that seq_aipe_R2() reports at each stage, with the package's pilot sample
# size and interval, for the two rows where the forms differ most.
r2_printed <- list(
  list(table = 13L, page = 506L, conf = 0.90,
       pub = list(N = 68.66, N_digits = 2L, N_se = 0.0963, p = 0.8678, w = 0.0733, w_se = NA)),
  list(table = 14L, page = 507L, conf = 0.95,
       pub = list(N = 84.17, N_digits = 2L, N_se = 0.1493, p = 0.9222, w = 0.0790, w_se = NA)))
conditions_r2_printed <- lapply(r2_printed, function(k) condition(
  id = sprintf("departure_r2_printed_%d", k$table), G = 5000L,
  run = function(g, seed) {
    set.seed(seed)
    z <- z_of(k$conf); cap <- 600L; b <- sqrt(0.9 / (0.1 * 2))
    start <- seq_aipe_R2_pilot(width = 0.10, p = 2, conf_level = k$conf)
    t(replicate(g, {
      x <- matrix(stats::rnorm(2 * cap), cap, 2)
      y <- drop(x %*% c(b, b)) + stats::rnorm(cap)
      n <- start
      repeat {
        chk <- seq_aipe_R2(x[seq_len(n), , drop = FALSE], y[seq_len(n)], width = 0.10,
                           conf_level = k$conf)
        r2 <- value_of(chk, "R2")
        if (n >= 16 * z^2 / 0.10^2 * (r2 + 1 / n) * (1 - (r2 + 1 / n))^2 || n >= cap) break
        n <- n + 1L
      }
      c(n = n, lower = value_of(chk, "lower_limit"), upper = value_of(chk, "upper_limit"),
        capped = as.numeric(n >= cap))
    }))
  },
  rows = function(parts) interval_rows(
    interval_summary_loop(parts, truth = 0.9), family_departures,
    "Equation 63 as printed (not the package's rule), on seq_aipe_R2() estimates",
    sprintf("%s, Table %d, p. %d", cite_2019, k$table, k$page),
    sprintf("2 predictors, population value 0.9; width 0.10, %g%%", 100 * k$conf),
    k$pub, kind = "rule as printed")))

# ---------------------------------------------------------------------------
# 10. Wide targets: coverage, or risk, where the guarantee is asymptotic
# ---------------------------------------------------------------------------
# The intervals for effect sizes: each condition has 10,000 replications, or
# 5,000 where a fifth element says so (the conditions with the largest final
# sample sizes).
family_wide <- "wide_targets"
wide_row <- function(s, procedure, condition) {
  if (!is.na(s$pct_capped) && s$pct_capped > 0) {
    stop("A replication reached the cap in: ", condition, call. = FALSE)
  }
  result_row(family_wide, procedure, "nominal level", condition, "coverage", 0.95,
             s$coverage, se = binom_se(s$coverage, s$G), replications = s$G,
             kind = "nominal level", mean_n = s$mean_n,
             note = if (!is.na(s$pct_unusable) && s$pct_unusable > 0)
               sprintf("variance estimate unusable at one or more checks in %.2f%% of studies",
                       s$pct_unusable) else "")
}
wide_sensitivity <- list(
  list("seq_aipe_smd()", seq_aipe_smd_sensitivity,
       "standardized mean difference 0.5, normal; width 0.8",
       list(true_value = 0.5, width = 0.8)),
  list("seq_aipe_smd()", seq_aipe_smd_sensitivity,
       "standardized mean difference 0.5, skewed; width 0.4",
       list(true_value = 0.5, width = 0.4, distribution = "skewed")),
  list("seq_aipe_sm()", seq_aipe_sm_sensitivity, "standardized mean 0.4, normal; width 0.5",
       list(true_value = 0.4, width = 0.5)),
  list("seq_aipe_c()", seq_aipe_c_sensitivity,
       "contrast 0.5 (weights 1, -0.5, -0.5), normal; width 1",
       list(true_means = c(0.5, 0, 0), contrast_weights = c(1, -0.5, -0.5), width = 1)),
  list("seq_aipe_r(), Pearson", seq_aipe_r_sensitivity, "correlation 0.4, normal; width 0.4",
       list(true_value = 0.4, width = 0.4)),
  list("seq_aipe_r(), Pearson", seq_aipe_r_sensitivity,
       "correlation 0.4, normal; width 0.4, m_0 = 40",
       list(true_value = 0.4, width = 0.4, m_0 = 40)),
  list("seq_aipe_r(), Pearson", seq_aipe_r_sensitivity, "correlation 0.4, normal; width 0.2",
       list(true_value = 0.4, width = 0.2), 5000L),
  list("seq_aipe_r(), Pearson", seq_aipe_r_sensitivity, "correlation 0.4, skewed; width 0.4",
       list(true_value = 0.4, width = 0.4, distribution = "skewed")),
  list("seq_aipe_r(), Pearson", seq_aipe_r_sensitivity, "correlation 0.8, normal; width 0.4",
       list(true_value = 0.8, width = 0.4)),
  list("seq_aipe_r(), Pearson", seq_aipe_r_sensitivity, "correlation 0.8, normal; width 0.2",
       list(true_value = 0.8, width = 0.2)),
  list("seq_aipe_r(), Pearson", seq_aipe_r_sensitivity, "correlation 0.9, normal; width 0.1",
       list(true_value = 0.9, width = 0.1)),
  list("seq_aipe_reg_coef()", seq_aipe_reg_coef_sensitivity, "slope 0.5, normal; width 0.4",
       list(true_value = 0.5, width = 0.4)),
  list("seq_aipe_reg_coef()", seq_aipe_reg_coef_sensitivity,
       "slope 0.5, normal; width 0.4, m_0 = 40",
       list(true_value = 0.5, width = 0.4, m_0 = 40)),
  list("seq_aipe_reg_coef()", seq_aipe_reg_coef_sensitivity, "slope 0.5, normal; width 0.2",
       list(true_value = 0.5, width = 0.2), 5000L),
  list("seq_aipe_reg_coef()", seq_aipe_reg_coef_sensitivity, "slope 0.5, skewed; width 0.4",
       list(true_value = 0.5, width = 0.4, distribution = "skewed")),
  list("seq_aipe_reg_coef()", seq_aipe_reg_coef_sensitivity, "slope 2, normal; width 0.2",
       list(true_value = 2, width = 0.2), 5000L),
  list("seq_aipe_cv()", seq_aipe_cv_sensitivity,
       "coefficient of variation 0.5, gamma; width 0.2",
       list(true_value = 0.5, width = 0.2, distribution = "gamma")),
  list("seq_aipe_cv()", seq_aipe_cv_sensitivity,
       "coefficient of variation 0.5, gamma; width 0.1",
       list(true_value = 0.5, width = 0.1, distribution = "gamma")),
  list("seq_aipe_cv()", seq_aipe_cv_sensitivity,
       "coefficient of variation 0.5, lognormal; width 0.2",
       list(true_value = 0.5, width = 0.2, distribution = "lognormal")),
  list("seq_aipe_cv()", seq_aipe_cv_sensitivity,
       "coefficient of variation 0.5, lognormal; width 0.1",
       list(true_value = 0.5, width = 0.1, distribution = "lognormal"), 5000L),
  list("seq_aipe_gini()", seq_aipe_gini_sensitivity, "Gini index 0.3, lognormal; width 0.2",
       list(true_value = 0.3, width = 0.2)),
  list("seq_aipe_gini()", seq_aipe_gini_sensitivity, "Gini index 0.3, lognormal; width 0.1",
       list(true_value = 0.3, width = 0.1)),
  list("seq_aipe_R2(), Wald interval", seq_aipe_R2_sensitivity,
       "squared multiple correlation 0.5, 5 predictors; width 0.3",
       list(true_value = 0.5, width = 0.3, p = 5)),
  list("seq_aipe_R2(), Wald interval", seq_aipe_R2_sensitivity,
       "squared multiple correlation 0.5, 10 predictors; width 0.2",
       list(true_value = 0.5, width = 0.2, p = 10)),
  list("seq_aipe_R2(), Wald interval", seq_aipe_R2_sensitivity,
       "squared multiple correlation 0.5, 10 predictors; width 0.3",
       list(true_value = 0.5, width = 0.3, p = 10)),
  list("seq_aipe_R2(), Wald interval", seq_aipe_R2_sensitivity,
       "squared multiple correlation 0.5, 15 predictors; width 0.3",
       list(true_value = 0.5, width = 0.3, p = 15)))
conditions_wide <- lapply(seq_along(wide_sensitivity), function(i) {
  k <- wide_sensitivity[[i]]
  condition(
    id = sprintf("wide_%02d", i), G = if (length(k) == 5L) k[[5L]] else 10000L,
    run = sensitivity_run(k[[2L]], c(k[[4L]], list(conf_level = 0.95))),
    rows = function(parts) wide_row(interval_summary_sensitivity(parts), k[[1L]], k[[3L]]))
})
# The rank correlations at a strong association, through the check function.
conditions_wide_rank <- lapply(c("kendall", "spearman"), function(method) condition(
  id = paste0("wide_", method), G = 10000L, run = rank_run(method, 0.8, 0.4, 1500L),
  rows = function(parts) wide_row(
    interval_summary_loop(parts, truth = rank_truth(method, 0.8)),
    sprintf("seq_aipe_r(), %s", rank_names[[method]]),
    "bivariate normal, rho = 0.8; width 0.4")))
# The distribution free two stage interval for a mean: first stage of 30,
# populations with standard deviation 2. The normal rows also have an exact
# value (the Stein calculation with a normal quantile in place of t).
dfr_populations <- list(
  list(id = "normal", label = "normal", mean = 50, draw = function(r) stats::rnorm(r, 50, 2)),
  list(id = "exponential", label = "exponential (skewness 2)", mean = 2,
       draw = function(r) stats::rexp(r, rate = 1 / 2)),
  list(id = "lognormal", label = "lognormal (skewness 4)", mean = 2,
       draw = function(r) stats::rlnorm(r, log(2) / 2, sqrt(log(2)))))
dfr_grid <- expand.grid(d = c(0.5, 0.2), population = seq_along(dfr_populations))
conditions_wide_dfr <- lapply(seq_len(nrow(dfr_grid)), function(i) {
  pop <- dfr_populations[[dfr_grid$population[i]]]
  d <- dfr_grid$d[i]
  cond <- sprintf("%s, standard deviation 2; half width %g, first stage of 30", pop$label, d)
  condition(
    id = sprintf("wide_dfr_%s_%g", pop$id, d), G = 20000L,
    key = if (pop$id == "normal") c(procedure = "ts_fwd_dfr_mean()", condition = cond),
    run = function(g, seed) {
      set.seed(seed)
      t(replicate(g, {
        pilot <- pop$draw(30)
        n <- attr(ts_fwd_dfr_mean(pilot, d = d), "criterion_n")
        c(n = n, ok = covers(ts_fwd_dfr_mean(pilot, second_stage(n - 30, pop$draw), d = d),
                             pop$mean))
      }))
    },
    rows = function(parts) {
      m <- pool_loop(parts)
      s <- list(G = nrow(m), coverage = mean(m[, "ok"]), mean_n = mean(m[, "n"]),
                pct_unusable = NA_real_, pct_capped = NA_real_)
      out <- wide_row(s, "ts_fwd_dfr_mean()", cond)
      if (pop$id == "normal") {
        exact <- exact_stein(30, z_of(0.95)^2 * 4 / d^2, quantile = "z")
        out <- rbind(out, result_row(
          family_wide, "ts_fwd_dfr_mean()", "exact calculation", cond, "coverage",
          exact[["second"]], s$coverage, se = binom_se(s$coverage, s$G),
          replications = s$G, kind = "exact calculation", mean_n = s$mean_n))
      }
      out
    })
})

# The Gini index at wide targets, with two intervals in the same studies.
# seq_aipe_gini() reports the estimate plus and minus z times its estimated
# standard error, which the stopping rule keeps below the target width; the
# interval of Chattopadhyay and De (2016) is the estimate plus and minus half
# the target width. The two are close when the final sample is large (section
# 7) and far apart when it is small. Lognormal data with a Gini index of 0.3,
# as in the rows of seq_aipe_gini_sensitivity() above, run here through the
# pilot and check functions so that both intervals can be formed.
conditions_wide_gini <- lapply(c(0.2, 0.1), function(width) condition(
  id = sprintf("wide_gini_intervals_%g", width), G = 10000L,
  run = function(g, seed) {
    set.seed(seed)
    cap <- 3000L
    sdlog <- sqrt(2) * stats::qnorm((0.3 + 1) / 2)   # Gini index 2 pnorm(sdlog / sqrt(2)) - 1
    start <- seq_aipe_gini_pilot(width = width, conf_level = 0.95)
    t(replicate(g, {
      x <- stats::rlnorm(cap, 0, sdlog)
      run_to_stopping(function(n) seq_aipe_gini(x[seq_len(n)], width = width,
                                                conf_level = 0.95),
                      start, cap)
    }))
  },
  rows = function(parts) {
    m <- pool_loop(parts)
    cond <- sprintf("Gini index 0.3, lognormal; width %g", width)
    if (any(m[, "capped"] > 0)) stop("A replication reached the cap in: ", cond, call. = FALSE)
    middle <- (m[, "lower"] + m[, "upper"]) / 2
    reported <- m[, "upper"] - m[, "lower"]
    stopifnot(all(reported <= width * (1 + 1e-12)))
    one <- function(quantity, covered, note) result_row(
      family_wide, "seq_aipe_gini_pilot(), seq_aipe_gini()", "nominal level", cond, quantity,
      0.95, mean(covered), se = binom_se(mean(covered), nrow(m)), replications = nrow(m),
      kind = "nominal level", mean_n = mean(m[, "n"]), note = note)
    rbind(
      one("coverage of the interval the package reports",
          m[, "lower"] <= 0.3 & 0.3 <= m[, "upper"],
          sprintf("mean width %.4f", mean(reported))),
      one("coverage of the article's interval (the estimate plus and minus half the width)",
          abs(middle - 0.3) <= width / 2, sprintf("width %g", width)))
  }))

# The purely sequential classical procedures (Mukhopadhyay and de Silva,
# 2009). Their confidence level, their risk bound, or their optimality holds
# in the limit as the target narrows. Each is run through its pilot function
# and its check function, one observation per stage, with the default pilot,
# at a target near the least favorable one that the function's help page
# reports and, for the two intervals for a mean, at a narrower target as
# well. n* is the sample size that would be needed if the population were
# known:
#   fixed width interval for a mean (half width d): z^2 sigma^2 / d^2;
#   interval for the location of a negative exponential distribution (width
#     d): log(1 / alpha) sigma / d;
#   interval for a median (length at most 2 d): z^2 / (4 d^2 f^2), with f the
#     density at the median (0.3989 for a standard normal population and 0.5
#     for a unit exponential one);
#   bounded risk for a normal mean, loss A |error|^k, bound w:
#     (B / w)^(2 / k) sigma^2 with B = A 2^(k / 2) gamma((k + 1) / 2) / sqrt(pi);
#   minimum risk, loss A error^2 + c n: sqrt(A / c) sigma for a normal mean
#     and sqrt(A / c) times the mean for an exponential mean, where the
#     smallest attainable risk is 2 c n*.
run_to_stopping_estimate <- function(check, start, cap, term = "mean") {
  n <- as.integer(start)
  repeat {
    chk <- check(n)
    if (isTRUE(attr(chk, "is_satisfied")) || n >= cap) break
    n <- n + 1L
  }
  c(n = n, estimate = value_of(chk, term),
    capped = as.numeric(!isTRUE(attr(chk, "is_satisfied"))))
}
classical_row <- function(m, procedure, condition, quantity, reference, kind, value, se) {
  if (any(m[, "capped"] > 0)) {
    stop("A replication reached the cap in: ", procedure, " ", condition, call. = FALSE)
  }
  result_row(family_wide, procedure, reference$source, condition, quantity,
             reference$value, value, se = se, replications = nrow(m), kind = kind,
             mean_n = mean(m[, "n"]))
}
# An interval procedure: coverage against the nominal 0.95.
classical_interval <- function(id, G, cap, procedure, label, pilot, check, draw, truth) {
  force(pilot); force(check); force(draw); force(truth); force(cap)
  condition(
    id = paste0("classical_", id), G = G,
    run = function(g, seed) {
      set.seed(seed)
      start <- pilot()
      t(replicate(g, {
        x <- draw(cap)
        run_to_stopping(function(n) check(x[seq_len(n)]), start, cap)
      }))
    },
    rows = function(parts) {
      m <- pool_loop(parts)
      coverage <- mean(m[, "lower"] <= truth & truth <= m[, "upper"])
      classical_row(m, procedure, label, "coverage",
                    list(source = "nominal level", value = 0.95), "nominal level",
                    coverage, binom_se(coverage, nrow(m)))
    })
}
# A risk procedure: the mean of loss(estimate, n) over the studies, where
# `loss` already divides by the bound or by the smallest attainable risk.
classical_risk <- function(id, G, cap, procedure, label, quantity, source, kind,
                           pilot, check, draw, loss) {
  force(pilot); force(check); force(draw); force(loss); force(cap)
  condition(
    id = paste0("classical_", id), G = G,
    run = function(g, seed) {
      set.seed(seed)
      start <- pilot()
      t(replicate(g, {
        x <- draw(cap)
        run_to_stopping_estimate(function(n) check(x[seq_len(n)]), start, cap)
      }))
    },
    rows = function(parts) {
      m <- pool_loop(parts)
      relative <- loss(m[, "estimate"], m[, "n"])
      classical_row(m, procedure, label, quantity, list(source = source, value = 1), kind,
                    mean(relative), batch_se(relative))
    })
}
n_star_mean <- function(d, sigma = 2) z_of(0.95)^2 * sigma^2 / d^2
n_star_median <- function(d, density) z_of(0.95)^2 / (4 * d^2 * density^2)
br_w <- c(squared = 0.65, absolute = 0.92)     # risk bounds: n* near 27 to 28
br_B <- c(squared = 2, absolute = 2 * sqrt(2 / pi))
conditions_classical <- list(
  classical_interval(
    "fwd_norm_0.77", 40000L, 500L, "seq_fwd_norm_mean()",
    sprintf("normal, standard deviation 2; half width 0.77 (n* = %.0f)", n_star_mean(0.77)),
    pilot = function() seq_fwd_norm_mean_pilot(d = 0.77),
    check = function(x) seq_fwd_norm_mean(x, d = 0.77),
    draw = function(r) stats::rnorm(r, 50, 2), truth = 50),
  classical_interval(
    "fwd_norm_0.25", 10000L, 2000L, "seq_fwd_norm_mean()",
    sprintf("normal, standard deviation 2; half width 0.25 (n* = %.0f)", n_star_mean(0.25)),
    pilot = function() seq_fwd_norm_mean_pilot(d = 0.25),
    check = function(x) seq_fwd_norm_mean(x, d = 0.25),
    draw = function(r) stats::rnorm(r, 50, 2), truth = 50),
  classical_interval(
    "fwd_dfr_normal_0.5", 20000L, 2000L, "seq_fwd_dfr_mean()",
    sprintf("normal, standard deviation 2; half width 0.5 (n* = %.0f)", n_star_mean(0.5)),
    pilot = function() seq_fwd_dfr_mean_pilot(d = 0.5),
    check = function(x) seq_fwd_dfr_mean(x, d = 0.5),
    draw = function(r) stats::rnorm(r, 50, 2), truth = 50),
  classical_interval(
    "fwd_dfr_exponential_0.5", 20000L, 2000L, "seq_fwd_dfr_mean()",
    sprintf("exponential, standard deviation 2; half width 0.5 (n* = %.0f)", n_star_mean(0.5)),
    pilot = function() seq_fwd_dfr_mean_pilot(d = 0.5),
    check = function(x) seq_fwd_dfr_mean(x, d = 0.5),
    draw = function(r) stats::rexp(r, rate = 1 / 2), truth = 2),
  classical_interval(
    "fwd_dfr_exponential_0.25", 10000L, 4000L, "seq_fwd_dfr_mean()",
    sprintf("exponential, standard deviation 2; half width 0.25 (n* = %.0f)", n_star_mean(0.25)),
    pilot = function() seq_fwd_dfr_mean_pilot(d = 0.25),
    check = function(x) seq_fwd_dfr_mean(x, d = 0.25),
    draw = function(r) stats::rexp(r, rate = 1 / 2), truth = 2),
  classical_interval(
    "fwd_nexp_0.24", 40000L, 500L, "seq_fwd_nexp_mean()",
    sprintf("negative exponential, scale 2; width 0.24 (n* = %.0f)", -log(0.05) * 2 / 0.24),
    pilot = function() seq_fwd_nexp_mean_pilot(d = 0.24),
    check = function(x) seq_fwd_nexp_mean(x, d = 0.24),
    draw = function(r) 5 + stats::rexp(r, rate = 1 / 2), truth = 5),
  classical_interval(
    "blci_median_normal_100", 10000L, 2000L, "seq_blci_dfr_median()",
    sprintf("standard normal; d = 0.2456 (n* = %.0f)",
            n_star_median(0.2456, stats::dnorm(0))),
    pilot = function() seq_blci_dfr_median_pilot(d = 0.2456),
    check = function(x) seq_blci_dfr_median(x, d = 0.2456),
    draw = function(r) stats::rnorm(r), truth = 0),
  classical_interval(
    "blci_median_exponential_100", 10000L, 2000L, "seq_blci_dfr_median()",
    sprintf("unit exponential; d = 0.196 (n* = %.0f)", n_star_median(0.196, 0.5)),
    pilot = function() seq_blci_dfr_median_pilot(d = 0.196),
    check = function(x) seq_blci_dfr_median(x, d = 0.196),
    draw = function(r) stats::rexp(r), truth = log(2)),
  classical_interval(
    "blci_median_normal_800", 5000L, 6000L, "seq_blci_dfr_median()",
    sprintf("standard normal; d = 0.08685 (n* = %.0f)",
            n_star_median(0.08685, stats::dnorm(0))),
    pilot = function() seq_blci_dfr_median_pilot(d = 0.08685),
    check = function(x) seq_blci_dfr_median(x, d = 0.08685),
    draw = function(r) stats::rnorm(r), truth = 0),
  classical_risk(
    "br_norm_squared", 50000L, 500L, "seq_br_norm_mean()",
    sprintf("normal, standard deviation 3; squared error loss, A = 2, w = %g (n* = %.0f)",
            br_w[["squared"]], (br_B[["squared"]] / br_w[["squared"]]) * 9),
    "risk relative to the bound", "risk bound", "asymptotic bound",
    pilot = function() seq_br_norm_mean_pilot(A = 2, k = 2, w = br_w[["squared"]]),
    check = function(x) seq_br_norm_mean(x, A = 2, k = 2, w = br_w[["squared"]]),
    draw = function(r) stats::rnorm(r, 50, 3),
    loss = function(estimate, n) 2 * (estimate - 50)^2 / br_w[["squared"]]),
  classical_risk(
    "br_norm_absolute", 50000L, 500L, "seq_br_norm_mean()",
    sprintf("normal, standard deviation 3; absolute error loss, A = 2, w = %g (n* = %.0f)",
            br_w[["absolute"]], (br_B[["absolute"]] / br_w[["absolute"]])^2 * 9),
    "risk relative to the bound", "risk bound", "asymptotic bound",
    pilot = function() seq_br_norm_mean_pilot(A = 2, k = 1, w = br_w[["absolute"]]),
    check = function(x) seq_br_norm_mean(x, A = 2, k = 1, w = br_w[["absolute"]]),
    draw = function(r) stats::rnorm(r, 50, 3),
    loss = function(estimate, n) 2 * abs(estimate - 50) / br_w[["absolute"]]),
  classical_risk(
    "mr_norm_10", 50000L, 500L, "seq_mr_norm_mean()",
    "normal, standard deviation 1; A = 100, c = 1 (n* = 10)",
    "risk relative to the smallest attainable", "smallest attainable risk",
    "asymptotic optimum",
    pilot = function() seq_mr_norm_mean_pilot(A = 100, c = 1),
    check = function(x) seq_mr_norm_mean(x, A = 100, c = 1),
    draw = function(r) stats::rnorm(r, 50, 1),
    loss = function(estimate, n) (100 * (estimate - 50)^2 + n) / (2 * 10)),
  classical_risk(
    "mr_exp_15", 50000L, 500L, "seq_mr_exp_mean()",
    "exponential, mean 1; A = 225, c = 1 (n* = 15)",
    "risk relative to the smallest attainable", "smallest attainable risk",
    "asymptotic optimum",
    pilot = function() seq_mr_exp_mean_pilot(A = 225, c = 1),
    check = function(x) seq_mr_exp_mean(x, A = 225, c = 1),
    draw = function(r) stats::rexp(r),
    loss = function(estimate, n) (225 * (estimate - 1)^2 + n) / (2 * 15)))

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
# The seeds of a condition depend on its position in this list, so a new
# condition is added at the end, whatever section defines it: the results of
# the conditions before it then stay as they were. (The rows are put in the
# order of the sections when the results are assembled.)
conditions <- c(conditions_smd_2018, conditions_cv_2018, conditions_pearson_2019,
                conditions_rank_2019, conditions_r2_2019, conditions_2016, conditions_2017,
                conditions_2025, conditions_gini, conditions_stein, conditions_nexp,
                conditions_br_exp, conditions_br_norm, conditions_unusable,
                conditions_cv_departure, conditions_eq26, conditions_r2_printed,
                conditions_wide, conditions_wide_rank, conditions_wide_dfr,
                conditions_2025_literal, conditions_gini_procedure, conditions_classical,
                conditions_gini_pareto, conditions_wide_gini)
ids <- vapply(conditions, function(k) k$id, character(1))
stopifnot(!anyDuplicated(ids))

# Run the conditions `which` (indices into `conditions`), each as `chunks`
# parts of `multiple` times its usual size, with seeds that depend only on the
# condition, the part, and `offset`. The parts are handed out one at a time as
# cores become free (mc.preschedule = FALSE), so the cores stay equally busy
# whatever a part costs. Returns the result rows of those conditions.
run_conditions <- function(which, multiple = 1L, offset = 0L) {
  jobs <- do.call(rbind, lapply(which, function(i) {
    g <- multiple * max(2L, as.integer(round(conditions[[i]]$G * scale / chunks)))
    data.frame(condition = i, part = seq_len(chunks), g = g,
               seed = master_seed + 1000L * i + offset + seq_len(chunks))
  }))
  done <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
    conditions[[jobs$condition[j]]]$run(jobs$g[j], jobs$seed[j])
  }, mc.cores = cores, mc.preschedule = FALSE)
  failed <- vapply(done, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(failed)) {
    stop("These conditions failed: ",
         paste(unique(ids[jobs$condition[failed]]), collapse = ", "),
         "\nFirst error: ", as.character(done[[which(failed)[1L]]]), call. = FALSE)
  }
  do.call(rbind, lapply(which, function(i) conditions[[i]]$rows(done[jobs$condition == i])))
}

# The standardized difference: the package's value minus the published value,
# divided by the square root of the sum of three variances: the Monte Carlo
# variance of the package's value, the Monte Carlo variance of the published
# value, and the variance of rounding the published value to the decimals
# printed (a uniform error of one unit in the last place). Every published
# value that is itself a simulation mean has a standard error: the one its
# article prints, or one derived as described at interval_rows() and in
# sections 1 and 6. A reference value that is not simulated (an exact value,
# a closed form, a population value, a nominal level, a bound) has none, and
# a missing standard error contributes nothing; when all three variances are
# absent the difference has no scale and z is NA.
with_z <- function(results) {
  or_zero <- function(x) ifelse(is.na(x), 0, x)
  rounding_var <- ifelse(is.na(results$published_digits), 0,
                         (10^-results$published_digits)^2 / 12)
  denominator <- sqrt(or_zero(results$se)^2 + or_zero(results$published_se)^2 + rounding_var)
  results$z <- ifelse(denominator > 0, (results$value - results$published) / denominator,
                      NA_real_)
  results
}
second_run_label <- "exact calculation, second run"

if (resume) {
  # The first runs are taken from the file already written.
  main <- readRDS(out_file)
  stopifnot(identical(attr(main, "seed"), master_seed), identical(attr(main, "scale"), scale),
            !nzchar(only))
  minutes <- attr(main, "minutes")
  processor <- attr(main, "processor_minutes")
  if (is.null(processor)) processor <- NA_real_
  main <- main[main$kind != second_run_label, setdiff(names(main), "z")]
  attr(main, "minutes") <- NULL
} else {
  selected <- if (nzchar(only)) grep(only, ids) else seq_along(conditions)
  if (!length(selected)) stop("No condition matches SMSD_VALIDATION_ONLY.", call. = FALSE)
  message(sprintf("%d conditions, %d parts, %d cores, scale %g", length(selected),
                  chunks * length(selected), cores, scale))
  main <- run_conditions(selected)
  if (!nzchar(only)) {
    main <- rbind(main, ck2016_example(), ck2017_example(), budget_closed_forms())
  }
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  processor <- processor_minutes() - processor_started
}
main <- with_z(main)

# Second runs. A simulated value more than three standard errors from an exact
# value is either chance or a defect, and a second run tells which. Every
# condition with such a row is run again with new seeds and ten times the
# replications, and its comparisons with the exact values are added to the
# results with the kind "exact calculation, second run".
second_started <- Sys.time()
second_processor_started <- processor_minutes()
far <- main[main$kind == "exact calculation" & !is.na(main$z) & abs(main$z) > 3, ]
again <- which(vapply(conditions, function(k) {
  !is.null(k$key) && any(far$procedure == k$key[["procedure"]] &
                           far$condition == k$key[["condition"]])
}, logical(1)))
second <- NULL
if (length(again)) {
  message("Second runs: ", paste(ids[again], collapse = ", "))
  second <- run_conditions(again, multiple = 10L, offset = 500L)
  second <- second[second$kind == "exact calculation", ]
  second$kind <- second_run_label
}

results <- with_z(rbind(main[, setdiff(names(main), "z")], second))
family_order <- c(family_2018, family_2019, family_r2, family_2016, family_2017,
                  family_2025, family_gini, family_classical, family_departures, family_wide)
results <- results[order(match(results$family, family_order)), ]
rownames(results) <- NULL

attr(results, "seed") <- master_seed
attr(results, "date") <- format(Sys.Date())
attr(results, "package_version") <- as.character(utils::packageVersion("SMSD"))
attr(results, "r_version") <- R.version.string
attr(results, "scale") <- scale
attr(results, "only") <- only
attr(results, "parts") <- chunks
attr(results, "cores") <- cores
attr(results, "minutes") <- minutes
attr(results, "processor_minutes") <- processor
attr(results, "second_run_minutes") <- as.numeric(difftime(Sys.time(), second_started,
                                                           units = "mins"))
attr(results, "second_run_processor_minutes") <- processor_minutes() - second_processor_started
saveRDS(results, out_file)
message(sprintf(paste("Wrote %s: %d rows; first runs %.1f minutes (%.0f of processor time),",
                      "second runs %.1f minutes (%.0f)."),
                out_file, nrow(results), attr(results, "minutes"),
                attr(results, "processor_minutes"), attr(results, "second_run_minutes"),
                attr(results, "second_run_processor_minutes")))
