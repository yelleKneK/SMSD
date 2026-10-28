# ---------------------------------------------------------------------------
# Monte Carlo verification of two sequential AIPE procedures in 'SMSD'.
#
# This script PRECOMPUTES the results displayed in the vignette
# "complete-tour". It is deliberately NOT run when the package is built or
# checked: the study takes hours of processor time (see "Run time" below).
# The results are written once to
#
#     inst/extdata/aipe_simulation_results.rds
#
# and the vignette loads that file with system.file(). Re-run this script
# whenever the procedures or the design change, and then read every table
# and sentence of the vignette's verification section against the new file.
#
# Procedures. The sequential accuracy in parameter estimation (AIPE)
# procedures for the standardized mean difference and for the Pearson
# correlation:
#   Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#     parameter estimation for a general class of effect sizes: A
#     sequential approach. Psychological Methods, 23(2), 226-243.
#     doi:10.1037/met0000127
#   Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential
#     accuracy in parameter estimation for population correlation
#     coefficients. Psychological Methods, 24(4), 492-515.
#     doi:10.1037/met0000203
#
# Design. Two questions are put to them:
#   * the DISTRIBUTION FREE property, by running conditions under both a
#     normal and a right skewed data generating mechanism; and
#   * the ASYMPTOTIC nature of the coverage guarantee, by sweeping the
#     target width for the Pearson correlation, where a wide target stops
#     the procedure in a small sample.
#
# The 18 conditions:
#   * standardized mean difference: delta = 0.2, 0.5, 0.8; normal and skewed
#     data; width 0.40 (6 conditions);
#   * Pearson correlation: rho = 0.3, 0.5, 0.7; normal data at widths 0.30,
#     0.15, and 0.08; skewed data at width 0.08 (12 conditions).
#
# Every replication runs the shipped procedure as a prospective study
# would: the pilot sample of seq_aipe_smd_pilot() or seq_aipe_r_pilot(),
# then seq_aipe_smd() or seq_aipe_r() after EVERY added observation (per
# group for the standardized mean difference, pair for the correlation),
# stopping at the first check whose is_satisfied attribute is TRUE. The
# loop and the data generating mechanisms are those of the package's own
# seq_aipe_smd_sensitivity() and seq_aipe_r_sensitivity(); the script
# asserts, before the long run, that it reproduces those functions exactly
# for the same seed. It keeps its own loop only to record, replication by
# replication, what the sensitivity summary does not return: the final
# sample size and width (for the vignette's figures) and whether the
# replication ended with the lower bound n^-3 in place of a smaller
# positive variance estimate.
#
# Summarized per condition:
#   * coverage: the proportion of final intervals that contain the true
#     value. The guarantee is asymptotic: it approaches the nominal level
#     as the target width narrows.
#   * final interval width: at most the target width, by construction.
#   * final sample size, beside the sample size the known asymptotic
#     variance would require, N* = ceiling((2 z / width)^2 * xi^2).
#   * how often the variance estimate was unusable (zero, negative, or not
#     finite) at one or more checks, how often a replication ended with
#     the bound n^-3 in place of a smaller positive estimate, and how
#     often it ended at the pilot sample size.
#
# The designs of the vignette's method sections. The vignette illustrates
# each procedure with one design. So that it can state how that very design
# behaves, each is also run here, with the package's exported sensitivity
# function, 5,000 times (see "The designs of the vignette's method sections"
# below). Its Methods 1 and 2 use conditions of the study itself. The
# budget design is run twice: with the package's default, and with the
# article's rule that a population is closed the first time it meets its
# criterion (procedure = "first_passage").
#
# Seeds. The study is split into chunks of `chunk_size` replications.
# Chunk j of condition i (in the order of `design` below) is run under
# set.seed(seed + 1000 * i + j), so the results do not depend on the
# number of cores or on the order in which chunks finish, and the first k
# chunks of a condition are the same whatever the total number of chunks.
# Chunk j of the kth design of a method section is the sensitivity function
# called with seed = seed + 100000 + 1000 * k + j.
#
# Reproduce with (from the package's source directory):
#   R CMD INSTALL .                          # the version to be verified
#   SMSD_SIM_CORES=6 Rscript data-raw/simulation_study.R
# SMSD_SIM_CORES is the number of processes (default 1; forked, so more
# than 1 needs macOS or Linux). Set SMSD_SIM_CACHE to a directory to keep
# each finished chunk there, so that an interrupted run resumes instead of
# starting over (use a fresh directory after any change to the package or
# to this script).
#
# Run time. The run that produced the cached file, 10,000 replications per
# condition (20 chunks of 500) in six processes, took about 18 hours of
# processor time and about 4 hours of elapsed time. More than four fifths
# of it went to the six correlation conditions at width 0.08: a replication
# there checks the rule several hundred to about 2,300 times, on samples
# that grow to that many pairs. The designs of the method sections took
# about half an hour of processor time more. The measured times are stored
# in the file's `meta` element and are quoted from there by the vignette:
# `timing` has the processor and elapsed seconds of every chunk,
# `cpu_hours` their total processor time, `wall_minutes` the time during
# which one or more chunks were being computed (so a resumed run counts
# only its working time), `cores` the number of processes, `cpu` and
# `platform` the computer, and `tour_cpu_minutes` the processor time of the
# designs of the method sections.
#
# Ken Kelley <kkelley@nd.edu>
# ---------------------------------------------------------------------------

if (!"SMSD" %in% loadedNamespaces()) library(SMSD)

# ---- Controls -------------------------------------------------------------

chunk_size <- 500L     # replications per chunk
chunks     <- 20L      # chunks per condition
G          <- chunks * chunk_size   # replications per condition: 10,000
conf_level <- 0.95     # nominal confidence level of every interval
seed       <- 113L     # fixed for the study (see "Seeds" above)
max_n      <- 100000L  # cap on the sample size of a replication

cores     <- as.integer(Sys.getenv("SMSD_SIM_CORES", "1"))
cache_dir <- Sys.getenv("SMSD_SIM_CACHE", "")
stopifnot(!is.na(cores), cores >= 1L)

z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)

# ---- Data generating mechanisms -------------------------------------------
# Standardized draws (mean 0, variance 1). "skewed" is a chi square
# distribution with 3 degrees of freedom, standardized: skewness
# sqrt(8 / 3), about 1.63, and fourth standardized moment 7 (3 for the
# normal distribution).
rstd <- function(distribution) {
  if (distribution == "normal") stats::rnorm
  else function(k) (stats::rchisq(k, df = 3) - 3) / sqrt(6)
}
kurtosis <- function(distribution) if (distribution == "normal") 3 else 7

# Standardized mean difference: two independent groups with unit variance
# and means delta and 0, so the population value is delta exactly. Both
# groups have the same shape, so the third moment term of the asymptotic
# variance (Kelley, Darku, & Chattopadhyay, 2018, Equation 26 as corrected)
# vanishes, and
#   xi^2 = 2 + (delta^2 / 4) * (kappa - 1) / 2,
# with kappa the fourth standardized moment: 2 + delta^2 / 4 for normal
# data and 2 + 3 delta^2 / 4 for the skewed data.
xi2_smd <- function(delta, distribution) {
  2 + (delta^2 / 4) * (kurtosis(distribution) - 1) / 2
}

# Pearson correlation: x = z1 and y = rho z1 + sqrt(1 - rho^2) z2, with z1
# and z2 independent standardized draws, so the population correlation is
# rho exactly for either mechanism. The asymptotic variance of sqrt(n) r is
# (Kelley, Darku, & Chattopadhyay, 2019, Equation 3), for standardized
# variables,
#   xi^2 = rho^2 (mu40 + mu04 + 2 mu22) / 4 + mu22 - rho (mu31 + mu13),
# and the joint central moments of (x, y) follow from the independence of
# z1 and z2. For normal data the result is (1 - rho^2)^2.
xi2_cor <- function(rho, distribution) {
  kappa <- kurtosis(distribution)
  s2 <- 1 - rho^2
  mu40 <- kappa
  mu04 <- rho^4 * kappa + 6 * rho^2 * s2 + s2^2 * kappa
  mu22 <- rho^2 * kappa + s2
  mu31 <- rho * kappa
  mu13 <- rho^3 * kappa + 3 * rho * s2
  rho^2 * (mu40 + mu04 + 2 * mu22) / 4 + mu22 - rho * (mu31 + mu13)
}

# ---- The design -----------------------------------------------------------

design <- rbind(
  # Standardized mean difference: one moderate target across effect size
  # and shape.
  expand.grid(method = "Standardized mean difference",
              true_value = c(0.2, 0.5, 0.8),
              distribution = c("normal", "skewed"),
              width = 0.40,
              stringsAsFactors = FALSE),
  # Pearson correlation: the target width swept from wide to narrow for
  # normal data, and the narrowest width repeated with skewed data.
  expand.grid(method = "Pearson correlation",
              true_value = c(0.3, 0.5, 0.7),
              distribution = "normal",
              width = c(0.30, 0.15, 0.08),
              stringsAsFactors = FALSE),
  expand.grid(method = "Pearson correlation",
              true_value = c(0.3, 0.5, 0.7),
              distribution = "skewed",
              width = 0.08,
              stringsAsFactors = FALSE)
)
design <- design[, c("method", "distribution", "true_value", "width")]
is_smd <- design$method == "Standardized mean difference"
design$xi_squared <- vapply(seq_len(nrow(design)), function(i) {
  if (is_smd[i]) xi2_smd(design$true_value[i], design$distribution[i])
  else xi2_cor(design$true_value[i], design$distribution[i])
}, numeric(1L))
design$theoretical_n <- ceiling((2 * z_crit / design$width)^2 *
                                  design$xi_squared)
design$pilot_n <- vapply(seq_len(nrow(design)), function(i) {
  pilot <- if (is_smd[i]) seq_aipe_smd_pilot else seq_aipe_r_pilot
  pilot(width = design$width[i], conf_level = conf_level)
}, integer(1L))

# ---- One replication ------------------------------------------------------
# The closures of a condition: draw(data, k) appends k observations (or
# returns a fresh sample of k when data is NULL) and check(data) is the
# package's check function on the data so far. The order of the random
# draws is that of the package's sensitivity functions.
condition_closures <- function(i) {
  value <- design$true_value[i]
  width <- design$width[i]
  rand <- rstd(design$distribution[i])
  if (is_smd[i]) {
    draw <- function(data, k) {
      g1 <- rand(k) + value
      g2 <- rand(k)
      if (is.null(data)) list(g1 = g1, g2 = g2)
      else list(g1 = c(data$g1, g1), g2 = c(data$g2, g2))
    }
    check <- function(data) {
      seq_aipe_smd(data$g1, data$g2, width = width, conf_level = conf_level)
    }
  } else {
    draw <- function(data, k) {
      z1 <- rand(k)
      z2 <- rand(k)
      x <- z1
      y <- value * z1 + sqrt(1 - value^2) * z2
      if (is.null(data)) list(x = x, y = y)
      else list(x = c(data$x, x), y = c(data$y, y))
    }
    check <- function(data) {
      seq_aipe_r(data$x, data$y, width = width, conf_level = conf_level,
                 method = "pearson")
    }
  }
  list(draw = draw, check = check)
}

# The purely sequential procedure: the rule is checked at the pilot sample
# size and after every observation added, and the replication ends at the
# first check that the package reports as met (never at a check whose
# variance estimate is unusable, which the package reports as not met).
# Returns the final sample size (per group for the standardized mean
# difference), the full width of the final interval, and 0/1 indicators:
# the interval contains the true value; the cap, not the rule, ended the
# replication; the variance estimate was unusable at one or more checks;
# the replication ended with the bound n^-3 in place of a smaller positive
# estimate (the xi_squared row of the final check is the bound itself).
run_one <- function(pilot_n, draw, check, true_value) {
  data <- draw(NULL, pilot_n)
  n <- pilot_n
  unusable <- FALSE
  repeat {
    res <- check(data)
    if (isTRUE(attr(res, "variance_unusable"))) unusable <- TRUE
    stopped <- isTRUE(attr(res, "is_satisfied"))
    if (stopped || n >= max_n) break
    data <- draw(data, 1L)
    n <- n + 1L
  }
  v <- stats::setNames(res$value, res$term)
  lower <- v[["lower_limit"]]
  upper <- v[["upper_limit"]]
  c(n = n,
    ci_width = upper - lower,
    covered = as.numeric(isTRUE(lower <= true_value && true_value <= upper)),
    capped = as.numeric(!stopped),
    unusable = as.numeric(unusable),
    bound = as.numeric(stopped && v[["xi_squared"]] <= n^(-3) * (1 + 1e-12)))
}

run_replications <- function(i, replications, chunk_seed) {
  f <- condition_closures(i)
  set.seed(chunk_seed)
  vapply(seq_len(replications), function(r) {
    run_one(design$pilot_n[i], f$draw, f$check, design$true_value[i])
  }, numeric(6L))
}

chunk_seed <- function(i, j) seed + 1000L * i + j

# ---- This loop is the package's own engine --------------------------------
# For one condition of each kind, the first replications of chunk 1 must
# equal what the exported sensitivity function returns for the same seed.
# If this fails, the sensitivity functions and this script have diverged
# (in the order of the random draws or in the stopping logic), and the
# script must be brought back in line before its results are described as
# those of the shipped procedures.
local({
  agrees <- function(i, sensitivity, replications = 25L) {
    mine <- run_replications(i, replications, chunk_seed(i, 1L))
    theirs <- sensitivity(true_value = design$true_value[i],
                          width = design$width[i], conf_level = conf_level,
                          distribution = design$distribution[i],
                          G = replications, seed = chunk_seed(i, 1L),
                          max_n = max_n)
    theirs <- stats::setNames(theirs$value, theirs$term)
    isTRUE(all.equal(
      c(mean(mine["n", ]), stats::sd(mine["n", ]), mean(mine["ci_width", ]),
        mean(mine["covered", ]), 100 * mean(mine["unusable", ])),
      unname(theirs[c("mean_n", "sd_n", "mean_width", "coverage",
                      "pct_variance_unusable")]),
      tolerance = 1e-12))
  }
  i_smd <- which(is_smd & design$distribution == "skewed" &
                   design$true_value == 0.5)
  i_cor <- which(!is_smd & design$distribution == "normal" &
                   design$true_value == 0.7 & design$width == 0.30)
  if (!agrees(i_smd, seq_aipe_smd_sensitivity)) {
    stop("The standardized mean difference loop no longer reproduces ",
         "seq_aipe_smd_sensitivity().")
  }
  if (!agrees(i_cor, seq_aipe_r_sensitivity)) {
    stop("The correlation loop no longer reproduces seq_aipe_r_sensitivity().")
  }
})

# ---- Run the design -------------------------------------------------------

tasks <- expand.grid(condition = seq_len(nrow(design)),
                     chunk = seq_len(chunks))
# Chunk by chunk across all conditions, the costliest condition first
# within a round, so that the processes stay evenly loaded.
tasks <- tasks[order(tasks$chunk, -design$theoretical_n[tasks$condition]), ]

run_task <- function(k) {
  i <- tasks$condition[k]
  j <- tasks$chunk[k]
  file <- if (nzchar(cache_dir)) {
    file.path(cache_dir, sprintf("condition%02d_chunk%02d.rds", i, j))
  } else ""
  if (nzchar(file) && file.exists(file)) return(readRDS(file))
  started <- Sys.time()
  cpu <- proc.time()
  mat <- run_replications(i, chunk_size, chunk_seed(i, j))
  cpu <- proc.time() - cpu
  out <- list(condition = i, chunk = j, seed = chunk_seed(i, j), mat = mat,
              started = started, finished = Sys.time(),
              cpu_seconds = sum(cpu[c("user.self", "sys.self")]),
              cores = cores)
  if (nzchar(file)) {
    tmp <- paste0(file, ".", Sys.getpid(), ".tmp")
    saveRDS(out, tmp)
    file.rename(tmp, file)
  }
  out
}

if (nzchar(cache_dir)) dir.create(cache_dir, recursive = TRUE,
                                  showWarnings = FALSE)
message("Running ", nrow(design), " conditions, ", G,
        " replications each, on ", cores, " core(s) ...")
done <- if (cores > 1L) {
  parallel::mclapply(seq_len(nrow(tasks)), run_task, mc.cores = cores,
                     mc.preschedule = FALSE)
} else {
  lapply(seq_len(nrow(tasks)), run_task)
}
failed <- !vapply(done, is.list, logical(1L))
if (any(failed)) {
  stop(sum(failed), " chunk(s) failed; the first error was: ",
       paste(format(done[[which(failed)[1L]]]), collapse = " "))
}

# ---- The designs of the vignette's method sections ------------------------
# One call of the matching sensitivity function per design, as the vignette
# makes it, without G and seed, which are added chunk by chunk. The designs
# of the vignette's Methods 1 and 2 are conditions of the study above.
tour_chunks <- 10L
tour_G <- tour_chunks * chunk_size   # replications per design: 5,000
tour_calls <- list(
  R2 = quote(seq_aipe_R2_sensitivity(true_value = 0.5, width = 0.10, p = 3)),
  cv = quote(seq_aipe_cv_sensitivity(true_value = 0.2, width = 0.04,
                                     distribution = "gamma")),
  reg_coef = quote(seq_aipe_reg_coef_sensitivity(true_value = 0.5,
                                                 width = 0.2)),
  sm = quote(seq_aipe_sm_sensitivity(true_value = 0.4, width = 0.40)),
  contrast = quote(seq_aipe_c_sensitivity(
    true_means = c(5.0, 4.6, 4.4), contrast_weights = c(1, -0.5, -0.5),
    width = 0.40)),
  gini = quote(seq_aipe_gini_sensitivity(true_value = 0.3, width = 0.10)),
  budget = quote(seq_budget_sensitivity(
    true_means = c(74, 75, 75), true_sds = c(8, 8, 8),
    contrast_weights = c(1, -0.5, -0.5), cost_per_group = c(20, 40, 40),
    budget = 6000, margin = 4, m_0 = 10)),
  mr_cv = quote(seq_mr_cv_sensitivity(true_value = 0.4, A = 5e6, c = 10,
                                      distribution = "gamma")),
  mr_smd = quote(seq_mr_smd_sensitivity(true_value = 0.5, A = 10000,
                                        c = 2.4)),
  # The budget design again, under the article's rule that a population is
  # closed the first time it meets its criterion.
  budget_first_passage = quote(seq_budget_sensitivity(
    true_means = c(74, 75, 75), true_sds = c(8, 8, 8),
    contrast_weights = c(1, -0.5, -0.5), cost_per_group = c(20, 40, 40),
    budget = 6000, margin = 4, m_0 = 10, procedure = "first_passage"))
)
tour_seed <- function(k, j) seed + 100000L + 1000L * k + j

tour_tasks <- expand.grid(design = seq_along(tour_calls),
                          chunk = seq_len(tour_chunks))
run_tour_task <- function(row) {
  k <- tour_tasks$design[row]
  j <- tour_tasks$chunk[row]
  file <- if (nzchar(cache_dir)) {
    file.path(cache_dir, sprintf("tour%02d_chunk%02d.rds", k, j))
  } else ""
  if (nzchar(file) && file.exists(file)) return(readRDS(file))
  cpu <- proc.time()
  call <- as.call(c(as.list(tour_calls[[k]]),
                    list(G = chunk_size, seed = tour_seed(k, j))))
  res <- eval(call)
  cpu <- proc.time() - cpu
  out <- list(design = k, chunk = j, seed = tour_seed(k, j),
              value = stats::setNames(res$value, res$term),
              cpu_seconds = sum(cpu[c("user.self", "sys.self")]))
  if (nzchar(file)) {
    tmp <- paste0(file, ".", Sys.getpid(), ".tmp")
    saveRDS(out, tmp)
    file.rename(tmp, file)
  }
  out
}
message("Running the ", length(tour_calls), " designs of the vignette's ",
        "method sections, ", tour_G, " replications each ...")
tour_done <- if (cores > 1L) {
  parallel::mclapply(seq_len(nrow(tour_tasks)), run_tour_task,
                     mc.cores = cores, mc.preschedule = FALSE)
} else {
  lapply(seq_len(nrow(tour_tasks)), run_tour_task)
}
failed <- !vapply(tour_done, is.list, logical(1L))
if (any(failed)) {
  stop(sum(failed), " chunk(s) of the method section designs failed; the ",
       "first error was: ",
       paste(format(tour_done[[which(failed)[1L]]]), collapse = " "))
}

# One row per design and summary term: the mean over the chunks, which for
# equal chunks is the value for all `tour_G` replications whenever the term
# is a mean, a proportion, a percentage, or a constant of the design.
# Medians, standard deviations, standard errors, and the variance of the
# estimates across replications are not, and are left out.
tour_summary <- do.call(rbind, lapply(seq_along(tour_calls), function(k) {
  mine <- tour_done[tour_tasks$design == k]
  values <- vapply(mine, function(d) d$value, mine[[1L]]$value)
  keep <- !grepl("^(median|sd|se|empirical)_", rownames(values))
  value <- rowMeans(values[keep, , drop = FALSE])
  if ("replications" %in% names(value)) value[["replications"]] <- tour_G
  data.frame(design = names(tour_calls)[k], term = names(value),
             value = unname(value), stringsAsFactors = FALSE)
}))
tour <- list(
  summary = tour_summary,
  calls = vapply(tour_calls, function(cl) {
    paste(trimws(deparse(cl, width.cutoff = 500L)), collapse = " ")
  }, character(1L)),
  G = tour_G,
  seed_rule = paste("chunk j of design k: the call with G = chunk_size and",
                    "seed = seed + 100000 + 1000 * k + j")
)

# ---- Summarize ------------------------------------------------------------

condition_matrix <- function(i) {
  mine <- done[tasks$condition == i]
  mine <- mine[order(vapply(mine, function(d) d$chunk, integer(1L)))]
  do.call(cbind, lapply(mine, function(d) d$mat))
}

summarize_condition <- function(i) {
  mat <- condition_matrix(i)
  stopifnot(ncol(mat) == G)
  n <- mat["n", ]
  w <- mat["ci_width", ]
  data.frame(
    method                = design$method[i],
    distribution          = design$distribution[i],
    true_value            = design$true_value[i],
    width                 = design$width[i],
    conf_level            = conf_level,
    G                     = G,
    pilot_n               = design$pilot_n[i],
    xi_squared            = design$xi_squared[i],
    theoretical_n         = design$theoretical_n[i],
    mean_n                = mean(n),
    median_n              = stats::median(n),
    sd_n                  = stats::sd(n),
    mean_width            = mean(w, na.rm = TRUE),
    median_width          = stats::median(w, na.rm = TRUE),
    pct_width_le_omega    = 100 * sum(w <= design$width[i], na.rm = TRUE) / G,
    coverage              = sum(mat["covered", ]) / G,
    pct_variance_unusable = 100 * sum(mat["unusable", ]) / G,
    pct_bound_stop        = 100 * sum(mat["bound", ]) / G,
    pct_stop_at_pilot     = 100 * sum(n == design$pilot_n[i] &
                                        mat["capped", ] == 0) / G,
    pct_hit_max_n         = 100 * sum(mat["capped", ]) / G,
    stringsAsFactors      = FALSE
  )
}

summary_df <- do.call(rbind, lapply(seq_len(nrow(design)),
                                    summarize_condition))
rownames(summary_df) <- NULL

# Replication by replication results for four conditions, for the figures.
spotlight_of <- function(i, label) {
  mat <- condition_matrix(i)
  list(n = as.integer(mat["n", ]), ci_width = mat["ci_width", ],
       covered = mat["covered", ] == 1, unusable = mat["unusable", ] == 1,
       bound = mat["bound", ] == 1,
       true_value = design$true_value[i], width = design$width[i],
       pilot_n = design$pilot_n[i],
       theoretical_n = design$theoretical_n[i], label = label)
}
find_condition <- function(smd, distribution, true_value, width) {
  which(is_smd == smd & design$distribution == distribution &
          design$true_value == true_value & design$width == width)
}
spotlight <- list(
  smd_normal = spotlight_of(find_condition(TRUE, "normal", 0.5, 0.40),
                            "SMD, delta = 0.5, normal"),
  smd_skewed = spotlight_of(find_condition(TRUE, "skewed", 0.5, 0.40),
                            "SMD, delta = 0.5, skewed"),
  cor_wide   = spotlight_of(find_condition(FALSE, "normal", 0.7, 0.30),
                            "Pearson r, rho = 0.7, width = 0.3"),
  cor_fine   = spotlight_of(find_condition(FALSE, "normal", 0.7, 0.08),
                            "Pearson r, rho = 0.7, width = 0.08")
)

started  <- do.call(c, lapply(done, function(d) d$started))
finished <- do.call(c, lapply(done, function(d) d$finished))

# Minutes during which one or more chunks were being computed: the length
# of the union of the chunks' intervals. For an uninterrupted run this is
# the time from the first start to the last finish; a resumed run does not
# count the time between its sessions.
busy_minutes <- function(started, finished) {
  from <- as.numeric(started)
  to <- as.numeric(finished)
  o <- order(from)
  from <- from[o]
  to <- to[o]
  total <- 0
  end <- from[1L]
  for (k in seq_along(from)) {
    if (to[k] > end) {
      total <- total + to[k] - max(from[k], end)
      end <- to[k]
    }
  }
  total / 60
}

# The processor, where the system reports one (macOS); NA elsewhere.
cpu_name <- function() {
  out <- tryCatch(
    suppressWarnings(system2("sysctl", c("-n", "machdep.cpu.brand_string"),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0L))
  if (length(out) == 1L && nzchar(out)) out else NA_character_
}

timing <- data.frame(
  condition       = vapply(done, function(d) d$condition, integer(1L)),
  chunk           = vapply(done, function(d) d$chunk, integer(1L)),
  cpu_seconds     = vapply(done, function(d) d$cpu_seconds, numeric(1L)),
  elapsed_seconds = as.numeric(finished) - as.numeric(started)
)
timing <- timing[order(timing$condition, timing$chunk), ]
rownames(timing) <- NULL
results <- list(
  summary   = summary_df,
  spotlight = spotlight,
  tour      = tour,
  meta = list(
    G            = G,
    chunk_size   = chunk_size,
    chunks       = chunks,
    conf_level   = conf_level,
    seed         = seed,
    seed_rule    = "chunk j of condition i: set.seed(seed + 1000 * i + j)",
    max_n        = max_n,
    generated_on = as.character(as.Date(max(finished))),
    r_version    = R.version.string,
    platform     = R.version$platform,
    cpu          = cpu_name(),
    smsd_version = as.character(utils::packageVersion("SMSD")),
    cores        = max(vapply(done, function(d) d$cores, integer(1L))),
    # For the 18 conditions: the time during which chunks were being
    # computed, the processor time summed over chunks, and both for every
    # chunk.
    wall_minutes = busy_minutes(started, finished),
    cpu_hours    = sum(timing$cpu_seconds) / 3600,
    timing       = timing,
    # For the designs of the method sections: processor time.
    tour_cpu_minutes = sum(vapply(tour_done, function(d) d$cpu_seconds,
                                  numeric(1L))) / 60
  )
)

out_dir <- file.path("inst", "extdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "aipe_simulation_results.rds")
saveRDS(results, out_file, compress = "xz")
message("Wrote ", out_file)
print(summary_df, digits = 4)
print(tour_summary[tour_summary$term %in%
                     c("mean_n", "coverage", "ratio_n", "risk_efficiency",
                       "noninferiority_rejection_rate"), ],
      digits = 4, row.names = FALSE)
str(results$meta)
