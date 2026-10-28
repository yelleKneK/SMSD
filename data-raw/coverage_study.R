# Operating characteristics of two difficult sequential interval designs.
# Run from the package root after loading the source with pkgload::load_all().
# Every replication calls the exported check on each additional observation.
# Builds read the cached result; they never run this simulation.
# SMSD_COVERAGE_G sets replications per condition (default 2000).
# SMSD_COVERAGE_OUT overrides the result file below.
if (!"SMSD" %in% loadedNamespaces()) library(SMSD)
G <- as.integer(Sys.getenv("SMSD_COVERAGE_G", "2000"))
stopifnot(length(G) == 1L, !is.na(G), G >= 1000L)
destination <- Sys.getenv("SMSD_COVERAGE_OUT",
                          "inst/extdata/coverage_study.rds")
started <- Sys.time()
term <- function(out, name) out$value[match(name, out$term)]
designs <- rbind(
  data.frame(procedure = "Pearson", truth = 0.8, p = 1L,
             width = c(0.4, 0.2, 0.4, 0.2),
             m_0 = c(4L, 4L, 50L, 50L)),
  data.frame(procedure = "Multiple R2", truth = 0, p = 20L,
             width = c(0.3, 0.2, 0.3, 0.2),
             m_0 = c(22L, 22L, 400L, 400L)),
  data.frame(procedure = "Multiple R2", truth = 0.5, p = 5L,
             width = c(0.3, 0.2), m_0 = c(7L, 7L))
)
designs$id <- seq_len(nrow(designs))
designs$seed <- 20261026L + 1000L * designs$id
designs$G <- G
designs$conf_level <- 0.95
designs$max_n <- 5000L

one_design <- function(a) {
  set.seed(a$seed)
  is_pearson <- a$procedure == "Pearson"
  pilot <- if (is_pearson) {
    seq_aipe_r_pilot(a$width, m_0 = a$m_0)
  } else {
    seq_aipe_R2_pilot(a$width, p = a$p, m_0 = a$m_0)
  }
  # Generate a whole stream once per study and reveal only the prefix at
  # a check. This changes no stopping rule and supplies reproducible data.
  draw <- function() {
    x <- matrix(stats::rnorm(a$max_n * a$p), ncol = a$p)
    beta <- if (is_pearson) a$truth else sqrt(a$truth / a$p)
    error_sd <- if (is_pearson) sqrt(1 - a$truth^2) else sqrt(1 - a$truth)
    y <- beta * rowSums(x) + error_sd * stats::rnorm(a$max_n)
    list(x = x, y = y)
  }
  runs <- lapply(seq_len(a$G), function(i) {
    data <- draw()
    unusable <- FALSE
    noted <- FALSE
    for (n in seq.int(pilot, a$max_n)) {
      index <- seq_len(n)
      out <- if (is_pearson) {
        seq_aipe_r(data$x[index, 1L], data$y[index], a$width)
      } else {
        seq_aipe_R2(data$x[index, , drop = FALSE], data$y[index], a$width)
      }
      unusable <- unusable || isTRUE(attr(out, "variance_unusable"))
      noted <- noted || length(attr(out, "notes")) > 0L
      stopped <- isTRUE(attr(out, "is_satisfied"))
      if (stopped) break
    }
    lower <- term(out, "lower_limit")
    upper <- term(out, "upper_limit")
    data.frame(design_id = a$id, replication = i, n = n,
      estimate = term(out, if (is_pearson) "r" else "R2"),
      lower = lower, upper = upper, width = upper - lower,
      covered = isTRUE(lower <= a$truth && a$truth <= upper),
      stopped = stopped, variance_unusable = unusable, noted = noted)
  })
  trials <- do.call(rbind, runs)
  coverage <- mean(trials$covered)
  ci <- stats::binom.test(sum(trials$covered), nrow(trials))$conf.int
  summary <- cbind(a, pilot = pilot, mean_n = mean(trials$n),
    se_mean_n = stats::sd(trials$n) / sqrt(a$G),
    median_n = stats::median(trials$n), p90_n = unname(stats::quantile(trials$n, 0.9)),
    coverage = coverage, coverage_mcse = sqrt(coverage * (1 - coverage) / a$G),
    coverage_lower = ci[1L], coverage_upper = ci[2L],
    mean_width = mean(trials$width, na.rm = TRUE),
    pct_width_within = 100 * mean(is.finite(trials$width) & trials$width <= a$width),
    pct_capped = 100 * mean(!trials$stopped),
    pct_variance_unusable = 100 * mean(trials$variance_unusable),
    pct_stop_at_pilot = 100 * mean(trials$stopped & trials$n == pilot))
  cat(sprintf("Design %d/%d: %s, truth %g, width %g, pilot %d: coverage %.4f, mean n %.1f\n",
    a$id, nrow(designs), a$procedure, a$truth, a$width, pilot, coverage, mean(trials$n)))
  flush.console()
  list(summary = summary, trials = trials)
}

results <- lapply(seq_len(nrow(designs)), function(i) one_design(designs[i, ]))
function_names <- c("seq_aipe_r", "seq_aipe_R2", ".var_r", ".R2_rule",
                    ".aipe_result", ".mu4", ".mu22", ".mu31", ".mu13")
fingerprint <- function(name) {
  f <- get(name, asNamespace("SMSD"))
  file <- tempfile()
  on.exit(unlink(file))
  writeLines(c(deparse(formals(f)), deparse(body(f), width.cutoff = 500L)), file)
  unname(tools::md5sum(file))
}
record <- list(summary = do.call(rbind, lapply(results, `[[`, "summary")),
  trials = do.call(rbind, lapply(results, `[[`, "trials")),
  meta = list(generated_on = as.character(Sys.Date()),
    package_version = as.character(utils::packageVersion("SMSD")),
    R = R.version.string, platform = R.version$platform,
    generator = "independent multivariate normal streams; one observation per stage",
    function_md5 = stats::setNames(vapply(function_names, fingerprint, ""), function_names),
    minutes = as.numeric(difftime(Sys.time(), started, units = "mins"))))
dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
saveRDS(record, destination, compress = "xz")
print(record$summary, row.names = FALSE)
cat("Saved", destination, "\n")
