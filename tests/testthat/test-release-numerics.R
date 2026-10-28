test_that("power planning is invariant to representable standardized units", {
  expect_identical(power_z(-1e308, 1e308, 1e308, beta = 0.2), 2L)
  expect_equal(power_z(-1e308, 1e308, 1e308, n = 1),
               stats::pnorm(stats::qnorm(0.95) - 2))
  expect_identical(power_z(-1e308, 1e308, 1e308, n_min = 1,
                           n_max = 5, n_values = 5),
                   power_z(-1, 1, 1, n_min = 1, n_max = 5, n_values = 5))
})

test_that("extreme valid confidence levels have finite quantiles and pilots", {
  level <- 1 - .Machine$double.eps / 2
  z <- stats::qnorm((1 - level) / 2, lower.tail = FALSE)
  expect_true(is.finite(z))
  expect_identical(seq_fwd_norm_mean_optimal(1, 1, level), 69L)
  expect_identical(seq_aipe_R2_pilot(0.2, p = 2, conf_level = level), 166L)
  expect_true(is.finite(seq_fwd_norm_mean_pilot(1, level)))
  expect_true(is.finite(seq_fwd_dfr_mean_pilot(1, level)))
  expect_true(is.finite(ts_fwd_norm_mean_pilot(1, level)))
  expect_true(is.finite(ts_fwd_dfr_mean_pilot(1, level)))
  expect_true(is.finite(seq_blci_dfr_median_pilot(1, level)))
})

test_that("rank variance calculations match the pairwise formulas with ties", {
  withr::local_seed(20261026)
  for (n in c(4L, 7L, 20L, 51L)) {
    for (tied in c(FALSE, TRUE)) {
      x <- stats::rnorm(n)
      y <- stats::rnorm(n)
      if (tied) {
        x <- round(x)
        y <- round(y)
      }
      rx <- rank(x)
      ry <- rank(y)
      count <- vapply(seq_len(n), function(i) {
        sum(rx <= rx[i] & ry <= ry[i])
      }, numeric(1))
      expect_identical(.rank_lower_count(rx, ry), count)
      upper <- vapply(seq_len(n), function(i) sum(ry[rx >= rx[i]]),
                       numeric(1))
      expect_identical(.rank_upper_sum(rx, ry), upper)
      w <- 2 * (n + 1) * count - n * (rx + ry)
      v <- n * rx * ry + (n + 1) *
        colSums(t(outer(rx, rx, "<=")) * ry +
                  t(outer(ry, ry, "<=")) * rx)
      expect_equal(.var_r(x, y, "kendall", floor = FALSE),
                   16 * stats::var(w) / (n * (n + 1))^2)
      expect_equal(.var_r(x, y, "spearman", floor = FALSE),
                   144 * stats::var(v) / (n * (n + 1)^2)^2)
    }
  }
})

test_that("AIPE results distinguish width satisfaction from coverage", {
  withr::local_seed(20261026)
  x <- stats::rnorm(50)
  y <- 0.8 * x + 0.6 * stats::rnorm(50)
  out <- seq_aipe_r(x, y, width = 0.4)
  expect_identical(attr(out, "coverage_guarantee"), "asymptotic")
  expect_true(any(grepl("Coverage        : asymptotic", format(out), fixed = TRUE)))
  expect_match(as.character(seq_report(out)),
                "meeting the width target does not establish nominal coverage", fixed = TRUE)
  r2 <- seq_aipe_R2(cbind(x, stats::rnorm(50)), y, width = 0.3)
  expect_identical(attr(r2, "coverage_guarantee"), "asymptotic")
})

test_that("coverage tables agree with the individual simulated studies", {
  study <- readRDS(system.file("extdata", "coverage_study.rds", package = "SMSD"))
  expect_equal(nrow(study$summary), 10)
  expect_true(all(study$summary$G >= 1000))
  expect_true(all(nzchar(study$meta$function_md5)))
  for (i in seq_len(nrow(study$summary))) {
    design <- study$summary[i, ]
    trials <- study$trials[study$trials$design_id == design$id, ]
    expect_equal(nrow(trials), design$G)
    expect_identical(trials$covered,
                     trials$lower <= design$truth & design$truth <= trials$upper)
    expect_equal(design$coverage, mean(trials$covered))
    expect_equal(design$mean_n, mean(trials$n))
    expect_equal(design$coverage_mcse,
                  sqrt(design$coverage * (1 - design$coverage) / design$G))
    expect_equal(design$pct_stop_at_pilot,
                  100 * mean(trials$stopped & trials$n == design$pilot))
    expect_equal(design$pct_capped, 100 * mean(!trials$stopped))
    expect_equal(design$pct_width_within,
                  100 * mean(is.finite(trials$width) & trials$width <= design$width))
  }
})

test_that("paper replication summaries agree with saved batch outputs", {
  study <- readRDS(system.file("extdata", "paper_replications.rds", package = "SMSD"))
  quantities <- c("mean_n", "mean_n1", "mean_n", "mean_n", "mean_n_1")
  expect_identical(study$summary$year, c(2016L, 2017L, 2018L, 2019L, 2025L))
  for (i in seq_len(nrow(study$summary))) {
    row <- study$summary[i, ]
    batch <- study$batches[[as.character(row$year)]]
    means <- vapply(batch$outputs, function(out) {
      out$value[match(quantities[i], out$term)]
    }, 0)
    expect_equal(length(unique(batch$seeds)), study$meta$batches)
    expect_equal(row$studies, study$meta$batches * study$meta$per_batch)
    expect_equal(row$reproduced, mean(means))
    expect_equal(row$mcse, stats::sd(means) / sqrt(length(means)))
    expect_equal(row$difference, row$reproduced - row$published)
  }
})
