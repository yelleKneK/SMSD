# Fresh Monte Carlo replications of one published simulation row per paper.
# Run from the package root after loading SMSD. Equal independent batches
# provide a Monte Carlo standard error of each reproduced mean.
if (!"SMSD" %in% loadedNamespaces()) library(SMSD)
batches <- 20L
per_batch <- 25L
term <- function(out, name) out$value[match(name, out$term)]
designs <- list(
  list(year = 2016L, label = "CV minimum risk", source = "Table 1, p. 642",
    published = 123.8528, published_se = 0.0870, quantity = "mean_n",
    run = function(seed) seq_mr_cv_sensitivity(0.2, A = 5e6, c = 10,
      distribution = "normal", variance = "published", G = per_batch, seed = seed)),
  list(year = 2017L, label = "SMD minimum risk", source = "Table 3, p. 104",
    published = 50.0210, published_se = 0.0020, quantity = "mean_n1",
    run = function(seed) seq_mr_smd_sensitivity(0, A = 1e6, c = 500,
      G = per_batch, seed = seed)),
  list(year = 2018L, label = "SMD sequential AIPE", source = "Table 2, p. 233",
    published = 793.6, published_se = NA_real_, quantity = "mean_n",
    run = function(seed) seq_aipe_smd_sensitivity(0.5, width = 0.2,
      G = per_batch, seed = seed)),
  list(year = 2019L, label = "Pearson sequential AIPE", source = "Table 2, p. 498",
    published = 204.0, published_se = 0.7435, quantity = "mean_n",
    run = function(seed) seq_aipe_r_sensitivity(0.5, width = 0.2,
      G = per_batch, seed = seed)),
  list(year = 2025L, label = "Budget contrast, group 1", source = "Table 1, p. 431",
    published = 95.868, published_se = 6.7458 / sqrt(5000), quantity = "mean_n_1",
    run = function(seed) seq_budget_sensitivity(c(1, 2, 3), c(4.5, 0.9, 1),
      contrast_weights = c(1, -2, 1), cost_per_group = c(38, 15, 10),
      budget = 5000, m_0 = 10, procedure = "first_passage",
      G = per_batch, seed = seed))
)
runs <- lapply(seq_along(designs), function(i) {
  d <- designs[[i]]
  seeds <- 20262026L + 1000L * i + seq_len(batches)
  outputs <- lapply(seeds, d$run)
  estimates <- vapply(outputs, term, 0, name = d$quantity)
  se <- stats::sd(estimates) / sqrt(batches)
  row <- data.frame(year = d$year, procedure = d$label, source = d$source,
    studies = batches * per_batch, published = d$published,
    reproduced = mean(estimates), mcse = se, published_se = d$published_se,
    difference = mean(estimates) - d$published)
  cat(sprintf("%d: article %.4f; reproduced %.4f, MC SE %.4f\n",
    d$year, d$published, row$reproduced, se))
  flush.console()
  list(summary = row, details = list(seeds = seeds, outputs = outputs))
})
details <- stats::setNames(lapply(runs, `[[`, "details"),
                           vapply(designs, function(d) as.character(d$year), ""))
result <- list(summary = do.call(rbind, lapply(runs, `[[`, "summary")), batches = details,
  meta = list(generated_on = as.character(Sys.Date()),
    package_version = as.character(utils::packageVersion("SMSD")),
    R = R.version.string, batches = batches, per_batch = per_batch,
    source_md5 = tools::md5sum(list.files("R", full.names = TRUE, pattern = "[.]R$"))))
saveRDS(result, "inst/extdata/paper_replications.rds", compress = "xz")
print(result$summary, row.names = FALSE)
