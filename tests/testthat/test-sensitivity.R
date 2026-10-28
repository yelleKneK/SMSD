# Tests for the seq_aipe_*_sensitivity Monte Carlo functions and the engine
# they share (R/sensitivity.R, R/sensitivity_internals.R).

# ---- helpers ---------------------------------------------------------------

sens_terms <- c("true_value", "target_width", "conf_level", "replications",
                "mean_n", "median_n", "sd_n", "mean_width", "median_width",
                "pct_width_within", "coverage", "pct_hit_max_n",
                "pct_variance_unusable")

sens_value <- function(out, term) out$value[match(term, out$term)]

# A check object for driving the engine with scripted answers. The "data" of
# a scripted study is just its sample size (see stub_draw).
stub_check <- function(n, satisfied, unusable = FALSE, lower = -1,
                       upper = 1) {
  out <- new_smsd_check(
    term = c("lower_limit", "estimate", "upper_limit"),
    value = c(-1, 0, 1),
    method = "scripted", is_satisfied = satisfied, current_n = n,
    criterion_n = n, conf_level = 0.95, width = 2,
    extra = list(variance_unusable = unusable)
  )
  # Set after construction: the constructor refuses Inf and NaN, and one
  # test scripts a check that holds NaN.
  out$value <- c(lower, (lower + upper) / 2, upper)
  out
}
stub_draw <- function(data, k) if (is.null(data)) k else data + k

# Run the engine on a scripted study: one replication unless G is given.
stub_run <- function(evaluate, pilot_n = 10L, max_n = 1000L, batch = 1L,
                     G = 1L, draw = stub_draw) {
  .sensitivity_run(G = G, seed = NULL, pilot_n = pilot_n, draw = draw,
                   evaluate = evaluate, true_value = 0, width = 2,
                   conf_level = 0.95, max_n = max_n, batch = batch)
}

# The sequential procedure as the help page describes it, written without the
# engine and without the functions' own generators: draw the pilot sample,
# check the stopping rule, add `batch` observations (to every group or
# variable), stop when the rule is met or at max_n, and summarize. gen(k)
# draws k new observations with base R generators and returns a list with one
# vector per group or variable; check(d) calls the check function that the
# sensitivity function is documented to run. Returns the nine computed rows
# of the summary, the final sample sizes, the sample size each final check
# reported, and the width of each final interval (NA where the final check
# gave no interval).
hand_study <- function(G, seed, pilot, gen, check, truth, width,
                       max_n = 100000L, batch = 1L) {
  one <- function(i) {
    d <- gen(pilot)
    n <- pilot
    unusable <- FALSE
    repeat {
      res <- check(d)
      unusable <- unusable || isTRUE(attr(res, "variance_unusable"))
      met <- attr(res, "is_satisfied")
      if (met || n >= max_n) break
      k <- min(batch, max_n - n)
      d <- Map(c, d, gen(k))
      n <- n + k
    }
    lower <- res$value[res$term == "lower_limit"]
    upper <- res$value[res$term == "upper_limit"]
    c(n = n, width = upper - lower,
      covered = isTRUE(lower <= truth && truth <= upper), capped = !met,
      unusable = unusable, current_n = attr(res, "current_n"))
  }
  m <- withr::with_seed(seed, vapply(seq_len(G), one, numeric(6)))
  w <- m["width", ]
  ok <- is.finite(w)
  list(
    pilot = pilot, n = unname(m["n", ]), current_n = unname(m["current_n", ]),
    widths = unname(w),
    rows = c(mean_n = mean(m["n", ]), median_n = stats::median(m["n", ]),
             sd_n = stats::sd(m["n", ]),
             mean_width = if (any(ok)) mean(w[ok]) else NA_real_,
             median_width = if (any(ok)) stats::median(w[ok]) else NA_real_,
             pct_width_within = 100 * sum(ok & w <= width) / G,
             coverage = sum(m["covered", ]) / G,
             pct_hit_max_n = 100 * sum(m["capped", ]) / G,
             pct_variance_unusable = 100 * sum(m["unusable", ]) / G)
  )
}

# The two unit variance generators of the help page, in base R: standard
# normal, and chi square with 3 degrees of freedom standardized to mean 0 and
# variance 1.
hand_rand <- list(
  normal = function(k) stats::rnorm(k),
  skewed = function(k) (stats::rchisq(k, df = 3) - 3) / sqrt(6)
)

# A sensitivity summary must equal the hand-written study row for row, and
# the comparison must be informative: sampling went past the pilot (so the
# batch size and the confidence level matter), the sample size reported is
# the one the final check counted, and, where a cap was set, it stopped some
# replications but not all.
#
# default_pilot is the pilot sample size at the default confidence level
# (same width and m_0), and batch the batch size of the study. A function
# whose pilot ignored conf_level would start at default_pilot; if that
# differed from the right pilot by a multiple of the batch, it would visit
# the same sample sizes over the same random stream and could give the same
# study. Where both are supplied, the difference must not be such a multiple.
expect_same_study <- function(out, ref, label, cap_binds = TRUE,
                              default_pilot = NULL, batch = NULL) {
  rows <- names(ref$rows)
  got <- sens_value(out, rows)
  names(got) <- rows
  expect_identical(got, ref$rows, info = label)
  expect_gt(sum(ref$n > ref$pilot), 0, label = paste(label, "past the pilot"))
  expect_identical(ref$n, ref$current_n, info = label)
  if (!is.null(default_pilot)) {
    expect_true((ref$pilot - default_pilot) %% batch != 0,
                label = paste(label, "pilot reveals the confidence level"))
  }
  capped <- ref$rows[["pct_hit_max_n"]]
  if (cap_binds) {
    expect_true(capped > 0 && capped < 100,
                label = paste(label, "cap stops some replications, not all"))
  } else {
    expect_identical(capped, 0, info = label)
  }
}

# Sample skewness (third standardized central moment).
sens_skewness <- function(x) {
  d <- x - mean(x)
  mean(d^3) / mean(d^2)^1.5
}

# The rows named in the \value section of the help page, in order.
sens_documented_rows <- function() {
  man <- testthat::test_path("..", "..", "man", "seq_aipe_sensitivity.Rd")
  rd <- if (file.exists(man)) {
    tools::parse_Rd(man)
  } else {
    tryCatch(tools::Rd_db("SMSD")[["seq_aipe_sensitivity.Rd"]],
             error = function(e) NULL)
  }
  if (is.null(rd)) return(NULL)
  value <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\value"), rd)
  rows <- character()
  walk <- function(x) {
    if (identical(attr(x, "Rd_tag"), "\\item") && length(x) == 2L) {
      codes <- Filter(function(y) identical(attr(y, "Rd_tag"), "\\code"),
                      x[[1L]])
      rows <<- c(rows, vapply(codes, function(y) paste(unlist(y),
                                                       collapse = ""),
                              character(1)))
    } else if (is.list(x)) {
      for (el in x) walk(el)
    }
  }
  walk(value)
  rows
}

# ---- the interface ---------------------------------------------------------

test_that("the sensitivity functions have the documented arguments", {
  shared <- c("true_value", "width", "conf_level", "distribution", "G",
              "seed", "m_0", "max_n", "batch")
  for (f in list(seq_aipe_smd_sensitivity, seq_aipe_r_sensitivity,
                 seq_aipe_cv_sensitivity, seq_aipe_reg_coef_sensitivity,
                 seq_aipe_sm_sensitivity)) {
    expect_identical(names(formals(f)), shared)
  }
  expect_identical(names(formals(seq_aipe_gini_sensitivity)),
                   setdiff(shared, "distribution"))
  expect_identical(
    names(formals(seq_aipe_c_sensitivity)),
    c("true_means", "contrast_weights", "width", "conf_level", "true_sds",
      "distribution", "G", "seed", "m_0", "max_n", "batch"))
  # Defaults: 95% confidence, purely sequential sampling, and no seed.
  for (f in list(seq_aipe_smd_sensitivity, seq_aipe_r_sensitivity,
                 seq_aipe_cv_sensitivity, seq_aipe_reg_coef_sensitivity,
                 seq_aipe_sm_sensitivity, seq_aipe_gini_sensitivity,
                 seq_aipe_c_sensitivity)) {
    fm <- formals(f)
    expect_identical(fm$conf_level, 0.95)
    expect_identical(fm$batch, 1L)
    expect_null(fm$seed)
    expect_identical(fm$G, 1000)
    expect_identical(fm$m_0, 4)
    expect_identical(fm$max_n, 100000L)
  }
  expect_identical(eval(formals(seq_aipe_cv_sensitivity)$distribution),
                   c("gamma", "lognormal"))
  expect_identical(eval(formals(seq_aipe_c_sensitivity)$distribution),
                   c("normal", "skewed"))
  expect_identical(formals(seq_aipe_c_sensitivity)$true_sds, 1)
})

# ---- what is returned ------------------------------------------------------

test_that("every sensitivity function returns the 13 documented rows", {
  outs <- list(
    smd = seq_aipe_smd_sensitivity(0.5, 1, G = 4, seed = 113),
    r = seq_aipe_r_sensitivity(0.4, 0.6, G = 4, seed = 113),
    cv = seq_aipe_cv_sensitivity(0.5, 0.5, G = 4, seed = 113),
    reg_coef = seq_aipe_reg_coef_sensitivity(0.5, 0.8, G = 4, seed = 113),
    sm = seq_aipe_sm_sensitivity(0.4, 1, G = 4, seed = 113),
    gini = seq_aipe_gini_sensitivity(0.3, 0.4, G = 4, seed = 113),
    c = seq_aipe_c_sensitivity(c(1, 0, 0), c(1, -0.5, -0.5), width = 1.5,
                               G = 4, seed = 113)
  )
  for (nm in names(outs)) {
    out <- outs[[nm]]
    # An smsd_sensitivity object: a term/value data frame with one more
    # class, which gives it its print method.
    expect_identical(class(out), c("smsd_sensitivity", "data.frame"),
                     info = nm)
    expect_identical(names(out), c("term", "value"), info = nm)
    expect_identical(dim(out), c(13L, 2L), info = nm)
    expect_identical(out$term, sens_terms, info = nm)
    expect_type(out$value, "double")
    expect_true(all(is.finite(out$value)), info = nm)
    # The reported number of replications is the number asked for.
    expect_identical(sens_value(out, "replications"), 4, info = nm)
    # Units: three percentages on 0 to 100, coverage a proportion.
    pct <- sens_value(out, c("pct_width_within", "pct_hit_max_n",
                             "pct_variance_unusable"))
    expect_true(all(pct >= 0 & pct <= 100), info = nm)
    cov <- sens_value(out, "coverage")
    expect_true(cov >= 0 && cov <= 1, info = nm)
  }
  # The echoed inputs.
  expect_identical(sens_value(outs$r, c("true_value", "target_width",
                                        "conf_level")), c(0.4, 0.6, 0.95))
  # For the contrast the true value is the contrast of the true means.
  expect_identical(sens_value(outs$c, "true_value"), 1)
})

test_that("a summary prints with the smsd_sensitivity method", {
  # A contrast of 1e20 (in units in which the standard deviation is 1e15):
  # the default data frame print would show the whole value column in
  # scientific notation (1e+20, 9.5e-01, ...). The class's method formats
  # each value on its own, in fixed notation.
  out <- seq_aipe_c_sensitivity(c(1e20, 0), c(1, -1), width = 1e16,
                                true_sds = 1e15, G = 2, seed = 113,
                                max_n = 8)
  expect_identical(sens_value(out, "true_value"), 1e20)
  shown <- capture.output(returned <- withVisible(print(out)))
  expect_identical(shown, format(out))
  expect_false(returned$visible)
  expect_identical(returned$value, out)
  expect_false(any(grepl("e[+-][0-9]", shown)))
  expect_match(shown[2], "^true_value +100000000000000000000$")
  expect_match(shown[4], "^conf_level +0\\.95$")
  expect_match(shown[5], "^replications +2$")
  # It is still an ordinary data frame to compute with.
  expect_identical(out[out$term == "replications", "value"], 2)
  expect_identical(as.data.frame(out)$term, sens_terms)
})

test_that("the help page lists every returned row, in order", {
  documented <- sens_documented_rows()
  skip_if(is.null(documented), "the help page source is not available")
  expect_identical(documented, sens_terms)
})

# ---- the whole procedure against a hand-written loop -----------------------
#
# Each sensitivity function is compared with hand_study() over the same random
# stream, and the two must agree exactly, row for row. The hand-written study
# draws the data itself, from the mechanism the help page states, and passes
# the confidence level, the minimum pilot size, the cap, and the batch size
# itself. A function that dropped one of those arguments, drew from the wrong
# distribution, ran another estimator, or judged coverage against another
# true value would therefore differ. Every case but the first uses a
# confidence level, a batch size, and a cap other than the defaults.

test_that("with every default, the summary equals a hand-written loop", {
  # The purely sequential procedure of the articles: one observation at a
  # time, 95% confidence, no cap in reach. Standardized mean 0.4, width 0.8.
  G <- 25
  theta <- 0.4
  width <- 0.8
  ref <- hand_study(
    G, 113, seq_aipe_sm_pilot(width = width),
    gen = function(k) list(stats::rnorm(k) + theta),
    check = function(d) seq_aipe_sm(d[[1]], width = width),
    truth = theta, width = width
  )
  out <- seq_aipe_sm_sensitivity(theta, width, G = G, seed = 113)
  expect_same_study(out, ref, "standardized mean, defaults", cap_binds = FALSE)
  expect_gt(stats::sd(ref$n), 0)
})

test_that("normal and skewed designs equal a hand-written loop", {
  G <- 20
  for (dist in c("normal", "skewed")) {
    rand <- hand_rand[[dist]]
    skew <- dist == "skewed"
    level <- if (skew) 0.99 else 0.90

    # Standardized mean difference 0.5: group 1 is shifted by delta.
    delta <- 0.5
    width <- if (skew) 1.2 else 0.9
    cap <- if (skew) 34L else 30L
    ref <- hand_study(
      G, 113, seq_aipe_smd_pilot(width = width, conf_level = level),
      gen = function(k) list(rand(k) + delta, rand(k)),
      check = function(d) {
        seq_aipe_smd(d[[1]], d[[2]], width = width, conf_level = level)
      },
      truth = delta, width = width, max_n = cap, batch = 3L
    )
    out <- seq_aipe_smd_sensitivity(delta, width, conf_level = level,
                                    distribution = dist, G = G, seed = 113,
                                    max_n = cap, batch = 3)
    expect_same_study(out, ref, paste("SMD,", dist),
                      default_pilot = seq_aipe_smd_pilot(width = width),
                      batch = 3)

    # Correlation -0.3: y = rho x + sqrt(1 - rho^2) e, Pearson. With the
    # normal data the minimum pilot size, 8, sets the pilot; with the skewed
    # data the confidence level does (9 at 0.99, 8 at 0.95).
    rho <- -0.3
    width <- if (skew) 0.6 else 0.5
    cap <- if (skew) 50L else 40L
    ref <- hand_study(
      G, 113, seq_aipe_r_pilot(width = width, conf_level = level, m_0 = 8),
      gen = function(k) {
        x <- rand(k)
        e <- rand(k)
        list(x, rho * x + sqrt(1 - rho^2) * e)
      },
      check = function(d) {
        seq_aipe_r(d[[1]], d[[2]], width = width, conf_level = level,
                   method = "pearson")
      },
      truth = rho, width = width, max_n = cap, batch = 2L
    )
    out <- seq_aipe_r_sensitivity(rho, width, conf_level = level,
                                  distribution = dist, G = G, seed = 113,
                                  m_0 = 8, max_n = cap, batch = 2)
    expect_same_study(
      out, ref, paste("correlation,", dist),
      default_pilot = if (skew) seq_aipe_r_pilot(width = width, m_0 = 8),
      batch = 2)
    if (!skew) {
      expect_gt(ref$pilot,
                seq_aipe_r_pilot(width = width, conf_level = level))
    }

    # Regression slope 0.5: y = beta x + e.
    beta <- 0.5
    width <- 0.7
    cap <- if (skew) 60L else 30L
    ref <- hand_study(
      G, 113, seq_aipe_reg_coef_pilot(width = width, conf_level = level),
      gen = function(k) {
        x <- rand(k)
        e <- rand(k)
        list(x, beta * x + e)
      },
      check = function(d) {
        seq_aipe_reg_coef(d[[1]], d[[2]], width = width, conf_level = level)
      },
      truth = beta, width = width, max_n = cap, batch = 4L
    )
    out <- seq_aipe_reg_coef_sensitivity(beta, width, conf_level = level,
                                         distribution = dist, G = G,
                                         seed = 113, max_n = cap, batch = 4)
    expect_same_study(out, ref, paste("slope,", dist),
                      default_pilot = seq_aipe_reg_coef_pilot(width = width),
                      batch = 4)

    # Standardized mean -0.6: one sample shifted by theta.
    theta <- -0.6
    width <- if (skew) 1.2 else 0.9
    cap <- if (skew) 45L else 16L
    ref <- hand_study(
      G, 113, seq_aipe_sm_pilot(width = width, conf_level = level),
      gen = function(k) list(rand(k) + theta),
      check = function(d) {
        seq_aipe_sm(d[[1]], width = width, conf_level = level)
      },
      truth = theta, width = width, max_n = cap, batch = 5L
    )
    out <- seq_aipe_sm_sensitivity(theta, width, conf_level = level,
                                   distribution = dist, G = G, seed = 113,
                                   max_n = cap, batch = 5)
    expect_same_study(out, ref, paste("standardized mean,", dist),
                      default_pilot = seq_aipe_sm_pilot(width = width),
                      batch = 5)

    # Contrast of four groups with unequal standard deviations, one group
    # with zero weight. True value 3.75 - (3.5 + 3.5) / 2 = 0.25.
    mu <- c(3.75, 3.5, 3.5, 3.1)
    sds <- c(0.5, 0.25, 1, 0.4)
    weights <- c(1, -0.5, -0.5, 0)
    width <- 0.7
    cap <- if (skew) 28L else 12L
    ref <- hand_study(
      G, 113, seq_aipe_c_pilot(width = width, conf_level = level),
      gen = function(k) lapply(1:4, function(j) mu[j] + sds[j] * rand(k)),
      check = function(d) {
        seq_aipe_c(d, contrast_weights = weights, width = width,
                   conf_level = level)
      },
      truth = 0.25, width = width, max_n = cap, batch = 3L
    )
    out <- seq_aipe_c_sensitivity(mu, weights, width, conf_level = level,
                                  true_sds = sds, distribution = dist, G = G,
                                  seed = 113, max_n = cap, batch = 3)
    expect_same_study(out, ref, paste("contrast,", dist),
                      default_pilot = seq_aipe_c_pilot(width = width),
                      batch = 3)
    expect_identical(sens_value(out, "true_value"), 0.25)
  }
})

test_that("the CV and Gini index functions equal a hand-written loop", {
  G <- 20
  # Coefficient of variation 0.7: gamma with shape 1 / 0.49, or lognormal
  # with log standard deviation sqrt(log(1.49)). (log1p(u) is log(1 + u)
  # without the rounding of 1 + u; the two can differ in the last digit.)
  kappa <- 0.7
  width <- 0.5
  gens <- list(
    gamma = function(k) list(stats::rgamma(k, shape = 1 / kappa^2)),
    lognormal = function(k) {
      list(stats::rlnorm(k, meanlog = 0, sdlog = sqrt(log1p(kappa^2))))
    }
  )
  for (dist in c("gamma", "lognormal")) {
    level <- if (dist == "gamma") 0.90 else 0.99
    cap <- if (dist == "gamma") 12L else 60L
    batch <- if (dist == "gamma") 3L else 2L
    ref <- hand_study(
      G, 113, seq_aipe_cv_pilot(width = width, conf_level = level, m_0 = 6),
      gen = gens[[dist]],
      check = function(d) {
        seq_aipe_cv(d[[1]], width = width, conf_level = level)
      },
      truth = kappa, width = width, max_n = cap, batch = batch
    )
    out <- seq_aipe_cv_sensitivity(kappa, width, conf_level = level,
                                   distribution = dist, G = G, seed = 113,
                                   m_0 = 6, max_n = cap, batch = batch)
    expect_same_study(out, ref, paste("coefficient of variation,", dist),
                      default_pilot = seq_aipe_cv_pilot(width = width,
                                                        m_0 = 6),
                      batch = batch)
  }

  # Gini index 0.45: lognormal with log standard deviation
  # sqrt(2) qnorm((0.45 + 1) / 2).
  gini <- 0.45
  width <- 0.3
  ref <- hand_study(
    G, 113, seq_aipe_gini_pilot(width = width, conf_level = 0.90, m_0 = 6),
    gen = function(k) {
      list(stats::rlnorm(k, meanlog = 0,
                         sdlog = sqrt(2) * stats::qnorm((gini + 1) / 2)))
    },
    check = function(d) seq_aipe_gini(d[[1]], width = width, conf_level = 0.90),
    truth = gini, width = width, max_n = 20L, batch = 2L
  )
  out <- seq_aipe_gini_sensitivity(gini, width, conf_level = 0.90, G = G,
                                   seed = 113, m_0 = 6, max_n = 20, batch = 2)
  expect_same_study(out, ref, "Gini index",
                    default_pilot = seq_aipe_gini_pilot(width = width,
                                                        m_0 = 6),
                    batch = 2)
})

# ---- the engine, driven by scripted checks --------------------------------

test_that("the engine stops at the first check that is met", {
  met_at_17 <- function(n) stub_check(n, satisfied = n >= 17)
  expect_identical(sens_value(stub_run(met_at_17), "mean_n"), 17)
  expect_identical(sens_value(stub_run(met_at_17, batch = 5), "mean_n"), 20)
  expect_identical(sens_value(stub_run(met_at_17, batch = 4), "mean_n"), 18)
  # Met at the pilot: no further sampling.
  expect_identical(sens_value(stub_run(met_at_17, pilot_n = 30L), "mean_n"),
                   30)
})

test_that("max_n is an exact cap for any batch size", {
  never <- function(n) stub_check(n, satisfied = FALSE)
  for (b in c(1L, 4L, 7L, 100L)) {
    out <- stub_run(never, pilot_n = 10L, max_n = 33L, batch = b)
    expect_identical(sens_value(out, "mean_n"), 33, info = paste("batch", b))
    expect_identical(sens_value(out, "pct_hit_max_n"), 100)
  }
  # The cap may equal the pilot: one check, no further sampling.
  out <- stub_run(never, pilot_n = 10L, max_n = 10L)
  expect_identical(sens_value(out, "mean_n"), 10)
  expect_identical(sens_value(out, "pct_hit_max_n"), 100)
})

test_that("pct_hit_max_n counts stops by the cap, not stops by the rule", {
  # The rule is first met exactly at the cap: stopped by the rule.
  met_at_30 <- function(n) stub_check(n, satisfied = n >= 30)
  out <- stub_run(met_at_30, pilot_n = 10L, max_n = 30L, batch = 100L)
  expect_identical(sens_value(out, "mean_n"), 30)
  expect_identical(sens_value(out, "pct_hit_max_n"), 0)
  # One more observation would have been needed: stopped by the cap.
  met_at_31 <- function(n) stub_check(n, satisfied = n >= 31)
  out <- stub_run(met_at_31, pilot_n = 10L, max_n = 30L, batch = 100L)
  expect_identical(sens_value(out, "mean_n"), 30)
  expect_identical(sens_value(out, "pct_hit_max_n"), 100)
})

test_that("max_n below the pilot sample size is an error that says so", {
  never <- function(n) stub_check(n, satisfied = FALSE)
  expect_error(stub_run(never, pilot_n = 10L, max_n = 9L),
               "'max_n' \\(9\\) is below the pilot sample size \\(10\\)")
  # Through an exported function: width 0.2 at 95% implies a pilot of 20.
  expect_identical(seq_aipe_smd_pilot(width = 0.2), 20L)
  expect_error(
    seq_aipe_smd_sensitivity(0.5, 0.2, G = 3, seed = 113, max_n = 10),
    "'max_n' \\(10\\) is below the pilot sample size \\(20\\)"
  )
  expect_error(
    seq_aipe_c_sensitivity(c(1, 0), c(1, -1), width = 0.2, G = 3,
                           max_n = 19),
    "below the pilot sample size \\(20\\)"
  )
  # The cap may equal the pilot.
  out <- seq_aipe_smd_sensitivity(0.5, 0.2, G = 3, seed = 113, max_n = 20)
  expect_identical(sens_value(out, "mean_n"), 20)
  # A pilot too large to represent is a clear error, not a missing value.
  expect_error(stub_run(never, pilot_n = NA_integer_),
               "pilot sample size.*too large to represent")
})

test_that("the engine runs exactly G replications", {
  pilots <- 0L
  counting_draw <- function(data, k) {
    if (is.null(data)) pilots <<- pilots + 1L
    stub_draw(data, k)
  }
  out <- stub_run(function(n) stub_check(n, satisfied = n >= 12), G = 7L,
                  draw = counting_draw)
  expect_identical(pilots, 7L)
  expect_identical(sens_value(out, "replications"), 7)
  # A single replication is allowed; the standard deviation is then NA.
  one <- stub_run(function(n) stub_check(n, satisfied = TRUE))
  expect_identical(sens_value(one, "replications"), 1)
  expect_identical(sens_value(one, "sd_n"), NA_real_)
})

test_that("pct_variance_unusable counts replications with such a check", {
  # The variance estimate is unusable at the pilot check of every other
  # replication; such a check is reported as not met, so sampling continues.
  replication <- 0L
  draw <- function(data, k) {
    if (is.null(data)) replication <<- replication + 1L
    stub_draw(data, k)
  }
  evaluate <- function(n) {
    unusable <- replication %% 2L == 1L && n == 10L
    stub_check(n, satisfied = !unusable && n >= 10L, unusable = unusable)
  }
  out <- stub_run(evaluate, G = 8L, draw = draw)
  expect_identical(sens_value(out, "pct_variance_unusable"), 50)
  # Those replications went on to 11; the others stopped at the pilot.
  expect_identical(sens_value(out, "mean_n"), 10.5)
  # A replication counts once however many of its checks were unusable: here
  # the first three checks of every replication are.
  three <- function(n) {
    stub_check(n, satisfied = n >= 13L, unusable = n < 13L)
  }
  out <- stub_run(three, G = 5L)
  expect_identical(sens_value(out, "pct_variance_unusable"), 100)
  expect_identical(sens_value(out, "mean_n"), 13)
  # The check at which the cap ends a replication counts as well.
  last_only <- function(n) {
    stub_check(n, satisfied = FALSE, unusable = n == 12L, lower = NA_real_,
               upper = NA_real_)
  }
  expect_identical(
    sens_value(stub_run(last_only, max_n = 12L), "pct_variance_unusable"), 100)
  expect_identical(
    sens_value(stub_run(last_only, max_n = 11L), "pct_variance_unusable"), 0)
  # A check function that does not report the attribute counts as usable.
  plain <- function(n) {
    res <- stub_check(n, satisfied = TRUE)
    attr(res, "variance_unusable") <- NULL
    res
  }
  expect_identical(sens_value(stub_run(plain), "pct_variance_unusable"), 0)
})

# ---- a replication that the cap ends without an interval -------------------
#
# A check at which no variance has been estimated (the raw estimate is zero,
# negative, or not finite) has no confidence limits and is never a stop, so a
# replication can end on one only when max_n ends it. The summary then counts
# the replication as not covering and as not within the target width, leaves
# it out of the two width rows, and counts it in every other row; no row
# becomes NA. With max_n equal to the pilot sample size each replication is
# exactly one check of its pilot sample, which makes every row computable by
# hand from the pilot samples, over the same random stream. A minimum pilot
# size of 10 at a width of 0.6 gives every kind of replication in numbers: of
# 400, about 60 to 70 without an interval, about 30 that meet the rule at
# once, and about 300 that the cap ends with an interval.

test_that("no interval at the cap: the slope summary, row by row", {
  # By hand, from the formulas of Kelley, Darku, and Chattopadhyay (2018):
  # the least squares slope b; V^2, the variance estimate of their Equation
  # 85 before its bound (from the package's estimator of it); and, when V^2
  # is usable, xi2 = max(V^2, n^-3) (Equation 84), the interval
  # b -/+ z sqrt(xi2 / n) (Equation 23), and the rule
  # n >= (2 z / width)^2 (xi2 + 1 / n) (Equation 22).
  G <- 400
  width <- 0.6
  beta <- 0.5
  n <- seq_aipe_reg_coef_pilot(width = width, m_0 = 10)
  expect_identical(n, 10L)
  z <- stats::qnorm(0.975)
  hand <- withr::with_seed(113, vapply(seq_len(G), function(i) {
    x <- stats::rnorm(n)
    e <- stats::rnorm(n)
    y <- beta * x + e
    b <- stats::cov(x, y) / stats::var(x)
    raw <- .var_reg_coef(x, y, floor = FALSE)
    usable <- !.var_unusable(raw)
    xi2 <- max(raw, n^-3)
    lower <- if (usable) b - z * sqrt(xi2 / n) else NA_real_
    upper <- if (usable) b + z * sqrt(xi2 / n) else NA_real_
    c(usable = usable,
      met = usable && n >= (2 * z / width)^2 * (xi2 + 1 / n),
      width = upper - lower,
      covered = usable && lower <= beta && beta <= upper)
  }, numeric(4)))
  usable <- hand["usable", ] == 1
  met <- hand["met", ] == 1
  # (No pilot in this stream has a positive estimate at or below the bound
  # 10^-3; that case has its own test below.)
  # The stream must hold all three kinds of replication, or the rows below
  # would not tell the readings apart: pilots with an unusable estimate,
  # pilots that meet the rule at once, and pilots that the cap ends with an
  # interval, many of which cover. (The counts themselves depend on the
  # estimator's arithmetic at borderline pilots and are not pinned.)
  expect_gt(sum(!usable), 30)
  expect_gt(sum(met), 10)
  expect_gt(sum(usable & !met), 100)
  expect_gt(sum(hand["covered", ]), 100)

  out <- seq_aipe_reg_coef_sensitivity(beta, width, G = G, seed = 113,
                                       m_0 = 10, max_n = n)
  expect_false(anyNA(out$value))
  expect_identical(sens_value(out, c("mean_n", "median_n", "sd_n")),
                   c(10, 10, 0))
  # Counted like any other replication.
  expect_identical(sens_value(out, "pct_variance_unusable"),
                   100 * sum(!usable) / G)
  expect_identical(sens_value(out, "pct_hit_max_n"), 100 * (G - sum(met)) / G)
  # Not covering and not within the width, as proportions of all G.
  expect_identical(sens_value(out, "coverage"), sum(hand["covered", ]) / G)
  expect_identical(
    sens_value(out, "pct_width_within"),
    100 * sum(usable & hand["width", ] <= width, na.rm = TRUE) / G)
  # Left out of the width rows. (The check function computes the slope on
  # rescaled data, so a width can differ from the hand value in its last
  # digits.)
  expect_equal(sens_value(out, "mean_width"), mean(hand["width", usable]),
               tolerance = 1e-12)
  expect_equal(sens_value(out, "median_width"),
               stats::median(hand["width", usable]), tolerance = 1e-12)
  # The other readings would give other numbers here, because some
  # replications have no interval and some of the others cover: coverage
  # among the replications that have an interval would be larger (its
  # divisor is sum(usable), not G), and a mean over every width would be NA.
})

test_that("no interval at the cap: the correlation, against the checks", {
  # The same design for the Pearson correlation, with the limits read from
  # the check function itself, one call per pilot sample.
  G <- 400
  width <- 0.6
  rho <- 0.4
  n <- seq_aipe_r_pilot(width = width, m_0 = 10)
  expect_identical(n, 10L)
  ref <- hand_study(
    G, 113, n,
    gen = function(k) {
      x <- stats::rnorm(k)
      e <- stats::rnorm(k)
      list(x, rho * x + sqrt(1 - rho^2) * e)
    },
    check = function(d) {
      seq_aipe_r(d[[1]], d[[2]], width = width, method = "pearson")
    },
    truth = rho, width = width, max_n = n
  )
  # The stream must hold replications of every kind: some that end without
  # an interval, some that meet the rule at the pilot, and many that the cap
  # ends with an interval. (The counts depend on the estimator's arithmetic
  # at borderline pilots and are not pinned.)
  no_interval <- is.na(ref$widths)
  capped <- round(G * ref$rows[["pct_hit_max_n"]] / 100)
  expect_gt(sum(no_interval), 30)
  expect_gt(G - capped, 10)
  expect_gt(capped - sum(no_interval), 100)
  # More intervals are within the target than replications met the rule
  # (the rule asks for a little more than the width alone, through its 1 / n
  # term), so pct_width_within and pct_hit_max_n are told apart.
  within <- sum(ref$widths <= width, na.rm = TRUE)
  expect_gt(within, G - capped)

  out <- seq_aipe_r_sensitivity(rho, width, G = G, seed = 113, m_0 = 10,
                                max_n = n)
  rows <- names(ref$rows)
  got <- sens_value(out, rows)
  names(got) <- rows
  expect_identical(got, ref$rows)
  expect_false(anyNA(out$value))
  # With normal data the correlation itself is always finite, so the
  # replications without an interval are those whose final (and only) check
  # reported an unusable variance estimate.
  expect_identical(sens_value(out, "pct_variance_unusable"),
                   100 * sum(no_interval) / G)
  expect_identical(sens_value(out, "pct_width_within"), 100 * within / G)
  expect_lte(sens_value(out, "coverage"), sum(!no_interval) / G)
  expect_identical(sens_value(out, "mean_width"),
                   mean(ref$widths[!no_interval]))
})

test_that("a positive estimate below the bound is used and is not counted", {
  # Correlation 0.4 at width 0.8: the pilot sample is 5 pairs and the bound
  # 5^-3 is 0.008, so a positive raw estimate below the bound is not rare.
  # With max_n equal to the pilot each replication is one check of its pilot
  # sample, and the rows follow by hand from Kelley, Darku, and
  # Chattopadhyay (2019): V^2, the raw variance estimate (from the package's
  # estimator of it, before the bound); and, when V^2 is positive,
  # xi2 = max(V^2, n^-3) (Equation 8, with gamma = 3), the rule
  # n >= (2 z / width)^2 (xi2 + 1 / n) (Equation 10), and the interval
  # r -/+ z sqrt(xi2 / n) (Equation 11). When V^2 is zero, negative, or not
  # finite there is no interval and the rule is not met; that, and only
  # that, is what pct_variance_unusable counts.
  G <- 2000
  width <- 0.8
  rho <- 0.4
  n <- seq_aipe_r_pilot(width = width)
  expect_identical(n, 5L)
  z <- stats::qnorm(0.975)
  hand <- withr::with_seed(113, vapply(seq_len(G), function(i) {
    x <- stats::rnorm(n)
    e <- stats::rnorm(n)
    y <- rho * x + sqrt(1 - rho^2) * e
    r <- stats::cor(x, y)
    raw <- .var_r(x, y, method = "pearson", floor = FALSE)
    unusable <- !is.finite(raw) || raw <= 0
    xi2 <- max(raw, n^-3)
    half <- if (unusable) NA_real_ else z * sqrt(xi2 / n)
    c(unusable = unusable,
      below = !unusable && raw <= n^-3,
      met = !unusable && n >= (2 * z / width)^2 * (xi2 + 1 / n),
      width = 2 * half,
      covered = !unusable && r - half <= rho && rho <= r + half)
  }, numeric(5)))
  unusable <- hand["unusable", ] == 1
  below <- hand["below", ] == 1
  # Not vacuous: the stream holds positive estimates below the bound, and
  # many estimates that are unusable.
  expect_gt(sum(below), 5)
  expect_gt(sum(unusable), 100)

  out <- seq_aipe_r_sensitivity(rho, width, G = G, seed = 113, max_n = n)
  # Only zero, negative, or non-finite estimates are counted.
  expect_identical(sens_value(out, "pct_variance_unusable"),
                   100 * sum(unusable) / G)
  # A positive estimate below the bound is used: the rule is evaluated with
  # the bound (and is met, since (2 z / 0.8)^2 (5^-3 + 1 / 5) = 4.99 is
  # below 5), so the cap did not end those replications, and their
  # intervals, of width 2 z sqrt(5^-3 / 5) = 0.157, are among those
  # summarized.
  expect_identical(sens_value(out, "pct_hit_max_n"),
                   100 * (G - sum(hand["met", ])) / G)
  expect_identical(sens_value(out, "pct_width_within"),
                   100 * sum(hand["width", ] <= width, na.rm = TRUE) / G)
  expect_equal(sens_value(out, "mean_width"),
               mean(hand["width", !unusable]), tolerance = 1e-10)
  expect_identical(sens_value(out, "coverage"), sum(hand["covered", ]) / G)
})

test_that("a pilot sample with a mean near zero does not abort the study", {
  # Gamma data with a coefficient of variation of 5 have shape 1 / 25, and
  # about 4 pilot samples of 4 in 1,000 then have a mean below the machine
  # epsilon (by simulation of the gamma distribution alone); with this seed
  # one occurs within the first 600. A check function that stopped with an
  # error there would abort the whole study. With max_n equal to the pilot,
  # each replication is one check of its pilot sample.
  expect_identical(seq_aipe_cv_pilot(width = 1), 4L)
  out <- seq_aipe_cv_sensitivity(5, 1, G = 600, seed = 113, max_n = 4)
  expect_identical(sens_value(out, "replications"), 600)
  expect_identical(sens_value(out, "mean_n"), 4)
  expect_true(all(is.finite(sens_value(out, c("coverage", "pct_width_within",
                                              "pct_hit_max_n")))))
})

test_that("a final interval that is not finite never counts as covering", {
  # Odd replications end at the cap with an estimate that is not finite; even
  # replications stop with the interval (-1, 1), which covers 0.
  replication <- 0L
  draw <- function(data, k) {
    if (is.null(data)) replication <<- replication + 1L
    stub_draw(data, k)
  }
  evaluate <- function(n) {
    if (replication %% 2L == 1L) {
      stub_check(n, satisfied = FALSE, lower = NaN, upper = NaN)
    } else {
      stub_check(n, satisfied = TRUE)
    }
  }
  out <- stub_run(evaluate, G = 6L, max_n = 10L, draw = draw)
  expect_identical(sens_value(out, "coverage"), 0.5)
  expect_identical(sens_value(out, "pct_width_within"), 50)
  expect_identical(sens_value(out, "pct_hit_max_n"), 50)
  expect_identical(sens_value(out, "mean_width"), 2)
  expect_identical(sens_value(out, "median_width"), 2)
})

test_that("when no replication ends with an interval, the widths are NA", {
  # Every replication is stopped by the cap at a check that gives no
  # interval (as the check functions do while the variance estimate is
  # unusable).
  no_interval <- function(n) {
    stub_check(n, satisfied = FALSE, unusable = TRUE, lower = NA_real_,
               upper = NA_real_)
  }
  out <- stub_run(no_interval, G = 4L, max_n = 12L)
  expect_identical(sens_value(out, "mean_width"), NA_real_)
  expect_identical(sens_value(out, "median_width"), NA_real_)
  expect_identical(sens_value(out, "pct_width_within"), 0)
  expect_identical(sens_value(out, "coverage"), 0)
  expect_identical(sens_value(out, "pct_hit_max_n"), 100)
  expect_identical(sens_value(out, "pct_variance_unusable"), 100)
  expect_identical(sens_value(out, "mean_n"), 12)
})

test_that("proportions and percentages are formed from counts", {
  # One replication in ten covers and is within the width, so coverage must
  # be exactly 0.1 and pct_width_within exactly 10; in a second study one
  # replication in ten has an unusable check and is stopped by the cap, so
  # those two rows must be exactly 10. Formed with mean() these values are
  # not exact where the long double type is no wider than a double (arm64,
  # for example, where mean(c(1, rep(0, 9))) is not 0.1): the test states the
  # intent on every platform and fails for mean() on those.
  replication <- 0L
  draw <- function(data, k) {
    if (is.null(data)) replication <<- replication + 1L
    stub_draw(data, k)
  }
  one_covers <- function(n) {
    if (replication == 1L) stub_check(n, satisfied = TRUE)
    else stub_check(n, satisfied = TRUE, lower = 1, upper = 4)
  }
  out <- stub_run(one_covers, G = 10L, draw = draw)
  expect_identical(sens_value(out, "coverage"), 0.1)
  expect_identical(sens_value(out, "pct_width_within"), 10)

  replication <- 0L
  one_capped <- function(n) {
    stub_check(n, satisfied = replication != 1L,
               unusable = replication == 1L)
  }
  out <- stub_run(one_capped, G = 10L, max_n = 12L, draw = draw)
  expect_identical(sens_value(out, "pct_hit_max_n"), 10)
  expect_identical(sens_value(out, "pct_variance_unusable"), 10)
  # Seven covering intervals in a thousand: exactly 0.007, and exactly 0.7
  # percent (100 * 7 / 1000; the product 100 * 0.007 is not 0.7 in double
  # precision on any platform).
  replication <- 0L
  seven_cover <- function(n) {
    if (replication <= 7L) stub_check(n, satisfied = TRUE)
    else stub_check(n, satisfied = TRUE, lower = 1, upper = 4)
  }
  out <- stub_run(seven_cover, G = 1000L, draw = draw)
  expect_identical(sens_value(out, "coverage"), 0.007)
  expect_identical(sens_value(out, "pct_width_within"), 0.7)
  # The same for the two rows about the cap and the unusable variance.
  replication <- 0L
  seven_capped <- function(n) {
    stub_check(n, satisfied = replication > 7L,
               unusable = replication <= 7L)
  }
  out <- stub_run(seven_capped, G = 1000L, max_n = 11L, draw = draw)
  expect_identical(sens_value(out, "pct_hit_max_n"), 0.7)
  expect_identical(sens_value(out, "pct_variance_unusable"), 0.7)
})

# ---- batch size ------------------------------------------------------------

test_that("with batch = 5 the final n is the pilot plus a multiple of 5", {
  # One replication per call, because only then is the final sample size
  # itself returned: twelve studies of one replication each (seeds 113 to
  # 124).
  width <- 0.8
  pilot <- seq_aipe_smd_pilot(width = width)
  final_n <- function(batch) {
    vapply(113:124, function(s) {
      out <- seq_aipe_smd_sensitivity(0.5, width, G = 1, seed = s,
                                      batch = batch)
      sens_value(out, "mean_n")
    }, numeric(1))
  }
  n <- final_n(5)
  expect_true(all((n - pilot) %% 5 == 0))
  expect_true(all(n > pilot))
  expect_gt(length(unique(n)), 1)
  # Not vacuous: checked after every observation, the same twelve studies do
  # not all end on that grid.
  expect_true(any((final_n(1) - pilot) %% 5 != 0))
  # Except at the cap, where the last batch is shortened: a width of 0.2
  # cannot be reached by 33 per group (pilot 20, then 25, 30, 33).
  capped <- seq_aipe_smd_sensitivity(0.5, 0.2, G = 3, seed = 113,
                                     max_n = 33, batch = 5)
  expect_identical(sens_value(capped, "mean_n"), 33)
  expect_identical(sens_value(capped, "sd_n"), 0)
  expect_identical(sens_value(capped, "pct_hit_max_n"), 100)
  expect_identical(sens_value(capped, "pct_width_within"), 0)
})

test_that("batch applies to every design", {
  # One replication per call, because only then is the final sample size
  # itself returned: eight studies of one replication each per function
  # (seeds 113 to 120), at targets the pilot sample does not usually meet.
  # With batch = 4 every final n must be the pilot plus a multiple of 4. The
  # same eight studies with batch = 1 do not all land on that grid, so a
  # function that ignored batch would fail. (The standardized mean
  # difference is in the test above.)
  designs <- list(
    r = list(
      run = function(...) seq_aipe_r_sensitivity(0.4, 0.4, ...),
      pilot = seq_aipe_r_pilot(width = 0.4)),
    cv = list(
      run = function(...) seq_aipe_cv_sensitivity(0.5, 0.3, ...),
      pilot = seq_aipe_cv_pilot(width = 0.3)),
    reg_coef = list(
      run = function(...) seq_aipe_reg_coef_sensitivity(0.5, 0.6, ...),
      pilot = seq_aipe_reg_coef_pilot(width = 0.6)),
    sm = list(
      run = function(...) seq_aipe_sm_sensitivity(0.4, 0.8, ...),
      pilot = seq_aipe_sm_pilot(width = 0.8)),
    gini = list(
      run = function(...) seq_aipe_gini_sensitivity(0.3, 0.3, ...),
      pilot = seq_aipe_gini_pilot(width = 0.3)),
    c = list(
      run = function(...) {
        seq_aipe_c_sensitivity(c(1, 0, 0), c(1, -0.5, -0.5), width = 1.2, ...)
      },
      pilot = seq_aipe_c_pilot(width = 1.2))
  )
  for (nm in names(designs)) {
    run <- designs[[nm]]$run
    pilot <- designs[[nm]]$pilot
    final_n <- function(batch) {
      vapply(113:120, function(s) {
        sens_value(run(G = 1, seed = s, batch = batch), "mean_n")
      }, numeric(1))
    }
    in_fours <- final_n(4)
    expect_true(all((in_fours - pilot) %% 4 == 0), info = nm)
    expect_true(any(in_fours > pilot), info = nm)
    expect_true(any((final_n(1) - pilot) %% 4 != 0), info = nm)
  }
})

# ---- the data generating mechanisms ----------------------------------------

test_that("the gamma and lognormal parameters give the requested CV", {
  # Population coefficient of variation by numerical integration of the
  # generating density (not the closed form the package uses).
  pop_cv <- function(dens) {
    m1 <- stats::integrate(function(x) x * dens(x), 0, Inf,
                           rel.tol = 1e-10)$value
    m2 <- stats::integrate(function(x) x^2 * dens(x), 0, Inf,
                           rel.tol = 1e-10)$value
    sqrt(m2 - m1^2) / m1
  }
  for (cv in c(0.2, 0.5, 1)) {
    shape <- .sensitivity_gamma_shape_cv(cv)
    expect_equal(pop_cv(function(x) stats::dgamma(x, shape = shape)), cv,
                 tolerance = 1e-6)
    sdlog <- .sensitivity_sdlog_cv(cv)
    expect_equal(pop_cv(function(x) stats::dlnorm(x, 0, sdlog)), cv,
                 tolerance = 1e-6)
  }
})

test_that("the lognormal parameter is accurate for a small CV", {
  # Lognormal data with log standard deviation s have a coefficient of
  # variation of sqrt(exp(s^2) - 1) (the closed form that the integrals above
  # confirm at ordinary values; at these values the density is too narrow to
  # integrate numerically). Computed as sqrt(log(1 + cv^2)), s loses the
  # digits of cv^2 that 1 + cv^2 cannot hold: the population value is then
  # off by 4e-5 of itself at 1e-6 and by a third at 1.1e-8.
  for (cv in c(1e-3, 1e-6, 1e-7, 3e-8, 1.1e-8)) {
    s <- .sensitivity_sdlog_cv(cv)
    expect_equal(sqrt(expm1(s^2)), cv, tolerance = 1e-12)
  }
})

test_that("the lognormal parameter gives the requested Gini index", {
  # Population Gini index from its definition, G = Delta / (2 mu) with
  # Delta = E|X1 - X2| = 2 * integral of F (1 - F), by numerical integration
  # (not the closed form the package uses).
  pop_gini <- function(sdlog) {
    half_delta <- stats::integrate(function(x) {
      p <- stats::plnorm(x, 0, sdlog)
      p * (1 - p)
    }, 0, Inf, rel.tol = 1e-10)$value
    half_delta / exp(sdlog^2 / 2)
  }
  for (g in c(0.1, 0.3, 0.6)) {
    expect_equal(pop_gini(.sensitivity_sdlog_gini(g)), g, tolerance = 1e-6)
  }
})

test_that("each function generates data with the stated value and shape", {
  # (with_mocked_bindings() needs testthat 3.2.0, which DESCRIPTION declares.)
  # The check function is replaced by one that records the pilot sample it is
  # handed and reports the rule as met, so each call below draws exactly one
  # sample of 100,000 (per group) from the function's own generator. The
  # sample statistic must then be within five standard errors of the
  # requested true value. The standard error (SE) beside each tolerance is
  # that of the statistic at n = 100,000, estimated by simulating the
  # mechanism as the help page describes it, independently of the package.
  #
  # The sample must also have the shape of the distribution asked for, judged
  # by its skewness: 0 for the normal mechanism (SE 0.008) and sqrt(8 / 3) =
  # 1.63 for the skewed one (SE 0.020); 2 * 0.5 = 1 for gamma data with a
  # coefficient of variation of 0.5 (SE 0.012) and (0.25 + 3) * 0.5 = 1.625
  # for lognormal data (SE 0.034). Each tolerance below is over five SEs, and
  # none reaches the other distribution.
  #
  # The recorded call also shows how the check function was called: the data
  # by position (the names of the data arguments differ across the check
  # functions and are not this family's to rely on) and every other argument
  # by name.
  n <- 100000
  seen <- NULL
  record <- function(...) {
    seen <<- list(...)
    stub_check(n, satisfied = TRUE)
  }
  expect_shape <- function(x, dist, what) {
    target <- c(normal = 0, skewed = sqrt(8 / 3), gamma = 1,
                lognormal = 1.625)[[dist]]
    tolerance <- c(normal = 0.05, skewed = 0.12, gamma = 0.07,
                   lognormal = 0.2)[[dist]]
    expect_lt(abs(sens_skewness(x) - target), tolerance,
              label = paste("skewness error,", what, dist))
  }

  for (dist in c("normal", "skewed")) {
    # Standardized mean difference 0.5: SE 0.0045 (normal), 0.0047 (skewed).
    with_mocked_bindings(
      seq_aipe_smd_sensitivity(0.5, 1, distribution = dist, G = 1,
                               seed = 113, m_0 = n, max_n = n),
      seq_aipe_smd = record
    )
    expect_named(seen, c("", "", "width", "conf_level"))
    g1 <- seen[[1]]
    g2 <- seen[[2]]
    expect_length(g1, n)
    pooled_sd <- sqrt((stats::var(g1) + stats::var(g2)) / 2)
    expect_lt(abs((mean(g1) - mean(g2)) / pooled_sd - 0.5), 0.025,
              label = paste("SMD error,", dist))
    expect_shape(g1, dist, "SMD group 1,")
    expect_shape(g2, dist, "SMD group 2,")

    # Correlation 0.4: SE 0.0026 (normal), 0.0031 (skewed).
    with_mocked_bindings(
      seq_aipe_r_sensitivity(0.4, 1, distribution = dist, G = 1, seed = 113,
                             m_0 = n, max_n = n),
      seq_aipe_r = record
    )
    expect_named(seen, c("", "", "width", "conf_level", "method"))
    expect_identical(seen$method, "pearson")
    x <- seen[[1]]
    y <- seen[[2]]
    expect_lt(abs(stats::cor(x, y) - 0.4), 0.016,
              label = paste("correlation error,", dist))
    # x, and the part of y that is independent of x.
    expect_shape(x, dist, "correlation x,")
    expect_shape(y - 0.4 * x, dist, "correlation e,")

    # Slope 0.5: SE 0.0030 (normal), 0.0032 (skewed).
    with_mocked_bindings(
      seq_aipe_reg_coef_sensitivity(0.5, 1, distribution = dist, G = 1,
                                    seed = 113, m_0 = n, max_n = n),
      seq_aipe_reg_coef = record
    )
    expect_named(seen, c("", "", "width", "conf_level"))
    x <- seen[[1]]
    y <- seen[[2]]
    expect_lt(abs(stats::cov(x, y) / stats::var(x) - 0.5), 0.016,
              label = paste("slope error,", dist))
    # The error term has unit variance: SE of its variance at most 0.008.
    expect_lt(abs(stats::var(y - 0.5 * x) - 1), 0.04,
              label = paste("slope error variance,", dist))
    expect_shape(x, dist, "slope x,")
    expect_shape(y - 0.5 * x, dist, "slope e,")

    # Standardized mean 0.4: SE 0.0033 (normal), 0.0022 (skewed).
    with_mocked_bindings(
      seq_aipe_sm_sensitivity(0.4, 1, distribution = dist, G = 1, seed = 113,
                              m_0 = n, max_n = n),
      seq_aipe_sm = record
    )
    expect_named(seen, c("", "width", "conf_level"))
    expect_lt(abs(mean(seen[[1]]) / stats::sd(seen[[1]]) - 0.4), 0.017,
              label = paste("standardized mean error,", dist))
    expect_shape(seen[[1]], dist, "standardized mean,")

    # Contrast (1, -1/2, -1/2) of means (75, 70, 70) with standard
    # deviations (10, 10, 20): true value 5; n times the variance of the
    # estimate is 100 + 25 + 100 = 225, so the SE is 0.047. The SE of a
    # group's standard deviation is at most 0.037 for sigma = 10 and 0.078
    # for sigma = 20 (the skewed data).
    with_mocked_bindings(
      seq_aipe_c_sensitivity(c(75, 70, 70), c(1, -0.5, -0.5), width = 1,
                             true_sds = c(10, 10, 20), distribution = dist,
                             G = 1, seed = 113, m_0 = n, max_n = n),
      seq_aipe_c = record
    )
    expect_named(seen, c("", "contrast_weights", "width", "conf_level"))
    groups <- seen[[1]]
    expect_length(groups, 3L)
    expect_identical(seen$contrast_weights, c(1, -0.5, -0.5))
    expect_identical(seen$width, 1)
    expect_identical(seen$conf_level, 0.95)
    means <- vapply(groups, mean, numeric(1))
    expect_lt(abs(sum(c(1, -0.5, -0.5) * means) - 5), 0.25,
              label = paste("contrast error,", dist))
    sds <- vapply(groups, stats::sd, numeric(1))
    expect_true(all(abs(sds - c(10, 10, 20)) < c(0.2, 0.2, 0.4)),
                label = paste("group standard deviations,", dist))
    for (j in 1:3) expect_shape(groups[[j]], dist, paste("contrast group", j))
  }

  # Coefficient of variation 0.5: SE 0.0012 (gamma), 0.0016 (lognormal).
  for (dist in c("gamma", "lognormal")) {
    with_mocked_bindings(
      seq_aipe_cv_sensitivity(0.5, 1, distribution = dist, G = 1, seed = 113,
                              m_0 = n, max_n = n),
      seq_aipe_cv = record
    )
    expect_named(seen, c("", "width", "conf_level"))
    x <- seen[[1]]
    expect_true(all(x > 0))
    expect_lt(abs(stats::sd(x) / mean(x) - 0.5), 0.009,
              label = paste("coefficient of variation error,", dist))
    expect_shape(x, dist, "coefficient of variation,")
  }

  # Gini index 0.3: SE 0.0007. (The width is below 1, the length of the
  # range of the index; the check is the recording stand-in, so the width
  # only has to be one the function accepts.)
  with_mocked_bindings(
    seq_aipe_gini_sensitivity(0.3, 0.9, G = 1, seed = 113, m_0 = n,
                              max_n = n),
    seq_aipe_gini = record
  )
  expect_named(seen, c("", "width", "conf_level"))
  x <- seen[[1]]
  expect_true(all(x > 0))
  expect_lt(abs(gini_index(x) - 0.3), 0.004)
})

test_that("the skewed mechanism is standardized", {
  # (chi square(3) - 3) / sqrt(6): mean 0, variance 1. With 200,000 draws the
  # SE of the mean is 0.0022 and of the variance sqrt(6 / n) = 0.0055 (the
  # kurtosis is 3 + 12 / 3 = 7), so the tolerances are over five SEs. The
  # third central moment is the skewness, sqrt(8 / 3) = 1.63: right skewed.
  z <- withr::with_seed(113, .skewed_std(200000))
  expect_lt(abs(mean(z)), 0.012)
  expect_lt(abs(stats::var(z) - 1), 0.03)
  expect_gt(mean((z - mean(z))^3), 1)
})

# ---- operating characteristics (seeded Monte Carlo) -------------------------
#
# Each test compares the mean final sample size with the optimal fixed sample
# size n_omega = (2 z / omega)^2 xi^2, computed from the population in closed
# form, and the coverage with the nominal 0.95. The sequential procedure does
# not hit n_omega exactly at a finite width; the size of that excess, and the
# coverage, were measured for each condition in a separate large run (5,000 or
# 10,000 replications, other seeds) and are quoted with their standard errors
# (SE). Every tolerance is that measured offset plus at least four Monte
# Carlo SEs of the quantity at the G used here.

test_that("SMD: the final n tracks the optimal n and the interval covers", {
  # Normal groups, delta = 0.5, width 0.8. Per group, xi^2 = 2 + delta^2 / 4.
  out <- seq_aipe_smd_sensitivity(0.5, 0.8, G = 120, seed = 113)
  n_omega <- (2 * stats::qnorm(0.975) / 0.8)^2 * (2 + 0.5^2 / 4)
  expect_equal(n_omega, 49.52, tolerance = 1e-3)
  # Large runs: mean n 50.58 (SE 0.03), 1.06 above n_omega, with an sd of
  # n of 3.35; and seq_aipe_smd_sensitivity(0.5, 0.8, G = 12000,
  # seed = 20261012) gives 50.55 (SE 0.03), 1.03 above, with an sd of 3.37.
  # The excess is about 1.05, and no more than 1.1. The SE of the mean at
  # G = 120 is 3.37 / sqrt(120) = 0.31, so the tolerance of 2.4 is the
  # excess plus at least 4.2 SE.
  expect_lt(abs(sens_value(out, "mean_n") - n_omega), 2.4)
  # Large run: coverage 0.9515 (SE 0.0021). At G = 120 the SE is 0.020, so
  # coverage below 0.9515 - 4 * 0.020 = 0.87 would be a failure.
  expect_gt(sens_value(out, "coverage"), 0.87)
  # The requested width is delivered in every replication.
  expect_identical(sens_value(out, "pct_width_within"), 100)
  expect_lte(sens_value(out, "mean_width"), 0.8)
  expect_identical(sens_value(out, "pct_hit_max_n"), 0)
})

test_that("standardized mean, skewed data: coverage and final n", {
  # theta = 0.4, width 0.5, chi square(3) data standardized to unit variance
  # (skewness g1 = sqrt(8 / 3), kurtosis k = 7). By the delta method,
  # xi^2 = 1 - theta g1 + theta^2 (k - 1) / 4.
  theta <- 0.4
  out <- seq_aipe_sm_sensitivity(theta, 0.5, distribution = "skewed",
                                 G = 300, seed = 113)
  xi2 <- 1 - theta * sqrt(8 / 3) + theta^2 * (7 - 1) / 4
  n_omega <- (2 * stats::qnorm(0.975) / 0.5)^2 * xi2
  expect_equal(n_omega, 36.07, tolerance = 1e-3)
  # Large runs: coverage 0.9597 (SE 0.0010) in
  # seq_aipe_sm_sensitivity(0.4, 0.5, distribution = "skewed", G = 40000,
  # seed = 20261012), and 0.9631 (SE 0.0019) in another run. At G = 300 the
  # SE is sqrt(0.9597 * 0.0403 / 300) = 0.0114, so coverage below
  # 0.9597 - 4 * 0.0114 = 0.914 would be a failure.
  expect_gt(sens_value(out, "coverage"), 0.914)
  # The same run of 40,000: mean n 40.40 (SE 0.04), 4.33 above n_omega at
  # this wide target; sd of n 7.61, so the SE of the mean at G = 300 is
  # 0.44. Tolerance 4.33 + 4 * 0.44 = 6.1.
  expect_lt(abs(sens_value(out, "mean_n") - n_omega), 6.1)
  expect_identical(sens_value(out, "pct_width_within"), 100)
})

test_that("contrast: the final n tracks the optimal n, the interval covers", {
  # Means (75, 70, 70), common standard deviation 10, weights (1, -1/2,
  # -1/2), width 6. Per group, xi^2 = sum c_k^2 sigma_k^2 = 150.
  w <- c(1, -0.5, -0.5)
  out <- seq_aipe_c_sensitivity(c(75, 70, 70), w, width = 6, true_sds = 10,
                                G = 150, seed = 113)
  n_omega <- (2 * stats::qnorm(0.975) / 6)^2 * sum(w^2 * 10^2)
  expect_equal(n_omega, 64.02, tolerance = 1e-3)
  # Large run: mean n 63.64 (SE 0.12), 0.38 below n_omega; sd of n 8.33, so
  # the SE of the mean at G = 150 is 0.68. Tolerance 0.38 + 4 * 0.68 = 3.1.
  expect_lt(abs(sens_value(out, "mean_n") - n_omega), 3.1)
  # Large run: coverage 0.9450 (SE 0.0032). At G = 150 the SE is 0.019, so
  # coverage below 0.9450 - 4 * 0.019 = 0.87 would be a failure.
  expect_gt(sens_value(out, "coverage"), 0.87)
  expect_identical(sens_value(out, "true_value"), 5)
  expect_identical(sens_value(out, "pct_width_within"), 100)
  expect_lte(sens_value(out, "mean_width"), 6)
})

# ---- seeds -----------------------------------------------------------------

test_that("a seeded call reproduces, and another seed gives another study", {
  a <- seq_aipe_sm_sensitivity(0.4, 0.8, G = 10, seed = 113)
  b <- seq_aipe_sm_sensitivity(0.4, 0.8, G = 10, seed = 113)
  # Any seed other than 113 will do for the study that must differ.
  other <- seq_aipe_sm_sensitivity(0.4, 0.8, G = 10, seed = 114)
  expect_identical(a, b)
  expect_false(identical(a$value, other$value))
  # A whole number seed may be stored as a double or an integer.
  expect_identical(seq_aipe_cv_sensitivity(0.5, 0.5, G = 5, seed = 113L),
                   seq_aipe_cv_sensitivity(0.5, 0.5, G = 5, seed = 113))
})

test_that("a seeded call leaves the caller's stream as it found it", {
  # The caller's stream is seeded with 113; the call inside it is given
  # another seed (2026) so that it cannot leave the caller's state unchanged
  # by coincidence.
  with_call <- withr::with_seed(113, {
    first <- stats::runif(1)
    seq_aipe_smd_sensitivity(0.5, 1, G = 3, seed = 2026)
    c(first, stats::runif(1))
  })
  without_call <- withr::with_seed(113, stats::runif(2))
  expect_identical(with_call, without_call)
})

test_that("seed = NULL draws from, and advances, the caller's stream", {
  from_stream <- withr::with_seed(
    113, seq_aipe_cv_sensitivity(0.5, 0.5, G = 5))
  expect_identical(from_stream,
                   seq_aipe_cv_sensitivity(0.5, 0.5, G = 5, seed = 113))
  two <- withr::with_seed(113, list(
    seq_aipe_cv_sensitivity(0.5, 0.5, G = 5),
    seq_aipe_cv_sensitivity(0.5, 0.5, G = 5)
  ))
  expect_false(identical(two[[1]]$value, two[[2]]$value))
})

test_that("a seeded call leaves a fresh session without a seed", {
  skip_on_cran()
  # A session that has not used the generator has no seed, and its first
  # draw is seeded from the clock. If a seeded call left its seed behind,
  # the next draw would instead be the same in every such session. Two fresh
  # R sessions each make the same seeded call and then draw one uniform
  # number: the two draws must differ. (No generator state is read.)
  ns_path <- getNamespaceInfo(asNamespace("SMSD"), "path")
  from_source <- file.exists(file.path(ns_path, "R", "sensitivity.R"))
  load_line <- if (from_source) {
    sprintf("suppressMessages(pkgload::load_all(%s, quiet = TRUE))",
            deparse(ns_path))
  } else {
    sprintf("suppressMessages(library(SMSD, lib.loc = %s))",
            deparse(dirname(ns_path)))
  }
  script <- withr::local_tempfile(fileext = ".R")
  writeLines(c(
    sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
    load_line,
    "out <- seq_aipe_smd_sensitivity(0.5, 1, G = 2, seed = 113)",
    "cat(format(stats::runif(1), digits = 17), '\\n')"
  ), script)
  fresh_draw <- function() {
    res <- suppressWarnings(system2(
      file.path(R.home("bin"), "Rscript"),
      c("--no-init-file", "--no-site-file", "--no-save", "--no-restore",
        shQuote(script)),
      stdout = TRUE, stderr = FALSE
    ))
    as.numeric(res[length(res)])
  }
  first <- fresh_draw()
  second <- fresh_draw()
  expect_length(first, 1L)
  expect_length(second, 1L)
  expect_true(is.finite(first) && is.finite(second))
  expect_false(identical(first, second))
})

# ---- argument checks -------------------------------------------------------

test_that("the shared arguments are checked, naming what was received", {
  f <- seq_aipe_smd_sensitivity
  # true_value
  expect_error(f(NA_real_, 0.8), "'true_value' must be a single finite number")
  expect_error(f(c(0.2, 0.5), 0.8), "'true_value'.*vector of length 2")
  expect_error(f(matrix(0.5), 0.8), "'true_value'.*dimensions 1 x 1")
  expect_error(f("a", 0.8), "'true_value' must be a single finite number")
  # width
  expect_error(f(0.5), "'width' must be supplied")
  expect_error(f(0.5, -1), "'width' must be a single finite number greater")
  expect_error(f(0.5, 0), "'width'.*got 0")
  expect_error(f(0.5, NA_real_), "'width'.*got NA")
  expect_error(f(0.5, Inf), "'width'.*got Inf")
  # conf_level
  expect_error(f(0.5, 0.8, conf_level = 1), "'conf_level'.*got 1")
  expect_error(f(0.5, 0.8, conf_level = 0), "'conf_level'.*got 0")
  expect_error(f(0.5, 0.8, conf_level = NA_real_), "'conf_level'.*got NA")
  expect_error(f(0.5, 0.8, conf_level = 95), "'conf_level'.*got 95")
  # G
  expect_error(f(0.5, 0.8, G = 2.5), "'G' must be a single whole number")
  expect_error(f(0.5, 0.8, G = 0), "'G'.*at least 1; got 0")
  expect_error(f(0.5, 0.8, G = NA_real_), "'G'.*got NA")
  expect_error(f(0.5, 0.8, G = "10"), "'G'.*got \"10\"")
  # m_0
  expect_error(f(0.5, 0.8, m_0 = 3), "'m_0'.*at least 4; got 3")
  expect_error(f(0.5, 0.8, m_0 = 4.5), "'m_0'.*got 4.5")
  expect_error(f(0.5, 0.8, m_0 = NA_real_), "'m_0'.*got NA")
  expect_error(f(0.5, 0.8, m_0 = Inf), "'m_0'.*got Inf")
  # max_n
  expect_error(f(0.5, 0.8, max_n = NA), "'max_n'.*got NA")
  expect_error(f(0.5, 0.8, max_n = "a"), "'max_n'.*got \"a\"")
  expect_error(f(0.5, 0.8, max_n = 30.5), "'max_n'.*got 30.5")
  expect_error(f(0.5, 0.8, max_n = 0), "'max_n'.*got 0")
  # batch
  expect_error(f(0.5, 0.8, batch = 0), "'batch'.*at least 1; got 0")
  expect_error(f(0.5, 0.8, batch = 2.5), "'batch'.*got 2.5")
  expect_error(f(0.5, 0.8, batch = NA), "'batch'.*got NA")
  expect_error(f(0.5, 0.8, batch = c(1, 2)), "'batch'.*vector of length 2")
  # seed
  expect_error(f(0.5, 0.8, seed = 1.7),
               "'seed' must be NULL or a single whole number.*got 1.7")
  expect_error(f(0.5, 0.8, seed = "a"), "'seed'.*got \"a\"")
  expect_error(f(0.5, 0.8, seed = NA_real_), "'seed'.*got NA")
  expect_error(f(0.5, 0.8, seed = c(1, 2)), "'seed'.*vector of length 2")
  expect_error(f(0.5, 0.8, seed = matrix(1)), "'seed'.*dimensions 1 x 1")
  expect_error(f(0.5, 0.8, seed = 1e10), "'seed'.*got 1e\\+10")
  # A negative whole number is a valid seed (hence -113, not 113), and it is
  # not read as its absolute value.
  expect_identical(f(0.5, 1, G = 2, seed = -113),
                   f(0.5, 1, G = 2, seed = -113))
  expect_false(identical(f(0.5, 1, G = 2, seed = -113)$value,
                         f(0.5, 1, G = 2, seed = 113)$value))
  # distribution
  expect_error(f(0.5, 0.8, distribution = "cauchy"),
               paste0("'distribution' must be one of \"normal\", ",
                      "\"skewed\"; got \"cauchy\""))
  expect_error(f(0.5, 0.8, distribution = 1), "'distribution'.*got 1")
  expect_error(f(0.5, 0.8, distribution = NA_character_), "'distribution'")
  expect_error(f(0.5, 0.8, distribution = c("skewed", "normal")),
               "'distribution'.*vector of length 2")
  expect_error(seq_aipe_cv_sensitivity(0.5, 0.5, distribution = "normal"),
               "'distribution' must be one of \"gamma\", \"lognormal\"")
})

test_that("every function checks the shared arguments", {
  calls <- list(
    smd = function(...) seq_aipe_smd_sensitivity(0.5, 0.8, ...),
    r = function(...) seq_aipe_r_sensitivity(0.4, 0.6, ...),
    cv = function(...) seq_aipe_cv_sensitivity(0.5, 0.5, ...),
    reg_coef = function(...) seq_aipe_reg_coef_sensitivity(0.5, 0.8, ...),
    sm = function(...) seq_aipe_sm_sensitivity(0.4, 1, ...),
    gini = function(...) seq_aipe_gini_sensitivity(0.3, 0.4, ...),
    c = function(...) seq_aipe_c_sensitivity(c(1, 0), c(1, -1), 1.5, ...)
  )
  for (nm in names(calls)) {
    f <- calls[[nm]]
    expect_error(f(G = 2.5), "'G'", info = nm)
    expect_error(f(m_0 = 3), "'m_0'", info = nm)
    expect_error(f(max_n = "a"), "'max_n'", info = nm)
    expect_error(f(batch = 0), "'batch'", info = nm)
    expect_error(f(seed = 1.7), "'seed'", info = nm)
    expect_error(f(conf_level = 1), "'conf_level'", info = nm)
    if (nm != "gini") {
      expect_error(f(distribution = "cauchy"),
                   "'distribution' must be one of .*; got \"cauchy\"",
                   info = nm)
    }
  }
  # The Gini index function has no distribution argument.
  expect_error(seq_aipe_gini_sensitivity(0.3, 0.4, distribution = "gamma"),
               "unused argument")
  expect_error(seq_aipe_r_sensitivity(0.4, -1), "'width'")
  expect_error(seq_aipe_cv_sensitivity(0.5, NA_real_), "'width'")
  expect_error(seq_aipe_reg_coef_sensitivity(0.5, 0), "'width'")
  expect_error(seq_aipe_sm_sensitivity(0.4), "'width' must be supplied")
  expect_error(seq_aipe_gini_sensitivity(0.3, Inf), "'width'")
  expect_error(seq_aipe_c_sensitivity(c(1, 0), c(1, -1)),
               "'width' must be supplied")
})

test_that("arguments are checked by the functions and again by the engine", {
  # The shared check returns the cleaned values: doubles for the three
  # numbers, integers for the four counts, an integer or NULL for the seed.
  a <- .sensitivity_check_args(0.5, 0.8, 0.9, 10, 6, 50, 2, 113)
  expect_identical(a, list(true_value = 0.5, width = 0.8, conf_level = 0.9,
                           G = 10L, m_0 = 6L, max_n = 50L, batch = 2L,
                           seed = 113L))
  expect_null(.sensitivity_check_args(0.5, 0.8, 0.9, 10, 6)$seed)
  expect_error(.sensitivity_check_args(NA_real_, 0.8, 0.95, 10, 4),
               "'true_value'")
  expect_error(.sensitivity_check_args(0.5, 0, 0.95, 10, 4), "'width'")
  expect_error(.sensitivity_check_args(0.5, 0.8, 1, 10, 4), "'conf_level'")
  expect_error(.sensitivity_check_args(0.5, 0.8, 0.95, 2.5, 4), "'G'")
  expect_error(.sensitivity_check_args(0.5, 0.8, 0.95, 10, 3), "'m_0'")
  expect_error(.sensitivity_check_args(0.5, 0.8, 0.95, 10, 4, max_n = 0),
               "'max_n'")
  expect_error(.sensitivity_check_args(0.5, 0.8, 0.95, 10, 4, batch = 0),
               "'batch'")
  expect_error(.sensitivity_check_args(0.5, 0.8, 0.95, 10, 4, seed = 1.5),
               "'seed'")
  # The engine checks its own arguments as well, for a caller that does not
  # go through the shared check.
  met <- function(n) stub_check(n, satisfied = TRUE)
  expect_error(stub_run(met, G = 0L), "'G' must be a single whole number")
  expect_error(stub_run(met, max_n = 2.5), "'max_n'.*got 2.5")
  expect_error(stub_run(met, batch = 0L), "'batch'.*at least 1; got 0")
  # The seed is refused before anything is drawn.
  draws <- 0L
  counting_draw <- function(data, k) {
    draws <<- draws + 1L
    stub_draw(data, k)
  }
  expect_error(
    .sensitivity_run(G = 1, seed = "a", pilot_n = 10L, draw = counting_draw,
                     evaluate = met, true_value = 0, width = 2,
                     conf_level = 0.95, max_n = 20),
    "'seed' must be NULL or a single whole number")
  expect_error(
    .sensitivity_run(G = 1, seed = 1e10, pilot_n = 10L, draw = counting_draw,
                     evaluate = met, true_value = 0, width = 2,
                     conf_level = 0.95, max_n = 20),
    "'seed' must be NULL or a single whole number")
  expect_identical(draws, 0L)
})

test_that("the true value is checked against the range of the effect size", {
  expect_error(seq_aipe_r_sensitivity(1.2, 0.3),
               "'true_value'.*greater than -1 and less than 1; got 1.2")
  expect_error(seq_aipe_r_sensitivity(1, 0.3), "'true_value'.*got 1")
  expect_error(seq_aipe_r_sensitivity(-1, 0.3), "'true_value'.*got -1")
  expect_error(seq_aipe_cv_sensitivity(-0.5, 0.3),
               "'true_value'.*greater than 0; got -0.5")
  expect_error(seq_aipe_cv_sensitivity(0, 0.3), "'true_value'.*got 0")
  expect_error(seq_aipe_gini_sensitivity(1.5, 0.1),
               "'true_value'.*greater than 0 and less than 1; got 1.5")
  expect_error(seq_aipe_gini_sensitivity(0, 0.1), "'true_value'.*got 0")
  expect_error(seq_aipe_gini_sensitivity(1, 0.1), "'true_value'.*got 1")
  expect_error(seq_aipe_reg_coef_sensitivity(Inf, 0.5), "'true_value'")
  expect_error(seq_aipe_sm_sensitivity(NA_real_, 0.5), "'true_value'")
  # A negative or zero value is allowed where the effect size allows it.
  expect_identical(
    sens_value(seq_aipe_smd_sensitivity(-0.5, 1, G = 2, seed = 113),
               "true_value"), -0.5)
  expect_identical(
    sens_value(seq_aipe_r_sensitivity(0, 0.6, G = 2, seed = 113),
               "true_value"), 0)
  # A standardized mean and a slope may be negative or zero, and the
  # procedure still stops by its rule (not by the cap).
  for (theta in c(-0.4, 0)) {
    out <- seq_aipe_sm_sensitivity(theta, 1, G = 3, seed = 113, max_n = 500)
    expect_identical(sens_value(out, "true_value"), theta)
    expect_identical(sens_value(out, "pct_hit_max_n"), 0)
  }
  out <- seq_aipe_reg_coef_sensitivity(-0.5, 0.8, G = 3, seed = 113,
                                       max_n = 500)
  expect_identical(sens_value(out, "pct_hit_max_n"), 0)
})

test_that("a true value too extreme to simulate is an error that says so", {
  # Each of these stops before any simulation is run. The single replication
  # and the small cap are there so that a missing check would show up at once
  # as a failed expectation, not as a study that cannot end.
  quick <- function(f) {
    force(f)
    function(...) f(..., G = 1, seed = 113, max_n = 8)
  }
  cv <- quick(seq_aipe_cv_sensitivity)
  gini <- quick(seq_aipe_gini_sensitivity)

  # A coefficient of variation so large that the gamma shape 1 / cv^2
  # underflows, or that the lognormal log standard deviation
  # sqrt(log(1 + cv^2)) is infinite.
  expect_error(cv(1e200, 0.5),
               "'true_value' is too extreme.*gamma.*of 0.*got 1e\\+200")
  expect_error(
    cv(1e200, 0.5, distribution = "lognormal"),
    "'true_value' is too extreme.*lognormal.*log standard deviation.*of Inf")
  # A Gini index so near 1 that sqrt(2) qnorm((g + 1) / 2) is infinite.
  expect_error(gini(1 - 1e-16, 0.5),
               "'true_value' is too extreme.*lognormal.*qnorm.*of Inf")
  # The check itself: a positive finite parameter is returned, and anything
  # else is an error that names the distribution, the parameter, and the
  # value received.
  expect_identical(.sensitivity_check_parameter(4, "gamma", "shape", 0.5), 4)
  expect_error(
    .sensitivity_check_parameter(Inf, "gamma", "shape, 1 / true_value^2",
                                 1e-170),
    paste0("'true_value' is too extreme to simulate: the gamma distribution ",
           "would need a shape, 1 / true_value\\^2, of Inf.*got 1e-170"))
  expect_error(.sensitivity_check_parameter(0, "lognormal", "log SD", 1e-9),
               "'true_value' is too extreme.*lognormal.*of 0.*got 1e-09")
  expect_error(.sensitivity_check_parameter(NaN, "gamma", "shape", -1),
               "'true_value' is too extreme.*of NaN")

  # The slope. The variables of y = beta x + e have mean zero, so they vary
  # whatever the slope; what can be lost is the error beside beta x. From
  # 1 / .Machine$double.eps on, adjacent doubles are at least 1 apart, so
  # unit variance errors cannot be represented. That is the bound, shown here
  # from the arithmetic itself: half a unit is lost beside the bound and kept
  # beside a quarter of it.
  bound <- 1 / .Machine$double.eps
  expect_false(identical((bound + 0.5) - bound, 0.5))
  expect_identical((bound / 4 + 0.5) - bound / 4, 0.5)
  slope <- quick(seq_aipe_reg_coef_sensitivity)
  expect_error(slope(bound, 0.8),
               paste0("'true_value' must be less than 1 / ",
                      ".Machine\\$double.eps .* in absolute value; ",
                      "got 4.5.*e\\+15.*unit variance of the errors"))
  expect_error(slope(-1e20, 0.8),
               "'true_value' must be less than.*got -1e\\+20")
  expect_error(slope(1e308, 0.8),
               "'true_value' must be less than.*got 1e\\+308")
  expect_identical(.sensitivity_check_slope(bound / 2), bound / 2)
  expect_identical(sens_value(slope(1e9, 0.8), "true_value"), 1e9)

  # The contrast: standard deviations whose squares overflow or underflow.
  f <- quick(seq_aipe_c_sensitivity)
  expect_error(f(c(1, 0), c(1, -1), 1, true_sds = 1e200),
               "'true_sds' must have squares .*; 1 of 1 would overflow")
  expect_error(f(c(1, 0), c(1, -1), 1, true_sds = c(1, 1e-200)),
               "'true_sds' must have squares .*; 1 of 2 would")
})

# ---- generated data must vary by the package's own standard ----------------
#
# The check functions take a sample whose values agree to ten significant
# digits to have no variation (.is_constant(): a range of no more than
# .smsd_constant_tol, 1e-10, times the largest absolute value). Such a sample
# gives no estimate, and in a contrast its group adds nothing to the variance.
# A population whose standard deviation is that small beside its mean would
# give a study in which no replication can stop, or one whose intervals leave
# a group's variance out. The sensitivity functions refuse it, with a margin
# of 100 on the tolerance: the mean of a generated variable must be less than
# 1e8 standard deviations in absolute value.

test_that("the limit on a generated mean is 100 times clear of the tolerance", {
  expect_identical(.sensitivity_mean_limit(), 1e8)
  # At the limit the tolerance is a range of a hundredth of a standard
  # deviation: beside a mean of 1e8, values 0.009 apart are taken to be
  # constant and values 0.011 apart are not.
  expect_true(.is_constant(1e8 + c(0, 0.009)))
  expect_false(.is_constant(1e8 + c(0, 0.011)))
  # How often a pilot sample of 4, the smallest there is, has a range that
  # small: P(range <= r) = n * integral of f(x) (F(x + r) - F(x))^(n - 1) dx
  # for n independent values with density f and distribution function F. It
  # is about 1.3 in ten million for normal data and 3.7 in ten million for
  # the skewed mechanism (chi square with 3 degrees of freedom, standardized),
  # so the tolerance plays no part in a study of any practical size. Without
  # the margin (r = 1) the same integral gives 0.11 for normal data.
  r <- .smsd_constant_tol * .sensitivity_mean_limit()
  expect_equal(r, 0.01)
  p_range <- function(r, dens, cdf, lower) {
    4 * stats::integrate(function(x) dens(x) * (cdf(x + r) - cdf(x))^3,
                         lower, Inf, rel.tol = 1e-10)$value
  }
  expect_lt(p_range(r, stats::dnorm, stats::pnorm, -Inf), 2e-7)
  expect_lt(
    p_range(r, function(x) sqrt(6) * stats::dchisq(3 + sqrt(6) * x, df = 3),
            function(x) stats::pchisq(3 + sqrt(6) * x, df = 3), -3 / sqrt(6)),
    5e-7)
  # The same, by drawing: of 20,000 normal pilot samples of 4 just inside
  # the limit, none is taken to be constant (0.003 are expected); at a mean
  # of 1e10 standard deviations, which the margin excludes, about one in ten
  # is.
  pilots <- withr::with_seed(113, matrix(stats::rnorm(80000), nrow = 4))
  expect_identical(sum(apply(9.9e7 + pilots, 2, .is_constant)), 0L)
  expect_gt(sum(apply(1e10 + pilots, 2, .is_constant)), 1000)
})

test_that("a population whose data would not vary enough is refused", {
  quick <- function(f) {
    force(f)
    function(...) f(..., G = 1, seed = 113, max_n = 8)
  }
  why <- "ten significant digits"
  # The standardized mean difference and the standardized mean: unit
  # variance data with a mean of true_value.
  for (g in list(seq_aipe_smd_sensitivity, seq_aipe_sm_sensitivity)) {
    f <- quick(g)
    expect_error(f(1e8, 0.8),
                 paste0("'true_value' must be less than 1e\\+08 in absolute ",
                        "value; got 1e\\+08.*", why))
    expect_error(f(-1e11, 0.8),
                 "'true_value' must be less than 1e\\+08.*got -1e\\+11")
    expect_error(f(1e308, 0.8), "'true_value' must be less than 1e\\+08")
    expect_identical(sens_value(f(9.9e7, 0.8), "true_value"), 9.9e7)
    expect_identical(sens_value(f(-9.9e7, 0.8), "true_value"), -9.9e7)
  }
  expect_identical(.sensitivity_check_location(5e7), 5e7)

  # The contrast: each mean against its own standard deviation.
  f <- quick(seq_aipe_c_sensitivity)
  expect_error(f(c(1e8, 0), c(1, -1), 1),
               paste0("'true_means' must each be less than 1e\\+08 standard ",
                      "deviations .*; 1 of 2 is not.*", why))
  expect_error(f(c(1e20, -1e20, 0), c(1, -1, 0), 1), "2 of 3 are not")
  expect_error(f(c(1e308, 1e308), c(1, 1), 1), "'true_means'.*2 of 2 are not")
  # (Of two means of 3e8, the one with a standard deviation of 4 is within
  # the limit and the one with a standard deviation of 1 is not, whichever
  # comes first.)
  expect_error(f(c(3e8, 3e8), c(1, -1), 1, true_sds = c(4, 1)),
               "'true_means' must each be less than .*; 1 of 2 is not")
  expect_error(f(c(3e8, 3e8), c(1, -1), 1, true_sds = c(1, 4)),
               "'true_means' must each be less than .*; 1 of 2 is not")
  expect_identical(
    sens_value(f(c(3e8, 3e8), c(1, -1), 30, true_sds = c(4, 4)), "true_value"),
    0)
  # A mean of 1e20 is refused with a standard deviation of 1e12 and accepted
  # with one of 1e13, and a mean of 1 is refused with one of 1e-8.
  expect_error(f(c(1e20, 0), c(1, -1), 1e13, true_sds = 1e12),
               "'true_means' must each be less than .*; 1 of 2 is not")
  out <- f(c(1e20, 0), c(1, -1), 1e14, true_sds = 1e13)
  expect_identical(sens_value(out, "true_value"), 1e20)
  expect_error(f(c(1, 0), c(1, -1), 1, true_sds = 1e-8),
               "'true_means' must each be less than .*; 1 of 2 is not")
  expect_identical(
    .sensitivity_check_location(c(1e20, 0), c(1e13, 1e13), "true_means"),
    c(1e20, 0))

  # The coefficient of variation is the standard deviation of the data over
  # their mean, so it must be greater than 1e-8, for either distribution.
  cv <- quick(seq_aipe_cv_sensitivity)
  for (dist in c("gamma", "lognormal")) {
    expect_error(cv(1e-8, 0.5, distribution = dist),
                 paste0("'true_value' must be greater than 1e-08; got 1e-08.*",
                        why))
    expect_error(cv(1e-170, 0.5, distribution = dist),
                 "'true_value' must be greater than 1e-08; got 1e-170")
    expect_identical(
      sens_value(cv(2e-8, 0.5, distribution = dist), "true_value"), 2e-8)
  }

  # The Gini index: lognormal data with log standard deviation s have a
  # standard deviation of sqrt(exp(s^2) - 1) times their mean, about s when
  # s is small, and a Gini index of 2 pnorm(s / sqrt(2)) - 1, about
  # s / sqrt(pi). The smallest Gini index is therefore about
  # 1e-8 / sqrt(pi) = 5.6e-9.
  expect_equal(.sensitivity_gini_lognormal(.sensitivity_sdlog_cv(1e-8)),
               1e-8 / sqrt(pi), tolerance = 1e-6)
  gini <- quick(seq_aipe_gini_sensitivity)
  expect_error(gini(5.6e-9, 0.5),
               paste0("'true_value' must be greater than about 5.6e-09 \\(the ",
                      "Gini index of lognormal data with a coefficient of ",
                      "variation of 1e-08\\); got 5.6e-09.*", why))
  expect_error(gini(1e-17, 0.5),
               "'true_value' must be greater than about 5.6e-09.*got 1e-17")
  expect_identical(sens_value(gini(5.7e-9, 0.5), "true_value"), 5.7e-9)
})

test_that("a study just inside the limit is the study it is anywhere else", {
  # A contrast, its variance, and its interval do not depend on where a
  # group's mean lies: adding 9.9e7 standard deviations to the first mean
  # moves the estimate and the true value alike. Over the same random
  # stream the two studies must therefore agree: the same sample sizes and
  # the same covering intervals, and widths equal up to the rounding of
  # unit variance values stored beside 9.9e7 (their last digit is 1.5e-8).
  # A check function that took such a group to have no variation would stop
  # earlier, with narrower intervals that cover less often.
  w <- c(1, -0.5, -0.5)
  study <- function(means) {
    seq_aipe_c_sensitivity(means, w, width = 2, G = 60, seed = 113,
                           max_n = 40)
  }
  far <- study(c(9.9e7, 0, 0))
  near <- study(c(0, 0, 0))
  expect_identical(sens_value(far, "true_value"), 9.9e7)
  rows <- c("mean_n", "median_n", "sd_n", "pct_width_within", "coverage",
            "pct_hit_max_n", "pct_variance_unusable")
  expect_identical(sens_value(far, rows), sens_value(near, rows))
  expect_equal(sens_value(far, "mean_width"), sens_value(near, "mean_width"),
               tolerance = 1e-6)
  # Not vacuous: sampling went past the pilot, and not every interval
  # covers.
  expect_gt(sens_value(near, "mean_n"), seq_aipe_c_pilot(width = 2))
  expect_gt(sens_value(near, "sd_n"), 0)
  expect_lt(sens_value(near, "coverage"), 1)
})

test_that("a distribution may be abbreviated; the default is the first", {
  full <- seq_aipe_sm_sensitivity(0.4, 1, distribution = "skewed", G = 3,
                                  seed = 113)
  expect_identical(
    seq_aipe_sm_sensitivity(0.4, 1, distribution = "sk", G = 3, seed = 113),
    full)
  expect_identical(
    seq_aipe_sm_sensitivity(0.4, 1, G = 3, seed = 113),
    seq_aipe_sm_sensitivity(0.4, 1, distribution = "normal", G = 3,
                            seed = 113))
  expect_false(identical(
    full$value,
    seq_aipe_sm_sensitivity(0.4, 1, G = 3, seed = 113)$value))
  expect_identical(
    seq_aipe_cv_sensitivity(0.5, 0.5, G = 3, seed = 113),
    seq_aipe_cv_sensitivity(0.5, 0.5, distribution = "gamma", G = 3,
                            seed = 113))
})

test_that("a required argument that is left out is asked for in the package's words", {
  # R's own message ('argument "true_value" is missing, with no default')
  # names the argument but not what it is.
  what <- c(
    seq_aipe_smd_sensitivity = "standardized mean difference",
    seq_aipe_r_sensitivity = "correlation",
    seq_aipe_cv_sensitivity = "coefficient of variation",
    seq_aipe_reg_coef_sensitivity = "regression slope",
    seq_aipe_sm_sensitivity = "standardized mean",
    seq_aipe_gini_sensitivity = "Gini index"
  )
  for (name in names(what)) {
    f <- get(name, envir = asNamespace("SMSD"))
    expect_error(
      f(width = 0.4, G = 1, seed = 113),
      paste0("'true_value' must be supplied: the population ", what[[name]],
             " to simulate."),
      fixed = TRUE, label = name
    )
    # Asked for before anything else is looked at.
    expect_error(f(width = -1, G = 0), "'true_value' must be supplied",
                 label = name)
    # The width is asked for as everywhere else.
    expect_error(f(0.3, G = 1, seed = 113), "'width' must be supplied",
                 label = name)
  }
  expect_error(
    seq_aipe_c_sensitivity(contrast_weights = c(1, -1), width = 1, G = 1),
    paste("'true_means' must be supplied: the population means to simulate,",
          "one per group."), fixed = TRUE
  )
  expect_error(
    seq_aipe_c_sensitivity(true_means = c(1, 0), width = 1, G = 1),
    "'contrast_weights' must be supplied: one contrast weight per group.",
    fixed = TRUE
  )
  expect_error(seq_aipe_c_sensitivity(c(1, 0), c(1, -1), G = 1),
               "'width' must be supplied")
  # Every required argument of every sensitivity function of this family,
  # left out one at a time with the others valid: none falls through to R.
  valid <- list(true_value = 0.4, width = 0.6, true_means = c(1, 0, 0),
                contrast_weights = c(1, -0.5, -0.5), p = 2)
  functions <- grep("^seq_aipe_.*_sensitivity$", getNamespaceExports("SMSD"),
                    value = TRUE)
  expect_length(functions, 8L)
  for (name in functions) {
    f <- get(name, envir = asNamespace("SMSD"))
    required <- names(formals(f))[vapply(formals(f), function(a) {
      is.symbol(a) && !nzchar(as.character(a))
    }, logical(1))]
    expect_true(all(required %in% names(valid)), label = name)
    for (omit in required) {
      expect_error(
        do.call(f, c(valid[setdiff(required, omit)], list(G = 1, seed = 113))),
        paste0("'", omit, "' must be supplied"), fixed = TRUE,
        label = paste(name, "without", omit)
      )
    }
  }
})

test_that("names on the contrast arguments must agree in order", {
  # Entries are matched to the groups by position. Weights named in another
  # order than the means would be applied to the wrong groups without a
  # word: here 75 * (-0.5) + 70 * (-0.5) + 70 * 1 = -2.5 in place of the
  # contrast the names describe, 75 - 35 - 35 = 5.
  f <- function(...) seq_aipe_c_sensitivity(..., width = 8, G = 2, seed = 113)
  mu <- c(in_person = 75, video = 70, hybrid = 70)
  expect_error(
    f(mu, c(video = -0.5, hybrid = -0.5, in_person = 1), true_sds = 10),
    paste0("'contrast_weights' is named \\(video, hybrid, in_person\\) but ",
           "'true_means' is named \\(in_person, video, hybrid\\)\\. Entries ",
           "are matched to groups by position, not by name")
  )
  # In the order of the means the names are accepted, and the true value
  # is the contrast they describe.
  w <- c(in_person = 1, video = -0.5, hybrid = -0.5)
  expect_identical(sens_value(f(mu, w, true_sds = 10), "true_value"), 5)
  # A single standard deviation is shared by every group, so its name is
  # not compared; one per group is compared like the rest.
  expect_identical(sens_value(f(mu, w, true_sds = c(sd = 10)), "true_value"),
                   5)
  expect_error(f(mu, w, true_sds = c(video = 10, in_person = 10, hybrid = 10)),
               "'true_sds' is named \\(video, in_person, hybrid\\) but 'true_means' is named")
  # Names on one argument only, or on none, are matched by position, as
  # documented.
  expect_identical(sens_value(f(mu, c(-0.5, -0.5, 1), true_sds = 10),
                              "true_value"), -2.5)
  expect_identical(sens_value(f(unname(mu), w, true_sds = 10), "true_value"),
                   5)
  # The names are compared before anything is drawn.
  next_draw <- withr::with_seed(113, {
    try(seq_aipe_c_sensitivity(mu, rev(w), width = 8, G = 1), silent = TRUE)
    stats::runif(1)
  })
  expect_identical(next_draw, withr::with_seed(113, stats::runif(1)))
})

test_that("the contrast arguments are checked", {
  # One replication and a small cap, so that a missing check would show up
  # at once as a failed expectation, not as a long study.
  f <- function(...) seq_aipe_c_sensitivity(..., G = 1, seed = 113, max_n = 8)
  w <- c(1, -0.5, -0.5)
  mu <- c(1, 0, 0)
  # true_means: a numeric vector of finite values, at least two groups.
  expect_error(f("a", w, 1),
               paste0("'true_means' must be a numeric vector \\(one ",
                      "population mean per group\\); got \"a\""))
  expect_error(f(matrix(mu, 1), w, 1), "'true_means'.*dimensions 1 x 3")
  expect_error(f(numeric(0), w, 1), "'true_means'.*length 0")
  expect_error(f(NULL, w, 1), "'true_means' must be a numeric vector.*NULL")
  expect_error(f(c(1, NA, 0), w, 1),
               "'true_means' must contain only finite numbers.*entry 2 is NA")
  expect_error(f(c(1, 0, Inf), w, 1),
               "'true_means' must contain only finite numbers.*entry 3 is Inf")
  expect_error(f(1, 1, 1), "'true_means'.*at least 2 groups")
  # contrast_weights: one finite weight per group, not all zero.
  expect_error(f(mu, c(1, -1), 1),
               paste0("'contrast_weights' must have 3 entries \\(one per ",
                      "group in 'true_means'\\); it has 2"))
  expect_error(f(mu, c(1, -1, 0, 0), 1),
               "'contrast_weights' must have 3 entries.*it has 4")
  # (The check function would refuse all-zero weights too, in other words,
  # but only once a pilot sample had been drawn: this function refuses them
  # itself, before it draws anything, so the caller's random number stream
  # has not moved.)
  expect_error(f(mu, c(0, 0, 0), 1),
               paste0("'contrast_weights' must not all be zero: a contrast ",
                      "needs at least one group with a nonzero weight"))
  next_draw <- withr::with_seed(113, {
    try(seq_aipe_c_sensitivity(mu, c(0, 0, 0), 1, G = 1, max_n = 8),
        silent = TRUE)
    stats::runif(1)
  })
  expect_identical(next_draw, withr::with_seed(113, stats::runif(1)))
  expect_error(
    f(mu, c(1, NA, -1), 1),
    "'contrast_weights' must contain only finite numbers.*entry 2 is NA")
  expect_error(
    f(mu, c(1, Inf, -1), 1),
    "'contrast_weights' must contain only finite numbers.*entry 2 is Inf")
  expect_error(f(mu, "w", 1), "'contrast_weights' must be a numeric vector")
  expect_error(f(mu, matrix(w, 1), 1),
               "'contrast_weights'.*dimensions 1 x 3")
  # true_sds: one positive finite value, or one per group.
  expect_error(f(mu, w, 1, true_sds = c(1, 2)),
               paste0("'true_sds' must have 3 entries \\(one per group in ",
                      "'true_means', or else a single entry shared by every ",
                      "group\\); it has 2"))
  expect_error(f(mu, w, 1, true_sds = c(1, 2, 3, 4)),
               "'true_sds' must have 3 entries.*it has 4")
  expect_error(f(mu, w, 1, true_sds = -1),
               "'true_sds' must contain only positive numbers.*entry 1 is -1")
  expect_error(f(mu, w, 1, true_sds = c(1, 0, 0)),
               "'true_sds' must contain only positive numbers.*entry 2 is 0")
  expect_error(f(mu, w, 1, true_sds = NA_real_),
               "'true_sds' must contain only finite numbers.*entry 1 is NA")
  expect_error(f(mu, w, 1, true_sds = "a"),
               "'true_sds' must be a numeric vector")
  expect_error(f(mu, w, 1, true_sds = matrix(1)),
               "'true_sds'.*dimensions 1 x 1")
  expect_error(f(mu, w, 1, true_sds = numeric(0)), "'true_sds'.*length 0")
  # One value is shared by every group; one per group is used as given. With
  # a cap equal to the pilot each call is a single check of 4 per group, so
  # the two calls draw the same standard normal numbers and differ only in
  # the standard deviations: the interval for (2, 2, 2) is twice as wide as
  # for the default of 1, and (2, 2, 2) is the same study as a single 2.
  g <- function(...) {
    seq_aipe_c_sensitivity(mu, w, width = 1, ..., G = 3, seed = 113,
                           max_n = 4)
  }
  expect_identical(g(true_sds = c(2, 2, 2)), g(true_sds = 2))
  expect_equal(sens_value(g(true_sds = 2), "mean_width"),
               2 * sens_value(g(), "mean_width"), tolerance = 1e-12)
  expect_false(identical(g(true_sds = c(2, 1, 1))$value,
                         g(true_sds = 2)$value))
  expect_error(f(mu, w, 1, distribution = "gamma"),
               "'distribution' must be one of \"normal\", \"skewed\"")
  # Means in range whose contrast overflows.
  expect_error(f(c(1, 1), c(1e308, 1e308), 1),
               "contrast implied by 'true_means' and 'contrast_weights' is not")
})
