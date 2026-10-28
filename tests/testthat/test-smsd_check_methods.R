# The smsd_check class: constructor contract and methods.

make_check <- function(...) {
  new_smsd_check(
    term = c("lower_limit", "smd", "upper_limit", "xi_squared",
             "criterion_n", "current_n"),
    value = c(0.1, 0.5, 0.9, 2.03, 120, 40),
    method = "standardized mean difference",
    is_satisfied = FALSE, current_n = 40, criterion_n = 120,
    conf_level = 0.95, width = 0.5, ...
  )
}

test_that("the constructor refuses an NA stopping flag", {
  expect_error(
    new_smsd_check("a", 1, "m", is_satisfied = NA, current_n = 1,
                   criterion_n = 1),
    "is.na"
  )
})

test_that("the constructor refuses a row that is Inf or NaN, and allows NA", {
  make <- function(value) {
    new_smsd_check(c("mean", "sd", "risk"), value, "m", is_satisfied = FALSE,
                   current_n = 3, criterion_n = NA_real_)
  }
  expect_identical(make(c(1, 0, NA))$value, c(1, 0, NA_real_))
  expect_identical(make(c(-1e308, 0, 1e308))$value, c(-1e308, 0, 1e308))
  expect_error(
    make(c(1, Inf, NA)),
    paste0("^The 'sd' row of the result is too large to represent \\(Inf\\)\\. ",
           "Express the data in larger units \\(divide them by a constant\\) ",
           "and state the precision target in the same units\\.$")
  )
  expect_error(make(c(-Inf, 2, Inf)),
               "^The 'mean' row of the result is too large to represent \\(-Inf\\)\\.")
  expect_error(make(c(1, 2, NaN)),
               "^The 'risk' row of the result is too large to represent \\(NaN\\)\\.")
})

test_that("notes default to empty and are printed when present", {
  out <- make_check()
  expect_identical(attr(out, "notes"), character(0))
  expect_false(any(grepl("Note:", format(out))))
  noted <- make_check(notes = .note_variance_unusable(),
                      extra = list(variance_unusable = TRUE))
  expect_true(attr(noted, "variance_unusable"))
  txt <- format(noted)
  expect_true(any(grepl("Note: The variance estimate is not yet usable",
                        txt)))
  expect_true(any(grepl("Stopping rule met: FALSE", txt)))
})

test_that("print writes exactly what format returns, invisibly", {
  out <- make_check()
  expect_identical(capture.output(res <- print(out)), format(out))
  expect_identical(res, out)
  invisible(capture.output(expect_invisible(print(out))))
})

test_that("digits affects display only", {
  out <- make_check()
  expect_true(any(grepl("2.03", format(out, digits = 3), fixed = TRUE)))
  expect_identical(out$value[4], 2.03)
})

test_that("summary and glance describe the check in one row", {
  out <- make_check()
  s <- summary(out)
  expect_s3_class(s, "data.frame")
  expect_identical(nrow(s), 1L)
  expect_identical(s$estimate_term, "smd")
  expect_identical(s$estimate, 0.5)
  expect_identical(s$lower_limit, 0.1)
  expect_identical(s$upper_limit, 0.9)
  expect_identical(s$current_n, 40)
  expect_identical(s$criterion_n, 120)
  expect_false(s$is_satisfied)
  expect_identical(s$conf_level, 0.95)
  expect_identical(s$width, 0.5)
  expect_identical(generics::glance(out), s)
})

test_that("summary reports NA for quantities a risk procedure lacks", {
  risk <- new_smsd_check(c("cv", "risk", "optimal_n", "current_n"),
                         c(0.4, 12, 30, 30), "minimum risk (cv)",
                         is_satisfied = TRUE, current_n = 30,
                         criterion_n = 30)
  s <- summary(risk)
  expect_identical(s$estimate_term, "cv")
  expect_true(is.na(s$lower_limit) && is.na(s$upper_limit))
  expect_true(is.na(s$conf_level) && is.na(s$width))
  expect_true(s$is_satisfied)
})

test_that("tidy uses the DMAR column layout", {
  td <- generics::tidy(make_check())
  expect_identical(names(td),
                   c("term", "estimate", "ci_lower", "ci_upper",
                     "conf_level"))
  expect_identical(td$term, "smd")
  expect_identical(c(td$ci_lower, td$estimate, td$ci_upper),
                   c(0.1, 0.5, 0.9))
})

test_that("as.data.frame returns the plain term and value frame", {
  df <- as.data.frame(make_check())
  expect_identical(class(df), "data.frame")
  expect_identical(names(df), c("term", "value"))
  expect_null(attr(df, "is_satisfied"))
})

test_that("a partial subset is a plain data frame; a full one keeps class", {
  out <- make_check()
  part <- out[1:2, ]
  expect_identical(class(part), "data.frame")
  expect_null(attr(part, "is_satisfied"))
  expect_null(attr(part, "method"))
  expect_s3_class(out[seq_len(nrow(out)), ], "smsd_check")
  expect_identical(out[out$term == "smd", "value"], 0.5)
  expect_identical(out$value[out$term == "smd"], 0.5)
  expect_identical(class(out[, "value", drop = FALSE]), "data.frame")
})

# The exact lines print() writes. Each expected vector below was written out
# from the documented layout (?smsd_check_methods) and the arguments given to
# the constructor, not copied from the method's own output. The display
# options are fixed so the lines do not depend on the session.
test_that("the printed lines of an interval check are exactly these", {
  withr::local_options(list(digits = 7, scipen = 0, OutDec = "."))
  expect_identical(
    format(make_check()),
    c("Sequential check (standardized mean difference)",
      "",
      "  Current n       : 40",
      "  Criterion n     : 120",
      "  Target width    : 0.5",
      "  Confidence level: 95%",
      "",
      "  lower_limit       : 0.1",
      "  smd               : 0.5",
      "  upper_limit       : 0.9",
      "  xi_squared        : 2.03",
      "  criterion_n       : 120",
      "  current_n         : 40",
      "",
      "  Stopping rule met: FALSE")
  )
})

test_that("a risk check prints no width or confidence line, and its note", {
  withr::local_options(list(digits = 7, scipen = 0, OutDec = "."))
  risk <- new_smsd_check(c("cv", "risk", "optimal_n", "current_n"),
                         c(0.4, 12.5, 28, 30), "minimum risk (cv)",
                         is_satisfied = TRUE, current_n = 30,
                         criterion_n = 28, notes = "A short note.")
  expect_identical(
    format(risk),
    c("Sequential check (minimum risk (cv))",
      "",
      "  Current n       : 30",
      "  Criterion n     : 28",
      "",
      "  cv                : 0.4",
      "  risk              : 12.5",
      "  optimal_n         : 28",
      "  current_n         : 30",
      "",
      "  Stopping rule met: TRUE",
      "",
      "  Note: A short note.")
  )
})

test_that("digits shortens the target width and every value shown", {
  withr::local_options(list(digits = 7, scipen = 0, OutDec = "."))
  long <- new_smsd_check(c("lower_limit", "r", "upper_limit", "xi_squared"),
                         c(0.123456789, 0.2345678, 0.3456789, 1.23456789),
                         "correlation", is_satisfied = FALSE, current_n = 25,
                         criterion_n = 81, conf_level = 0.9,
                         width = 0.123456)
  expect_identical(
    format(long, digits = 3),
    c("Sequential check (correlation)",
      "",
      "  Current n       : 25",
      "  Criterion n     : 81",
      "  Target width    : 0.123",
      "  Confidence level: 90%",
      "",
      "  lower_limit       : 0.123",
      "  r                 : 0.235",
      "  upper_limit       : 0.346",
      "  xi_squared        : 1.23",
      "",
      "  Stopping rule met: FALSE")
  )
  expect_identical(capture.output(print(long, digits = 3)),
                   format(long, digits = 3))
  # The default is the session's 'digits' option (7 here).
  expect_true("  lower_limit       : 0.1234568" %in% format(long))
  expect_true("  Target width    : 0.123456" %in% format(long))
  expect_identical(long$value[1], 0.123456789)
})

test_that("a budget check ended by exhaustion says so on the stopping line", {
  # The budget constrained procedure ends in two ways, and is_satisfied is
  # TRUE for both. When the budget ended the sampling with a group still
  # below its criterion (budget_exhausted TRUE, rule_met FALSE), the line
  # must not read "TRUE" above a criterion the sample has not reached. The
  # objects are built from the documented attributes of ?seq_budget, so the
  # expected lines are independent of seq_budget() itself.
  withr::local_options(list(digits = 7, scipen = 0, OutDec = "."))
  budget <- function(rule_met, budget_exhausted, is_satisfied) {
    new_smsd_check(
      c("lower_limit", "contrast", "upper_limit", "variance",
        "current_total_n", "current_total_cost"),
      c(-1.5, 2, 5.5, 3.2, 399, 2984),
      "budget constrained K-population contrast",
      is_satisfied = is_satisfied, current_n = 399, criterion_n = 400,
      conf_level = 0.95,
      extra = list(rule_met = rule_met, budget_exhausted = budget_exhausted)
    )
  }
  exhausted <- budget(FALSE, TRUE, TRUE)
  expect_identical(
    format(exhausted),
    c("Sequential check (budget constrained K-population contrast)",
      "",
      "  Current n       : 399",
      "  Criterion n     : 400",
      "  Confidence level: 95%",
      "",
      "  lower_limit       : -1.5",
      "  contrast          : 2",
      "  upper_limit       : 5.5",
      "  variance          : 3.2",
      "  current_total_n   : 399",
      "  current_total_cost: 2984",
      "",
      "  Stopping rule met: FALSE (sampling has ended: the budget is exhausted)")
  )
  # Display only: the flag that drives a sampling loop is untouched.
  expect_true(attr(exhausted, "is_satisfied"))
  expect_true(summary(exhausted)$is_satisfied)
  last_line <- function(x) {
    lines <- format(x)
    lines[length(lines)]
  }
  # Every group at its criterion; sampling still under way; and the budget
  # exhausted with a variance that cannot be used (the flag is then FALSE).
  expect_identical(last_line(budget(TRUE, FALSE, TRUE)),
                   "  Stopping rule met: TRUE")
  expect_identical(last_line(budget(FALSE, FALSE, FALSE)),
                   "  Stopping rule met: FALSE")
  # The wording is for the ending by exhaustion alone: rule_met decides.
  expect_identical(last_line(budget(TRUE, TRUE, TRUE)),
                   "  Stopping rule met: TRUE")
  expect_identical(
    last_line(budget(FALSE, TRUE, FALSE)),
    "  Stopping rule met: FALSE (sampling has ended: the budget is exhausted)"
  )
})

test_that("checks bound together are refused by every method that reads one", {
  # rbind() of checks keeps the class and the attributes of the first, so a
  # log of stages kept that way would print, summarize, and be reported as
  # its first stage alone. Two checks with known contents: not met at n = 40,
  # met at n = 120.
  first <- make_check()
  second <- new_smsd_check(
    term = first$term, value = c(0.25, 0.5, 0.75, 2.03, 120, 120),
    method = "standardized mean difference", is_satisfied = TRUE,
    current_n = 120, criterion_n = 120, conf_level = 0.95, width = 0.5
  )
  refused <- paste0(
    "^This object has class 'smsd_check' but is no longer a check as the ",
    "check function returned it \\(the term lower_limit appears %d times, ",
    "so it holds more than one check\\)\\. It was changed afterward, for ",
    "example by binding checks together or by changing an attribute\\. Use ",
    "as\\.data\\.frame\\(\\) for the plain table of terms; to keep a ",
    "record of the stages of a study, bind the one row summary\\(\\) of ",
    "each check, or use seq_path\\(\\)\\.$"
  )
  bound <- list(
    two = rbind(first, second),
    reversed = rbind(second, first),
    # The loop a user writes: log <- NULL; log <- rbind(log, check).
    logged = Reduce(rbind, list(first, second, first), NULL)
  )
  copies <- c(two = 2L, reversed = 2L, logged = 3L)
  for (nm in names(bound)) {
    msg <- sprintf(refused, copies[[nm]])
    expect_error(format(bound[[nm]]), msg)
    expect_error(print(bound[[nm]]), msg)
    expect_error(summary(bound[[nm]]), msg)
    expect_error(generics::tidy(bound[[nm]]), msg)
    expect_error(generics::glance(bound[[nm]]), msg)
    expect_error(seq_report(bound[[nm]]), msg)
    # The plain table is always available, with every row.
    plain <- as.data.frame(bound[[nm]])
    expect_identical(class(plain), "data.frame")
    expect_identical(nrow(plain), 6L * copies[[nm]])
    expect_null(attr(plain, "is_satisfied"))
  }
  expect_identical(as.data.frame(bound$two)$value,
                   c(first$value, second$value))
  # One check bound to nothing is that check.
  expect_identical(format(rbind(NULL, first)), format(first))
  # The documented way to keep a log: one summary() row per check.
  log <- rbind(summary(first), summary(second))
  expect_identical(log$current_n, c(40, 120))
  expect_identical(log$is_satisfied, c(FALSE, TRUE))
})

test_that("a check whose attributes or columns were changed is refused", {
  intact <- make_check()
  expect_null(.smsd_check_problem(intact))
  expect_null(.smsd_check_problem(intact[seq_len(nrow(intact)), ]))
  # A criterion that is not defined is NA by design, and is still a check.
  no_criterion <- new_smsd_check("median", 3, "m", is_satisfied = FALSE,
                                 current_n = 5, criterion_n = NA_real_)
  expect_null(.smsd_check_problem(no_criterion))
  expect_identical(summary(no_criterion)$criterion_n, NA_real_)

  no_flag <- intact
  attr(no_flag, "is_satisfied") <- NULL
  expect_error(
    summary(no_flag),
    "no longer a check .* \\(the attribute is_satisfied is missing or not as the check function set it\\)"
  )
  expect_error(print(no_flag), "the attribute is_satisfied is missing")
  na_flag <- intact
  attr(na_flag, "is_satisfied") <- NA
  expect_error(format(na_flag), "the attribute is_satisfied is missing or not")
  two_gone <- intact
  attr(two_gone, "current_n") <- NULL
  attr(two_gone, "criterion_n") <- "120"
  expect_error(
    generics::glance(two_gone),
    "\\(the attributes current_n, criterion_n are missing or not as the check function set them\\)"
  )
  no_method <- intact
  attr(no_method, "method") <- NULL
  expect_error(seq_report(no_method), "the attribute method is missing")
  text_values <- intact
  text_values$value <- as.character(text_values$value)
  expect_error(
    generics::tidy(text_values),
    "\\(its columns are not 'term' \\(character\\) and 'value' \\(numeric\\)\\)"
  )
  # A value changed in place leaves one check, which the methods still read.
  rounded <- intact
  rounded$value <- round(rounded$value, 1)
  expect_identical(summary(rounded)$estimate, 0.5)
})

test_that("what the limit rows hold when no variance can be estimated", {
  # The statement on ?smsd_check_methods, family by family, on data with no
  # variation (or, for the correlation, the slope, and R2, collinear data,
  # whose variance estimate is zero): every check is flagged and not met;
  # the accuracy in parameter estimation functions report NA limits; the
  # fixed width functions report limits that are the estimate and d; and
  # seq_budget() reports limits computed from the understated variance.
  flat <- rep(5, 40)
  line <- 1:40 + 0
  withr::local_seed(113)
  noise <- rnorm(40)
  three <- list(rep(1, 10), rep(2, 10), rep(3, 10))
  w <- c(1, -0.5, -0.5)
  limits <- function(out) {
    out$value[match(c("lower_limit", "upper_limit"), out$term)]
  }
  aipe <- list(
    seq_aipe_sm = seq_aipe_sm(flat, width = 0.5),
    seq_aipe_smd = seq_aipe_smd(flat, flat + 1, width = 0.5),
    seq_aipe_cv = seq_aipe_cv(flat, width = 0.5),
    seq_aipe_gini = seq_aipe_gini(flat, width = 0.5),
    seq_aipe_r = seq_aipe_r(line, 2 * line + 1, width = 0.5),
    seq_aipe_reg_coef = seq_aipe_reg_coef(line, 2 * line + 1, width = 0.5),
    seq_aipe_R2 = seq_aipe_R2(cbind(line, noise), 2 * line + 1, width = 0.2),
    seq_aipe_c = seq_aipe_c(three, w, width = 1)
  )
  for (nm in names(aipe)) {
    # For a sample correlation of exactly 1, whether the cause is recorded
    # as the variance or as the estimate is the correlation function's
    # business; the limits are NA and the rule is not met either way.
    if (nm != "seq_aipe_r") {
      expect_true(attr(aipe[[nm]], "variance_unusable"), info = nm)
    }
    expect_false(attr(aipe[[nm]], "is_satisfied"), info = nm)
    expect_identical(limits(aipe[[nm]]), c(NA_real_, NA_real_), info = nm)
    td <- generics::tidy(aipe[[nm]])
    expect_identical(c(td$ci_lower, td$ci_upper), c(NA_real_, NA_real_),
                     info = nm)
  }
  fwd <- list(
    seq_fwd_norm_mean = list(seq_fwd_norm_mean(flat, d = 0.5), c(4.5, 5.5)),
    seq_fwd_dfr_mean = list(seq_fwd_dfr_mean(flat, d = 0.5), c(4.5, 5.5)),
    # One sided: the sample minimum minus d, and the sample minimum.
    seq_fwd_nexp_mean = list(seq_fwd_nexp_mean(flat, d = 0.5), c(4.5, 5))
  )
  for (nm in names(fwd)) {
    out <- fwd[[nm]][[1L]]
    expect_true(attr(out, "variance_unusable"), info = nm)
    expect_false(attr(out, "is_satisfied"), info = nm)
    expect_identical(limits(out), fwd[[nm]][[2L]], info = nm)
    # tidy() carries the limits and nothing that marks them: the flag is
    # the thing to test, as the help page says.
    td <- generics::tidy(out)
    expect_identical(c(td$ci_lower, td$ci_upper), fwd[[nm]][[2L]], info = nm)
    expect_false("variance_unusable" %in% names(td), info = nm)
  }
  # Every group constant: the contrast is 1 - (2 + 3) / 2 = -1.5 and the
  # variance row is zero, so both limits equal the estimate.
  all_flat <- seq_budget(three, w, cost_per_group = c(38, 15, 10),
                         budget = 1500)
  expect_true(attr(all_flat, "variance_unusable"))
  expect_false(attr(all_flat, "is_satisfied"))
  expect_identical(limits(all_flat), c(-1.5, -1.5))
  # One group constant: the variance leaves that group out. By hand, with
  # the 97.5 percent normal quantile, from the two groups that vary.
  g1 <- c(4, 6, 5, 7, 3, 5, 6, 4, 5, 5)
  g3 <- c(2, 4, 3, 3, 1, 5, 3, 2, 4, 3)
  one_flat <- seq_budget(list(g1, rep(2, 10), g3), w,
                         cost_per_group = c(38, 15, 10), budget = 1500)
  expect_true(attr(one_flat, "variance_unusable"))
  expect_false(attr(one_flat, "is_satisfied"))
  estimate <- mean(g1) - 0.5 * 2 - 0.5 * mean(g3)
  half <- stats::qnorm(0.975) * sqrt(var(g1) / 10 + 0.25 * var(g3) / 10)
  expect_equal(limits(one_flat), c(estimate - half, estimate + half),
               tolerance = 1e-12)
})

test_that("no check function returns Inf or NaN in the value column", {
  # The contract of the check object: a row is a finite number, or NA where
  # the quantity is documented as not defined. Data at the edges of double
  # precision either give such an object or stop with one of the package's
  # own messages; they never put Inf or NaN in a row. Every check function
  # is run on ordinary data and on the same data moved to those edges.
  withr::local_seed(113)
  x <- rnorm(40, 10, 2)
  pos <- rgamma(40, shape = 16, rate = 20)
  y <- 0.5 * x + rnorm(40)
  z <- rnorm(40)
  shapes <- list(
    as_is = function(v) v,
    times_1e150 = function(v) v * 1e150,
    times_1e200 = function(v) v * 1e200,
    times_1e300 = function(v) v * 1e300,
    times_1e_200 = function(v) v * 1e-200,
    times_1e_300 = function(v) v * 1e-300,
    plus_1e15 = function(v) v + 1e15,
    one_at_1e300 = function(v) replace(v, 1L, 1e300),
    largest_double = function(v) {
      rep_len(.Machine$double.xmax * c(1, 0.5), length(v))
    },
    # From one end of the range to the other: a difference of two of these
    # (the length of an order statistic interval) is not a double.
    both_ends = function(v) {
      rep_len(c(-1.7e308, -1.6e308, 1.5e308, 1.7e308, 1.6e308), length(v))
    },
    constant = function(v) rep(5, length(v))
  )
  checks <- list(
    seq_aipe_sm = function(f) seq_aipe_sm(f(x), width = 0.6),
    seq_aipe_cv = function(f) seq_aipe_cv(f(pos), width = 0.1),
    seq_aipe_gini = function(f) seq_aipe_gini(f(pos), width = 0.1),
    seq_aipe_smd = function(f) seq_aipe_smd(f(x), f(x + z), width = 0.9),
    seq_aipe_r = function(f) seq_aipe_r(f(x), f(y), width = 0.4),
    seq_aipe_r_kendall = function(f) {
      seq_aipe_r(f(x), f(y), width = 0.4, method = "kendall")
    },
    seq_aipe_reg_coef = function(f) seq_aipe_reg_coef(f(x), f(y), width = 0.4),
    seq_aipe_R2 = function(f) {
      seq_aipe_R2(cbind(f(x), f(z + 3)), f(y), width = 0.4)
    },
    seq_aipe_c = function(f) {
      seq_aipe_c(list(f(x), f(y), f(x + z)), c(1, -0.5, -0.5), width = 8)
    },
    seq_mr_cv = function(f) seq_mr_cv(f(pos), A = 4e5, c = 75),
    seq_mr_smd = function(f) seq_mr_smd(f(x), f(x + z), A = 1e4, c = 5),
    seq_mr_norm_mean = function(f) seq_mr_norm_mean(f(x), A = 100, c = 1),
    seq_mr_exp_mean = function(f) seq_mr_exp_mean(f(pos), A = 100, c = 0.01),
    seq_mr_nexp_mean = function(f) {
      seq_mr_nexp_mean(f(pos), A = 100, c = 0.01, k = 1, t = 2)
    },
    seq_mr_dfr_mean_arb = function(f) {
      seq_mr_dfr_mean_arb(f(x), A = 100, c = 1)
    },
    seq_br_norm_mean = function(f) seq_br_norm_mean(f(x), A = 100, k = 2, w = 1),
    seq_blci_dfr_median = function(f) seq_blci_dfr_median(f(x), d = 1),
    seq_fwd_norm_mean = function(f) seq_fwd_norm_mean(f(x), d = 0.6),
    seq_fwd_dfr_mean = function(f) seq_fwd_dfr_mean(f(x), d = 0.6),
    seq_fwd_nexp_mean = function(f) seq_fwd_nexp_mean(f(pos), d = 0.1),
    seq_budget = function(f) {
      seq_budget(list(f(x), f(y), f(x + z)), c(1, -0.5, -0.5),
                 cost_per_group = c(38, 15, 10), budget = 3000)
    },
    ts_fwd_norm_mean = function(f) {
      ts_fwd_norm_mean(f(x)[1:15], f(x)[16:40], d = 0.6)
    },
    ts_fwd_dfr_mean = function(f) {
      ts_fwd_dfr_mean(f(x)[1:30], f(x)[31:40], d = 0.6)
    },
    ts_fwd_nexp_mean = function(f) {
      ts_fwd_nexp_mean(f(pos)[1:15], f(pos)[16:40], d = 1)
    },
    ts_br_norm_mean = function(f) {
      ts_br_norm_mean(f(x)[1:15], f(x)[16:40], A = 100, k = 2, w = 1)
    },
    ts_br_exp_mean = function(f) {
      ts_br_exp_mean(f(pos)[1:20], f(pos)[21:40], A = 100, w = 1)
    }
  )
  returned <- 0L
  for (fn in names(checks)) {
    for (shape in names(shapes)) {
      info <- paste(fn, shape)
      out <- suppressWarnings(
        tryCatch(checks[[fn]](shapes[[shape]]), error = identity)
      )
      if (inherits(out, "error")) {
        # One of the package's own messages: no call attached, and not R's
        # text for a condition that was NA or of the wrong length.
        expect_null(conditionCall(out), info = info)
        expect_false(grepl("missing value where|argument is of length|non-numeric",
                           conditionMessage(out)), info = info)
      } else {
        returned <- returned + 1L
        expect_s3_class(out, "smsd_check")
        expect_false(any(is.infinite(out$value) | is.nan(out$value)),
                     info = info)
        expect_true(attr(out, "is_satisfied") %in% c(TRUE, FALSE),
                    info = info)
      }
    }
  }
  # Every function returns an object on the ordinary data, and most calls
  # at the edges return one too, so the loop above is not a list of errors.
  expect_gte(returned, length(checks) + 100L)
})
