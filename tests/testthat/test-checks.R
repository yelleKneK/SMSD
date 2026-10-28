# The shared argument checks: each must fire with a message that names the
# argument, and must never surface R's own condition errors.

test_that(".check_number accepts a valid value and rejects each bad kind", {
  expect_identical(.check_number(0.3, "width"), 0.3)
  expect_identical(.check_number(2L, "k"), 2)
  for (bad in list(NA_real_, NaN, Inf, -Inf, "a", NULL, c(1, 2), TRUE,
                   matrix(1, 1, 1))) {
    expect_error(.check_number(bad, "width"),
                 "'width' must be a single finite number")
  }
})

test_that(".check_number honors open and closed bounds", {
  expect_error(.check_number(0, "w", lower = 0, lower_open = TRUE),
               "greater than 0")
  expect_identical(.check_number(0, "w", lower = 0), 0)
  expect_error(.check_number(1, "p", lower = 0, upper = 1,
                             lower_open = TRUE, upper_open = TRUE),
               "greater than 0 and less than 1")
  expect_error(.check_number(2, "p", upper = 1), "at most 1")
})

test_that(".check_width, .check_positive, and .check_conf_level fire", {
  expect_error(.check_width(), "'width' must be supplied")
  expect_error(.check_width(0), "'width' must be a single finite number")
  expect_error(.check_width(NA_real_), "got NA")
  expect_error(.check_positive(-1, "A"), "'A' must be a single finite number")
  expect_error(.check_conf_level(1), "'conf_level'")
  expect_error(.check_conf_level(NA_real_), "'conf_level'")
  expect_identical(.check_conf_level(0.9), 0.9)
})

test_that(".check_count requires a whole number at or above the minimum", {
  expect_identical(.check_count(10, "G"), 10L)
  expect_error(.check_count(2.5, "G"), "'G' must be a single whole number")
  expect_error(.check_count(3, "m_0", min = 4L), "at least 4")
  expect_error(.check_count(NA_real_, "G"), "'G'")
  expect_error(.check_count(1e12, "G"), "'G'")
})

test_that(".check_flag requires TRUE or FALSE", {
  expect_true(.check_flag(TRUE, "na.rm"))
  expect_error(.check_flag(NA, "na.rm"), "'na.rm' must be TRUE or FALSE")
  expect_error(.check_flag("yes", "na.rm"), "'na.rm' must be TRUE or FALSE")
})

test_that(".check_data cleans valid data and returns doubles", {
  expect_identical(.check_data(1:5, "x"), as.numeric(1:5))
  expect_identical(.check_data(c(1, NA, 3), "x", na.rm = TRUE), c(1, 3))
})

test_that(".check_data rejects each kind of unusable data with advice", {
  expect_error(.check_data(matrix(1:8, ncol = 2), "x"),
               "must be a numeric vector, not a matrix")
  expect_error(.check_data(data.frame(a = 1:4), "x"), "not a data frame")
  expect_error(.check_data(letters, "x"), "must be a numeric vector")
  expect_error(.check_data(c(1, NA, 3), "x"), "Remove or impute")
  expect_error(.check_data(c(1, NA, 3), "x", na.rm = FALSE),
               "set na.rm = TRUE")
  expect_error(.check_data(c(1, Inf, 3), "x"), "infinite values")
  expect_error(.check_data(c(1, 2), "x", min_n = 4L),
               "needs at least 4 observations.*it has 2")
  expect_error(.check_data(c(1, NA, NA, 2), "x", min_n = 4L, na.rm = TRUE),
               "needs at least 4 observations.*it has 2")
  expect_error(.check_data(c(1, 0, 2), "x", positive = TRUE),
               "only positive values")
  expect_error(.check_data(c(1, -1, 2), "x", nonnegative = TRUE),
               "no negative values")
  expect_identical(.check_data(c(1, 0, 2), "x", nonnegative = TRUE),
                   c(1, 0, 2))
})

test_that(".as_sample_size rounds up and refuses unrepresentable values", {
  expect_identical(.as_sample_size(10.2), 11L)
  expect_identical(.as_sample_size(10), 10L)
  expect_error(.as_sample_size(1e24), "too large to represent")
  expect_error(.as_sample_size(Inf), "too large to represent")
  expect_error(.as_sample_size(NaN), "too large to represent")
})

test_that(".var_unusable is TRUE for zero, negative, and non-finite estimates", {
  expect_true(.var_unusable(-0.2))
  expect_true(.var_unusable(0))
  expect_true(.var_unusable(NaN))
  expect_true(.var_unusable(NA_real_))
  expect_true(.var_unusable(Inf))
  # A positive estimate is usable however small: the bound, not this
  # predicate, deals with it.
  expect_false(.var_unusable(1e-12))
  expect_false(.var_unusable(0.5))
  expect_identical(.var_unusable(c(0.5, 0, -1, NA, 1e-12)),
                   c(FALSE, TRUE, TRUE, TRUE, FALSE))
})

test_that(".var_floor applies the bound n^-3 and can return the raw value", {
  expect_identical(.var_floor(0.5, 10), 0.5)
  expect_identical(.var_floor(1e-12, 10), 10^-3)
  expect_identical(.var_floor(-0.2, 10), 10^-3)
  expect_identical(.var_floor(NaN, 10), 10^-3)
  expect_identical(.var_floor(-0.2, 10, floor = FALSE), -0.2)
  expect_identical(.var_floor(1e-12, 10, floor = FALSE), 1e-12)
  expect_true(is.nan(.var_floor(NaN, 10, floor = FALSE)))
})

test_that(".check_choice names the argument and matches partially", {
  ch <- c("pearson", "kendall", "spearman")
  expect_identical(.check_choice(ch, ch, "method"), "pearson")
  expect_identical(.check_choice("kend", ch, "method"), "kendall")
  expect_error(.check_choice("tau", ch, "method"),
               "'method' must be one of \"pearson\", \"kendall\", \"spearman\"")
  expect_error(.check_choice(NA_character_, ch, "method"), "'method'")
  expect_error(.check_choice(1, ch, "method"), "'method'")
  expect_error(.check_choice(c("pearson", "kendall"), ch, "method"),
               "'method'")
})

test_that(".check_seed accepts NULL or a whole number set.seed() accepts", {
  expect_null(.check_seed(NULL))
  expect_identical(.check_seed(113), 113L)
  expect_identical(.check_seed(-5), -5L)
  for (bad in list(1.5, NA_real_, Inf, 1e10, "a", c(1, 2), matrix(1, 1, 1))) {
    expect_error(.check_seed(bad), "'seed' must be NULL or a single whole")
  }
})

test_that(".check_vector checks type, finiteness, length, and sign", {
  expect_identical(.check_vector(1:3, "w", "one weight per group"),
                   c(1, 2, 3))
  expect_error(.check_vector("a", "w", "one weight per group"),
               "'w' must be a numeric vector \\(one weight per group\\)")
  expect_error(.check_vector(matrix(1:4, 2), "w", "weights"), "'w'")
  expect_error(.check_vector(numeric(0), "w", "weights"), "length 0")
  expect_error(.check_vector(c(1, NA, 3), "w", "weights"),
               "entry 2 is NA")
  expect_error(.check_vector(c(1, NA, Inf), "w", "weights"),
               "2 entries in all")
  expect_error(.check_vector(c(1, 2), "w", "weights", len = 3,
                             len_what = "one per group"),
               "must have 3 entries \\(one per group\\); it has 2")
  expect_error(.check_vector(c(1, 0), "cost", "costs", positive = TRUE),
               "entry 2 is 0")
})

test_that(".is_constant uses one relative tolerance of 1e-10", {
  expect_true(.is_constant(rep(5, 10)))
  expect_true(.is_constant(c(0, 0, 0)))
  # Rounding noise on equal values: 0.3 computed as differences near 1e5.
  gains <- (1e5 + 0.3) - rep(1e5, 12) + c(0, 1e-12)
  expect_true(.is_constant(gains))
  expect_false(.is_constant(c(1, 1 + 1e-8)))
  expect_false(.is_constant(c(-1, 1)))
  expect_identical(.smsd_constant_tol, 1e-10)
})

test_that(".as_sample_size ignores representation error but not real excess", {
  expect_identical(.as_sample_size(56.00000000000001), 56L)
  expect_identical(.as_sample_size(56.001), 57L)
  expect_identical(.as_sample_size(1e8 + 0.5), 100000001L)
  expect_identical(.as_sample_size(c(10.2, 3)), c(11L, 3L))
  expect_error(.as_sample_size(c(10, Inf)), "too large to represent")
  expect_error(.as_sample_size(1e24, "pilot sample size",
                               advice = "Lower 'A'."),
               "implied pilot sample size .* Lower 'A'\\.")
})

test_that("a one dimensional numeric array is taken as the vector it is", {
  # tapply() returns one (a mean per participant, a mean per group), and
  # mean(), sd(), and t.test() take it as a vector.
  means <- tapply(c(1, 2, 3, 4), c(1, 1, 2, 2), mean)
  expect_identical(dim(means), 2L)
  expect_identical(.check_data(means, "x"), c(1.5, 3.5))
  expect_identical(.check_data(array(c(2L, 4L, 9L), 3), "x"), c(2, 4, 9))
  expect_identical(.check_vector(means, "true_means", "one mean per group"),
                   c(1.5, 3.5))
  expect_identical(.check_vector(means, "true_means", "one mean per group",
                                 len = 2, positive = TRUE), c(1.5, 3.5))
  # The values are still checked.
  expect_error(.check_data(array(c(1, NA, 3), 3), "x"),
               "'x' contains missing values")
  expect_error(.check_data(array(c(1, Inf, 3), 3), "x"), "infinite values")
  expect_error(.check_data(array(numeric(0), 0), "x", min_n = 4L),
               "needs at least 4 observations; it has 0")
  expect_error(.check_vector(array(c(1, 0), 2), "cost", "costs",
                             positive = TRUE), "entry 2 is 0")
  expect_error(.check_vector(means, "w", "weights", len = 3),
               "must have 3 entries; it has 2")
})

test_that("a matrix or higher array is refused, with advice that can work", {
  # A matrix: the advice names a column, and x[, 1] is a vector.
  expect_error(
    .check_data(matrix(1:8, ncol = 2), "x"),
    paste0("^'x' must be a numeric vector, not a matrix or array with ",
           "dimensions 4 x 2\\. Pass a single column, for example ",
           "x\\[, 1\\] or x\\[\\[\"column\"\\]\\]\\.$")
  )
  expect_error(.check_data(matrix(1:8, ncol = 1), "x"),
               "dimensions 8 x 1\\. Pass a single column")
  expect_error(.check_data(data.frame(a = 1:4), "x"),
               "not a data frame with dimensions 4 x 1\\. Pass a single column")
  # Three dimensions: x[, 1] is an error there, so it is not advised.
  expect_identical(
    tryCatch(.check_data(array(1, c(2, 2, 2)), "x"),
             error = conditionMessage),
    paste0("'x' must be a numeric vector, not a matrix or array with ",
           "dimensions 2 x 2 x 2.")
  )
  # Nothing is unclassed into acceptance. A one dimensional table() is
  # numeric, but it holds the counts of the values (2, 1, 3 here), not the
  # six values; a by() result is classed too. A Date is not numeric, with or
  # without a dim attribute, and neither is text.
  counts <- table(c(3, 3, 5, 7, 7, 7))
  expect_true(is.numeric(counts) && length(dim(counts)) == 1L)
  expect_identical(
    tryCatch(.check_data(counts, "x"), error = conditionMessage),
    "'x' must be a numeric vector, not a matrix or array with dimensions 3."
  )
  expect_error(.check_vector(counts, "true_means", "one mean per group"),
               "got a matrix or array with dimensions 3")
  by_means <- by(c(1, 2, 3, 4), c(1, 1, 2, 2), mean)
  expect_error(.check_data(by_means, "x"),
               "not a matrix or array with dimensions 2\\.$")
  dates <- structure(as.Date("2026-01-01") + 1:4, dim = 4L)
  expect_identical(
    tryCatch(.check_data(dates, "x"), error = conditionMessage),
    "'x' must be a numeric vector, not a matrix or array with dimensions 4."
  )
  expect_error(.check_data(array(letters[1:3], 3), "x"),
               "not a matrix or array with dimensions 3\\.$")
  expect_error(.check_vector(dates, "w", "weights"),
               "'w' must be a numeric vector \\(weights\\)")
  expect_error(.check_vector(matrix(1:4, 2), "w", "weights"),
               "got a matrix or array with dimensions 2 x 2")
})

test_that("the check functions take tapply() output as they take a vector", {
  # Reaction times, four per participant; the analysis is of the 30 per
  # participant means, which tapply() returns as a one dimensional array.
  withr::local_seed(113)
  d <- data.frame(id = rep(1:30, each = 4), rt = rgamma(120, 16, 20))
  person_mean <- tapply(d$rt, d$id, mean)
  expect_identical(dim(person_mean), 30L)
  plain <- as.vector(person_mean)
  expect_identical(seq_aipe_cv(person_mean, width = 0.1),
                   seq_aipe_cv(plain, width = 0.1))
  expect_identical(seq_fwd_norm_mean(person_mean, d = 0.1),
                   seq_fwd_norm_mean(plain, d = 0.1))
  expect_identical(seq_mr_cv(person_mean, A = 400, c = 1),
                   seq_mr_cv(plain, A = 400, c = 1))
  half <- tapply(d$rt, list(d$id, rep(1:2, 60)), mean)
  expect_identical(seq_aipe_smd(person_mean, array(half[, 2], 30),
                                width = 1),
                   seq_aipe_smd(plain, half[, 2], width = 1))
  expect_identical(
    ts_fwd_norm_mean(person_mean[1:12], person_mean[13:30], d = 0.2),
    ts_fwd_norm_mean(plain[1:12], plain[13:30], d = 0.2)
  )
  # Group means and standard deviations from tapply() as population values.
  score <- c(48, 52, 50, 46, 54, 49, 51, 50, 47, 53, 50, 50)
  group <- rep(1:3, each = 4)
  args <- list(contrast_weights = c(1, -0.5, -0.5),
               cost_per_group = c(38, 15, 10), budget = 600, G = 4,
               seed = 113)
  expect_identical(
    do.call(seq_budget_sensitivity,
            c(list(true_means = tapply(score, group, mean),
                   true_sds = tapply(score, group, sd)), args)),
    do.call(seq_budget_sensitivity,
            c(list(true_means = as.vector(tapply(score, group, mean)),
                   true_sds = as.vector(tapply(score, group, sd))), args))
  )
})

# A class of number for which every comparison and every arithmetic
# operation is NA, standing in for a 64 bit integer of package bit64 (what
# data.table::fread() returns for large whole numbers), whose comparison with
# Inf is NA and whose arithmetic truncates to whole numbers. bit64 is not a
# dependency of the package. The Ops method is registered for the calling
# test only (an Ops method is found through the registry, not the calling
# environment) and removed when that test ends. Returns the constructor.
local_incomparable <- function(env = parent.frame()) {
  registerS3method("Ops", "smsd_test_incomparable", function(e1, e2) NA)
  withr::defer(rm(list = "Ops.smsd_test_incomparable",
                  envir = get(".__S3MethodsTable__.", envir = baseenv())),
               envir = env)
  function(x) structure(x, class = "smsd_test_incomparable")
}

test_that("a number with a class is refused by the scalar checks", {
  odd <- local_incomparable()
  expect_true(is.na(odd(1) <= Inf))
  expect_true(is.na(odd(4) == round(odd(4))))
  expect_true(is.numeric(odd(1)) && is.finite(odd(1)))
  plain <- function(name, what = "a plain number") {
    paste0("^'", name, "' must be ", what, ", not an object of class ",
           "'smsd_test_incomparable'\\. Pass as\\.numeric\\(", name,
           "\\)\\.$")
  }
  expect_error(.check_number(odd(1), "width"), plain("width"))
  expect_error(.check_positive(odd(0.5), "A"), plain("A"))
  expect_error(.check_width(odd(0.5)), plain("width"))
  expect_error(.check_conf_level(odd(0.9)), plain("conf_level"))
  expect_error(.check_count(odd(4), "m_0", min = 4L),
               plain("m_0", "a plain whole number"))
  expect_error(.check_seed(odd(113)),
               plain("seed", "NULL or a plain whole number"))
  # Also when the value it holds is out of range: the class is named first.
  expect_error(.check_positive(odd(-1), "A"), plain("A"))
  # What is not a number is described as before.
  expect_error(.check_number(factor("a"), "A"),
               "^'A' must be a single finite number; got an object of class")
  expect_error(.check_count(as.Date("2026-01-01"), "G"),
               "^'G' must be a single whole number of at least 1; got ")
  # Vectors are returned as doubles, so a class does not reach the caller.
  expect_identical(.check_vector(odd(c(38, 15, 10)), "cost_per_group",
                                 "costs", positive = TRUE), c(38, 15, 10))
  expect_error(.check_vector(odd(c(1, 0)), "cost", "costs", positive = TRUE),
               "entry 2 is 0")
  expect_identical(.check_data(odd(c(2, 4, 9)), "x"), c(2, 4, 9))
})

test_that("every exported function refuses a number with a class", {
  # Each numeric scalar argument of each exported function in turn. Before
  # this check, such a number either stopped with R's "missing value where
  # TRUE/FALSE needed" or was used in arithmetic as it stood: a bit64 'A' in
  # ts_br_exp_mean() gave a final sample size of 20 against 2756.
  odd <- local_incomparable()
  withr::local_seed(113)
  x <- rgamma(40, 16, 20)
  y <- 0.5 * x + rnorm(40, 0, 0.1)
  g <- list(rnorm(12, 5), rnorm(12, 6), rnorm(12, 7))
  w <- c(1, -0.5, -0.5)
  cst <- c(38, 15, 10)
  bud <- seq_budget(g, w, cst, budget = 3000)
  calls <- list(
    list(decide_seq_budget_equivalence, list(check = bud, margin_lower = -1,
                                             margin_upper = 1,
                                             alpha_level = 0.05)),
    list(decide_seq_budget_noninferiority, list(check = bud, margin = 1,
                                                alpha_level = 0.05)),
    list(power_z, list(mu_0 = 100, mu_1 = 104, sigma = 15, conf_level = 0.95,
                       n = 50)),
    list(power_z, list(mu_0 = 100, mu_1 = 104, sigma = 15, n_min = 10,
                       n_max = 20, n_values = 5)),
    list(power_z, list(mu_0 = 100, mu_1 = 104, sigma = 15, beta = 0.2)),
    list(seq_aipe_c, list(groups = g, contrast_weights = w, width = 2,
                          conf_level = 0.95)),
    list(seq_aipe_c_sensitivity, list(true_means = c(5, 6, 7),
                                      contrast_weights = w, width = 2,
                                      conf_level = 0.95, G = 2, seed = 1,
                                      m_0 = 4, max_n = 1000, batch = 1)),
    list(seq_aipe_cv, list(x = x, width = 0.2, conf_level = 0.95)),
    list(seq_aipe_cv_pilot, list(width = 0.2, conf_level = 0.95, m_0 = 4)),
    list(seq_aipe_cv_sensitivity, list(true_value = 0.2, width = 0.2,
                                       conf_level = 0.95, G = 2, seed = 1,
                                       m_0 = 4, max_n = 1000, batch = 1)),
    list(seq_aipe_gini, list(x = x, width = 0.2, conf_level = 0.95)),
    list(seq_aipe_r, list(x = x, y = y, width = 0.4, conf_level = 0.95)),
    list(seq_aipe_R2_pilot, list(width = 0.3, p = 2, conf_level = 0.95,
                                 m_0 = 4)),
    list(seq_aipe_reg_coef, list(x = x, y = y, width = 0.5,
                                 conf_level = 0.95)),
    list(seq_aipe_sm, list(x = x, width = 0.5, conf_level = 0.95)),
    list(seq_aipe_smd, list(group_1 = x, group_2 = y, width = 1,
                            conf_level = 0.95)),
    list(seq_aipe_smd_sensitivity, list(true_value = 0.5, width = 1,
                                        conf_level = 0.95, G = 2, seed = 1,
                                        m_0 = 4, max_n = 1000, batch = 1)),
    list(seq_blci_dfr_median, list(x = x, d = 0.3, conf_level = 0.95)),
    list(seq_blci_dfr_median_pilot, list(d = 0.3, conf_level = 0.95,
                                         m_0 = 10, gamma = 1)),
    list(seq_br_norm_mean, list(x = x, A = 100, k = 2, w = 1)),
    list(seq_br_norm_mean_pilot, list(A = 100, k = 2, w = 1, m_0 = 10,
                                      gamma = 1)),
    list(seq_budget, list(groups = g, contrast_weights = w,
                          cost_per_group = cst, budget = 3000,
                          conf_level = 0.95)),
    list(seq_budget_optimal, list(contrast_weights = w, sigmas = c(1, 2, 3),
                                  cost_per_group = cst, budget = 3000)),
    list(seq_budget_pilot, list(m_0 = 5)),
    list(seq_budget_sensitivity, list(true_means = c(50, 50, 50),
                                      true_sds = c(9, 6, 6),
                                      contrast_weights = w,
                                      cost_per_group = cst, budget = 1500,
                                      conf_level = 0.95, margin = 1,
                                      alpha_level = 0.05, G = 2, seed = 1,
                                      m_0 = 5)),
    list(seq_fwd_dfr_mean, list(x = x, d = 0.3, conf_level = 0.95)),
    list(seq_fwd_dfr_mean_pilot, list(d = 0.3, conf_level = 0.95, m_0 = 10,
                                      gamma = 1)),
    list(seq_fwd_nexp_mean, list(x = x, d = 0.3, conf_level = 0.95)),
    list(seq_fwd_nexp_mean_pilot, list(d = 0.3, conf_level = 0.95, m_0 = 10,
                                       gamma = 1)),
    list(seq_fwd_norm_mean, list(x = x, d = 0.3, conf_level = 0.95)),
    list(seq_fwd_norm_mean_optimal, list(variance = 30, d = 2,
                                         conf_level = 0.95)),
    list(seq_fwd_norm_mean_pilot, list(d = 0.3, conf_level = 0.95, m_0 = 10,
                                       gamma = 1)),
    list(seq_mr_cv, list(x = x, A = 400, c = 1, gamma = 0.49)),
    list(seq_mr_cv_pilot, list(A = 400, c = 1, m_0 = 4, gamma = 0.49)),
    list(seq_mr_cv_sensitivity, list(true_value = 0.3, A = 20000, c = 75,
                                     G = 2, seed = 1, m_0 = 4, gamma = 0.49,
                                     max_n = 1000, batch = 1)),
    list(seq_mr_dfr_mean_arb, list(x = x, A = 100, c = 1, gamma = 1)),
    list(seq_mr_dfr_mean_arb_pilot, list(A = 100, c = 1, m_0 = 4,
                                         gamma = 1)),
    list(seq_mr_exp_mean, list(x = x, A = 100, c = 1)),
    list(seq_mr_exp_mean_pilot, list(A = 100, c = 1, m_0 = 4, gamma = 1)),
    list(seq_mr_nexp_mean, list(x = x, A = 100, c = 1, k = 1, t = 2)),
    list(seq_mr_nexp_mean_pilot, list(A = 100, c = 1, k = 1, t = 2, m_0 = 4,
                                      gamma = 1)),
    list(seq_mr_norm_mean, list(x = x, A = 100, c = 1)),
    list(seq_mr_norm_mean, list(x = x, c = 1, price = 10, epsilon = 1)),
    list(seq_mr_norm_mean_pilot, list(A = 100, c = 1, m_0 = 4, gamma = 1)),
    list(seq_mr_norm_mean_pilot, list(c = 1, price = 10, epsilon = 1)),
    list(seq_mr_smd, list(group_1 = x, group_2 = y, A = 1000, c = 1,
                          gamma = 0.49)),
    list(seq_mr_smd_pilot, list(A = 1000, c = 1, m_0 = 4, gamma = 0.49)),
    list(seq_mr_smd_sensitivity, list(true_value = 0.5, A = 2000, c = 5,
                                      G = 2, seed = 1, m_0 = 4, gamma = 0.49,
                                      max_n = 1000, batch = 1)),
    list(ts_br_exp_mean, list(pilot_data = x[1:20],
                              additional_data = x[21:40], A = 100, w = 3)),
    list(ts_br_exp_mean_pilot, list(A = 100, w = 3, m_0 = 20, gamma = 1)),
    list(ts_br_norm_mean, list(pilot_data = x[1:15],
                               additional_data = x[16:40], A = 100, k = 2,
                               w = 3)),
    list(ts_br_norm_mean_pilot, list(A = 100, k = 2, w = 3, m_0 = 10,
                                     gamma = 1)),
    list(ts_fwd_dfr_mean, list(pilot_data = x[1:30],
                               additional_data = x[31:40], d = 0.3,
                               conf_level = 0.95)),
    list(ts_fwd_dfr_mean_pilot, list(d = 0.3, conf_level = 0.95, m_0 = 30,
                                     gamma = 1)),
    list(ts_fwd_nexp_mean, list(pilot_data = x[1:15],
                                additional_data = x[16:40], d = 0.3,
                                conf_level = 0.95)),
    list(ts_fwd_norm_mean, list(pilot_data = x[1:15],
                                additional_data = x[16:40], d = 0.3,
                                conf_level = 0.95)),
    list(ts_fwd_norm_mean_pilot, list(d = 0.3, conf_level = 0.95, m_0 = 10,
                                      gamma = 1))
  )
  checked <- 0L
  for (call in calls) {
    args <- call[[2L]]
    scalar <- vapply(args, function(v) is.numeric(v) && length(v) == 1L &&
                       !is.object(v), logical(1))
    for (name in names(args)[scalar]) {
      with_class <- args
      with_class[[name]] <- odd(args[[name]])
      expect_error(
        do.call(call[[1L]], with_class),
        paste0("'", name, "' must be [^.]*plain[^.]*, not an object of ",
               "class 'smsd_test_incomparable'\\. Pass as\\.numeric\\(",
               name, "\\)\\.")
      )
      checked <- checked + 1L
    }
  }
  expect_gte(checked, 200L)
  # The two costs of the standardized mean difference procedures as well.
  expect_error(seq_mr_smd(x, y, A = 1000, c = odd(c(1, 3))),
               "^'c' must be one or two plain numbers, not an object of class")
  expect_error(seq_mr_smd_pilot(A = 1000, c = odd(c(1, 3))),
               "Pass as\\.numeric\\(c\\)\\.$")
  # Vector arguments are converted to doubles: the same result as the plain
  # vector.
  expect_identical(seq_budget(g, odd(w), odd(cst), budget = 3000), bud)
  expect_identical(seq_budget_optimal(odd(w), odd(c(1, 2, 3)), odd(cst),
                                      budget = 3000),
                   seq_budget_optimal(w, c(1, 2, 3), cst, budget = 3000))
  expect_identical(seq_aipe_c(g, odd(w), width = 2),
                   seq_aipe_c(g, w, width = 2))
})
