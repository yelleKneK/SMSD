# Tests for seq_aipe_c() and seq_aipe_c_pilot(): the sequential AIPE
# procedure of Kelley, Darku, and Chattopadhyay (2018) applied to a contrast
# of group means (their Example 6).

row_value <- function(out, term) out$value[out$term == term]

test_that("seq_aipe_c_pilot is max{m_0, ceiling(2 z / width)} (2018, Equation 21)", {
  expect_identical(seq_aipe_c_pilot(width = 0.40), 10L)
  expect_identical(seq_aipe_c_pilot(width = 6), 4L)
  expect_identical(seq_aipe_c_pilot(width = 0.40, conf_level = 0.90), 9L)
  expect_identical(seq_aipe_c_pilot(width = 0.40, m_0 = 20), 20L)
})

test_that("seq_aipe_c_pilot checks its arguments", {
  expect_pilot_checks(seq_aipe_c_pilot)
})

test_that("seq_aipe_c returns the documented rows and attributes", {
  withr::local_seed(113)
  groups <- list(rnorm(30, mean = 0.5), rnorm(30), rnorm(30))
  cw <- c(1, -0.5, -0.5)
  out <- seq_aipe_c(groups, contrast_weights = cw, width = 0.50)
  expect_s3_class(out, c("smsd_check", "data.frame"), exact = TRUE)
  expect_identical(out$term, c("lower_limit", "contrast", "upper_limit",
                               "xi_squared", "criterion_n", "current_n"))
  expect_true(all(is.finite(out$value)))
  expect_identical(attr(out, "method"), "contrast")
  expect_identical(attr(out, "n_per_group"), c(30L, 30L, 30L))
  expect_identical(attr(out, "contrast_weights"), cw)
  expect_identical(attr(out, "variance_unusable"), FALSE)
  expect_identical(attr(out, "notes"), character(0))
  expect_identical(attr(out, "conf_level"), 0.95)
  expect_identical(attr(out, "width"), 0.50)
  # criterion_n and current_n are on the same (per group) scale.
  expect_equal(attr(out, "current_n"), 30)
  expect_equal(row_value(out, "current_n"), 30)
  expect_silent(seq_aipe_c(groups, cw, width = 0.50))
  # Integer weights are stored as doubles.
  expect_identical(
    attr(seq_aipe_c(groups, c(1L, -1L, 0L), width = 0.5), "contrast_weights"),
    c(1, -1, 0))
})

test_that("seq_aipe_c keeps group names on n_per_group and nowhere else", {
  # A named list (what split(y, g) returns) used to pass the name of its
  # first group to current_n and criterion_n, and through them to the row
  # name of summary().
  withr::local_seed(113)
  groups <- list(rnorm(30, mean = 0.5), rnorm(30), rnorm(30))
  cw <- c(1, -0.5, -0.5)
  plain <- seq_aipe_c(groups, cw, width = 0.5)
  named <- seq_aipe_c(list(new = groups[[1]], a = groups[[2]],
                           b = groups[[3]]), cw, width = 0.5)
  expect_identical(attr(named, "n_per_group"), c(new = 30L, a = 30L, b = 30L))
  expect_identical(named$value, plain$value)
  expect_identical(attr(named, "current_n"), 30L)
  expect_identical(attr(named, "current_n"), attr(plain, "current_n"))
  expect_null(names(attr(named, "criterion_n")))
  expect_identical(attr(named, "criterion_n"), attr(plain, "criterion_n"))
  expect_identical(attr(named, "is_satisfied"),
                   attr(named, "current_n") >= attr(named, "criterion_n"))
  expect_identical(rownames(summary(named)), "1")
  expect_identical(summary(named), summary(plain))
  # The same in the "not yet usable" state.
  k <- rep(3, 10)
  flat <- seq_aipe_c(list(new = k, old = k + 1), c(1, -1), width = 0.5)
  expect_null(names(attr(flat, "current_n")))
  expect_null(names(attr(flat, "criterion_n")))
  expect_identical(attr(flat, "is_satisfied"),
                   attr(flat, "current_n") >= attr(flat, "criterion_n"))
})

test_that("seq_aipe_c computes the contrast, its variance, the interval, and the criterion", {
  # Hand example: groups 1..4, 2 4 6 8, and 3 3 5 9 with weights 1, -1/2,
  # -1/2. Means 2.5, 5, 5; variances 5/3, 20/3, 8. Contrast
  # 2.5 - 2.5 - 2.5 = -2.5; xi^2 = 5/3 + 20/12 + 8/4 = 16/3.
  groups <- list(c(1, 2, 3, 4), c(2, 4, 6, 8), c(3, 3, 5, 9))
  cw <- c(1, -0.5, -0.5)
  for (cl in c(0.90, 0.95)) {
    out <- seq_aipe_c(groups, cw, width = 4, conf_level = cl)
    z <- qnorm(1 - (1 - cl) / 2)
    expect_equal(row_value(out, "contrast"), -2.5)
    expect_equal(row_value(out, "xi_squared"), 16 / 3)
    expect_equal(row_value(out, "lower_limit"), -2.5 - z * sqrt(16 / 3 / 4))
    expect_equal(row_value(out, "upper_limit"), -2.5 + z * sqrt(16 / 3 / 4))
    expect_equal(row_value(out, "criterion_n"),
                 ceiling((2 * z / 4)^2 * (16 / 3 + 1 / 4)))
    expect_identical(attr(out, "is_satisfied"),
                     4 >= row_value(out, "criterion_n"))
  }
  # Weights need not sum to zero or to one: a single group mean and a sum.
  expect_equal(row_value(seq_aipe_c(groups, c(1, 0, 0), width = 4),
                         "contrast"), 2.5)
  expect_equal(row_value(seq_aipe_c(groups, c(1, 1, 1), width = 4),
                         "xi_squared"), 5 / 3 + 20 / 3 + 8)
})

test_that("a positive variance estimate below n^-3 is used through the bound, with a note about units", {
  # The contrast has units; the bound n^-3 has none. Scores in meters with
  # a spread of a centimeter: each group has variance 0.00025 (the variance
  # of 1, ..., 5 is 2.5, times 0.01^2), so with weights 1, -1/2, -1/2 the
  # variance estimate is 1.5 * 0.00025 = 0.000375, positive but below
  # 5^-3 = 0.008. The same scores in centimeters give 3.75. The bound
  # stands in for the estimate, the rule is evaluated, and the interval is
  # built; a note says so. (An earlier version reported such a check as
  # "not yet usable", with NA limits.)
  meters <- list(c(0.48, 0.49, 0.50, 0.51, 0.52),
                 c(0.58, 0.59, 0.60, 0.61, 0.62),
                 c(0.38, 0.39, 0.40, 0.41, 0.42))
  cw <- c(1, -0.5, -0.5)
  z <- qnorm(0.975)
  expect_equal(.var_c(meters, cw, floor = FALSE), 0.000375)
  out <- seq_aipe_c(meters, cw, width = 0.02)
  expect_false(attr(out, "variance_unusable"))
  expect_identical(row_value(out, "xi_squared"), 5^(-3))
  expect_equal(row_value(out, "contrast"), 0.50 - 0.30 - 0.20)
  expect_equal(row_value(out, "upper_limit") - row_value(out, "lower_limit"),
               2 * z * sqrt(5^(-3) / 5))
  # (2 z / 0.02)^2 (0.008 + 1/5) = 7990.2
  expect_equal(row_value(out, "criterion_n"), 7991)
  expect_false(attr(out, "is_satisfied"))
  note <- attr(out, "notes")
  expect_length(note, 1L)
  expect_match(note, "positive but below its lower bound")
  expect_match(note, "the bound was used in its place")
  expect_match(note, "Rescale the outcome")
  expect_match(note, "seq_aipe_c")
  expect_false(grepl("not yet usable", note))

  # The note informs; it does not change the decision. With a target of 2
  # meters the rule with the bound is met at once,
  # (2 z / 2)^2 (0.008 + 1/5) = 0.8, and the note is still there.
  wide <- seq_aipe_c(meters, cw, width = 2)
  expect_equal(row_value(wide, "criterion_n"), 1)
  expect_true(attr(wide, "is_satisfied"))
  expect_match(attr(wide, "notes"), "the bound was used in its place")

  centimeters <- lapply(meters, function(g) 100 * g)
  rescaled <- seq_aipe_c(centimeters, cw, width = 2)
  expect_false(attr(rescaled, "variance_unusable"))
  expect_equal(row_value(rescaled, "xi_squared"), 3.75)
  expect_identical(attr(rescaled, "notes"), character(0))

  # Groups with no variation at all give no estimate of the variance: "not
  # yet usable", with the standard note and one that names the groups.
  k <- rep(3, 5)
  expect_true(is.nan(.var_c(list(k, k + 2, k - 1), cw, floor = FALSE)))
  flat <- seq_aipe_c(list(k, k + 2, k - 1), cw, width = 0.4)
  expect_true(attr(flat, "variance_unusable"))
  expect_false(attr(flat, "is_satisfied"))
  expect_true(is.na(row_value(flat, "lower_limit")))
  expect_length(attr(flat, "notes"), 2L)
  expect_identical(attr(flat, "notes")[1L], .note_variance_unusable())
  expect_match(attr(flat, "notes")[2L],
               "^The observations in groups 1, 2, and 3 are all equal")
})

test_that("a group with no variation makes the variance of the contrast not yet usable", {
  # Ratings on a 1 to 5 scale from three programs of 8; every rating in the
  # second program is 3. The variance of that program has not been
  # estimated: summing the other two terms would leave it out of the
  # variance of the contrast, and at this width the rule would be reported
  # as met with an interval built on two programs. (seq_budget() reads the
  # same data the same way.)
  ratings <- list(new = c(4, 5, 3, 4, 4, 5, 3, 4),
                  existing_a = rep(3, 8),
                  existing_b = c(2, 3, 4, 3, 2, 3, 4, 3))
  cw <- c(1, -0.5, -0.5)
  # With the constant program left out: variances 4/7 and 4/7, so
  # xi^2 = 4/7 + 0.25 * 4/7 = 5/7, and the criterion at width 3 would be
  # ceiling((2 z / 3)^2 (5/7 + 1/8)) = 2, far below 8.
  expect_equal(stats::var(ratings$new), 4 / 7)
  expect_equal(stats::var(ratings$existing_b), 4 / 7)
  expect_true(aipe_rule_met(8, 5 / 7, 3))

  out <- seq_aipe_c(ratings, cw, width = 3)
  expect_false(attr(out, "is_satisfied"))
  expect_true(attr(out, "variance_unusable"))
  expect_true(is.na(row_value(out, "lower_limit")))
  expect_true(is.na(row_value(out, "upper_limit")))
  # The contrast itself is estimated: 4 - 0.5 * 3 - 0.5 * 3 = 1.
  expect_equal(row_value(out, "contrast"), 1)
  expect_identical(row_value(out, "xi_squared"), 8^(-3))
  expect_identical(attr(out, "criterion_n"), 9)
  expect_identical(attr(out, "notes")[1L], .note_variance_unusable())
  expect_identical(
    attr(out, "notes")[2L],
    paste("The observations in group 2 are all equal, or differ too little",
          "relative to their size for a variance to be computed, so the",
          "variance of that group has not been estimated. If a group shows",
          "no variation well past the pilot sample, reconsider the measure",
          "(a ceiling or a floor).")
  )
  expect_true(attr(seq_budget(ratings, cw, c(1, 1, 1), budget = 1000),
                   "variance_unusable"))
  # Two such groups are both named.
  two <- seq_aipe_c(list(ratings$new, rep(3, 8), rep(2, 8)), cw, width = 3)
  expect_match(attr(two, "notes")[2L],
               paste("^The observations in groups 2 and 3 are all equal.*so",
                     "the variance of those groups has not been estimated"))
  # A constant group whose weight is zero takes no part in the contrast:
  # the rule is evaluated on the groups that do, as before. xi^2 = 8/7 and
  # the criterion is ceiling((2 z / 3)^2 (8/7 + 1/8)) = 3.
  pair <- seq_aipe_c(ratings, c(1, 0, -1), width = 3)
  expect_false(attr(pair, "variance_unusable"))
  expect_equal(row_value(pair, "xi_squared"), 8 / 7)
  expect_identical(attr(pair, "criterion_n"), 3)
  expect_true(attr(pair, "is_satisfied"))
  expect_identical(attr(pair, "notes"), character(0))
  # One more rating that differs ends it: every program then varies, and
  # the estimate is the weighted sum of the three sample variances.
  ratings$new <- c(ratings$new, 4)
  ratings$existing_a <- c(ratings$existing_a, 4)
  ratings$existing_b <- c(ratings$existing_b, 3)
  later <- seq_aipe_c(ratings, cw, width = 3)
  expect_false(attr(later, "variance_unusable"))
  expect_equal(row_value(later, "xi_squared"),
               sum(cw^2 * vapply(ratings, stats::var, numeric(1))))
})

test_that("seq_aipe_c refuses weights whose names disagree with the groups", {
  # split() returns the groups in alphabetical order; weights typed in the
  # order of the design, with names, are not in that order. Matched by
  # position they would give another contrast without a word: here
  # 1 * 3 - 0.5 * 7 - 0.5 * 4 = -2.5 in place of 7 - 0.5 * 4 - 0.5 * 3 = 3.5.
  score <- c(6, 7, 8, 7, 3, 4, 5, 4, 2, 3, 4, 3)
  format <- rep(c("in_person", "video", "hybrid"), each = 4)
  groups <- split(score, format)
  expect_identical(names(groups), c("hybrid", "in_person", "video"))
  w <- c(in_person = 1, video = -0.5, hybrid = -0.5)
  expect_error(
    seq_aipe_c(groups, w, width = 8),
    paste0("'contrast_weights' is named \\(in_person, video, hybrid\\) but ",
           "'groups' is named \\(hybrid, in_person, video\\)\\. Entries are ",
           "matched to groups by position, not by name")
  )
  # In the order of the groups the names are accepted and the contrast is
  # the intended one.
  expect_equal(row_value(seq_aipe_c(groups, w[names(groups)], width = 8),
                         "contrast"), 3.5)
  expect_equal(row_value(seq_aipe_c(groups[names(w)], w, width = 8),
                         "contrast"), 3.5)
  # Names on one argument only, or on neither, are matched by position, as
  # documented: hybrid, in_person, video get 1, -0.5, -0.5.
  by_position <- 3 - 0.5 * 7 - 0.5 * 4
  expect_equal(row_value(seq_aipe_c(groups, unname(w), width = 8),
                         "contrast"), by_position)
  expect_equal(row_value(seq_aipe_c(unname(groups), w, width = 8),
                         "contrast"), by_position)
  expect_equal(row_value(seq_aipe_c(unname(groups), unname(w), width = 8),
                         "contrast"), by_position)
  # The weights are checked before their names are compared.
  expect_error(seq_aipe_c(groups, width = 8),
               "'contrast_weights' must be supplied")
  expect_error(seq_aipe_c(groups, c(a = 1, b = -1), width = 8),
               "'contrast_weights' must have 3 entries")
})

test_that("seq_aipe_c rejects unbalanced groups", {
  withr::local_seed(113)
  groups <- list(rnorm(40), rnorm(20), rnorm(25))
  expect_error(
    seq_aipe_c(groups, contrast_weights = c(1, -1, 0), width = 0.50),
    "same per group sample size.*40, 20, 25"
  )
})

test_that("seq_aipe_c satisfied implies delivered width within the target", {
  # Regression for a removed unbalanced branch, which could report
  # is_satisfied = TRUE while delivering an interval up to sqrt(K) times
  # wider than requested. Under the balanced aggregate criterion,
  # n >= (2z/w)^2 (xi^2 + 1/n) implies 2 z sqrt(xi^2 / n) <= w.
  withr::local_seed(113)
  scale_to <- function(x, m, s2) (x - mean(x)) / stats::sd(x) * sqrt(s2) + m
  w <- 0.5
  z <- stats::qnorm(0.975)
  n <- 1000
  cw <- c(1, -2, 1)
  # Variances at 98% of the PER GROUP budget: a per group rule would say
  # stop here even though the aggregate interval is ~sqrt(3) too wide.
  s2 <- 0.98 * (n * (w / (2 * z))^2 - 1 / n) / cw^2
  groups <- lapply(1:3, function(k) scale_to(rnorm(n), k, s2[k]))
  out <- seq_aipe_c(groups, cw, width = w)
  expect_false(attr(out, "is_satisfied"))

  # And whenever the rule IS satisfied, the delivered width obeys the
  # guarantee.
  withr::local_seed(113)
  groups_ok <- lapply(c(0.5, 0, 0), function(m) rnorm(400, mean = m, sd = 1))
  out_ok <- seq_aipe_c(groups_ok, c(1, -0.5, -0.5), width = 0.50)
  expect_true(attr(out_ok, "is_satisfied"))
  achieved <- row_value(out_ok, "upper_limit") - row_value(out_ok, "lower_limit")
  expect_lte(achieved, 0.50)
})

test_that("seq_aipe_c checks 'groups'", {
  withr::local_seed(113)
  x <- rnorm(30)
  y <- rnorm(30)
  expect_error(seq_aipe_c(contrast_weights = c(1, -1), width = 0.5),
               "'groups' must be supplied")
  # A matrix or a data frame is not accepted: a data frame in long format
  # (outcome and numeric group code) would otherwise be read as two groups.
  expect_error(seq_aipe_c(cbind(x, y), c(1, -1), width = 0.5),
               "'groups' must be a list.*matrix")
  expect_error(seq_aipe_c(data.frame(x, y), c(1, -1), width = 0.5),
               "'groups' must be a list.*data frame")
  expect_error(seq_aipe_c(x, c(1, -1), width = 0.5),
               "'groups' must be a list")
  expect_error(seq_aipe_c(list(x), 1, width = 0.5), "at least 2 groups")
  expect_error(seq_aipe_c(list(x, letters), c(1, -1), width = 0.5),
               "'groups\\[\\[2\\]\\]' must be a numeric vector")
  expect_error(seq_aipe_c(list(matrix(x, ncol = 2), y), c(1, -1), width = 0.5),
               "'groups\\[\\[1\\]\\]'.*not a matrix")
  expect_error(seq_aipe_c(list(x, replace(y, 3, NA)), c(1, -1), width = 0.5),
               "'groups\\[\\[2\\]\\]'.*missing values")
  expect_error(seq_aipe_c(list(x, replace(y, 3, Inf)), c(1, -1), width = 0.5),
               "'groups\\[\\[2\\]\\]'.*infinite")
  expect_error(seq_aipe_c(list(x, y[1:3]), c(1, -1), width = 0.5),
               "'groups\\[\\[2\\]\\]'.*at least 4")
})

test_that("seq_aipe_c checks 'contrast_weights'", {
  withr::local_seed(113)
  g <- list(rnorm(30), rnorm(30), rnorm(30))
  expect_error(seq_aipe_c(g, width = 0.5),
               "'contrast_weights' must be supplied")
  # The messages are those of the shared vector check.
  expect_error(seq_aipe_c(g, c(1, -0.5), width = 0.5),
               paste0("'contrast_weights' must have 3 entries (one per ",
                      "group in 'groups'); it has 2."), fixed = TRUE)
  # Previously an NA weight reached the stopping flag, and all zero weights
  # gave a zero width interval with the rule "met".
  expect_error(seq_aipe_c(g, c(1, NA, -1), width = 0.5),
               "'contrast_weights' must contain only finite.*entry 2 is NA")
  expect_error(seq_aipe_c(g, c(1, Inf, -1), width = 0.5),
               "'contrast_weights' must contain only finite.*entry 2 is Inf")
  expect_error(seq_aipe_c(g, c(0, 0, 0), width = 0.5), "not all be zero")
  expect_error(seq_aipe_c(g, c("a", "b", "c"), width = 0.5),
               "'contrast_weights' must be a numeric vector")
  expect_error(seq_aipe_c(g, matrix(c(1, 0, -1), 1), width = 0.5),
               "'contrast_weights' must be a numeric vector")
  expect_error(seq_aipe_c(g, numeric(0), width = 0.5),
               "'contrast_weights' must be a numeric vector.*length 0")
  expect_error(seq_aipe_c(g, NULL, width = 0.5),
               "'contrast_weights' must be a numeric vector.*NULL")
})

test_that("seq_aipe_c checks 'width' and 'conf_level'", {
  withr::local_seed(113)
  g <- list(rnorm(30), rnorm(30))
  expect_scalar_checks(function(...) seq_aipe_c(g, c(1, -1), ...))
})

test_that("operating characteristics: the contrast procedure reaches its target", {
  # No published table exists for the contrast, so the targets are the
  # theoretical sample size and the nominal coverage. Three normal groups
  # with unit variance and weights 1, -1/2, -1/2: xi^2 = 1.5, so at width
  # 0.5 the theoretical per group sample size is
  # ceiling((2 * 1.96 / 0.5)^2 * 1.5) = 93. In 12,000 replications of this
  # condition (200 other seeds of this test) the mean final per group
  # sample size was 92.6 and coverage 0.946; an earlier run of 3,000 gave
  # 92.4 and 0.949. The mean of the 60 replications here had a standard
  # deviation of at most 1.38 over 200 and 300 other seeds, so its window
  # is 4 * 1.38 + 0.6 = 6.2. The standard error of the coverage is
  # sqrt(0.946 * 0.054 / 60) = 0.0292, so its window is
  # 4 * 0.0292 + 0.005 = 0.122.
  withr::local_seed(113)
  w <- 0.5
  cw <- c(1, -0.5, -0.5)
  z <- qnorm(0.975)
  expect_equal(ceiling((2 * z / w)^2 * 1.5), 93)
  take <- function(d, n) lapply(d, function(v) v[1:n])
  stat <- function(d) {
    c(sum(cw * vapply(d, mean, numeric(1))), .var_c(d, cw, floor = FALSE))
  }
  check <- function(d) seq_aipe_c(d, cw, width = w)
  runs <- replicate(60, {
    d <- list(rnorm(400, 0.5), rnorm(400), rnorm(400))
    aipe_run(d, take, stat, check, width = w, pilot = seq_aipe_c_pilot(w))
  })
  s <- aipe_summary(runs, truth = 0.5)
  expect_lt(abs(s[["mean_n"]] - 93), 6.2)
  expect_lt(abs(s[["coverage"]] - 0.95), 0.122)
  expect_lte(s[["max_width"]], w)
})
