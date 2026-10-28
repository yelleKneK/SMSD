# seq_path(): the stopping path of a sequential procedure for data in hand,
# and the methods of the smsd_path object it returns.
#
# The independent value throughout is a direct call of the check function on
# the first n observations, read from the term rows and the attributes of
# the check object (not through seq_path() or summary()).

# One row of what the path should hold at sample size n, from a direct call.
direct_row <- function(res, estimate_term) {
  pick <- function(term) {
    v <- res$value[res$term == term]
    if (length(v)) v else NA_real_
  }
  data.frame(
    estimate = pick(estimate_term),
    lower_limit = pick("lower_limit"),
    upper_limit = pick("upper_limit"),
    current_n = attr(res, "current_n"),
    criterion_n = attr(res, "criterion_n"),
    is_satisfied = attr(res, "is_satisfied")
  )
}

expect_path_equals_direct <- function(path, direct, estimate_term) {
  rows <- do.call(rbind, lapply(direct, direct_row,
                                estimate_term = estimate_term))
  got <- as.data.frame(path)
  expect_identical(got$estimate, rows$estimate)
  expect_identical(got$lower_limit, rows$lower_limit)
  expect_identical(got$upper_limit, rows$upper_limit)
  expect_equal(got$current_n, rows$current_n)
  expect_equal(got$criterion_n, rows$criterion_n)
  expect_identical(got$is_satisfied, rows$is_satisfied)
}

# The smallest sample size a check accepts, found by calling it directly.
first_accepted <- function(call_at) {
  for (n in 1:50) {
    ok <- tryCatch({
      call_at(n)
      TRUE
    }, error = function(e) FALSE)
    if (ok) return(n)
  }
  NA_integer_
}

# A check function with a variance flag, a note, and no interval, built on
# the object contract alone: the rule is met once n reaches `target`, except
# that below `usable_from` the variance estimate is flagged as unusable (as
# a check function does when the estimate is zero, negative, or not finite).
toy_check <- function(x, target, usable_from = 0) {
  n <- length(x)
  if (n < 2L) stop("need at least 2 observations", call. = FALSE)
  unusable <- n < usable_from
  new_smsd_check(
    term = c("mean", "criterion_n", "current_n"),
    value = c(mean(x), target, n),
    method = "toy procedure",
    is_satisfied = !unusable && n >= target,
    current_n = n, criterion_n = target,
    extra = list(variance_unusable = unusable),
    notes = if (unusable) .note_variance_unusable() else NULL
  )
}

# A two group risk procedure that reports the total sample size, as
# seq_mr_smd() does: met once the total reaches `target`.
toy_two_group <- function(group_1, group_2, target) {
  n <- length(group_1) + length(group_2)
  new_smsd_check(
    term = c("difference", "optimal_n", "current_n"),
    value = c(mean(group_1) - mean(group_2), target, n),
    method = "toy two group procedure", is_satisfied = n >= target,
    current_n = n, criterion_n = target
  )
}

test_that("one sample: each row equals a direct call on the first n", {
  withr::local_seed(113)
  x <- rgamma(60, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.2, start = "earliest")
  expect_s3_class(path, c("smsd_path", "data.frame"), exact = TRUE)
  expect_identical(names(path),
                   c("n", "estimate", "lower_limit", "upper_limit",
                     "current_n", "criterion_n", "is_satisfied",
                     "variance_unusable"))
  # The path starts at the smallest n the check accepts and runs to the
  # last observation.
  first <- first_accepted(function(n) seq_aipe_cv(x[seq_len(n)], width = 0.2))
  expect_gt(first, 1L)
  expect_identical(path$n, first:60)
  direct <- lapply(first:60, function(n) {
    seq_aipe_cv(x[seq_len(n)], width = 0.2)
  })
  expect_path_equals_direct(path, direct, "cv")
  # The same quantities summary() gives for each direct call.
  via_summary <- do.call(rbind, lapply(direct, summary))
  expect_identical(path$estimate, via_summary$estimate)
  expect_identical(path$lower_limit, via_summary$lower_limit)
  expect_identical(path$upper_limit, via_summary$upper_limit)
  expect_equal(path$current_n, via_summary$current_n)
  expect_equal(path$criterion_n, via_summary$criterion_n)
  expect_identical(path$is_satisfied, via_summary$is_satisfied)

  # By default the path starts at the pilot sample size and is otherwise
  # the same path.
  from_pilot <- seq_path(seq_aipe_cv, x, width = 0.2)
  pilot <- seq_aipe_cv_pilot(width = 0.2)
  expect_identical(from_pilot$n, pilot:60)
  later <- path[path$n >= pilot, ]
  rownames(later) <- NULL
  expect_identical(as.data.frame(from_pilot), later)
})

test_that("two groups: named and unnamed lists equal direct calls", {
  withr::local_seed(113)
  g1 <- rnorm(50, mean = 0.5)
  g2 <- rnorm(50)
  named <- seq_path(seq_aipe_smd, list(group_1 = g1, group_2 = g2),
                    width = 0.9)
  expect_identical(max(named$n), 50L)
  direct <- lapply(named$n, function(n) {
    seq_aipe_smd(g1[seq_len(n)], g2[seq_len(n)], width = 0.9)
  })
  expect_path_equals_direct(named, direct, "smd")
  unnamed <- seq_path(seq_aipe_smd, list(g1, g2), width = 0.9)
  expect_identical(as.data.frame(unnamed), as.data.frame(named))
  # Names decide which argument an element is, whatever their order.
  swapped <- seq_path(seq_aipe_smd, list(group_2 = g2, group_1 = g1),
                      width = 0.9)
  expect_identical(as.data.frame(swapped), as.data.frame(named))
})

test_that("paired data: a list and a data frame equal direct calls", {
  withr::local_seed(113)
  x <- rnorm(45)
  y <- 0.6 * x + rnorm(45, sd = 0.8)
  path_r <- seq_path(seq_aipe_r, list(x, y), width = 0.5,
                     method = "kendall")
  expect_identical(max(path_r$n), 45L)
  direct_r <- lapply(path_r$n, function(n) {
    seq_aipe_r(x[seq_len(n)], y[seq_len(n)], width = 0.5,
               method = "kendall")
  })
  expect_path_equals_direct(path_r, direct_r, "r")
  expect_identical(attr(path_r, "method"), "Kendall correlation")
  # The proper name is capitalized wherever the label is shown.
  expect_identical(format(path_r)[1], "Sequential path: Kendall correlation")
  # Named after the function's own two data arguments, whatever they are
  # called (they are being renamed to x and y in this release), in either
  # order: the same path.
  data_args <- names(formals(seq_aipe_r))[1:2]
  expect_false(any(data_args %in% c("width", "conf_level", "method")))
  named_r <- seq_path(seq_aipe_r, stats::setNames(list(y, x), rev(data_args)),
                      width = 0.5, method = "kendall")
  expect_identical(named_r, path_r)

  # A data frame is a named list of columns: the names, not the order,
  # decide which is x and which is y.
  path_b <- seq_path(seq_aipe_reg_coef, data.frame(y = y, x = x),
                     width = 0.6)
  direct_b <- lapply(path_b$n, function(n) {
    seq_aipe_reg_coef(x[seq_len(n)], y[seq_len(n)], width = 0.6)
  })
  expect_path_equals_direct(path_b, direct_b, "reg_coef")
})

test_that("a list of groups: every group is cut to its first n", {
  withr::local_seed(113)
  groups <- list(rnorm(40, 0.5), rnorm(40), rnorm(40))
  w <- c(1, -0.5, -0.5)
  path <- seq_path(seq_aipe_c, list(groups = groups), contrast_weights = w,
                   width = 1)
  expect_identical(max(path$n), 40L)
  direct <- lapply(path$n, function(n) {
    seq_aipe_c(lapply(groups, function(g) g[seq_len(n)]),
               contrast_weights = w, width = 1)
  })
  expect_path_equals_direct(path, direct, "contrast")
  expect_identical(attr(path, "arguments"),
                   list(contrast_weights = w, width = 1))
})

test_that("a matrix element is cut by rows", {
  withr::local_seed(113)
  m <- cbind(rnorm(30), rnorm(30))
  by_rows <- function(m, width) seq_aipe_r(m[, 1], m[, 2], width = width)
  path <- seq_path(by_rows, m, width = 0.6)
  # A wrapper has no pilot function: the path starts at the smallest n the
  # wrapper accepts.
  expect_identical(attr(path, "start_rule"), "earliest")
  expect_identical(max(path$n), 30L)
  direct <- lapply(path$n, function(n) {
    seq_aipe_r(m[seq_len(n), 1], m[seq_len(n), 2], width = 0.6)
  })
  expect_path_equals_direct(path, direct, "r")
  # A data frame inside a list is cut by rows too.
  frame_rows <- function(d, width) seq_aipe_r(d$a, d$b, width = width)
  framed <- seq_path(frame_rows, list(data.frame(a = m[, 1], b = m[, 2])),
                     width = 0.6)
  expect_identical(as.data.frame(framed), as.data.frame(path))
})

test_that("risk procedures: no interval, and c and d reach the check", {
  withr::local_seed(113)
  x <- rgamma(50, shape = 9, rate = 3)
  # `c` is a prefix of `check` and `d` of `data`: R's partial matching would
  # take them for seq_path()'s own arguments.
  path <- seq_path(seq_mr_cv, x, A = 400000, c = 75)
  expect_identical(max(path$n), 50L)
  direct <- lapply(path$n, function(n) {
    seq_mr_cv(x[seq_len(n)], A = 400000, c = 75)
  })
  expect_path_equals_direct(path, direct, "cv")
  expect_true(all(is.na(path$lower_limit)) && all(is.na(path$upper_limit)))
  expect_null(attr(path, "conf_level"))
  expect_null(attr(path, "width"))
  expect_identical(attr(path, "arguments"), list(A = 400000, c = 75))

  path_d <- seq_path(seq_fwd_norm_mean, x, d = 0.4)
  direct_d <- lapply(path_d$n, function(n) {
    seq_fwd_norm_mean(x[seq_len(n)], d = 0.4)
  })
  expect_path_equals_direct(path_d, direct_d, "mean")
  # The one sample argument by its own name, whatever it is called (the
  # classical functions' `data` is being renamed to `x` in this release).
  data_arg <- names(formals(seq_fwd_norm_mean))[1]
  expect_identical(seq_path(seq_fwd_norm_mean,
                            stats::setNames(list(x), data_arg), d = 0.4),
                   path_d)
})

test_that("current_n is the check's own scale, n the stage", {
  withr::local_seed(113)
  g1 <- rnorm(40, 0.5)
  g2 <- rnorm(40)
  path <- seq_path(seq_mr_smd, list(g1, g2), A = 2000, c = 1)
  direct <- lapply(path$n, function(n) {
    seq_mr_smd(g1[seq_len(n)], g2[seq_len(n)], A = 2000, c = 1)
  })
  expect_path_equals_direct(path, direct, "smd")
  # seq_mr_smd() reports the total over both groups. (seq_report() words
  # its sample size as a total for this method; if this changes, change the
  # n_scale of its entry in .report_methods() with it.)
  expect_equal(path$current_n, 2 * path$n)
  toy <- seq_path(toy_two_group, list(g1, g2), target = 30)
  expect_equal(toy$current_n, 2 * toy$n)
  expect_identical(attr(toy, "stop_n"), 15L)
})

test_that("the ways of calling seq_path() agree", {
  withr::local_seed(113)
  x <- rgamma(40, shape = 9, rate = 3)
  ref <- seq_path(seq_mr_cv, x, A = 400000, c = 75, by = 3)
  expect_identical(seq_path(data = x, check = seq_mr_cv, A = 400000, c = 75,
                            by = 3), ref)
  expect_identical(seq_path(by = 3, A = 400000, x, c = 75,
                            check = seq_mr_cv), ref)
  expect_identical(do.call(seq_path, list(seq_mr_cv, x, A = 400000, c = 75,
                                          by = 3)), ref)
  wrapper <- function(...) seq_path(seq_mr_cv, x, ...)
  expect_identical(wrapper(A = 400000, c = 75, by = 3), ref)
  expect_identical(lapply(list(x), seq_path, check = seq_mr_cv, A = 400000,
                          c = 75, by = 3)[[1]], ref)
  expect_identical(Map(seq_path, list(seq_mr_cv), list(x), A = 400000,
                       c = 75, by = 3)[[1]], ref)
  # Each argument is evaluated once, `start` and `by` included.
  counts <- new.env()
  counts$n <- 0L
  tick <- function(value) {
    counts$n <- counts$n + 1L
    value
  }
  expect_identical(seq_path(tick(seq_mr_cv), tick(x), A = tick(400000),
                            c = tick(75), by = tick(3), start = tick(18)),
                   seq_path(seq_mr_cv, x, A = 400000, c = 75, by = 3,
                            start = 18))
  expect_identical(counts$n, 6L)
  # What the native pipe makes of x |> seq_path(check = seq_mr_cv, ...).
  expect_identical(seq_path(x, check = seq_mr_cv, A = 400000, c = 75,
                            by = 3), ref)
  # Unnamed arguments after the data go to the check in order.
  expect_identical(as.data.frame(seq_path(seq_aipe_cv, x, 0.3)),
                   as.data.frame(seq_path(seq_aipe_cv, x, width = 0.3)))
})

test_that("stop_n is the first n at which a hand written loop stops", {
  withr::local_seed(113)
  x <- rgamma(120, shape = 9, rate = 3)
  # The loop every vignette writes: start at the pilot, add one observation,
  # check, and stop the first time the rule is met.
  n <- seq_aipe_cv_pilot(width = 0.15)
  repeat {
    out <- seq_aipe_cv(x[seq_len(n)], width = 0.15)
    if (attr(out, "is_satisfied")) break
    n <- n + 1L
  }
  path <- seq_path(seq_aipe_cv, x, width = 0.15)
  expect_identical(attr(path, "start"), seq_aipe_cv_pilot(width = 0.15))
  expect_identical(attr(path, "stop_n"), n)
  expect_identical(attr(path, "final_check"), out)
  expect_false(any(path$is_satisfied[path$n < n]))
  expect_true(path$is_satisfied[path$n == n])
  # The path goes on to the last observation, for illustration.
  expect_identical(max(path$n), 120L)
  expect_gt(max(path$n), n)
})

test_that("a rule never met gives stop_n NA and the last check", {
  withr::local_seed(113)
  x <- rgamma(30, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.02, start = "earliest")
  expect_identical(attr(path, "stop_n"), NA_integer_)
  expect_false(any(path$is_satisfied))
  expect_identical(attr(path, "final_check"), seq_aipe_cv(x, width = 0.02))
})

test_that("start: the pilot by default, and the other choices", {
  withr::local_seed(113)
  x <- rgamma(80, shape = 9, rate = 3)
  # The pilot of Kelley, Darku, and Chattopadhyay (2018) is
  # max(4, ceiling(2 z / width)): 20 at 95% and 26 at 99% for width 0.2.
  expect_identical(ceiling(2 * stats::qnorm(0.975) / 0.2), 20)
  expect_identical(ceiling(2 * stats::qnorm(0.995) / 0.2), 26)
  default <- seq_path(seq_aipe_cv, x, width = 0.2)
  expect_identical(attr(default, "start"), 20L)
  expect_identical(attr(default, "start_rule"), "pilot")
  expect_identical(default$n, 20:80)
  # The pilot function gets the arguments it shares with the check, however
  # they were given: by name, by position, or by a partial name.
  at_99 <- seq_path(seq_aipe_cv, x, width = 0.2, conf_level = 0.99)
  expect_identical(attr(at_99, "start"), 26L)
  expect_identical(attr(at_99, "conf_level"), 0.99)
  expect_identical(attr(seq_path(seq_aipe_cv, x, 0.2, 0.99), "start"), 26L)
  expect_identical(attr(seq_path(seq_aipe_cv, x, wid = 0.2, conf = 0.99),
                        "start"), 26L)
  expect_identical(seq_path(seq_aipe_cv, x, width = 0.2, start = "pilot"),
                   default)
  by_function <- seq_path(seq_aipe_cv, x, width = 0.2,
                          start = seq_aipe_cv_pilot)
  expect_identical(by_function, default)
  own_pilot <- seq_path(seq_aipe_cv, x, width = 0.2,
                        start = function(width, conf_level) 31)
  expect_identical(attr(own_pilot, "start"), 31L)
  expect_identical(attr(own_pilot, "start_rule"), "pilot")

  by_number <- seq_path(seq_aipe_cv, x, width = 0.2, start = 25)
  expect_identical(by_number$n, 25:80)
  expect_identical(attr(by_number, "start"), 25L)
  expect_identical(attr(by_number, "start_rule"), "supplied")

  # "earliest": the smallest n the check accepts, here found by direct calls.
  first <- first_accepted(function(n) seq_aipe_cv(x[seq_len(n)], width = 0.2))
  earliest <- seq_path(seq_aipe_cv, x, width = 0.2, start = "earliest")
  expect_identical(attr(earliest, "start"), first)
  expect_identical(attr(earliest, "start_rule"), "earliest")
  expect_lt(first, 20L)

  # A function that is not one of the package's has no pilot: the default
  # is the earliest start, and asking for "pilot" is an error.
  toy <- seq_path(toy_check, x, target = 9)
  expect_identical(attr(toy, "start"), 2L)
  expect_identical(attr(toy, "start_rule"), "earliest")
  expect_error(seq_path(toy_check, x, target = 9, start = "pilot"),
               "No pilot function is known for this check function")
  wrapper <- function(x, width) seq_aipe_cv(x, width = width)
  expect_identical(attr(seq_path(wrapper, x, width = 0.2), "start"), first)
})

test_that("every check function can be given its pilot by default", {
  # The default start calls <check>_pilot() with the arguments it shares
  # with <check>. That works only while every argument a pilot function
  # requires is also an argument of its check function. Where it is not,
  # seq_path(<check>, ...) stops and says how to supply 'start' (tested
  # below). The pairs for which that is so are written down here, so that a
  # new one is noticed: say on the help page of that procedure that 'start'
  # must be supplied, then add the pair.
  needs_start <- list(
    # The squared multiple correlation: the pilot needs the number of
    # predictors, which the check function reads from the data.
    seq_aipe_R2 = "p"
  )
  ns <- asNamespace("SMSD")
  pilots <- ls(ns, pattern = "^seq_.*_pilot$")
  expect_gt(length(pilots), 10L)
  for (pilot_name in pilots) {
    check_name <- sub("_pilot$", "", pilot_name)
    check <- get0(check_name, envir = ns, inherits = FALSE)
    if (!is.function(check) || check_name == "seq_budget") next
    companion <- .path_pilot_companion(check)
    expect_identical(companion$name, pilot_name)
    pilot_formals <- formals(get(pilot_name, envir = ns))
    required <- names(pilot_formals)[vapply(pilot_formals, function(a) {
      is.symbol(a) && !nzchar(as.character(a))
    }, logical(1))]
    expected <- needs_start[[check_name]]
    if (is.null(expected)) expected <- character(0)
    expect_identical(setdiff(required, c(names(formals(check)), "...")),
                     expected, label = pilot_name)
  }
  # Nothing but the package's own function has a companion.
  expect_null(.path_pilot_companion(function(x, width) NULL))
  expect_null(.path_pilot_companion(mean))
})

test_that("a pilot function that cannot be called from '...' says what to do", {
  withr::local_seed(113)
  x <- rgamma(40, shape = 9, rate = 3)
  # The pilot needs an argument the check function does not take.
  expect_error(
    seq_path(seq_aipe_cv, x, width = 0.3, start = function(width, p) p + 9),
    paste0("The pilot function given as 'start' needs 'p', which the check ",
           "function does not take, so the pilot sample size cannot be ",
           "worked out from the arguments in '...'")
  )
  expect_error(
    seq_path(seq_aipe_cv, x, width = 0.3, start = function(width, p) p + 9),
    "supply the number it returns as 'start'"
  )
  # Called by the user it is a number like any other.
  pilot <- function(width, p) p + 9
  expect_identical(attr(seq_path(seq_aipe_cv, x, width = 0.3,
                                 start = pilot(0.3, p = 3)), "start"), 12L)
  # An argument of the package's pilot function given in '...': R's own
  # message would be "unused argument".
  expect_error(
    seq_path(seq_aipe_cv, x, width = 0.3, m_0 = 30),
    paste0("'m_0' is an argument of seq_aipe_cv_pilot\\(\\), not of ",
           "seq_aipe_cv\\(\\), so it cannot be passed through '...'")
  )
  expect_error(seq_path(seq_aipe_cv, x, width = 0.3, m_0 = 30),
               "start = seq_aipe_cv_pilot(..., m_0 = )", fixed = TRUE)
  expect_identical(attr(seq_path(seq_aipe_cv, x, width = 0.3,
                                 start = seq_aipe_cv_pilot(width = 0.3,
                                                           m_0 = 30)),
                        "start"), 30L)
  # An argument of neither is left to R.
  expect_error(seq_path(seq_aipe_cv, x, width = 0.3, zzz = 1),
               "unused argument")
  # A pilot function that fails: its message is closed with a period before
  # the advice follows.
  expect_error(
    seq_path(seq_aipe_cv, x, width = 0.3,
             start = function(width) stop("no pilot today")),
    "stopped with an error: no pilot today\\. It is called with"
  )
  # More than one sample size, equal or not: the reason given is true of
  # both.
  for (sizes in list(c(12, 9), c(18, 18))) {
    expect_error(
      seq_path(seq_aipe_cv, x, width = 0.3, start = function(width) sizes),
      paste0("returned more than one sample size \\(", sizes[1], ", ",
             sizes[2], "\\), one for each group, as it does for a design ",
             "whose groups are sampled at different rates")
    )
  }
})

test_that("the default start matters for a rule with no correction term", {
  # The classical minimum risk rule for a normal mean, n >= sqrt(A / c) s,
  # has no correction term, so from the smallest sample size it can be met
  # by chance where the procedure, which starts at its pilot, would not
  # stop. The path starts at the pilot by default.
  withr::local_seed(113)
  pilot <- seq_mr_norm_mean_pilot(A = 100, c = 1)
  before_pilot <- 0L
  for (i in 1:40) {
    x <- stats::rnorm(40)
    default <- seq_path(seq_mr_norm_mean, x, A = 100, c = 1)
    earliest <- seq_path(seq_mr_norm_mean, x, A = 100, c = 1,
                         start = "earliest")
    expect_identical(attr(default, "start"), pilot)
    expect_identical(attr(default, "start_rule"), "pilot")
    expect_gte(attr(default, "stop_n"), pilot)
    # An earlier start checks every stage the later one does, and more, so
    # it can only stop at the same sample size or sooner.
    expect_lte(attr(earliest, "stop_n"), attr(default, "stop_n"))
    if (attr(earliest, "stop_n") < pilot) before_pilot <- before_pilot + 1L
  }
  # Started at the smallest sample size, 2, the path stops before the pilot
  # of 4 in about a fifth of data sets: 57 and 60 of 300 drawn as here with
  # seeds 1 and 2. The chance that none of 40 does is then about
  # 0.805^40 = 0.0002. (A pilot of 2 would leave nothing before it.)
  if (pilot > 2L) expect_gt(before_pilot, 0L)
})

test_that("the correction term makes stop_n the same from either start", {
  # n >= (2 z / w)^2 (xi^2 + 1 / n) cannot hold for n < 2 z / w, the pilot
  # size, so starting below the pilot cannot stop the path earlier.
  withr::local_seed(113)
  for (i in 1:8) {
    x <- rgamma(70, shape = 4, rate = 2)
    early <- seq_path(seq_aipe_cv, x, width = 0.3, start = "earliest")
    pilot <- seq_path(seq_aipe_cv, x, width = 0.3)
    expect_identical(attr(early, "stop_n"), attr(pilot, "stop_n"))
    expect_lt(attr(early, "start"), attr(pilot, "start"))
  }
})

test_that("by: stages are start, start + by, ..., and the tail is unused", {
  withr::local_seed(113)
  x <- rgamma(60, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.2, start = 10, by = 7)
  expect_identical(path$n, seq.int(10L, 60L, by = 7L))
  expect_identical(max(path$n), 59L)
  expect_identical(attr(path, "by"), 7L)
  expect_identical(attr(path, "n_observations"), 60L)
  direct <- lapply(path$n, function(n) seq_aipe_cv(x[seq_len(n)], width = 0.2))
  expect_path_equals_direct(path, direct, "cv")
  # With batches the stop is the first checked n at which the rule is met.
  every <- seq_path(seq_aipe_cv, x, width = 0.2, start = 10)
  expect_gte(attr(path, "stop_n"), attr(every, "stop_n"))
  expect_identical(attr(path, "stop_n"),
                   path$n[which(path$is_satisfied)[1]])

  one <- seq_path(seq_aipe_cv, x, width = 0.2, start = 60)
  expect_identical(nrow(one), 1L)
  expect_identical(one$n, 60L)
})

test_that("the attributes describe the path", {
  withr::local_seed(113)
  x <- rgamma(60, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.2, conf_level = 0.9)
  expect_identical(attr(path, "method"), "coefficient of variation")
  expect_identical(attr(path, "estimate_term"), "cv")
  expect_identical(attr(path, "conf_level"), 0.9)
  expect_identical(attr(path, "width"), 0.2)
  expect_identical(attr(path, "start"), path$n[1])
  expect_identical(attr(path, "start_rule"), "pilot")
  expect_identical(attr(path, "by"), 1L)
  expect_identical(attr(path, "n_observations"), 60L)
  expect_identical(attr(path, "arguments"),
                   list(width = 0.2, conf_level = 0.9))
  expect_s3_class(attr(path, "final_check"), "smsd_check")
  expect_type(path$n, "integer")
  expect_type(path$is_satisfied, "logical")
})

test_that("variance_unusable follows the check, and such a stage is unmet", {
  x <- c(2, 4, 6, 8, 10, 12, 14, 16)
  path <- seq_path(toy_check, x, target = 3, usable_from = 6)
  expect_identical(path$n, 2:8)
  expect_identical(path$variance_unusable, path$n < 6)
  expect_identical(path$is_satisfied, path$n >= 6)
  expect_identical(attr(path, "stop_n"), 6L)
  expect_identical(path$estimate, cumsum(x)[2:8] / 2:8)
  expect_identical(attr(attr(path, "final_check"), "notes"), character(0))

  # The print: the header counts the stages, the table shows the column,
  # and a line under the table names the stages (n = 2 to 5) and says what
  # the column means, in the words of the package's rule (zero, negative, or
  # not finite).
  txt <- format(path)
  expect_true("  Variance        : not usable at 4 of 7 stages" %in% txt)
  expect_true(grepl("variance_unusable$", txt[grepl("criterion_n", txt)]))
  expect_identical(
    txt[seq.int(length(txt) - 7L, length(txt))],
    c("",
      "  At n = 2 to 5 the variance estimate was zero, negative, or not finite",
      "  (column variance_unusable), so the stopping rule was reported as not",
      "  met there.",
      "",
      "  The estimate to report is the one at n = 6. Later rows show what the",
      "  rule would have said had sampling continued; they are not grounds for",
      "  choosing a different stopping point.")
  )
  expect_false(any(grepl("floor", txt)))
  # One such stage: n = 2.
  one <- format(seq_path(toy_check, x, target = 3, usable_from = 3))
  expect_true("  Variance        : not usable at 1 of 7 stages" %in% one)
  expect_true(paste("  At n = 2 the variance estimate was zero, negative, or",
                    "not finite") %in% one)
  only <- format(seq_path(toy_check, x[1:2], target = 3, usable_from = 3))
  expect_true("  Variance        : not usable at 1 of 1 stage" %in% only)
  # Checked after every second observation the stages are n = 2, 4, 6, and
  # 8. With no usable variance below n = 8 the first three are named as a
  # run of stages two observations apart.
  batches <- seq_path(toy_check, x, target = 3, usable_from = 8, by = 2)
  expect_identical(batches$n[batches$variance_unusable], c(2L, 4L, 6L))
  expect_match(paste(trimws(format(batches)), collapse = " "),
               paste("At n = 2 to 6 by 2 the variance estimate was zero,",
                     "negative, or not finite (column variance_unusable)"),
               fixed = TRUE)
  # None: no header line, no column, no line under the table.
  none <- format(seq_path(toy_check, x, target = 3))
  expect_false(any(grepl("Variance|variance", none)))
  expect_identical(summary(path)$variance_unusable, FALSE)
  expect_identical(
    summary(seq_path(toy_check, x, target = 3, usable_from = 99))
    $variance_unusable,
    TRUE
  )

  # A check that reports no variance flag gives NA, not FALSE.
  no_flag <- function(x) {
    new_smsd_check("mean", mean(x), "no flag", is_satisfied = FALSE,
                   current_n = length(x), criterion_n = 99)
  }
  expect_identical(seq_path(no_flag, x)$variance_unusable, rep(NA, 8))
  expect_false(any(grepl("Variance|variance", format(seq_path(no_flag, x)))))
})

test_that("the print names unusable stages that the table leaves out", {
  # Thirty observations checked from n = 2: the table shows n = 2, 3, 4 and
  # 28, 29, 30. The variance estimate is unusable at n = 10, 11, and 12
  # only, so none of the rows shown is marked, and the reader learns which
  # stages are meant from the line under the table.
  middle <- function(x, target) {
    n <- length(x)
    unusable <- n >= 10L && n <= 12L
    new_smsd_check(c("mean", "criterion_n", "current_n"),
                   c(mean(x), target, n), "toy procedure",
                   is_satisfied = !unusable && n >= target, current_n = n,
                   criterion_n = target,
                   extra = list(variance_unusable = unusable))
  }
  path <- seq_path(middle, as.numeric(1:30), target = 11, start = 2)
  expect_identical(path$n[path$variance_unusable], 10:12)
  expect_identical(attr(path, "stop_n"), 13L)
  txt <- format(path)
  expect_true("  Variance        : not usable at 3 of 29 stages" %in% txt)
  table_rows <- txt[grepl("^ +[0-9]+ +[0-9.]+ ", txt)]
  expect_identical(as.integer(sub("^ +([0-9]+) .*$", "\\1", table_rows)),
                   c(2L, 3L, 4L, 28L, 29L, 30L))
  expect_true(all(grepl(" FALSE +FALSE$", table_rows[1:3])))
  expect_true(all(grepl(" TRUE +FALSE$", table_rows[4:6])))
  expect_match(
    paste(trimws(txt), collapse = " "),
    paste("At n = 10 to 12 the variance estimate was zero, negative, or not",
          "finite (column variance_unusable), so the stopping rule was",
          "reported as not met there."), fixed = TRUE
  )

  # How the stages are worded. A run of three or more stages in a row is
  # given by its ends; one or two are written out.
  n <- 1:40
  expect_identical(.path_stages_text(n, 10L, by = 1L), "n = 10")
  expect_identical(.path_stages_text(n, 10:11, by = 1L), "n = 10 and 11")
  expect_identical(.path_stages_text(n, 10:12, by = 1L), "n = 10 to 12")
  expect_identical(.path_stages_text(n, c(3:6, 10L, 14:16), by = 1L),
                   "n = 3 to 6, 10, and 14 to 16")
  expect_identical(.path_stages_text(n, c(3L, 4L, 8L), by = 1L),
                   "n = 3, 4, and 8")
  expect_identical(.path_stages_text(n, c(3L, 8L), by = 1L), "n = 3 and 8")
  # The rows are positions among the stages, not sample sizes: for a path
  # that starts at n = 7, rows 4 and 5 are n = 10 and 11.
  expect_identical(.path_stages_text(7:60, 4:5, by = 1L), "n = 10 and 11")
  # Stages five observations apart: rows 1 to 3 are n = 10, 15, and 20, a
  # run of stages although not of sample sizes, and "by 5" says so.
  fives <- seq(10L, 60L, by = 5L)
  expect_identical(.path_stages_text(fives, 1:3, by = 5L), "n = 10 to 20 by 5")
  expect_identical(.path_stages_text(fives, c(1L, 3L), by = 5L),
                   "n = 10 and 20")
  expect_identical(.path_stages_text(fives, c(1:3, 7L, 9:11), by = 5L),
                   "n = 10 to 20 by 5, 40, and 50 to 60 by 5")
  # Eight pieces are still listed; more are counted, with where to look.
  expect_identical(.path_stages_text(n, seq(1L, 15L, by = 2L), by = 1L),
                   "n = 1, 3, 5, 7, 9, 11, 13, and 15")
  expect_identical(
    .path_stages_text(n, seq(1L, 17L, by = 2L), by = 1L),
    "9 stages between n = 1 and n = 17 (as.data.frame() shows which)"
  )
  scattered <- function(x) {
    n <- length(x)
    new_smsd_check("mean", mean(x), "toy procedure", is_satisfied = FALSE,
                   current_n = n, criterion_n = 99,
                   extra = list(variance_unusable = n %% 2L == 0L))
  }
  many <- paste(trimws(format(seq_path(scattered, as.numeric(1:30),
                                       start = 2))), collapse = " ")
  expect_match(many,
               paste("At 15 stages between n = 2 and n = 30 (as.data.frame()",
                     "shows which) the variance estimate was zero, negative,",
                     "or not finite (column variance_unusable), so the",
                     "stopping rule was reported as not met there."),
               fixed = TRUE)
})

test_that("the header gives no interval for a stage with no usable variance", {
  # A check that holds limits at every stage, as the fixed width procedures
  # do, and flags its variance estimate below `usable_from`. The rule is
  # never met here, so the header describes the last stage: 8 observations
  # with mean 9, limits 8 and 10.
  with_limits <- function(x, usable_from) {
    n <- length(x)
    new_smsd_check(c("lower_limit", "mean", "upper_limit"),
                   c(mean(x) - 1, mean(x), mean(x) + 1), "toy procedure",
                   is_satisfied = FALSE, current_n = n, criterion_n = 99,
                   conf_level = 0.9, width = 2,
                   extra = list(variance_unusable = n < usable_from))
  }
  x <- c(2, 4, 6, 8, 10, 12, 14, 16)
  usable <- format(seq_path(with_limits, x, usable_from = 5, start = 2))
  expect_true("  Latest estimate : mean = 9 [8, 10]" %in% usable)
  flagged <- seq_path(with_limits, x, usable_from = 99, start = 2)
  txt <- format(flagged)
  expect_true("  Latest estimate : mean = 9" %in% txt)
  expect_false(any(grepl("mean = 9 [", txt, fixed = TRUE)))
  # The table still shows what the check returned at that stage.
  expect_identical(flagged$upper_limit[flagged$n == 8L], 10)
  expect_true(any(grepl("^ +8 +9 +8 +10 ", txt)))
})

test_that("a bound on the interval's length is not called a target width", {
  # The bounded length interval for a median stops when the interval is no
  # longer than 2 d; the interval reported is usually shorter. A check with
  # that method label, built on the object contract: width 1.2 is the bound,
  # and the interval at the stop (n = 6) runs from 6.6 to 7.4, a length of
  # 0.8.
  bounded <- function(x, d) {
    n <- length(x)
    new_smsd_check(c("lower_limit", "median", "upper_limit"),
                   c(stats::median(x) - 0.4, stats::median(x),
                     stats::median(x) + 0.4),
                   "bounded length confidence interval (distribution free)",
                   is_satisfied = n >= 6L, current_n = n,
                   criterion_n = NA_real_, conf_level = 0.9, width = 2 * d)
  }
  x <- c(2, 4, 6, 8, 10, 12, 14, 16)
  path <- seq_path(bounded, x, d = 0.6, start = 2)
  expect_identical(attr(path, "width"), 1.2)
  txt <- format(path)
  expect_true("  Length at most  : 1.2" %in% txt)
  expect_false(any(grepl("Target width", txt, fixed = TRUE)))
  expect_true("  Estimate there  : median = 7 [6.6, 7.4]" %in% txt)
  expect_identical(.path_subtitle(path, as.data.frame(path)),
                   "Length at most 1.2, 90% confidence; rule first met at n = 6")
  # The three outputs for one path say the same thing of the number.
  expect_match(as.character(seq_report(path)),
               "with a 90% confidence interval no longer than 1.2.",
               fixed = TRUE)
  # The package's own procedure carries that label.
  withr::local_seed(113)
  y <- stats::rgamma(60, shape = 9, rate = 3)
  own <- seq_path(seq_blci_dfr_median, y, d = 0.6)
  expect_true("  Length at most  : 1.2" %in% format(own))
  expect_match(.path_subtitle(own, as.data.frame(own)),
               "^Length at most 1\\.2, 95% confidence; rule first met at n = ")
  # Every other width is a target: a procedure aimed at a width, and a
  # method label the package does not know.
  cv <- seq_path(seq_aipe_cv, y, width = 0.2)
  expect_true("  Target width    : 0.2" %in% format(cv))
  expect_identical(
    .path_subtitle(cv, as.data.frame(cv)),
    paste0("Target width 0.2, 95% confidence; rule first met at n = ",
           attr(cv, "stop_n"))
  )
  other <- function(x, d) {
    out <- bounded(x, d)
    attr(out, "method") <- "a method of the future"
    out
  }
  expect_true("  Target width    : 1.2" %in%
                format(seq_path(other, x, d = 0.6, start = 2)))
})

test_that("a package check: stages with a variance of zero are unusable", {
  # Ratings that begin with six equal scores. The sample variance of the
  # first n is zero exactly for n = 3 to 6, so by the package's rule no
  # variance has been estimated at those stages and the stopping rule is
  # not met there; from n = 7 on a variance exists.
  ratings <- c(5, 5, 5, 5, 5, 5, 4, 6, 5, 7, 3, 5, 6, 4, 5, 6, 4, 5, 7, 3)
  running_var <- vapply(3:20, function(n) stats::var(ratings[seq_len(n)]),
                        numeric(1))
  expect_identical(which(running_var == 0), 1:4)
  path <- seq_path(seq_fwd_norm_mean, ratings, d = 1.2, start = 3)
  expect_identical(path$n, 3:20)
  expect_identical(path$variance_unusable, running_var == 0)
  expect_false(any(path$is_satisfied[path$variance_unusable]))
  expect_true(is.na(attr(path, "stop_n")) || attr(path, "stop_n") >= 7L)
  txt <- format(path)
  expect_true("  Variance        : not usable at 4 of 18 stages" %in% txt)
  # The check at such a stage says so in the package's words.
  at_five <- seq_fwd_norm_mean(ratings[1:5], d = 1.2)
  expect_true(attr(at_five, "variance_unusable"))
  expect_false(attr(at_five, "is_satisfied"))
  expect_true(.note_variance_unusable() %in% attr(at_five, "notes"))
})

test_that("a pilot function gets its arguments by name, in any order", {
  x <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
  # A check with a correction exponent gamma, and two pilot functions that
  # differ only in the order of m_0 and gamma (the order of those two is
  # being changed in the package's minimum risk pilot functions). The pilot
  # is m_0 + 10 * gamma, so the two arguments cannot be exchanged unseen:
  # with m_0 = 4 it is 4 + 10 * 0.2 = 6 for gamma = 0.2 and 4 + 10 * 0.5 = 9
  # for the default gamma.
  with_gamma <- function(x, target, gamma = 0.5) toy_check(x, target = target)
  gamma_first <- function(target, gamma = 0.5, m_0 = 4) m_0 + 10 * gamma
  m_0_first <- function(target, m_0 = 4, gamma = 0.5) m_0 + 10 * gamma
  for (pilot in list(gamma_first, m_0_first)) {
    # gamma by name, by position (the check's third argument), and left at
    # the pilot function's default.
    expect_identical(attr(seq_path(with_gamma, x, target = 3, gamma = 0.2,
                                   start = pilot), "start"), 6L)
    expect_identical(attr(seq_path(with_gamma, x, 3, 0.2, start = pilot),
                          "start"), 6L)
    expect_identical(attr(seq_path(with_gamma, x, target = 3, start = pilot),
                          "start"), 9L)
  }
  # The package's two minimum risk procedures: the default start is the
  # pilot function called by name with the gamma given to the check.
  withr::local_seed(113)
  y <- stats::rgamma(60, shape = 9, rate = 3)
  expect_identical(
    attr(seq_path(seq_mr_cv, y, A = 400000, c = 75, gamma = 0.3), "start"),
    seq_mr_cv_pilot(A = 400000, c = 75, gamma = 0.3)
  )
  expect_false(seq_mr_cv_pilot(A = 400000, c = 75, gamma = 0.3) ==
                 seq_mr_cv_pilot(A = 400000, c = 75))
  g1 <- stats::rnorm(60, mean = 0.5)
  g2 <- stats::rnorm(60)
  expect_identical(
    attr(seq_path(seq_mr_smd, list(g1, g2), A = 5000, c = 1, gamma = 0.3),
         "start"),
    seq_mr_smd_pilot(A = 5000, c = 1, gamma = 0.3)
  )
  expect_false(seq_mr_smd_pilot(A = 5000, c = 1, gamma = 0.3) ==
                 seq_mr_smd_pilot(A = 5000, c = 1))
})

test_that("a warning from the check is reported once, with its stages", {
  x <- c(2, 4, 6, 8, 10, 12, 14, 16)
  noisy <- function(x, target) {
    if (length(x) %% 2L == 0L) warning("an even sample size", call. = FALSE)
    toy_check(x, target = target)
  }
  expect_warning(
    path <- seq_path(noisy, x, target = 5),
    "gave this warning at 4 of 7 stages (first at n = 2): an even sample size",
    fixed = TRUE
  )
  expect_identical(attr(path, "stop_n"), 5L)
  expect_silent(seq_path(toy_check, x, target = 5))
})

test_that("seq_path() refuses what is not a check function", {
  x <- c(2, 4, 6, 8, 10, 12)
  expect_error(seq_path("seq_aipe_cv", x, width = 0.2),
               "'check' must be a sequential check function")
  expect_error(seq_path(function(v) mean(v), x),
               "returns an object of class 'smsd_check'")
  expect_error(seq_path(), "'check' must be supplied")
  expect_error(seq_path(seq_aipe_cv, width = 0.2), "'data' must be supplied")
  expect_error(seq_path(seq_aipe_cv, zzz_no_such_object, width = 0.2),
               "zzz_no_such_object")
})

test_that("two stage and budget procedures are refused with advice", {
  withr::local_seed(113)
  x <- rnorm(30)
  expect_error(seq_path(ts_fwd_norm_mean, list(x[1:10], x[11:30]), d = 0.3),
               "does not apply to the two stage procedures")
  expect_error(seq_path(ts_fwd_norm_mean, list(x[1:10], x[11:30]), d = 0.3),
               "once on the pilot sample")
  groups <- list(rnorm(20), rnorm(20), rnorm(20))
  expect_error(
    seq_path(seq_budget, list(groups = groups),
             contrast_weights = c(1, -0.5, -0.5),
             cost_per_group = c(1, 1, 1), budget = 300),
    "does not apply to the budget constrained procedure"
  )
  # A wrapper hides the arguments; the label of the result still says so.
  wrapped <- function(v, d) ts_fwd_norm_mean(v[1:10], v[-(1:10)], d = d)
  expect_error(seq_path(wrapped, x, d = 0.3),
               "does not apply to the two stage procedures")
  wrapped_budget <- function(groups) {
    seq_budget(groups, contrast_weights = c(1, -0.5, -0.5),
               cost_per_group = c(1, 1, 1), budget = 300)
  }
  expect_error(seq_path(wrapped_budget, list(groups = groups)),
               "does not apply to the budget constrained procedure")

  # The package's own functions are known by what they are, whatever their
  # arguments are called: in a stand-in for the namespace, a two stage
  # function and a budget function with arguments of other names.
  stand_in <- new.env()
  stand_in$ts_other <- function(first, second, d) NULL
  stand_in$ts_other_pilot <- function(d) NULL
  stand_in$seq_budget <- function(g, w, costs, total) NULL
  stand_in$seq_other <- function(x, width) NULL
  expect_identical(.path_family(stand_in$ts_other, ns = stand_in),
                   "two_stage")
  expect_identical(.path_family(stand_in$seq_budget, ns = stand_in), "budget")
  expect_identical(.path_family(stand_in$seq_other, ns = stand_in), "")
  # A pilot function is not a check function of either family.
  expect_identical(.path_family(stand_in$ts_other_pilot, ns = stand_in), "")
  # A function of the same shape that is not the package's (it lives in
  # another environment) is not taken for one, and the argument names of
  # the two families still are.
  expect_identical(.path_family(local(function(first, second, d) NULL),
                                ns = stand_in), "")
  expect_identical(.path_family(function(pilot_data, more, d) NULL,
                                ns = stand_in), "two_stage")
  expect_identical(.path_family(function(groups, budget) NULL, ns = stand_in),
                   "budget")
  # Every two stage check function of the package and the budget function.
  ns <- asNamespace("SMSD")
  two_stage <- sub("_pilot$", "", ls(ns, pattern = "^ts_.*_pilot$"))
  expect_gt(length(two_stage), 4L)
  for (name in two_stage) {
    expect_identical(.path_family(get(name, envir = ns)), "two_stage",
                     label = name)
  }
  expect_identical(.path_family(seq_budget), "budget")
  expect_identical(.path_family(seq_aipe_cv), "")
})

test_that("data that do not fit are refused with a message that says why", {
  withr::local_seed(113)
  g1 <- rnorm(30)
  g2 <- rnorm(30)
  expect_error(seq_path(seq_aipe_smd, list(g1, g2[1:25]), width = 0.9),
               "same number of observations.*: 30, .*: 25")
  expect_error(
    seq_path(seq_aipe_c, list(groups = list(g1, g2, g2[1:9])),
             contrast_weights = c(1, -0.5, -0.5), width = 1),
    "'groups'\\[\\[3\\]\\]: 9"
  )
  expect_error(seq_path(seq_aipe_smd, list(a = g1, b = g2), width = 0.9),
               "named a, b, which are not arguments of the check function")
  expect_error(seq_path(seq_aipe_reg_coef, data.frame(u = g1, v = g2),
                        width = 0.9),
               "A data frame is read as a named list of its columns")
  expect_error(seq_path(seq_aipe_smd, list(group_1 = g1, g2), width = 0.9),
               "Name every element of 'data'")
  expect_error(seq_path(seq_aipe_smd, list(group_1 = g1, group_1 = g2),
                        width = 0.9), "distinct names")
  expect_error(seq_path(seq_aipe_cv, list(x = g1), x = 3, width = 0.9),
               "supplied both in 'data' and in '...'")
  expect_error(seq_path(seq_aipe_cv, list(), width = 0.9), "'data' is empty")
  expect_error(seq_path(seq_aipe_cv, NULL, width = 0.9), "got NULL")
  expect_error(seq_path(seq_aipe_cv, numeric(0), width = 0.9),
               "holds no observations")
  expect_error(seq_path(seq_aipe_cv, list(list()), width = 0.9),
               "is an empty list")
  expect_error(seq_path(seq_aipe_cv, list(mean), width = 0.9),
               "must be a vector, or a matrix or data frame")
  expect_error(
    seq_path(seq_aipe_c, list(groups = list(g1, g2, "a", NULL)),
             contrast_weights = c(1, -1, 0, 0), width = 1),
    "Member 4 of 'groups' must be a vector"
  )
  expect_error(seq_path(seq_aipe_cv, array(1, c(2, 2, 2)), width = 0.9),
               "must be a vector, or a matrix or data frame")
  # An unnamed list of groups reads as three data arguments; the message
  # names the argument, not the data, and says how to pass a list of groups.
  err <- tryCatch(
    seq_path(seq_aipe_c, list(g1, g2, g2),
             contrast_weights = c(1, -0.5, -0.5), width = 1),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "list(groups = list(...))", fixed = TRUE)
  expect_lt(nchar(err), 500)
})

test_that("data written after 'data' are refused, not passed whole", {
  withr::local_seed(113)
  g1 <- rnorm(60, mean = 0.5)
  g2 <- rnorm(60)
  # Everything in '...' is passed whole at every stage. seq_mr_smd() accepts
  # groups of unequal size, so the second group written after `data` would
  # be used in full from the first stage on, with no error: this direct call
  # is the first stage such a path would hold, 14 observations from group 1
  # and all 60 from group 2.
  expect_equal(attr(seq_mr_smd(g1[1:14], g2, A = 5000, c = 1), "current_n"),
               74)
  expect_error(
    seq_path(seq_mr_smd, g1, g2, A = 5000, c = 1),
    paste0("Unnamed argument 1 after 'data' \\(the check function's ",
           "'group_2'\\) holds 60 observations, and the check function ",
           "also runs when it is cut along with the data, so it looks like ",
           "data")
  )
  expect_error(seq_path(seq_mr_smd, g1, g2, A = 5000, c = 1),
               "seq_path(check, list(group_1, group_2), ...)", fixed = TRUE)
  expect_error(seq_path(seq_mr_smd, g1, g2, A = 5000, c = 1),
               "wrap it in I() to pass it whole", fixed = TRUE)
  # Given by name it is the same mistake.
  expect_error(seq_path(seq_mr_smd, g1, group_2 = g2, A = 5000, c = 1),
               "The argument 'group_2' holds 60 observations")
  expect_error(seq_path(seq_mr_smd, list(group_1 = g1), g2, A = 5000, c = 1),
               "list(group_1, group_2)", fixed = TRUE)
  # A second group of another length is caught too, and the advice says to
  # trim.
  expect_error(
    seq_path(seq_mr_smd, g1, g2[1:45], A = 5000, c = 1),
    "holds 45 observations \\('data' holds 60\\).*trimmed to a common number"
  )
  # The paired and two group functions, for which the mistake used to end in
  # the advice that a larger 'start' may be needed.
  expect_error(seq_path(seq_aipe_smd, g1, g2, width = 0.9), "looks like data")
  expect_error(seq_path(seq_aipe_r, g1, g2, width = 0.5), "looks like data")
  expect_error(seq_path(seq_aipe_reg_coef, g1, g2, width = 0.6),
               "seq_path(check, list(x, y), ...)", fixed = TRUE)
  expect_error(seq_path(toy_two_group, g1, g2, target = 30),
               "the check function's 'group_2'")
  # A matrix or a data frame of observations is data as well.
  by_rows <- function(m, other, width) {
    seq_aipe_r(m[, 1], other[, 1], width = width)
  }
  expect_error(seq_path(by_rows, cbind(g1), cbind(g2), width = 0.5),
               "holds 60 observations")
  expect_error(seq_path(by_rows, cbind(g1), data.frame(g2), width = 0.5),
               "looks like data")

  # With all of the data in `data` every stage holds n per group.
  right <- seq_path(seq_mr_smd, list(g1, g2), A = 5000, c = 1)
  expect_equal(right$current_n, 2 * right$n)
})

test_that("a vector argument that is not data still reaches the check", {
  withr::local_seed(113)
  # Contrast weights: one for each group, however many observations there
  # are, including as many weights as observations (five groups of five).
  groups <- lapply(1:5, function(i) rnorm(5))
  w <- c(1, -0.25, -0.25, -0.25, -0.25)
  path <- seq_path(seq_aipe_c, list(groups = groups), contrast_weights = w,
                   width = 3, start = 4)
  expect_identical(path$n, 4:5)
  direct <- lapply(4:5, function(n) {
    seq_aipe_c(lapply(groups, function(g) g[seq_len(n)]),
               contrast_weights = w, width = 3)
  })
  expect_path_equals_direct(path, direct, "contrast")
  expect_identical(attr(path, "arguments")$contrast_weights, w)

  # A pair of costs.
  g1 <- rnorm(40, mean = 0.5)
  g2 <- rnorm(40)
  costs <- seq_path(seq_mr_smd, list(g1, g2), A = 2000, c = c(1, 1),
                    start = 12)
  expect_identical(attr(costs, "arguments")$c, c(1, 1))
  expect_identical(costs$estimate[1],
                   summary(seq_mr_smd(g1[1:12], g2[1:12], A = 2000,
                                      c = c(1, 1)))$estimate)

  # An argument that the test takes for data although it is not: a check
  # that runs on two observations and takes a vector of three. I() says it
  # is not data; the check receives it as it was, without the wrapper.
  with_probs <- function(x, probs) {
    new_smsd_check("mean", mean(x), "toy with a vector argument",
                   is_satisfied = length(x) >= 6, current_n = length(x),
                   criterion_n = 6, extra = list(seen = probs))
  }
  probs <- c(0.1, 0.5, 0.9)
  expect_error(
    seq_path(with_probs, g1, probs = probs),
    "The argument 'probs' holds 3 observations \\('data' holds 40\\)"
  )
  kept <- seq_path(with_probs, g1, probs = I(probs))
  expect_identical(attr(kept, "stop_n"), 6L)
  expect_identical(attr(kept, "arguments"), list(probs = probs))
  expect_identical(attr(attr(kept, "final_check"), "seen"), probs)
  # I() takes off only its own class.
  dated <- structure(probs, class = "my_probs")
  kept <- seq_path(with_probs, g1, probs = I(dated))
  expect_identical(attr(attr(kept, "final_check"), "seen"), dated)
  # A single number, a string, and a function are never taken for data.
  scalar_args <- function(x, target, label = "a", f = mean) {
    toy_check(x, target = target)
  }
  expect_identical(attr(seq_path(scalar_args, g1, target = 5, label = "b",
                                 f = stats::median), "stop_n"), 5L)
})

test_that("missing values in the data are refused", {
  withr::local_seed(113)
  x <- rgamma(40, shape = 9, rate = 3)
  y <- rgamma(40, shape = 9, rate = 3)
  holes <- x
  holes[c(3, 7, 9)] <- NA
  holes[12] <- NaN
  # Why: a check that drops missing values reports fewer observations than
  # the stage holds, so the path's n and the check's n would disagree. This
  # is seq_mr_cv() at what would be stage 18.
  expect_equal(attr(seq_mr_cv(holes[1:18], A = 400000, c = 75), "current_n"),
               14)
  expect_error(
    seq_path(seq_mr_cv, holes, A = 400000, c = 75),
    paste0("'data' holds missing values \\(NA or NaN\\): 'data', 4 of 40 ",
           "observations\\. A path is indexed by the number of observations ",
           "collected")
  )
  expect_error(seq_path(seq_mr_cv, holes, A = 400000, c = 75),
               "Remove the incomplete observations first")
  # Asking the check to drop them does not make a path of them.
  expect_error(seq_path(seq_mr_cv, holes, A = 400000, c = 75, na.rm = TRUE),
               "'data' holds missing values")
  # Whatever the start, and for a check of the user's own.
  expect_error(seq_path(toy_check, holes, target = 5, start = 20),
               "'data', 4 of 40 observations")
  # Each element is named with its count.
  expect_error(seq_path(seq_aipe_smd, list(x, holes), width = 0.9),
               "element 2 of 'data', 4 of 40 observations")
  expect_error(
    seq_path(seq_aipe_smd, list(group_1 = holes, group_2 = holes[40:1]),
             width = 0.9),
    "'group_1', 4 of 40 observations; 'group_2', 4 of 40 observations"
  )
  expect_error(
    seq_path(seq_aipe_c, list(groups = list(x, holes, y)),
             contrast_weights = c(1, -0.5, -0.5), width = 1),
    "'groups'\\[\\[2\\]\\], 4 of 40 observations"
  )
  # In a matrix or a data frame an observation is a row: three missing
  # values in two rows are two incomplete observations.
  m <- cbind(x, y)
  m[5, ] <- NA
  m[9, 2] <- NA
  by_rows <- function(m, width) seq_aipe_r(m[, 1], m[, 2], width = width)
  expect_error(seq_path(by_rows, m, width = 0.6),
               "'data', 2 of 40 observations")
  expect_error(seq_path(by_rows, list(as.data.frame(m)), width = 0.6),
               "'data', 2 of 40 observations")
  expect_error(seq_path(seq_aipe_reg_coef, as.data.frame(m), width = 0.6),
               "'y', 2 of 40 observations")
  # Complete data pass.
  expect_s3_class(seq_path(seq_mr_cv, x, A = 400000, c = 75), "smsd_path")
})

test_that("an error of the check is reported with the sample size", {
  withr::local_seed(113)
  x <- rgamma(30, shape = 9, rate = 3)
  # A missing argument of the check fails on all the data, at once. R's
  # message has no period; one is added.
  expect_error(seq_path(seq_aipe_cv, x),
               "stopped with an error when given all 30 observations")
  expect_error(seq_path(toy_check, x), "with no default\\.$")
  expect_error(seq_path(seq_aipe_cv, 3, width = 0.3),
               "when given the one observation in hand: ")
  accepts_one <- function(x) {
    new_smsd_check("mean", mean(x), "accepts one", is_satisfied = FALSE,
                   current_n = length(x), criterion_n = 9)
  }
  expect_error(seq_path(accepts_one, 3, start = 2),
               "'start' is 2, but only 1 observation is in hand")
  expect_identical(nrow(seq_path(accepts_one, 3)), 1L)
  # A start below what the check accepts.
  expect_error(seq_path(seq_aipe_cv, x, width = 0.3, start = 2),
               "stopped with an error at n = 2: .*larger 'start'")
  # An error part of the way along the path.
  fragile <- function(x, target) {
    if (length(x) == 7L) stop("seven is unlucky", call. = FALSE)
    toy_check(x, target = target)
  }
  expect_error(seq_path(fragile, x, target = 5, start = 2),
               "stopped with an error at n = 7: seven is unlucky")
  # Found by trying, the default start steps over the failure at n = 7 only
  # if it comes first; here n = 2 already works, so the path fails at 7.
  expect_error(seq_path(fragile, x, target = 5), "at n = 7")
})

test_that("start and by are validated", {
  withr::local_seed(113)
  x <- rgamma(30, shape = 9, rate = 3)
  expect_error(seq_path(seq_aipe_cv, x, width = 0.3, start = 2.5),
               "'start' must be a single whole number")
  expect_error(seq_path(seq_aipe_cv, x, width = 0.3, start = NA),
               "'start' must be a single whole number")
  expect_error(seq_path(seq_aipe_cv, x, width = 0.3, start = 31),
               "'start' is 31, but only 30 observations are in hand")
  # The pilot for a width of 0.01 is ceiling(2 z / 0.01) = 392 observations:
  # with 30 in hand the procedure has not reached its first check.
  expect_identical(ceiling(2 * stats::qnorm(0.975) / 0.01), 392)
  expect_error(seq_path(seq_aipe_cv, x, width = 0.01),
               "The pilot sample size is 392, but only 30 observations are")
  expect_error(seq_path(seq_aipe_cv, x, width = 0.01),
               "start = \"earliest\"", fixed = TRUE)
  expect_s3_class(seq_path(seq_aipe_cv, x, width = 0.01, start = "earliest"),
                  "smsd_path")
  expect_error(seq_path(seq_aipe_cv, x, width = 0.3,
                        start = function(width) stop("no pilot")),
               "pilot function given as 'start' stopped with an error")
  expect_error(seq_path(seq_aipe_cv, x, width = 0.3,
                        start = function(width) "ten"),
               "the value the pilot function returned")
  expect_error(seq_path(seq_aipe_cv, x, width = 0.3,
                        start = function(width) c(12, 9)),
               "returned more than one sample size (12, 9)", fixed = TRUE)
  # A string goes through the package's shared choice check, with what else
  # 'start' accepts added to its message.
  for (bad in list("soon", NA_character_, "")) {
    expect_error(
      seq_path(seq_aipe_cv, x, width = 0.3, start = bad),
      paste0("'start' must be one of \"pilot\", \"earliest\"; got .*\\. It ",
             "may also be NULL, a single whole number, or a pilot function\\.")
    )
  }
  expect_error(seq_path(seq_aipe_cv, x, width = 0.3, start = "soon"),
               "got \"soon\"", fixed = TRUE)
  # Anything else must be a whole number, and the message says what else
  # 'start' accepts: several strings are not a choice.
  expect_error(
    seq_path(seq_aipe_cv, x, width = 0.3, start = c("pilot", "earliest")),
    paste0("'start' must be a single whole number of at least 1; got a ",
           "character vector of length 2. It may also be NULL, a pilot ",
           "function, \"pilot\", or \"earliest\"."), fixed = TRUE
  )
  for (bad in list(2.5, 0, NA, TRUE, factor("pilot"), list(),
                   character(0))) {
    expect_error(
      seq_path(seq_aipe_cv, x, width = 0.3, start = bad),
      paste0("'start' must be a single whole number of at least 1; got .*\\. ",
             "It may also be NULL, a pilot function, \"pilot\", or ",
             "\"earliest\"\\.$")
    )
  }
  # A number that a pilot function returned is not 'start' as the user gave
  # it, and its message says nothing of the kind.
  returned <- tryCatch(
    seq_path(seq_aipe_cv, x, width = 0.3, start = function(width) "ten"),
    error = function(e) conditionMessage(e)
  )
  expect_match(returned, "the value the pilot function returned", fixed = TRUE)
  expect_false(grepl("It may also be", returned, fixed = TRUE))
  # As with the other choices in the package, a string may be abbreviated.
  expect_identical(seq_path(seq_aipe_cv, x, width = 0.3, start = "earl"),
                   seq_path(seq_aipe_cv, x, width = 0.3, start = "earliest"))
  expect_identical(seq_path(seq_aipe_cv, x, width = 0.3, start = "p"),
                   seq_path(seq_aipe_cv, x, width = 0.3))
  for (bad in list(0, 1.5, NA, "a", c(1, 2))) {
    expect_error(seq_path(seq_aipe_cv, x, width = 0.3, by = bad),
                 "'by' must be a single whole number of at least 1")
  }
})

test_that("an argument given twice is refused for what it is", {
  withr::local_seed(113)
  x <- stats::rgamma(60, shape = 16, rate = 20)
  g1 <- stats::rnorm(60, mean = 0.5)
  g2 <- stats::rnorm(60)
  # The same argument twice in '...': said plainly, and not described as an
  # argument of the pilot function (which the second copy used to be taken
  # for, because R matches only the first of two equal names).
  err <- tryCatch(seq_path(seq_aipe_cv, x, width = 0.1, width = 0.2),
                  error = function(e) conditionMessage(e))
  expect_match(err, "'width' is supplied more than once in '...'",
               fixed = TRUE)
  expect_false(grepl("_pilot", err, fixed = TRUE))
  expect_error(.path_refuse_repeated(list(a = 1, b = 2, a = 3, b = 4, 5, 6)),
               "'a', 'b' are supplied more than once in '...'", fixed = TRUE)
  expect_silent(.path_refuse_repeated(list(1, 2, a = 3)))
  expect_silent(.path_refuse_repeated(list()))

  # A data argument given again in '...' while `data` is a bare vector. The
  # vector goes to the first argument of the check function, so `x = ` in
  # '...' is the same argument twice; the message used to be about
  # conf_level, which the user never wrote.
  err <- tryCatch(seq_path(seq_aipe_cv, x, x = x, width = 0.1),
                  error = function(e) conditionMessage(e))
  expect_identical(err, paste0(
    "The argument x of the check function is supplied both in 'data' and ",
    "in '...'. Data belong in 'data'; everything else in '...'."
  ))
  # The same for an unnamed list, element by element, with the reminder of
  # how such a list is read.
  expect_error(seq_path(seq_aipe_smd, list(g1, g2), group_2 = g2,
                        width = 0.9),
               paste0("The argument group_2 of the check function is ",
                      "supplied both in 'data' and in '...'\\..*An unnamed ",
                      "list is read as one element for each data argument"))
  # An argument that is not one of those the data go to still reaches the
  # check, whatever its name.
  expect_s3_class(seq_path(seq_mr_cv, x, A = 400000, c = 75), "smsd_path")
  # A name that is exactly a later argument's goes to that argument (R
  # matches exact names first), even when it is the start of the first
  # one's: here `x` is the second argument, and the data go to `xs`.
  second <- function(xs, x, width) seq_aipe_cv(xs, width = width)
  expect_identical(attr(seq_path(second, x, x = 1, width = 0.1), "stop_n"),
                   attr(seq_path(seq_aipe_cv, x, width = 0.1), "stop_n"))
  # Partial matching, which R applies after exact matching, is a collision
  # too: `xs` would be matched by `xs` exactly and by `xs_data` partially.
  partial <- function(xs_data, width) seq_aipe_cv(xs_data, width = width)
  expect_error(seq_path(partial, x, xs = x, width = 0.1),
               "The argument xs_data of the check function is supplied both")
})

test_that("the Arguments line is short and breaks only between arguments", {
  # A number is written as it would be typed while that is short, and in
  # R's own notation beyond: fixed notation for 1e150 is 151 digits, and
  # for 2^70 it shows digits a double does not hold.
  expect_identical(.path_arg_text(400000), "400000")
  expect_identical(.path_arg_text(0), "0")
  expect_identical(.path_arg_text(0.00001), "0.00001")
  expect_identical(.path_arg_text(999999999999999), "999999999999999")
  expect_identical(.path_arg_text(1e150), "1e+150")
  expect_identical(.path_arg_text(2^70), "1.180592e+21")
  expect_identical(.path_arg_text(-1e15), "-1e+15")
  expect_identical(.path_arg_text(1e-7), "1e-07")
  # The same whatever the session's preference for fixed notation.
  withr::with_options(list(scipen = 100), {
    expect_identical(.path_arg_text(1e150), "1e+150")
    expect_identical(.path_arg_text(400000), "400000")
  })
  withr::with_options(list(scipen = -9), {
    expect_identical(.path_arg_text(400000), "400000")
    expect_identical(.path_arg_text(0.00001), "0.00001")
    expect_identical(.path_arg_text(c(1, -0.5, -0.5)), "c(1, -0.5, -0.5)")
  })
  expect_identical(.path_arg_text(c(2L, 5L)), "c(2, 5)")
  expect_identical(.path_arg_text(c(1, 1e20)), "c(1, 1e+20)")
  expect_identical(.path_arg_text(c(1, NA)), "c(1, NA)")
  # Anything that is not a single plain number.
  expect_identical(.path_arg_text(mean), "a function")
  expect_identical(.path_arg_text(c(1, -0.5, -0.5)), "c(1, -0.5, -0.5)")
  expect_identical(.path_arg_text("kendall"), "\"kendall\"")
  expect_identical(.path_arg_text(seq(0.1, 1.2, by = 0.1)),
                   "a numeric vector of length 12")
  expect_identical(.path_arg_text(matrix(0, 2, 2)),
                   "a matrix or array with dimensions 2 x 2")
  long <- .path_arg_text(strrep("abcdefghij", 5))
  expect_identical(nchar(long), 40L)
  expect_match(long, "\\.\\.\\.$")

  withr::local_seed(113)
  huge <- format(seq_path(seq_fwd_norm_mean, stats::rnorm(40, 10, 2) * 1e150,
                          d = 1e150, start = 10))
  expect_true("  Arguments       : d = 1e+150" %in% huge)
  expect_lte(max(nchar(huge)), 80L)

  # A vector argument of a contrast path.
  groups <- list(stats::rnorm(30, 1), stats::rnorm(30), stats::rnorm(30))
  contrast <- format(seq_path(seq_aipe_c, list(groups),
                              contrast_weights = c(1, -0.5, -0.5),
                              width = 2))
  expect_true("  Arguments       : contrast_weights = c(1, -0.5, -0.5)" %in%
                contrast)

  # Two arguments that do not fit on one line of 58 characters used to be
  # broken inside the second one ("f" on one line, "= a function" on the
  # next). They fit the line as it is now filled.
  toy <- function(x, w, f = NULL) {
    new_smsd_check(c("mean", "current_n"), c(mean(x), length(x)),
                   "toy procedure", is_satisfied = length(x) >= 6,
                   current_n = length(x), criterion_n = 6)
  }
  txt <- format(seq_path(toy, 1:10 + 0.5, w = I(seq(0.1, 1.2, by = 0.1)),
                         f = mean, start = 2))
  expect_true(paste("  Arguments       : w = a numeric vector of length 12,",
                    "f = a function") %in% txt)
  # When they do not fit, the break falls between arguments. At a width of
  # 60: the first piece and its comma end at column 58; the second would
  # reach 79, so it starts a line; the third would reach 66 on that line;
  # the fourth fits after the third (20 + 25 + 1 + 10 = 56).
  pieces <- c("first = a numeric vector of length 12", "second = a function",
              "third = c(1, -0.5, -0.5)", "fourth = 3")
  expect_identical(
    .path_args_lines(pieces, width = 60L),
    c("  Arguments       : first = a numeric vector of length 12,",
      "                    second = a function,",
      "                    third = c(1, -0.5, -0.5), fourth = 3")
  )
  # A piece longer than the line has a line to itself.
  expect_identical(
    .path_args_lines(c("a = 1", strrep("b", 70), "c = 3"), width = 60L),
    c("  Arguments       : a = 1,",
      paste0("                    ", strrep("b", 70), ","),
      "                    c = 3")
  )
  expect_identical(.path_args_pieces(list(width = 0.2, A = 5, 7),
                                     omit = "width"), c("A = 5", "7"))
  expect_identical(.path_args_text(list(width = 0.2), omit = "width"), "")
})

test_that("the notes of the check at the stage reported are printed", {
  caution <- paste("A caution that belongs to the stage at which the rule is",
                   "first met, long enough to need a second line when it is",
                   "printed under the estimate.")
  noted <- function(x, target, at = target) {
    n <- length(x)
    new_smsd_check(c("mean", "current_n"), c(mean(x), n), "toy procedure",
                   is_satisfied = n >= target, current_n = n,
                   criterion_n = target,
                   notes = if (n == at) caution else NULL)
  }
  txt <- format(seq_path(noted, 1:8, target = 5, start = 2))
  i <- which(startsWith(txt, "  Estimate there  : "))
  expect_length(i, 1L)
  note_lines <- txt[seq.int(i + 1L, which(txt == "")[2L] - 1L)]
  expect_gt(length(note_lines), 1L)
  expect_true(startsWith(note_lines[1L], "  Note            : A caution"))
  expect_true(all(startsWith(note_lines[-1L], strrep(" ", 20L))))
  expect_lte(max(nchar(note_lines)), 78L)
  # The lines hold the note, whole and once.
  expect_identical(
    paste(trimws(sub("^  Note            : ", "", note_lines)),
          collapse = " "),
    caution
  )
  # A note that belongs to another stage is not shown: the path reports the
  # stage at which the rule is first met.
  expect_false(any(grepl("Note   ", format(seq_path(noted, 1:8, target = 5,
                                                    at = 7, start = 2)),
                         fixed = TRUE)))

  # The package's own case: the slope of data on scales that make its
  # variance estimate positive but below n^-3. The rule is met at the pilot
  # through the bound, and the path now says so where the estimate is.
  withr::local_seed(113)
  x <- stats::rnorm(40)
  y <- 1e-6 * (2 * x + stats::rnorm(40))
  path <- seq_path(seq_aipe_reg_coef, list(x = x, y = y), width = 0.5)
  expect_identical(attr(path, "stop_n"), path$n[1L])
  final_notes <- attr(attr(path, "final_check"), "notes")
  expect_length(final_notes, 1L)
  expect_match(final_notes, "the bound was used in its place")
  txt <- format(path)
  expect_true(any(startsWith(txt, paste("  Note            : The variance",
                                        "estimate is positive but below"))))
  expect_true(any(grepl("the bound was used in its place", txt,
                        fixed = TRUE)))

  # A rule never met whose last stage has no usable variance: the line
  # under the table already says so for every such stage, so the standard
  # note is not printed a second time.
  never <- format(seq_path(toy_check, 1:6, target = 4, usable_from = 10))
  expect_false(any(startsWith(never, "  Note            : ")))
  expect_true(any(grepl("the variance estimate was zero, negative, or not",
                        never, fixed = TRUE)))
  # Any other note of such a stage is printed: one group of a standardized
  # mean difference with no variation, which the check names.
  flat <- format(seq_path(seq_aipe_smd, list(c(4, 5, 3, 4, 4, 5, 3, 4),
                                             rep(3, 8)), width = 1.5))
  expect_identical(sum(startsWith(flat, "  Note            : ")), 1L)
  expect_true(any(startsWith(flat, paste("  Note            : The",
                                         "observations in 'group_2' are"))))
})

test_that("the target and the confidence level are written as typed, whatever the display options", {
  # The table and the estimate follow digits and scipen, as print() does
  # everywhere; the constants the user chose do not, so that 95% never
  # reads 9.5e+01%.
  withr::local_seed(113)
  x <- stats::rgamma(60, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.2, conf_level = 0.9)
  withr::local_options(list(scipen = -5, digits = 3))
  txt <- format(path)
  expect_true("  Target width    : 0.2" %in% txt)
  expect_true("  Confidence level: 90%" %in% txt)
  expect_match(.path_subtitle(path, as.data.frame(path)),
               "^Target width 0\\.2, 90% confidence; rule ")
  # A pilot function that returns one size per group is named with its
  # sizes as numbers are typed.
  two <- function(...) c(12, 9)
  expect_error(seq_path(seq_aipe_cv, x, width = 0.2, start = two),
               "returned more than one sample size (12, 9)", fixed = TRUE)
})

test_that("with a decimal comma the limits are separated by a semicolon", {
  withr::local_seed(113)
  x <- stats::rgamma(60, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.2)
  i <- which(path$n == attr(path, "stop_n"))
  usual <- sprintf("  Estimate there  : cv = %s [%s, %s]",
                   format(path$estimate[i]), format(path$lower_limit[i]),
                   format(path$upper_limit[i]))
  expect_true(usual %in% format(path))
  withr::local_options(OutDec = ",")
  comma <- sprintf("  Estimate there  : cv = %s [%s; %s]",
                   format(path$estimate[i]), format(path$lower_limit[i]),
                   format(path$upper_limit[i]))
  expect_match(comma, "cv = 0,[0-9]+ \\[0,[0-9]+; 0,[0-9]+\\]$")
  expect_true(comma %in% format(path))
})

test_that("the line under the plot title names the last stage of a rule never met", {
  withr::local_seed(113)
  x <- stats::rgamma(30, shape = 9, rate = 3)
  never <- seq_path(seq_aipe_cv, x, width = 0.05, start = 10)
  expect_identical(attr(never, "stop_n"), NA_integer_)
  # The first stage is n = 10 and the last n = 30: the rule was not met by
  # the last.
  expect_identical(.path_subtitle(never, as.data.frame(never)),
                   "Target width 0.05, 95% confidence; rule not met by n = 30")
})

test_that("a function with a budget argument is refused through seq_path()", {
  # seq_budget() itself is refused by what it is and by the label of its
  # result; a function of the user's is known by the arguments of that
  # family, before any stage is run.
  withr::local_seed(113)
  x <- stats::rgamma(40, shape = 9, rate = 3)
  expect_error(
    seq_path(function(x, budget) seq_aipe_cv(x, 0.2), x, budget = 1),
    "does not apply to the budget constrained procedure"
  )
  expect_error(
    seq_path(function(x, cost_per_group) seq_aipe_cv(x, 0.2), x,
             cost_per_group = 1),
    "does not apply to the budget constrained procedure"
  )
})

test_that("a list of groups written after 'data' is taken for data", {
  withr::local_seed(113)
  x <- stats::rgamma(40, shape = 9, rate = 3)
  g1 <- stats::rgamma(40, shape = 9, rate = 3)
  g2 <- stats::rgamma(40, shape = 9, rate = 3)
  pooled <- function(x, groups, width) {
    seq_aipe_cv(c(x, unlist(groups)), width = width)
  }
  expect_error(seq_path(pooled, x, groups = list(g1, g2), width = 0.3),
               paste0("The argument 'groups' holds 40 observations, and the ",
                      "check function also runs when it is cut along with ",
                      "the data, so it looks like data"))
  # Wrapped in I() it is passed whole, on the user's word.
  expect_s3_class(seq_path(pooled, x, groups = I(list(g1, g2)), width = 0.3),
                  "smsd_path")
  # What counts as a number of observations.
  expect_identical(.path_extra_count(data.frame(a = 1:7, b = 7:1)), 7L)
  expect_identical(.path_extra_count(matrix(0, 5, 2)), 5L)
  expect_identical(.path_extra_count(1:9), 9L)
  expect_identical(.path_extra_count(list(1:4, 5:8)), 4L)
  expect_identical(.path_extra_count(list(1:4, matrix(0, 4, 2))), 4L)
  expect_identical(.path_extra_count(list()), NA_integer_)
  expect_identical(.path_extra_count(list(1:4, 1:5)), NA_integer_)
  expect_identical(.path_extra_count(list(1:4, list(1:4))), NA_integer_)
  expect_identical(.path_extra_count(NULL), NA_integer_)
  expect_identical(.path_extra_count(mean), NA_integer_)
  expect_identical(.path_extra_count(array(0, c(2, 2, 2))), NA_integer_)
})

test_that("a stage whose rows lack the estimate term is read through summary()", {
  # The name of the estimate row is taken from the check on all the data.
  # A check that names the row differently at some stages (here "average"
  # at odd sample sizes) is still read: the estimate is the row summary()
  # takes for it.
  renaming <- function(x, target) {
    n <- length(x)
    new_smsd_check(c(if (n %% 2L == 0L) "mean" else "average", "current_n"),
                   c(mean(x), n), "toy procedure", is_satisfied = n >= target,
                   current_n = n, criterion_n = target)
  }
  x <- c(2, 4, 9, 1, 6, 8)
  path <- seq_path(renaming, x, target = 5, start = 2)
  expect_identical(attr(path, "estimate_term"), "mean")
  expect_equal(path$estimate, cumsum(x)[2:6] / 2:6)
  expect_identical(path$is_satisfied, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_identical(path$criterion_n, rep(5, 5))
  expect_identical(path$current_n, as.numeric(2:6))
  expect_identical(attr(path, "stop_n"), 5L)
})

test_that("print shows the header, the ends of the table, and the reminder", {
  withr::local_seed(113)
  x <- rgamma(60, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.2)
  stop_n <- attr(path, "stop_n")
  at_stop <- seq_aipe_cv(x[seq_len(stop_n)], width = 0.2)
  txt <- format(path)
  expect_identical(txt[1], "Sequential path: coefficient of variation")
  expect_true(sprintf("  Stages checked  : %d (n = %d to 60)", nrow(path),
                      path$n[1]) %in% txt)
  expect_true("  First check at  : the pilot sample size" %in% txt)
  expect_true("  Target width    : 0.2" %in% txt)
  expect_true("  Confidence level: 95%" %in% txt)
  expect_true(paste("  Rule first met  : n =", stop_n) %in% txt)
  at <- function(term) at_stop$value[at_stop$term == term]
  est_line <- sprintf("  Estimate there  : cv = %s [%s, %s]",
                      format(at("cv")), format(at("lower_limit")),
                      format(at("upper_limit")))
  expect_true(est_line %in% txt)
  # Three rows from each end, with a break between.
  expect_true("  ..." %in% txt)
  expect_identical(sum(grepl("^ +[0-9]+ +[0-9.]+ ", txt)), 6L)
  expect_true(any(grepl(paste("The estimate to report is the one at n =",
                              stop_n), txt, fixed = TRUE)))
  # Columns with nothing to show are left out of the printed table only.
  header <- txt[grepl("criterion_n", txt)]
  expect_false(grepl("variance_unusable|current_n", header))
  expect_true(grepl("lower_limit", header))

  expect_identical(capture.output(res <- print(path)), txt)
  expect_identical(res, path)
  expect_identical(sum(grepl("^ +[0-9]+ +[0-9.]+ ",
                             format(path, rows = 1))), 2L)
  expect_error(format(path, rows = 0), "'rows' must be")
  # digits is checked before it reaches base::format().
  for (bad in list(0, "a", NA, 2.5, c(3, 4))) {
    expect_error(format(path, digits = bad),
                 "'digits' must be a single whole number of at least 1")
    expect_error(print(path, digits = bad), "'digits' must be")
  }
  expect_error(format(path, digits = 23),
               "'digits' must be a single whole number from 1 to 22; got 23")
  expect_error(format(path, digits = 22), NA)

  # digits changes the display, never the stored values.
  short <- format(path, digits = 3)
  expect_true(any(grepl(sprintf("cv = %s ", format(at("cv"), digits = 3)),
                        short, fixed = TRUE)))
  expect_identical(path$estimate[path$n == stop_n], at("cv"))
})

test_that("print adapts to risk procedures, batches, and a rule never met", {
  withr::local_seed(113)
  g1 <- rnorm(40, 0.5)
  g2 <- rnorm(40)
  risk <- format(seq_path(toy_two_group, list(g1, g2), target = 30, by = 2,
                          start = 2))
  expect_true("  Arguments       : target = 30" %in% risk)
  expect_true("  Stages checked  : 20 (n = 2 to 40 by 2)" %in% risk)
  expect_true("  First check at  : the sample size supplied as 'start'" %in%
                risk)
  expect_true("  Rule first met  : n = 16" %in% risk)
  expect_false(any(grepl("Target width|Confidence level", risk)))
  header <- risk[grepl("criterion_n", risk)]
  expect_true(grepl("current_n", header))
  expect_false(grepl("lower_limit", header))
  # The loss constants of a package procedure are shown the same way.
  expect_true("  Arguments       : A = 2000, c = 1" %in%
                format(seq_path(seq_mr_smd, list(g1, g2), A = 2000, c = 1)))

  never <- format(seq_path(seq_aipe_smd, list(g1, g2), width = 0.05,
                           start = "earliest"))
  expect_true("  First check at  : the smallest sample size the check accepts"
              %in% never)
  expect_true("  Rule first met  : at no stage checked" %in% never)
  expect_true(any(grepl("^  Latest estimate : smd = ", never)))
  expect_true(any(grepl("was not met at any stage checked", never)))

  # A short path is printed in full, with no break.
  short <- format(seq_path(toy_check, c(1, 2, 3, 4, 5), target = 4))
  expect_false("  ..." %in% short)
  expect_identical(sum(grepl("^ +[0-9]+ +[0-9.]+ ", short)), 4L)
  # The documented boundary: a path of no more than 2 * rows + 1 stages is
  # shown in full. With rows = 3, seven stages (n = 2 to 8) are all shown
  # and eight (n = 2 to 9) are cut to the first and last three.
  stage_rows <- function(lines) sum(grepl("^ +[0-9]+ +[0-9.]+ ", lines))
  seven <- format(seq_path(toy_check, 1:8, target = 4))
  expect_false("  ..." %in% seven)
  expect_identical(stage_rows(seven), 7L)
  eight <- format(seq_path(toy_check, 1:9, target = 4))
  expect_true("  ..." %in% eight)
  expect_identical(stage_rows(eight), 6L)
  # With rows = 1 the boundary is three stages.
  expect_false("  ..." %in% format(seq_path(toy_check, 1:4, target = 4),
                                   rows = 1))
  expect_true("  ..." %in% format(seq_path(toy_check, 1:5, target = 4),
                                  rows = 1))
  # A loss constant is shown as it would be typed, not as 4e+05.
  expect_true("  Arguments       : A = 400000, c = 75" %in%
                format(seq_path(seq_mr_cv, stats::rgamma(40, 9, 3),
                                A = 4e5, c = 75)))
  # Met at the last stage: nothing comes later to warn about.
  last <- format(seq_path(toy_check, c(1, 2, 3, 4, 5), target = 5))
  expect_identical(last[length(last)],
                   "  The estimate to report is the one at n = 5.")
})

test_that("summary is the stopping stage in one row", {
  withr::local_seed(113)
  x <- rgamma(60, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.2)
  stop_n <- attr(path, "stop_n")
  at_stop <- seq_aipe_cv(x[seq_len(stop_n)], width = 0.2)
  s <- summary(path)
  expect_s3_class(s, "data.frame", exact = TRUE)
  expect_identical(nrow(s), 1L)
  expect_identical(names(s),
                   c("method", "estimate_term", "stop_n", "n", "estimate",
                     "lower_limit", "upper_limit", "current_n",
                     "criterion_n", "is_satisfied", "variance_unusable",
                     "conf_level", "width", "stages"))
  expect_identical(s$method, "coefficient of variation")
  expect_identical(s$estimate_term, "cv")
  expect_identical(s$stop_n, stop_n)
  expect_identical(s$n, stop_n)
  expect_identical(s$estimate, at_stop$value[at_stop$term == "cv"])
  expect_identical(s$lower_limit,
                   at_stop$value[at_stop$term == "lower_limit"])
  expect_identical(s$upper_limit,
                   at_stop$value[at_stop$term == "upper_limit"])
  expect_equal(s$criterion_n, attr(at_stop, "criterion_n"))
  expect_true(s$is_satisfied)
  expect_identical(s$conf_level, 0.95)
  expect_identical(s$width, 0.2)
  expect_identical(s$stages, nrow(path))

  # Never met: the last stage, with stop_n NA; no interval for a risk rule.
  never <- summary(seq_path(toy_check, c(3, 1, 4, 1, 5), target = 99))
  expect_identical(never$stop_n, NA_integer_)
  expect_identical(never$n, 5L)
  expect_false(never$is_satisfied)
  expect_identical(never$estimate, 14 / 5)
  expect_identical(never$criterion_n, 99)
  expect_identical(never$stages, 4L)
  expect_false(never$variance_unusable)
  expect_true(is.na(never$conf_level) && is.na(never$width))
  expect_true(is.na(never$lower_limit) && is.na(never$upper_limit))
})

test_that("as.data.frame and [ return plain data for anything partial", {
  withr::local_seed(113)
  x <- rgamma(40, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.2, start = "earliest")
  plain <- as.data.frame(path)
  expect_identical(class(plain), "data.frame")
  expect_identical(sort(names(attributes(plain))),
                   c("class", "names", "row.names"))
  expect_identical(dim(plain), c(length(path$n), 8L))
  expect_identical(plain$estimate, path$estimate)

  part <- path[1:5, ]
  expect_identical(class(part), "data.frame")
  expect_null(attr(part, "stop_n"))
  expect_null(attr(part, "final_check"))
  expect_identical(class(head(path, 3)), "data.frame")
  expect_identical(class(path[, c("n", "estimate")]), "data.frame")
  expect_identical(path[path$n == 10L, "estimate"],
                   path$estimate[match(10L, path$n)])
  expect_type(path[["n"]], "integer")
  # The whole path keeps its class and its attributes.
  whole <- path[, ]
  expect_s3_class(whole, "smsd_path")
  expect_identical(attr(whole, "stop_n"), attr(path, "stop_n"))
})

test_that("a subset that keeps the whole path is still the path", {
  withr::local_seed(113)
  x <- rgamma(40, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.2)
  # Forms of subsetting for which the data frame method keeps the class but
  # drops the attributes.
  expect_identical(path[names(path)], path)
  expect_identical(path[1:8], path)
  expect_identical(path[, 1:8], path)
  expect_identical(path[seq_len(nrow(path)), ], path)
  expect_identical(subset(path, n > 0), path)
  expect_identical(format(path[names(path)]), format(path))
  # Anything less is plain, attributes and all.
  for (part in list(path[1:7], path[, 1:7], subset(path, n > 25),
                    path[rev(seq_len(nrow(path))), ],
                    path[c("estimate", "n")])) {
    expect_identical(class(part), "data.frame")
    expect_identical(sort(names(attributes(part))),
                     c("class", "names", "row.names"))
  }
})

test_that("an object that is no longer a path is refused by every method", {
  withr::local_seed(113)
  x <- rgamma(40, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.2)
  expect_null(.path_problem(path))
  # What a data frame verb of another package can leave: the class and the
  # attributes of the path on other rows.
  dressed <- function(rows) {
    out <- as.data.frame(path)[rows, , drop = FALSE]
    for (a in .path_attrs) attr(out, a) <- attr(path, a)
    class(out) <- class(path)
    out
  }
  stop_n <- attr(path, "stop_n")
  expect_lt(stop_n, 40L)
  later <- dressed(path$n > stop_n)          # the stopping stage is gone
  earlier <- dressed(path$n < 40L)           # the tail is gone
  doubled <- rbind(path, path)
  expect_s3_class(doubled, "smsd_path")
  stripped <- path
  attr(stripped, "stop_n") <- NULL
  moved <- path
  attr(moved, "stop_n") <- 40L               # a later, more favorable stage
  no_columns <- structure(list(n = 1:3), class = c("smsd_path", "data.frame"),
                          row.names = 1:3)
  cases <- list(
    list(later,
         "its rows are not the stages seq_path\\(\\) checked, n = 20 to 40"),
    list(earlier, "its rows are not the stages"),
    list(doubled, "its rows are not the stages"),
    list(stripped, "the attribute stop_n is missing"),
    list(moved, "'stop_n' is not the first stage at which the rule is met"),
    list(no_columns, "columns are missing")
  )
  withr::local_pdf(NULL)
  for (case in cases) {
    broken <- case[[1]]
    expect_match(.path_problem(broken), case[[2]])
    for (method in list(format, print, summary, plot, seq_report)) {
      expect_error(method(broken),
                   "has class 'smsd_path' but is no longer a path")
      expect_error(method(broken), case[[2]])
      expect_error(method(broken),
                   "Use as.data.frame\\(\\) for the plain table")
    }
  }
  # as.data.frame() is the way out, and gives what the object holds.
  expect_identical(class(as.data.frame(later)), "data.frame")
  expect_identical(as.data.frame(later)$n, path$n[path$n > stop_n])
  expect_identical(nrow(as.data.frame(doubled)), 2L * nrow(path))
  expect_identical(as.data.frame(no_columns), data.frame(n = 1:3))
  # A path with batches and a start of its own is a path.
  expect_null(.path_problem(seq_path(seq_aipe_cv, x, width = 0.2, start = 7,
                                     by = 4)))
  expect_null(.path_problem(seq_path(seq_aipe_cv, x, width = 0.01,
                                     start = 30)))
})

test_that("plot draws every kind of path and leaves par() as it found it", {
  withr::local_seed(113)
  x <- rgamma(60, shape = 9, rate = 3)
  g1 <- rnorm(40, 0.5)
  g2 <- rnorm(40)
  interval <- seq_path(seq_aipe_cv, x, width = 0.2)
  risk <- seq_path(seq_mr_cv, x, A = 400000, c = 75)
  never <- seq_path(seq_aipe_smd, list(g1, g2), width = 0.05,
                    start = "earliest")
  batches <- seq_path(seq_aipe_smd, list(g1, g2), width = 0.9, by = 6)
  one_stage <- seq_path(seq_aipe_cv, x, width = 0.2, start = 60)
  expect_false(is.na(attr(interval, "stop_n")))
  expect_true(all(is.na(risk$lower_limit)))
  expect_true(is.na(attr(never, "stop_n")))

  withr::local_pdf(NULL)
  graphics::par(mfrow = c(1, 2), mar = c(3, 3, 2, 1), cex = 0.9)
  graphics::plot(1)
  before <- graphics::par(no.readonly = TRUE)
  # What any plot changes: the coordinates of the last panel drawn and the
  # position in the layout (the position is pinned in the next test).
  moves <- c("usr", "xaxp", "yaxp", "fig", "mfg", "plt", "pin")
  same <- setdiff(names(before), moves)
  for (path in list(interval, risk, never, batches, one_stage)) {
    expect_error(res <- plot(path), NA)
    expect_identical(res, path)
    expect_identical(graphics::par(no.readonly = TRUE)[same], before[same])
    for (panel in c("estimate", "criterion")) {
      expect_error(plot(path, which = panel), NA)
      expect_identical(graphics::par(no.readonly = TRUE)[same], before[same])
    }
  }
  expect_invisible(plot(interval))
  expect_error(plot(interval, legend = FALSE, main = "A title",
                    xlab = "Participants", ylab = c("CV", "N"),
                    col = "black", lwd = 1, lty = 2, las = 0), NA)
  expect_error(plot(interval, main = ""), NA)
  expect_error(plot(interval, which = "criterion", log = "y",
                    ylim = c(1, 1000), cex.main = 0.8), NA)
  expect_warning(plot(interval, which = "criterion", log = "xy"), NA)
  expect_warning(plot(risk), NA)
  expect_warning(plot(never, which = "estimate", ylim = c(-2, 2)), NA)
  expect_identical(graphics::par(no.readonly = TRUE)[same], before[same])
})

test_that("a single panel takes one cell of the user's layout", {
  withr::local_seed(113)
  x <- rgamma(60, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.2)
  withr::local_pdf(NULL)
  # In a 2 by 2 layout filled by rows, three plots in a row occupy cells
  # (1, 1), (1, 2), and (2, 1). A panel that reset the layout would start
  # again at the first cell.
  graphics::par(mfrow = c(2, 2))
  graphics::plot(1)
  expect_identical(graphics::par("mfg"), c(1L, 1L, 2L, 2L))
  plot(path, which = "estimate")
  expect_identical(graphics::par("mfg"), c(1L, 2L, 2L, 2L))
  plot(path, which = "criterion")
  expect_identical(graphics::par("mfg"), c(2L, 1L, 2L, 2L))
  expect_identical(graphics::par("mfrow"), c(2L, 2L))
  # Both panels are stacked in a layout of their own, and the user's layout
  # is put back afterward.
  plot(path)
  expect_identical(graphics::par("mfrow"), c(2L, 2L))
})

test_that("each panel shows the series and the stopping point it should", {
  withr::local_seed(113)
  x <- rgamma(60, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.2)
  stop_n <- attr(path, "stop_n")
  # Independent values: direct calls of the check on the first n.
  direct <- lapply(path$n, function(n) seq_aipe_cv(x[seq_len(n)], width = 0.2))
  value <- function(term) {
    vapply(direct, function(d) d$value[d$term == term], numeric(1))
  }
  first_met <- path$n[which(vapply(direct, attr, logical(1),
                                   which = "is_satisfied"))[1]]
  expect_identical(stop_n, first_met)
  expect_lt(stop_n, 60L)
  i_stop <- match(first_met, path$n)
  k <- nrow(path)
  tab <- as.data.frame(path)

  est <- .path_panel_data(tab, stop_n, "estimate")
  expect_identical(est$n, path$n)
  expect_identical(est$y, value("cv"))
  expect_identical(est$lower, value("lower_limit"))
  expect_identical(est$upper, value("upper_limit"))
  expect_null(est$criterion)
  expect_true(est$has_band)
  # The path proper runs up to the stopping stage; the gray part starts
  # there (so that the lines join) and runs to the last stage.
  expect_identical(est$before, seq_len(i_stop))
  expect_identical(est$after, i_stop:k)
  expect_identical(est$marker,
                   c(n = first_met, y = value("cv")[i_stop]))
  expect_identical(est$y_range,
                   range(c(value("lower_limit"), value("upper_limit"))))

  crit <- .path_panel_data(tab, stop_n, "criterion")
  current <- vapply(direct, attr, numeric(1), which = "current_n")
  criterion <- vapply(direct, attr, numeric(1), which = "criterion_n")
  expect_equal(crit$y, current)
  expect_equal(crit$criterion, criterion)
  expect_null(crit$lower)
  expect_false(crit$has_band)
  expect_equal(crit$marker, c(n = first_met, y = current[i_stop]))
  expect_equal(crit$y_range, range(c(current, criterion)))

  # What reaches the device: with the default axis style the vertical axis
  # spans the range of the data plus 4 percent at each end (no legend, so
  # no room is added for one).
  withr::local_pdf(NULL)
  plot(path, which = "criterion", legend = FALSE)
  span <- range(c(current, criterion))
  expect_equal(graphics::par("usr")[3:4], span + c(-0.04, 0.04) * diff(span))
  expect_equal(graphics::par("usr")[1:2], c(20, 60) + c(-0.04, 0.04) * 40)
  plot(path, which = "estimate", legend = FALSE)
  span <- range(c(value("lower_limit"), value("upper_limit")))
  expect_equal(graphics::par("usr")[3:4], span + c(-0.04, 0.04) * diff(span))
  # With a legend, room is left above the data: the axis reaches higher.
  plot(path, which = "estimate")
  expect_gt(graphics::par("usr")[4], span[2] + 0.2 * diff(span))
  expect_equal(graphics::par("usr")[3], span[1] - 0.04 * 1.28 * diff(span))

  # A rule never met: every stage is the path proper, and nothing is marked.
  never <- seq_path(seq_aipe_cv, x, width = 0.02, start = 50)
  shown <- .path_panel_data(as.data.frame(never), attr(never, "stop_n"),
                            "estimate")
  expect_identical(shown$before, 1:11)
  expect_identical(shown$after, integer(0))
  expect_null(shown$marker)
  expect_identical(shown$i_stop, NA_integer_)
  # Met at the last stage: nothing comes after.
  last <- seq_path(toy_check, c(1, 2, 3, 4, 5), target = 5)
  shown <- .path_panel_data(as.data.frame(last), attr(last, "stop_n"),
                            "estimate")
  expect_identical(shown$before, 1:4)
  expect_identical(shown$after, integer(0))
  expect_identical(shown$marker, c(n = 5, y = 3))
  # Values that are not finite are not drawn, and a flat series still gets
  # an axis of positive length.
  flat <- data.frame(n = 2:4, estimate = c(Inf, 7, 7), lower_limit = NA_real_,
                     upper_limit = NA_real_, current_n = 2:4,
                     criterion_n = NA_real_, is_satisfied = FALSE,
                     variance_unusable = NA)
  shown <- .path_panel_data(flat, NA_integer_, "estimate")
  expect_identical(shown$y, c(NA, 7, 7))
  expect_false(shown$has_band)
  expect_identical(shown$y_range, c(3.5, 10.5))
})

test_that("the stopping line ends below the legend, not through it", {
  # Coordinates as par("usr") gives them. A legend box with its left edge at
  # 80, its top at 1, and a height of 0.2 has its bottom at 0.8.
  box <- list(w = 40, h = 0.2, left = 80, top = 1)
  # Under the box the line stops at the bottom of the box.
  expect_identical(.path_marker_top(limit = 1, point = 0.5, at = 100,
                                    box = box), 0.8)
  # Left of the box it runs to the limit.
  expect_identical(.path_marker_top(limit = 1, point = 0.5, at = 60,
                                    box = box), 1)
  # A limit below the box (room was left above the data) is kept.
  expect_identical(.path_marker_top(limit = 0.7, point = 0.5, at = 100,
                                    box = box), 0.7)
  # The line always reaches the point it marks, even inside the box.
  expect_identical(.path_marker_top(limit = 1, point = 0.9, at = 100,
                                    box = box), 0.9)
  # No legend: the limit; a point that cannot be drawn does not matter.
  expect_identical(.path_marker_top(limit = 1, point = 0.5, at = 100), 1)
  expect_identical(.path_marker_top(limit = 1, point = NA_real_, at = 100,
                                    box = box), 0.8)
})

test_that("plot handles a procedure with no criterion and non-finite stages", {
  x <- c(3, 1, 4, 1, 5, 9, 2, 6)
  # A procedure that reports no criterion sample size (as the bounded length
  # interval for a median does) and whose estimate is not finite at first.
  no_criterion <- function(x, target) {
    n <- length(x)
    est <- if (n < 4L) NaN else stats::median(x)
    out <- new_smsd_check(
      term = c("lower_limit", "median", "upper_limit"),
      value = c(0, 0, 0),
      method = "no criterion", is_satisfied = n >= target,
      current_n = n, criterion_n = NA_real_, conf_level = 0.9, width = 0.5
    )
    # Set after construction: the constructor refuses NaN.
    out$value <- c(est - 1 / n, est, est + 1 / n)
    out
  }
  path <- seq_path(no_criterion, x, target = 6)
  expect_true(all(is.na(path$criterion_n)))
  expect_true(all(is.nan(path$estimate[path$n < 4])))
  tab <- as.data.frame(path)
  # By default the criterion panel is left out; asked for by name it is an
  # error.
  expect_identical(.path_plot_panels(tab, c("estimate", "criterion"),
                                     default_which = TRUE), "estimate")
  expect_error(.path_plot_panels(tab, c("estimate", "criterion"),
                                 default_which = FALSE),
               "reports no criterion sample size")
  expect_identical(
    .path_plot_panels(data.frame(estimate = 1, criterion_n = 5),
                      c("criterion", "estimate"), default_which = TRUE),
    c("criterion", "estimate")
  )
  withr::local_pdf(NULL)
  # One panel only is drawn: in a layout of two cells the default plot takes
  # the next cell and leaves the layout alone.
  graphics::par(mfrow = c(2, 2))
  graphics::plot(1)
  expect_error(plot(path), NA)
  expect_identical(graphics::par("mfg"), c(1L, 2L, 2L, 2L))
  expect_warning(plot(path), NA)
  expect_identical(graphics::par("mfg"), c(2L, 1L, 2L, 2L))
  expect_error(plot(path, which = "criterion"),
               "reports no criterion sample size")
  # 'which' goes through the package's shared choice check, name by name.
  expect_error(
    plot(path, which = "nonsense"),
    paste0("'which' must be one of \"estimate\", \"criterion\"; got ",
           "\"nonsense\". Give both, or leave 'which' out, to draw both ",
           "panels."), fixed = TRUE
  )
  expect_error(plot(path, which = c("estimate", "zzz")), "got \"zzz\"",
               fixed = TRUE)
  for (bad in list(2, NA_character_, character(0), NULL, TRUE)) {
    expect_error(plot(path, which = bad), "'which' must be one of")
  }
  # One panel or both, in the order given; a name may be abbreviated.
  expect_identical(.path_which(c("estimate", "criterion")),
                   c("estimate", "criterion"))
  expect_identical(.path_which(c("criterion", "estimate")),
                   c("criterion", "estimate"))
  expect_identical(.path_which("criterion"), "criterion")
  expect_identical(.path_which("est"), "estimate")
  expect_identical(.path_which(c("crit", "criterion", "crit")), "criterion")
  expect_error(plot(path, legend = NA), "'legend' must be TRUE or FALSE")
  never_finite <- seq_path(function(x) {
    out <- new_smsd_check("mean", NA_real_, "never finite",
                          is_satisfied = FALSE, current_n = length(x),
                          criterion_n = 10)
    out$value <- NaN   # set after construction: the constructor refuses NaN
    out
  }, x)
  expect_error(plot(never_finite), "not finite at any stage")
  expect_error(plot(never_finite, which = "criterion"), NA)
})

test_that("titles and axis labels are checked, and may be expressions", {
  withr::local_seed(113)
  x <- rgamma(60, shape = 9, rate = 3)
  path <- seq_path(seq_aipe_cv, x, width = 0.2)
  withr::local_pdf(NULL)
  expect_error(plot(path, main = c("a", "b")),
               "'main' must be one string or expression; got a character")
  expect_error(plot(path, main = 3), "'main' must be one string or expression")
  expect_error(plot(path, main = NA_character_), "'main' must be")
  expect_error(plot(path, xlab = c("a", "b")),
               "'xlab' must be one string or expression")
  expect_error(plot(path, ylab = 1:3),
               "'ylab' must be one string or expression, or one for each panel")
  expect_error(plot(path, ylab = c("a", "b", "c")), "panel drawn \\(2\\)")
  expect_error(plot(path, which = "estimate", ylab = c("a", "b")),
               "'ylab' must be one string or expression; got")
  expect_error(plot(path, ylab = list("a", 2)), "'ylab' must be")
  # One label for each panel, in any of the forms; an expression is plotted
  # as mathematics.
  for (ylab in list("CV", c("CV", "N"), expression(hat(kappa)),
                    expression(hat(kappa), italic(N)), quote(hat(kappa)),
                    list(quote(hat(kappa)), "N"))) {
    expect_warning(plot(path, ylab = ylab), NA)
  }
  expect_warning(plot(path, which = "estimate", ylab = expression(hat(kappa)),
                      xlab = quote(italic(n)),
                      main = expression(paste("Path for ", kappa))), NA)
  # The labels as the panels receive them: a list with one string or one
  # expression for each panel.
  expect_identical(.path_plot_label("CV", "ylab", 2L), list("CV", "CV"))
  expect_identical(.path_plot_label(c("CV", "N"), "ylab", 2L),
                   list("CV", "N"))
  expect_identical(.path_plot_label(quote(hat(kappa)), "ylab", 1L),
                   list(expression(hat(kappa))))
  expect_identical(.path_plot_label(expression(a, b), "ylab", 2L),
                   list(expression(a), expression(b)))
})

test_that("operating characteristic: stop_n of an exactly known rule", {
  # Inverse sampling: stop at the first n with r successes among Bernoulli(p)
  # trials. The stopping size is negative binomial, with mean r / p and
  # variance r (1 - p) / p^2 exactly, so the first passage that seq_path()
  # reports can be checked against a closed form.
  r <- 4
  p <- 0.5
  inverse_sampling <- function(x, r) {
    n <- length(x)
    new_smsd_check(
      term = c("proportion", "successes", "current_n"),
      value = c(mean(x), sum(x), n), method = "inverse sampling",
      is_satisfied = sum(x) >= r, current_n = n, criterion_n = NA_real_
    )
  }
  reps <- 200
  stop_n <- withr::with_seed(113, vapply(seq_len(reps), function(i) {
    x <- stats::rbinom(32, size = 1, prob = p)
    attr(seq_path(inverse_sampling, x, r = r, start = r), "stop_n")
  }, integer(1)))
  # P(N > 32) = P(fewer than 4 successes in 32 trials) = 1.3e-6, so no path
  # ends with the rule unmet.
  expect_equal(stats::pbinom(r - 1, 32, p), 1.28e-6, tolerance = 0.01)
  expect_false(anyNA(stop_n))
  # Tolerance: four Monte Carlo standard errors of the mean stopping size,
  # sqrt(r (1 - p) / p^2 / reps) = 0.2 each.
  se <- sqrt(r * (1 - p) / p^2 / reps)
  expect_equal(se, 0.2)
  expect_lt(abs(mean(stop_n) - r / p), 4 * se)
  # The spread is the negative binomial's too: variance 8. The standard
  # error of the sample variance of 200 such stopping sizes is 1.07 (found
  # by simulating the negative binomial), so 4.5 is more than four of them.
  expect_lt(abs(stats::var(stop_n) - r * (1 - p) / p^2), 4.5)
  # The rule cannot be met before r observations.
  expect_true(all(stop_n >= r))
})

test_that("operating characteristic: the full procedure, stage by stage", {
  # A small Monte Carlo of the sequential AIPE procedure for the coefficient
  # of variation, run through seq_path() and through the loop a user would
  # write. Normal data with mean 10 and SD 2; 95% interval of full width
  # 0.12; the path starts at the pilot, as the loop does.
  width <- 0.12
  pilot <- seq_aipe_cv_pilot(width = width)
  reps <- 25
  met <- 0L
  withr::with_seed(113, {
    for (i in seq_len(reps)) {
      x <- stats::rnorm(70, mean = 10, sd = 2)
      n <- pilot
      repeat {
        out <- seq_aipe_cv(x[seq_len(n)], width = width)
        if (attr(out, "is_satisfied") || n == 70L) break
        n <- n + 1L
      }
      s <- summary(seq_path(seq_aipe_cv, x, width = width))
      if (attr(out, "is_satisfied")) {
        met <- met + 1L
        # The path stops exactly where the loop stops, and the interval it
        # reports is no wider than the target: the rule
        # n >= (2 z / w)^2 (xi^2 + 1 / n) implies 2 z sqrt(xi^2 / n) < w.
        expect_identical(s$stop_n, n)
        expect_lt(s$upper_limit - s$lower_limit, width)
      } else {
        expect_identical(s$stop_n, NA_integer_)
      }
    }
  })
  # Nearly every replication stops within 70 observations: in 2,000
  # replications of this loop with seed 1 and no cap, the stopping size
  # averaged 45.5 with a standard deviation of 5.7, and 7 went beyond 70.
  # Three or more of 25 beyond 70 then has a chance of about 1e-4.
  expect_gte(met, reps - 2L)
})
