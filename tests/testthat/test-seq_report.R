# seq_report(): the paragraph that reports a sequential check or path.
#
# The objects are built with the constructor from known numbers, so every
# number expected in the text is written down here, not computed by the
# function under test.

known_check <- function(is_satisfied = TRUE, ...) {
  new_smsd_check(
    term = c("lower_limit", "smd", "upper_limit", "xi_squared",
             "criterion_n", "current_n"),
    value = c(0.20456, 0.50044, 0.79632, 2.03, 88, 88),
    method = "standardized mean difference",
    is_satisfied = is_satisfied, current_n = 88,
    criterion_n = if (is_satisfied) 88 else 121,
    conf_level = 0.95, width = 0.6, ...
  )
}

known_risk <- function(...) {
  new_smsd_check(
    term = c("cv", "risk", "xi_squared", "optimal_n", "current_n"),
    value = c(0.21987, 5964.345, 0.031, 35, 35),
    method = "minimum risk point estimation (coefficient of variation)",
    is_satisfied = TRUE, current_n = 35L, criterion_n = 35L,
    extra = list(loss_A = 400000, cost_c = 75, gamma = 0.49), ...
  )
}

test_that("the paragraph holds the right numbers for a known check", {
  rep <- seq_report(known_check())
  expect_s3_class(rep, c("smsd_report", "character"), exact = TRUE)
  expect_length(rep, 1L)
  txt <- as.character(rep)
  expect_identical(
    txt,
    paste(
      "Sample size was determined with the sequential accuracy in parameter",
      "estimation procedure for the standardized mean difference (Kelley,",
      "Darku, & Chattopadhyay, 2018), with a target full width of 0.6 for",
      "the 95% confidence interval. The stopping rule was met at n = 88 per",
      "group. At that sample size the estimate of the standardized mean",
      "difference was 0.500, 95% CI [0.205, 0.796]."
    )
  )
  expect_identical(attr(rep, "method"), "standardized mean difference")
  expect_identical(
    attr(rep, "reference"),
    paste(
      "Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in",
      "parameter estimation for a general class of effect sizes: A",
      "sequential approach. Psychological Methods, 23(2), 226\u2013243.",
      "doi:10.1037/met0000127"
    )
  )
})

test_that("a rule not yet met is reported as an interim result", {
  txt <- as.character(seq_report(known_check(is_satisfied = FALSE)))
  expect_match(txt, "Sample size is being determined with", fixed = TRUE)
  expect_match(txt, paste("At n = 88 per group the stopping rule has not yet",
                          "been met (the criterion sample size at this stage",
                          "is 121), so sampling continues."), fixed = TRUE)
  expect_match(txt, paste("The interim estimate of the standardized mean",
                          "difference is 0.500, 95% CI [0.205, 0.796]."),
               fixed = TRUE)
  expect_false(grepl("was met", txt, fixed = TRUE))
  expect_false(grepl("was determined", txt, fixed = TRUE))
})

test_that("digits changes the text only", {
  chk <- known_check()
  expect_match(as.character(seq_report(chk, digits = 2)),
               "was 0.50, 95% CI [0.20, 0.80].", fixed = TRUE)
  expect_match(as.character(seq_report(chk, digits = 4)),
               "was 0.5004, 95% CI [0.2046, 0.7963].", fixed = TRUE)
  expect_match(as.character(seq_report(chk, digits = 0)),
               "was 1, 95% CI [0, 1].", fixed = TRUE)
  expect_identical(chk$value[2], 0.50044)
  for (bad in list(-1, 2.5, NA, "3", c(2, 3))) {
    expect_error(seq_report(chk, digits = bad), "'digits' must be")
  }
  # A double holds 15 significant decimal digits; more are refused, not
  # printed as noise with a warning from formatC().
  expect_error(seq_report(chk, digits = 16),
               "'digits' must be a single whole number from 0 to 15")
  expect_error(seq_report(chk, digits = 400), "got 400")
  expect_warning(seq_report(chk, digits = 15), NA)
  expect_match(as.character(seq_report(chk, digits = 15)),
               "was 0.500440000000000, 95% CI", fixed = TRUE)
  # A target the user gave as a fraction is shown to 7 significant digits,
  # not to 15; a whole number is shown in full.
  third <- new_smsd_check(
    c("lower_limit", "smd", "upper_limit"), c(0.1, 0.25, 0.4),
    "standardized mean difference", is_satisfied = TRUE, current_n = 1234567,
    criterion_n = 1234567, conf_level = 0.955, width = 1 / 3
  )
  txt <- as.character(seq_report(third))
  expect_match(txt, "a target full width of 0.3333333 for the 95.5% confidence",
               fixed = TRUE)
  expect_match(txt, "was met at n = 1,234,567 per group.", fixed = TRUE)
  expect_identical(attr(third, "width"), 1 / 3)
  # A value that rounds to zero is not shown with a minus sign.
  near_zero <- new_smsd_check(
    c("lower_limit", "smd", "upper_limit"), c(-0.3, -0.0002, 0.3),
    "standardized mean difference", is_satisfied = TRUE, current_n = 50,
    criterion_n = 50, conf_level = 0.9, width = 0.6
  )
  expect_match(as.character(seq_report(near_zero, digits = 2)),
               "was 0.00, 90% CI [-0.30, 0.30].", fixed = TRUE)
})

test_that("a risk procedure reports its loss constants and no interval", {
  rep <- seq_report(known_risk())
  expect_identical(
    as.character(rep),
    paste(
      "Sample size was determined with the sequential minimum risk point",
      "estimation procedure for the coefficient of variation (Chattopadhyay",
      "& Kelley, 2016), with a loss constant of A = 400,000 and a sampling",
      "cost of c = 75 per observation. The stopping rule was met at n = 35.",
      "At that sample size the estimate of the coefficient of variation was",
      "0.220."
    )
  )
  expect_match(attr(rep, "reference"),
               "627\u2013648. doi:10.1080/00273171.2016.1203279", fixed = TRUE)

  two_costs <- new_smsd_check(
    c("smd", "current_n1", "current_n2"), c(0.4, 60, 90),
    "minimum risk point estimation (standardized mean difference)",
    is_satisfied = TRUE, current_n = 150, criterion_n = 150,
    extra = list(loss_A = 50000, cost_1 = 2, cost_2 = 8)
  )
  txt <- as.character(seq_report(two_costs))
  expect_match(txt, paste("with a loss constant of A = 50,000 and sampling",
                          "costs of 2 and 8 per observation in the first and",
                          "second groups."), fixed = TRUE)
  expect_match(txt, "was met at a total sample size of N = 150.",
               fixed = TRUE)

  bounded <- new_smsd_check(
    c("mean", "sd", "risk", "optimal_n", "current_n"), c(2.5, 1, 0.4, 40, 40),
    "bounded risk estimation (normal mean)", is_satisfied = TRUE,
    current_n = 40, criterion_n = 40,
    extra = list(loss_A = 100, k = 1, risk_bound = 0.5)
  )
  # The object carries the power of the error in its loss (k), so the loss
  # is stated with it.
  expect_match(as.character(seq_report(bounded)),
               paste("with a loss of A |e|^k for an estimation error of e",
                     "(A = 100, k = 1) and a risk bound of w = 0.5."),
               fixed = TRUE)

  # Without the constants on the object the sentence simply ends.
  bare <- new_smsd_check(
    c("cv", "current_n"), c(0.2, 35),
    "minimum risk point estimation (coefficient of variation)",
    is_satisfied = TRUE, current_n = 35, criterion_n = 35
  )
  expect_match(as.character(seq_report(bare)),
               "(Chattopadhyay & Kelley, 2016). The stopping rule",
               fixed = TRUE)
})

test_that("notes and an unusable variance are carried into the paragraph", {
  # The package's standard note for such a check, word for word (it is
  # written in R/var_floor.R; if its wording is changed there, change it
  # here too): the variance estimate was zero, negative, or not finite, so
  # no variance has been estimated and the stopping rule is not met.
  standard <- paste(
    "The variance estimate is not yet usable at this sample size (it is",
    "zero, negative, or not finite), so the stopping rule is reported as not",
    "met. Collect more data and check again."
  )
  expect_identical(.note_variance_unusable(), standard)
  unusable <- known_check(is_satisfied = FALSE,
                          extra = list(variance_unusable = TRUE),
                          notes = .note_variance_unusable())
  txt <- as.character(seq_report(unusable))
  expect_match(txt, paste("Note:", standard), fixed = TRUE)
  expect_identical(lengths(regmatches(txt, gregexpr("Note:", txt))), 1L)
  expect_false(grepl("floor", txt, fixed = TRUE))
  # The flag alone is enough: an object that is flagged and carries no note
  # gets the standard one. The object holds a criterion sample size (121)
  # and limits (0.205, 0.796), as known_check() built it; neither rests on
  # an estimated variance, so the paragraph gives neither.
  flag_only <- known_check(is_satisfied = FALSE,
                           extra = list(variance_unusable = TRUE))
  expect_identical(attr(flag_only, "criterion_n"), 121)
  expect_identical(flag_only$value[flag_only$term == "upper_limit"], 0.79632)
  expect_identical(
    as.character(seq_report(flag_only)),
    paste(
      "Sample size is being determined with the sequential accuracy in",
      "parameter estimation procedure for the standardized mean difference",
      "(Kelley, Darku, & Chattopadhyay, 2018), with a target full width of",
      "0.6 for the 95% confidence interval. At n = 88 per group the stopping",
      "rule has not yet been met, so sampling continues. The interim",
      "estimate of the standardized mean difference is 0.500. Note:",
      standard
    )
  )
  # A check function that explains the unusable variance in a note of its
  # own (as several do) is not given the standard note as well.
  own_note <- known_check(
    is_satisfied = FALSE, extra = list(variance_unusable = TRUE),
    notes = "The scores in group 2 are all equal, so no variance exists yet."
  )
  txt <- as.character(seq_report(own_note))
  expect_match(txt, paste("is 0.500. Note: The scores in group 2 are all",
                          "equal, so no variance exists yet."), fixed = TRUE)
  expect_identical(lengths(regmatches(txt, gregexpr("Note:", txt))), 1L)
  expect_false(grepl("Collect more data", txt, fixed = TRUE))
  expect_false(grepl("criterion sample size|CI \\[|121|0\\.796", txt))
  # Not flagged: no note is added, and the criterion sample size and the
  # interval are given. The flag is what takes them out.
  plain <- as.character(seq_report(
    known_check(is_satisfied = FALSE, extra = list(variance_unusable = FALSE))
  ))
  expect_false(grepl("Note:", plain))
  expect_match(plain, paste("(the criterion sample size at this stage is",
                            "121), so sampling continues. The interim",
                            "estimate of the standardized mean difference is",
                            "0.500, 95% CI [0.205, 0.796]."), fixed = TRUE)

  # A check function of the package whose object holds limits at such a
  # stage: eight equal scores have no variance, and the fixed width interval
  # is the mean plus and minus d whatever the data, 5 - 0.3 and 5 + 0.3.
  flat <- seq_fwd_norm_mean(rep(5, 8), d = 0.3)
  expect_true(attr(flat, "variance_unusable"))
  expect_equal(flat$value[flat$term %in% c("lower_limit", "upper_limit")],
               c(4.7, 5.3))
  txt <- as.character(seq_report(flat))
  expect_match(txt, paste("At n = 8 the stopping rule has not yet been met,",
                          "so sampling continues. The interim estimate of the",
                          "mean is 5.000. Note:"), fixed = TRUE)
  expect_false(grepl("CI [", txt, fixed = TRUE))

  two_notes <- known_check(notes = c("First note.", "Second note."))
  expect_match(as.character(seq_report(two_notes)),
               "Note: First note. Note: Second note.", fixed = TRUE)
  expect_false(grepl("Note:", as.character(seq_report(known_check()))))

  not_finite <- new_smsd_check(
    c("lower_limit", "sm", "upper_limit"), c(NA_real_, NA_real_, NA_real_),
    "standardized mean", is_satisfied = FALSE, current_n = 10,
    criterion_n = 12, conf_level = 0.95, width = 0.4,
    notes = .note_not_estimable()
  )
  # The constructor refuses NaN; a check that holds it all the same (changed
  # after it was returned) must still be reported without printing it.
  not_finite$value <- c(NaN, NaN, NaN)
  txt <- as.character(seq_report(not_finite))
  expect_match(txt, paste("The estimate of the standardized mean is not",
                          "finite at this sample size."), fixed = TRUE)
  expect_match(txt, paste("Note:", .note_not_estimable()), fixed = TRUE)
  expect_false(grepl("NaN", txt, fixed = TRUE))
})

test_that("a path's report counts the earlier stages with no usable variance", {
  # A check built on the object contract alone: the rule is met once n
  # reaches `target`, except that below `usable_from` the variance estimate
  # is flagged as unusable and the rule reported as not met.
  flagged <- function(x, target, usable_from = 0) {
    n <- length(x)
    unusable <- n < usable_from
    new_smsd_check(
      term = c("mean", "criterion_n", "current_n"),
      value = c(mean(x), target, n), method = "toy procedure",
      is_satisfied = !unusable && n >= target, current_n = n,
      criterion_n = target, extra = list(variance_unusable = unusable),
      notes = if (unusable) .note_variance_unusable() else NULL
    )
  }
  x <- c(2, 4, 6, 8, 10, 12, 14, 16)
  sentence <- function(stages) {
    paste("At", stages, "the variance estimate was zero, negative, or not",
          "finite, so the stopping rule was reported as not met there.")
  }
  # Stages n = 2, ..., 8; unusable at n = 2 to 5; the rule (n >= 3) is first
  # met at n = 6, the fifth stage. All four earlier stages were unusable.
  all_four <- seq_path(flagged, x, target = 3, usable_from = 6, start = 2)
  expect_identical(attr(all_four, "stop_n"), 6L)
  expect_identical(
    as.character(seq_report(all_four, digits = 1)),
    paste(
      "Sample size was determined with a sequential procedure (toy",
      "procedure). The stopping rule was evaluated from n = 2 onward, after",
      "each additional observation, and was first met at n = 6. At that",
      "sample size the estimate of mean was 7.0.",
      sentence("all 4 earlier stages")
    )
  )
  # Unusable at n = 2 and 3 only; the rule (n >= 7) is first met at n = 7,
  # the sixth stage: 2 of the 5 earlier stages.
  some <- seq_path(flagged, x, target = 7, usable_from = 4, start = 2)
  expect_identical(attr(some, "stop_n"), 7L)
  expect_match(as.character(seq_report(some)),
               paste("was 8.000.", sentence("2 of the 5 earlier stages")),
               fixed = TRUE)
  # One earlier stage, and it was unusable.
  one <- seq_path(flagged, x, target = 3, usable_from = 6, start = 5)
  expect_identical(attr(one, "stop_n"), 6L)
  expect_match(as.character(seq_report(one)),
               paste("was 7.000.", sentence("the one earlier stage")),
               fixed = TRUE)
  # Two earlier stages (n = 4 and 5). Both unusable: "both", not "all 2".
  both <- seq_path(flagged, x, target = 3, usable_from = 6, start = 4)
  expect_identical(attr(both, "stop_n"), 6L)
  expect_match(as.character(seq_report(both)),
               paste("was 7.000.", sentence("both earlier stages")),
               fixed = TRUE)
  # One of the two (n = 4; at n = 5 the variance is usable and the rule,
  # n >= 6, is not yet met).
  half <- seq_path(flagged, x, target = 6, usable_from = 5, start = 4)
  expect_identical(attr(half, "stop_n"), 6L)
  expect_match(as.character(seq_report(half)),
               paste("was 7.000.", sentence("1 of the 2 earlier stages")),
               fixed = TRUE)
  # Stages after the stopping stage do not count, and a path with no such
  # stage says nothing: here every stage is usable up to the stop.
  none <- seq_path(flagged, x, target = 3, start = 2)
  expect_false(grepl("variance", as.character(seq_report(none))))
  first <- seq_path(flagged, x, target = 3, usable_from = 6, start = 6)
  expect_false(grepl("variance", as.character(seq_report(first))))
  # Never met, and still unusable at the last stage: the earlier stages are
  # counted and the check's own note about the last stage follows. The
  # criterion sample size the last check holds (3, which the 8 observations
  # in hand exceed) is not given, because the stage is flagged.
  never <- as.character(seq_report(
    seq_path(flagged, x, target = 3, usable_from = 99, start = 2)
  ))
  expect_match(
    never,
    paste("after each additional observation. At n = 8 the stopping rule",
          "has not yet been met, so sampling continues. The interim estimate",
          "of mean is 9.000.", sentence("all 6 earlier stages"), "Note:",
          .note_variance_unusable()),
    fixed = TRUE
  )
  # The same path with a usable last stage gives the criterion.
  waiting <- as.character(seq_report(
    seq_path(flagged, x, target = 30, usable_from = 4, start = 2)
  ))
  expect_match(waiting,
               paste("At n = 8 the stopping rule has not yet been met (the",
                     "criterion sample size at this stage is 30), so sampling",
                     "continues."), fixed = TRUE)
  # A check object has no earlier stages.
  expect_false(grepl("earlier", as.character(seq_report(
    flagged(x, target = 3, usable_from = 99)
  ))))
  # A check that reports no flag says nothing about its variance estimate
  # (the path holds NA), so its stages are not counted as unusable: here the
  # rule (n >= 5) is first met at the fourth stage, after three such stages.
  no_flag <- function(x, target) {
    new_smsd_check(c("mean", "criterion_n", "current_n"),
                   c(mean(x), target, length(x)), "toy procedure",
                   is_satisfied = length(x) >= target,
                   current_n = length(x), criterion_n = target)
  }
  silent <- seq_path(no_flag, x, target = 5, start = 2)
  expect_identical(attr(silent, "stop_n"), 5L)
  expect_identical(silent$variance_unusable, rep(NA, 7))
  expect_false(grepl("earlier|variance", as.character(seq_report(silent))))
})

test_that("an unknown method is reported without a reference", {
  unknown <- new_smsd_check(c("lower_limit", "omega", "upper_limit"),
                            c(0.7, 0.8, 0.9), "a method of the future",
                            is_satisfied = TRUE, current_n = 140,
                            criterion_n = 131, conf_level = 0.9, width = 0.2)
  expect_error(rep <- seq_report(unknown), NA)
  expect_identical(
    as.character(rep),
    paste(
      "Sample size was determined with a sequential procedure (a method of",
      "the future), with a target full width of 0.2 for the 90% confidence",
      "interval. The stopping rule was met at n = 140. At that sample size",
      "the estimate of omega was 0.800, 90% CI [0.700, 0.900]."
    )
  )
  expect_identical(attr(rep, "reference"), character(0))
  expect_false(any(grepl("Reference", format(rep))))
})

test_that("a path is reported at its stopping stage, with how it was checked", {
  withr::local_seed(113)
  g1 <- rnorm(80, mean = 0.5)
  g2 <- rnorm(80)
  path <- seq_path(seq_aipe_smd, list(g1, g2), width = 0.8, by = 5)
  stop_n <- attr(path, "stop_n")
  expect_false(is.na(stop_n))
  expect_lt(stop_n, 80L)
  at_stop <- seq_aipe_smd(g1[seq_len(stop_n)], g2[seq_len(stop_n)],
                          width = 0.8)
  rep <- seq_report(path, digits = 3)
  # The path's report is the report of the check at the stopping stage, with
  # the sentence about the rule saying where checking began and how often.
  expect_identical(attr(path, "final_check"), at_stop)
  pilot <- seq_aipe_smd_pilot(width = 0.8)
  expect_identical(
    as.character(rep),
    sub("The stopping rule was met at",
        paste0("The stopping rule was evaluated from the pilot sample size ",
               "(n = ", pilot, " per group) onward, after every 5 ",
               "additional observations per group, and was first met at"),
        as.character(seq_report(at_stop, digits = 3)), fixed = TRUE)
  )
  expect_match(as.character(rep), paste0("first met at n = ", stop_n,
                                         " per group."), fixed = TRUE)
  expect_identical(attr(rep, "reference"),
                   attr(seq_report(at_stop), "reference"))

  # One observation at a time, from a start the user chose.
  every <- as.character(seq_report(
    seq_path(seq_aipe_smd, list(g1, g2), width = 0.8, start = 12)
  ))
  expect_match(every, paste("evaluated from n = 12 per group onward, after",
                            "each additional observation per group, and was",
                            "first met at"), fixed = TRUE)

  # Never met: the last stage, as an interim result.
  never <- as.character(seq_report(
    seq_path(seq_aipe_smd, list(g1, g2), width = 0.1, start = 70)
  ))
  expect_match(never, "Sample size is being determined with", fixed = TRUE)
  expect_match(never, paste("additional observation per group. At n = 80 per",
                            "group the stopping rule has not yet been met"),
               fixed = TRUE)
  expect_match(never, "The interim estimate of", fixed = TRUE)

  # A single stage.
  one <- as.character(seq_report(
    seq_path(seq_aipe_smd, list(g1, g2), width = 0.8, start = 80)
  ))
  expect_match(one, paste("The stopping rule was evaluated at n = 80 per",
                          "group, and was first met at n = 80 per group."),
               fixed = TRUE)
})

test_that("two stage and budget objects are worded for what they are", {
  two_stage <- new_smsd_check(
    c("lower_limit", "mean", "upper_limit", "optimal_n", "current_n"),
    c(2.6, 2.9, 3.2, 39, 60), "two stage fixed width interval (normal mean)",
    is_satisfied = TRUE, current_n = 60, criterion_n = 39,
    conf_level = 0.95, width = 0.6
  )
  txt <- as.character(seq_report(two_stage))
  expect_match(txt, paste("the two stage fixed width confidence interval",
                          "procedure for the mean of a normal distribution",
                          "(Mukhopadhyay & de Silva, 2009)"), fixed = TRUE)
  expect_match(txt, paste("With n = 60, the total sample size the procedure",
                          "requires (39) had been reached."), fixed = TRUE)
  pilot_only <- new_smsd_check(
    c("mean", "risk", "optimal_n", "current_n"), c(2.9, 0.4, 28, 10),
    "two stage bounded risk estimation (exponential mean)",
    is_satisfied = FALSE, current_n = 10, criterion_n = 28,
    extra = list(loss_A = 1, risk_bound = 0.5)
  )
  expect_match(as.character(seq_report(pilot_only)),
               paste("With n = 10, the total sample size the procedure",
                     "requires (28) has not yet been reached, so sampling",
                     "continues."), fixed = TRUE)
  # A two stage object that flags its variance estimate: the total it holds
  # rests on no estimated variance and is not given.
  flagged_pilot <- new_smsd_check(
    c("mean", "risk", "optimal_n", "current_n"), c(2.9, 0.4, 11, 10),
    "two stage bounded risk estimation (normal mean)",
    is_satisfied = FALSE, current_n = 10, criterion_n = 11,
    extra = list(loss_A = 1, risk_bound = 0.5, variance_unusable = TRUE)
  )
  expect_match(as.character(seq_report(flagged_pilot)),
               paste("With n = 10, the total sample size the procedure",
                     "requires has not yet been reached, so sampling",
                     "continues."), fixed = TRUE)

  budget <- new_smsd_check(
    c("lower_limit", "contrast", "upper_limit", "variance"),
    c(0.2, 0.5, 0.8, 0.02), "budget constrained K-population contrast",
    is_satisfied = TRUE, current_n = 300, criterion_n = 300,
    conf_level = 0.9, extra = list(budget = 3500),
    notes = "The budget is exhausted."
  )
  rep <- seq_report(budget)
  expect_match(as.character(rep),
               paste("(Chattopadhyay, Bandyopadhyay, Kelley, & Padalunkal,",
                     "2025), with a 90% confidence interval and a total",
                     "budget of 3,500."), fixed = TRUE)
  # An object that does not say how it ended: "sampling ended" is true
  # whether every group reached its criterion or the budget ran out, and
  # "the stopping rule was met" would be false in the second case.
  expect_match(as.character(rep),
               paste("total budget of 3,500. Sampling ended at a total",
                     "sample size of N = 300. At that sample size"),
               fixed = TRUE)
  expect_false(grepl("was met", as.character(rep), fixed = TRUE))
  expect_match(as.character(rep), "Note: The budget is exhausted.",
               fixed = TRUE)
  expect_match(attr(rep, "reference"),
               paste("Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., &",
                     "Padalunkal, J. J. (2025). A sequential approach for",
                     "noninferiority or equivalence of a linear contrast",
                     "under cost constraints. Psychological Methods, 30(2),",
                     "425\u2013439. doi:10.1037/met0000570"),
               fixed = TRUE)
})

test_that("a budget that ran out is not reported as a rule that was met", {
  # The budget constrained procedure ends when every group has reached its
  # criterion or when the budget cannot pay for another needed observation.
  # In both cases is_satisfied is TRUE (sampling has ended); the attributes
  # rule_met and budget_exhausted say which.
  budget_check <- function(rule_met, budget_exhausted, is_satisfied,
                           current_n, criterion_n, notes = NULL, ...) {
    new_smsd_check(
      c("lower_limit", "contrast", "upper_limit", "variance"),
      c(0.339, 0.673, 1.007, 0.029),
      "budget constrained K-population contrast",
      is_satisfied = is_satisfied, current_n = current_n,
      criterion_n = criterion_n, conf_level = 0.95,
      extra = list(budget = 360, rule_met = rule_met,
                   budget_exhausted = budget_exhausted, ...),
      notes = notes
    )
  }
  opening <- paste(
    "Sample size was determined with the sequential budget constrained",
    "procedure for a linear contrast of population means (Chattopadhyay,",
    "Bandyopadhyay, Kelley, & Padalunkal, 2025), with a 95% confidence",
    "interval and a total budget of 360."
  )
  closing <- paste("At that sample size the estimate of the contrast was",
                   "0.673, 95% CI [0.339, 1.007].")
  exhausted <- as.character(seq_report(budget_check(
    rule_met = FALSE, budget_exhausted = TRUE, is_satisfied = TRUE,
    current_n = 180, criterion_n = 236
  )))
  expect_identical(
    exhausted,
    paste(
      opening,
      "Sampling ended at a total sample size of N = 180 because the budget",
      "could not pay for another observation from a group still below its",
      "criterion sample size (the criterion sample sizes at that point came",
      "to 236 in total).",
      closing
    )
  )
  expect_false(grepl("stopping rule was met", exhausted, fixed = TRUE))

  met <- as.character(seq_report(budget_check(
    rule_met = TRUE, budget_exhausted = FALSE, is_satisfied = TRUE,
    current_n = 240, criterion_n = 236
  )))
  expect_identical(
    met,
    paste(opening,
          "Sampling ended at a total sample size of N = 240, when every",
          "group had reached its criterion sample size.", closing)
  )

  ongoing <- as.character(seq_report(budget_check(
    rule_met = FALSE, budget_exhausted = FALSE, is_satisfied = FALSE,
    current_n = 90, criterion_n = 236
  )))
  expect_match(ongoing, "Sample size is being determined with", fixed = TRUE)
  expect_match(ongoing, paste("At a total sample size of N = 90 the stopping",
                              "rule has not yet been met (the criterion",
                              "sample size at this stage is 236), so sampling",
                              "continues. The interim estimate of the",
                              "contrast is 0.673"), fixed = TRUE)
  expect_false(grepl("Sampling ended", ongoing, fixed = TRUE))

  # Under way, with the variance of a group not estimable: the total the
  # object holds (320) includes a placeholder for that group, and its
  # interval leaves that group's variance out, so neither is given.
  own_note <- paste("The observations in group 2 are all equal, so the",
                    "variance of the contrast cannot be estimated.")
  waiting <- as.character(seq_report(budget_check(
    rule_met = FALSE, budget_exhausted = FALSE, is_satisfied = FALSE,
    current_n = 90, criterion_n = 320, variance_unusable = TRUE,
    notes = own_note
  )))
  expect_identical(
    waiting,
    paste(sub("was determined", "is being determined", opening, fixed = TRUE),
          "At a total sample size of N = 90 the stopping rule has not yet",
          "been met, so sampling continues. The interim estimate of the",
          "contrast is 0.673. Note:", own_note)
  )

  # The budget ran out while the variance of a group could not be estimated:
  # the check reports the rule as not met (is_satisfied is FALSE), but
  # sampling has ended, and the paragraph must not say that it continues.
  # The total of the criterion sample sizes is left out, because the
  # criterion of such a group is a placeholder, and so is the interval.
  closing_no_interval <- paste("At that sample size the estimate of the",
                               "contrast was 0.673.")
  stuck <- as.character(seq_report(budget_check(
    rule_met = FALSE, budget_exhausted = TRUE, is_satisfied = FALSE,
    current_n = 180, criterion_n = 320, variance_unusable = TRUE,
    notes = own_note
  )))
  expect_identical(
    stuck,
    paste(opening,
          "Sampling ended at a total sample size of N = 180 because the",
          "budget could not pay for another observation from a group still",
          "below its criterion sample size.", closing_no_interval, "Note:",
          own_note)
  )
  expect_false(grepl("continues|interim|is being|Collect more data", stuck))
  # The same object without a note: what is added does not advise collecting
  # data that the budget cannot pay for.
  bare <- as.character(seq_report(budget_check(
    rule_met = FALSE, budget_exhausted = TRUE, is_satisfied = FALSE,
    current_n = 180, criterion_n = 320, variance_unusable = TRUE
  )))
  expect_identical(
    bare,
    paste(opening,
          "Sampling ended at a total sample size of N = 180 because the",
          "budget could not pay for another observation from a group still",
          "below its criterion sample size.", closing_no_interval,
          "Note: The variance estimate is not usable at this sample size (it",
          "is zero, negative, or not finite), so the stopping rule is",
          "reported as not met and no confidence interval is given.")
  )

  # Data that already cost more than the budget are reported as that,
  # whichever of the other two attributes is set.
  for (flags in list(c(TRUE, FALSE), c(FALSE, TRUE))) {
    over <- as.character(seq_report(budget_check(
      rule_met = flags[1], budget_exhausted = flags[2], is_satisfied = TRUE,
      current_n = 180, criterion_n = 198, budget_exceeded = TRUE
    )))
    expect_identical(
      over,
      paste(opening, "Sampling ended at a total sample size of N = 180,",
            "which cost more than the budget.", closing)
    )
  }
  # The attribute FALSE changes nothing.
  expect_identical(
    as.character(seq_report(budget_check(
      rule_met = TRUE, budget_exhausted = FALSE, is_satisfied = TRUE,
      current_n = 240, criterion_n = 236, budget_exceeded = FALSE
    ))),
    met
  )
  # The wording belongs to the budget constrained procedure: another
  # procedure whose object happens to carry such attributes is worded as
  # any other.
  other <- known_check(extra = list(rule_met = FALSE,
                                    budget_exhausted = TRUE))
  expect_match(as.character(seq_report(other)),
               "The stopping rule was met at n = 88 per group.", fixed = TRUE)
})

test_that("seq_budget() objects carry the attributes the report reads", {
  # Three groups of 60 at costs 1, 2, and 3: the data cost
  # 60 * (1 + 2 + 3) = 360.
  withr::local_seed(113)
  groups <- list(stats::rnorm(60, 0.4), stats::rnorm(60),
                 stats::rnorm(60, -0.2))
  w <- c(1, -0.5, -0.5)
  run <- function(data = groups, budget, closed = NULL) {
    seq_budget(data, contrast_weights = w, cost_per_group = c(1, 2, 3),
               budget = budget, closed = closed)
  }
  flags <- function(chk) {
    out <- attributes(chk)[c("is_satisfied", "rule_met", "budget_exhausted",
                             "budget_exceeded", "variance_unusable")]
    # Present, and each one TRUE or FALSE.
    expect_identical(lengths(out),
                     c(is_satisfied = 1L, rule_met = 1L, budget_exhausted = 1L,
                       budget_exceeded = 1L, variance_unusable = 1L))
    unlist(out)
  }
  state <- function(is_satisfied, rule_met, budget_exhausted, budget_exceeded,
                    variance_unusable) {
    c(is_satisfied = is_satisfied, rule_met = rule_met,
      budget_exhausted = budget_exhausted, budget_exceeded = budget_exceeded,
      variance_unusable = variance_unusable)
  }
  ended_at_180 <- "Sampling ended at a total sample size of N = 180"

  # The budget is spent exactly, with groups still below their criteria.
  spent <- run(budget = 360)
  expect_identical(flags(spent), state(TRUE, FALSE, TRUE, FALSE, FALSE))
  expect_identical(attr(spent, "budget"), 360)
  expect_gt(attr(spent, "criterion_n"), 180)
  txt <- as.character(seq_report(spent))
  expect_match(txt, "a 95% confidence interval and a total budget of 360.",
               fixed = TRUE)
  # The data cost 60 * (1 + 2 + 3) = 360, all of the budget, so nothing was
  # left; the sentence says so with the amounts.
  expect_identical(attr(spent, "remaining_budget"), 0)
  expect_match(
    txt,
    paste0(ended_at_180, " because the budget could not pay for another ",
           "observation from a group still below its criterion sample ",
           "size: 0 of the 360 was left (the criterion sample sizes at ",
           "that point came to ", attr(spent, "criterion_n"), " in total)."),
    fixed = TRUE
  )
  expect_false(grepl("stopping rule was met", txt, fixed = TRUE))
  # seq_budget() says in a note of its own that the budget is spent. The
  # sentence above says the same, so that note is not repeated: the ending
  # is stated once, and the budget is written one way throughout.
  expect_match(attr(spent, "notes")[1L],
               "^Sampling has ended because the budget cannot ")
  expect_false(grepl("Note: Sampling has ended", txt, fixed = TRUE))
  expect_identical(lengths(regmatches(txt, gregexpr("budget could not pay",
                                                    txt))), 1L)
  expect_identical(lengths(regmatches(txt, gregexpr("Note:", txt))),
                   length(attr(spent, "notes")) - 1L)

  # A budget of 2,000 leaves room for more: sampling continues.
  going <- run(budget = 2000)
  expect_identical(flags(going), state(FALSE, FALSE, FALSE, FALSE, FALSE))
  txt <- as.character(seq_report(going))
  expect_match(txt, "Sample size is being determined with", fixed = TRUE)
  # Every variance is estimable, so the total of the criterion sample sizes
  # and the interval are given.
  expect_gt(attr(going, "criterion_n"), 180)
  expect_match(
    txt,
    paste0("a total budget of 2,000. At a total sample size of N = 180 the ",
           "stopping rule has not yet been met (the criterion sample size ",
           "at this stage is ",
           format(attr(going, "criterion_n"), big.mark = ","),
           "), so sampling continues. The interim estimate of the contrast ",
           "is "),
    fixed = TRUE
  )
  expect_match(txt, ", 95% CI [", fixed = TRUE)

  # Every group closed under the first passage rule: the rule is met.
  closed <- run(budget = 400, closed = c(TRUE, TRUE, TRUE))
  expect_identical(flags(closed), state(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_match(as.character(seq_report(closed)),
               paste0(ended_at_180, ", when every group had reached its ",
                      "criterion sample size."), fixed = TRUE)

  # A budget of 300 is less than the 360 the data cost.
  over <- run(budget = 300)
  expect_identical(flags(over), state(TRUE, FALSE, TRUE, TRUE, FALSE))
  # The cost and the budget are in the sentence (360 is 60 more than 300),
  # and the note of seq_budget() that says the same is not repeated.
  txt <- as.character(seq_report(over))
  expect_match(txt,
               paste0(ended_at_180, ", which cost 360, more than the budget ",
                      "of 300, so the budget constraint does not hold for ",
                      "this result."),
               fixed = TRUE)
  expect_match(attr(over, "notes")[1L], "^The data already cost ")
  expect_false(grepl("Note: The data already cost", txt, fixed = TRUE))
  expect_identical(lengths(regmatches(txt, gregexpr("Note:", txt))),
                   length(attr(over, "notes")) - 1L)

  # One group constant: its variance cannot be estimated. With the budget
  # spent, sampling has ended although the rule is reported as not met.
  flat <- groups
  flat[[2]] <- rep(3, 60)
  stuck <- run(flat, budget = 360)
  expect_identical(flags(stuck), state(FALSE, FALSE, TRUE, FALSE, TRUE))
  txt <- as.character(seq_report(stuck))
  expect_match(txt, "Sample size was determined with", fixed = TRUE)
  # The contrast of the three sample means, by hand; the object's interval
  # (which its own note says is too narrow to rely on) is not given.
  contrast <- sum(w * vapply(flat, mean, numeric(1)))
  expect_true(all(is.finite(
    stuck$value[stuck$term %in% c("lower_limit", "upper_limit")]
  )))
  expect_match(txt,
               paste0(ended_at_180, " because the budget could not pay for ",
                      "another observation from a group still below its ",
                      "criterion sample size: 0 of the 360 was left. At ",
                      "that sample size the estimate of the contrast was ",
                      formatC(contrast, digits = 3, format = "f"),
                      ". Note: "), fixed = TRUE)
  expect_false(grepl("CI [", txt, fixed = TRUE))
  # Every note of the object but the one the sentence above words (that the
  # budget is spent) is carried over, and nothing in the paragraph says
  # that sampling continues or that more data should be collected.
  spent_note <- startsWith(attr(stuck, "notes"),
                           "Sampling has ended because the budget cannot ")
  expect_identical(sum(spent_note), 1L)
  expect_gt(sum(!spent_note), 0L)
  for (note in attr(stuck, "notes")[!spent_note]) {
    expect_match(txt, paste("Note:", note), fixed = TRUE)
  }
  expect_identical(lengths(regmatches(txt, gregexpr("Note:", txt))),
                   sum(!spent_note))
  expect_false(grepl("sampling continues|interim|Collect more data", txt))
  # With budget to spare the same data are an interim result.
  waiting <- run(flat, budget = 2000)
  expect_identical(flags(waiting), state(FALSE, FALSE, FALSE, FALSE, TRUE))
  txt <- as.character(seq_report(waiting))
  # The total the object holds includes a placeholder for the constant
  # group, so it is not given, and neither is the interval.
  expect_true(is.finite(attr(waiting, "criterion_n")))
  expect_match(txt,
               paste0("At a total sample size of N = 180 the stopping rule ",
                      "has not yet been met, so sampling continues. The ",
                      "interim estimate of the contrast is ",
                      formatC(contrast, digits = 3, format = "f"),
                      ". Note: "), fixed = TRUE)
  expect_false(grepl("criterion sample size at this stage|CI \\[", txt))
  expect_identical(lengths(regmatches(txt, gregexpr("Note:", txt))),
                   length(attr(waiting, "notes")))
})

test_that("check objects from the package are reported with their source", {
  withr::local_seed(113)
  x <- rgamma(60, shape = 9, rate = 3)
  y <- 0.5 * x + rnorm(60)
  cases <- list(
    list(seq_aipe_cv(x, width = 0.2), "Kelley, Darku, & Chattopadhyay, 2018",
         "the coefficient of variation"),
    list(seq_aipe_r(x, y, width = 0.4, method = "kendall"),
         "Kelley, Darku, & Chattopadhyay, 2019", "Kendall's tau"),
    list(seq_aipe_r(x, y, width = 0.4, method = "spearman"),
         "Kelley, Darku, & Chattopadhyay, 2019", "Spearman's rho"),
    list(seq_aipe_r(x, y, width = 0.4),
         "Kelley, Darku, & Chattopadhyay, 2019",
         "the Pearson correlation coefficient"),
    list(seq_aipe_gini(x, width = 0.1), "Chattopadhyay & De, 2016",
         "the Gini index"),
    list(seq_mr_cv(x, A = 400000, c = 75), "Chattopadhyay & Kelley, 2016",
         "the coefficient of variation"),
    list(seq_fwd_norm_mean(x, d = 0.3), "Mukhopadhyay & de Silva, 2009",
         "the mean")
  )
  for (case in cases) {
    chk <- case[[1]]
    rep <- seq_report(chk)
    s <- summary(chk)
    txt <- as.character(rep)
    expect_match(txt, paste0("(", case[[2]], ")"), fixed = TRUE)
    expect_match(txt, paste("of", case[[3]]), fixed = TRUE)
    expect_match(txt, formatC(s$estimate, digits = 3, format = "f"),
                 fixed = TRUE)
    expect_length(attr(rep, "reference"), 1L)
  }
})

# What each method label is credited to, written out label by label. The
# credits were read from the sources: Kelley, Darku, and Chattopadhyay
# (2018) treats the standardized mean difference, the coefficient of
# variation, the standardized mean, the regression coefficient in simple
# linear regression, and contrasts (its Examples 1 to 4 and 6) and has no
# procedure for the Gini index, which is Chattopadhyay and De (2016); the
# 2019 article treats the three correlation coefficients and extends the
# method to the squared multiple correlation coefficient; the classical
# procedures are those of Mukhopadhyay and de Silva (2009).
expected_source <- c(
  "coefficient of variation" = "kdc2018",
  "standardized mean" = "kdc2018",
  "standardized mean difference" = "kdc2018",
  "simple linear regression slope" = "kdc2018",
  "contrast" = "kdc2018",
  "Gini index" = "cd2016",
  "Pearson correlation" = "kdc2019",
  "Kendall correlation" = "kdc2019",
  "Spearman correlation" = "kdc2019",
  "squared multiple correlation, Wald interval" = "kdc2019",
  "squared multiple correlation, log transformed interval" = "kdc2019",
  "minimum risk point estimation (coefficient of variation)" = "ck2016",
  "minimum risk point estimation (standardized mean difference)" = "ck2017",
  "minimum risk point estimation (normal mean)" = "mds2009",
  "minimum risk point estimation (exponential mean)" = "mds2009",
  "minimum risk point estimation (negative exponential location)" = "mds2009",
  "minimum risk point estimation (distribution free mean, arbitrary variance)" =
    "mds2009",
  "bounded risk estimation (normal mean)" = "mds2009",
  "fixed width interval (normal mean)" = "mds2009",
  "fixed width interval (distribution free mean)" = "mds2009",
  "fixed width interval (negative exponential location)" = "mds2009",
  "bounded length confidence interval (distribution free)" = "mds2009",
  "two stage fixed width interval (normal mean)" = "mds2009",
  "two stage fixed width interval (distribution free mean)" = "mds2009",
  "two stage fixed width interval (negative exponential location)" = "mds2009",
  "two stage bounded risk estimation (normal mean)" = "mds2009",
  "two stage bounded risk estimation (exponential mean)" = "mds2009",
  "budget constrained K-population contrast" = "cbkp2025"
)

# The full references, as the first pages of the articles give them (title,
# volume, issue, pages, and DOI) and, for the two sources not on hand, as
# the package's list of citations gives them. They are written as the
# reference lists of the help pages show them: an en dash in a page range,
# the DOI as doi:..., and the ISBN of the book.
expected_reference <- c(
  kdc2018 = paste(
    "Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in",
    "parameter estimation for a general class of effect sizes: A",
    "sequential approach. Psychological Methods, 23(2), 226\u2013243.",
    "doi:10.1037/met0000127"
  ),
  kdc2019 = paste(
    "Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential",
    "accuracy in parameter estimation for population correlation",
    "coefficients. Psychological Methods, 24(4), 492\u2013515.",
    "doi:10.1037/met0000203"
  ),
  ck2016 = paste(
    "Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient",
    "of variation with minimum risk: A sequential method for minimizing",
    "sampling error and study cost. Multivariate Behavioral Research,",
    "51(5), 627\u2013648. doi:10.1080/00273171.2016.1203279"
  ),
  ck2017 = paste(
    "Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized",
    "mean difference with minimum risk: Maximizing accuracy and minimizing",
    "cost with sequential estimation. Psychological Methods, 22(1),",
    "94\u2013113. doi:10.1037/met0000089"
  ),
  cbkp2025 = paste(
    "Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.",
    "(2025). A sequential approach for noninferiority or equivalence of a",
    "linear contrast under cost constraints. Psychological Methods, 30(2),",
    "425\u2013439. doi:10.1037/met0000570"
  ),
  cd2016 = paste(
    "Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within",
    "pre-specified error bound. Econometrics, 4(3), 30.",
    "doi:10.3390/econometrics4030030"
  ),
  mds2009 = paste(
    "Mukhopadhyay, N., & de Silva, B. M. (2009). Sequential methods and",
    "their applications. Chapman & Hall/CRC. ISBN:9780367386535"
  )
)

expected_cite <- c(
  kdc2018 = "Kelley, Darku, & Chattopadhyay, 2018",
  kdc2019 = "Kelley, Darku, & Chattopadhyay, 2019",
  ck2016 = "Chattopadhyay & Kelley, 2016",
  ck2017 = "Chattopadhyay & Kelley, 2017",
  cbkp2025 = "Chattopadhyay, Bandyopadhyay, Kelley, & Padalunkal, 2025",
  cd2016 = "Chattopadhyay & De, 2016",
  mds2009 = "Mukhopadhyay & de Silva, 2009"
)

# One call (or more) of every check function the package exports, with the
# labels it may produce and the target the paragraph must state: the
# constants are those passed in the call, written out here, so the test
# fails if the object stops carrying them under the names seq_report()
# reads. A check function is an exported function that has an exported
# _pilot companion. Data are passed by position and everything else by name.
report_cases <- function() {
  withr::local_seed(113)
  x <- stats::rgamma(60, shape = 9, rate = 3)
  y <- 0.5 * x + stats::rnorm(60)
  x2 <- cbind(x, stats::rnorm(60))
  g1 <- stats::rnorm(60, mean = 0.5)
  g2 <- stats::rnorm(60)
  g3 <- stats::rnorm(60, mean = -0.2)
  e <- stats::rexp(60, rate = 1 / 3)
  w <- c(1, -0.5, -0.5)
  case <- function(fn, labels, target, ...) {
    list(fn = fn, labels = labels, target = target, checks = list(...))
  }
  width <- function(value) {
    paste0(", with a target full width of ", value,
           " for the 95% confidence interval.")
  }
  loss_and_cost <- paste(", with a loss constant of A = 100 and a sampling",
                         "cost of c = 1 per observation.")
  # The three procedures that let the user choose the power of the error in
  # the loss (k) state the loss with it.
  loss_with_power <- paste(", with a loss of A |e|^k for an estimation error",
                           "of e (A = 100, k = 1)")
  loss_and_bound <- paste(loss_with_power, "and a risk bound of w = 5.")
  list(
    case("seq_aipe_cv", "coefficient of variation", width("0.2"),
         seq_aipe_cv(x, width = 0.2)),
    case("seq_aipe_sm", "standardized mean", width("0.6"),
         seq_aipe_sm(x, width = 0.6)),
    case("seq_aipe_gini", "Gini index", width("0.1"),
         seq_aipe_gini(x, width = 0.1)),
    case("seq_aipe_smd", "standardized mean difference", width("0.8"),
         seq_aipe_smd(g1, g2, width = 0.8)),
    case("seq_aipe_r",
         c("Pearson correlation", "Kendall correlation",
           "Spearman correlation"), width("0.4"),
         seq_aipe_r(x, y, width = 0.4),
         seq_aipe_r(x, y, width = 0.4, method = "kendall"),
         seq_aipe_r(x, y, width = 0.4, method = "spearman")),
    case("seq_aipe_reg_coef", "simple linear regression slope", width("0.5"),
         seq_aipe_reg_coef(x, y, width = 0.5)),
    case("seq_aipe_c", "contrast", width("1"),
         seq_aipe_c(list(g1, g2, g3), contrast_weights = w, width = 1)),
    case("seq_aipe_R2",
         c("squared multiple correlation, Wald interval",
           "squared multiple correlation, log transformed interval"),
         width("0.3"),
         seq_aipe_R2(x2, y, width = 0.3, ci_method = "wald"),
         seq_aipe_R2(x2, y, width = 0.3, ci_method = "log_transformed")),
    case("seq_mr_cv",
         "minimum risk point estimation (coefficient of variation)",
         paste(", with a loss constant of A = 400,000 and a sampling cost of",
               "c = 75 per observation."),
         seq_mr_cv(x, A = 400000, c = 75)),
    case("seq_mr_smd",
         "minimum risk point estimation (standardized mean difference)",
         paste(", with a loss constant of A = 5,000 and a sampling cost of",
               "c = 1 per observation."),
         seq_mr_smd(g1, g2, A = 5000, c = 1)),
    case("seq_mr_norm_mean", "minimum risk point estimation (normal mean)",
         loss_and_cost, seq_mr_norm_mean(x, A = 100, c = 1)),
    case("seq_mr_exp_mean",
         "minimum risk point estimation (exponential mean)",
         loss_and_cost, seq_mr_exp_mean(e, A = 100, c = 1)),
    case("seq_mr_nexp_mean",
         "minimum risk point estimation (negative exponential location)",
         paste(loss_with_power, "and a sampling cost of c = 1 per",
               "observation."),
         seq_mr_nexp_mean(e + 2, A = 100, c = 1, k = 1, t = 1)),
    case("seq_mr_dfr_mean_arb",
         paste("minimum risk point estimation (distribution free mean,",
               "arbitrary variance)"),
         loss_and_cost,
         seq_mr_dfr_mean_arb(x, A = 100, c = 1, gamma = 0.49)),
    case("seq_br_norm_mean", "bounded risk estimation (normal mean)",
         loss_and_bound, seq_br_norm_mean(x, A = 100, k = 1, w = 5)),
    # The fixed width procedures take the half width d; the paragraph gives
    # the full width, 2 d = 0.6. For the negative exponential location the
    # interval runs from the estimate down, so its whole length is d.
    case("seq_fwd_norm_mean", "fixed width interval (normal mean)",
         width("0.6"), seq_fwd_norm_mean(x, d = 0.3)),
    case("seq_fwd_dfr_mean", "fixed width interval (distribution free mean)",
         width("0.6"), seq_fwd_dfr_mean(x, d = 0.3)),
    case("seq_fwd_nexp_mean",
         "fixed width interval (negative exponential location)",
         width("0.3"), seq_fwd_nexp_mean(e + 2, d = 0.3)),
    # The bounded length interval: 2 d = 2 is the largest length allowed.
    case("seq_blci_dfr_median",
         "bounded length confidence interval (distribution free)",
         ", with a 95% confidence interval no longer than 2.",
         seq_blci_dfr_median(x, d = 1)),
    case("ts_fwd_norm_mean", "two stage fixed width interval (normal mean)",
         width("0.6"), ts_fwd_norm_mean(x[1:15], x[16:60], d = 0.3)),
    case("ts_fwd_dfr_mean",
         "two stage fixed width interval (distribution free mean)",
         width("0.6"), ts_fwd_dfr_mean(x[1:15], x[16:60], d = 0.3)),
    case("ts_fwd_nexp_mean",
         "two stage fixed width interval (negative exponential location)",
         width("0.3"), ts_fwd_nexp_mean(e[1:15] + 2, e[16:60] + 2, d = 0.3)),
    case("ts_br_norm_mean", "two stage bounded risk estimation (normal mean)",
         loss_and_bound,
         ts_br_norm_mean(x[1:15], x[16:60], A = 100, k = 1, w = 5)),
    case("ts_br_exp_mean",
         "two stage bounded risk estimation (exponential mean)",
         ", with a loss constant of A = 1 and a risk bound of w = 0.5.",
         ts_br_exp_mean(e[1:15], e[16:60], A = 1, w = 0.5)),
    case("seq_budget", "budget constrained K-population contrast",
         ", with a 95% confidence interval and a total budget of 360.",
         seq_budget(list(g1, g2, g3), contrast_weights = w,
                    cost_per_group = c(1, 2, 3), budget = 360))
  )
}

test_that("every label is credited to the source it comes from", {
  methods <- .report_methods()
  sources <- .report_sources()
  credit <- vapply(methods, function(m) m$source, character(1))
  # Label for label: nothing missing, nothing extra, nothing credited to
  # another source.
  expect_identical(credit[order(names(credit))],
                   expected_source[order(names(expected_source))])
  expect_identical(anyDuplicated(names(methods)), 0L)
  # The full reference and the citation of every source, word for word.
  expect_identical(sort(names(sources)), sort(names(expected_reference)))
  for (key in names(expected_reference)) {
    expect_identical(sources[[key]]$reference, expected_reference[[key]])
    expect_identical(sources[[key]]$cite, expected_cite[[key]])
    # Plain ASCII apart from the en dash of a page range, and no bare URL:
    # the DOI is written as on the help pages.
    expect_false(grepl("[^ -~]", gsub("\u2013", "-", sources[[key]]$reference,
                                      fixed = TRUE)))
    expect_false(grepl("[0-9]-[0-9]|https?://", sources[[key]]$reference))
  }
  for (label in names(methods)) {
    entry <- methods[[label]]
    expect_match(entry$procedure, "^the ", label = label)
    expect_true(nzchar(entry$parameter), label = label)
    # A two stage procedure is worded as reaching a required total; nothing
    # else is.
    expect_identical(entry$two_stage, grepl("^two stage", label),
                     label = label)
  }
  # How the sample size of each procedure is worded: per group where the
  # check reports the size of one group, a total where it reports the sum
  # over the groups, and plain otherwise.
  n_scale <- vapply(methods, function(m) m$n_scale, character(1))
  expect_identical(
    sort(names(n_scale)[n_scale == "per_group"]),
    c("contrast", "standardized mean difference")
  )
  expect_identical(
    sort(names(n_scale)[n_scale == "total"]),
    c("budget constrained K-population contrast",
      "minimum risk point estimation (standardized mean difference)")
  )
  expect_true(all(n_scale %in% c("", "per_group", "total")))
  stage_unit <- vapply(methods, function(m) m$stage_unit, character(1))
  expect_identical(
    sort(names(stage_unit)[stage_unit == "per_group"]),
    c("contrast",
      "minimum risk point estimation (standardized mean difference)",
      "standardized mean difference")
  )
  expect_true(all(stage_unit %in% c("", "per_group")))
  # Only the budget constrained procedure can end with its rule not met,
  # and only the bounded length interval has a width that is a bound.
  flagged <- function(field) {
    names(methods)[vapply(methods, function(m) isTRUE(m[[field]]),
                          logical(1))]
  }
  expect_identical(flagged("budget"),
                   "budget constrained K-population contrast")
  expect_identical(flagged("width_is_bound"),
                   "bounded length confidence interval (distribution free)")
})

test_that("every check function in the package is reported with its source", {
  cases <- report_cases()
  methods <- .report_methods()
  sources <- .report_sources()
  ns <- asNamespace("SMSD")
  exported <- getNamespaceExports(ns)
  check_functions <- exported[paste0(exported, "_pilot") %in% exported]
  expect_gt(length(check_functions), 20L)
  # A check function without a case here has no tested entry: add its call
  # above and its label to .report_methods(). A case for a function the
  # package does not have is a leftover: remove it and its entries.
  covered <- vapply(cases, function(case) case$fn, character(1))
  expect_identical(sort(covered), sort(check_functions))
  # Every entry of the table belongs to a case, so none is left behind when
  # a label changes.
  listed <- unlist(lapply(cases, function(case) case$labels))
  expect_identical(sort(listed), sort(names(methods)))

  for (case in cases) {
    produced <- vapply(case$checks, function(chk) attr(chk, "method"),
                       character(1))
    # The labels the function produces are the labels listed for it.
    expect_identical(sort(produced), sort(case$labels), label = case$fn)
    for (chk in case$checks) {
      label <- attr(chk, "method")
      rep <- seq_report(chk)
      source <- sources[[expected_source[[label]]]]
      expect_identical(attr(rep, "reference"), source$reference,
                       label = label)
      # The first sentence ends with the citation and the target, with the
      # constants of the call: the width and confidence level, the loss
      # constants, or the budget.
      expect_match(as.character(rep),
                   paste0(" (", source$cite, ")", case$target, " "),
                   fixed = TRUE, label = label)
      expect_identical(attr(rep, "method"), label)
    }
  }
})

test_that("the cost, the loss, and the width are worded as what they are", {
  # seq_mr_nexp_mean() prices n observations at c n^t. Only with t = 1 is c
  # the cost of one observation; with t = 2 and c = 0.5, four observations
  # cost 0.5 * 4^2 = 8, not 2.
  risk <- function(t, k = 2) {
    new_smsd_check(
      c("location", "risk", "optimal_n", "current_n"), c(2.01, 30, 8, 60),
      "minimum risk point estimation (negative exponential location)",
      is_satisfied = TRUE, current_n = 60, criterion_n = 8,
      extra = c(list(loss_A = 100, cost_c = 0.5), if (!is.null(k)) list(k = k),
                list(t = t))
    )
  }
  expect_match(as.character(seq_report(risk(2))),
               paste("with a loss of A |e|^k for an estimation error of e",
                     "(A = 100, k = 2) and a sampling cost of c n^t for n",
                     "observations (c = 0.5, t = 2). The"), fixed = TRUE)
  expect_false(grepl("per observation", as.character(seq_report(risk(2)))))
  expect_match(as.character(seq_report(risk(1))),
               paste("with a loss of A |e|^k for an estimation error of e",
                     "(A = 100, k = 2) and a sampling cost of c = 0.5 per",
                     "observation. The"), fixed = TRUE)
  # The loss is stated with the power the object carries, whatever it is,
  # and as a plain loss constant when the object carries none.
  expect_match(as.character(seq_report(risk(1, k = 1.5))),
               "(A = 100, k = 1.5) and a sampling cost", fixed = TRUE)
  expect_match(as.character(seq_report(risk(1, k = NULL))),
               paste("with a loss constant of A = 100 and a sampling cost of",
                     "c = 0.5 per observation. The"), fixed = TRUE)
  # The package's own function carries t under that name.
  withr::local_seed(113)
  e <- stats::rexp(40, rate = 1 / 3) + 2
  expect_match(
    as.character(seq_report(seq_mr_nexp_mean(e, A = 100, c = 0.5, k = 2,
                                             t = 1.5))),
    "a sampling cost of c n^t for n observations (c = 0.5, t = 1.5).",
    fixed = TRUE
  )
  # Two analyses that differ only in the power of the error in the loss are
  # different designs (the criterion sample size depends on it), and their
  # paragraphs say which is which. The three package functions that take k
  # carry it under that name.
  # The target as the first sentence states it: what follows the citation,
  # up to the period that ends the sentence.
  target <- function(chk) {
    txt <- as.character(seq_report(chk))
    rest <- substring(txt, regexpr("), with ", txt, fixed = TRUE) + 8L)
    substring(rest, 1L, regexpr(". ", rest, fixed = TRUE) - 1L)
  }
  expect_identical(
    target(seq_mr_nexp_mean(e, A = 100, c = 0.5, k = 1, t = 1)),
    paste("a loss of A |e|^k for an estimation error of e (A = 100, k = 1)",
          "and a sampling cost of c = 0.5 per observation")
  )
  expect_identical(
    target(seq_mr_nexp_mean(e, A = 100, c = 0.5, k = 3, t = 1)),
    paste("a loss of A |e|^k for an estimation error of e (A = 100, k = 3)",
          "and a sampling cost of c = 0.5 per observation")
  )
  expect_identical(
    target(seq_br_norm_mean(e, A = 100, k = 2, w = 5)),
    paste("a loss of A |e|^k for an estimation error of e (A = 100, k = 2)",
          "and a risk bound of w = 5")
  )
  expect_identical(
    target(ts_br_norm_mean(e[1:15], e[16:40], A = 100, k = 2, w = 5)),
    paste("a loss of A |e|^k for an estimation error of e (A = 100, k = 2)",
          "and a risk bound of w = 5")
  )

  # The bounded length interval for a median stops when the interval is no
  # longer than 2 d; the interval reported is usually shorter, so 2 d is not
  # a width the interval was aimed at.
  bounded <- new_smsd_check(
    c("lower_limit", "median", "upper_limit", "interval_length",
      "current_n"), c(2.354, 2.72, 2.964, 0.61, 60),
    "bounded length confidence interval (distribution free)",
    is_satisfied = TRUE, current_n = 60, criterion_n = NA_real_,
    conf_level = 0.9, width = 2
  )
  expect_identical(
    as.character(seq_report(bounded)),
    paste(
      "Sample size was determined with the distribution free sequential",
      "bounded length confidence interval procedure for a median",
      "(Mukhopadhyay & de Silva, 2009), with a 90% confidence interval no",
      "longer than 2. The stopping rule was met at n = 60. At that sample",
      "size the estimate of the median was 2.720, 90% CI [2.354, 2.964]."
    )
  )
  # The flag the report, the print method of a path, and the line under the
  # title of its plot all read: the bounded length label only.
  expect_true(.report_width_is_bound(
    "bounded length confidence interval (distribution free)"
  ))
  for (label in list("fixed width interval (normal mean)",
                     "coefficient of variation", "a method of the future",
                     NULL, NA_character_, character(0), 3)) {
    expect_false(.report_width_is_bound(label))
  }
})

test_that("the sample size is worded on the scale the check reports", {
  withr::local_seed(113)
  g1 <- stats::rnorm(40, mean = 0.5)
  g2 <- stats::rnorm(40)
  g3 <- stats::rnorm(40, mean = -0.2)
  # The contrast and the standardized mean difference report the size of one
  # group: 40 observations in each of three groups is n = 40 per group.
  contrast <- seq_aipe_c(list(g1, g2, g3),
                         contrast_weights = c(1, -0.5, -0.5), width = 5)
  expect_true(attr(contrast, "is_satisfied"))
  expect_match(as.character(seq_report(contrast)),
               "The stopping rule was met at n = 40 per group.", fixed = TRUE)
  smd <- seq_aipe_smd(g1, g2, width = 5)
  expect_true(attr(smd, "is_satisfied"))
  expect_match(as.character(seq_report(smd)),
               "The stopping rule was met at n = 40 per group.", fixed = TRUE)
  # seq_mr_smd() reports the two groups together: 40 and 40 is N = 80.
  risk <- seq_mr_smd(g1, g2, A = 20, c = 1)
  expect_true(attr(risk, "is_satisfied"))
  expect_match(as.character(seq_report(risk)),
               "was met at a total sample size of N = 80.", fixed = TRUE)
  # One cost for both groups is worded as one cost, whether the object
  # carries it once or once for each group.
  expect_match(as.character(seq_report(risk)),
               paste("with a loss constant of A = 20 and a sampling cost of",
                     "c = 1 per observation."), fixed = TRUE)
  equal_costs <- new_smsd_check(
    c("smd", "current_n1", "current_n2"), c(0.4, 75, 75),
    "minimum risk point estimation (standardized mean difference)",
    is_satisfied = TRUE, current_n = 150, criterion_n = 150,
    extra = list(loss_A = 50000, cost_1 = 3, cost_2 = 3)
  )
  expect_match(as.character(seq_report(equal_costs)),
               paste("with a loss constant of A = 50,000 and a sampling cost",
                     "of c = 3 per observation."), fixed = TRUE)
  # A single sample and paired data: plain n.
  x <- stats::rgamma(40, shape = 9, rate = 3)
  expect_match(as.character(seq_report(seq_aipe_cv(x, width = 2))),
               "The stopping rule was met at n = 40.", fixed = TRUE)
})

test_that("print wraps the paragraph and shows the reference", {
  rep <- seq_report(known_check())
  lines <- capture.output(res <- print(rep, width = 60))
  expect_identical(res, rep)
  expect_true(all(nchar(lines) <= 60))
  expect_identical(paste(lines[seq_len(which(lines == "")[1] - 1L)],
                         collapse = " "), as.character(rep))
  expect_true("Reference" %in% lines)
  expect_identical(
    paste(trimws(lines[seq.int(which(lines == "Reference") + 1L,
                               length(lines))]), collapse = " "),
    attr(rep, "reference")
  )
  expect_gt(length(capture.output(print(rep, width = 40))), length(lines))
  capture.output(expect_invisible(print(rep)))
  expect_error(print(rep, width = 5), "'width' must be")
  expect_error(print(rep, width = NA), "'width' must be")
  # as.character() gives the bare string.
  expect_identical(class(as.character(rep)), "character")
  expect_null(attributes(as.character(rep)))
})

test_that("format gives one plain string, so a report fits in a data frame", {
  rep <- seq_report(known_check())
  # format() methods return one element for each element of their input.
  expect_identical(format(rep), as.character(rep))
  expect_length(format(rep), 1L)
  # A report as a column of a data frame: built and printed without the
  # warning of a corrupt data frame that a longer format() result causes.
  expect_warning(d <- data.frame(id = 1L, report = rep), NA)
  expect_identical(dim(d), c(1L, 2L))
  expect_warning(shown <- capture.output(print(d)), NA)
  expect_true(any(grepl("Sample size was determined with", shown,
                        fixed = TRUE)))
  d2 <- data.frame(id = 1L)
  d2$report <- rep
  expect_warning(capture.output(print(d2)), NA)
  expect_identical(as.character(format(d2)$report), as.character(rep))
})

test_that("a path checked every 2 observations says so", {
  # The wording has three cases: each observation (by = 1), every `by`
  # observations, and a single stage. The first and last are tested above
  # with by = 1 and by = 5; by = 2 is the value nearest to 1.
  withr::local_seed(113)
  g1 <- rnorm(80, mean = 0.5)
  g2 <- rnorm(80)
  txt <- as.character(seq_report(
    seq_path(seq_aipe_smd, list(g1, g2), width = 0.8, start = 12, by = 2)
  ))
  expect_match(txt, paste("evaluated from n = 12 per group onward, after",
                          "every 2 additional observations per group, and",
                          "was first met at"), fixed = TRUE)
  expect_false(grepl("each additional observation", txt, fixed = TRUE))
})

test_that("numbers are written for an English text whatever the session's decimal mark", {
  # An estimate in the thousands takes a comma between the thousands and a
  # period before the decimals: 12345.678 is "12,345.678".
  large <- new_smsd_check(
    c("lower_limit", "mean", "upper_limit", "current_n"),
    c(12340.1234, 12345.678, 12351.2326, 40),
    "fixed width interval (normal mean)", is_satisfied = TRUE,
    current_n = 40, criterion_n = 40, conf_level = 0.95, width = 11.1092
  )
  expected <- paste("At that sample size the estimate of the mean was",
                    "12,345.678, 95% CI [12,340.123, 12,351.233].")
  expect_match(as.character(seq_report(large)), expected, fixed = TRUE)
  expect_match(as.character(seq_report(large)),
               "with a target full width of 11.1092 for the 95% confidence",
               fixed = TRUE)
  # With a decimal comma in the session the paragraph is the same string,
  # and no warning about the two marks being the same character.
  usual <- as.character(seq_report(large))
  known <- as.character(seq_report(known_check()))
  withr::local_options(OutDec = ",")
  expect_no_warning(comma <- as.character(seq_report(large)))
  expect_identical(comma, usual)
  expect_no_warning(known_comma <- as.character(seq_report(known_check())))
  expect_identical(known_comma, known)
  expect_match(known_comma, "was 0.500, 95% CI [0.205, 0.796].", fixed = TRUE)
})

test_that("the target is worded from whatever the object carries", {
  # A width with no confidence level: the width is given, and the limits
  # are called an interval, not a confidence interval of a stated level.
  no_level <- new_smsd_check(
    c("lower_limit", "smd", "upper_limit", "current_n"),
    c(0.20456, 0.50044, 0.79632, 88), "standardized mean difference",
    is_satisfied = TRUE, current_n = 88, criterion_n = 88, width = 0.6
  )
  txt <- as.character(seq_report(no_level))
  expect_match(txt, paste("2018), with a target full width of 0.6 for the",
                          "confidence interval. The stopping rule"),
               fixed = TRUE)
  expect_match(txt, "was 0.500, interval [0.205, 0.796].", fixed = TRUE)
  expect_false(grepl("% CI", txt, fixed = TRUE))

  # One design constant, and three: the list is joined as a sentence.
  risk <- function(...) {
    new_smsd_check(
      c("mean", "current_n"), c(2.5, 40),
      "minimum risk point estimation (normal mean)", is_satisfied = TRUE,
      current_n = 40, criterion_n = 40, extra = list(...)
    )
  }
  expect_match(as.character(seq_report(risk(loss_A = 400000))),
               "2009), with a loss constant of A = 400,000. The stopping",
               fixed = TRUE)
  expect_match(
    as.character(seq_report(risk(loss_A = 100, cost_c = 2,
                                 risk_bound = 0.5))),
    paste("2009), with a loss constant of A = 100, a sampling cost of c = 2",
          "per observation, and a risk bound of w = 0.5. The stopping"),
    fixed = TRUE
  )

  # An object whose notes attribute is gone is reported without notes.
  bare <- known_check()
  attr(bare, "notes") <- NULL
  expect_identical(as.character(seq_report(bare)),
                   as.character(seq_report(known_check())))
})

test_that("a budget ending is stated once, with the amounts the object carries", {
  spent_note <- paste(
    "Sampling has ended because the budget cannot pay for another",
    "observation: 12.5 of 360 is left, and one more observation from a",
    "group still below its criterion (group 1) costs at least 40."
  )
  other_note <- "Another note that the report has no sentence for."
  budget_check <- function(...) {
    new_smsd_check(
      c("lower_limit", "contrast", "upper_limit", "variance"),
      c(0.339, 0.673, 1.007, 0.029),
      "budget constrained K-population contrast",
      is_satisfied = TRUE, current_n = 180, criterion_n = 236,
      conf_level = 0.95,
      extra = list(budget = 360, rule_met = FALSE, budget_exhausted = TRUE,
                   ...),
      notes = c(spent_note, other_note)
    )
  }
  # With the amount left on the object, the sentence gives it and the note
  # that says the same is left out; any other note is kept.
  txt <- as.character(seq_report(budget_check(remaining_budget = 12.5)))
  expect_match(txt, paste("below its criterion sample size: 12.5 of the 360",
                          "was left (the criterion sample sizes at that",
                          "point came to 236 in total)."), fixed = TRUE)
  expect_false(grepl("Note: Sampling has ended", txt, fixed = TRUE))
  expect_match(txt, paste("Note:", other_note), fixed = TRUE)
  expect_identical(lengths(regmatches(txt, gregexpr("Note:", txt))), 1L)
  # An amount that is 12 up to the rounding of sums of costs reads 12.
  txt <- as.character(seq_report(budget_check(
    remaining_budget = 360 - (348 + 1e-12)
  )))
  expect_match(txt, "criterion sample size: 12 of the 360 was left",
               fixed = TRUE)
  # Without the amount the sentence cannot give it, so the note stays.
  txt <- as.character(seq_report(budget_check()))
  expect_match(txt, paste("below its criterion sample size (the criterion",
                          "sample sizes"), fixed = TRUE)
  expect_match(txt, paste("Note:", spent_note), fixed = TRUE)
  expect_identical(lengths(regmatches(txt, gregexpr("Note:", txt))), 2L)

  # The helper leaves the notes of every other ending alone.
  amounts <- list(budget = "360", remaining = "12.5", cost = "347.5")
  notes <- c(spent_note, other_note)
  expect_identical(.report_budget_notes(notes, "budget", amounts), other_note)
  expect_identical(.report_budget_notes(notes, "criterion", amounts), notes)
  expect_identical(.report_budget_notes(notes, "ended", amounts), notes)
  expect_identical(.report_budget_notes(notes, "exceeded", amounts), notes)
  expect_identical(.report_budget_notes(notes, "budget", NULL), notes)
})

test_that("seq_report() refuses anything that is not a check or a path", {
  expect_error(seq_report(),
               "'x' must be supplied: a sequential check object")
  expect_error(seq_report(1:3), "'x' must be a sequential check object")
  expect_error(seq_report(as.data.frame(known_check())),
               "'x' must be a sequential check object")
  expect_error(seq_report("text"), "class 'smsd_check'")
})
