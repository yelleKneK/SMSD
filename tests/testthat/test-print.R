# Tests for the print.smsd_check S3 method's conditional output. The width and
# confidence-level lines appear for interval (AIPE) methods; both are omitted
# for risk methods, whose conf_level attribute is NULL, per the documented
# contract. The exact lines are pinned in test-smsd_check_methods.R.

test_that("print.smsd_check shows width and confidence for an AIPE object", {
  set.seed(113)
  out <- seq_aipe_smd(rnorm(40, 0.5), rnorm(40, 0.0), width = 0.8)
  txt <- capture.output(print(out))
  expect_true(any(grepl("Target width", txt)))
  expect_true(any(grepl("Confidence level", txt)))
})

test_that("print.smsd_check omits the confidence line for a risk method", {
  set.seed(113)
  out <- seq_mr_smd(rnorm(40, 0.5), rnorm(40, 0.0), A = 100, c = 4)
  txt <- capture.output(print(out))
  expect_false(any(grepl("Confidence level", txt)))
})

test_that("print.smsd_check returns its input invisibly", {
  set.seed(113)
  out <- seq_mr_cv(rgamma(40, shape = 4, scale = 2), A = 400, c = 1)
  capture.output(expect_invisible(print(out)))
})
