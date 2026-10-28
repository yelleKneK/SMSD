# Regression tests that do not belong to one function's own test file:
# missing values are refused by every seq_aipe_* check function, max_n is
# an exact cap of the sensitivity engine when one observation is added at a
# time, and slr_slope() checks its input instead of recycling it.

test_that("seq_aipe_* reject NA input instead of silently flooring xi^2", {
  withr::local_seed(113)
  x_na <- c(rnorm(29), NA)
  ok <- rnorm(30)
  expect_error(seq_aipe_smd(x_na, ok, width = 0.7), "NA")
  expect_error(seq_aipe_r(x_na, ok, width = 0.5), "NA")
  expect_error(seq_aipe_cv(abs(x_na) + 3, width = 0.3), "NA")
  expect_error(seq_aipe_sm(x_na, width = 0.4), "NA")
  expect_error(seq_aipe_gini(abs(x_na) + 1, width = 0.2), "NA")
  expect_error(seq_aipe_reg_coef(x_na, ok, width = 0.5), "NA")
  expect_error(
    seq_aipe_c(list(x_na, ok), c(1, -1), width = 0.5),
    "NA"
  )
})

test_that("sensitivity max_n is an exact cap under one-at-a-time growth", {
  out <- seq_aipe_smd_sensitivity(true_value = 0.5, width = 0.2,
                                  G = 3, seed = 113, max_n = 25)
  # With a narrow width the rule cannot be satisfied by n = 25, so every
  # replication must stop exactly at the cap.
  expect_equal(out$value[out$term == "mean_n"], 25)
})

test_that("slr_slope validates its input instead of recycling", {
  expect_error(slr_slope(1:10, 1:9), "same length")
  expect_error(slr_slope(letters[1:4], 1:4), "numeric")
})
