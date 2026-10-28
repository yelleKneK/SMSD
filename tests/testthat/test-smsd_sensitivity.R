# The smsd_sensitivity class: a term / value data frame that prints in fixed
# notation.

test_that("the constructor returns a classed term and value data frame", {
  out <- new_smsd_sensitivity(c("replications", "coverage"), c(1000, 0.95))
  expect_s3_class(out, "smsd_sensitivity")
  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("term", "value"))
  expect_identical(out$value, c(1000, 0.95))
  expect_identical(out$value[out$term == "coverage"], 0.95)
  expect_error(new_smsd_sensitivity("a", c(1, 2)))
  expect_error(new_smsd_sensitivity("a", "b"))
})

test_that("values print in fixed notation, each formatted on its own", {
  out <- new_smsd_sensitivity(
    c("true_value", "replications", "mean_n", "coverage", "sd_n"),
    c(0.5, 1e+05, 1234.5678, 0.9412, NA)
  )
  txt <- format(out)
  expect_length(txt, 6L)
  expect_false(any(grepl("e[+-]", txt)))
  expect_match(txt[2], "true_value +0\\.5$")
  expect_match(txt[3], "replications +100000$")
  expect_match(txt[6], "sd_n +NA$")
})

test_that("digits affects display only, and print returns its input", {
  out <- new_smsd_sensitivity(c("mean_n", "coverage"), c(1234.5678, 0.94125))
  expect_match(format(out, digits = 3)[2], "mean_n +1235$")
  expect_identical(out$value[1], 1234.5678)
  expect_identical(capture.output(res <- print(out)), format(out))
  expect_identical(res, out)
})
