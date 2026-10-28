# Values pinned against sources outside the package's own code: the
# stopping criterion of every seq_aipe_* check function against the formula
# of Kelley, Darku, and Chattopadhyay (2018, Equation 22; 2019, Equation
# 10), and the three correlation variance estimators against their known
# asymptotic values under bivariate normality.

test_that("every seq_aipe_* criterion equals ceiling((2z/w)^2 (xi^2 + 1/n))", {
  withr::local_seed(113)
  z <- stats::qnorm(0.975)
  pin <- function(out, w) {
    xi2 <- out$value[out$term == "xi_squared"]
    n <- attr(out, "current_n")
    expect_equal(out$value[out$term == "criterion_n"],
                 ceiling((2 * z / w)^2 * (xi2 + 1 / n)))
    expect_identical(attr(out, "is_satisfied"),
                     n >= out$value[out$term == "criterion_n"])
  }
  x1 <- rnorm(40, 0.5); x2 <- rnorm(40); xp <- rgamma(40, 4, scale = 2)
  pin(seq_aipe_smd(x1, x2, width = 0.7), 0.7)
  pin(seq_aipe_r(x1, x2, width = 0.5), 0.5)
  pin(seq_aipe_cv(xp, width = 0.3), 0.3)
  pin(seq_aipe_sm(x1, width = 0.5), 0.5)
  pin(seq_aipe_reg_coef(x1, x2, width = 0.6), 0.6)
  pin(seq_aipe_gini(xp, width = 0.2), 0.2)
  pin(seq_aipe_c(list(x1, x2, rnorm(40)), c(1, -0.5, -0.5), width = 0.6), 0.6)
})

test_that(".var_r Pearson approaches (1 - rho^2)^2 under bivariate normality", {
  # The standard deviation of the estimate is 0.0070 at n = 30,000 (400
  # samples), so the relative tolerance of 0.08 (0.033 in absolute terms)
  # is 4.7 of them. (At the n of 6,000 this test used before, the same
  # tolerance was two standard deviations.)
  withr::local_seed(113)
  n <- 30000; rho <- 0.6
  x <- rnorm(n); y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
  expect_equal(.var_r(x, y, method = "pearson"), (1 - rho^2)^2,
               tolerance = 0.08)
})

test_that(".var_r Kendall approaches 4/9 at independence", {
  withr::local_seed(113)
  n <- 2000
  expect_equal(.var_r(rnorm(n), rnorm(n), method = "kendall"), 4 / 9,
               tolerance = 0.15)
})

test_that(".var_r Spearman approaches 1 at independence", {
  withr::local_seed(113)
  n <- 1500
  expect_equal(.var_r(rnorm(n), rnorm(n), method = "spearman"), 1,
               tolerance = 0.15)
})
