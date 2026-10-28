# Internal machinery for the coefficient of variation minimum risk point
# estimation procedure of Chattopadhyay and Kelley (2016). Not exported.

# The correction term tuning constant shared by the two minimum risk point
# estimation articles: strictly between 0 and 1/2 (Chattopadhyay & Kelley,
# 2016, p. 636; 2017, p. 100).
#' @noRd
.mrpe_gamma <- function(gamma) {
  .check_number(gamma, "gamma", lower = 0, upper = 0.5,
                lower_open = TRUE, upper_open = TRUE)
}

# Round a criterion up to a whole number exactly as .as_sample_size() does
# (its ceiling ignores representation error), but without its error: a
# criterion that is not finite or too large to represent is returned as it
# is. The rule functions below decide "met" with this value and the check
# functions report .as_sample_size() of the same criterion, so the decision
# and the reported criterion cannot disagree; and because nothing here can
# fail, a simulated replication cannot abort.
#' @noRd
.mrpe_ceiling <- function(criterion) {
  if (is.finite(criterion) && criterion <= .Machine$integer.max) {
    as.numeric(.as_sample_size(criterion))
  } else {
    criterion
  }
}

# The power of two by which the rule functions of this family divide the
# data before any sum of squares is formed, so that nothing they compute
# depends on the units of the data. `r` holds the range (or the ranges) of
# the data. When the largest absolute value lies between 2^-200 and 2^200
# (about 6e-61 and 2e60) no sum of squares can overflow or underflow, and
# the scale is 1: the caller leaves the data as they are and no copy is
# made. Outside that band it is .aipe_scale(), the power of two of the
# order of the largest absolute value, which seq_aipe_cv() and
# seq_aipe_smd() use in the same way. Dividing by a power of two is exact,
# so data on any scale give the same bits as the same data on an ordinary
# one.
#' @noRd
.mrpe_scale <- function(r) {
  largest <- max(abs(r))
  if (largest >= 2^-200 && largest <= 2^200) 1 else .aipe_scale(r)
}

# Chattopadhyay and Kelley (2016), Equation 22, as printed:
#   V_n^2 = s^4 / (4 xbar^4) + mu4_hat / (4 xbar^4)
#           + s^2 / (2 xbar^2) - mu3_hat / xbar^3,
# with the unbiased U-statistic estimators of the third and fourth central
# moments (their Equations 20 and 21). It equals the asymptotic variance of
# the sample coefficient of variation only when the population has no excess
# kurtosis (see inst/cv_variance_note.md).
#
# x: numeric vector (at least 4 values) with a positive mean.
# floor: TRUE returns max{V_n^2, n^-3}, the bound the other variance
#   estimators in the package carry; FALSE returns the raw estimate, which
#   is what .mr_cv_rule() needs in order to tell an unusable estimate (zero,
#   negative, or not finite) from a small positive one.
#
# The expression does not depend on the units of the data, so it is
# evaluated on the data divided by their mean. Evaluated on the data as
# given, the fourth power of the mean overflows for data of magnitude 1e77
# and underflows at 1e-81, where .var_cv() (which rescales in the same
# spirit) still returns the estimate.
#' @noRd
.var_cv_published <- function(x, floor = TRUE) {
  n <- length(x)
  z <- x / mean(x)
  mz <- mean(z)
  vz <- stats::var(z)
  value <- vz^2 / (4 * mz^4) + .mu4(z) / (4 * mz^4) + vz / (2 * mz^2) -
    .mu3(z) / mz^3
  .var_floor(value, n, floor)
}

# One evaluation of the stopping rule on a cleaned data vector. Shared by
# seq_mr_cv() and seq_mr_cv_sensitivity(), so the Monte Carlo wrapper runs
# exactly the rule the check function reports.
#
# Returns a list: estimable (is the sample mean positive beyond rounding
# error? If so there is a finite estimate), cv, xi2 (the variance estimate
# in use: the raw estimate bounded below by n^-3), unusable (was the raw
# estimate zero, negative, or not finite?), criterion (the right side of the
# stopping rule on the scale of n, not rounded; it can be very large and is
# never converted with an error here, so this function cannot fail), risk
# (the estimated risk at the current n), and satisfied.
#
# The unusable case follows the package wide rule in R/var_floor.R: no
# variance has been estimated, so the rule is not met. A positive raw
# estimate, however small, is used as the articles prescribe.
#' @noRd
.mr_cv_rule <- function(x, A, c, gamma, variance) {
  n <- length(x)
  # Nothing computed here depends on the units of the data, which are
  # divided by a power of two when they are on an extreme scale (see
  # .mrpe_scale()). The sum of squares behind the standard deviation can
  # then neither overflow (computed on the data as given, it does for data
  # of magnitude 1e154 and beyond, which would leave no estimate) nor
  # underflow (at 1e-160 the estimate would be wrong in its fifth digit, and
  # at 1e-200 it would be 0, with the rule reported as met).
  r <- c(min(x), max(x))
  scale <- .mrpe_scale(r)
  if (scale != 1) {
    x <- x / scale
    r <- r / scale
  }
  mx <- mean(x)
  # The mean must be positive by more than rounding error. Centered or
  # standardized scores have a mean of about 1e-17 with either sign; treating
  # a positive one as a real mean would give a coefficient of variation near
  # 1e16. The tolerance has the form of the one stats::t.test() uses to call
  # data essentially constant (max(abs(r)) is the largest absolute value).
  if (!(mx > 10 * .Machine$double.eps * max(abs(r)))) {
    return(list(estimable = FALSE, cv = NA_real_, xi2 = NA_real_,
                unusable = FALSE, criterion = NA_real_, risk = NA_real_,
                satisfied = FALSE))
  }
  # Finite: the sum of squares is in range and the mean is more than 2e-15
  # times the largest absolute value.
  cv <- stats::sd(x) / mx
  # .is_constant() needs only the range.
  raw <- if (.is_constant(r)) {
    # No variation beyond rounding error: nothing to estimate a variance
    # from (the moment ratios would be rounding noise).
    NaN
  } else if (variance == "published") {
    .var_cv_published(x, floor = FALSE)
  } else {
    .var_cv(x, floor = FALSE)
  }
  unusable <- .var_unusable(raw)
  xi2 <- .var_floor(raw, n)
  # Chattopadhyay and Kelley (2016), Equation 23: stop at the first n with
  # n^2 >= (A / c) (V_n^2 + n^(-2 gamma)); `criterion` is the same condition
  # on the scale of n.
  criterion <- sqrt((A / c) * (xi2 + n^(-2 * gamma)))
  list(estimable = TRUE, cv = cv, xi2 = xi2,
       unusable = unusable, criterion = criterion,
       risk = A * xi2 / n + c * n,
       satisfied = !unusable && isTRUE(n >= .mrpe_ceiling(criterion)))
}

# The note attached when the sample mean is zero or negative, or positive
# only by rounding error.
#' @noRd
.note_cv_mean_not_positive <- function() {
  paste0("The sample mean is not positive (it is zero, negative, or not ",
         "distinguishable from zero at machine precision), so the ",
         "coefficient of variation cannot be estimated from these data and ",
         "the stopping rule is reported as not met. The procedure is for a ",
         "population with a positive mean.")
}

# The note attached when the published (2016) variance is requested.
#' @noRd
.note_cv_published <- function() {
  paste0("The variance estimate is Equation 22 of Chattopadhyay and Kelley ",
         "(2016), as requested with variance = \"published\". It agrees ",
         "with the default estimator for normal data and can differ ",
         "otherwise; see ?seq_mr_cv.")
}


#' Pilot Sample Size for Minimum Risk Point Estimation of the Coefficient of
#' Variation
#'
#' Returns the pilot (first stage) sample size for the purely sequential
#' minimum risk point estimation (MRPE) procedure for the population
#' coefficient of variation \eqn{\kappa = \sigma / \mu} of Chattopadhyay and
#' Kelley (2016). It needs only the two cost constants, so it can be computed
#' before any data are collected.
#'
#' @details
#' The pilot size is Equation 24 of Chattopadhyay and Kelley (2016),
#' \deqn{m = \max\left\{m_0,\ \left\lceil (A/c)^{1/(2 + 2\gamma)}
#'   \right\rceil\right\},}
#' where \eqn{\lceil\cdot\rceil} rounds up to the next whole number. The
#' stopping rule of \code{\link{seq_mr_cv}} can never be met below
#' \eqn{(A/c)^{1/(2 + 2\gamma)}} observations, whatever the data, which is why
#' the pilot starts there. In the article's worked example (\eqn{A} = 400,000,
#' \eqn{c} = 75, \eqn{\gamma} = 0.49; p. 637) the pilot size is 18.
#'
#' @param A The loss constant: the price of one unit of squared estimation
#'   error, in the same currency as \code{c}. The article builds it as the
#'   amount one is willing to pay to estimate the coefficient of variation to
#'   within \eqn{\epsilon}, divided by \eqn{\epsilon^2} (for example
#'   \eqn{1000 / 0.05^2 = 400000}). A single positive number.
#' @param c The cost of sampling one observation. A single positive number.
#' @param m_0 The smallest pilot size allowed: a whole number of at least 4
#'   (the default), because the variance estimator uses the fourth central
#'   moment, which needs four observations.
#' @param gamma The tuning constant of the stopping rule's correction term,
#'   strictly between 0 and 1/2. The default, 0.49, is the value
#'   Chattopadhyay and Kelley (2016) suggest. Use the same value in
#'   \code{\link{seq_mr_cv}}.
#'
#' @return A length one integer, the pilot sample size.
#'
#' @references
#' Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient
#' of variation with minimum risk: A sequential method for minimizing
#' sampling error and study cost. \emph{Multivariate Behavioral Research},
#' \emph{51}(5), 627--648. \doi{10.1080/00273171.2016.1203279}
#'
#' @author
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in},
#' Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # The worked example of Chattopadhyay and Kelley (2016, p. 637): a home
#' # visit costs 75 dollars, and the team would pay 1,000 dollars to
#' # estimate the coefficient of variation to within 0.05.
#' seq_mr_cv_pilot(A = 1000 / 0.05^2, c = 75)
#'
#' @seealso \code{\link{seq_mr_cv}} for the stopping rule and
#'   \code{\link{seq_mr_cv_sensitivity}} to study the procedure before
#'   collecting data; the accuracy in parameter estimation companion
#'   \code{\link{seq_aipe_cv}}; in \pkg{DMAR}, \code{\link[DMAR]{mr_cv}}.
#' @family seq_mr_cv
#' @export
seq_mr_cv_pilot <- function(A, c, m_0 = 4, gamma = 0.49) {
  A <- .check_positive(A, "A")
  c <- .check_positive(c, "c")
  m_0 <- .check_count(m_0, "m_0", min = 4L)
  gamma <- .mrpe_gamma(gamma)
  max(m_0, .as_sample_size((A / c)^(1 / (2 + 2 * gamma)),
                           "pilot sample size",
                           advice = "Check 'A' and 'c'."))
}


#' Minimum Risk Point Estimation of the Coefficient of Variation
#'
#' Given the data collected so far, evaluates the purely sequential minimum
#' risk point estimation (MRPE) stopping rule of Chattopadhyay and Kelley
#' (2016) for the population coefficient of variation
#' \eqn{\kappa = \sigma / \mu}. The procedure is for a study in which each
#' observation costs money and a less accurate estimate also has a cost: it
#' stops sampling when the estimated total of the two is as small as it can
#' be. It targets a risk, not a confidence interval, so there is no
#' confidence level.
#'
#' The function does not collect data or keep state. Call
#' \code{\link{seq_mr_cv_pilot}} for the starting sample size, collect that
#' many observations, then call \code{seq_mr_cv()} after each new observation
#' (or batch) until \code{attr(out, "is_satisfied")} is \code{TRUE}.
#'
#' @details
#' \strong{The risk and the stopping rule.} The loss of estimating
#' \eqn{\kappa} by the sample coefficient of variation \eqn{k_n = s_n /
#' \bar{X}_n} from \eqn{n} observations is the squared error priced at
#' \eqn{A} plus the sampling cost, and its expectation is the risk
#' \eqn{R_n(\kappa) = A\,E[(k_n - \kappa)^2] + c\,n} (Equation 11 of
#' Chattopadhyay and Kelley, 2016). For samples that are not too small this
#' is approximately \eqn{A\,\xi^2 / n + c\,n} (their Equation 14), where
#' \eqn{\xi^2} is the asymptotic variance of \eqn{\sqrt{n}\,(k_n - \kappa)}.
#' The approximate risk is smallest at \eqn{n_c = (A/c)^{1/2}\,\xi} (their
#' Equation 15), where it equals \eqn{2 c\, n_c} (their Equation 16). Because
#' \eqn{\xi} is unknown, the purely sequential rule (their Equation 23) stops
#' at the first \eqn{n \ge m} with
#' \deqn{n^2 \ge \frac{A}{c}\left(\hat{\xi}_n^2 + n^{-2\gamma}\right),}
#' where \eqn{m} is the pilot size, \eqn{\hat{\xi}_n^2} estimates
#' \eqn{\xi^2} from the data, and \eqn{n^{-2\gamma}} with
#' \eqn{0 < \gamma < 1/2} keeps the rule from stopping too early while the
#' estimate is still noisy. The reported \code{optimal_n} is the same
#' criterion on the scale of \eqn{n},
#' \eqn{\lceil (A/c)^{1/2} (\hat{\xi}_n^2 + n^{-2\gamma})^{1/2} \rceil}, and
#' the reported \code{risk} is \eqn{A\,\hat{\xi}_n^2 / n + c\,n}. The
#' coefficient of variation and the estimate \eqn{\hat{\xi}_n^2} are ratios
#' in which the units of the data cancel, so the result is the same
#' whatever the units (centimeters or meters, dollars or thousands of
#' dollars), up to rounding error.
#'
#' \strong{Two variance estimators.} The argument \code{variance} chooses
#' \eqn{\hat{\xi}_n^2}.
#' \itemize{
#'   \item \code{"corrected"} (the default; use this one) is the estimator of
#'     Kelley, Darku, and Chattopadhyay (2018, Equations 31 and 32), the
#'     sample version of the asymptotic variance in their Equation 30,
#'     \deqn{\xi^2 = \frac{\mu_4}{4\mu^2\sigma^2} - \frac{\sigma^2}{4\mu^2}
#'       - \frac{\mu_3}{\mu^3} + \frac{\sigma^4}{\mu^4},}
#'     where \eqn{\mu_3} and \eqn{\mu_4} are the third and fourth central
#'     moments. It is consistent for any population with a positive mean and
#'     a finite fourth moment, and it is the estimator
#'     \code{\link{seq_aipe_cv}} uses.
#'   \item \code{"published"} is Equation 22 of Chattopadhyay and Kelley
#'     (2016), the sample version of their Equation 9,
#'     \deqn{\frac{\mu_4}{4\mu^4} + \frac{\sigma^4}{4\mu^4}
#'       - \frac{\mu_3}{\mu^3} + \frac{\sigma^2}{2\mu^2}.}
#'     It is offered so that the article's worked example and tables, and
#'     the results of \code{\link[DMAR]{mr_cv}} in \pkg{DMAR} and
#'     \code{mr.cv()} in \pkg{MBESS}, can be reproduced.
#' }
#' The two expressions differ by
#' \eqn{\frac{1}{4}\kappa^2 (1 - \kappa^2)(\mu_4/\sigma^4 - 3)}. They are
#' therefore equal for a normal population, and for any population with no
#' excess kurtosis, and they differ otherwise. For gamma data with
#' \eqn{\kappa = 0.8}, for example, \eqn{n} times the mean squared error of
#' \eqn{k_n} at \eqn{n = 400} was 0.516 in a simulation of 200,000 samples
#' (standard error 0.002); the corrected expression gives 0.525 and the
#' published one 0.304. In the full procedure for that population (\eqn{A} =
#' 5,000,000, \eqn{c} = 10; 5,000 replications) the risk is smallest at 512
#' observations; the mean final sample size was 479 with the corrected
#' estimator and about 374 with the published one (the article's Table 1
#' reports 374.0). The derivation, with code that
#' a reader can run, is in the file \code{cv_variance_note.md} installed
#' with the package
#' (\code{system.file("cv_variance_note.md", package = "SMSD")}). Both
#' estimators use the unbiased estimators of \eqn{\mu_3} and \eqn{\mu_4}
#' (Equations 20 and 21 of the 2016 article), and a positive estimate from
#' either is bounded below by \eqn{n^{-3}}, as in Equation 31 of the 2018
#' article (Equation 22 of the 2016 article has no such bound; the package
#' applies it to the published estimator as well).
#'
#' \strong{When the variance cannot be estimated.} If the variance estimate
#' is zero, negative, or not finite (every observation so far is the same,
#' or the same to ten significant digits, or the sample gives a negative
#' estimate), the function reports the stopping rule as not met, sets
#' \code{attr(out, "variance_unusable")} to \code{TRUE}, and adds a note;
#' see \code{\link{smsd_check_methods}}. This is a deliberate departure from
#' the article's rule, which read literally would stop at once on such an
#' estimate. A positive estimate, however small, is used as the articles
#' prescribe, so a small coefficient of variation does not by itself make
#' an estimate unusable.
#'
#' What the departure costs depends on how often the estimate is negative.
#' For normal data that is seldom. With a coefficient of variation of 0.10,
#' \eqn{A} = 400, and \eqn{c} = 1 (a pilot size of 8), the mean final sample
#' size in 2,000 simulated runs was 8.05, where the rule as printed, with
#' the bound in place of a negative estimate, gave 8.00. It is frequent for
#' a quantity that takes only two or three distinct values and has a small
#' coefficient of variation: the asymptotic variance of such a population
#' is close to zero (it is \eqn{\kappa^4} for two equally likely values),
#' the default estimate is a difference of terms that nearly cancel, and it
#' is negative at many stages in a row. For observations equal to 9 or 11
#' with equal probability (\eqn{\kappa} = 0.10) and the same constants, the
#' rule as printed stopped at 8 in every one of 2,000 runs. This function
#' reports the rule as met at 8 with probability exactly one half (when 1,
#' 2, 5, 6, or 7 of the 8 observations are 11) and keeps sampling
#' otherwise: in the same 2,000 runs it stopped at 8 in 48.5 percent, the
#' mean final sample size was 25, the 90th percentile 51, the 99th
#' percentile 244, and the largest 693. (With
#' \code{variance = "published"} the mean was 8.01, because
#' Equation 22 seldom turns negative for such data; for this population,
#' however, its large sample value is about fifty times the asymptotic
#' variance.) If the quantity you measure takes very few distinct values,
#' decide on a largest sample size before sampling begins.
#'
#' \strong{Agreement with \code{mr_cv()} in \pkg{DMAR}.} With
#' \code{variance = "published"} this function and
#' \code{\link[DMAR]{mr_cv}} apply the same rule to the same number whenever
#' Equation 22 is at least \eqn{n^{-3}}, as it is in the article's example.
#' \code{mr_cv()} applies Equation 22 as printed, with no lower bound, so
#' for a smaller value, which occurs in small samples, the two can differ,
#' in two ways.
#' \itemize{
#'   \item Equation 22 is zero or negative. This function reports the rule
#'     as not met; \code{mr_cv()} uses the value and can stop. In 4,000
#'     samples of 8 observations from a lognormal population (log scale
#'     standard deviation 1, \eqn{A} = 400, \eqn{c} = 1), Equation 22 was
#'     zero or negative in 223 (5.6 percent), and \code{mr_cv()} stopped on
#'     each of them.
#'   \item Equation 22 is positive but below \eqn{n^{-3}}. This function
#'     uses \eqn{n^{-3}} in its place, so its \code{risk} row is the larger
#'     of the two, and so is its criterion: the decisions differ when
#'     \eqn{n^2} lies between \eqn{(A/c)(V_n^2 + n^{-2\gamma})} and
#'     \eqn{(A/c)(n^{-3} + n^{-2\gamma})}, with \eqn{V_n^2} for Equation 22.
#'     The band is narrow (with the default \eqn{\gamma}, at most 1.5
#'     percent of \eqn{n^2} at \eqn{n} = 8 and 0.6 percent at \eqn{n} = 13),
#'     but when \eqn{A/c} puts \eqn{n^2} inside it the two functions
#'     disagree for every sample with so small an estimate. For the eight
#'     heights 168.2, 171.9, 166.4, 173.5, 170.8, 169.1, 174.6, and 167.3
#'     (Equation 22 is 0.00015, and \eqn{8^{-3}} is 0.00195) with \eqn{c} =
#'     1, both functions report the rule as met at \eqn{A} = 400 and as not
#'     met at \eqn{A} = 500. At \eqn{A} = 485, \code{mr_cv()} reports it as
#'     met (63.3 against \eqn{n^2} = 64) and this function as not met (64.1
#'     against 64). In 4,000 samples of 8 normal observations with a
#'     coefficient of variation of 0.02, the two decisions were the same in
#'     every sample at \eqn{A} = 400 and different in every sample at
#'     \eqn{A} = 485.
#' }
#'
#' \strong{When the mean is not positive.} If the sample mean is zero or
#' negative, the coefficient of variation is not estimable: the estimate,
#' risk, variance, and criterion are \code{NA}, the rule is reported as not
#' met, and a note says why. A mean that differs from zero only by rounding
#' error (no larger than about 2e-15 times the largest absolute observation,
#' as for centered or standardized scores) is treated as zero, whatever its
#' sign. A positive mean that is very close to zero gives a very large
#' coefficient of variation, and the function stops with an error if the
#' criterion that follows is too large to represent.
#'
#' \strong{What to expect.} The article proves that the expected final
#' sample size and the risk approach their optimal values as \eqn{A/c} grows;
#' at practical sizes the agreement is approximate. In simulations of the
#' default procedure with \eqn{A} = 5,000,000 and \eqn{c} = 10 (5,000
#' replications each), the ratio of the mean final sample size to
#' \eqn{n_c} was 1.18 for normal and for gamma data with \eqn{\kappa = 0.2}
#' (there \eqn{n_c} is near the pilot size, so the correction term matters
#' most), 1.01 for normal data with \eqn{\kappa = 0.4}, and about 0.99,
#' 0.95, and 0.94 for gamma data with \eqn{\kappa} = 0.4, 0.6, and 0.8; the
#' ratio of the realized risk to the minimum risk was within about 0.03 of 1
#' in each case (its Monte Carlo standard error is about 0.01). When
#' \eqn{A/c} is small and the population is strongly skewed the
#' procedure stops early: the variance estimate rests on the third and
#' fourth sample moments, which are very variable in small samples, and
#' sampling ends the first time the criterion is met. For gamma data with
#' \eqn{\kappa = 1.3} and \eqn{c = 1} the ratio of sample sizes was about
#' 0.58, 0.74, 0.85, and 0.93 at \eqn{A} = 900, 9,000, 90,000, and 900,000
#' (0.86 and 0.94 for the last two in a second set of runs). Use
#' \code{\link{seq_mr_cv_sensitivity}} to see the behavior for your own
#' constants before collecting data.
#'
#' @param x A numeric vector: the observations collected so far, at least
#'   four. The coefficient of variation is meaningful for a quantity
#'   measured on a ratio scale with a positive mean (the article assumes
#'   nonnegative data); it is not meaningful for centered or standardized
#'   scores or for differences, whose mean can be zero.
#' @inheritParams seq_mr_cv_pilot
#' @param gamma The tuning constant of the correction term, strictly between
#'   0 and 1/2; default 0.49. Use the same value as in
#'   \code{\link{seq_mr_cv_pilot}}.
#' @param variance Which estimator of the asymptotic variance \eqn{\xi^2} to
#'   use: \code{"corrected"} (the default), the consistent estimator of
#'   Kelley, Darku, and Chattopadhyay (2018); or \code{"published"},
#'   Equation 22 of Chattopadhyay and Kelley (2016), for reproducing that
#'   article. See Details.
#' @param na.rm Should missing values be removed before the computation?
#'   \code{TRUE} (the default) removes them; \code{FALSE} stops with an
#'   error if any are present.
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} with
#'   columns \code{term} and \code{value} and one row for each of
#'   \describe{
#'     \item{\code{cv}}{the sample coefficient of variation, \eqn{s_n /
#'       \bar{X}_n};}
#'     \item{\code{risk}}{the estimated risk at the current sample size,
#'       \eqn{A\,\hat{\xi}_n^2 / n + c\,n};}
#'     \item{\code{xi_squared}}{the variance estimate \eqn{\hat{\xi}_n^2}
#'       from the estimator chosen by \code{variance}, bounded below by
#'       \eqn{n^{-3}} (it equals \eqn{n^{-3}} when the variance cannot be
#'       estimated);}
#'     \item{\code{optimal_n}}{the criterion sample size, rounded up: the
#'       stopping rule is met exactly when \code{current_n} is at least
#'       \code{optimal_n}. When the variance cannot be estimated it is at
#'       least \code{current_n + 1};}
#'     \item{\code{current_n}}{the number of observations supplied.}
#'   }
#'   Attributes: \code{is_satisfied} (whether to stop sampling: \code{TRUE}
#'   or \code{FALSE}, never \code{NA}), \code{current_n},
#'   \code{criterion_n} (equal to \code{optimal_n}), \code{method},
#'   \code{loss_A}, \code{cost_c}, \code{gamma}, \code{variance} (the
#'   estimator used), \code{variance_unusable}, and \code{notes} (a
#'   character vector, printed under the result). The attribute
#'   \code{variance_unusable} is \code{TRUE} when the variance estimate is
#'   zero, negative, or not finite (see Details).
#'   \code{conf_level} and \code{width} are \code{NULL}. When the sample
#'   mean is not positive, \code{cv}, \code{risk}, \code{xi_squared},
#'   \code{optimal_n}, and \code{criterion_n} are \code{NA}. See
#'   \code{\link{smsd_check_methods}} for the print, summary, and conversion
#'   methods.
#'
#' @references
#' Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient
#' of variation with minimum risk: A sequential method for minimizing
#' sampling error and study cost. \emph{Multivariate Behavioral Research},
#' \emph{51}(5), 627--648. \doi{10.1080/00273171.2016.1203279}
#'
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#' parameter estimation for a general class of effect sizes: A sequential
#' approach. \emph{Psychological Methods}, \emph{23}(2), 226--243.
#' \doi{10.1037/met0000127}
#'
#' @author
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in},
#' Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # The worked example of Chattopadhyay and Kelley (2016, pp. 637-638): the
#' # number of age appropriate books in the homes of third graders. A home
#' # visit costs 75 dollars and A = 1000 / 0.05^2.
#' books <- c(36, 53, 19, 11, 10, 24, 14, 65, 18, 48, 25, 35, 13, 18, 3, 41,
#'            5, 3, 44, 26, 13, 39, 2, 3, 26, 22, 8, 15, 12, 22, 5, 21, 23,
#'            40, 18)
#' m <- seq_mr_cv_pilot(A = 400000, c = 75)
#' m
#'
#' # After the pilot sample the rule is not met, so sampling continues.
#' seq_mr_cv(books[1:m], A = 400000, c = 75)
#'
#' # With all 35 observations the default (corrected) variance asks for one
#' # more; the article's own variance reproduces its printed result (risk
#' # 4891.284, stop at 35).
#' seq_mr_cv(books, A = 400000, c = 75)
#' seq_mr_cv(books, A = 400000, c = 75, variance = "published")
#'
#' @seealso \code{\link{seq_mr_cv_pilot}} for the pilot size and
#'   \code{\link{seq_mr_cv_sensitivity}} to study the procedure before
#'   collecting data; the accuracy in parameter estimation companion
#'   \code{\link{seq_aipe_cv}}; in \pkg{DMAR}, \code{\link[DMAR]{mr_cv}}
#'   (the same procedure with the published variance), the estimator
#'   \code{\link[DMAR]{cv}}, and the fixed sample size planner
#'   \code{\link[DMAR]{ss_aipe_cv}}.
#' @family seq_mr_cv
#' @export
seq_mr_cv <- function(x, A, c, gamma = 0.49,
                      variance = c("corrected", "published"), na.rm = TRUE) {
  A <- .check_positive(A, "A")
  c <- .check_positive(c, "c")
  gamma <- .mrpe_gamma(gamma)
  variance <- .check_choice(variance, c("corrected", "published"),
                            "variance")
  x <- .mr_check_data(x, na.rm, min_n = 4L, arg = "x")
  n <- length(x)

  rule <- .mr_cv_rule(x, A, c, gamma, variance)
  notes <- character(0)
  optimal_n <- NA_integer_
  if (!rule$estimable) {
    # The only way to have no estimate: a sample mean that is not positive.
    notes <- .note_cv_mean_not_positive()
  } else {
    optimal_n <- .as_sample_size(
      rule$criterion, "optimal sample size",
      advice = paste0("Check 'A', 'c', and the data: a sample mean close ",
                      "to zero makes the coefficient of variation, and ",
                      "with it the criterion, very large."))
    if (rule$unusable) {
      # No variance has been estimated, so the rule is not met whatever the
      # sample size. The criterion reported is at least one more
      # observation, so that "current_n >= optimal_n" says the same as
      # is_satisfied on every path.
      optimal_n <- max(optimal_n, n + 1L)
      notes <- .note_variance_unusable()
    }
    if (variance == "published") {
      notes <- c(notes, .note_cv_published())
    }
  }

  new_smsd_check(
    term = c("cv", "risk", "xi_squared", "optimal_n", "current_n"),
    value = c(rule$cv, rule$risk, rule$xi2, optimal_n, n),
    method = "minimum risk point estimation (coefficient of variation)",
    is_satisfied = rule$satisfied,
    current_n = n,
    criterion_n = as.numeric(optimal_n),
    extra = list(loss_A = A, cost_c = c, gamma = gamma, variance = variance,
                 variance_unusable = rule$unusable),
    notes = notes
  )
}
