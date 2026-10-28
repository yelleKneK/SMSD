# The display options the tests were written under.
#
# Many tests pin text the package writes with format(): a note, an error
# message, a printed check. That text follows the session's `digits` and
# `scipen`, and so does deparse(), which the hygiene detectors use to read
# every function. testthat sets OutDec, the width, the language, and the
# collation for each test, but not these two, so in a session where either
# was changed (a profile with options(scipen = 999) is common) correct code
# fails tests that compare text. The defaults are set here for the whole
# suite and restored when it ends. A test that needs another value sets it
# for itself with withr::local_options(), which takes precedence.
withr::local_options(list(digits = 7, scipen = 0),
                     .local_envir = testthat::teardown_env())
