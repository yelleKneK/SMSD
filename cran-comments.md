# cran-comments

## Submission

This is the first release of SMSD ("Sequential Methods for Study Design"),
version 1.0.0.

SMSD provides sequential procedures for planning a study when the population
values that a sample size formula needs are unknown. Data are collected in
stages, and a stopping rule computed from the data in hand says whether
sampling can stop. The package implements the methods of the articles that
`citation("SMSD")` lists and of the book cited in DESCRIPTION. It is the
sequential companion of DMAR, by the same maintainer.

## An Earlier Submission Under This Name

An earlier submission of a package with this name was not accepted, because
its sequential functions were interactive. The package has been redesigned
since. Each point raised then is answered as follows.

* **The functions wrote to the console with `cat()`.** No function writes to
  the console now, apart from the print methods registered in NAMESPACE
  (`print.smsd_check()`, `print.smsd_path()`, `print.smsd_report()`, and
  `print.smsd_sensitivity()`). What used to be printed is returned in the
  object. A test, `tests/testthat/test-rd_hygiene.R`, inspects every function
  in the namespace for console output outside a registered print method. It
  runs wherever the tests run, CRAN included.

* **The results were lists with title cased, punctuated names** (for example
  `"Is.Satisfied?"`). A check function returns a data frame of class
  `smsd_check` with the columns `term` and `value`. Whether the stopping rule
  is met is the attribute `is_satisfied`, which is always `TRUE` or `FALSE`.
  A pilot function returns the pilot sample size as an integer.

* **A function switched between modes and expected to be called again as
  the data arrived.** Each procedure is now a pair of functions with one job
  each: a `*_pilot()` function returns the pilot sample size, and a check
  function evaluates the stopping rule once on the data it is given. A check
  function does not prompt, reads nothing from the console, and keeps
  nothing between calls. The user's own code collects the data and decides
  when to call again; the examples show that loop in ordinary R.

The same test file enforces other patterns that `R CMD check` does not
report: no Rd markup inside `\eqn{}` or `\deqn{}`, no `\dontrun{}` or
`\donttest{}` (every example runs), no example guarded by
`requireNamespace()`, no commented out code in examples, no function that
sets a seed without restoring the caller's random number stream or that
touches the global environment, no option, graphical parameter, or working
directory changed without being restored, no default that writes a file,
and no function that takes more than two cores.

## Dependency on DMAR

SMSD lists DMAR (>= 1.0.0) in `Depends`. DMAR provides the fixed sample size
counterparts of the sequential procedures, and SMSD uses its `smd()`
estimator. DMAR 1.0.0 has been on CRAN since 2026-09-21.

## Test Environments

* Local, 2026-10-26: R 4.6.1, aarch64-apple-darwin23, macOS 27.0.1;
  `R CMD check --as-cran`, including remote incoming checks and both
  manuals. Strong dependencies were checked against CRAN; DMAR 1.0.0
  was used in a private library. HTML validation used HTML Tidy 5.8.0;
  DESCRIPTION spelling was explicitly enabled with Aspell 0.60.8.1 and
  its English dictionaries.

Windows R-release, Windows R-devel, Linux, and the declared minimum R
version have not yet been tested for this release candidate. Obtain and
add those results for the final tarball before submission.

## R CMD Check Results

The completed local check reports 0 errors, 0 warnings, and 1 note:

* CRAN incoming feasibility: new submission, with the author surnames
  and acronyms listed below reported as possibly misspelled.

Examples, the test suite, all twelve vignettes, and both manual builds
completed successfully. Seven longer tests are deliberately skipped
under CRAN; the release gate also runs the full suite separately.

CRAN's incoming checks add the NOTE for a new submission. It may list the
words below as possibly misspelled in DESCRIPTION. Each is spelled as
intended.

* AIPE: the acronym for accuracy in parameter estimation, defined in the
  Description.
* MRPE: the acronym for minimum risk point estimation, defined in the
  Description.
* Bandyopadhyay: the surname of an author of a cited article.
* Chattopadhyay: the surname of a package author, who is also an author of
  cited articles.
* Darku: the surname of a package author, who is also an author of cited
  articles.
* De: the surname of an author of a cited article (Chattopadhyay and De,
  2016).
* Kelley: the surname of the maintainer.
* Mukhopadhyay: the surname of an author of the cited book.
* Padalunkal: the surname of an author of a cited article.
* de Silva: the surname of an author of the cited book.
  The prefix `de` is also spelled as intended.

## Downstream Dependencies

None. No package on CRAN depends on SMSD.
