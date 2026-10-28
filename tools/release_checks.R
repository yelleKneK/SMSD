# Static release checks for SMSD.
#
# Every invariant a release of this package must satisfy that can be read
# from the source without running the package: the patterns CRAN's pretest
# and reviewers have returned packages for, version and citation lockstep,
# the house conventions that are easy to state and easy to break, and what
# may and may not be in the tarball and the repository. R CMD check reports
# none of these.
#
# Usage, from the package root:
#   Rscript tools/release_checks.R                      # the source tree
#   Rscript tools/release_checks.R --artifact=TARBALL   # a built tarball and
#                                                       # the history behind it
#   Rscript tools/release_checks.R --history            # the history alone:
#                                                       # run before any push,
#                                                       # the first one included
#
# One line is printed for each check, [PASS] or [FAIL], with the offending
# files and lines under a failure. The exit status is nonzero on any failure,
# and tools/release_gate.R runs this script in both forms and reads those
# lines. Every detector is first shown a known bad input, which it must flag,
# and a clean one, which it must pass; a detector that cannot do both fails
# its check without looking at the package, so no check passes because its
# detector is broken.
#
# The patterns CRAN forbids in package code are defined once, in
# tests/testthat/test-rd_hygiene.R, which runs them on the installed
# namespace wherever the tests run. This script evaluates that file's
# detectors and applies them to the source files, so the two cannot disagree.
#
# The script never contains the private strings it searches for: they are
# assembled from character codes below. It writes only under tempdir().
# tools/ is excluded from the tarball by .Rbuildignore.

args <- commandArgs(trailingOnly = TRUE)
is_artifact <- grepl("^--artifact=", args)
is_history <- args == "--history"
# An argument this script does not know stops it: a mistyped --artifact=
# would otherwise run the source tree checks and report them as the answer.
if (any(!is_artifact & !is_history))
  stop("Unknown argument: ", paste(args[!is_artifact & !is_history], collapse = " "),
       " (the arguments are --artifact=TARBALL and --history)", call. = FALSE)
if (sum(is_artifact) > 1L) stop("--artifact was given more than once.", call. = FALSE)
if (any(is_artifact) && any(is_history))
  stop("--artifact sweeps the history as well; give one of the two.", call. = FALSE)
# The history and the identities behind it, without a tarball. Nothing may be
# pushed that this sweep has not passed, and the first push comes before any
# tarball can be certified, so the sweep has a form that needs none.
history_only <- any(is_history)
artifact <- sub("^--artifact=", "", args[is_artifact])
if (length(artifact)) artifact <- normalizePath(artifact, mustWork = TRUE)

if (!file.exists("DESCRIPTION")) stop("Run from the package root.", call. = FALSE)
root <- normalizePath(".")

state <- new.env(parent = emptyenv())
state$failed <- 0L
state$run <- 0L

report <- function(ok, name, detail = character()) {
  state$run <- state$run + 1L
  if (!ok) state$failed <- state$failed + 1L
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", name))
  if (!ok) for (d in detail) cat("       ", d, "\n", sep = "")
  invisible(ok)
}

# Prove the detector on a known bad and a known clean input, then apply it.
# Where the bad input plants several defects, `n_bad` is how many findings it
# must give: a detector that sees one of them and misses the rest is broken.
check <- function(name, detector, target, bad, good, n_bad = NULL) {
  found_bad <- tryCatch(length(detector(bad)), error = function(e) 0L)
  flags_bad <- if (is.null(n_bad)) found_bad > 0L else found_bad == n_bad
  passes_good <- tryCatch(length(detector(good)) == 0L, error = function(e) FALSE)
  if (!flags_bad || !passes_good) {
    return(report(FALSE, name, c(
      if (!flags_bad) sprintf("the detector gave %d finding(s) on its known bad input; %s expected",
                              found_bad, if (is.null(n_bad)) "at least one was" else paste(n_bad, "were")),
      if (!passes_good) "the detector flagged its known clean input",
      "the detector is broken, so the package was not examined")))
  }
  hits <- detector(target)
  report(length(hits) == 0L, name, hits)
}

chr <- function(...) intToUtf8(c(...))

desc <- read.dcf("DESCRIPTION")
field <- function(f) if (f %in% colnames(desc)) unname(desc[1L, f]) else NA_character_
pkg <- field("Package")

# ---- Private Strings ------------------------------------------------------
# Strings that belong to the author's working environment and not to the
# package. None may reach the repository, a tarball, or the site. They are
# assembled from character codes so that this file never matches itself, and
# a finding is reported by class and location, never by content.
private_patterns <- c(
  tool      = chr(99, 108, 97, 117, 100, 101),
  tool      = chr(97, 110, 116, 104, 114, 111, 112, 105, 99),
  tool      = chr(99, 104, 97, 116, 103, 112, 116),
  tool      = chr(111, 112, 101, 110, 97, 105),
  tool      = chr(99, 111, 112, 105, 108, 111, 116),
  trailer   = chr(99, 111, 45, 97, 117, 116, 104, 111, 114, 101, 100, 45, 98, 121),
  tool_file = chr(97, 103, 101, 110, 116, 115, 91, 46, 93, 109, 100),
  machine_path = chr(47, 85, 115, 101, 114, 115, 47),
  machine_path = chr(91, 92, 92, 93, 85, 115, 101, 114, 115, 91, 92, 92, 93),
  machine_path = chr(47, 104, 111, 109, 101, 47, 91, 94, 47, 93, 43, 47),
  machine_path = chr(47, 118, 97, 114, 47, 102, 111, 108, 100, 101, 114, 115, 47),
  machine_path = chr(47, 112, 114, 105, 118, 97, 116, 101, 47, 40, 116, 109, 112,
                     124, 118, 97, 114, 41, 47),
  folder = chr(68, 114, 111, 112, 98, 111, 120),
  folder = chr(91, 47, 92, 92, 93, 81, 117, 97, 108, 105, 116, 121, 91,
               32, 95, 93, 67, 111, 110, 116, 114, 111, 108, 124, 81,
               117, 97, 108, 105, 116, 121, 91, 32, 95, 93, 67, 111,
               110, 116, 114, 111, 108, 91, 47, 92, 92, 95, 93),
  working_document = chr(65, 85, 84, 72, 79, 82, 83, 95, 84, 79, 68, 79),
  working_document = chr(115, 117, 103, 103, 101, 115, 116, 105, 111, 110, 115,
                         95, 102, 111, 114, 95, 115, 109, 115, 100),
  working_document = chr(119, 111, 114, 107, 105, 110, 103, 32, 104, 101, 108,
                         112, 32, 102, 105, 108, 101, 115),
  working_document = chr(81, 67, 95, 114, 101, 112, 111, 114, 116),
  working_document = chr(114, 101, 118, 105, 101, 119, 95, 40, 97, 105, 112,
                         101, 124, 99, 108, 97, 115, 115, 105, 99, 97, 108, 124,
                         100, 111, 99, 115, 124, 109, 114, 112, 101, 41)
)
private_bad_line <- paste("see", chr(47, 85, 115, 101, 114, 115, 47), "x and",
                          chr(99, 108, 97, 117, 100, 101))
# One sample for each pattern, in the order of `private_patterns`, entered a
# second time from character codes and matched by its own pattern only. The
# patterns cannot be proofread, since they are character codes: a mistyped one
# is found here, because its sample is then not flagged and the control that
# counts the samples fails.
private_samples <- c(
  chr(67, 108, 97, 117, 100, 101),
  chr(65, 110, 116, 104, 114, 111, 112, 105, 99),
  chr(67, 104, 97, 116, 71, 80, 84),
  chr(79, 112, 101, 110, 65, 73),
  chr(67, 111, 112, 105, 108, 111, 116),
  chr(67, 111, 45, 65, 117, 116, 104, 111, 114, 101, 100, 45, 66, 121, 58, 32, 120),
  chr(65, 71, 69, 78, 84, 83, 46, 109, 100),
  chr(47, 85, 115, 101, 114, 115, 47, 120),
  chr(67, 58, 92, 85, 115, 101, 114, 115, 92, 120),
  chr(47, 104, 111, 109, 101, 47, 120, 47, 121),
  chr(47, 118, 97, 114, 47, 102, 111, 108, 100, 101, 114, 115, 47, 120),
  chr(47, 112, 114, 105, 118, 97, 116, 101, 47, 116, 109, 112, 47, 120),
  chr(68, 114, 111, 112, 98, 111, 120),
  chr(82, 47, 81, 117, 97, 108, 105, 116, 121, 32, 67, 111, 110, 116, 114, 111, 108),
  chr(65, 85, 84, 72, 79, 82, 83, 95, 84, 79, 68, 79),
  chr(115, 117, 103, 103, 101, 115, 116, 105, 111, 110, 115, 95, 102, 111, 114,
      95, 115, 109, 115, 100),
  chr(119, 111, 114, 107, 105, 110, 103, 32, 104, 101, 108, 112, 32, 102, 105,
      108, 101, 115),
  chr(81, 67, 95, 114, 101, 112, 111, 114, 116),
  chr(114, 101, 118, 105, 101, 119, 95, 100, 111, 99, 115)
)
stopifnot(length(private_samples) == length(private_patterns))

# Classes of private string on each line of a character vector, as
# "label:line (class, class)".
private_in_lines <- function(lines, label) {
  if (!length(lines)) return(character())
  classes <- vapply(seq_along(private_patterns), function(k)
    grepl(private_patterns[[k]], lines, ignore.case = TRUE, perl = TRUE,
          useBytes = TRUE), logical(length(lines)))
  if (is.null(dim(classes))) classes <- matrix(classes, nrow = length(lines))
  hit <- which(rowSums(classes) > 0L)
  vapply(hit, function(i)
    sprintf("%s:%d (%s)", label, i,
            paste(unique(names(private_patterns)[classes[i, ]]), collapse = ", ")),
    character(1))
}

# The bytes of a file as lines; a compressed R data file is read expanded.
file_lines <- function(path) {
  size <- file.size(path)
  if (is.na(size) || size == 0) return(character())
  raw <- readBin(path, "raw", size)
  if (grepl("[.](rds|rda|RData|gz)$", path, ignore.case = TRUE)) {
    for (type in c("gzip", "bzip2", "xz")) {
      expanded <- tryCatch(memDecompress(raw, type = type), error = function(e) NULL,
                           warning = function(w) NULL)
      if (!is.null(expanded)) { raw <- expanded; break }
    }
  }
  txt <- rawToChar(raw[raw != as.raw(0L)])
  strsplit(txt, "\n", fixed = TRUE, useBytes = TRUE)[[1L]]
}

private_in_files <- function(files, root = ".") {
  out <- character()
  for (f in files) {
    path <- file.path(root, f)
    if (dir.exists(path) || !file.exists(path)) next
    out <- c(out, private_in_lines(f, paste0(f, " (file name)")),
             private_in_lines(file_lines(path), f))
  }
  sub(" \\(file name\\):1 ", " (file name) ", out)
}

control_dir <- file.path(tempdir(), "release_checks_controls")
unlink(control_dir, recursive = TRUE)
dir.create(control_dir, recursive = TRUE)
writeLines(c("A clean line.", private_samples), file.path(control_dir, "bad.txt"))
writeLines(c("A clean line.", "Quality control charts are sequential too."),
           file.path(control_dir, "good.txt"))

# ---- What May Be in a Tarball ----------------------------------------------
standard_surfaces <- c("DESCRIPTION", "NAMESPACE", "NEWS.md", "README.md", "R",
                       "build", "inst", "man", "tests", "vignettes", "MD5")
nonstandard_surfaces <- function(contents) {
  rel <- sub("^[^/]+/", "", contents[grepl("/", contents)])
  top <- unique(sub("/.*", "", rel[nzchar(rel)]))
  working <- grep(paste0("(^|/)(tools|dev|data-raw|docs|pkgdown|revdep)/|[.]Rcheck|",
                         "_problems|[.]orig$|(^|/)[.][^/]+$|[.]tar[.]gz$|",
                         "^inst/(?!(eq26_correction|cv_variance_note)[.]md$)[^/]+[.]md$"),
                  rel, value = TRUE, perl = TRUE)
  unique(c(setdiff(top, standard_surfaces), working))
}

inspect_tarball <- function(tarball, what) {
  contents <- utils::untar(tarball, list = TRUE)
  check(sprintf("%s holds only the standard package surfaces (%d files)", what,
                sum(!grepl("/$", contents))),
        nonstandard_surfaces, contents,
        bad = c("pkg/DESCRIPTION", "pkg/R/a.R", "pkg/tools/release_gate.R",
                "pkg/cran-comments.md", "pkg/inst/notes.md"),
        good = c("pkg/DESCRIPTION", "pkg/R/a.R", "pkg/inst/doc/a.html",
                 "pkg/inst/eq26_correction.md", "pkg/build/vignette.rds",
                 "pkg/NEWS.md", "pkg/man/figures/logo.png"))
  ex <- file.path(tempdir(), "release_checks_untar")
  unlink(ex, recursive = TRUE)
  dir.create(ex)
  utils::untar(tarball, exdir = ex)
  files <- list.files(ex, recursive = TRUE, all.files = TRUE)
  check(sprintf("%s carries no private strings (names and content)", what),
        function(x) private_in_files(x$files, x$root),
        list(files = files, root = ex),
        bad = list(files = "bad.txt", root = control_dir),
        good = list(files = "good.txt", root = control_dir),
        n_bad = length(private_patterns))
  unlink(ex, recursive = TRUE)
  size_mb <- file.size(tarball) / 1048576
  report(size_mb < 5, sprintf("%s is %.2f MB, under the 5 MB guideline", what, size_mb))
}

# ---- Artifact Mode: A Built Tarball and the History Behind It --------------
# With --history the tarball is left out and the history is swept alone.
if (length(artifact) || history_only) {
  if (length(artifact)) inspect_tarball(artifact, "the tarball")

  # Every ref that can be pushed: the current branch, every local branch with
  # an upstream, and the tags. A private branch with no upstream and the
  # generated site branch are left out.
  git <- function(...) suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = FALSE))
  inside <- identical(git("rev-parse", "--is-inside-work-tree"), "true")
  if (!inside) {
    report(FALSE, "the history that can be pushed carries no private strings",
           "git did not report a work tree here, so the history was not read")
  } else {
    heads <- git("for-each-ref", shQuote("--format=%(refname:short)|%(upstream:short)"), "refs/heads")
    tracked <- sub("[|].*", "", heads[!grepl("[|]$", heads)])
    refs <- unique(c(git("rev-parse", "--abbrev-ref", "HEAD"), tracked,
                     git("tag", "--list")))
    refs <- setdiff(refs, c("gh-pages", "HEAD", ""))
    if (!length(refs)) refs <- "HEAD"
    log <- git("log", shQuote(refs), "-p", "--no-color", "--no-ext-diff",
               shQuote("--format=commit %H%nauthor %an <%ae>%ncommitter %cn <%ce>%n%B"))
    history_hits <- function(lines) {
      if (!length(lines)) return("git log returned nothing: the history was not read")
      at <- grep("^commit [0-9a-f]{40}$", lines)
      if (!length(at)) return("no commit found in the git log output")
      lines <- lines[at[1L]:length(lines)]
      at <- at - at[1L] + 1L
      of <- findInterval(seq_along(lines), at)
      commit <- substr(lines[at], 8L, 14L)[of]
      # The two lines after "commit" name the author and the committer, with
      # their addresses. A line belongs to a commit's diff once that commit's
      # first "diff --git" line has passed; before it, to the message.
      diffs <- grep("^diff --git ", lines)
      in_diff <- findInterval(seq_along(lines), diffs) > 0L &
        c(0L, diffs)[findInterval(seq_along(lines), diffs) + 1L] > at[of]
      where <- ifelse((seq_along(lines) - at[of]) %in% 1:2, "the author or committer line",
                      ifelse(in_diff, "a diff", "the message"))
      hits <- private_in_lines(lines, "line")
      n <- as.integer(sub("^line:([0-9]+) .*", "\\1", hits))
      unique(sprintf("commit %s, %s %s", commit[n], where[n], sub("^line:[0-9]+ ", "", hits)))
    }
    control_commit <- c(paste("commit", strrep("a", 40)), "author A B <a@b.org>", "committer A B <a@b.org>",
                        "A change", "diff --git a/x b/x")
    check(sprintf("the history that can be pushed carries no private strings (%s)",
                  paste(refs, collapse = ", ")),
          history_hits, log,
          bad = c(control_commit, paste0("+", private_bad_line),
                  paste("commit", strrep("b", 40)),
                  paste0("author A B <noreply@", private_patterns[[2L]], ".com>"),
                  "committer A B <a@b.org>", "Another change"),
          good = c(control_commit, "+x <- 1"), n_bad = 2L)

    # Who made the commits: the name and the address of every author, and of
    # every committer. An author is a person of Authors@R, under an address
    # given there. A committer is one of them too, or the hosting service's
    # own identity, which it uses for a merge made on its web pages. Anyone
    # else is printed and fails the check, and the address is part of the
    # test: a commit under an author's name with another address is public
    # under that address once it is pushed.
    people <- eval(parse(text = field("Authors@R")))
    known <- vapply(people, function(p) paste(c(p$given, p$family), collapse = " "), "")
    known_mail <- tolower(unlist(lapply(people, function(p) p$email)))
    host_identity <- "GitHub <noreply@github.com>"
    foreign <- function(x) {
      who <- unique(grep(paste0("^", x$role, " "), x$lines, value = TRUE))
      if (!length(who)) return(sprintf("no %s line found in the git log output", x$role))
      who <- sub(paste0("^", x$role, " "), "", who)
      name <- trimws(sub(" *<[^<>]*>$", "", who))
      mail <- tolower(ifelse(grepl("<[^<>]*>$", who), sub("^.*<([^<>]*)>$", "\\1", who), ""))
      ok <- (name %in% x$names & mail %in% x$mail) | (x$role == "committer" & who == host_identity)
      sprintf("%s %s: %s", x$role, who[!ok],
              ifelse(!name[!ok] %in% x$names, "not a person of Authors@R",
                     "the address is not one given in Authors@R"))
    }
    for (role in c("author", "committer")) {
      seen <- unique(sub(paste0("^", role, " "), "", grep(paste0("^", role, " "), log, value = TRUE)))
      check(sprintf("every %s is a package author under an address of Authors@R%s (%s)", role,
                    if (role == "committer") ", or the hosting service" else "",
                    paste(seen, collapse = "; ")),
            foreign, list(role = role, lines = log, names = known, mail = known_mail),
            bad = list(role = role, names = "Ada Bell", mail = "ada@example.org",
                       lines = paste(role, c("Somebody Else <ada@example.org>",
                                             "Ada Bell <noreply@elsewhere.example>", "Ada Bell"))),
            good = list(role = role, names = "Ada Bell", mail = "ada@example.org",
                        lines = paste(role, c("Ada Bell <Ada@Example.org>",
                                              if (role == "committer") host_identity))),
            n_bad = 3L)
    }
  }

  cat(sprintf("\n%d check(s), %d failed.\n", state$run, state$failed))
  quit(status = if (state$failed > 0L) 1L else 0L)
}

# ---- Reading the Sources ----------------------------------------------------
r_files    <- list.files("R", pattern = "[.][Rr]$", full.names = TRUE)
rd_files   <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
# Every R Markdown source under vignettes/, the site-only articles in its
# subfolders included: the text detectors read them all. The vignettes that
# are installed are the ones at the top of vignettes/ (the build leaves the
# subfolders out); the checks of titles, index entries, and vignette("name")
# are about those.
vig_files  <- list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE, recursive = TRUE)
vig_installed <- list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE)
test_files <- list.files("tests", pattern = "[.][Rr]$", full.names = TRUE, recursive = TRUE)
prose_files <- c(r_files, rd_files, vig_files, test_files,
                 Filter(file.exists, c("DESCRIPTION", "NEWS.md", "README.md",
                                       "inst/CITATION", "CITATION.cff",
                                       ".zenodo.json", "cran-comments.md")))

read_all <- function(files) {
  x <- lapply(files, readLines, warn = FALSE, encoding = "UTF-8")
  names(x) <- files
  x
}
lines_of <- read_all(prose_files)

# Lines of a named list of files that match a pattern, as "file:line".
grep_lines <- function(pattern, files, ...) {
  out <- character()
  for (f in names(files)) {
    m <- grep(pattern, files[[f]], ...)
    if (length(m)) out <- c(out, sprintf("%s:%d", f, m))
  }
  out
}

# A phrase that wraps across lines is still found: each file is searched as
# one string with its white space collapsed and its roxygen markers removed,
# and a match is reported as "file:line" for the line it begins on. A fixed
# string in `except` (a cited title, say) is blanked before the search.
grep_collapsed <- function(pattern, files, except = character(), ...) {
  out <- character()
  for (f in names(files)) {
    l <- gsub("\\s+", " ", trimws(sub("^#' ?", "", files[[f]])))
    kept <- which(nzchar(l))
    if (!length(kept)) next
    starts <- cumsum(c(1L, nchar(l[kept])[-length(kept)] + 1L))
    text <- paste(l[kept], collapse = " ")
    for (x in except) text <- gsub(x, strrep("_", nchar(x)), text, fixed = TRUE)
    m <- gregexpr(pattern, text, ...)[[1L]]
    if (m[1L] > 0L) out <- c(out, sprintf("%s:%d", f, unique(kept[findInterval(m, starts)])))
  }
  out
}

# The parse data of an R source file (or of text, for a control): terminal
# tokens in order, each with the name of the top-level function it sits in.
# Comments and roxygen never reach a detector.
read_source <- function(path = NULL, text = NULL, label = path) {
  exprs <- if (is.null(text)) parse(path, keep.source = TRUE) else
    parse(text = text, keep.source = TRUE)
  pd <- utils::getParseData(exprs)
  if (is.null(pd) || !nrow(pd)) return(NULL)
  refs <- attr(exprs, "srcref")
  first <- vapply(refs, function(r) r[1L], integer(1))
  last  <- vapply(refs, function(r) r[3L], integer(1))
  is_named <- function(e) is.call(e) && length(e) == 3L && as.character(e[[1L]])[1L] %in% c("<-", "=") &&
    is.name(e[[2L]])
  top_name <- vapply(seq_along(exprs), function(k)
    if (is_named(exprs[[k]])) as.character(exprs[[k]][[2L]]) else "", character(1))
  # The argument names of each top-level function, and the line it begins on.
  is_fn <- vapply(seq_along(exprs), function(k)
    is_named(exprs[[k]]) && is.call(exprs[[k]][[3L]]) && identical(exprs[[k]][[3L]][[1L]], as.name("function")),
    logical(1))
  fns <- lapply(which(is_fn), function(k)
    list(name = top_name[k], line = first[k], formals = names(exprs[[k]][[3L]][[2L]])))
  tk <- pd[pd$terminal & pd$token != "COMMENT", ]
  tk <- tk[order(tk$line1, tk$col1), ]
  tk$top <- vapply(tk$line1, function(l) {
    k <- which(first <= l & l <= last)
    if (length(k)) k[1L] else NA_integer_
  }, integer(1))
  tk$fn <- ifelse(is.na(tk$top), "", top_name[tk$top])
  rownames(tk) <- NULL
  list(label = label, tk = tk, pd = pd, fns = fns)
}
sources <- Filter(Negate(is.null), lapply(r_files, read_source))
control_source <- function(text) list(read_source(text = text, label = "control"))

# A finding as "file:line: what in function".
at <- function(src, rows, what) {
  if (!length(rows)) return(character())
  fn <- src$tk$fn[rows]
  sprintf("%s:%d: %s in %s", src$label, src$tk$line1[rows], what, ifelse(nzchar(fn), fn, "top-level code"))
}
over_sources <- function(per_source) function(srcs) unlist(lapply(srcs, per_source))

# What the NAMESPACE file registers. A print or format method is known by its
# registration and not by its name: a function called print_progress is not
# a print method.
namespace <- tryCatch(parseNamespaceFile(basename(root), dirname(root)), error = function(e) NULL)
s3_methods <- if (is.null(namespace)) matrix(NA_character_, 0L, 4L) else namespace$S3methods
s3_function <- ifelse(is.na(s3_methods[, 3L]), paste(s3_methods[, 1L], s3_methods[, 2L], sep = "."),
                      s3_methods[, 3L])
print_methods <- unname(s3_function[s3_methods[, 1L] == "print"])
display_methods <- unname(s3_function[s3_methods[, 1L] %in% c("print", "format", "toString")])
if (is.null(namespace))
  report(FALSE, "the NAMESPACE file can be read",
         "it could not be parsed, so no function is known to be a print method")

# ---- What CRAN Forbids in Package Code -------------------------------------
# One definition of each pattern: the detectors are the ones that
# tests/testthat/test-rd_hygiene.R runs on the installed namespace, on CRAN
# as well as locally. That file is evaluated here down to its marker line.
# Here the detectors read the source files, one top-level expression at a
# time, so they name the file and the line, they see code that is outside
# any function, and they run before anything is installed.
hygiene_file <- file.path("tests", "testthat", "test-rd_hygiene.R")
hygiene <- local({
  if (!file.exists(hygiene_file)) return(NULL)
  lines <- readLines(hygiene_file, warn = FALSE)
  marker <- grep("^# ---- Reading the Package", lines)
  if (length(marker) != 1L) return(NULL)
  env <- new.env(parent = globalenv())
  tryCatch({ eval(parse(text = lines[seq_len(marker - 1L)], keep.source = FALSE), env); env },
           error = function(e) NULL)
})
if (is.null(hygiene))
  report(FALSE, "the detectors of tests/testthat/test-rd_hygiene.R can be read",
         c("the file is missing, does not parse, or has lost the line that begins \"# ---- Reading the Package\";",
           "every check that uses its detectors fails below"))
shared <- function(name) {
  x <- if (is.null(hygiene)) NULL else get0(name, envir = hygiene, inherits = FALSE)
  if (is.null(x)) function(...) stop("the detector could not be read") else x
}
unquote <- shared("unquote")
is_call_to <- shared("is_call_to")
call_args <- shared("call_args")
is_named_arg <- shared("is_named_arg")

# A detector on token tables, applied to each top-level expression of each
# source file. A function named in `exempt` is not examined.
in_sources <- function(detector, exempt = character()) function(srcs) {
  unique(unlist(lapply(srcs, function(s) lapply(split(s$tk, ifelse(is.na(s$tk$top), 0L, s$tk$top)), function(tk) {
    rownames(tk) <- NULL
    fn <- tk$fn[1L]
    if (fn %in% exempt) return(character())
    found <- detector(tk)
    if (!nrow(found)) return(character())
    sprintf("%s:%d: %s in %s", s$label, tk$line1[found$row], found$what,
            if (nzchar(fn)) fn else "top-level code")
  }))))
}

# The generator's state lives in the global environment, and CRAN policy is
# that a package does not touch it. Seeds go through one helper built on the
# withr package; a function that sets the seed itself leaves the caller's
# stream changed.
check("no code touches the generator's state or the global environment",
      in_sources(shared("global_state_hits")), sources,
      bad = control_source(paste0(
        "f <- function(seed) {\n  old <- get(\".Random.seed\", envir = globalenv())\n",
        "  set.seed(seed)\n  stats::rnorm(1)\n}\n",
        "g <- function(x) assign(\"x\", x, pos = 1)\n",
        "RNGkind(\"Mersenne-Twister\")")),
      good = control_source(paste0(
        "# .Random.seed is restored by withr\n",
        "f <- function(seed) withr::with_seed(seed, stats::rnorm(1))\n",
        "g <- function(x, y) graphics::text(x, y, labels = \"a\", pos = 1)")),
      n_bad = 5L)

check("no code uses the superassignment operator",
      in_sources(shared("superassign_hits")), sources,
      bad = control_source(paste0("f <- function() { n <- 0; g <- function() n <<", "- n + 1; g(); n }\n",
                                  "h <- function() { n <- 0; g <- function() assign(\"n\", 1, inherits = TRUE); g(); n }")),
      good = control_source("f <- function() { s <- new.env(); s$n <- 0; s$n <- s$n + 1; s$n }"),
      n_bad = 2L)

check("no code sets the warn option",
      in_sources(shared("warn_option_hits")), sources,
      bad = control_source(paste0(
        "f <- function(x) { old <- options(digits = 3,\n  warn = -1); on.exit(options(old)); as.numeric(x) }\n",
        "g <- function(x) { quiet <- list(warn = -1); do.call(options, quiet); x }")),
      good = control_source(paste0(
        "f <- function(x) { suppressWarnings(as.numeric(x)); getOption(\"warn\") }\n",
        "g <- function(con) readLines(con, warn = FALSE)")),
      n_bad = 2L)

check("no function has a default file path, a save switch, or a path written into a call",
      in_sources(shared("default_path_hits")), sources,
      bad = control_source(paste0(
        "f <- function(x, save = FALSE, filename = \"results.csv\") {\n",
        "  utils::write.csv(x, \"copy.csv\")\n  cat(x, file = \"log.txt\")\n}\n",
        "g <- function(x, path = file.path(tempdir(), \"x.rds\")) x")),
      good = control_source(paste0(
        "f <- function(x, filename = NULL, file, path = TRUE, dir = c(\"upper\", \"lower\"), con = stdout()) {\n",
        "  if (!is.null(filename)) utils::write.csv(x, filename)\n",
        "  cat(\"done\", file = \"\")\n",
        "  system.file(\"extdata\", \"results.rds\", package = \"stats\")\n}")),
      n_bad = 5L)

# CRAN policy on the user's session. par(), options(), and the working
# directory are restored by the function that changes them; a package does
# not install packages, quit, or start other software.
check("no function changes the session without restoring it",
      in_sources(shared("session_state_hits")), sources,
      bad = control_source(paste0(
        "f <- function(x) { old <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(old)); plot(x) }\n",
        "g <- function(x) { graphics::par(mfrow = c(1, 2)); plot(x) }\n",
        "h <- function(x, con) { on.exit(close(con)); options(digits = 3); x }\n",
        "k <- function() quit(save = \"no\")")),
      good = control_source(paste0(
        "f <- function(x) { old <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(old)); plot(x) }\n",
        "g <- function() { q <- function(u) stats::qnorm(u); c(q(0.975), graphics::par(\"usr\"), getOption(\"digits\")) }")),
      n_bad = 3L)

# SMSD's own return from CRAN: an earlier submission was not accepted because
# the sequential functions wrote to the console from inside their bodies. A
# function returns its result; only a print method that the NAMESPACE file
# registers writes.
check("no function writes to the console outside a registered print method",
      function(x) in_sources(shared("console_hits"), exempt = x$print_methods)(x$sources),
      list(sources = sources, print_methods = print_methods),
      bad = list(print_methods = "print.smsd_x", sources = control_source(paste0(
        "seq_x <- function(x) { cat(\"Stopping rule met\\n\"); x }\n",
        "print_progress <- function(x) { message(\"stage \", x); invisible(x) }\n",
        "show_all <- function(x) invisible(lapply(x, print))\n",
        "cat(\"loaded\\n\")"))),
      good = list(print_methods = "print.smsd_x", sources = control_source(paste0(
        "print.smsd_x <- function(x, ...) { cat(format(x), sep = \"\\n\"); invisible(x) }\n",
        "seq_x <- function(x) { if (x < 0) warning(\"negative\"); utils::capture.output(print(x)) }"))),
      n_bad = 4L)

# ---- House Conventions in the Code ------------------------------------------
# A seed is the caller's to give. No function has one as a default (under
# any argument name that says seed), assigns one to `seed`, or hands a number
# to the seed helpers.
seeded <- c("with_seed", "local_seed", ".with_seed", "with_preserve_seed", "set.seed")
is_number_token <- function(tk, rows)
  length(rows) == 1L && tk$token[rows] == "NUM_CONST" && grepl("^[0-9.]", tk$text[rows])
check("no seed is a constant, as a default or in a function body",
      over_sources(function(s) {
        tk <- s$tk; n <- nrow(tk)
        i <- which(tk$token == "SYMBOL_FORMALS" & grepl("seed", tk$text, ignore.case = TRUE))
        i <- i[i + 2L <= n]
        i <- i[tk$token[i + 1L] == "EQ_FORMALS" & tk$token[i + 2L] != "NULL_CONST"]
        j <- which(tk$token == "SYMBOL" & grepl("seed", tk$text, ignore.case = TRUE))
        j <- j[j + 2L <= n]
        j <- j[tk$token[j + 1L] %in% c("LEFT_ASSIGN", "EQ_ASSIGN") &
                 vapply(j + 2L, function(r) is_number_token(tk, r), logical(1))]
        out <- c(at(s, i, sprintf("%s has a default other than NULL", tk$text[i])),
                 at(s, j, sprintf("%s is assigned a number", tk$text[j])))
        for (k in which(is_call_to(tk, seeded))) {
          args <- call_args(tk, k)
          for (a in seq_along(args)) {
            named <- is_named_arg(tk, args[[a]])
            value <- if (named) args[[a]][-(1:2)] else args[[a]]
            if ((if (named) tk$text[args[[a]][1L]] == "seed" else a == 1L) && is_number_token(tk, value))
              out <- c(out, at(s, k, sprintf("%s() is given the number %s", tk$text[k], tk$text[value])))
          }
        }
        out
      }),
      sources,
      bad = control_source(paste0(
        "f <- function(G = 1000, seed = 113) G\n",
        "g <- function(G = 10, seed = NULL) { if (is.null(seed)) seed <- 113; withr::with_seed(seed, stats::rnorm(G)) }\n",
        "h <- function(G = 10) withr::with_seed(113, stats::rnorm(G))\n",
        "k <- function(G = 10) .with_seed(function() stats::rnorm(G), seed = 2026L)\n",
        "m <- function(G = 10, random_seed = 113) G")),
      good = control_source(paste0(
        "f <- function(G = 1000, seed = NULL) { set <- 113; .with_seed(seed, function() stats::rnorm(G)) }")),
      n_bad = 5L)

# A function the user calls takes conf_level and never alpha (inside a
# function, alpha is 1 - conf_level and may be handed to a helper), and no
# argument name of any function has a dot other than na.rm. A method keeps
# the argument names of its generic, as as.data.frame() has row.names; where
# the generic cannot be found, the method's dotted names are not judged.
generic_formals <- function(fn) {
  i <- match(fn, s3_function)
  if (is.na(i)) return(character())
  generic <- tryCatch(if (is.na(s3_methods[i, 4L])) get(s3_methods[i, 1L], mode = "function")
                      else getExportedValue(s3_methods[i, 4L], s3_methods[i, 1L]),
                      error = function(e) NULL)
  if (is.null(generic)) NA_character_ else as.character(names(formals(generic)))
}
user_facing <- c(if (!is.null(namespace)) namespace$exports, unname(s3_function))
control_generic <- function(fn) if (fn == "as.data.frame.smsd_x") names(formals(as.data.frame)) else character()
check("no exported function has an alpha argument, and no argument name has a dot other than na.rm",
      function(x) unlist(lapply(x$sources, function(s) lapply(s$fns, function(f) {
        allowed <- x$generic_formals(f$name)
        dotted <- f$formals[grepl(".", f$formals, fixed = TRUE) & !f$formals %in% c("...", "na.rm")]
        dotted <- if (anyNA(allowed)) character() else setdiff(dotted, allowed)
        named <- c(if (f$name %in% x$user_facing) intersect(f$formals, "alpha"), dotted)
        if (length(named)) sprintf("%s:%d: %s has the argument %s", s$label, f$line, f$name, named)
      }))),
      list(sources = sources, generic_formals = generic_formals, user_facing = user_facing),
      bad = list(generic_formals = control_generic, user_facing = c("f", "as.data.frame.smsd_x"),
                 sources = control_source(paste0(
                   "f <- function(x, alpha = 0.05) x\n",
                   ".g <- function(x, conf.level = 0.95) x\n",
                   "as.data.frame.smsd_x <- function(x, row.names = NULL, optional = FALSE, keep.attributes = TRUE, ...) x"))),
      good = list(generic_formals = control_generic, user_facing = c("f", "as.data.frame.smsd_x"),
                  sources = control_source(paste0(
                    "f <- function(x, conf_level = 0.95, na.rm = FALSE, ...) .z(1 - conf_level)\n",
                    ".z <- function(alpha) stats::qnorm(1 - alpha / 2)\n",
                    "as.data.frame.smsd_x <- function(x, row.names = NULL, optional = FALSE, ...) x"))),
      n_bad = 3L)

# Stored values are never rounded; rounding belongs to display. What is
# watched is a number rounded for show: round() with a digits argument,
# signif(), a number turned into text and back (as.numeric(format(x, digits = 3))), and
# trunc(x * 1000) / 1000 and its relatives. round(x) to a whole number is
# arithmetic on a count, as ceiling() is, and formatC() returns text; neither
# is watched.
#
# Each such call is followed to where its value goes. It is display code,
# and passes, when the value ends in text, a condition message, a plot label,
# or a comparison (the whole number test x == round(x, 0) gives TRUE or
# FALSE, not a rounded number), when it is discarded, or when the function is
# a registered print or format method. It fails when the rounded number is
# returned by a function or kept in an object that is. A value given to a
# local variable is followed through that variable's later uses, and a value
# returned by a local function through that function's calls.
absorbers <- c("stop", "warning", "message", "sprintf", "gettextf", "paste", "paste0",
               "format", "formatC", "prettyNum", "toString", "cat", "print", "nchar",
               "all.equal", "identical", "isTRUE", "isFALSE", "match",
               "text", "mtext", "title", "legend", "axis",
               ".smsd_stop", ".smsd_describe")
comparisons <- c("EQ", "NE", "LT", "GT", "LE", "GE")
rounding_hits <- function(s, display_fns) {
  pd <- s$pd
  pd <- pd[order(pd$line1, pd$col1, -pd$line2, -pd$col2), ]
  ids <- pd$id; tok <- pd$token; txt <- pd$text
  row_of <- function(id) match(id, ids)
  parent_of <- function(id) pd$parent[row_of(id)]
  kids_by_parent <- split(ids, pd$parent)
  kids <- function(id) { k <- kids_by_parent[[as.character(id)]]; if (is.null(k)) integer() else k }
  kid_tokens <- function(id) tok[row_of(kids(id))]
  expr_kids <- function(id) { k <- kids(id); k[!pd$terminal[row_of(k)]] }
  start <- pd$line1 * 1e6 + pd$col1
  end <- pd$line2 * 1e6 + pd$col2
  inside <- function(id) which(start >= start[row_of(id)] & end <= end[row_of(id)])
  fn_nodes <- unique(pd$parent[tok %in% c("FUNCTION", "'\\\\'")])
  assign_nodes <- unique(pd$parent[tok %in% c("LEFT_ASSIGN", "RIGHT_ASSIGN", "EQ_ASSIGN")])
  compare_nodes <- unique(pd$parent[tok %in% comparisons | (tok == "SPECIAL" & txt == "%in%")])
  loop_nodes <- unique(pd$parent[tok %in% c("FOR", "WHILE", "REPEAT")])
  # The function a call node calls, or NA when the node is not a call.
  call_name <- function(id) {
    k <- kids(id)
    if (length(k) < 2L || tok[row_of(k[2L])] != "'('") return(NA_character_)
    head <- row_of(kids(k[1L]))
    name <- txt[head][tok[head] == "SYMBOL_FUNCTION_CALL"]
    if (length(name) == 1L) name else NA_character_
  }
  n_args <- function(id) sum(kid_tokens(id) == "','") + (length(kids(id)) > 3L)
  enclosing <- function(id, nodes) {
    p <- parent_of(id)
    while (p != 0L && !p %in% nodes) p <- parent_of(p)
    p
  }
  # The value and the target of an assignment node.
  sides <- function(a) {
    k <- expr_kids(a)
    if ("RIGHT_ASSIGN" %in% kid_tokens(a)) list(value = k[1L], target = k[length(k)])
    else list(value = k[length(k)], target = k[1L])
  }
  first_symbol <- function(e) {
    r <- inside(e)
    r <- r[tok[r] == "SYMBOL"]
    if (length(r)) txt[r[1L]] else NA_character_
  }
  seen <- new.env(parent = emptyenv())

  # Where the value of `node` ends up: character(0) when it goes no further
  # than display, otherwise what keeps it.
  flows <- function(node, depth = 0L) {
    if (depth > 40L) return("goes further than this check can follow")
    out <- character()
    repeat {
      p <- parent_of(node)
      if (p == 0L) return(c(out, "is kept at the top level of the file"))
      if (p %in% compare_nodes) return(out)
      called <- call_name(p)
      if (!is.na(called)) {
        if (called %in% absorbers) return(out)
        if (called == "return") return(c(out, result_of(enclosing(p, fn_nodes), depth + 1L)))
      }
      if (p %in% assign_nodes) {
        name <- first_symbol(sides(p)$target)
        scope <- enclosing(p, fn_nodes)
        if (scope == 0L)
          return(c(out, sprintf("is kept in %s at the top level of the file", name)))
        if (!is.na(name)) out <- c(out, uses_of(name, scope, p, depth + 1L))
        node <- p
        next
      }
      kt <- kid_tokens(p)
      if ("'{'" %in% kt) {
        ek <- expr_kids(p)
        if (node != ek[length(ek)]) return(out)
      } else if ("IF" %in% kt) {
        if (node == expr_kids(p)[1L]) return(out)
      } else if (any(c("FOR", "WHILE", "REPEAT") %in% kt)) {
        return(out)
      } else if (tok[row_of(p)] == "forcond") {
        r <- row_of(kids(p))
        return(c(out, uses_of(txt[r][tok[r] == "SYMBOL"][1L], enclosing(p, fn_nodes), p, depth + 1L)))
      } else if (p %in% fn_nodes) {
        ek <- expr_kids(p)
        if (node == ek[length(ek)]) return(c(out, result_of(p, depth + 1L)))
        r <- row_of(kids(p))
        formal <- txt[r][tok[r] == "SYMBOL_FORMALS" & start[r] < start[row_of(node)]]
        return(c(out, uses_of(formal[length(formal)], p, p, depth + 1L)))
      }
      node <- p
    }
  }
  # The later uses of a variable that holds the value, within its function.
  # Inside a loop an earlier line runs again, so every use counts there.
  uses_of <- function(name, scope, from, depth) {
    if (is.na(name) || scope == 0L) return(character())
    key <- paste(name, scope)
    if (!is.null(seen[[key]])) return(character())
    seen[[key]] <- TRUE
    r <- inside(scope)
    r <- r[txt[r] == name & tok[r] %in% c("SYMBOL", "SYMBOL_FUNCTION_CALL")]
    looped <- enclosing(from, c(loop_nodes, scope)) != scope
    if (!looped) r <- r[start[r] >= start[row_of(from)]]
    out <- character()
    for (i in r) {
      e <- pd$parent[i]
      if (tok[i] == "SYMBOL_FUNCTION_CALL") {
        out <- c(out, flows(parent_of(e), depth))
        next
      }
      up <- parent_of(e)
      if (up %in% assign_nodes && identical(sides(up)$target, e)) next
      out <- c(out, flows(e, depth))
    }
    out
  }
  # Where the value of a function goes when the function returns the value.
  result_of <- function(f, depth) {
    if (f == 0L) return("is kept at the top level of the file")
    p <- parent_of(f)
    if (p %in% assign_nodes && identical(sides(p)$value, f)) {
      k <- kids(sides(p)$target)
      name <- if (length(k) == 1L && tok[row_of(k)] == "SYMBOL") txt[row_of(k)] else NA_character_
      if (parent_of(p) == 0L) {
        if (!is.na(name) && name %in% display_fns) return(character())
        return(sprintf("is returned by %s", if (is.na(name)) "a function" else name))
      }
      if (!is.na(name)) return(c(uses_of(name, enclosing(p, fn_nodes), p, depth), flows(p, depth)))
    }
    flows(f, depth)
  }
  power_of_ten <- function(e) {
    r <- row_of(kids(e))
    if (length(r) == 1L && tok[r] == "NUM_CONST") {
      v <- suppressWarnings(as.numeric(sub("L$", "", txt[r])))
      return(!is.na(v) && v >= 10 && abs(log10(v) - round(log10(v))) < 1e-9)
    }
    length(r) == 3L && tok[r[2L]] == "'^'" && identical(txt[row_of(kids(kids(e)[1L]))], "10")
  }

  out <- character()
  for (i in which(tok == "SYMBOL_FUNCTION_CALL")) {
    name <- txt[i]
    call <- parent_of(pd$parent[i])
    label <- NULL
    if (name == "signif" || (name == "round" && n_args(call) >= 2L)) {
      label <- paste0(name, "()")
    } else if (name %in% c("as.numeric", "as.double", "as.integer")) {
      arg <- expr_kids(call)[2L]
      inner <- call_name(arg)
      if (!is.na(inner) && inner %in% c("format", "formatC", "sprintf", "prettyNum", "toString")) {
        # A part of a date (format(x, "%Y")) and a whole number written out
        # (sprintf("%d", k)) lose nothing; a precision does.
        r <- inside(arg)
        strings <- txt[r][tok[r] == "STR_CONST"]
        harmless <- (inner == "format" && any(grepl("%[A-Za-z]", strings))) ||
          (inner == "sprintf" && length(strings) > 0L && !any(grepl("%[-+ 0#]*[0-9]*[.][0-9*]", strings)))
        if (!harmless) label <- sprintf("%s(%s())", name, inner)
      }
    } else if (name %in% c("trunc", "floor", "ceiling", "round")) {
      p <- parent_of(call)
      k <- kids(p)
      if (p != 0L && "'/'" %in% kid_tokens(p) && identical(k[1L], call) && power_of_ten(k[length(k)]) &&
          "'*'" %in% kid_tokens(expr_kids(call)[2L])) {
        label <- sprintf("%s() of a multiple, divided back", name)
        call <- p
      }
    }
    if (is.null(label)) next
    rm(list = ls(seen), envir = seen)
    why <- unique(flows(call))
    if (length(why)) out <- c(out, sprintf("%s:%d: the number from %s %s", s$label, pd$line1[i], label, why[1L]))
  }
  out
}
check("no stored value is rounded (rounding only where the result is displayed or compared)",
      function(x) unlist(lapply(x$sources, rounding_hits, display_fns = x$display)),
      list(sources = sources, display = display_methods),
      bad = list(display = "format.smsd_x", sources = control_source(paste0(
        "a <- function(x) { out <- data.frame(term = \"est\", value = round(mean(x), 3)); out }\n",
        "b <- function(x) signif(mean(x), 4)\n",
        "d <- function(x) round(x, digits = 0)\n",
        "e <- function(x) { est <- as.numeric(format(mean(x), digits = 3)); list(est = est) }\n",
        "f <- function(x) as.numeric(sprintf(\"%.3f\", mean(x)))\n",
        "g <- function(x) trunc(mean(x) * 1000) / 1000\n",
        "h <- function(x) { est <- round(mean(x), 3); if (est == 0) 0 else est }\n",
        "summary.smsd_x <- function(object, ...) { object$value <- signif(object$value, 3); object }\n",
        "k <- function(x) { near <- function(v) round(v, 2); list(est = near(mean(x))) }\n",
        "shown <- round(pi, 2)"))),
      good = list(display = "format.smsd_x", sources = control_source(paste0(
        "a <- function(x) { if (x != round(x, 0) || abs(x - round(x)) > 0) stop(\"got \", round(x, 2)); x }\n",
        "format.smsd_x <- function(x, digits = 3, ...) { shown <- signif(unclass(x), digits); shown }\n",
        # The call sites of this package that an earlier form of this check
        # flagged, each of them correct: a grid of whole sample sizes, text
        # padded for a table, a number formatted inside a sentence, and an
        # amount of money rounded and then formatted.
        "grid <- function(a, b, k) unique(as.integer(round(seq(a, b, length.out = k))))\n",
        "pad <- function(nm, txt) formatC(c(nm, txt), width = max(nchar(c(nm, txt))))\n",
        "report <- function(v, digits) {\n",
        "  num <- function(v) formatC(round(v, digits) + 0, digits = digits, format = \"f\", big.mark = \",\")\n",
        "  paste0(\"The estimate was \", num(v), \".\")\n}\n",
        "money <- function(x, budget) {\n",
        "  places <- max(0, floor(13 - log10(budget)))\n  x <- round(x, places)\n",
        "  if (x == 0) return(\"0\")\n  format(x, digits = 15L)\n}\n",
        "keys <- function(a, b) match(round(a, 2), round(b, 3))\n",
        "share <- function(n, p) round(n * p) / n\n",
        "even <- function(n) 2 * ceiling(n / 2)\n",
        "label <- function(x) { out <- data.frame(value = mean(x)); attr(out, \"label\") <- paste0(round(mean(x), 1)); out }\n",
        "year <- function(day, k) c(as.numeric(format(day, \"%Y\")), as.integer(sprintf(\"%d\", k)))\n",
        "draw <- function(x, y) { graphics::text(x, y, labels = round(y, 1)); invisible(NULL) }"))),
      n_bad = 10L)

# ---- Examples and Help Pages ----------------------------------------------
check("no \\donttest or \\dontrun in the sources or the help pages",
      function(files) grep_lines("\\\\donttest|\\\\dontrun", files),
      lines_of[c(r_files, rd_files)],
      bad = list("control.R" = c("#' @examples", "#' \\donttest{", "#' mean(1:3)", "#' }")),
      good = list("control.R" = c("#' @examples", "#' mean(1:3)")))

# The @examples block that follows each "#' @examples" line of a source file,
# as the example lines and their line numbers in the file.
example_blocks <- function(lines) {
  blocks <- list()
  for (s in grep("^#' @examples", lines)) {
    i <- s + 1L; idx <- integer()
    while (i <= length(lines) && grepl("^#'", lines[i]) && !grepl("^#' @[a-zA-Z]", lines[i])) {
      idx <- c(idx, i); i <- i + 1L
    }
    blocks[[length(blocks) + 1L]] <- list(text = sub("^#' ?", "", lines[idx]), at = idx)
  }
  blocks
}
is_code <- shared("is_code")
# A reviewer's words: "Some code lines in examples are commented out. Please
# never do that." Every window of consecutive comment lines in an example is
# parsed with the markers removed, and a window that parses to a call or an
# assignment is flagged. The rule is mechanical: a line of prose that happens
# to be valid R ("Kelley (2007)", "N - 1") is flagged too and is reworded as
# a sentence, as the top of tests/testthat/test-rd_hygiene.R explains.
check("no commented-out code in any @examples block",
      function(files) {
        out <- character()
        for (f in names(files)) for (b in example_blocks(files[[f]])) {
          cm <- grep("^\\s*#", b$text)
          if (!length(cm)) next
          for (run in split(cm, cumsum(c(1L, diff(cm) != 1L)))) {
            txt <- sub("^\\s*#+\\s?", "", b$text[run]); n <- length(run); flagged <- integer()
            for (a in seq_len(n)) for (z in a:n) {
              w <- paste(txt[a:z], collapse = "\n")
              if (nzchar(trimws(w)) && is_code(w)) flagged <- c(flagged, run[a:z])
            }
            if (length(flagged)) out <- c(out, sprintf("%s:%d", f, b$at[sort(unique(flagged))]))
          }
        }
        out
      },
      lines_of[r_files],
      bad = list("control.R" = c("#' @examples", "#' # The pilot sample size comes first.",
                                 "#' # seq_x_pilot(width = 0.2)", "#' seq_x(1:10)",
                                 "#' # Kelley (2007)", "#' @export")),
      good = list("control.R" = c("#' @examples", "#' # The pilot sample size comes first.",
                                  "#' # Kelley, 2007, Equation 15.",
                                  "#' seq_x(1:10)  # then the check", "#' @export",
                                  "# old <- seq_x(1:10)")),
      n_bad = 2L)

check("no requireNamespace() guard in any @examples block",
      function(files) {
        out <- character()
        for (f in names(files)) for (b in example_blocks(files[[f]])) {
          m <- which(grepl("requireNamespace", b$text, fixed = TRUE) & !grepl("^\\s*#", b$text))
          if (length(m)) out <- c(out, sprintf("%s:%d", f, b$at[m]))
        }
        out
      },
      lines_of[r_files],
      bad = list("control.R" = c("#' @examples", "#' if (requireNamespace(\"withr\")) mean(1:3)")),
      good = list("control.R" = c("#' @examples", "#' # withr is imported, so no requireNamespace() is needed",
                                  "#' mean(1:3)", "f <- function() requireNamespace(\"withr\")")))

# CRAN policy: a package uses at most two cores unless the user asks for
# more, and that covers what CRAN runs as well as what it installs. The
# detector is the test file's, which reads the functions of the namespace.
# Here it reads the source files, and with them the code the test cannot
# see: the test files, the examples, and the R chunks of the vignettes and
# of the README. Each is read with its other lines blanked, so a finding
# names the line of the file. A chunk or an example that does not parse is
# code that is shown and not run, and is left out.
code_only <- function(lines, keep, label) {
  code <- ifelse(keep, lines, "")
  if (!any(nzchar(trimws(code)))) return(NULL)
  whole <- tryCatch(read_source(text = paste(code, collapse = "\n"), label = label), error = function(e) NULL)
  if (!is.null(whole)) return(list(whole))
  runs <- split(which(keep), cumsum(c(1L, diff(which(keep)) != 1L)))
  Filter(Negate(is.null), lapply(runs, function(r) {
    part <- rep("", length(lines)); part[r] <- lines[r]
    tryCatch(read_source(text = paste(part, collapse = "\n"), label = label), error = function(e) NULL)
  }))
}
chunk_sources <- function(lines, label) {
  open <- grep("^\\s*```+\\s*\\{\\s*[rR]\\b", lines)
  close <- grep("^\\s*```+\\s*$", lines)
  keep <- rep(FALSE, length(lines))
  for (o in open) {
    end <- close[close > o][1L]
    if (!is.na(end) && end > o + 1L) keep[(o + 1L):(end - 1L)] <- TRUE
  }
  code_only(lines, keep, label)
}
example_sources <- function(lines, label) {
  keep <- rep(FALSE, length(lines))
  for (b in example_blocks(lines)) keep[b$at] <- TRUE
  code_only(gsub("\\%", "%", sub("^#' ?", "", lines), fixed = TRUE), keep, paste(label, "(an example)"))
}
run_code <- function(files) unlist(lapply(names(files), function(f) {
  if (grepl("[.]Rmd$", f)) chunk_sources(files[[f]], f)
  else if (grepl("^R/", f)) example_sources(files[[f]], f)
  else code_only(files[[f]], rep(TRUE, length(files[[f]])), f)
}), recursive = FALSE)
rmd_files <- c(vig_files, Filter(file.exists, "README.Rmd"))
check("no function, example, test, or vignette takes more than two cores, by default or from the machine",
      function(x) in_sources(shared("core_count_hits"))(c(x$sources, run_code(x$files))),
      list(sources = sources, files = c(lines_of[c(r_files, test_files)], read_all(rmd_files))),
      bad = list(sources = control_source(paste0(
                   "f <- function(x, cores = parallel::detectCores() - 1) parallel::mclapply(x, sqrt, mc.cores = cores)\n",
                   "g <- function(x, workers = 4L) { cl <- parallel::makeCluster(8); on.exit(parallel::stopCluster(cl)); x }")),
                 files = list("vignettes/control.Rmd" = c("The study runs on every core.", "", "```{r study}",
                                                          "out <- parallel::mclapply(1:8, sqrt, mc.cores = 6)", "```"),
                              "R/control.R" = c("#' @examples", "#' options(mc.cores = 4)", "#' f(1:3)", "f <- function(x) x"),
                              "tests/testthat/test-control.R" = "cl <- parallel::makePSOCKcluster(3L)")),
      good = list(sources = control_source(paste0(
                    "f <- function(x, cores = 1L) parallel::mclapply(x, sqrt, mc.cores = cores)\n",
                    "g <- function(x, G = 1000) { cl <- parallel::makeCluster(2L); on.exit(parallel::stopCluster(cl)); x }")),
                  files = list("vignettes/control.Rmd" = c("With mc.cores = 6 the study takes an hour; detectCores() is not used.", "",
                                                           "```{r study}", "out <- parallel::mclapply(1:8, sqrt, mc.cores = 2)", "```", "",
                                                           "```", "mclapply(x, f, mc.cores = 6) and so on, with <your cores>", "```"),
                               "R/control.R" = c("#' @examples", "#' f(1:3)  # two cores at most", "f <- function(x) x"),
                               "tests/testthat/test-control.R" = "x <- rep(1, 8)")),
      n_bad = 6L)

# Parsed, not searched: the markup and the \eqn{ it sits in can be on
# different lines. \code{} inside math becomes \texttt in LaTeX math and stops
# the PDF manual. The detector is the test file's, which also sees a macro
# written with an option, as in \link[stats]{qnorm}.
rd_nodes <- function(x, tags) {
  tag <- attr(x, "Rd_tag")
  hit <- if (!is.null(tag) && tag %in% tags) list(x) else list()
  if (is.list(x)) for (el in x) hit <- c(hit, rd_nodes(el, tags))
  hit
}
rd_text <- function(...) {
  con <- textConnection(paste0(...)); on.exit(close(con)); tools::parse_Rd(con)
}
rd_db <- tools::Rd_db(dir = ".")
math_markup_hits <- shared("math_markup_hits")
check(sprintf("no Rd markup inside \\eqn{} or \\deqn{} (parsed, %d pages)", length(rd_db)),
      function(db) names(db)[vapply(db, function(rd) length(math_markup_hits(rd)) > 0L, logical(1))],
      rd_db,
      bad = list(a.Rd = rd_text("\\name{x}\\alias{x}\\title{x}\\description{The width is \\eqn{2 z\n",
                                "\\code{sd} / \\sqrt{n}}.}"),
                 b.Rd = rd_text("\\name{x}\\alias{x}\\title{x}\\description{The quantile is ",
                                "\\eqn{\\link[stats]{qnorm}(0.975)}.}")),
      good = list(control.Rd = rd_text("\\name{x}\\alias{x}\\title{x}\\description{The width is ",
                                       "\\eqn{2 z s / \\sqrt{n}}, with \\code{sd} as \\eqn{s} and ",
                                       "\\eqn{z} from \\code{\\link[stats]{qnorm}}.}")),
      n_bad = 2L)

# Help page titles are AP title case with no period at the end: the first and
# last words and every principal word are capitalized, and only articles,
# conjunctions, and prepositions of three letters or fewer may be lowercase.
# Code, math, single symbols, and identifiers keep their own case.
may_be_lower <- c("a", "an", "the", "and", "but", "or", "nor", "for", "yet", "so",
                  "as", "if", "at", "by", "in", "of", "off", "on", "per", "to", "up",
                  "via", "vs",
                  # particles of surnames, lowercase by convention
                  "de", "van", "von", "der", "den", "la", "le", "du", "di", "da", "del")
title_problems <- function(title) {
  t <- trimws(gsub("\\s+", " ", title))
  out <- character()
  if (grepl("[.]$", t)) out <- c(out, "ends with a period")
  words <- strsplit(t, " ", fixed = TRUE)[[1L]]
  words <- gsub("^[\\[(\"'`*_]+|[])\"'`*_,:;.!?]+$", "", words)
  words <- words[nzchar(words)]
  wordlike <- grepl("^[A-Za-z][A-Za-z'-]*$", words) & nchar(sub("-.*", "", words)) > 1L
  first_part <- sub("-.*", "", words)
  # A word with a capital anywhere (AIPE, dMACS) keeps the case it has.
  capital <- grepl("[A-Z]", first_part)
  edge <- seq_along(words) %in% c(1L, length(words))
  wrong <- wordlike & !capital & (edge | !tolower(first_part) %in% may_be_lower)
  if (any(wrong)) out <- c(out, sprintf("not capitalized: %s", paste(words[wrong], collapse = ", ")))
  out
}
# The text of a help page's title with marked-up code and math replaced by a
# placeholder, so their case is not judged.
rd_title <- function(rd) {
  flat <- function(x) {
    tag <- attr(x, "Rd_tag")
    if (!is.null(tag) && tag %in% c("\\code", "\\eqn", "\\deqn", "\\pkg", "\\link", "\\verb",
                                   "\\samp", "\\var", "\\file", "\\env", "\\url", "\\emph"))
      return(" CODE ")
    if (is.list(x)) paste(vapply(x, flat, character(1)), collapse = "") else paste(x, collapse = "")
  }
  node <- rd_nodes(rd, "\\title")
  if (length(node)) flat(node[[1L]]) else ""
}
check("help page titles are AP title case with no period at the end",
      function(db) unlist(lapply(names(db), function(p) {
        pr <- title_problems(rd_title(db[[p]]))
        if (length(pr)) sprintf("%s: \"%s\": %s", p, trimws(gsub("\\s+", " ", rd_title(db[[p]]))),
                                paste(pr, collapse = "; "))
      })),
      rd_db,
      bad = list(a.Rd = rd_text("\\name{a}\\alias{a}\\title{Sequential AIPE check for the standardized mean difference}"),
                 b.Rd = rd_text("\\name{b}\\alias{b}\\title{Pilot Sample Size.}")),
      good = list(a.Rd = rd_text("\\name{a}\\alias{a}\\title{Sequential AIPE Check for the Standardized Mean ",
                                 "Difference, With \\code{conf_level} and Cohen's \\emph{d} (Mukhopadhyay and de Silva)}"),
                  b.Rd = rd_text("\\name{b}\\alias{b}\\title{Methods for smsd_check Objects: What a Check Is Made Of}")))

# ---- Versions and Citation Files -------------------------------------------
version <- field("Version")
check(sprintf("the version %s has exactly three components", version),
      function(v) if (grepl("^[0-9]+[.][0-9]+[.][0-9]+$", v)) character() else v,
      version, bad = "1.0.0.9000", good = "1.0.0")

news_version <- function(lines) sub("^# +[A-Za-z]+ +", "", grep("^# ", lines, value = TRUE)[1L])
lockstep <- function(x) {
  off <- x$versions[is.na(x$versions) | x$versions != x$description]
  sprintf("%s has %s; DESCRIPTION has %s", names(off), off, x$description)
}
yaml_ok <- requireNamespace("yaml", quietly = TRUE)
json_ok <- requireNamespace("jsonlite", quietly = TRUE)
cff <- if (yaml_ok && file.exists("CITATION.cff")) yaml::read_yaml("CITATION.cff") else NULL
cit_meta <- as.list(desc[1L, ])
cit <- tryCatch(utils::readCitationFile("inst/CITATION", meta = cit_meta), error = function(e) e)
cit_manual <- if (inherits(cit, "error")) NULL else
  Filter(function(b) identical(tolower(b$bibtype), "manual"), cit)
cit_note <- if (length(cit_manual)) cit_manual[[1L]]$note else NA_character_
# The letter to CRAN names the version it submits, in its Submission section
# ("version 1.0.0"). A letter left over from the release before would name
# the old one.
letter_version <- function(lines) {
  at <- grep("^## +Submission *$", lines)
  if (!length(at)) return(NA_character_)
  rest <- lines[-seq_len(at[1L])]
  ends <- grep("^## ", rest)
  if (length(ends)) rest <- rest[seq_len(ends[1L] - 1L)]
  text <- paste(rest, collapse = " ")
  m <- regmatches(text, regexpr("\\bversion +[0-9]+([.][0-9]+)+", text))
  if (length(m)) sub("^version +", "", m) else NA_character_
}
letter_ok <- identical(letter_version(c("## Submission", "", "This is the first release of pkg,", "version 1.0.1.",
                                        "", "## Results", "It replaces version 0.9.0.")), "1.0.1") &&
  is.na(letter_version(c("## Submission", "", "A first release.", "", "## Results", "version 1.0.1")))
if (!letter_ok)
  report(FALSE, "the version in the letter to CRAN can be read",
         "the reader of cran-comments.md failed its own control, so the letter's version is not compared")
check("the version is the same in DESCRIPTION, NEWS.md, CITATION.cff, inst/CITATION, and cran-comments.md",
      lockstep,
      list(description = version,
           versions = c("the top heading of NEWS.md" = news_version(lines_of[["NEWS.md"]]),
                        "CITATION.cff" = if (is.null(cff)) NA_character_ else as.character(cff$version),
                        "the package entry of inst/CITATION" =
                          if (is.na(cit_note)) NA_character_ else sub("^.*version +", "", cit_note),
                        if (file.exists("cran-comments.md"))
                          c("the Submission section of cran-comments.md" =
                              letter_version(lines_of[["cran-comments.md"]])))),
      bad = list(description = "1.0.1", versions = c(NEWS = "1.0.0", cff = "1.0.1", letter = NA_character_)),
      good = list(description = "1.0.1", versions = c(NEWS = "1.0.1", cff = "1.0.1", letter = "1.0.1")),
      n_bad = 2L)

people <- eval(parse(text = field("Authors@R")))
desc_family <- vapply(people, function(p) paste(p$family, collapse = " "), "")
desc_given  <- vapply(people, function(p) paste(p$given, collapse = " "), "")
desc_orcid  <- vapply(people, function(p) {
  o <- p$comment[["ORCID"]]
  if (is.null(o)) "" else sub("^https?://orcid.org/", "", o)
}, "")
same_people <- function(x) {
  out <- character()
  if (!identical(unname(x$family), unname(desc_family)) || !identical(unname(x$given), unname(desc_given)))
    out <- c(out, sprintf("authors are [%s]; DESCRIPTION has [%s]",
                          paste(x$given, x$family, collapse = "; "),
                          paste(desc_given, desc_family, collapse = "; ")))
  if (!is.null(x$orcid) && !identical(unname(x$orcid), unname(desc_orcid)))
    out <- c(out, sprintf("ORCIDs are [%s]; DESCRIPTION has [%s]",
                          paste(x$orcid, collapse = "; "), paste(desc_orcid, collapse = "; ")))
  out
}
people_bad <- list(family = c("Else", desc_family), given = c("Somebody", desc_given),
                   orcid = c("", desc_orcid))
people_good <- list(family = desc_family, given = desc_given, orcid = desc_orcid)
if (is.null(cff)) {
  report(FALSE, "CITATION.cff is present and can be read",
         if (!yaml_ok) "the yaml package is not installed, so the file was not read" else "no CITATION.cff")
} else {
  check("CITATION.cff names the authors of DESCRIPTION, in its order, with their ORCIDs",
        same_people,
        list(family = vapply(cff$authors, function(a) a[["family-names"]], ""),
             given = vapply(cff$authors, function(a) a[["given-names"]], ""),
             orcid = vapply(cff$authors, function(a)
               if (is.null(a$orcid)) "" else sub("^https?://orcid.org/", "", a$orcid), "")),
        bad = people_bad, good = people_good)
  urls <- trimws(strsplit(field("URL"), ",")[[1L]])
  check("the repository in CITATION.cff is a URL of DESCRIPTION",
        function(u) if (isTRUE(u %in% urls)) character() else sprintf("repository-code is %s", u),
        cff[["repository-code"]], bad = "https://example.org/elsewhere", good = urls[1L])
}
if (inherits(cit, "error") || !length(cit_manual)) {
  report(FALSE, "inst/CITATION can be read by utils::readCitationFile() alone",
         if (inherits(cit, "error")) conditionMessage(cit) else "no entry of type Manual")
} else {
  m <- cit_manual[[1L]]
  check("inst/CITATION names the authors of DESCRIPTION, in its order",
        same_people,
        list(family = vapply(m$author, function(p) paste(p$family, collapse = " "), ""),
             given = vapply(m$author, function(p) paste(p$given, collapse = " "), "")),
        bad = people_bad[c("family", "given")], good = people_good[c("family", "given")])
  date_year <- substr(field("Date"), 1L, 4L)
  check(sprintf("the year in inst/CITATION is the year of the DESCRIPTION Date (%s)", date_year),
        function(y) if (identical(as.character(y), date_year)) character() else sprintf("year is %s", y),
        m$year, bad = "1999", good = date_year)
}
zen <- if (json_ok && file.exists(".zenodo.json"))
  tryCatch(jsonlite::fromJSON(".zenodo.json", simplifyVector = FALSE), error = function(e) NULL)
if (is.null(zen)) {
  report(FALSE, ".zenodo.json is present and can be read",
         if (!json_ok) "the jsonlite package is not installed, so the file was not read" else "missing or not valid JSON")
} else {
  zn <- strsplit(vapply(zen$creators, function(c) c$name, ""), ", *")
  check(".zenodo.json names the authors of DESCRIPTION, in its order, with their ORCIDs",
        same_people,
        list(family = vapply(zn, `[`, "", 1L), given = vapply(zn, `[`, "", 2L),
             orcid = vapply(zen$creators, function(c) if (is.null(c$orcid)) "" else c$orcid, "")),
        bad = people_bad, good = people_good)
  check(".zenodo.json and CITATION.cff carry the same title",
        function(x) if (identical(x[[1L]], x[[2L]])) character() else sprintf("\"%s\" against \"%s\"", x[[1L]], x[[2L]]),
        list(zen$title, if (is.null(cff)) NA_character_ else cff$title),
        bad = list("A", "B"), good = list("A", "A"))
}

# ---- DESCRIPTION ------------------------------------------------------------
# The Title and Description are read by a person at CRAN, and R CMD check
# enforces little of what that person wants.
description_problems <- function(d) {
  out <- character()
  title <- gsub("\\s+", " ", d$title); text <- gsub("\\s+", " ", d$description)
  if (startsWith(tolower(title), tolower(d$package))) out <- c(out, "the Title starts with the package name")
  if (grepl("\\b(in|for|with) R\\b|\\bR [Pp]ackage\\b", title)) out <- c(out, "the Title says 'in R', 'for R', or 'R package'")
  out <- c(out, sprintf("the Title: %s", title_problems(title)))
  if (grepl(sprintf("^(%s\\b|This package|A package|The package)", d$package), text, ignore.case = TRUE))
    out <- c(out, "the Description starts with the package name or 'This package'")
  if (grepl("doi: |<doi:[^>]* [^>]*>|https?://(dx[.])?doi[.]org", text))
    out <- c(out, "a DOI in the Description is not written <doi:10.prefix/suffix>")
  if (is.na(d$language) || d$language != "en-US") out <- c(out, "Language is not en-US")
  if (is.na(d$url) || is.na(d$bugs)) out <- c(out, "URL or BugReports is missing")
  if (!any(nzchar(d$orcid))) out <- c(out, "no author in Authors@R has an ORCID")
  out
}
check("DESCRIPTION: Title, Description, Language, URL, BugReports, and ORCID as CRAN expects",
      description_problems,
      list(package = pkg, title = field("Title"), description = field("Description"),
           language = field("Language"), url = field("URL"), bugs = field("BugReports"),
           orcid = desc_orcid),
      bad = list(package = "pkg", title = "pkg: tools for R.", description = "This package does things, doi: 10.1/x",
                 language = NA_character_, url = NA_character_, bugs = NA_character_, orcid = ""),
      good = list(package = "pkg", title = "Tools for Study Design", description = "Methods of Kelley (2007) <doi:10.3758/BF03192993>.",
                  language = "en-US", url = "https://example.org", bugs = "https://example.org/issues", orcid = "0000"))

roxygen_pinned <- c(field("Config/roxygen2/version"), field("RoxygenNote"))
roxygen_pinned <- roxygen_pinned[!is.na(roxygen_pinned)][1L]
check(sprintf("the roxygen2 version pinned in DESCRIPTION (%s) is the one installed", roxygen_pinned),
      function(x) if (identical(x[[1L]], x[[2L]])) character() else
        sprintf("DESCRIPTION pins %s; %s is installed, and a different version rewrites every help page", x[[1L]], x[[2L]]),
      list(roxygen_pinned, as.character(utils::packageVersion("roxygen2"))),
      bad = list("7.0.0", "8.1.0"), good = list("8.1.0", "8.1.0"))

site <- if (yaml_ok && file.exists("_pkgdown.yml")) yaml::read_yaml("_pkgdown.yml")$url else NULL
desc_urls <- sub("/$", "", trimws(strsplit(field("URL"), ",")[[1L]]))
check("the site address in _pkgdown.yml is a URL of DESCRIPTION",
      function(u) if (length(u) == 1L && sub("/$", "", u) %in% desc_urls) character() else
        sprintf("_pkgdown.yml has url %s; DESCRIPTION lists %s", if (is.null(u)) "(none)" else u,
                paste(desc_urls, collapse = ", ")),
      site, bad = "https://example.org/nowhere/", good = paste0(desc_urls[1L], "/"))

# ---- Citations in the Text --------------------------------------------------
# A work cited as submitted or in press is swept to its published form
# everywhere at once. The status may follow a parenthesis, a comma, or a
# semicolon, may close a reference in the journal's own form ("Manuscript
# submitted for publication."), or may be a sentence ("is in press at"),
# and it may wrap from one line to the next.
stale <- paste0("[(,;] *(submitted|under review|in preparation|in press|forthcoming|to appear)\\b",
                "|\\bmanuscript (submitted|in preparation|under review)\\b",
                "|\\badvance online publication\\b",
                "|\\b(is|are|was|were) (in press|forthcoming|under review)\\b")
check("no citation carries a stale status (submitted, in press, and so on)",
      function(files) grep_collapsed(stale, files, ignore.case = TRUE),
      lines_of,
      bad = list("control.Rmd" = c("Kelley and Darku (2031, in press) develop the method.",
                                   "Kelley, K. (2031). *A title*. Manuscript submitted for publication.",
                                   "Kelley, K. (2031). A title. *A Journal*. Advance online",
                                   "publication.",
                                   "The article is in press at *A Journal*.")),
      good = list("control.Rmd" = c("The data are submitted to the stopping rule in stages.",
                                    "The interval is in preparation for the next stage, and the press is ignored.",
                                    "Kelley, K. (2031). A title. *A Journal, 1*(2), 3--4.")),
      n_bad = 4L)

# The same for a package: a sentence that says a package is on its way to
# CRAN is false from the day the package is accepted.
pending <- "\\bheaded to CRAN\\b|\\bwill be submitted\\b|\\b(has|have) been submitted to CRAN\\b|\\bnot yet on CRAN\\b"
check("no sentence says a package is still on its way to CRAN",
      function(files) grep_collapsed(pending, files, ignore.case = TRUE),
      lines_of,
      bad = list("control.md" = c("Both packages are headed to CRAN; the second will be", "submitted after the first is accepted.",
                                  "", "The first has been submitted", "to CRAN.")),
      good = list("control.md" = c("Both packages are on CRAN.", "The data are submitted to the stopping rule in stages.")),
      n_bad = 2L)

# Chattopadhyay, Bandyopadhyay, Kelley, and Padalunkal (2025) was once cited
# here under a working title and with wrong initials. The published form is
# the only one: "A sequential approach for noninferiority or equivalence of a
# linear contrast under cost constraints", Psychological Methods, 30(2),
# 425--439. The old form is assembled in pieces so this file does not carry it.
old_2025 <- paste0("budget constrai", "ned sample size plan", "ning,? +under +K +popul", "ations",
                   "|Bandyopadhyay, +", "U[.]|Padalunkal, +", "T[.]|",
                   "given *= *\"U[.]\", *family *= *\"Bandyo|given *= *\"T[.]\", *family *= *\"Padal")
check("the 2025 article is cited in its published form everywhere",
      function(files) grep_collapsed(old_2025, files, ignore.case = TRUE),
      lines_of,
      bad = list("control.R" = c(paste0("#' Chattopadhyay, B., Bandyopadhyay, ", "U., Kelley, K., & Padalunkal,"),
                                 "#' T. (2025). Sequential, budget", paste0("#' constrai", "ned sample size planning under K"),
                                 "#' populations.")),
      good = list("control.R" = c("#' Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal,",
                                  "#' J. J. (2025). A sequential approach for noninferiority or equivalence",
                                  "#' of a linear contrast under cost constraints.")))

# Maxwell, Delaney, and Kelley, fourth edition, is dated 2027 here and in the
# sibling package.
check("Maxwell, Delaney, and Kelley (4th ed.) is dated 2027 everywhere",
      function(files) grep_collapsed("Maxwell,( S[.] E[.],)? Delaney,( H[.] D[.],)? (&|and) Kelley,?( K[.])? ?[(]?20(2[0-689]|[013-9][0-9])",
                                     files),
      lines_of,
      bad = list("control.Rmd" = c("Maxwell, S. E., Delaney, H. D., & Kelley, K.", "(2026). *Designing experiments*")),
      good = list("control.Rmd" = c("Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing", "experiments*",
                                    "Maxwell, Kelley, and Rausch (2008)")))

# ---- House Style -------------------------------------------------------------
em_dash <- intToUtf8(8212L)
check("no em dash in the package text",
      function(files) grep_lines(em_dash, files, fixed = TRUE),
      lines_of[setdiff(names(lines_of), "cran-comments.md")],
      bad = list("control.Rmd" = paste0("The rule ", em_dash, " as stated ", em_dash, " stops early.")),
      good = list("control.Rmd" = "The rule, as stated, stops early (pages 305--324)."))

check("R/ and man/ are ASCII",
      function(files) grep_lines("[^\001-\177]", files, perl = TRUE, useBytes = TRUE),
      lines_of[c(r_files, rd_files)],
      bad = list("control.R" = paste0("# the Gini index of Cram", intToUtf8(233L), "r")),
      good = list("control.R" = "# the Gini index of Cramer"))

# The methodological noun phrases of the house style are written open. A
# cited title keeps its own spelling: a title that hyphenates one of them is
# entered in `cited_titles`, in full and as printed, with the work it
# belongs to. None is cited today.
open_phrases <- paste0("accuracy-in-parameter|in-parameter-estimation|parameter-estimation|minimum-risk|",
                       "risk-point-estimation|point-estimation|effect-sizes?\\b|sample-sizes?\\b|",
                       "confidence-intervals?\\b|chi-squared?\\b|data-generating|factor-analys[ie]s|",
                       "latent-variables?\\b|eta-squared|repeated-measures|model-comparisons?\\b")
cited_titles <- character()
control_title <- "A Sample-Size Table for the Chi-Square Test"
check("the open noun phrases of the house style are never hyphenated (effect size, sample size, and so on)",
      function(x) grep_collapsed(open_phrases, x$files, except = x$titles, ignore.case = TRUE),
      list(files = lines_of, titles = cited_titles),
      bad = list(titles = control_title,
                 files = list("control.Rmd" = c("the accuracy-in-parameter-estimation approach",
                                                "a sample-size plan for a narrow confidence-interval",
                                                "Bell, A. (2031). A Sample-Size Table for the Chi-Squared Test."))),
      good = list(titles = control_title,
                  files = list("control.Rmd" = c("the accuracy in parameter estimation (AIPE) approach",
                                                 "a sample size plan, a two-stage procedure, and a risk-minimizing rule",
                                                 "the mixed-effects model and the effect sizes",
                                                 "Bell, A. (2031). A Sample-Size Table for the",
                                                 "Chi-Square Test. *A Journal, 1*, 2--3."))),
      n_bad = 3L)

# ---- Tests and Vignettes -----------------------------------------------------
# The hygiene detectors were promised to run on CRAN as well as locally, so
# their file exists, has no skip, and the suite is run whole.
check("tests/testthat/test-rd_hygiene.R exists and has no skip, and the suite runs unfiltered",
      function(x) {
        out <- character()
        if (is.null(x$hygiene)) return("tests/testthat/test-rd_hygiene.R is missing")
        src <- read_source(text = x$hygiene, label = "test-rd_hygiene.R")
        rows <- which(src$tk$token == "SYMBOL_FUNCTION_CALL" & grepl("^skip", src$tk$text))
        out <- c(out, at(src, rows, paste0(src$tk$text[rows], "()")))
        if (any(grepl("filter *=", x$runner))) out <- c(out, "tests/testthat.R filters the suite")
        out
      },
      list(hygiene = if (file.exists("tests/testthat/test-rd_hygiene.R"))
             readLines("tests/testthat/test-rd_hygiene.R", warn = FALSE),
           runner = readLines("tests/testthat.R", warn = FALSE)),
      bad = list(hygiene = c("test_that(\"x\", {", "  skip_on_cran()", "  expect_true(TRUE)", "})"),
                 runner = "test_check(\"SMSD\")"),
      good = list(hygiene = c("# no skip here", "test_that(\"x\", expect_true(TRUE))"),
                  runner = "test_check(\"SMSD\")"))

vig_names <- sub("[.]Rmd$", "", basename(vig_installed))
check("every vignette(\"name\") in the package names a vignette that ships",
      function(files) {
        out <- character()
        for (f in names(files)) {
          m <- regmatches(files[[f]], gregexpr("vignette\\(\"[^\"]+\"", files[[f]]))
          for (i in which(lengths(m) > 0L)) {
            nm <- sub("vignette\\(\"([^\"]+)\"", "\\1", m[[i]])
            bad <- setdiff(nm, vig_names)
            if (length(bad)) out <- c(out, sprintf("%s:%d -> %s", f, i, paste(bad, collapse = ", ")))
          }
        }
        out
      },
      lines_of[c(r_files, rd_files, vig_files, "README.md", "NEWS.md")],
      bad = list("control.R" = "#' See vignette(\"no-such-vignette-anywhere\")."),
      good = list("control.R" = sprintf("#' See vignette(\"%s\").", vig_names[1L])))

check("each vignette's title and index entry agree and are AP title case",
      function(files) {
        out <- character()
        for (f in names(files)) {
          l <- files[[f]]
          title <- unquote(trimws(sub("^title: *", "", grep("^title:", l, value = TRUE)[1L])))
          entry <- sub("^.*VignetteIndexEntry\\{(.*)\\} *$", "\\1", grep("VignetteIndexEntry", l, value = TRUE)[1L])
          if (is.na(title) || is.na(entry) || !identical(title, entry))
            out <- c(out, sprintf("%s: the title and the index entry differ", f))
          pr <- title_problems(title)
          if (length(pr)) out <- c(out, sprintf("%s: %s", f, paste(pr, collapse = "; ")))
        }
        out
      },
      lines_of[vig_installed],
      bad = list("control.Rmd" = c("title: \"Sequential AIPE for the standardized mean difference\"",
                                   "  %\\VignetteIndexEntry{Sequential AIPE}")),
      good = list("control.Rmd" = c("title: \"Sequential AIPE for the Standardized Mean Difference\"",
                                    "  %\\VignetteIndexEntry{Sequential AIPE for the Standardized Mean Difference}")))

# ---- Figures Quoted From the Help Pages --------------------------------------
# The vignettes and the README quote measured figures from the help pages:
# coverage, sample sizes, the share of unusable variance estimates. When a
# help page is remeasured, R CMD check says nothing about the quotations it
# leaves behind. The rule that keeps them together:
#
#   * A paragraph that names a help page (as `?name`, as "the help page of
#     `name()`", or as help("name")) has every decimal figure of its prose on
#     one of the pages it names. The name may sit in an HTML comment, which
#     the reader does not see: <!-- figures from `?name` -->.
#   * A sentence that names a page in its prose is held to the pages it
#     names, not to every page of its paragraph, and its whole numbers count
#     too ("about 52 observations" against a page that says 53). The
#     sentences after it in the paragraph are held to the pages it names and
#     to any the paragraph names in a comment; a sentence before the first
#     that names a page is held to all the pages of the paragraph.
#   * A row of a table is held to the pages named in the row or, failing
#     that, in the paragraph before the table. In such a row every number
#     counts, whole numbers too, and the numbers of one row are together on
#     one line of the page (a row of its table) or in one paragraph of it.
#     A row that names no page and whose first cell is a function, as
#     `name()`, describes that function: its decimal figures are on that
#     function's page, as for a paragraph.
#   * A figure that is the paragraph's own, such as a design value of the
#     example beside it, is declared in a comment in that paragraph:
#     <!-- not from a help page: 0.80, 0.30 -->. A declaration in the
#     paragraph before a table covers the table.
#
# Figures in code and in inline R are computed where they stand and are not
# read. A leading zero and a thousands comma are not part of a figure (0.668
# is .668, 1,079 is 1079). A decimal figure in prose may be the page's figure
# rounded to fewer decimals (0.89 for 0.891, 0.67 for .668), so a remeasured
# page whose figure no longer rounds to the quotation is found; in a table
# row every figure is the one the page prints.
rd_aliases <- unlist(lapply(names(rd_db), function(f) {
  al <- vapply(rd_nodes(rd_db[[f]], "\\alias"), function(n) paste(unlist(n), collapse = ""), character(1))
  stats::setNames(rep(f, length(al)), al)
}))
help_text <- local({
  cache <- new.env(parent = emptyenv())
  function(file) {
    if (is.null(cache[[file]])) {
      out <- tempfile(fileext = ".txt")
      on.exit(unlink(out))
      tools::Rd2txt(rd_db[[file]], out = out, options = list(underline_titles = FALSE))
      cache[[file]] <- readLines(out, warn = FALSE, encoding = "UTF-8")
    }
    cache[[file]]
  }
})
figures <- function(text, whole = FALSE) {
  text <- gsub("(?<=[0-9]),(?=[0-9]{3}(?![0-9]))", "", text, perl = TRUE)
  m <- regmatches(text, gregexpr("(?<![A-Za-z0-9_.])[0-9]*[.]?[0-9]+(?![0-9]*[.][0-9]|[A-Za-z0-9_])",
                                 text, perl = TRUE))[[1L]]
  if (!whole) m <- m[grepl(".", m, fixed = TRUE)]
  unique(sub("^0+(?=[.])", "", m, perl = TRUE))
}
# The pages a piece of text names, as the keys of `pages` (a named list of
# rendered pages, each a character vector of lines).
named_pages <- function(text) {
  nm <- c(regmatches(text, gregexpr("`[?][A-Za-z0-9_.]+`", text))[[1L]],
          regmatches(text, gregexpr("help pages? of `[A-Za-z0-9_.]+\\(\\)`", text))[[1L]],
          regmatches(text, gregexpr("help\\(\"[^\"]+\"\\)", text))[[1L]])
  unique(gsub("^`[?]|`$|^help pages? of `|\\(\\)`$|^help\\(\"|\"\\)$", "", nm))
}
# The decimal figures of `figs` that are not on a page whose figures are
# `printed`: neither printed as they are nor the rounding of a figure printed
# with more decimals.
not_quoted <- function(figs, printed) {
  decimals <- function(x) nchar(sub("^[^.]*[.]?", "", x))
  printed_n <- suppressWarnings(as.numeric(printed))
  figs[!vapply(figs, function(fig) {
    if (fig %in% printed) return(TRUE)
    longer <- decimals(printed) > decimals(fig) & !is.na(printed_n)
    any(abs(printed_n[longer] - as.numeric(fig)) <= 0.5 * 10^-decimals(fig) + 1e-9)
  }, logical(1))]
}
# The figures of a sentence. In a sentence that names a page every number
# counts, except a locator (p. 103, Table 1, Equation 8) and a year, which
# cite an article and not the page.
sentence_figures <- function(text, whole) {
  if (whole) {
    text <- gsub("\\b(pp?[.]|Tables?|Equations?|Eqs?[.]|Sections?|Theorems?|Examples?|Chapters?|Figures?|Appendix)\\s*[0-9]+(\\s*(-+|\u2013|to|and)\\s*[0-9]+)?",
                 " ", text, perl = TRUE)
    text <- gsub("(?<![0-9,.])(19|20)[0-9]{2}(?![0-9]|,[0-9]|[.][0-9])", " ", text, perl = TRUE)
  }
  figures(text, whole = whole)
}
# The figures of `figs` that are a percentage in `text` (followed by % or
# "percent") and that a page gives as a proportion instead (95% as 0.95).
as_proportion_on <- function(figs, text, printed) {
  m <- regmatches(text, gregexpr("[0-9][0-9.,]*(?=\\s*(%|percent))", text, perl = TRUE))[[1L]]
  pct <- intersect(figs, figures(paste(m, collapse = " "), whole = TRUE))
  share <- function(x) sub("^0+(?=[.])", "", format(as.numeric(x) / 100, scientific = FALSE, digits = 15), perl = TRUE)
  pct[vapply(pct, function(x) !length(not_quoted(share(x), printed)), logical(1))]
}
declared_own <- function(text) {
  m <- regmatches(text, gregexpr("<!-- *not from a help page:[^>]*-->", text))[[1L]]
  figures(paste(sub("^<!-- *not from a help page:", "", m), collapse = " "), whole = TRUE)
}
# The lines of a paragraph with their comments, code, inline R, and link
# targets blanked. A piece of inline R may run from one line to the next, so
# the paragraph is read as a whole, and each line keeps its place.
prose_only <- function(lines) {
  text <- paste(lines, collapse = "\n")
  for (pattern in c("<!--.*?-->", "`[^`]*`", "(?<=\\])\\([^)]*\\)")) {
    m <- gregexpr(pattern, text, perl = TRUE)
    regmatches(text, m) <- lapply(regmatches(text, m), function(x) gsub("[^\n]", " ", x))
  }
  strsplit(paste0(text, "\n"), "\n", fixed = TRUE)[[1L]]
}
# The sentences of a paragraph's text, as start and end positions. A sentence
# ends at a period, a question mark, or an exclamation mark (with a closing
# parenthesis, quotation mark, or emphasis mark after it) that is followed by
# white space and a capital letter or a piece of code; a decimal point and
# the period of "p. 234" are not. Comments are not read for the ends.
sentence_spans <- function(text) {
  m <- gregexpr("<!--.*?-->", text, perl = TRUE)
  regmatches(text, m) <- lapply(regmatches(text, m), function(x) gsub("[^\n]", " ", x))
  m <- gregexpr("`[^`]*`", text, perl = TRUE)
  regmatches(text, m) <- lapply(regmatches(text, m), function(x) gsub("[^\n]", "X", x))
  e <- gregexpr("[.!?][)\"'*_]*(?=\\s+[A-Z])", text, perl = TRUE)[[1L]]
  ends <- if (e[1L] > 0L) e + attr(e, "match.length") - 1L else integer()
  ends <- c(ends, nchar(text))
  data.frame(start = c(1L, utils::head(ends, -1L) + 1L), end = ends)
}
# The numbers of one table row, all on one line of a page or all in one of
# its paragraphs.
together <- function(tk, page) {
  if (!length(tk)) return(TRUE)
  has_all <- function(text) all(tk %in% figures(text, whole = TRUE))
  if (any(vapply(page, has_all, logical(1)))) return(TRUE)
  blank <- !nzchar(trimws(page))
  paragraphs <- split(page[!blank], cumsum(blank)[!blank])
  # A block of the page that is a table (columns two or more spaces apart on
  # every line) is read a line at a time, above, and not as a whole: the
  # figures of two of its rows are not a quotation of either.
  tabular <- vapply(paragraphs, function(p) all(lengths(gregexpr("\\S {2,}\\S", trimws(p))) >= 2L &
                                                  grepl("\\S {2,}\\S", trimws(p))), logical(1))
  any(vapply(paragraphs[!tabular], function(p) has_all(paste(p, collapse = " ")), logical(1)))
}
quoted_figure_hits <- function(x) {
  pages <- x$pages; aliases <- x$aliases
  resolve <- function(names) unique(unname(aliases[intersect(names, names(aliases))]))
  out <- character()
  for (f in names(x$files)) {
    lines <- x$files[[f]]
    if (length(lines) && lines[1L] == "---") {
      end <- which(lines == "---")
      if (length(end) >= 2L) lines[seq_len(end[2L])] <- ""
    }
    fence <- grep("^\\s*(```|~~~)", lines)
    for (k in seq_len(length(fence) %/% 2L)) lines[fence[2L * k - 1L]:fence[2L * k]] <- ""
    blank <- !nzchar(trimws(lines))
    blocks <- split(which(!blank), cumsum(blank)[!blank])
    before <- list(pages = character(), own = character())
    for (b in blocks) {
      text <- paste(lines[b], collapse = " ")
      named <- named_pages(text)
      unknown <- setdiff(named, names(aliases))
      if (length(unknown))
        out <- c(out, sprintf("%s:%d: no help page is named %s", f, b[1L], paste(unknown, collapse = ", ")))
      here <- resolve(named)
      own <- declared_own(text)
      if (all(grepl("^\\s*[|]", lines[b]))) {
        for (i in b) {
          if (grepl("^\\s*[|][-: |]+$", lines[i])) next
          first <- trimws(strsplit(sub("^\\s*[|]", "", lines[i]), "|", fixed = TRUE)[[1L]][1L])
          row_pages <- resolve(named_pages(lines[i]))
          if (!length(row_pages)) row_pages <- before$pages
          if (!length(row_pages)) {
            described <- if (grepl("^`[A-Za-z0-9_.]+\\(\\)`$", first)) resolve(gsub("^`|\\(\\)`$", "", first))
            if (length(described)) {
              missing <- not_quoted(setdiff(figures(prose_only(lines[i])), before$own),
                                    figures(paste(unlist(pages[described]), collapse = "\n"), whole = TRUE))
              if (length(missing))
                out <- c(out, sprintf("%s:%d: %s is not on %s", f, i, paste(missing, collapse = ", "),
                                      paste(sub("[.]Rd$", "", described), collapse = ", ")))
            }
            next
          }
          tk <- setdiff(figures(prose_only(lines[i]), whole = TRUE), before$own)
          if (!any(vapply(pages[row_pages], function(p) together(tk, p), logical(1))))
            out <- c(out, sprintf("%s:%d: the figures of this row (%s) are not together on %s", f, i,
                                  paste(tk, collapse = ", "), paste(sub("[.]Rd$", "", row_pages), collapse = ", ")))
        }
        # The paragraph before a table speaks for that table and no further.
        before <- list(pages = character(), own = character())
        next
      }
      before <- list(pages = here, own = own)
      if (!length(here)) next
      raw <- paste(lines[b], collapse = "\n")
      prose <- paste(prose_only(lines[b]), collapse = "\n")
      comments <- regmatches(raw, gregexpr("<!--.*?-->", raw, perl = TRUE))[[1L]]
      in_comments <- resolve(named_pages(gsub("\\s+", " ", paste(comments, collapse = " "))))
      line_end <- cumsum(nchar(lines[b]) + 1L) - 1L
      line_start <- line_end - nchar(lines[b]) + 1L
      last_named <- character()
      for (s in seq_len(nrow(spans <- sentence_spans(raw)))) {
        sentence <- substr(raw, spans$start[s], spans$end[s])
        named_here <- resolve(named_pages(gsub("\\s+", " ", gsub("<!--.*?-->", " ", sentence, perl = TRUE))))
        whole <- length(named_here) > 0L
        scope <- if (whole) named_here else if (length(last_named)) union(last_named, in_comments) else here
        if (whole) last_named <- named_here
        on_pages <- figures(paste(unlist(pages[scope]), collapse = "\n"), whole = TRUE)
        for (k in seq_along(b)) {
          from <- max(line_start[k], spans$start[s]); to <- min(line_end[k], spans$end[s])
          if (from > to) next
          text <- substr(prose, from, to)
          missing <- not_quoted(setdiff(sentence_figures(text, whole), own), on_pages)
          missing <- setdiff(missing, as_proportion_on(missing, text, on_pages))
          if (length(missing))
            out <- c(out, sprintf("%s:%d: %s is not on %s", f, b[k], paste(missing, collapse = ", "),
                                  paste(sub("[.]Rd$", "", scope), collapse = ", ")))
        }
      }
    }
  }
  out
}
control_page <- c("Details:", "",
                  "     With normal data coverage was 0.952 at 'width = 0.4' (6,000",
                  "     studies) and .668 with 10 predictors.", "",
                  "       Effect size and true value   Width  Mean n  Coverage ",
                  "       Correlation, 0.4               0.4      55      0.78 ",
                  "       Correlation, 0.4               0.2     265      0.93 ",
                  "       Correlation, 0.9               0.1   1,079      0.45 ")
control_page_z <- c("Details:", "",
                    "     At 'conf_level = 0.95' coverage was 0.31 at 'width = 0.2' and 0.33",
                    "     at 'width = 0.1'; the estimate was unusable in 16.4% of samples.")
control_aliases <- c(seq_x = "seq_x.Rd", seq_x_sensitivity = "seq_x.Rd", seq_z = "seq_z.Rd")
control_pages <- list(seq_x.Rd = control_page, seq_z.Rd = control_page_z)
# The README is read in its source, where a figure its own code produces is
# inline R and a quoted one is prose; README.md alone is read when there is
# no source.
quoted_files <- c(vig_files, if (file.exists("README.Rmd")) "README.Rmd" else "README.md")
check(sprintf("every figure quoted from a help page is on the page it names (%d files)", length(quoted_files)),
      quoted_figure_hits,
      list(files = read_all(quoted_files), aliases = rd_aliases,
           pages = stats::setNames(lapply(unique(unname(rd_aliases)), help_text), unique(unname(rd_aliases)))),
      bad = list(aliases = control_aliases, pages = control_pages,
                 files = list("control.Rmd" = c(
                   "Coverage was 0.951 at a width of 0.4 (`?seq_x`).", "",
                   "The help page of `seq_z()` reports 0.952 at a width of 0.4. The help page of",
                   "`seq_x()` reports 0.668 with 10 predictors.", "",
                   "Coverage was 0.952 at a width of 0.4 in about 5,900 studies (`?seq_x`).", "",
                   "The help page of `seq_x()` reports 0.952. The help page of `seq_z()` reports",
                   "0.31. It also gives 0.668.", "",
                   "The help page of `seq_x_sensitivity()` gives 0.66 with 10 predictors.", "",
                   "The interval was 0.31 wide here; the page has 0.952 (`?seq_y`).", "",
                   "These rows are from `?seq_x`.", "",
                   "| Condition | Width | Mean n | Coverage |", "|---|---|---|---|",
                   "| Correlation, 0.4 | 0.4 | 55 | 0.93 |",
                   "| Correlation, 0.9 | 0.1 | 1,080 | 0.45 |", "",
                   "| Function | Measured |", "|---|---|", "| `seq_x()` | 0.953 at a width of 0.4 in 6,001 studies |"))),
      good = list(aliases = control_aliases, pages = control_pages,
                  files = list("control.Rmd" = c(
                    "---", "title: \"A width of 0.123\"", "---", "",
                    "*Coverage.* The help page of `seq_z()` reports 0.31 at a width of 0.2 for 95%",
                    "intervals (Kelley, 2019, p. 234; Table 3), with 16% of samples unusable. It",
                    "also gives 0.33. The help page of `seq_x()` reports 0.952 in 6,000 studies.", "",
                    "Coverage was 0.33 and 0.952. Both are on the help pages (`?seq_z`, `?seq_x`).", "",
                    "Coverage was 0.952 at a width of 0.4, with 6,000 studies and 0.668",
                    "with 10 predictors (`?seq_x`); the code gives `r 0.5 + 0.25`, `r c(0.61,",
                    "0.62)`, and `0.77`.", "",
                    "```{r}", "0.999  # see ?seq_x", "```", "",
                    "The interval was 0.31 wide here, and the rule is on `?seq_x`.",
                    "<!-- not from a help page: 0.31 -->", "",
                    "The second decimal can differ: 0.93 at 0.2. <!-- figures from `?seq_x` -->", "",
                    "Rounded, the page has 0.95 and 0.67 (`?seq_x`).", "",
                    "A paragraph that names no page may say 0.123.", "",
                    "These rows are from `?seq_x`.", "",
                    "| Condition | Width | Mean n | Coverage |", "|---|---|---|---|",
                    "| Correlation, 0.4 | 0.4 | 55 | 0.78 |",
                    "| Correlation, 0.9 | 0.1 | 1,079 | 0.45 |", "",
                    "| Function | Measured |", "|---|---|",
                    "| `seq_x()` | 0.952 at a width of 0.4, .668 with 12 predictors (Stein, 1945) |",
                    "| `mean()` | 0.123 |"))),
      n_bad = 9L)

# ---- The Repository and the Tarball ------------------------------------------
tracked <- suppressWarnings(system2("git", c("ls-files", "-co", "--exclude-standard"),
                                    stdout = TRUE, stderr = FALSE))
if (!length(tracked)) {
  report(FALSE, "no tracked or trackable file carries a private string",
         "git listed no files here, so nothing was read")
} else {
  check(sprintf("no tracked or trackable file carries a private string (%d files)", length(tracked)),
        function(x) private_in_files(x$files, x$root),
        list(files = tracked, root = "."),
        bad = list(files = "bad.txt", root = control_dir),
        good = list(files = "good.txt", root = control_dir),
        n_bad = length(private_patterns))
}

# A tarball built from the working tree, to see what .Rbuildignore lets
# through. Built without the vignettes, so its size is a floor.
build_dir <- file.path(tempdir(), "release_checks_build")
unlink(build_dir, recursive = TRUE)
dir.create(build_dir)
built <- local({
  old <- setwd(build_dir); on.exit(setwd(old))
  system2("R", c("CMD", "build", "--no-manual", "--no-build-vignettes", shQuote(root)),
          stdout = FALSE, stderr = FALSE) == 0L
})
tb <- list.files(build_dir, pattern = "[.]tar[.]gz$", full.names = TRUE)
if (!built || length(tb) != 1L) {
  report(FALSE, "R CMD build of the working tree succeeds", "run R CMD build by hand for the error")
} else {
  inspect_tarball(tb, "a tarball of the working tree")
}
unlink(c(build_dir, control_dir), recursive = TRUE)

cat(sprintf("\n%d check(s), %d failed.\n", state$run, state$failed))
quit(status = if (state$failed > 0L) 1L else 0L)
