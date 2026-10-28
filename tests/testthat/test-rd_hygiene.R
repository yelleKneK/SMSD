# Documentation and run time hygiene that R CMD check reports only with the
# manuals built, or not at all.
#
# Each pattern below has caused a return from CRAN's pretest or from a CRAN
# reviewer for a package whose check was clean, or is CRAN policy:
#
#   help pages    Rd markup inside \eqn{} or \deqn{} (the PDF manual stops on
#                 it); \donttest{} and \dontrun{} (one \donttest{} makes
#                 R CMD check --as-cran run every example twice); an example
#                 gated on requireNamespace(); comment lines in an example
#                 that are R code.
#   functions     reading or writing the random number generator's state or
#                 the global environment; the superassignment operator;
#                 setting the warn option; a file path supplied by default or
#                 written into a call; changing par(), options(), or the
#                 working directory without restoring them on exit;
#                 installing packages, quitting, or starting other software;
#                 taking more than two cores, by default or by asking the
#                 machine how many it has; and writing to the console from a
#                 function body, which is what an earlier submission of this
#                 package was returned for.
#
# The help pages are parsed and every function in the namespace is deparsed
# and tokenized, so nothing here depends on how a line happens to be wrapped
# or on a comment that mentions a pattern. Every detector is first shown a
# known bad input, which it must flag, and a clean twin, which it must pass:
# a detector that cannot fire never passes by default. The file runs on CRAN
# as well as locally and has no skip; if the help pages or the namespace
# cannot be read, that is a failure and not a pass.
#
# The known bad inputs are built from text and are never run.
#
# The rule for comment lines in an example is the reviewer's, applied
# mechanically: no comment line may parse as R code. Most prose does not
# parse, but a short line can. A word followed by a parenthesis is a call
# ("Kelley (2007)", "Equation (15)"), and a hyphenated pair or bare
# arithmetic is an expression ("two-stage", "N - 1", "sigma^2"). When the
# test names such a line, reword it as a sentence: "Kelley, 2007,
# Equation 15." does not parse.
#
# tools/release_checks.R evaluates this file down to the line that begins
# "# ---- Reading the Package" and applies the same detectors to the source
# files, where it can name the file and the line. Nothing above that line
# may read the package or run a test.

# ---- The Package Under Test and Its Help Pages --------------------------

package_under_test <- function() {
  pkg <- tryCatch(testthat::testing_package(), error = function(e) "")
  if (!nzchar(pkg)) {
    desc <- file.path("..", "..", "DESCRIPTION")
    if (file.exists(desc)) pkg <- unname(read.dcf(desc, "Package")[1L])
  }
  pkg
}

# The Rd sources when the tests run from the source tree, and the installed
# help database when they run from an installed package (as on CRAN).
rd_sources <- function() {
  man <- file.path("..", "..", "man")
  tryCatch(
    if (dir.exists(man)) tools::Rd_db(dir = file.path("..", ".."))
    else tools::Rd_db(package_under_test()),
    error = function(e) list()
  )
}

rd_from_text <- function(...) {
  con <- textConnection(paste0(...))
  on.exit(close(con))
  tools::parse_Rd(con)
}

# Every node of a parsed help page whose tag is one of `tags`.
rd_nodes <- function(x, tags) {
  tag <- attr(x, "Rd_tag")
  hit <- if (!is.null(tag) && tag %in% tags) list(x) else list()
  if (is.list(x)) for (el in x) hit <- c(hit, rd_nodes(el, tags))
  hit
}

math_text <- function(rd) {
  vapply(rd_nodes(rd, c("\\eqn", "\\deqn")),
         function(node) paste(unlist(node), collapse = ""), character(1))
}

examples_text <- function(rd) {
  paste(unlist(rd_nodes(rd, "\\examples")), collapse = "")
}

sweep_pages <- function(db, detector) {
  hits <- lapply(db, detector)
  hits <- hits[lengths(hits) > 0L]
  as.character(unlist(lapply(names(hits), function(page)
    paste(page, hits[[page]], sep = ": "))))
}

# ---- Detectors on Help Pages --------------------------------------------

# Rd markup inside math. \code{} inside \eqn{} becomes \texttt in LaTeX math
# and stops the PDF manual; the HTML manual reports a math rendering problem.
# A macro may carry an option before its brace, as in \link[stats]{qnorm}.
markup_in_math <- paste0(
  "\\\\(code|emph|link|linkS4class|strong|pkg|bold|verb|var|samp|file|env|",
  "dQuote|sQuote|href|url|doi|acronym|cite|dfn|kbd|option|command|email|",
  "CRANpkg|abbr|enc|out|Sexpr)[[:space:]]*(\\[[^]]*\\])?[[:space:]]*\\{")

math_markup_hits <- function(rd) {
  math <- math_text(rd)
  math[grepl(markup_in_math, math)]
}

wrapper_hits <- function(rd) {
  vapply(rd_nodes(rd, c("\\donttest", "\\dontrun")),
         function(node) attr(node, "Rd_tag"), character(1))
}

# A comment that mentions requireNamespace() is not a guard.
guard_hits <- function(rd) {
  ex <- strsplit(examples_text(rd), "\n", fixed = TRUE)[[1L]]
  ex <- ex[!grepl("^\\s*#", ex)]
  trimws(ex[grepl("requireNamespace", ex, fixed = TRUE)])
}

# Whether text parses to a call or an assignment. A formula alone does not
# count: "# y ~ x" in an example names a model and runs nothing.
is_code <- function(txt) {
  parsed <- tryCatch(parse(text = txt, keep.source = FALSE),
                     error = function(e) NULL)
  if (is.null(parsed) || !length(parsed)) return(FALSE)
  is_call <- vapply(parsed, is.call, logical(1))
  is_formula <- vapply(parsed, function(e)
    is.call(e) && identical(e[[1L]], as.name("~")), logical(1))
  any(is_call & !is_formula)
}

# Comment lines in an example that are R code. A CRAN reviewer's words: "Some
# code lines in examples are commented out. Please never do that." Every
# window of consecutive comment lines is parsed with the comment markers
# removed; a window that parses to a call or an assignment is flagged,
# whatever its intent. That includes a line of prose that happens to be
# valid R (see the top of this file), which is reworded.
commented_code <- function(rd) {
  lines <- strsplit(examples_text(rd), "\n", fixed = TRUE)[[1L]]
  cm <- grep("^\\s*#", lines)
  if (!length(cm)) return(character())
  hits <- integer()
  for (run in split(cm, cumsum(c(1L, diff(cm) != 1L)))) {
    txt <- sub("^\\s*#+\\s?", "", lines[run])
    n <- length(run)
    for (a in seq_len(n)) for (b in a:n) {
      window <- paste(txt[a:b], collapse = "\n")
      if (nzchar(trimws(window)) && is_code(window)) hits <- c(hits, run[a:b])
    }
  }
  trimws(lines[sort(unique(hits))])
}

# ---- Reading Functions --------------------------------------------------

fn_from_text <- function(text) {
  eval(parse(text = text, keep.source = FALSE), envir = baseenv())
}

# The terminal tokens of a function, in order. Source references are dropped
# first, so a comment never reaches a detector and the result is the same
# for an installed package and for one loaded from source. deparse() writes
# a number as the session's scipen option says (pos = 1 becomes
# pos = 1e+00 when scipen is negative), and the detectors read numbers, so
# the default is used here whatever the session has. NULL when the function
# cannot be read this way.
function_tokens <- function(f) {
  withr::local_options(list(scipen = 0L))
  txt <- deparse(utils::removeSource(f), width.cutoff = 500L)
  parsed <- tryCatch(parse(text = txt, keep.source = TRUE),
                     error = function(e) NULL)
  if (is.null(parsed)) return(NULL)
  d <- utils::getParseData(parsed)
  if (is.null(d)) return(NULL)
  d <- d[d$terminal, c("line1", "col1", "token", "text")]
  d <- d[order(d$line1, d$col1), ]
  rownames(d) <- NULL
  d
}

# What a detector on functions returns: one row for each finding, with the
# row of the token table it was found at and what was found.
hit <- function(rows = integer(), what = character()) {
  rows <- as.integer(rows)
  data.frame(row = rows, what = rep_len(as.character(what), length(rows)),
             stringsAsFactors = FALSE)
}

unquote <- function(x) gsub("^[\"']|[\"']$", "", x)

is_call_to <- function(tk, names) {
  tk$token == "SYMBOL_FUNCTION_CALL" & tk$text %in% names
}

token_before <- function(tk) c("", tk$token)[seq_len(nrow(tk))]
token_after <- function(tk) c(tk$token, "")[seq_len(nrow(tk)) + 1L]
is_qualified <- function(tk) token_before(tk) %in% c("NS_GET", "NS_GET_INT")

# The arguments of the call whose name is token `i`: a list with one vector
# of token rows for each argument.
call_args <- function(tk, i) {
  n <- nrow(tk)
  if (i >= n || tk$token[i + 1L] != "'('") return(list())
  opening <- c("'('" = 1L, "'['" = 1L, "LBB" = 2L, "'{'" = 1L)
  closing <- c("')'" = 1L, "']'" = 1L, "'}'" = 1L)
  depth <- 1L
  args <- list()
  current <- integer()
  j <- i + 2L
  while (j <= n) {
    tok <- tk$token[j]
    if (tok %in% names(opening)) depth <- depth + opening[[tok]]
    if (tok %in% names(closing)) depth <- depth - closing[[tok]]
    if (depth == 0L) break
    if (depth == 1L && tok == "','") {
      args[[length(args) + 1L]] <- current
      current <- integer()
    } else {
      current <- c(current, j)
    }
    j <- j + 1L
  }
  if (length(current) || length(args)) args[[length(args) + 1L]] <- current
  args
}

is_named_arg <- function(tk, arg) {
  length(arg) >= 2L && tk$token[arg[2L]] == "EQ_SUB"
}

is_string_literal <- function(tk, rows) {
  length(rows) == 1L && tk$token[rows] == "STR_CONST" &&
    nzchar(unquote(tk$text[rows]))
}

# The tokens of a formal's default, which begins at row `j`: up to the comma
# or the parenthesis that ends it.
default_tokens <- function(tk, j) {
  n <- nrow(tk)
  depth <- 0L
  rows <- integer()
  while (j <= n) {
    tok <- tk$token[j]
    if (tok %in% c("'('", "'['", "'{'")) {
      depth <- depth + 1L
    } else if (tok == "LBB") {
      depth <- depth + 2L
    } else if (tok %in% c("')'", "']'", "'}'")) {
      if (depth == 0L) break
      depth <- depth - 1L
    } else if (tok == "','" && depth == 0L) {
      break
    }
    rows <- c(rows, j)
    j <- j + 1L
  }
  rows
}

# The names a function binds itself: its formals, its local variables, and
# its loop variables. A call to q() in a function that has defined a local
# q is a call to that local and not to base::q().
local_names <- function(tk) {
  assigned <- tk$token == "SYMBOL" &
    !token_before(tk) %in% c("'$'", "'@'") &
    (token_after(tk) %in% c("LEFT_ASSIGN", "EQ_ASSIGN", "IN") |
       token_before(tk) == "RIGHT_ASSIGN")
  unique(tk$text[tk$token == "SYMBOL_FORMALS" | assigned])
}

# The rows where one of the functions in `names` is used without being
# called on the spot: passed to another function (lapply(x, print)), stored
# (show <- message), or named in a string that is handed to a function that
# calls its argument (do.call("cat", x)). A name the function has bound
# itself is that local object and is left alone, as is a list element
# (e$message).
functionals <- c("do.call", "match.fun", "lapply", "sapply", "vapply", "Map",
                 "mapply", "Reduce", "Filter", "apply", "tapply", "outer",
                 "get", "get0", "mget", "getExportedValue",
                 "getFromNamespace")

value_uses <- function(tk, names) {
  symbols <- which(tk$token == "SYMBOL" & tk$text %in% names &
                     !token_before(tk) %in% c("'$'", "'@'") &
                     (!tk$text %in% local_names(tk) | is_qualified(tk)))
  strings <- integer()
  for (i in which(is_call_to(tk, functionals))) {
    rows <- unlist(call_args(tk, i))
    strings <- c(strings, rows[tk$token[rows] == "STR_CONST" &
                                 unquote(tk$text[rows]) %in% names])
  }
  sort(unique(c(symbols, strings)))
}

# ---- Detectors on Functions ---------------------------------------------

# CRAN policy: a package does not touch the global environment, and that
# includes the generator's state, which lives there. Seeds go through one
# helper built on the withr package, which sets the seed and restores the
# caller's stream; a function that sets the seed itself leaves the caller's
# stream changed. The global environment can also be reached by position:
# as.environment(1), sys.frame(0), and pos = 1 in assign() and its
# relatives.
env_accessors <- c("assign", "get", "get0", "mget", "exists", "rm", "remove",
                   "eval", "evalq", "delayedAssign", "makeActiveBinding")

global_state_hits <- function(tk) {
  watched <- c(".Random.seed", ".GlobalEnv")
  fns <- c("globalenv", "set.seed", "RNGkind")
  named <- which((tk$token == "SYMBOL" & tk$text %in% watched) |
                   (tk$token == "STR_CONST" & unquote(tk$text) %in% watched))
  called <- which(is_call_to(tk, fns))
  passed <- value_uses(tk, fns)
  out <- rbind(
    hit(named, unquote(tk$text[named])),
    hit(called, paste0(tk$text[called], "()")),
    hit(passed, paste(unquote(tk$text[passed]), "passed as a function")))
  is_number <- function(rows, values) {
    length(rows) == 1L && tk$token[rows] == "NUM_CONST" &&
      tk$text[rows] %in% values
  }
  for (i in which(is_call_to(tk, c("as.environment", "sys.frame")))) {
    args <- call_args(tk, i)
    global <- if (tk$text[i] == "sys.frame") c("0", "0L") else c("1", "1L")
    if (length(args) && is_number(args[[1L]], global))
      out <- rbind(out, hit(i, sprintf("%s(%s)", tk$text[i],
                                       tk$text[args[[1L]]])))
  }
  for (i in which(is_call_to(tk, env_accessors))) {
    for (a in call_args(tk, i)) {
      if (is_named_arg(tk, a) && tk$text[a[1L]] %in% c("pos", "envir") &&
          is_number(a[-(1:2)], c("1", "1L")))
        out <- rbind(out, hit(i, sprintf("%s(%s = 1)", tk$text[i],
                                         tk$text[a[1L]])))
    }
  }
  out
}

# Legal when it modifies a local of an enclosing function, but a reviewer
# reads the operator as a write to the global environment.
# assign(inherits = TRUE) is the same write by another name.
superassign_hits <- function(tk) {
  ops <- which(tk$token %in% c("LEFT_ASSIGN", "RIGHT_ASSIGN") &
                 tk$text %in% c("<<-", "->>"))
  out <- hit(ops, tk$text[ops])
  for (i in which(is_call_to(tk, "assign"))) {
    for (a in call_args(tk, i)) {
      if (is_named_arg(tk, a) && tk$text[a[1L]] == "inherits" &&
          length(a) == 3L && tk$text[a[3L]] %in% c("TRUE", "T"))
        out <- rbind(out, hit(i, "assign(inherits = TRUE)"))
    }
  }
  out
}

# CRAN rejects a negative warn outright, restored on exit or not. A known
# warning is handled at the call that raises it, so no function sets warn:
# not in options() or its withr relatives, and not in a list or a vector
# built to be handed to them.
warn_option_hits <- function(tk) {
  out <- hit()
  for (i in which(is_call_to(tk, c("options", "local_options",
                                   "with_options")))) {
    rows <- unlist(call_args(tk, i))
    named <- rows[tk$token[rows] %in% c("SYMBOL_SUB", "STR_CONST") &
                    unquote(tk$text[rows]) == "warn"]
    named <- named[named < nrow(tk)]
    if (any(tk$token[named + 1L] == "EQ_SUB"))
      out <- rbind(out, hit(i, paste0(tk$text[i], "(warn = )")))
  }
  for (i in which(is_call_to(tk, c("list", "c")))) {
    for (a in call_args(tk, i)) {
      value <- a[-(1:2)]
      if (is_named_arg(tk, a) && unquote(tk$text[a[1L]]) == "warn" &&
          length(value) &&
          all(tk$token[value] %in% c("NUM_CONST", "'-'", "'+'")) &&
          any(grepl("^[0-9.]", tk$text[value])))
        out <- rbind(out, hit(i, sprintf("%s(warn = %s)", tk$text[i],
                                         paste(tk$text[value], collapse = ""))))
    }
  }
  out
}

# CRAN policy: nothing is written outside the session's temporary directory
# by default. A function that can write takes `filename = NULL` and writes
# only when the user supplies a path. An argument named for a path has a
# default path when the default is a string or a call that builds a path;
# a flag (path = TRUE), a choice (dir = c("upper", "lower")), a connection
# (con = stdout()), and the empty string that means the console are not
# paths. A single string is read as one whatever it says, because that is
# how a reviewer reads it: name the argument for what it is.
path_args <- c("file", "filename", "file_name", "filepath", "file_path",
               "path", "dir", "directory", "destfile", "outfile",
               "output_file", "con")
path_builders <- c("file.path", "tempfile", "tempdir", "getwd", "path.expand",
                   "normalizePath", "paste", "paste0", "sprintf", "Sys.getenv",
                   "R.home", "system.file", "basename", "dirname", "here")
path_first_writers <- c("sink", "pdf", "png", "jpeg", "tiff", "bmp", "svg",
                        "postscript", "file.create", "dir.create")
path_second_writers <- c("write", "write.csv", "write.csv2", "write.table",
                         "saveRDS", "dump", "writeLines")

default_path_hits <- function(tk) {
  n <- nrow(tk)
  out <- hit()
  # The default of a formal, of the function or of one defined inside it.
  for (i in which(tk$token == "SYMBOL_FORMALS" & token_after(tk) == "EQ_FORMALS")) {
    value <- default_tokens(tk, i + 2L)
    if (tk$text[i] %in% path_args) {
      head <- value[!tk$token[value] %in% c("SYMBOL_PACKAGE", "NS_GET",
                                            "NS_GET_INT")][1L]
      builds <- !is.na(head) && tk$token[head] == "SYMBOL_FUNCTION_CALL" &&
        tk$text[head] %in% path_builders
      if (is_string_literal(tk, value) || builds)
        out <- rbind(out, hit(i, sprintf("argument %s has the default %s",
                                         tk$text[i],
                                         paste(tk$text[value], collapse = ""))))
    }
    if (tk$text[i] == "save" && length(value) == 1L &&
        tk$text[value] %in% c("TRUE", "FALSE", "T", "F"))
      out <- rbind(out, hit(i, "a save switch"))
  }
  # A quoted path as a named argument of any call.
  named <- which(tk$token == "SYMBOL_SUB" & tk$text %in% path_args)
  for (i in named[named + 2L <= n]) {
    if (tk$token[i + 1L] == "EQ_SUB" && is_string_literal(tk, i + 2L) &&
        (i + 3L > n || tk$token[i + 3L] %in% c("','", "')'")))
      out <- rbind(out, hit(i, sprintf("%s = %s", tk$text[i],
                                       tk$text[i + 2L])))
  }
  # A quoted path in the position where a writer takes its destination.
  for (i in which(is_call_to(tk, c(path_first_writers, path_second_writers)))) {
    args <- call_args(tk, i)
    unnamed <- args[!vapply(args, function(a) is_named_arg(tk, a), logical(1))]
    k <- if (tk$text[i] %in% path_first_writers) 1L else 2L
    if (length(unnamed) >= k && is_string_literal(tk, unnamed[[k]]))
      out <- rbind(out, hit(i, sprintf("%s(%s)", tk$text[i],
                                       tk$text[unnamed[[k]]])))
  }
  out
}

# CRAN policy on the user's session: par(), options(), and the working
# directory are restored by the function that changes them, and a package
# does not install packages, quit, or start other software. A change counts
# as restored only when the same function is called again inside an
# on.exit() (or a withr::defer()): an on.exit() that closes a connection
# restores nothing.
session_state_hits <- function(tk) {
  banned <- c("install.packages", "q", "quit", "system", "system2", "shell",
              "shell.exec", "browseURL")
  rows <- which(is_call_to(tk, banned) &
                  (!tk$text %in% local_names(tk) | is_qualified(tk)))
  out <- hit(rows, paste0(tk$text[rows], "()"))
  restored <- unlist(lapply(which(is_call_to(tk, c("on.exit", "defer"))),
                            function(i) {
                              inner <- unlist(call_args(tk, i))
                              tk$text[inner][tk$token[inner] == "SYMBOL_FUNCTION_CALL"]
                            }))
  only_reads <- function(args) {
    all(vapply(args, function(a) {
      is_string_literal(tk, a) ||
        (is_named_arg(tk, a) && tk$text[a[1L]] == "no.readonly") ||
        (length(a) >= 4L && tk$token[a[1L]] == "SYMBOL_FUNCTION_CALL" &&
           tk$text[a[1L]] == "c" &&
           all(tk$token[a[-1L]] %in% c("'('", "')'", "','", "STR_CONST")))
    }, logical(1)))
  }
  for (i in which(is_call_to(tk, c("par", "options", "setwd")))) {
    if (tk$text[i] != "setwd" && only_reads(call_args(tk, i))) next
    if (!tk$text[i] %in% restored)
      out <- rbind(out, hit(i, sprintf(
        "%s() changed and not restored by on.exit() in the function",
        tk$text[i])))
  }
  out
}

# This package's own return from CRAN: the sequential functions once wrote
# their progress to the console and expected to be called again. A function
# returns its result, and only a print method writes. A writer handed to
# another function to call (lapply(x, print)) writes just the same; a
# writer whose output is captured as text does not.
console_writers <- c("cat", "print", "print.default", "message", "writeLines",
                     "str", "txtProgressBar", "setTxtProgressBar")

console_hits <- function(tk) {
  captured <- unlist(lapply(which(is_call_to(tk, "capture.output")),
                            function(i) unlist(call_args(tk, i))))
  called <- setdiff(which(is_call_to(tk, console_writers)), captured)
  passed <- setdiff(value_uses(tk, console_writers), captured)
  rbind(hit(called, paste0(tk$text[called], "()")),
        hit(passed, paste(unquote(tk$text[passed]), "passed as a function")))
}

# CRAN policy: a package uses at most two cores unless the user asks for
# more, in its examples, tests, and vignettes as in its functions. No
# function asks the machine how many cores it has (detectCores() and its
# relatives: the answer on a check server is not two), has more than two as
# the default of an argument that counts cores, workers, or threads, or
# writes more than two into a call. A count the user supplies passes.
core_counters <- c("detectCores", "availableCores", "availableWorkers",
                   "defaultNumThreads")
core_args <- c("mc.cores", "cores", "ncores", "n_cores", "num_cores",
               "ncpus", "Ncpus", "n_cpus", "workers", "n_workers", "nworkers",
               "threads", "nthreads", "n_threads", "nthread", "num_threads",
               "OMP_NUM_THREADS", "MC_CORES")
worker_starters <- c("makeCluster", "makePSOCKcluster", "makeForkCluster",
                     "registerDoParallel", "registerDoMC", "setDTthreads",
                     "setThreadOptions")

core_count_hits <- function(tk) {
  n <- nrow(tk)
  # A number above two, however it is written (4, 4L, 4e+00).
  over_two <- function(rows) {
    length(rows) == 1L && tk$token[rows] == "NUM_CONST" &&
      isTRUE(suppressWarnings(as.numeric(sub("L$", "", tk$text[rows]))) > 2)
  }
  called <- which(is_call_to(tk, core_counters))
  passed <- value_uses(tk, core_counters)
  out <- rbind(
    hit(called, paste0(tk$text[called], "()")),
    hit(passed, paste(unquote(tk$text[passed]), "passed as a function")))
  for (i in which(tk$token == "SYMBOL_FORMALS" & tk$text %in% core_args &
                  token_after(tk) == "EQ_FORMALS")) {
    value <- default_tokens(tk, i + 2L)
    if (over_two(value))
      out <- rbind(out, hit(i, sprintf("argument %s has the default %s",
                                       tk$text[i], tk$text[value])))
  }
  named <- which(tk$token %in% c("SYMBOL_SUB", "STR_CONST") &
                   unquote(tk$text) %in% core_args)
  for (i in named[named + 2L <= n]) {
    if (tk$token[i + 1L] == "EQ_SUB" && over_two(i + 2L) &&
        (i + 3L > n || tk$token[i + 3L] %in% c("','", "')'")))
      out <- rbind(out, hit(i, sprintf("%s = %s", unquote(tk$text[i]),
                                       tk$text[i + 2L])))
  }
  for (i in which(is_call_to(tk, worker_starters))) {
    args <- call_args(tk, i)
    unnamed <- args[!vapply(args, function(a) is_named_arg(tk, a), logical(1))]
    if (length(unnamed) && over_two(unnamed[[1L]]))
      out <- rbind(out, hit(i, sprintf("%s(%s)", tk$text[i],
                                       tk$text[unnamed[[1L]]])))
  }
  out
}

# Only a print method that the package registers may write: a function
# that is merely named print_something is not one, and a format method
# returns its text.
may_write <- function(name, print_methods) name %in% print_methods

flags <- function(detector, text) {
  nrow(detector(function_tokens(fn_from_text(text)))) > 0L
}

sweep_namespace <- function(tokens, detector) {
  hits <- lapply(tokens, function(tk) unique(detector(tk)$what))
  hits <- hits[lengths(hits) > 0L]
  as.character(unlist(lapply(names(hits), function(nm)
    paste(nm, hits[[nm]], sep = ": "))))
}

namespace_functions <- function() {
  ns <- asNamespace(package_under_test())
  objs <- mget(ls(ns, all.names = TRUE), envir = ns, inherits = FALSE)
  fns <- Filter(is.function, objs)
  fns[!vapply(fns, is.primitive, logical(1))]
}

# The print methods the package registers, by the name of the function: the
# S3method(print, class) directives of its NAMESPACE file, which is beside
# DESCRIPTION in the source tree and in the installed package alike.
registered_print_methods <- function() {
  file <- system.file("NAMESPACE", package = package_under_test())
  directives <- tryCatch(parse(file, keep.source = FALSE),
                         error = function(e) NULL)
  word <- function(x) unquote(paste(deparse(x), collapse = ""))
  methods <- character()
  for (d in as.list(directives)) {
    if (!is.call(d) || !identical(d[[1L]], as.name("S3method")) ||
        length(d) < 3L || word(d[[2L]]) != "print") next
    methods <- c(methods, if (length(d) >= 4L) word(d[[4L]])
                 else paste("print", word(d[[3L]]), sep = "."))
  }
  methods
}

# ---- Reading the Package ------------------------------------------------

rd_db <- rd_sources()
ns_fns <- namespace_functions()
ns_tokens <- lapply(ns_fns, function_tokens)
print_methods <- registered_print_methods()

# ---- The Tests ----------------------------------------------------------

test_that("the help pages and every function in the namespace can be read", {
  expect_true(nzchar(package_under_test()))
  expect_gt(length(rd_db), 0L)
  expect_gt(length(ns_fns), 0L)
  unread <- names(ns_tokens)[vapply(ns_tokens, is.null, logical(1))]
  expect_identical(unread, character(0))
  # The registered print methods are found, and each is a function of the
  # namespace; were they not, no function would be allowed to write.
  expect_true("print.smsd_check" %in% print_methods)
  expect_identical(setdiff(print_methods, names(ns_fns)), character(0))
})

test_that("no help page carries Rd markup inside \\eqn{} or \\deqn{}", {
  bad <- rd_from_text(
    "\\name{x}\\alias{x}\\title{x}\\description{The region is ",
    "\\eqn{(-\\code{delta_lower}, \\code{delta_upper})}.}")
  split_lines <- rd_from_text(
    "\\name{x}\\alias{x}\\title{x}\\description{The bound is \\deqn{n \\ge\n",
    "\\emph{z}^2 / d^2}{n >= z^2 / d^2}.}")
  with_option <- rd_from_text(
    "\\name{x}\\alias{x}\\title{x}\\description{The quantile is ",
    "\\eqn{\\link[stats]{qnorm}(0.975)}, from \\eqn{\\CRANpkg{stats}}.}")
  good <- rd_from_text(
    "\\name{x}\\alias{x}\\title{x}\\description{The region is ",
    "\\eqn{(-\\delta_L, +\\delta_U)}, with \\code{delta_lower} as ",
    "\\eqn{\\delta_L} and \\eqn{\\text{conf\\_level}} as given; the ",
    "quantile \\eqn{z} is from \\code{\\link[stats]{qnorm}}, and ",
    "\\eqn{\\varepsilon + \\overline{x} + \\mathrm{Var}(x)} is plain LaTeX.}")
  expect_length(math_markup_hits(bad), 1L)
  expect_length(math_markup_hits(split_lines), 1L)
  expect_length(math_markup_hits(with_option), 2L)
  expect_length(math_markup_hits(good), 0L)

  expect_identical(sweep_pages(rd_db, math_markup_hits), character(0))
})

test_that("no help page uses \\donttest{} or \\dontrun{}", {
  bad <- rd_from_text(
    "\\name{x}\\alias{x}\\title{x}\\examples{\nmean(1:3)\n",
    "\\donttest{\nmean(1:1e6)\n}\n\\dontrun{\nmean(1:1e9)\n}\n}")
  good <- rd_from_text(
    "\\name{x}\\alias{x}\\title{x}\\examples{\nmean(1:3)\n}")
  expect_setequal(wrapper_hits(bad), c("\\donttest", "\\dontrun"))
  expect_length(wrapper_hits(good), 0L)

  expect_identical(sweep_pages(rd_db, wrapper_hits), character(0))
})

test_that("no example is gated on requireNamespace()", {
  bad <- rd_from_text(
    "\\name{x}\\alias{x}\\title{x}\\examples{\n",
    "if (requireNamespace(\"stats\", quietly = TRUE)) {\n",
    "  stats::median(1:3)\n}\n}")
  good <- rd_from_text(
    "\\name{x}\\alias{x}\\title{x}\\examples{\n",
    "# stats is attached, so no requireNamespace() is needed here.\n",
    "stats::median(1:3)\n}")
  expect_length(guard_hits(bad), 1L)
  expect_length(guard_hits(good), 0L)

  expect_identical(sweep_pages(rd_db, guard_hits), character(0))
})

test_that("no example carries commented-out code", {
  bad <- rd_from_text(
    "\\name{x}\\alias{x}\\title{x}\\examples{\n",
    "# The fit is shown for comparison.\n",
    "# fit_model(X, Y, M, outcome = \"continuous\")\n",
    "# start <- Sys.time()\n",
    "# Sys.sleep(0.2)\n",
    "default_grid(300, 500)\n",
    "# A call can be commented out across lines:\n",
    "# fit_model(X, Y,\n",
    "#           outcome = \"continuous\")\n}")
  hits <- commented_code(bad)
  expect_length(hits, 5L)
  expect_false(any(grepl("shown for comparison|across lines", hits)))
  good <- rd_from_text(
    "\\name{x}\\alias{x}\\title{x}\\examples{\n",
    "# Twenty observations per group, then the stopping rule.\n",
    "# The model is y ~ x, as in Kelley, 2007, Equation 15.\n",
    "# y ~ x\n",
    "default_grid(300, 500)  # the pilot sample size\n}")
  expect_length(commented_code(good), 0L)
  # The rule is mechanical: prose that happens to be valid R is flagged
  # too, and is reworded (see the top of this file).
  prose_that_parses <- rd_from_text(
    "\\name{x}\\alias{x}\\title{x}\\examples{\n",
    "# Kelley (2007)\n",
    "default_grid(300, 500)\n",
    "# N - 1\n}")
  expect_identical(commented_code(prose_that_parses),
                   c("# Kelley (2007)", "# N - 1"))

  expect_identical(sweep_pages(rd_db, commented_code), character(0))
})

test_that("no function touches the generator's state or the global environment", {
  detector <- global_state_hits
  expect_true(flags(detector, paste0(
    "function(seed) { old <- get(\".Random.seed\", envir = globalenv()); ",
    "on.exit(assign(\".Random.seed\", old, envir = globalenv())); ",
    "stats::rnorm(1) }")))
  expect_true(flags(detector, "function() exists(\".Random.seed\")"))
  expect_true(flags(detector, "function(x) assign(\"x\", x, envir = .GlobalEnv)"))
  expect_true(flags(detector, "function(x) rm(x, envir = globalenv())"))
  expect_true(flags(detector, "function(seed) { set.seed(seed); stats::rnorm(1) }"))
  # The same things reached indirectly.
  expect_true(flags(detector,
                    "function(seed) { do.call(set.seed, list(seed)); stats::rnorm(1) }"))
  expect_true(flags(detector,
                    "function(seed) { do.call(\"set.seed\", list(seed)); stats::rnorm(1) }"))
  expect_true(flags(detector, "function(x) assign(\"x\", x, pos = 1)"))
  expect_true(flags(detector,
                    "function(x) assign(\"x\", x, envir = as.environment(1L))"))
  expect_true(flags(detector, "function(x) assign(\"x\", x, envir = sys.frame(0))"))
  expect_false(flags(detector, paste0(
    "function(seed) withr::with_seed(seed, stats::rnorm(1))")))
  expect_false(flags(detector, paste0(
    "function(x) { e <- new.env(parent = emptyenv()); ",
    "assign(\"seed\", x, envir = e); e }")))
  # pos = 1 is "below" for a label and the first frame up for a caller.
  expect_false(flags(detector, paste0(
    "function(x, y) { graphics::text(x, y, labels = \"a\", pos = 1); ",
    "get(\"x\", envir = parent.frame(1)); ",
    "stop(\"call set.seed() before the simulation\") }")))

  expect_identical(sweep_namespace(ns_tokens, detector), character(0))
})

test_that("functions are read the same way under any scipen option", {
  # A session that prefers scientific notation must not hide pos = 1, nor
  # turn a count of cores into another number.
  withr::local_options(list(scipen = -5L, digits = 3L))
  expect_true(flags(global_state_hits, "function(x) assign(\"x\", x, pos = 1)"))
  expect_true(flags(core_count_hits, "function(x, cores = 4L) x"))
  expect_false(flags(global_state_hits, paste0(
    "function(x, y) graphics::text(x, y, labels = \"a\", pos = 1)")))
  tk <- function_tokens(fn_from_text("function(x) x + 100000"))
  expect_identical(tk$text[tk$token == "NUM_CONST"], "1e+05")
  withr::local_options(list(scipen = 100L))
  tk <- function_tokens(fn_from_text("function(x) x + 100000"))
  expect_identical(tk$text[tk$token == "NUM_CONST"], "1e+05")
})

test_that("no function uses the superassignment operator", {
  detector <- superassign_hits
  expect_true(flags(detector, paste0(
    "function() { n <- 0; f <- function() n <<", "- n + 1; f(); n }")))
  expect_true(flags(detector, paste0(
    "function() { n <- 0; f <- function() n + 1 -", ">> n; f(); n }")))
  expect_true(flags(detector, paste0(
    "function() { n <- 0; f <- function() assign(\"n\", n + 1, ",
    "inherits = TRUE); f(); n }")))
  expect_false(flags(detector, paste0(
    "function() { state <- new.env(parent = emptyenv()); state$n <- 0; ",
    "f <- function() state$n <- state$n + 1; f(); ",
    "get(\"n\", envir = state, inherits = TRUE) }")))

  expect_identical(sweep_namespace(ns_tokens, detector), character(0))
})

test_that("no function sets the warn option", {
  detector <- warn_option_hits
  expect_true(flags(detector, paste0(
    "function(x) { old <- options(warn = -1); on.exit(options(old)); ",
    "as.numeric(x) }")))
  expect_true(flags(detector,
                    "function(x) { base::options(digits = 3, warn = -1); x }"))
  expect_true(flags(detector,
                    "function(x) withr::with_options(list(warn = -1), as.numeric(x))"))
  expect_true(flags(detector, paste0(
    "function(x) { quiet <- list(warn = -1); old <- options(quiet); ",
    "on.exit(options(old)); as.numeric(x) }")))
  expect_true(flags(detector,
                    "function(x) { do.call(options, list(warn = 2)); x }"))
  expect_false(flags(detector, "function(x) suppressWarnings(as.numeric(x))"))
  expect_false(flags(detector, "function() getOption(\"warn\")"))
  expect_false(flags(detector, paste0(
    "function(x) { old <- options(digits = 3); on.exit(options(old)); ",
    "warning(\"warn = \", x) }")))
  expect_false(flags(detector, paste0(
    "function(con) { control <- list(warn = TRUE, n = 3); ",
    "readLines(con, warn = FALSE) }")))

  expect_identical(sweep_namespace(ns_tokens, detector), character(0))
})

test_that("no function has a default file path or writes to a path of its own", {
  detector <- default_path_hits
  expect_true(flags(detector, "function(x, filename = \"results.csv\") x"))
  expect_true(flags(detector, "function(x, save = FALSE, filename) x"))
  expect_true(flags(detector,
                    "function(x, path = file.path(tempdir(), \"x.rds\")) x"))
  expect_true(flags(detector,
                    "function(x, file = base::tempfile(fileext = \".csv\")) x"))
  expect_true(flags(detector, "function(x, dir = \"output\") x"))
  expect_true(flags(detector,
                    "function(x) utils::write.csv(x, \"results.csv\")"))
  expect_true(flags(detector,
                    "function(x) cat(x, file = \"log.txt\", append = TRUE)"))
  expect_true(flags(detector,
                    "function(x) { grDevices::pdf(\"plot.pdf\"); plot(x) }"))
  expect_true(flags(detector,
                    "function(x) lapply(x, function(el, file = \"a.txt\") el)"))
  expect_false(flags(detector, paste0(
    "function(x, filename = NULL, file) { ",
    "if (!is.null(filename)) utils::write.csv(x, filename); ",
    "cat(\"done\", \"\\n\", file = \"\", sep = \"\"); ",
    "system.file(\"extdata\", \"results.rds\", package = \"stats\") }")))
  # An argument named like a path whose default is not one.
  expect_false(flags(detector, paste0(
    "function(x, path = TRUE, dir = c(\"upper\", \"lower\"), ",
    "con = stdout(), file = \"\", save) x")))

  expect_identical(sweep_namespace(ns_tokens, detector), character(0))
})

test_that("no function changes the session without restoring it", {
  detector <- session_state_hits
  expect_true(flags(detector,
                    "function(x) { graphics::par(mfrow = c(1, 2)); plot(x) }"))
  expect_true(flags(detector, "function(x) { options(digits = 3); print(x) }"))
  expect_true(flags(detector, "function(d) { setwd(d); list.files() }"))
  expect_true(flags(detector, "function(p) utils::install.packages(p)"))
  expect_true(flags(detector, "function() quit(save = \"no\")"))
  expect_true(flags(detector, "function(cmd) system(cmd)"))
  # An on.exit() that restores something else, or nothing.
  expect_true(flags(detector, paste0(
    "function(x, con) { on.exit(close(con)); ",
    "graphics::par(mfrow = c(1, 2)); plot(x) }")))
  expect_true(flags(detector, paste0(
    "function(x) { on.exit(NULL); graphics::par(mfrow = c(1, 2)); plot(x) }")))
  expect_true(flags(detector, paste0(
    "function(x) { old <- graphics::par(\"mar\"); on.exit(graphics::par(old)); ",
    "options(digits = 3); x }")))
  expect_true(flags(detector, "function(p) { q <- 1; base::q(\"no\") }"))
  expect_false(flags(detector, paste0(
    "function(x) { old <- graphics::par(mfrow = c(1, 2)); ",
    "on.exit(graphics::par(old), add = TRUE); plot(x) }")))
  expect_false(flags(detector, paste0(
    "function(d) { old <- setwd(d); on.exit({ setwd(old); invisible() }); ",
    "list.files() }")))
  expect_false(flags(detector, paste0(
    "function() { usr <- graphics::par(\"usr\"); all <- graphics::par(); ",
    "keep <- graphics::par(no.readonly = TRUE); ",
    "two <- graphics::par(c(\"usr\", \"mfrow\")); ",
    "c(usr, getOption(\"digits\"), options(\"digits\")[[1L]]) }")))
  # A local function that happens to be named q is not base::q().
  expect_false(flags(detector, paste0(
    "function(p, system = NULL) { q <- function(u) stats::qnorm(u); q(p) }")))

  expect_identical(sweep_namespace(ns_tokens, detector), character(0))
})

test_that("no function writes to the console outside a print method", {
  detector <- console_hits
  expect_true(flags(detector,
                    "function(x) { cat(\"Stopping rule met\\n\"); x }"))
  expect_true(flags(detector,
                    "function(x) { message(\"Collect more data.\"); invisible(x) }"))
  expect_true(flags(detector, "function(x) { print(x); invisible(x) }"))
  expect_true(flags(detector, "function(x) { utils::str(x); invisible(x) }"))
  expect_true(flags(detector,
                    "function(x) { print.default(unclass(x)); invisible(x) }"))
  # A writer handed to another function to call.
  expect_true(flags(detector, "function(x) invisible(lapply(x, print))"))
  expect_true(flags(detector, "function(x) do.call(cat, as.list(x))"))
  expect_true(flags(detector, "function(x) do.call(\"cat\", as.list(x))"))
  expect_true(flags(detector, "function(x) { show <- message; show(x) }"))
  expect_false(flags(detector, paste0(
    "function(x) { if (x < 0) warning(\"negative\"); ",
    "structure(list(x), class = \"result\") }")))
  # Text that is captured, a flag named print, and a condition's message
  # are not console output.
  expect_false(flags(detector, paste0(
    "function(x, print = FALSE) { lines <- utils::capture.output(print(x)); ",
    "str <- paste(lines, collapse = \"\"); ",
    "tryCatch(if (print) str, error = function(e) e$message, ",
    "message = function(m) NULL) }")))

  registered <- c("print.smsd_check", "print.summary.smsd_check")
  expect_true(may_write("print.smsd_check", registered))
  expect_false(may_write("format.smsd_check", registered))
  expect_false(may_write("summary.smsd_check", registered))
  expect_false(may_write("print_progress", registered))
  expect_false(may_write(".print_notes", registered))
  expect_false(may_write("seq_aipe_smd", registered))

  quiet <- ns_tokens[!may_write(names(ns_tokens), print_methods)]
  expect_identical(sweep_namespace(quiet, detector), character(0))
})

test_that("no function takes more than two cores, by default or from the machine", {
  detector <- core_count_hits
  expect_true(flags(detector, paste0(
    "function(x, f) parallel::mclapply(x, f, ",
    "mc.cores = parallel::detectCores())")))
  expect_true(flags(detector,
                    "function(x, f) parallel::mclapply(x, f, mc.cores = 8)"))
  expect_true(flags(detector, "function(x, cores = 4L) x"))
  expect_true(flags(detector, "function() parallel::makeCluster(6)"))
  expect_true(flags(detector, "function() parallel::makePSOCKcluster(4L)"))
  expect_true(flags(detector, "function() doParallel::registerDoParallel(cores = 3)"))
  expect_true(flags(detector, paste0(
    "function() { old <- options(mc.cores = 12); on.exit(options(old)); 1 }")))
  expect_true(flags(detector,
                    "function() Sys.setenv(\"OMP_NUM_THREADS\" = 16)"))
  # The machine's count reached indirectly.
  expect_true(flags(detector, paste0(
    "function() { count <- parallel::detectCores; count() - 1 }")))
  expect_true(flags(detector,
                    "function() do.call(\"detectCores\", list()) - 1"))
  # Two at most, and a count that is the caller's to give.
  expect_false(flags(detector, paste0(
    "function(x, f, cores = 1L) parallel::mclapply(x, f, mc.cores = cores)")))
  expect_false(flags(detector, "function(x, cores = 2) parallel::makeCluster(2L)"))
  expect_false(flags(detector, paste0(
    "function(x, workers = NULL) { threads <- 1; ",
    "stop(\"detectCores() is not used; pass cores = 4 yourself\") }")))
  # A number above two that counts something else.
  expect_false(flags(detector, "function(x, G = 1000, batch = 5L) rep(x, 4)"))

  expect_identical(sweep_namespace(ns_tokens, detector), character(0))
})
