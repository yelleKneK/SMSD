# Build a site from the committed source, never from ignored local notes.
# Usage: Rscript tools/build_site.R --output=DIR
# DIR must lie outside the repository. Nothing is published by this script.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !grepl("^--output=.+", args)) {
  stop("Supply exactly one --output=DIR outside the repository.", call. = FALSE)
}
root <- normalizePath(".", mustWork = TRUE)
if (!file.exists("DESCRIPTION")) stop("Run from the package root.", call. = FALSE)
changes <- system2("git", c("status", "--porcelain"), stdout = TRUE)
if (!is.null(attr(changes, "status")) || length(changes)) {
  stop("Commit the source changes before building the site.", call. = FALSE)
}
destination <- sub("^--output=", "", args)
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
destination <- normalizePath(destination, mustWork = TRUE)
if (startsWith(paste0(destination, "/"), paste0(root, "/"))) {
  stop("The output directory must lie outside the repository.", call. = FALSE)
}
if (length(list.files(destination, all.files = TRUE, no.. = TRUE))) {
  stop("The output directory must be empty; choose a new directory.", call. = FALSE)
}
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
if (!is.null(attr(commit, "status"))) stop("Cannot identify the source commit.")
archive <- tempfile(fileext = ".tar")
source <- tempfile(pattern = "site-source-")
dir.create(source)
status <- system2("git", c("archive", "--format=tar",
                           paste0("--output=", shQuote(archive)), "HEAD"))
if (status != 0L) stop("Cannot archive the source commit.")
utils::untar(archive, exdir = source)
pkgdown::build_site(source, install = TRUE, new_process = TRUE, preview = FALSE,
                    override = list(destination = destination))
pages <- list.files(destination, pattern = "[.]html$")
allowed <- c("index.html", "404.html", "authors.html", "LICENSE.html", "LICENSE-text.html")
unexpected <- setdiff(pages, allowed)
if (length(unexpected)) {
  stop("Unexpected top-level pages: ", paste(unexpected, collapse = ", "), call. = FALSE)
}
writeLines(c(paste("Source commit:", commit),
             paste("Built:", format(Sys.time(), tz = "UTC", usetz = TRUE))),
            file.path(destination, "source-record.txt"))
unlink(c(archive, source), recursive = TRUE)
cat("Site built from", commit, "in", destination, "\n")
