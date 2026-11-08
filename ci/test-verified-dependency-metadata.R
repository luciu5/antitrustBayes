#!/usr/bin/env Rscript

# Exercise the exact R heredocs embedded in qa.yml. This catches drift between
# the hosted verification code and the regression checks below.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste(
    "usage: Rscript --vanilla ci/test-verified-dependency-metadata.R",
    "<antitrust-source> <coordination-source> <trade-source>",
    "<antitrust-sha> <coordination-sha>"
  ), call. = FALSE)
}

packages <- c("antitrust", "coordination", "trade")
source_dirs <- args[seq_along(packages)]
names(source_dirs) <- packages
antitrust_sha <- tolower(args[[4L]])
coordination_sha <- tolower(args[[5L]])
if (!grepl("^[[:xdigit:]]{40}$", antitrust_sha) ||
    !grepl("^[[:xdigit:]]{40}$", coordination_sha)) {
  stop("The selected structural dependency refs must be full commit SHAs.", call. = FALSE)
}

workflow <- readLines(".github/workflows/qa.yml", warn = FALSE)
heredoc_starts <- grep(
  'Rscript --vanilla - "${checkout_dir}" "${dependency}" "$RELEASE_PIN_BUILD"',
  workflow,
  fixed = TRUE
)
if (length(heredoc_starts) != 2L) {
  stop("Expected both structural dependency verification heredocs in qa.yml.", call. = FALSE)
}

extract_heredoc <- function(start) {
  marker <- which(grepl("<<'RSCRIPT'", workflow[seq.int(start, length(workflow))], fixed = TRUE))[1L]
  if (is.na(marker)) stop("Could not find R heredoc marker after workflow invocation.", call. = FALSE)
  first <- start + marker
  end <- which(trimws(workflow[seq.int(first, length(workflow))]) == "RSCRIPT")[1L]
  if (is.na(end)) stop("Could not find the end of a workflow R heredoc.", call. = FALSE)
  end <- first + end - 1L
  code <- workflow[seq.int(first, end - 1L)]
  sub("^ {10}", "", code)
}

workflow_code <- lapply(heredoc_starts, extract_heredoc)
if (!identical(workflow_code[[1L]], workflow_code[[2L]])) {
  stop("The two hosted dependency verification heredocs have drifted.", call. = FALSE)
}

for (package in packages) {
  if (!file.exists(file.path(source_dirs[[package]], "DESCRIPTION"))) {
    stop("Missing selected source DESCRIPTION for ", package, call. = FALSE)
  }
}

installed_description <- lapply(packages, utils::packageDescription)
names(installed_description) <- packages
source_description <- lapply(source_dirs, function(path) {
  read.dcf(file.path(path, "DESCRIPTION"))[1L, ]
})

# The selected antitrust source and installation intentionally omit Remotes.
# Confirm that the historical source DCF expression fails before exercising
# the fixed workflow code against both missing-field representations.
if ("Remotes" %in% names(source_description$antitrust) ||
    "Remotes" %in% names(installed_description$antitrust)) {
  stop("The antitrust regression fixture must have no Remotes field.", call. = FALSE)
}
old_source_access_fails <- tryCatch({
  source_description$antitrust[["Remotes"]]
  FALSE
}, error = function(error) TRUE)
if (!old_source_access_fails) {
  stop("The regression fixture no longer reproduces the old missing-Dcf-field failure.", call. = FALSE)
}

script_file <- tempfile(fileext = ".R")
on.exit(unlink(script_file), add = TRUE)

run_workflow_verifier <- function(code, source_dir, package, release_pin,
                                  anti_sha = antitrust_sha,
                                  coord_sha = coordination_sha) {
  writeLines(code, script_file, useBytes = TRUE)
  command_args <- c(
    "--vanilla", shQuote(script_file), shQuote(source_dir), shQuote(package),
    shQuote(if (release_pin) "true" else "false"), shQuote(anti_sha), shQuote(coord_sha)
  )
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    args = command_args,
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(status = as.integer(status), output = output)
}

expect_success <- function(result, label) {
  if (result$status != 0L) {
    stop(label, " unexpectedly failed:\n", paste(result$output, collapse = "\n"), call. = FALSE)
  }
  cat("PASS:", label, "\n")
}

expect_rejection <- function(result, label) {
  if (result$status == 0L) {
    stop(label, " was accepted by the workflow verifier.", call. = FALSE)
  }
  cat("PASS:", label, "\n")
}

wrong_sha <- function(sha) {
  paste0(if (substr(sha, 1L, 1L) == "0") "1" else "0", substr(sha, 2L, 40L))
}

for (copy_index in seq_along(workflow_code)) {
  code <- workflow_code[[copy_index]]
  label_prefix <- paste0("qa.yml verifier ", copy_index)

  for (package in packages) {
    result <- run_workflow_verifier(code, source_dirs[[package]], package, TRUE)
    expect_success(result, paste(label_prefix, "accepts selected", package, "source/install metadata"))
  }

  # Wrong exact pins must still fail for both packages whose Remotes are
  # expected to identify the selected structural dependencies.
  result <- run_workflow_verifier(
    code, source_dirs[["coordination"]], "coordination", TRUE,
    anti_sha = wrong_sha(antitrust_sha)
  )
  expect_rejection(result, paste(label_prefix, "rejects a wrong antitrust pin in coordination"))

  result <- run_workflow_verifier(
    code, source_dirs[["trade"]], "trade", TRUE,
    coord_sha = wrong_sha(coordination_sha)
  )
  expect_rejection(result, paste(label_prefix, "rejects a wrong coordination pin in trade"))

  # A source Version that differs from the installed package must not pass.
  for (package in packages) {
    mismatch_dir <- tempfile(pattern = paste0("source-version-mismatch-", package, "-"))
    dir.create(mismatch_dir)
    mismatch_dcf <- matrix(
      source_description[[package]],
      nrow = 1L,
      dimnames = list(NULL, names(source_description[[package]]))
    )
    mismatch_dcf[1L, "Version"] <- paste0(mismatch_dcf[1L, "Version"], ".ci-mismatch")
    write.dcf(mismatch_dcf, file = file.path(mismatch_dir, "DESCRIPTION"))
    result <- run_workflow_verifier(code, mismatch_dir, package, FALSE)
    expect_rejection(result, paste(label_prefix, "rejects an installed/source version mismatch for", package))
    unlink(mismatch_dir, recursive = TRUE)
  }
}

cat("Both qa.yml verification snippets passed the metadata regressions.\n")
