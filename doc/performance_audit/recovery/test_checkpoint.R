suppressPackageStartupMessages(library(testthat))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
source(file.path(root, "doc/performance_audit/recovery/checkpoint.R"))
base <- file.path(root, "doc/statistical_toolkit_runs/20260916/performance-audit/recovery/checkpoint-tests")
dir.create(base, recursive = TRUE, showWarnings = FALSE)
artifacts <- if (length(args) >= 2L) args[[2L]] else tempfile("test-", tmpdir = base)
dir.create(artifacts, recursive = TRUE, showWarnings = FALSE)

crash_job <- function() {
  if (!file.exists(file.path(artifacts, "allow_completion"))) quit(save = "no", status = 89L)
  42L
}
if ("--child" %in% args) {
  audit_checkpoint_run("crash", list(fixture = "intentional-process-exit"),
                       crash_job, artifacts)
  quit(status = 1L)
}

test_that("completed jobs retain warnings, cache data, and distinguish specifications", {
  count <- 0L
  run <- function() { count <<- count + 1L; warning("retained warning"); list(answer = 17L) }
  a <- audit_checkpoint_run("success", list(data = 1:3, seed = 123L), run, artifacts)
  b <- audit_checkpoint_run("success", list(data = 1:3, seed = 123L), run, artifacts)
  expect_identical(a$status, "COMPLETE")
  expect_identical(a$warnings, "retained warning")
  expect_identical(b$result, a$result)
  expect_true(b$cached)
  expect_identical(count, 1L)
  c <- audit_checkpoint_run("success", list(data = 1:3, seed = 124L), run, artifacts)
  expect_false(c$cached)
  expect_identical(count, 2L)
  expect_false(identical(a$fingerprint, c$fingerprint))
  audit_atomic_rds("corrupt fixture payload", file.path(a$attempt, "result.rds"))
  expect_error(audit_checkpoint_run("success", list(data = 1:3, seed = 123L),
                                  run, artifacts), "absent or corrupt")
  expect_error(audit_checkpoint_run("../escape", list(), run, artifacts), "directory name")
})

test_that("failures persist and are retried only explicitly", {
  count <- 0L
  run <- function() { count <<- count + 1L; warning("before failure"); stop("fixture failure") }
  a <- audit_checkpoint_run("failed", list(seed = 5L), run, artifacts)
  b <- audit_checkpoint_run("failed", list(seed = 5L), run, artifacts)
  expect_identical(a$status, "FAILED")
  expect_identical(a$error, "fixture failure")
  expect_identical(a$warnings, "before failure")
  expect_true(b$cached)
  expect_identical(count, 1L)
  c <- audit_checkpoint_run("failed", list(seed = 5L), run, artifacts, retry_failed = TRUE)
  expect_identical(c$status, "FAILED")
  expect_identical(count, 2L)
  expect_false(identical(a$attempt, c$attempt))
  expect_true(file.exists(file.path(a$attempt, "finished.rds")))
})

test_that("an exited process resumes only explicitly and retains its interruption", {
  script <- file.path(root, "doc/performance_audit/recovery/test_checkpoint.R")
  exit <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    args = shQuote(c(script, root, artifacts, "--child")),
    stdout = file.path(artifacts, "child.stdout"), stderr = file.path(artifacts, "child.stderr")))
  expect_identical(as.integer(exit), 89L)
  expect_error(audit_checkpoint_run("crash", list(fixture = "intentional-process-exit"),
    crash_job, artifacts), "Unfinished attempt")
  file.create(file.path(artifacts, "allow_completion"))
  recovered <- audit_checkpoint_run("crash", list(fixture = "intentional-process-exit"),
    crash_job, artifacts, resume_incomplete = TRUE)
  expect_identical(recovered$status, "COMPLETE")
  expect_identical(recovered$result, 42L)
  cached <- audit_checkpoint_run("crash", list(fixture = "intentional-process-exit"),
    crash_job, artifacts)
  expect_true(cached$cached)
  records <- list.files(file.path(artifacts, "crash"), pattern = "finished.rds$",
                        recursive = TRUE, full.names = TRUE)
  statuses <- vapply(records, function(p) readRDS(p)$status, character(1))
  expect_setequal(statuses, c("INTERRUPTED", "COMPLETE"))
  expect_identical(length(records), 2L)
})

cat("Checkpoint tests completed; fixtures retained at:", artifacts, "\n")
