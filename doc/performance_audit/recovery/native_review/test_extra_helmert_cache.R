suppressPackageStartupMessages(library(testthat))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  root <- normalizePath(args[[1L]], mustWork = TRUE)
} else {
  root <- normalizePath("doc/performance_audit/recovery/native_review",
                       mustWork = TRUE)
}
source(file.path(root, "cache_helpers.R"))
artifact_root <- tempfile("native-review-cache-")
dir.create(artifact_root, recursive = TRUE)

test_that("cache identity separates generations and reuses exact completion", {
  calls <- 0L
  run <- function() {
    calls <<- calls + 1L
    list(status = "ok", value = calls)
  }
  spec <- list(schema = 2L, input_hash = "input-a", kernel_hash = "kernel-a",
               raw_points_hash = NA_character_)
  first <- native_review_cache(file.path(artifact_root, "case"), spec, run)
  cached <- native_review_cache(file.path(artifact_root, "case"), spec, run)
  expect_false(first$reused)
  expect_true(cached$reused)
  expect_identical(calls, 1L)
  changed <- native_review_cache(
    file.path(artifact_root, "case"),
    modifyList(spec, list(kernel_hash = "kernel-b")), run)
  expect_false(changed$reused)
  expect_false(identical(first$generation, changed$generation))
  expect_identical(calls, 2L)
})

test_that("failed attempts persist and force creates a new attempt", {
  calls <- 0L
  run <- function() {
    calls <<- calls + 1L
    list(status = "error", error = "fixture failure", value = calls)
  }
  spec <- list(schema = 2L, fixture = "failed")
  first <- native_review_cache(file.path(artifact_root, "failed"), spec, run)
  cached <- native_review_cache(file.path(artifact_root, "failed"), spec, run)
  retried <- native_review_cache(file.path(artifact_root, "failed"), spec, run,
                                 force = TRUE)
  expect_identical(first$state, "FAILED")
  expect_true(cached$reused)
  expect_false(retried$reused)
  expect_false(identical(first$attempt, retried$attempt))
  expect_identical(calls, 2L)
})

test_that("candidate raw replay requires intact provenance", {
  u <- c(-.2, .3, .4)
  source_hashes <- list(runner = "r", fixture = "f", cache_helper = "c")
  settings <- list(iter = 1L, warmup = 0L)
  seeds <- list(target = 58003L, generated_quantities = 58005L)
  entry <- list(case = "fixture", u = u,
    u_hash = native_review_hash(u), input_hash = "data",
    kernel_code_hash = "candidate-kernel",
    source_hashes = source_hashes, settings = settings,
    schema = 2L, role = "candidate", seeds = seeds)
  expect_silent(native_review_validate_raw(entry, "fixture", "data",
                                           source_hashes, settings, seeds))
  expect_error(native_review_validate_raw(
    modifyList(entry, list(u = c(u, .5))), "fixture", "data",
    source_hashes, settings, seeds), "provenance mismatch")
  expect_error(native_review_validate_raw(
    modifyList(entry, list(input_hash = "changed")), "fixture", "data",
    source_hashes, settings, seeds), "provenance mismatch")
  expect_error(native_review_validate_raw(
    modifyList(entry, list(settings = list(iter = 2L, warmup = 0L))),
    "fixture", "data", source_hashes, settings, seeds),
    "provenance mismatch")
})

cat("Extra native-review cache tests completed; artifacts retained at:",
    artifact_root, "\n")
