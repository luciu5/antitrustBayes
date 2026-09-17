# Run the pure R migration-policy checks without an installed Bayes native library.
# Full native target restoration remains a separate gate.
suppressPackageStartupMessages(library(testthat))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
test_env <- new.env(parent = globalenv())
sys.source(file.path(root, "R", "bma.R"), envir = test_env)
test_file(file.path(root, "tests", "testthat", "test-bridge-performance-migration.R"),
          env = test_env, reporter = "summary", stop_on_failure = TRUE)
