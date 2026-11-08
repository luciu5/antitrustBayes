#!/usr/bin/env Rscript
# Permanent A: these checks guard scientific/provenance failures that ordinary
# policy-equilibrium tests cannot detect. They require no fits or native sampler.
source("doc/validation/saved_gate_policy_replay.R")

testthat::test_that("A: replay rejects changed likelihood or reconstruction code despite matching package version", {
  paths <- c("DESCRIPTION", sort(list.files("R", "[.]R$", full.names = TRUE)),
             sort(list.files("inst/stan", "[.]stan$", full.names = TRUE)))
  hashes <- stats::setNames(vapply(paths, saved_gate_hash, character(1)), paths)
  for (file in c("R/posterior_policy_models.R", "R/priors.R", "inst/stan/ces_supply.stan")) {
    changed <- hashes; changed[file] <- paste(rep("0", 64), collapse = "")
    testthat::expect_error(saved_gate_check_sources(changed, TRUE), "Sealed inference source changed")
  }
  changed <- hashes; changed["DESCRIPTION"] <- paste(rep("0", 64), collapse = "")
  testthat::expect_error(saved_gate_check_sources(changed), "DESCRIPTION changed")
  testthat::expect_false(saved_gate_check_sources(changed, TRUE)$identical[1L])
})

testthat::test_that("A: replay cannot silently reuse a refit or drop a failed positive-support case", {
  path <- tempfile(); on.exit(unlink(path))
  saveRDS(list(draws = 1:3), path); seal <- saved_gate_hash(path)
  saveRDS(list(draws = 4:6), path)
  testthat::expect_error(saved_gate_check_artifact(path, seal), "artifact seal differs")
  calls <- list(logit_C_M4T4 = list(priors = list(alpha = 1)), ces_C_M4T4_finite_noise = list())
  m <- list(status = "sampling_complete", model_specification = 14L,
            calls = calls, calls_hash = digest::digest(calls, algo = "sha256"))
  m$calls$logit_C_M4T4$priors$alpha <- 2
  testthat::expect_error(saved_gate_check_calls(m), "fixture/prior/sampler call identity changed")
  m$calls <- calls[1L]; m$calls_hash <- digest::digest(m$calls, algo = "sha256")
  testthat::expect_error(saved_gate_check_calls(m), "incomplete fits cannot be omitted")
})

testthat::test_that("A: replay refuses paths that could overwrite inference evidence or a prior replay", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  saved <- file.path(root, "saved"); dir.create(saved)
  testthat::expect_error(saved_gate_output(saved, saved), "separate from the read-only")
  testthat::expect_error(saved_gate_output(saved, file.path(saved, "replay")), "separate from the read-only")
  output <- file.path(root, "replay"); dir.create(output); writeLines("retained", file.path(output, "failed.log"))
  testthat::expect_error(saved_gate_output(saved, output), "Preserve existing replay output")
  testthat::expect_identical(readLines(file.path(output, "failed.log")), "retained")
})
