.nutpie_pilot_env <- new.env(parent = globalenv())
.nutpie_pilot_path <- testthat::test_path("..", "..", "doc", "validation", "nutpie_certification_pilot.R")
testthat::skip_if_not(file.exists(.nutpie_pilot_path), "Repository-only certification experiment.")
sys.source(.nutpie_pilot_path, envir = .nutpie_pilot_env)

test_that("certification pilot cannot start inference without explicit opt-in", {
  withr::local_envvar(AB_NUTPIE_PILOT_EXECUTE = "false")
  expect_error(.nutpie_pilot_env$nutpie_pilot_worker("absent-task.rds"), "AB_NUTPIE_PILOT_EXECUTE=true")
})

test_that("one bad chain or unavailable diagnostic prevents pilot eligibility", {
  thresholds <- .nutpie_pilot_env$nutpie_pilot_design()$diagnostic_thresholds
  d <- data.frame(divergences = rep(0, 4), depth_hits = 0, bfmi = .8,
    max_rhat = 1.005, min_bulk_ess = 800, min_tail_ess = 600, unavailable_posterior_diagnostics = 0)
  check <- function(x) .nutpie_pilot_env$nutpie_pilot_diagnostic_pass(x, thresholds)
  expect_true(check(d))
  bad_values <- c(divergences = 1, depth_hits = 1, bfmi = .2, max_rhat = 1.02,
    min_bulk_ess = 399, min_tail_ess = 399, unavailable_posterior_diagnostics = 1)
  for (field in names(bad_values)) {
    bad <- d; bad[4L, field] <- bad_values[[field]]
    expect_false(check(bad), info = field)
    bad[4L, field] <- NA_real_
    expect_false(check(bad), info = paste("missing", field))
  }
  expect_false(check(d[-4L, ]))
})

test_that("pilot collection preserves failures and never certifies a successful subset", {
  directory <- withr::local_tempdir()
  design <- .nutpie_pilot_env$nutpie_pilot_design()
  manifest <- list(design = design)
  for (i in which(design$tasks$variant == "nutpie_diag")) {
    path <- file.path(directory, design$tasks$id[i]); dir.create(path)
    saveRDS(list(status = "completed", identity_equal = TRUE, diagnostics_passed = TRUE),
      file.path(path, "result.rds"))
  }
  value <- .nutpie_pilot_env$nutpie_pilot_collect(manifest, directory, directory)
  expect_identical(value$status, "incomplete")
  expect_true(value$eligible[["nutpie_diag"]])
  expect_false(value$eligible[["cmdstanr"]])
  path <- file.path(directory, "CF_pooled__nutpie_diag", "result.rds")
  failed <- list(status = "failed", reason = "failed posterior", identity_equal = TRUE)
  saveRDS(failed, path)
  value <- .nutpie_pilot_env$nutpie_pilot_collect(manifest, directory, directory)
  expect_false(value$eligible[["nutpie_diag"]])
  expect_equal(nrow(value$tasks), 9L)
  expect_equal(sum(value$tasks$status == "not_run"), 6L)
  expect_identical(readRDS(path), failed)
  expect_false(value$public_evidence_gate_changed)
})

test_that("prepared pilot rejects runtime changes before inspecting or mutating fits", {
  manifest <- list(design = .nutpie_pilot_env$nutpie_pilot_design(), runtime = "sealed")
  original <- .nutpie_pilot_env$nutpie_pilot_runtime
  withr::defer(assign("nutpie_pilot_runtime", original, envir = .nutpie_pilot_env))
  .nutpie_pilot_env$nutpie_pilot_runtime <- function(repository) "changed"
  expect_error(.nutpie_pilot_env$nutpie_pilot_assert_manifest(manifest, "."), "changed after preparation")
})

test_that("constant fields are allowed but cannot hide unavailable nonconstant diagnostics", {
  draws <- array(c(rep(1, 8), 1:8), dim = c(4, 2, 2),
    dimnames = list(NULL, NULL, c("constant", "varying")))
  summary <- data.frame(variable = c("constant", "varying"), rhat = c(NA, 1.005),
    ess_bulk = c(NA, 800), ess_tail = c(NA, 700))
  constants <- .nutpie_pilot_env$nutpie_pilot_constant_reference(draws)
  check <- function(x, s) .nutpie_pilot_env$nutpie_pilot_posterior_availability(x, s, constants)
  expect_true(all(check(draws, summary)$diagnostics_available))
  summary$rhat[2L] <- NA_real_
  expect_identical(check(draws, summary)$diagnostics_available, c(TRUE, FALSE))
  summary$rhat[2L] <- Inf
  expect_identical(check(draws, summary)$diagnostics_available, c(TRUE, FALSE))
  draws[1L, 1L, 1L] <- NA_real_
  expect_false(check(draws, summary)$diagnostics_available[1L])
  draws[,,1L] <- 2
  expect_false(check(draws, summary)$diagnostics_available[1L])
  draws[,,1L] <- 1
  draws[,,2L] <- 0
  summary[2L, c("rhat", "ess_bulk", "ess_tail")] <- c(1, 800, 700)
  expect_false(check(draws, summary)$diagnostics_available[2L])
})
