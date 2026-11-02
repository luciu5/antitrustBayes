.pilot_audit_path <- testthat::test_path("..", "..", "doc", "validation", "nutpie_pilot_diagnostic_audit.R")
testthat::skip_if_not(file.exists(.pilot_audit_path), "Repository-only no-HMC pilot audit.")
.pilot_audit_env <- new.env(parent = globalenv())
sys.source(.pilot_audit_path, envir = .pilot_audit_env)

test_that("both Nutpie adaptations use native reference constants without changing CmdStan", {
  manifest <- list(original_root = "original", reference_fits = list(
    CES_CS = list(path = "sealed-cmdstan-fit.rds")))
  reference <- .pilot_audit_env$nutpie_pilot_audit_reference
  expect_identical(reference("CES_CS", "cmdstanr", manifest), "sealed-cmdstan-fit.rds")
  expected <- file.path("original", "CES_CS__nutpie_diag", "fit.rds")
  expect_identical(reference("CES_CS", "nutpie_diag", manifest), expected)
  expect_identical(reference("CES_CS", "nutpie_low_rank", manifest), expected)
})

test_that("precision correction preserves exact constant identity and missing-diagnostic failures", {
  helper <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "doc", "validation", "nutpie_certification_pilot.R"), envir = helper)
  draws <- array(c(rep(1 / 3, 8), 1:8), c(4, 2, 2),
    dimnames = list(NULL, NULL, c("fixed", "varying")))
  summary <- data.frame(variable = c("fixed", "varying"), rhat = c(NA, 1.005),
    ess_bulk = c(NA, 700), ess_tail = c(NA, 600))
  old <- helper$nutpie_pilot_posterior_availability(draws, summary, list(fixed = .33333333))
  corrected <- helper$nutpie_pilot_posterior_availability(draws, summary, list(fixed = 1 / 3))
  expect_identical(old$diagnostics_available, c(FALSE, TRUE))
  expect_identical(corrected$diagnostics_available, c(TRUE, TRUE))
  expect_identical(.pilot_audit_env$nutpie_pilot_audit_reason(draws, old),
    c("constant absent from reference or not exactly equal", ""))
  summary$rhat[2L] <- NA_real_
  unavailable <- helper$nutpie_pilot_posterior_availability(draws, summary, list(fixed = 1 / 3))
  expect_identical(unavailable$diagnostics_available, c(TRUE, FALSE))
  draws[, , "fixed"] <- 1 / 3 + 1e-14
  changed <- helper$nutpie_pilot_posterior_availability(draws, summary, list(fixed = 1 / 3))
  expect_identical(changed$diagnostics_available, c(FALSE, FALSE))
})

test_that("diagnostic audit cannot run accidentally when sourced", {
  withr::local_envvar(AB_NUTPIE_DIAGNOSTIC_AUDIT_EXECUTE = "false")
  expect_error(.pilot_audit_env$nutpie_pilot_diagnostic_audit("absent-repository"),
    "AB_NUTPIE_DIAGNOSTIC_AUDIT_EXECUTE=true")
})
