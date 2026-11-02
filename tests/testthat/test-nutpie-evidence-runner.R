.nutpie_evidence_runner_path <- testthat::test_path(
  "..", "..", "doc", "validation", "nutpie_evidence_validation.R")
testthat::skip_if_not(file.exists(.nutpie_evidence_runner_path),
  "Research runner checks require the repository doc/ tree (excluded from the installed package).")

.nutpie_evidence_runner <- function() {
  env <- new.env(parent = globalenv())
  source(.nutpie_evidence_runner_path, local = env)
  env
}

test_that("evidence runner has a fixed full kernel/model plan without inference on source", {
  env <- .nutpie_evidence_runner()
  tasks <- env$nutpie_evidence_tasks()
  expect_identical(nrow(tasks), 15L)
  expect_identical(sum(tasks$stage == "sample"), 4L)
  expect_identical(sum(tasks$stage == "bridge"), 6L)
  expect_identical(sum(tasks$stage == "target"), 4L)
  expect_setequal(tasks$kernel[tasks$stage == "target"],
    c("logit_supply", "logit_stackelberg", "ces_supply", "ces_stackelberg"))
  expect_identical(anyDuplicated(tasks$id), 0L)
})

test_that("evidence budgets are durable and separate from historical sampling clocks", {
  env <- .nutpie_evidence_runner()
  now <- as.POSIXct("2026-10-19 12:00:00", tz = "UTC")
  expect_error(env$nutpie_evidence_budget(list(), now), "invocation clock is required")
  expect_identical(env$nutpie_evidence_budget(list(started_at = now - 100, deadline = now + 60), now), 55)
  expect_error(env$nutpie_evidence_budget(list(started_at = now - 14400, deadline = now), now), "exhausted")
  directory <- tempfile("evidence-clock-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  path <- env$nutpie_evidence_invocation(directory, seconds = 60, now = now)
  invocation <- readRDS(path)
  expect_identical(invocation$started_at, now)
  expect_identical(invocation$deadline, now + 60)
  expect_identical(env$nutpie_evidence_budget(invocation, now), 55)
  expect_error(env$nutpie_evidence_invocation(directory, seconds = 14401), "between 10 and 14400")
})

test_that("native target checks detect omitted constants and Jacobians", {
  env <- .nutpie_evidence_runner()
  normalized <- c(-10, -12)
  without_jacobian <- c(-9, -10)
  checks <- env$nutpie_evidence_point_checks(normalized, normalized, normalized,
    without_jacobian, without_jacobian, without_jacobian)
  expect_true(all(checks$passed))
  expect_equal(checks$jacobian, c(-1, -2))
  bad <- env$nutpie_evidence_point_checks(normalized, normalized + log(2), normalized,
    without_jacobian, without_jacobian, without_jacobian)
  expect_false(any(bad$passed))
  bad <- env$nutpie_evidence_point_checks(normalized, normalized, normalized,
    without_jacobian, without_jacobian, without_jacobian + .2)
  expect_false(any(bad$passed))
  expect_error(env$nutpie_evidence_point_checks(NA_real_, 1, 1, 1, 1, 1), "finite paired")
})

test_that("shared evidence arguments retain B observations and priors across candidates", {
  env <- .nutpie_evidence_runner()
  original <- list(prices = c(1, 2), quantities = c(3, 4), margins = c(.2, .3),
    priors = list(prior_alpha_sd = .7), supply_model = "bertrand", use_rho = TRUE)
  manifest <- list(fixtures = list(B = list(arguments = original)))
  for (backend in c("cmdstanr", "nutpie")) for (model in c("MC", "C")) {
    got <- env$nutpie_evidence_shared_arguments(manifest, model, backend, tempdir())
    expect_identical(got[c("prices", "quantities", "margins", "priors", "use_rho")],
                     original[c("prices", "quantities", "margins", "priors", "use_rho")])
    expect_identical(got$stan_options$chains, 4L)
    expect_identical(got$stan_options$iter, 2000L)
    expect_identical(got$stan_options$warmup, 1000L)
  }
})

test_that("bridge equivalence retains every repetition and rejects low precision", {
  env <- .nutpie_evidence_runner()
  reps <- lapply(c(-10.01, -10, -9.99), function(x) list(logml = x, mcse_logml = .01))
  summary <- env$nutpie_evidence_summary(reps)
  expect_equal(summary$logml, -10)
  expect_error(env$nutpie_evidence_summary(reps[1:2]), "three planned")
  expect_error(env$nutpie_evidence_summary(c(reps[1:2], list(list(logml = NA_real_)))), "must succeed")
  reference <- list(model = "B", identity = list(data = "same"), summary = summary)
  alternative <- reference
  alternative$summary$logml <- -9.99
  expect_true(env$nutpie_evidence_compare(reference, alternative)$passed)
  alternative$summary$uncertainty <- 1
  expect_false(env$nutpie_evidence_compare(reference, alternative)$passed)
  alternative <- reference
  alternative$summary$mcse_available <- FALSE
  expect_false(env$nutpie_evidence_compare(reference, alternative)$passed)
  alternative$identity$data <- "different"
  expect_error(env$nutpie_evidence_compare(reference, alternative), "different source/data")
})

test_that("analytic evidence reference has unit mass including both Jacobians", {
  env <- .nutpie_evidence_runner()
  scale_integral <- integrate(function(x) exp(-exp(x) + x), -Inf, Inf)$value
  share_integral <- integrate(function(x) plogis(x) * plogis(x, lower.tail = FALSE), -Inf, Inf)$value
  expect_equal(scale_integral * share_integral, 1, tolerance = 1e-7)
  expect_equal(env$nutpie_evidence_reference_density(c(0, 0)), -1 - log(4))
})

test_that("evidence collection with missing fits retains the entire planned support", {
  env <- .nutpie_evidence_runner()
  root <- tempfile("nutpie-evidence-")
  dir.create(file.path(root, "evidence"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  report <- file.path(root, "reports")
  dir.create(report)
  protocol <- env$nutpie_evidence_collect(list(tasks = env$nutpie_evidence_tasks()), root, report)
  expect_identical(protocol$status, "incomplete")
  expect_false(protocol$public_gate_changed)
  status <- read.csv(file.path(report, "nutpie_evidence_status.csv"))
  expect_identical(nrow(status), 15L)
  expect_true(all(status$status == "not_run"))
  expect_error(env$nutpie_evidence_require_targets(root), "All four native target")
})

test_that("evidence precision cannot certify poor or missing HMC diagnostics", {
  env <- .nutpie_evidence_runner()
  d <- data.frame(divergences = rep(0, 4), depth_hits = rep(0, 4), bfmi = rep(.5, 4),
    max_rhat = rep(1, 4), min_bulk_ess = rep(200, 4), min_tail_ess = rep(200, 4))
  expect_true(env$nutpie_evidence_diagnostic_pass(list(diagnostics = list(summary = d))))
  d$max_rhat[[2L]] <- 1.1
  expect_false(env$nutpie_evidence_diagnostic_pass(list(diagnostics = list(summary = d))))
  d$max_rhat[[2L]] <- NA_real_
  expect_false(env$nutpie_evidence_diagnostic_pass(list(diagnostics = list(summary = d))))
})
