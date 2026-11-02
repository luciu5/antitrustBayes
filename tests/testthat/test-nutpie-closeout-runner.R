closeout_runner <- function() {
  path <- testthat::test_path("..", "..", "doc", "validation", "nutpie_closeout_validation.R")
  if (!file.exists(path)) testthat::skip("Research runner is excluded from the built package.")
  environment <- new.env(parent = globalenv())
  sys.source(path, envir = environment)
  environment
}

test_that("closeout schedules every feasible original k state without retuning", {
  runner <- closeout_runner()
  design <- runner$nutpie_closeout_design()
  expect_identical(design$k_support, 0:4)
  expect_equal(unname(design$k_prior), c(.10, .15, .25, .30, .20))
  expect_equal(length(design$tasks), 42L)
  samples <- Filter(function(x) x$stage == "sample", design$tasks)
  expect_length(samples, 8L)
  expect_setequal(vapply(samples, `[[`, integer(1), "k"), c(0L, 1L, 3L, 4L))
  expect_identical(design$controls$control, list(adapt_delta = .95, max_treedepth = 13L))
  policies <- Filter(function(x) x$stage == "policy", design$tasks)
  expect_length(policies, 10L)
  expect_equal(table(vapply(policies, `[[`, integer(1), "k")), rep(2L, 5L), ignore_attr = TRUE)
})

test_that("closeout uncertainty retains failed repetitions and precision limits", {
  runner <- closeout_runner()
  partial <- runner$nutpie_closeout_bridge_summary(list(list(logml = 1, mcse_logml = .01),
    NULL, list(logml = 1.02, mcse_logml = .02)))
  expect_equal(partial$logml, 1.01)
  expect_equal(partial$completed, 2L)
  expect_equal(partial$requested, 3L)
  expect_true(is.na(partial$repetitions[2L]))
  reference <- alternative <- partial
  comparison <- runner$nutpie_closeout_compare(reference, alternative, "partial")
  expect_false(comparison$precision_passed)
  expect_true(comparison$agreement_passed)
  reference$completed <- alternative$completed <- 3L
  reference$uncertainty <- alternative$uncertainty <- .1
  comparison <- runner$nutpie_closeout_compare(reference, alternative, "imprecise")
  expect_false(comparison$precision_passed)
  expect_true(comparison$agreement_passed)
  comparison <- runner$nutpie_closeout_compare(reference, alternative, "mixed", shared_sources = TRUE)
  expect_equal(comparison$combined_mcse, .2)
  expect_equal(comparison$uncertainty_combination, "sum_shared_sources")
})

test_that("every public closeout comparison receives a distinct prespecified RNG seed", {
  runner <- closeout_runner()
  design <- runner$nutpie_closeout_design()
  paths <- expand.grid(model = c("B", "MC", "C", "CF", "CS"), backend = c("cmdstanr", "nutpie"))
  manifest <- list(design = design, original_paths = stats::setNames(rep("fake-fit.rds", nrow(paths)),
    paste(paths$model, paths$backend, sep = "_")))
  runner$readRDS <- function(path) NULL
  runner$nutpie_closeout_require <- function(...) NULL
  runner$nutpie_closeout_component_path <- function(...) "fake-component.rds"
  runner$nutpie_closeout_component_fit <- function(...) NULL
  runner$nutpie_closeout_public <- function(fits, directory, universe = FALSE, seed) seed
  index <- which(vapply(design$tasks, function(x) x$stage %in% c("public", "fixed", "component"), logical(1)))
  seeds <- vapply(design$tasks[index], function(task) runner$nutpie_closeout_task(task, manifest, "unused"), integer(1))
  expect_length(seeds, 18L)
  expect_equal(seeds, design$bridge_seed + 1000L * index)
  expect_false(anyDuplicated(seeds) > 0L)
})

test_that("conditional policy mixing retains scientific weights after attrition", {
  runner <- closeout_runner()
  complete <- runner$nutpie_closeout_weighted_summary(list(c(0, 0), c(10, 10)), c(.25, .75))
  failed <- runner$nutpie_closeout_weighted_summary(list(c(0, NA_real_), c(10, 10)), c(.25, .75))
  expect_equal(complete, failed)
  expect_equal(complete[["mean"]], 7.5)
  expect_equal(complete[["sd"]], sqrt(18.75))
  expect_equal(complete[["prob_gt_5"]], .75)
  expect_equal(complete[["median"]], 10)
  expect_error(runner$nutpie_closeout_weighted_summary(list(c(NA_real_), c(10, 10)), c(.25, .75)),
    "positive-weight k component")
})

test_that("collection records all unexecuted tasks without fitting or compiling", {
  runner <- closeout_runner()
  root <- tempfile("closeout-raw-"); report <- tempfile("closeout-report-")
  dir.create(root); dir.create(report)
  manifest <- list(root = root, design = runner$nutpie_closeout_design())
  status <- runner$nutpie_closeout_collect(manifest, report)
  expect_equal(nrow(status), 42L)
  expect_true(all(status$status == "not_run"))
  expect_length(list.files(root, all.files = FALSE), 0L)
  expect_true(file.exists(file.path(report, "nutpie_closeout_status.csv")))
  expect_true(file.exists(file.path(report, "nutpie_closeout_results.md")))
})

test_that("policy collection distinguishes exact identities from different constants", {
  runner <- closeout_runner()
  root <- tempfile("closeout-raw-"); report <- tempfile("closeout-report-")
  dir.create(root); dir.create(report); dir.create(file.path(root, "conditional_policy"))
  manifest <- list(root = root, design = runner$nutpie_closeout_design())
  mixture <- list()
  for (backend in c("cmdstanr", "nutpie", "mixed")) for (regime in c("merger", "tariff25")) {
    mixture[[paste(backend, regime, sep = "_")]] <- list(backend = backend, regime = regime,
      estimate = c(mean = if (backend == "mixed") 2 else 1, sd = 0),
      mcse = c(mean = NA_real_, sd = NA_real_), state_mcse = c(mean = 0, sd = 0),
      evidence_mcse = c(mean = 0, sd = 0), constant = c(mean = TRUE, sd = TRUE),
      failure_mass = 0, status = "all_states_valid")
  }
  saveRDS(list(status = "completed", reason = "", value = mixture,
    task = list(stage = "mixture")), file.path(root, "conditional_policy", "result.rds"))
  runner$nutpie_closeout_collect(manifest, report)
  checks <- read.csv(file.path(report, "nutpie_closeout_policy_comparisons.csv"))
  identity <- checks[checks$route == "nutpie", ]
  expect_true(all(identity$exact_constant_identity))
  expect_true(all(identity$agreement_passed))
  expect_true(all(identity$precision_passed))
  expect_true(all(identity$precision_status == "exempt_exact_identity"))
  expect_true(all(is.na(identity$combined_mcse)))
  different <- checks[checks$route == "mixed" & checks$metric == "mean", ]
  expect_true(all(!different$exact_constant_identity))
  expect_true(all(!different$agreement_passed))
  expect_true(all(!different$precision_passed))
})
