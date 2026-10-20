.validation_runner_env <- new.env(parent = globalenv())
.validation_runner_path <- testthat::test_path("..", "..", "doc", "validation", "nutpie_backend_validation.R")
testthat::skip_if_not(file.exists(.validation_runner_path),
  "Repository-only research runner is excluded from the source package.")
sys.source(.validation_runner_path, envir = .validation_runner_env)

test_that("validation workers cannot bypass the explicit inference opt-in", {
  withr::local_envvar(AB_NUTPIE_EXECUTE = "false")
  expect_error(.validation_runner_env$nutpie_validation_worker("nonexistent-task.rds"),
               "require AB_NUTPIE_EXECUTE=true")
})

test_that("worker failures are durably recorded before any inference", {
  directory <- tempfile("nutpie-worker-", tmpdir = tempdir())
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  task_path <- file.path(directory, "task.rds")
  saveRDS(list(root = directory, directory = directory, repository = ".",
    library_paths = .libPaths(), variant = "nutpie_diag",
    fixture = list(id = "B", conduct = "B", arguments = list(demand_model = "logit"))), task_path)
  withr::local_envvar(AB_NUTPIE_EXECUTE = "true")
  withr::local_options(list(antitrustBayes.cmdstan_cache = getOption("antitrustBayes.cmdstan_cache"),
                           antitrustBayes.nutpie_cache = getOption("antitrustBayes.nutpie_cache")))
  testthat::local_mocked_bindings(load_all = function(...) invisible(NULL), .package = "pkgload")
  testthat::local_mocked_bindings(.nutpie_model_for_kernel = function(...) stop("fixture compiler failure"),
                                  .package = "antitrustBayes")
  expect_message(result <- .validation_runner_env$nutpie_validation_worker(task_path), "fixture compiler failure")
  saved <- readRDS(file.path(directory, "status.rds"))
  expect_identical(result$status, "failed")
  expect_identical(saved$status, "failed")
  expect_identical(saved$reason, "fixture compiler failure")
  expect_false(file.exists(file.path(directory, "fit.rds")))
})

test_that("paired validation summaries detect shifts without normalizing away failed pairs", {
  reference <- data.frame(variable = c("alpha[1]", "rho"), mean = c(2, .3),
    median = c(2, .3), q025 = c(1, -.1), q975 = c(3, .7), sd = c(.5, .2),
    mcse_mean = .01, mcse_median = .01, mcse_q025 = .02, mcse_q975 = .02)
  alternative <- reference
  alternative$mean[1L] <- 2.01
  alternative$mean[2L] <- .8
  got <- .validation_runner_env$nutpie_validation_compare(reference, alternative, TRUE)
  expect_identical(got$status[got$metric == "mean"],
                   c("within_mcse_tolerance", "outside_mcse_tolerance"))
  expect_equal(got$combined_mcse[got$metric == "mean"], rep(sqrt(2) * .01, 2L))
  expect_equal(got$cmdstan_minus_nutpie[got$metric == "mean"], c(-.01, -.5))
  expect_true(all(.validation_runner_env$nutpie_validation_compare(reference, reference, FALSE)$status == "identity_mismatch"))
  alternative$mcse_q025[1L] <- NA_real_
  got <- .validation_runner_env$nutpie_validation_compare(reference, alternative, TRUE)
  expect_identical(got$status[got$metric == "q025" & got$variable == "alpha[1]"], "unavailable_mcse")
  reordered <- .validation_runner_env$nutpie_validation_compare(reference, reference[2:1, ], TRUE)
  expect_true(all(reordered$status == "within_mcse_tolerance"))
  different <- reference; different$variable[1L] <- "other_parameter"
  expect_error(.validation_runner_env$nutpie_validation_compare(reference, different, TRUE),
               "different validation parameter")
  reference$mcse_sd <- alternative$mcse_sd <- .01
  alternative$sd[2L] <- .4
  got <- .validation_runner_env$nutpie_validation_compare(reference, alternative, TRUE)
  expect_identical(got$status[got$metric == "sd"],
    c("within_mcse_tolerance", "outside_mcse_tolerance"))
})

test_that("validation fixtures preserve fixed conduct and active hierarchical routes", {
  skip_if_not_installed("coordination", "0.0.0.9004")
  fixed <- .validation_runner_env$nutpie_validation_fixture("BS")
  expect_identical(fixed$arguments$structure, "stackelberg")
  expect_identical(fixed$arguments$k_support, 2L)
  expect_equal(fixed$arguments$k_prior, 1)
  expect_true(all(fixed$equilibrium$passed))
  pooled <- .validation_runner_env$nutpie_validation_fixture("C", pooled = TRUE)
  expect_equal(length(unique(pooled$arguments$market)), 3L)
  expect_equal(length(unique(pooled$arguments$year)), 2L)
  expect_true(pooled$arguments$use_rho)
  expect_true(pooled$arguments$use_alpha_market_pool)
  expect_identical(pooled$arguments$s0_pooling, "market")
  expect_identical(pooled$arguments$quality_parameterization, "collapsed")
  expect_identical(pooled$arguments$market_intercept_parameterization, "collapsed")
  expect_identical(pooled$arguments$priors, bayes_prior_manifest())
  expect_true(all(pooled$equilibrium$passed))
})

test_that("validation collection retains every absent backend pair explicitly", {
  directory <- tempfile("nutpie-runner-", tmpdir = tempdir())
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  manifest <- list(run_id = "fixture-run", fixtures = list(C = list()),
    tasks = data.frame(fixture = rep("C", 3L), variant = c("cmdstanr", "nutpie_diag", "nutpie_low_rank"),
      id = c("C__cmdstanr", "C__nutpie_diag", "C__nutpie_low_rank")))
  .validation_runner_env$nutpie_validation_collect(manifest, directory, directory)
  status <- read.csv(file.path(directory, "nutpie_backend_status.csv"))
  comparisons <- read.csv(file.path(directory, "nutpie_backend_comparisons.csv"))
  expect_equal(nrow(status), 3L)
  expect_true(all(status$status == "not_run"))
  expect_equal(nrow(comparisons), 2L)
  expect_true(all(comparisons$status == "unavailable_pair"))
})

test_that("policy MCSE restores iteration and chain axes and withholds holes", {
  records <- data.frame(draw_index = c(5L, 2L, 1L, 6L, 3L, 4L),
                        price_pct = c(50, 20, 10, 60, 30, 40))
  got <- .validation_runner_env$nutpie_validation_policy_chains(records, c(3L, 2L))
  expect_equal(got, matrix(c(10, 20, 30, 40, 50, 60), 3L, 2L))
  expect_null(.validation_runner_env$nutpie_validation_policy_chains(records[-1L, ], c(3L, 2L)))
  records$price_pct[1L] <- NA_real_
  expect_null(.validation_runner_env$nutpie_validation_policy_chains(records, c(3L, 2L)))
  records$draw_index[1L] <- 2L
  expect_error(.validation_runner_env$nutpie_validation_policy_chains(records, c(3L, 2L)), "Duplicate posterior rows")
})

test_that("validation budgets are fresh invocation controls, with matched retained draws", {
  controls <- .validation_runner_env$nutpie_validation_controls()
  expect_identical(controls$chains, 4L)
  expect_identical(controls$iter - controls$warmup, 1000L)
  expect_identical(controls$cores, 2L)
  withr::local_envvar(AB_NUTPIE_MAX_SECONDS = "90")
  expect_equal(.validation_runner_env$nutpie_validation_budget(), 90)
  Sys.setenv(AB_NUTPIE_MAX_SECONDS = "yesterday")
  expect_error(.validation_runner_env$nutpie_validation_budget(), "integer >= 30")
})

test_that("postfit validation reuses identical saved draws without compiling or sampling", {
  directory <- tempfile("nutpie-postfit-", tmpdir = tempdir())
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  controls <- .validation_runner_env$nutpie_validation_controls()
  fixture <- list(id = "B", conduct = "B", arguments = list(demand_model = "logit"))
  args <- c(fixture$arguments, list(backend = "nutpie", stan_options = controls,
                                  nutpie_options = list(adaptation = "diag")))
  saveRDS(args, file.path(directory, "arguments.rds"))
  values <- array(c(rep(-1, 4000L), seq_len(4000L)), c(1000L, 4L, 2L),
    dimnames = list(NULL, NULL, c("lp__", "theta")))
  native <- structure(list(draws_array = values, metadata = list(sampling_time_native = 2,
    num_samples = 1000L, num_chains = 4L)), class = "NutpieFit")
  fit <- methods::new("BayesLogit", stanfit = native, stan_data = list(),
    metadata = list(stan_kernel = "logit_supply", model_spec_version = "test-fixture",
      stan_source_hash = digest::digest(file = .cmdstan_kernel_file("logit_supply"), algo = "sha256"),
      stan_data_hash = "fixture-data"))
  saveRDS(fit, file.path(directory, "fit.rds"))
  saveRDS(list(status = "postfit_failed", estimation_seconds = 3,
    sampling_seconds = NA_real_, total_seconds = 5), file.path(directory, "status.rds"))
  task_path <- file.path(directory, "task.rds")
  saveRDS(list(root = directory, directory = directory, repository = ".", fixture = fixture,
    variant = "nutpie_diag", library_paths = .libPaths(), sampler_controls = controls), task_path)
  withr::local_envvar(AB_NUTPIE_EXECUTE = "true")
  withr::local_options(list(antitrustBayes.cmdstan_cache = getOption("antitrustBayes.cmdstan_cache"),
                           antitrustBayes.nutpie_cache = getOption("antitrustBayes.nutpie_cache")))
  testthat::local_mocked_bindings(load_all = function(...) invisible(NULL), .package = "pkgload")
  testthat::local_mocked_bindings(
    .nutpie_model_for_kernel = function(...) stop("Unexpected compilation"),
    bayes_est = function(...) stop("Unexpected resampling"), .package = "antitrustBayes")
  worker <- .validation_runner_env$nutpie_validation_worker
  sandbox <- new.env(parent = environment(worker))
  sandbox$nutpie_validation_parameters <- function(...) list(summary = data.frame(ess_bulk = 10))
  sandbox$nutpie_validation_diagnostics <- function(...) list(summary = data.frame(min_bulk_ess = 10))
  sandbox$nutpie_validation_policy <- function(...) list(summary = data.frame(mean = 0))
  environment(worker) <- sandbox
  result <- worker(task_path)
  expect_identical(result$status, "completed")
  expect_true(result$fit_reused)
  expect_identical(readRDS(file.path(directory, "fit.rds")), fit)
  saved <- readRDS(file.path(directory, "status.rds"))
  expect_true(saved$total_seconds >= 5)
  expect_identical(saved$reason, "")
  expect_true(file.exists(file.path(directory, "validation.rds")))
})

test_that("runtime seal detects source changes before validation can reuse scientific output", {
  directory <- tempfile("nutpie-runtime-", tmpdir = tempdir())
  dir.create(file.path(directory, "R"), recursive = TRUE)
  dir.create(file.path(directory, "doc", "validation"), recursive = TRUE)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  runner <- file.path(directory, "doc", "validation", "nutpie_backend_validation.R")
  writeLines("# test runner", runner)
  source <- file.path(directory, "R", "policy.R")
  writeLines("# first policy", source)
  first <- .validation_runner_env$nutpie_validation_runtime(directory)
  expect_identical(first, .validation_runner_env$nutpie_validation_runtime(directory))
  writeLines("# changed policy", source)
  second <- .validation_runner_env$nutpie_validation_runtime(directory)
  expect_false(identical(first, second))
  expect_identical(first$packages, second$packages)
  expect_false(identical(first$source_hashes[["R/policy.R"]], second$source_hashes[["R/policy.R"]]))
})

test_that("CES Stackelberg validation runs merger and tariff policies with inherited leaders", {
  skip_if_not_installed("trade")
  fixture <- .validation_runner_env$nutpie_validation_ces_fixture("CS")
  a <- fixture$arguments
  market <- bayes_market(prices = a$prices, quantities = a$quantities,
    margins = a$margins, group = a$group, outside_price = 1,
    insideSize = 1000, market_side = "output")
  specified <- bayes_specify(bayes_params(gamma = 6, s0 = .2,
    supply_model = 3L, structure = "stackelberg", partition_prior = "ordered_k",
    k_support = 2L, k_prior = 1), market)
  plan <- posterior_policy_draws(specified, n_draws = 1L, seed = 491871L)
  prepared <- suppressWarnings(posterior_policy_models(specified, draw_plan = plan,
    insideSize = 1000, demand_anchoring = "latent_quality"))
  expect_s4_class(prepared$models[[1L]][[1L]], "StackelbergCES")
  fit <- methods::new("BayesLogit", stanfit = NULL,
    data = as.data.frame(a[c("prices", "quantities", "margins", "market", "year", "group")]))
  ## A specified economic state replaces sampling only. Run the complete
  ## validation policy function, including native promotion and both solves.
  local_mocked_bindings(
    .bayes_backend_draws = function(...) array(0, c(1L, 1L, 1L)),
    posterior_policy_draws = function(...) plan,
    posterior_policy_models = function(...) prepared,
    .package = "antitrustBayes")
  native_simulate <- trade::simulate
  solved <- list()
  local_mocked_bindings(simulate = function(...) {
    value <- native_simulate(...)
    solved[[length(solved) + 1L]] <<- value
    value
  }, .package = "trade")
  result <- suppressWarnings(.validation_runner_env$nutpie_validation_policy(fit))
  expect_equal(nrow(result$draws), 2L)
  expect_setequal(result$draws$regime, c("merger", "tariff25"))
  expect_true(all(result$draws$status == "ok"), info = paste(result$draws$reason, collapse = "; "))
  expect_true(all(is.finite(result$draws$price_pct)))
  expect_true(all(result$draws$cost_identity_error < 1e-8))
  expect_true(all(result$summary$n_valid == 1L))
  expect_identical(result$plan, plan)
  expect_length(solved, 2L)
  if (length(solved) == 2L) {
    firms <- aggregate(a$quantities, list(firm = a$group), sum)
    pair <- as.character(firms$firm[order(-firms$x)[1:2]])
    inherited <- as.character(solved[[1L]]@tariff_state$leadersPre)
    inherited[inherited == pair[[2L]]] <- pair[[1L]]
    expect_identical(as.character(solved[[1L]]@tariff_state$leadersPost), unique(inherited))
    expect_identical(as.character(solved[[2L]]@tariff_state$leadersPost),
      as.character(solved[[2L]]@tariff_state$leadersPre))
  }
})
