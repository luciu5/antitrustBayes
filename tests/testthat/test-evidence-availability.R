.availability_fit <- function(backend) {
  data <- data.frame(prices = c(1, 2), quantities = c(4, 6), margins = c(.25, .2),
    market = "m", year = "1", group = c("a", "b"), loo_cluster = "case",
    margin_observed = TRUE, outside_price = 1)
  outcomes <- list(shareIn = c(.4, .6), clr_shareIn = log(c(.4, .6)),
    marginNormInv = c(4, 5), log_marginNormInv = log(c(4, 5)), margin_obs = c(1L, 1L),
    market_year_idx = c(1L, 1L), cluster = c(1L, 1L), likelihood_mode = 1L)
  native <- if (backend == "nutpie") structure(list(metadata = list(
    provenance = list(nutpieR = "fixture"), resolved_options = list(adaptation = "diag"))),
    class = "NutpieFit") else NULL
  new("BayesLogit", stanfit = native, data = data, stan_data = outcomes,
    supply_model = 1L, demand_model = "logit", market_side = "output", price_sd = 1,
    metadata = list(backend = backend, inference = "hmc", quality_parameterization = "collapsed",
      prior_manifest = "baseline", resolved_priors = bayes_prior_manifest("baseline"),
      sampler_controls = list(adaptation = if (backend == "nutpie") "diag" else "diag_e")))
}

.availability_mock_target <- function(env) {
  # These tests start after scientific/target validation and inject numerical
  # bridge outcomes. Native integration validates the actual target separately.
  testthat::local_mocked_bindings(
    validate_fit_spec = function(...) invisible(TRUE),
    .bayes_validate_bridge = function(...) invisible(TRUE),
    .bayes_evidence_contract = function(...) list(kernel = "fixture"),
    .bayes_bridge_restore_context = function(...) list(samples = matrix(0, 20, 1),
      bridge_args = list(), source = "fixture", kernel = "fixture", saved_hash = "same",
      current_hash = "same", migration = FALSE,
      cleanup = function() { env$cleaned <- env$cleaned + 1L }),
    .bayes_backend_diagnostic_summary = function(...) list(chains = data.frame(divergences = 20)),
    .package = "antitrustBayes", .env = parent.frame())
}

test_that("public mixed evidence retains partial attempts without a sampler permission gate", {
  state <- new.env(parent = emptyenv()); state$attempt <- state$cleaned <- 0L
  .availability_mock_target(state)
  testthat::local_mocked_bindings(bridge_sampler = function(...) {
    state$attempt <- state$attempt + 1L
    if (state$attempt == 2L) stop("retained numerical failure")
    if (state$attempt == 4L) warning("retained finite-result warning")
    list(logml = c(-2, NA, -2.2, -3, -3.2, -3.4)[state$attempt], mcse_logml = .02)
  }, .package = "bridgesampling")
  expect_warning(expect_warning(comp <- bayes_compare(
    native = .availability_fit("nutpie"), cmd = .availability_fit("cmdstanr"),
    prior_weights = c(.4, .6), bridge_repetitions = 3L), "retained numerical failure"),
    "retained finite-result warning")
  expect_equal(state$attempt, 6L); expect_equal(state$cleaned, 2L)
  expect_equal(comp$logml, c(-2.1, -3.2))
  expect_equal(comp$logml_sd, c(sd(c(-2, -2.2)), .2))
  expect_identical(comp$numerical_status, c("partial", "complete"))
  expect_equal(comp$completed_repetitions, c(2, 3))
  expect_equal(comp$requested_repetitions, c(3, 3))
  info <- attr(comp, "bridge_context_info")
  expect_identical(info$native$repetitions$status, c("available", "failed", "available"))
  expect_match(info$native$repetitions$error[2], "retained numerical failure")
  expect_match(info$cmd$repetitions$warnings[1], "retained finite-result warning")
  expect_equal(info$native$repetitions$mcse_logml, c(.02, NA, .02))
  expect_true(is.null(attr(comp, "bridge_list")[[1]][[2]]))
  expect_true(all(is.finite(bayes_weights(comp, "bridge"))))
  expect_equal(sum(comp$post_prob), 1)
  metadata <- attr(comp, "backend_metadata")
  expect_identical(metadata$native$adaptation, "diag")
  expect_equal(metadata$native$diagnostics$chains$divergences, 20)
})

test_that("unavailable positive-prior evidence withholds weights for the complete set", {
  state <- new.env(parent = emptyenv()); state$attempt <- state$cleaned <- 0L
  .availability_mock_target(state)
  testthat::local_mocked_bindings(bridge_sampler = function(...) {
    state$attempt <- state$attempt + 1L
    list(logml = if (state$attempt <= 3L) -2 else Inf)
  }, .package = "bridgesampling")
  comp <- bayes_compare(a = .availability_fit("nutpie"), b = .availability_fit("cmdstanr"))
  expect_equal(comp$logml[1], -2)
  expect_true(is.na(comp$logml[2])); expect_true(all(is.na(comp$post_prob)))
  expect_identical(comp$numerical_status, c("complete", "unavailable"))
  expect_true(all(attr(comp, "bridge_context_info")$b$repetitions$status == "failed"))
  expect_true(all(is.infinite(attr(comp, "bridge_context_info")$b$repetitions$logml)))
  expect_error(bayes_weights(comp, "bridge"), "Every positive-prior")
  state$attempt <- 0L
  zero_prior <- bayes_compare(a = .availability_fit("nutpie"), b = .availability_fit("cmdstanr"),
    prior_weights = c(1, 0))
  expect_equal(zero_prior$post_prob, c(1, 0))
  expect_equal(bayes_weights(zero_prior, "bridge"), c(a = 1, b = 0))
  expect_identical(zero_prior$evidence_status, c("available", "unavailable"))
})

test_that("zero-prior evidence does not suppress required candidates or cause numerical underflow", {
  state <- new.env(parent = emptyenv()); state$attempt <- state$cleaned <- 0L
  .availability_mock_target(state)
  testthat::local_mocked_bindings(bridge_sampler = function(...) {
    state$attempt <- state$attempt + 1L
    list(logml = if (state$attempt <= 3L) -10000 else 10000)
  }, .package = "bridgesampling")
  comp <- bayes_compare(a = .availability_fit("nutpie"), b = .availability_fit("cmdstanr"),
    prior_weights = c(1, 0))
  expect_equal(comp$post_prob, c(1, 0))
  expect_equal(bayes_weights(comp, "bridge"), c(a = 1, b = 0))
})
