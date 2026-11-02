.policy_models_fixture <- function(demand = "logit", conduct = 1L,
                                   structure = "simultaneous", n = 3L) {
  args <- list(s0 = .2, supply_model = conduct, structure = structure)
  if (demand == "logit") args$alpha <- rep(6, n) else args$gamma <- rep(10, n)
  if (structure != "simultaneous") {
    args$partition_prior <- "ordered_k"
    args$k_support <- if (structure == "fringe") 0:2 else 1:2
    args$k_prior <- rep(1 / length(args$k_support), length(args$k_support))
  }
  params <- do.call(bayes_params, args)
  market <- bayes_market(
    prices = c(1.1, 1.2, 1.3), quantities = c(40, 30, 20),
    margins = c(.5, .4, .3), group = c("a", "a", "b"),
    outside_price = 1, market_side = "output", insideSize = 90,
    labels = c("one", "two", "three"))
  bayes_specify(params, market)
}

test_that("finite negative fringe costs retain draw identity and a typed failure", {
  params <- bayes_params(
    s0 = .2, supply_model = 3L, structure = "fringe", alpha = rep(1, 3),
    partition_prior = "ordered_k", k_support = 1L, k_prior = 1)
  market <- bayes_market(
    prices = c(1.1, 1.2, 1.3), quantities = c(40, 30, 20),
    margins = c(.5, .4, .3), group = c("a", "a", "b"),
    outside_price = 1, market_side = "output", insideSize = 90,
    labels = c("one", "two", "three"))
  fit <- bayes_specify(params, market)
  plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 11L)
  realized <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan))
  expect_identical(realized$draw_indices, plan$draw_indices)
  expect_true(all(vapply(realized$models[[1L]], is.null, logical(1))))
  expect_identical(realized$failures$draw_index, plan$draw_indices)
  expect_true(all(realized$failures$category == "bayes_finite_nonpositive_cost"))
})

test_that("finite negative baseline costs are typed for every candidate model", {
  candidates <- c(MC = 2L, B = 1L, C = 3L, A = 4L,
                  BF = 1L, CF = 3L, BS = 1L, CS = 3L)
  for (candidate in names(candidates)) {
    structure <- if (candidate %in% c("BF", "CF")) "fringe" else
      if (candidate %in% c("BS", "CS")) "stackelberg" else "simultaneous"
    args <- list(s0 = .2, supply_model = candidates[[candidate]],
                 structure = structure, alpha = rep(.5, 3L))
    if (structure != "simultaneous") {
      args$partition_prior <- "ordered_k"
      args$k_support <- 1L
      args$k_prior <- 1
    }
    fit <- bayes_specify(do.call(bayes_params, args), bayes_market(
      prices = c(1.1, 1.2, 1.3), quantities = c(40, 30, 20),
      margins = c(.5, .4, .3), group = c("a", "b", "c"),
      outside_price = 1, market_side = "output", insideSize = 90,
      labels = c("one", "two", "three")))
    plan <- posterior_policy_draws(fit, n_draws = 1L, seed = 11L)
    realized <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan))
    expect_null(realized$models[[1L]][[1L]], info = candidate)
    expect_identical(realized$failures$draw_index, plan$draw_indices,
                     info = candidate)
    expect_identical(realized$failures$category,
                     "bayes_finite_nonpositive_cost", info = candidate)
  }
})

test_that("public realization preserves canonical baseline state without post simulation", {
  testthat::local_mocked_bindings(
    .simulate_antitrust_fit = function(...) stop("unexpected post simulation"),
    .package = "antitrust")
  for (demand in c("logit", "ces")) for (conduct in c(1L, 2L, 3L)) {
    fit <- .policy_models_fixture(demand, conduct)
    before <- serialize(fit, NULL)
    plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 11L)
    prepared <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan))
    expect_s3_class(prepared, "bayes_policy_models")
    expect_identical(prepared$draw_plan, plan)
    expect_identical(prepared$draw_indices, plan$draw_indices)
    expect_equal(nrow(prepared$failures), 0L)
    for (pos in seq_along(plan$draw_indices)) {
      model <- prepared$models$unit_1[[pos]]
      expect_s4_class(model, "AntitrustFit")
      expect_equal(model@model@pricePre, fit@market$prices, ignore_attr = TRUE)
      expect_equal(unname(model@model@pricePost), unname(model@model@pricePre))
      expect_identical(attr(model, "bayes_policy_draw")$draw_index, plan$draw_indices[[pos]])
    }
    expect_identical(serialize(fit, NULL), before)
  }
})

test_that("native realizations preserve planned k and demand across anchoring modes", {
  skip_if_not_installed("coordination", "0.0.0.9004")
  for (demand in c("logit", "ces")) for (conduct in c(1L, 3L)) {
    for (structure in c("fringe", "stackelberg")) {
      fit <- .policy_models_fixture(demand, conduct, structure)
      plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 18L)
      observed <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan))
      latent <- suppressWarnings(posterior_policy_models(
        fit, draw_plan = plan, demand_anchoring = "latent_quality"))
      expect_equal(nrow(observed$failures), 0L)
      expect_equal(nrow(latent$failures), 0L)
      for (pos in seq_along(plan$draw_indices)) {
        one <- observed$models$unit_1[[pos]]
        two <- latent$models$unit_1[[pos]]
        expect_s4_class(one, paste0(if (structure == "fringe") "CoreFringe" else "Stackelberg",
                                    if (demand == "logit") "Logit" else "CES"))
        expect_equal(attr(one, "bayes_policy_draw")$k, plan$k_cache[[pos]][["unit_1"]])
        expect_equal(attr(two, "bayes_policy_draw")$k, plan$k_cache[[pos]][["unit_1"]])
        expect_equal(attr(two, "bayes_policy_draw")$realized_k,
                     plan$k_cache[[pos]][["unit_1"]])
        expect_equal(one@mcPre, two@mcPre)
        expect_equal(unname(one@pricePre), unname(two@pricePre), tolerance = 1e-9)
        expect_equal(antitrust::calcShares(one, TRUE, revenue = FALSE),
                     antitrust::calcShares(two, TRUE, revenue = FALSE), tolerance = 1e-9)
        expect_identical(attr(two, "bayes_policy_draw")$plan_fingerprint,
                         plan$fit_signature$latent_state_fingerprint)
      }
    }
  }
})

test_that("specified realization handles missing margins and absent structure metadata", {
  fit <- .policy_models_fixture()
  fit@metadata$structure <- NULL
  fit@market$margins[[2L]] <- NA_real_
  fit@market$margin_observed[[2L]] <- FALSE
  plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 10L)
  prepared <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan))
  expect_equal(nrow(prepared$failures), 0L)
  expect_identical(prepared$draw_plan, plan)
  expect_identical(prepared$draw_indices, plan$draw_indices)
  expect_true(all(vapply(prepared$models$unit_1, methods::is, logical(1), "AntitrustFit")))
})

test_that("default product labels remain usable for multiproduct firms", {
  fit <- .policy_models_fixture(structure = "stackelberg")
  fit@market$labels <- fit@market$group
  plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 10L)
  prepared <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan))
  expect_equal(nrow(prepared$failures), 0L)
  for (model in prepared$models$unit_1) {
    expect_equal(anyDuplicated(model@labels), 0L)
    expect_equal(model@firmOwnerPre, fit@market$group)
    expect_equal(unname(model@pricePre), fit@market$prices, tolerance = 1e-9)
  }
})

test_that("realization failures keep draw positions and posterior identities", {
  fit <- .policy_models_fixture(n = 3L)
  plan <- posterior_policy_draws(fit, n_draws = 3L, seed = 8L)
  original <- antitrustBayes:::build_conduct_object
  bad_draw <- plan$draw_indices[[2L]]
  testthat::local_mocked_bindings(build_conduct_object = function(params, draw_idx, ...) {
    if (draw_idx == bad_draw) stop("fixture construction failure")
    original(params, draw_idx, ...)
  }, .package = "antitrustBayes")
  expect_warning(prepared <- posterior_policy_models(fit, draw_plan = plan),
                 "1 structural draw realizations failed")
  expect_length(prepared$models$unit_1, 3L)
  expect_null(prepared$models$unit_1[[2L]])
  expect_s4_class(prepared$models$unit_1[[1L]], "AntitrustFit")
  expect_s4_class(prepared$models$unit_1[[3L]], "AntitrustFit")
  expect_equal(prepared$failures$draw_index, bad_draw)
  expect_equal(prepared$failures$stage, "baseline_construction")
  expect_equal(names(prepared$models$unit_1), as.character(plan$draw_indices))
})

test_that("realization rejects incompatible plans, selectors and models", {
  fit <- .policy_models_fixture()
  plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 3L)
  expect_error(posterior_policy_models(fit, draw_plan = plan, seed = 3), "new draw_plan")
  expect_error(posterior_policy_models(fit, draw_plan = plan, n_draws = 1L), "conflicts")
  expect_error(posterior_policy_models(fit, draw_plan = plan, market = "absent"), "markets")
  changed <- plan
  changed$quality_draws[1L, 1L] <- 1
  expect_error(posterior_policy_models(fit, draw_plan = changed), "latent states")
  fit@market_side <- "input"
  expect_error(posterior_policy_models(fit), "output Logit/CES")
})
