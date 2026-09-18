.latent_policy_fixture <- function(demand = "logit", supply = 1L,
                                   structure = NULL,
                                   partition_prior = "cutoff",
                                   k_support = NULL, k_prior = NULL,
                                   n_draws = 1L, market_side = "output",
                                   group = c("a", "b", "c")) {
  prices <- c(1.1, 1.2, 1.3)
  quantities <- c(40, 30, 20)
  margins <- c(.50, .40, .30)
  if (demand == "logit") {
    params <- bayes_params(
      alpha = rep(6, n_draws), s0 = .2, supply_model = supply,
      structure = structure, partition_prior = partition_prior,
      k_support = k_support, k_prior = k_prior
    )
  } else {
    params <- bayes_params(
      gamma = rep(10, n_draws), s0 = .2, supply_model = supply,
      structure = structure, partition_prior = partition_prior,
      k_support = k_support, k_prior = k_prior
    )
  }
  market <- bayes_market(
    prices = prices, quantities = quantities, margins = margins,
    group = group, outside_price = 1, market_side = market_side,
    insideSize = 90, labels = paste0("p", seq_along(prices))
  )
  list(
    fit = bayes_specify(params, market), market = market,
    owner = {
      owner <- group
      firms <- unique(as.character(group))
      if (length(firms) > 1L) owner[group == firms[[2L]]] <- firms[[1L]]
      owner
    }, params = params
  )
}

test_that("posterior policy plans are reusable and strict about conflicts", {
  fx <- .latent_policy_fixture(n_draws = 4L)
  plan <- posterior_policy_draws(fx$fit, n_draws = 3L, seed = 12L)
  expect_s3_class(plan, "bayes_policy_draw_plan")
  expect_equal(length(plan$draw_indices), 3L)
  expect_equal(colnames(plan$quality_draws), c("a", "b", "c"))

  one <- suppressWarnings(simulate(
    fx$fit, fx$owner, draw_plan = plan, mc.cores = 1L,
    demand_anchoring = "latent_quality"
  ))
  two <- suppressWarnings(bayes_sim(
    fx$fit, fx$owner, draw_plan = plan, mc.cores = 1L,
    demand_anchoring = "latent_quality"
  ))
  expect_equal(one@metadata$draw_indices, plan$draw_indices)
  expect_equal(two@metadata$draw_indices, plan$draw_indices)
  expect_equal(one@metadata$demand_anchoring, "latent_quality")
  expect_equal(one@metadata$draw_plan$draw_indices, plan$draw_indices)
  expect_error(
    simulate(fx$fit, fx$owner, n_draws = 2L, draw_plan = plan,
             mc.cores = 1L),
    "conflicts"
  )
  expect_error(
    simulate(fx$fit, fx$owner, draw_plan = plan, role_mapping = "map",
             mc.cores = 1L),
    "conflicts"
  )
  changed <- fx$fit
  changed@metadata$latent_test_marker <- TRUE
  expect_error(
    simulate(changed, fx$owner, draw_plan = plan, mc.cores = 1L),
    "different fit"
  )
})

test_that("draw plans bind quality and role states to the selected posterior rows", {
  fx <- .latent_policy_fixture(n_draws = 4L)
  plan <- posterior_policy_draws(fx$fit, n_draws = 3L, seed = 12L)
  expect_identical(rownames(plan$quality_draws), as.character(plan$draw_indices))
  changed <- plan
  changed$draw_indices <- rev(changed$draw_indices)
  expect_error(.lq_validate_plan(fx$fit, changed), "posterior rows were changed")
  changed$metadata$draw_indices <- changed$draw_indices
  rownames(changed$quality_draws) <- as.character(changed$draw_indices)
  expect_error(.lq_validate_plan(fx$fit, changed), "posterior rows were changed")
  changed <- plan
  changed$quality_draws[1L, ] <- c(.01, -.01, 0)
  expect_error(.lq_validate_plan(fx$fit, changed), "latent states")
  changed <- plan
  changed$k_cache <- list(changed = TRUE)
  expect_error(.lq_validate_plan(fx$fit, changed), "latent states")
})

test_that("sampled and collapsed quality maps retain fitted group labels", {
  make_fit <- function(mode) {
    stan_data <- list(
      N = 4L, N_group = 2L, group = c(1L, 2L, 1L, 2L),
      log_shareIn = c(.2, .3, .1, .4),
      quality_parameterization = if (mode == "sampled") 1L else 2L,
      share_likelihood = 2L, likelihood_mode = 1L
    )
    attr(stan_data, "group_levels") <- c("z", "a")
    methods::new(
      "BayesLogit", stanfit = NULL, stan_data = stan_data,
      data = data.frame(group = c("z", "a", "z", "a")),
      supply_model = 1L, demand_model = "logit", market_side = "output",
      metadata = list(
        quality_parameterization = mode, posterior_draws = 2L,
        share_likelihood = "gaussian"
      ), price_sd = 1
    )
  }
  sampled <- make_fit("sampled")
  draws <- list(
    a_market = matrix(c(0, 0), 2L, 1L),
    b_quality_all = matrix(c(.2, -.2, .3, -.3), 2L, 2L, byrow = TRUE),
    share_mean = matrix(0, 2L, 4L), sigma_b_strat = rep(.5, 2L),
    sigma_share = rep(.4, 2L)
  )
  got <- .lq_quality_draws(sampled, draws, c(2L, 1L))
  expect_equal(colnames(got$draws), c("z", "a"))
  expect_equal(rownames(got$draws), c("2", "1"))
  expect_equal(unname(got$draws[, "z"]), c(.3, .2))
  expect_equal(unname(got$draws[, "a"]), c(-.3, -.2))

  collapsed <- make_fit("collapsed")
  got_collapsed <- .lq_quality_draws(collapsed, draws, c(2L, 1L), seed = 3L)
  expect_equal(colnames(got_collapsed$draws), c("z", "a"))
  expect_true(all(abs(rowSums(got_collapsed$draws)) < 1e-8))

  collapsed@stan_data$use_rho <- 1L
  collapsed@stan_data$likelihood_mode <- 1L
  collapsed@stan_data$margin_idx <- 1L
  collapsed@metadata$partition_prior <- "ordered_k"
  expect_error(
    .lq_quality_draws(collapsed, draws, c(1L, 2L), seed = 3L),
    "component-specific"
  )
})

test_that("latent demand replacement keeps the observed cost primitive and solves pre FOCs", {
  for (demand in c("logit", "ces")) {
    fx <- .latent_policy_fixture(demand = demand)
    baseline <- build_conduct_object(
      fx$params, 1L, fx$market, ownerPre = fx$market$group,
      ownerPost = fx$market$group
    )
    antitrust_fit <- suppressWarnings(as_antitrust_fit(baseline))
    model <- antitrust_fit@model
    old_mc <- model@mcPre
    quality_score <- log(c(.60, .25, .15))
    info <- .lq_replace_demand(
      model, quality_score, s0 = .2, prices_reference = fx$market$prices,
      price_outside = 1, quality = c(.3, -.1, -.2), draw_index = 1L
    )
    expect_equal(info$mc_pre, unname(old_mc))
    expect_equal(info$reference_price, fx$market$prices)
    expect_true(max(abs(info$foc_solved)) < 1e-6)
    expect_equal(sum(info$reference_quantity_shares), .8, tolerance = 1e-10)
    expect_true(any(abs(info$price_discrepancy) > 1e-8))
    if (demand == "ces") {
      expect_equal(info$model@mktSize, model@mktSize)
    }
  }
  fx <- .latent_policy_fixture()
  baseline <- build_conduct_object(
    fx$params, 1L, fx$market, ownerPre = fx$market$group,
    ownerPost = fx$market$group
  )
  baseline@mcPre[] <- NA_real_
  expect_error(
    .lq_replace_demand(baseline, log(c(.5, .3, .2)), .2,
                       fx$market$prices, 1, rep(0, 3), 1L),
    "numeric mcPre"
  )
})

test_that("latent mode supports output Logit and CES policy paths", {
  for (demand in c("logit", "ces")) {
    fx <- .latent_policy_fixture(demand = demand)
    result <- suppressWarnings(simulate(
      fx$fit, fx$owner, mc.cores = 1L, demand_anchoring = "latent_quality"
    ))
    expect_equal(as.integer(result@n_valid), 1L)
    expect_identical(result@metadata$demand_anchoring, "latent_quality")
    expect_true(is.list(result@metadata$latent_diagnostics$unit_1[[1L]]))
    expect_true(result@metadata$latent_diagnostics$unit_1[[1L]]$
                  max_abs_foc_at_latent_pre < 1e-6)
  }
  fx <- .latent_policy_fixture(n_draws = 2L)
  plan <- posterior_policy_draws(fx$fit, n_draws = 2L, seed = 17L)
  merger_cost <- antitrust::counterfactual(
    ownership = fx$owner, costs = rep(.10, fx$market$n_prod)
  )
  result <- suppressWarnings(simulate(
    fx$fit, merger_cost, draw_plan = plan, mc.cores = 1L,
    demand_anchoring = "latent_quality"
  ))
  expect_equal(as.integer(result@n_valid), 2L)
  expect_equal(result@metadata$draw_indices, plan$draw_indices)
})

test_that("latent input Logit preserves value primitives across B/C policies", {
  for (supply in c(1L, 3L)) {
    fx <- .latent_policy_fixture(supply = supply, market_side = "input")
    baseline <- build_conduct_object(
      fx$params, 1L, fx$market, ownerPre = fx$market$group,
      ownerPost = fx$market$group
    )
    antitrust_fit <- suppressWarnings(as_antitrust_fit(baseline))
    value_pre <- unname(antitrust_fit@model@mcPre)
    result <- suppressWarnings(simulate(
      fx$fit, fx$owner, mc.cores = 1L, demand_anchoring = "latent_quality"
    ))
    expect_equal(as.integer(result@n_valid), 1L)
    simulated <- result@sims$unit_1[[1L]]
    expect_false(isTRUE(simulated@output))
    expect_equal(unname(simulated@mcPre), value_pre, tolerance = 1e-10)
    expect_equal(unname(simulated@mcPost), value_pre, tolerance = 1e-10)
    expect_lt(max(abs(unname(simulated@mcPre) -
                      unname(simulated@pricePre) -
                      antitrust::calcMargins(simulated,
                                             preMerger = TRUE,
                                             level = TRUE))), 1e-6)
    expect_lt(max(abs(unname(simulated@mcPost) -
                      unname(simulated@pricePost) -
                      antitrust::calcMargins(simulated,
                                             preMerger = FALSE,
                                             level = TRUE))), 1e-6)
    diag <- result@metadata$latent_diagnostics$unit_1[[1L]]
    expect_equal(sum(diag$reference_quantity_shares), .8, tolerance = 1e-10)
    expect_equal(sum(diag$latent_pre_quantity_shares), .8, tolerance = 1e-10)
    expect_identical(diag$foc_definition,
                     "mc - price - positive dollar markdown")
    expect_lt(diag$max_abs_foc_at_latent_pre, 1e-6)
  }
})

test_that("nonzero input quality moves the baseline without changing values", {
  for (supply in c(1L, 3L)) {
    fx <- .latent_policy_fixture(supply = supply, market_side = "input")
    baseline <- build_conduct_object(
      fx$params, 1L, fx$market, ownerPre = fx$market$group,
      ownerPost = fx$market$group
    )
    model <- suppressWarnings(as_antitrust_fit(baseline))@model
    original_values <- unname(model@mcPre)
    changed <- .lq_replace_demand(
      model, log(c(.60, .25, .15)), s0 = .2,
      prices_reference = fx$market$prices, price_outside = 1,
      quality = c(.3, -.1, -.2), draw_index = 1L
    )
    expect_equal(unname(changed$model@mcPre), original_values)
    expect_true(any(abs(changed$price_discrepancy) > 1e-8))
    expect_equal(sum(changed$reference_quantity_shares), .8,
                 tolerance = 1e-10)
    markdown <- antitrust::calcMargins(changed$model, preMerger = TRUE,
                                       level = TRUE)
    expect_true(all(is.finite(markdown) & markdown > 0))
    expect_lt(max(abs(original_values - changed$model@pricePre - markdown)),
              1e-6)
  }
})

test_that("latent input Logit supports native Stackelberg and core-fringe paths", {
  skip_if_not_installed("coordination")
  for (structure in c("stackelberg", "fringe")) {
    fx <- .latent_policy_fixture(
      supply = 1L, structure = structure, partition_prior = "ordered_k",
      k_support = 1:2, k_prior = c(.5, .5), market_side = "input"
    )
    plan <- posterior_policy_draws(fx$fit, seed = 28L)
    result <- suppressWarnings(simulate(
      fx$fit, fx$owner, draw_plan = plan, mc.cores = 1L,
      demand_anchoring = "latent_quality"
    ))
    expect_equal(as.integer(result@n_valid), 1L)
    simulated <- result@sims$unit_1[[1L]]
    expect_false(isTRUE(simulated@output))
    expect_equal(unname(simulated@mcPost), unname(simulated@mcPre),
                 tolerance = 1e-10)
    residual <- unname(simulated@mcPost) - unname(simulated@pricePost) -
      antitrust::calcMargins(simulated, preMerger = FALSE, level = TRUE)
    expect_true(all(is.finite(residual)))
    expect_lt(max(abs(residual)), 1e-6)
    diag <- result@metadata$latent_diagnostics$unit_1[[1L]]
    expect_identical(diag$foc_definition,
                     "mc - price - positive dollar markdown")
    expect_identical(as.integer(diag$k),
                     as.integer(plan$k_cache[[as.character(plan$draw_indices[[1L]])]]
                                      [["unit_1"]]))
  }
})

test_that("a repeated firm quality is kept fixed through a multiproduct policy", {
  fx <- .latent_policy_fixture(group = c("a", "a", "b"))
  plan <- posterior_policy_draws(fx$fit, seed = 4L)
  result <- suppressWarnings(simulate(
    fx$fit, fx$owner, draw_plan = plan, mc.cores = 1L,
    demand_anchoring = "latent_quality"
  ))
  expect_equal(as.integer(result@n_valid), 1L)
  diag <- result@metadata$latent_diagnostics$unit_1[[1L]]
  expect_equal(length(diag$quality), 2L)
  expect_equal(sum(diag$quality), 0, tolerance = 1e-10)
})

test_that("exact-k plans share one role draw across native policy simulations", {
  skip_if_not_installed("coordination")
  for (demand in c("logit", "ces")) {
    for (supply in c(1L, 3L)) {
      for (structure in c("stackelberg", "fringe")) {
        fx <- .latent_policy_fixture(
          demand = demand, supply = supply, structure = structure,
          partition_prior = "ordered_k", k_support = 1:2,
          k_prior = c(.5, .5)
        )
        plan <- posterior_policy_draws(fx$fit, seed = 8L)
        result <- suppressWarnings(simulate(
          fx$fit, fx$owner, draw_plan = plan, mc.cores = 1L,
          demand_anchoring = "latent_quality"
        ))
        expect_equal(as.integer(result@n_valid), 1L)
        expect_true(isTRUE(result@metadata$draw_plan$exact_k))
        expect_equal(result@metadata$latent_diagnostics$unit_1[[1L]]$k,
                     plan$k_cache[[as.character(plan$draw_indices[[1L]])]][["unit_1"]])
      }
    }
  }
})

test_that("latent mode rejects unsupported normalization and conduct combinations", {
  input <- .latent_policy_fixture(market_side = "input")
  expect_silent(suppressWarnings(simulate(
    input$fit, input$owner, demand_anchoring = "latent_quality", mc.cores = 1L
  )))
  input_ces <- .latent_policy_fixture(demand = "ces", market_side = "input")
  expect_error(
    simulate(input_ces$fit, input_ces$owner,
             demand_anchoring = "latent_quality", mc.cores = 1L), "input Logit"
  )
  auction <- .latent_policy_fixture(supply = 4L)
  expect_error(
    simulate(auction$fit, auction$owner, demand_anchoring = "latent_quality",
             mc.cores = 1L), "auction or bargaining"
  )
})

test_that("the observed-share route remains the default and numerically unchanged", {
  fx <- .latent_policy_fixture(n_draws = 3L)
  set.seed(22L)
  implicit <- suppressWarnings(simulate(
    fx$fit, fx$owner, n_draws = 3L, mc.cores = 1L
  ))
  set.seed(22L)
  explicit <- suppressWarnings(simulate(
    fx$fit, fx$owner, n_draws = 3L, mc.cores = 1L,
    demand_anchoring = "observed"
  ))
  expect_equal(implicit@price_pre, explicit@price_pre)
  expect_equal(implicit@price_post, explicit@price_post)
  expect_identical(implicit@metadata$demand_anchoring, "observed")
})


test_that("historical policy-plan state hashes retain their version-1 key set", {
  fx <- .latent_policy_fixture(n_draws = 4L)
  plan <- posterior_policy_draws(fx$fit, n_draws = 3L, seed = 12L)
  v2_keys <- c("draw_indices", "quality_draws", "group_labels",
               "quality_scope", "quality_units", "quality_lookup",
               "quality_source", "quality_conditioning", "quality_mode",
               "k_cache", "role_mapping", "fixed_k", "role_offset")
  expect_identical(.lq_plan_state_fingerprint(plan),
                   digest::digest(plan[v2_keys], algo = "xxhash64", serialize = TRUE))

  historical_keys <- c("draw_indices", "quality_draws", "group_labels",
                       "quality_source", "quality_conditioning", "quality_mode",
                       "k_cache", "role_mapping", "fixed_k", "role_offset")
  historical <- plan
  historical$fit_signature$fit_signature_version <- NULL
  historical$fit_signature$quality_scope <- NULL
  historical$fit_signature$quality_units <- NULL
  historical$fit_signature$quality_lookup <- NULL
  historical$quality_scope <- NULL
  historical$quality_units <- NULL
  historical$quality_lookup <- NULL
  historical$fit_signature$fit_fingerprint <-
    .lq_fit_fingerprint(fx$fit, version = 1L)
  historical$fit_signature$latent_state_fingerprint <-
    digest::digest(historical[historical_keys], algo = "xxhash64", serialize = TRUE)
  expect_identical(.lq_plan_state_fingerprint(historical),
                   historical$fit_signature$latent_state_fingerprint)
  expect_no_error(.lq_validate_plan(fx$fit, historical))
  historical$fit_signature$fit_signature_version <- 1L
  expect_no_error(.lq_validate_plan(fx$fit, historical))
  changed <- historical
  changed$quality_draws[1L, 1L] <- changed$quality_draws[1L, 1L] + .01
  expect_error(.lq_validate_plan(fx$fit, changed), "latent states")
  changed <- plan
  changed$quality_units <- rev(changed$quality_units)
  expect_error(.lq_validate_plan(fx$fit, changed), "latent states")
})
