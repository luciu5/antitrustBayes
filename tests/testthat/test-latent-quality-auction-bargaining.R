## Material failures guarded here: applying posted-price Logit inversion to
## auction win probabilities, substituting antitrust's different multiproduct
## bargaining equation, and recovering costs again after quality is realized.
## Posterior extraction is mocked; these are economic policy tests, not HMC.
.family_quality_fixture <- function(demand, supply, side = "output") {
  group <- c("a", "a", "b")
  data <- data.frame(prices = c(1.1, 1.2, 1.3), quantities = c(40, 30, 20),
    margins = c(.5, .4, .3), market = "m1", year = "1", group = group,
    loo_cluster = "m1", margin_observed = TRUE, outside_price = 1)
  sd <- list(quality_parameterization = 1L, N_group = 2L,
             group = c(1L, 1L, 2L))
  attr(sd, "model_spec_version") <- 14L
  attr(sd, "quality_scope") <- "global"
  attr(sd, "group_levels") <- c("a", "b")
  fit <- methods::new("BayesLogit",
    stanfit = NULL, stan_data = sd, data = data,
    supply_model = as.integer(supply), demand_model = demand,
    market_side = side, price_sd = 1,
    metadata = list(quality_parameterization = "sampled", posterior_draws = 2L,
      structure = "simultaneous", partition_prior = "cutoff",
      model_spec_version = 14L))
  q <- rbind(c(.20, -.20), c(-.15, .15))
  draws <- list(a_market = matrix(if (demand == "logit") 6 else 9, 2L, 1L),
    s0_prob = matrix(.2, 2L, 1L), b_quality_all = q,
    share_mean = matrix(rep(log(c(40, 30, 20) / 90), 2L), 2L, byrow = TRUE) +
      q[, c(1L, 1L, 2L)], bargaining_power = rep(.3, 2L))
  list(fit = fit, draws = draws, owners = group)
}

## Independent conduct identities use solved shares and ownership, rather than
## calcMargins(), to catch a wrong family, sign, or bargaining reconstruction.
.family_expected_markup <- function(model, supply, pre) {
  shares <- antitrust::calcShares(model, preMerger = pre,
                                  revenue = methods::is(model, "CES"))
  owner <- if (pre) model@ownerPre else model@ownerPost
  firm <- as.numeric(owner %*% shares)
  if (methods::is(model, "CES")) {
    gamma <- model@slopes$gamma
    proportion <- if (supply == 4L) 1 - (1 - firm)^(1 / (gamma - 1)) else
      (1 - .3) / (gamma - (gamma - 1) * firm)
    return(proportion * if (pre) model@pricePre else model@pricePost)
  }
  alpha <- abs(model@slopes$alpha)
  if (supply == 4L) return(-log1p(-firm) / (alpha * firm))
  kappa <- .3 / (1 - .3)
  -log1p(-firm) /
    (alpha * (kappa * firm - (1 - firm) * log1p(-firm)))
}

test_that("nonzero posterior auction and bargaining states retain primitives and solve policies", {
  for (demand in c("logit", "ces")) for (supply in c(4L, 5L)) {
    for (side in if (demand == "logit") c("output", "input") else "output") {
      fixture <- .family_quality_fixture(demand, supply, side)
      testthat::local_mocked_bindings(
        .cmdstan_extract = function(...) fixture$draws,
        .package = "antitrustBayes")
      fit <- fixture$fit
      plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 29L)
      no_policy <- suppressWarnings(simulate(fit, fixture$owners,
        draw_plan = plan, mc.cores = 1L))
      policy <- antitrust::counterfactual(ownership = rep("merged", 3L),
                                         costs = rep(.05, 3L))
      changed <- suppressWarnings(simulate(fit, policy,
        draw_plan = plan, mc.cores = 1L))
      prepared <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan))
      observed <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan,
                                                           demand_anchoring = "observed"))
      case <- paste(demand, supply, side)
      expect_equal(as.integer(no_policy@n_valid), 2L, info = case)
      expect_equal(as.integer(changed@n_valid), 2L, info = case)
      expect_equal(nrow(prepared$failures), 0L, info = case)
      expect_equal(nrow(observed$failures), 0L, info = case)
      expect_identical(no_policy@metadata$demand_anchoring, "latent_quality")
      expect_identical(prepared$draw_plan, plan)
      for (i in seq_len(2L)) {
        base <- no_policy@sims$m1[[i]]
        post <- changed@sims$m1[[i]]
        realized <- prepared$models$m1[[i]]@model
        original <- observed$models$m1[[i]]@model
        base_info <- no_policy@metadata$latent_diagnostics$m1[[i]]
        post_info <- changed@metadata$latent_diagnostics$m1[[i]]
        expect_equal(base_info$quality, as.numeric(plan$quality_draws[i, ]))
        expect_equal(base_info$quality, post_info$quality)
        expect_equal(unname(base@mcPre), unname(original@mcPre))
        expect_equal(post@mcPre, base@mcPre)
        expect_equal(unname(realized@mcPre), unname(base@mcPre))
        expect_equal(unname(realized@pricePre), unname(base@pricePre), tolerance = 1e-8)
        ## Downstream consumers may reconstruct demand from the public fit
        ## parameters. A stale observed-share intercept would discard q.
        parameters <- prepared$models$m1[[i]]@parameters
        if (demand == "ces") {
          odds <- parameters$meanval * realized@pricePre^(-parameters$gamma)
          implied_share <- odds / (1 + sum(odds))
        } else {
          eta <- parameters$meanval
          if (supply != 4L)
            eta <- eta + parameters$alpha * (realized@pricePre - realized@priceOutside)
          implied_share <- exp(eta) / (1 + sum(exp(eta)))
        }
        expect_equal(as.numeric(implied_share),
          as.numeric(antitrust::calcShares(realized, TRUE, revenue = FALSE)),
          tolerance = 1e-10)
        replay <- suppressWarnings(antitrust::simulate(prepared$models$m1[[i]],
                                                       fixture$owners))
        expect_equal(unname(replay@pricePost), unname(base@pricePre), tolerance = 1e-7)
        expect_gt(max(abs(base_info$reference_quantity_shares - c(40, 30, 20) / 90 * .8)),
                  1e-4)
        expect_equal(unname(base@pricePre), unname(base@pricePost), tolerance = 1e-7)
        expect_equal(no_policy@price_change$m1[i, ], rep(0, 3L),
                     ignore_attr = TRUE, tolerance = 1e-7)
        expected_cost <- if (demand == "logit" && supply == 4L)
          post@mcPre + .05 else post@mcPre * 1.05
        expect_equal(unname(post@mcPost), unname(expected_cost))
        for (pre in c(TRUE, FALSE)) {
          price <- if (pre) post@pricePre else post@pricePost
          cost <- if (pre) post@mcPre else post@mcPost
          residual <- if (side == "output") price - cost else cost - price
          expect_lt(max(abs(residual - .family_expected_markup(post, supply, pre))),
                    1e-6)
        }
        if (demand == "logit" && supply == 4L) {
          expect_identical(base_info$demand_representation,
                           "auction_cost_based_value_index")
          alpha <- post@slopes$alpha
          eta <- post@slopes$meanval + alpha * (post@mcPost - post@mcPre)
          expected_share <- exp(eta) / (1 + sum(exp(eta)))
          expect_equal(as.numeric(antitrust::calcShares(post, FALSE)),
                       as.numeric(expected_share), tolerance = 1e-10)
          expect_equal(base@slopes$meanval, post@slopes$meanval)
          expect_equal(base_info$auction_value,
                       as.numeric(post@slopes$meanval - alpha * post@mcPre))
          perturbed <- post
          perturbed@pricePost <- post@pricePost * 1.2
          expect_equal(antitrust::calcShares(perturbed, FALSE),
                       antitrust::calcShares(post, FALSE))
        }
        if (demand == "logit" && supply == 5L) {
          expect_s4_class(realized, "BayesConductLogit")
          expect_equal(realized@bargaining_power, .3)
        }
      }
    }
  }
})

test_that("smooth-cutoff Logit auctions never use posted-price latent reconstruction", {
  params <- bayes_params(alpha = 6, s0 = .2, supply_model = 4L,
    structure = "fringe", partition_prior = "cutoff", use_cutoff = TRUE,
    cutoff_share = .2)
  market <- bayes_market(prices = c(1.1, 1.2, 1.3), quantities = c(40, 30, 20),
    margins = c(.5, .4, .3), group = c("a", "a", "b"),
    outside_price = 1, market_side = "output")
  observed <- suppressWarnings(as_antitrust_fit(build_conduct_object(
    params, 1L, market, ownerPre = market$group, ownerPost = market$group)))
  expect_error(.lq_replace_demand(observed@model, log(c(.5, .3, .2)), .2,
    market$prices, 1, c(.1, -.1), 1L), "no implemented cost-based value-index")
})

test_that("posterior auction paths carry fixed quality, value indices, and exited choices", {
  for (side in c("output", "input")) {
    fixture <- .family_quality_fixture("logit", 4L, side)
    testthat::local_mocked_bindings(.cmdstan_extract = function(...) fixture$draws,
                                   .package = "antitrustBayes")
    plan <- posterior_policy_draws(fixture$fit, n_draws = 2L, seed = 29L)
    scenarios <- list(
      antitrust::add_step(antitrust::add_step(
        antitrust::counterfactual(costs = rep(.05, 3L)), ownership = rep("merged", 3L))),
      antitrust::add_step(antitrust::counterfactual(exit = "b"),
                          ownership = rep("merged", 3L)))
    for (scenario in scenarios) {
      result <- suppressWarnings(simulate(fixture$fit, scenario, draw_plan = plan,
        labels = c("a", "b", "c"), mc.cores = 1L))
      expect_equal(as.integer(result@n_valid), 2L)
      for (i in 1:2) {
        path <- result@sims$m1[[i]]
        expect_s4_class(path, "CounterfactualPath")
        quality <- attr(path, "bayes_latent_quality")
        expect_equal(quality$quality, as.numeric(plan$quality_draws[i, ]))
        value <- quality$auction_value
        for (j in seq_along(path@results)) {
          model <- antitrust::result_at(path, j)
          expect_equal(as.numeric(model@slopes$meanval -
            model@slopes$alpha * model@mcPre), value, tolerance = 1e-10)
          for (pre in c(TRUE, FALSE)) {
            active <- if (pre) attr(model, "antitrust_subset_pre", exact = TRUE) else
              model@subset
            if (is.null(active)) active <- rep(TRUE, 3L)
            cost <- if (pre) model@mcPre else model@mcPost
            price <- if (pre) model@pricePre else model@pricePost
            owner <- if (pre) model@ownerPre else model@ownerPost
            weights <- exp(value[active] + model@slopes$alpha * cost[active])
            expected <- weights / (1 + sum(weights))
            shares <- antitrust::calcShares(model, preMerger = pre)
            expect_equal(as.numeric(shares[active]), as.numeric(expected), tolerance = 1e-10)
            expect_true(all(is.na(shares[!active])))
            firm <- as.numeric(owner[active, active, drop = FALSE] %*% expected)
            markup <- -log1p(-firm) / (abs(model@slopes$alpha) * firm)
            actual <- if (side == "output") price[active] - cost[active] else
              cost[active] - price[active]
            expect_equal(as.numeric(actual), as.numeric(markup), tolerance = 1e-8)
          }
          if (j > 1L) expect_equal(antitrust::calcShares(model, TRUE),
            antitrust::calcShares(antitrust::result_at(path, j - 1L), FALSE))
        }
        if (length(path@results) == 3L) {
          terminal <- antitrust::final_result(path)
          expect_equal(terminal@pricePre, terminal@pricePost, tolerance = 1e-10)
        }
      }
    }
  }
})

## Material failure: adding a period demand state must preserve cost-based
## auction value inversion and bargaining primitives across sequential policies.
test_that("demand shocks reach auction and bargaining sequential policy routes", {
  for (demand in c("logit", "ces")) for (supply in c(4L, 5L))
    for (side in if (demand == "logit") c("output", "input") else "output") {
      fixture <- .family_quality_fixture(demand, supply, side)
      fit <- fixture$fit
      d <- fit@stan_data
      d$N <- 3L; d$N_market <- 1L; d$market <- rep.int(1L, 3L)
      d$N_market_year <- 1L; d$market_year_idx <- rep.int(1L, 3L)
      d$mky_to_market <- 1L; d$share_basis <- 1L
      d$share_residual_model <- 2L; d$likelihood_mode <- 1L
      d$share_likelihood <- 2L; d$margin_likelihood <- 2L; d$use_rho <- 1L
      d$margin_obs <- rep.int(1L, 3L); d$margin_idx <- 1:3
      d$log_shareIn <- log(c(40, 30, 20) / 90)
      d$log_marginNormInv <- log(c(2, 2.5, 3))
      attr(d, "model_spec_version") <- .antitrust_bayes_model_spec_version
      fit@stan_data <- d
      fit@metadata$model_spec_version <- .antitrust_bayes_model_spec_version
      fit@metadata$share_residual_model <- "demand_shock"
      fit@metadata$share_basis <- "relative"
      fit@metadata$share_likelihood <- "gaussian"
      fit@metadata$margin_likelihood <- "gaussian"
      fit@metadata$market_intercept_parameterization_effective <- "derived_market_year"
      draws <- fixture$draws
      draws$relative_share_predictor <- draws$share_mean
      draws$sigma_demand_shock <- c(.3, .35)
      draws$sigma_measurement <- c(.2, .25)
      draws$rho_demand_shock <- c(.4, .4)
      draws$sigma_share <- sqrt(draws$sigma_demand_shock^2 + draws$sigma_measurement^2)
      draws$sigma_margin <- c(.4, .5)
      draws$rho <- draws$rho_demand_shock * draws$sigma_demand_shock / draws$sigma_share
      draws$pred_marginNormInv <- matrix(2, 2L, 3L)
      testthat::local_mocked_bindings(.cmdstan_extract = function(...) draws,
                                    .package = "antitrustBayes")
      plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 420L)
      path <- antitrust::add_step(antitrust::counterfactual(costs = rep(.02, 3L)),
                                  ownership = rep("merged", 3L))
      result <- suppressWarnings(simulate(fit, path, draw_plan = plan, mc.cores = 1L))
      prepared <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan))
      expect_equal(as.integer(result@n_valid), 2L)
      expect_equal(nrow(prepared$failures), 0L)
      for (r in 1:2) {
        state <- attr(prepared$models[[1L]][[r]], "bayes_policy_draw")
        expect_equal(state$demand_shock, as.numeric(plan$demand_shock_draws[r, ]))
        expect_equal(result@metadata$latent_diagnostics[[1L]][[r]]$demand_shock,
                     state$demand_shock)
        post <- antitrust::final_result(result@sims[[1L]][[r]])
        for (pre in c(TRUE, FALSE)) {
          price <- if (pre) post@pricePre else post@pricePost
          cost <- if (pre) post@mcPre else post@mcPost
          residual <- if (side == "output") price - cost else cost - price
          expect_lt(max(abs(residual - .family_expected_markup(post, supply, pre))), 1e-6)
        }
      }
    }
})
