.relative_posterior_native <- function(draws) {
  n_draw <- nrow(draws$relative_share_predictor)
  columns <- unlist(lapply(names(draws), function(name) {
    x <- draws[[name]]
    if (is.null(dim(x))) name else paste0(name, "[", seq_len(ncol(x)), "]")
  }))
  flat <- do.call(cbind, draws)
  structure(list(draws_array = array(as.numeric(flat), c(n_draw, 1L, ncol(flat)),
    dimnames = list(NULL, NULL, columns))), class = "NutpieFit")
}

.relative_posterior_fixture <- function(quality = "collapsed", rho = .6,
                                       demand = "logit") {
  n_draw <- 2L
  group <- rep(1:3, 2L)
  cell <- rep(1:2, each = 3L)
  data <- data.frame(prices = c(1, 1.2, 1.6, 1.1, 1.3, 1.7),
    quantities = c(50, 30, 20, 20, 50, 30), margins = .35,
    group = letters[group], market = "m1", year = as.character(cell),
    outside_price = 1, margin_observed = c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE))
  d <- list(N = 6L, N_group = 3L, group = group, N_market = 1L,
    market = rep.int(1L, 6L), N_market_year = 2L, market_year_idx = cell,
    mky_to_market = rep.int(1L, 2L), log_shareIn = log(data$quantities / 100),
    likelihood_mode = 1L, use_rho = as.integer(rho != 0),
    share_basis = 1L, quality_parameterization = if (quality == "sampled") 1L else 2L,
    share_likelihood = 2L, margin_likelihood = 2L, margin_obs = as.integer(data$margin_observed),
    margin_idx = which(data$margin_observed),
    log_marginNormInv = log(2) + c(.12, 0, -.08, .1, -.07, 0))
  attr(d, "group_levels") <- letters[1:3]
  attr(d, "market_levels") <- "m1"
  attr(d, "quality_scope") <- "global"
  value <- rbind(c(.1, -.2, .3, .5, .2, .1), c(-.1, .4, .2, .3, -.2, .1))
  q <- rbind(c(.2, -.1, -.1), c(-.1, .15, -.05))
  draws <- list(a_market = matrix(4, n_draw, 1L),
    s0_prob = matrix(c(.2, .25, .35, .4), n_draw, 2L),
    relative_share_predictor = value, sigma_b_strat = c(.7, .8),
    sigma_share = c(.4, .45), sigma_margin = c(.5, .55),
    rho = rep(rho, n_draw), pred_marginNormInv = matrix(2, n_draw, 6L))
  if (quality == "sampled") {
    draws$b_quality_all <- q
    draws$relative_share_predictor <- value + q[, group]
  }
  native <- .relative_posterior_native(draws)
  fit <- methods::new("BayesLogit", stanfit = native, stan_data = d,
    data = data, price_sd = 1, supply_model = 1L,
    demand_model = demand, market_side = "output", metadata = list(
      quality_parameterization = quality, quality_scope = "global",
      share_likelihood = "gaussian", margin_likelihood = "gaussian",
      share_basis = "relative", observation_basis = "ilr_quantity_shares_v1",
      market_intercept_parameterization_effective = "derived_market_year",
      structure = "simultaneous", partition_prior = "cutoff",
      market_year_labels = c("m1_1", "m1_2"), posterior_draws = n_draw,
      model_spec_version = .antitrust_bayes_model_spec_version))
  list(fit = fit, draws = draws)
}

## Material failure: a public default or wrapper can silently condition on
## observed shares, or pair recovered quality with a different slope/role row.
.policy_default_fixture <- function(quality = "collapsed", demand = "logit",
                                    structure = "simultaneous", k = 2L) {
  fixture <- .relative_posterior_fixture(quality, demand = demand)
  fit <- fixture$fit
  data <- fit@data
  native <- structure != "simultaneous"
  prepare <- if (demand == "ces") prepare_ces_stan_data else prepare_stan_data
  fit@stan_data <- prepare(
    prices = data$prices, quantities = data$quantities, margins = data$margins,
    market = data$market, year = data$year, group = data$group,
    market_side = "output", supply_model = 1L,
    quality_scope = if (native) "market" else "global",
    quality_parameterization = quality, share_basis = "relative",
    share_residual_model = "combined", use_rho = TRUE,
    use_cutoff = FALSE, use_hmt = FALSE, priors = list(),
    structure = structure, partition_prior = if (native) "ordered_k" else "cutoff",
    k_support = if (native) k else NULL,
    k_prior = if (native) rep(1 / length(k), length(k)) else NULL)
  fit@price_sd <- if (demand == "ces") fit@stan_data$log_price_sd else
    fit@stan_data$price_sd
  fit@metadata$market_year_labels <- attr(fit@stan_data, "market_year_size")$market_year_labels
  fixture$draws$a_market[, 1L] <- c(4, 5)
  if (native) {
    fit@metadata$quality_scope <- "market"
    fit@metadata$structure <- structure
    fit@metadata$partition_prior <- "ordered_k"
    fit@metadata$k_support <- k
    fit@metadata$k_prior <- rep(1 / length(k), length(k))
    fit@metadata$firm_market_ranking <- attr(fit@stan_data, "firm_market_ranking")
    fixture$draws$posterior_component_prob <- matrix(
      1 / length(k), 2L, fit@stan_data$N_partition_component)
    fixture$draws$prior_component_prob <- fixture$draws$posterior_component_prob
  }
  fit@stanfit <- .relative_posterior_native(fixture$draws)
  fit
}

test_that("estimated defaults reuse posterior quality while observed anchoring stays conditional", {
  for (quality in c("sampled", "collapsed")) for (demand in c("logit", "ces")) {
    fit <- .policy_default_fixture(quality, demand)
    saved <- serialize(fit, NULL)
    plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 291L)
    realize <- function(...) suppressWarnings(posterior_policy_models(
      fit, draw_plan = plan, year = "1", ...))
    implicit <- realize()
    latent <- realize(demand_anchoring = "latent_quality")
    observed <- realize(demand_anchoring = "observed")
    expect_identical(implicit$demand_anchoring, "latent_quality")
    expect_equal(nrow(implicit$failures), 0L)
    expect_equal(implicit$models, latent$models)
    scenario <- antitrust::counterfactual(ownership = c("a", "b", "c"))
    run <- function(...) suppressWarnings(simulate(
      fit, scenario, draw_plan = plan, year = "1", mc.cores = 1L, ...))
    default_sim <- run()
    latent_sim <- run(demand_anchoring = "latent_quality")
    observed_sim <- run(demand_anchoring = "observed")
    wrapped <- suppressWarnings(bayes_sim(
      fit, scenario, draw_plan = plan, year = "1", mc.cores = 1L))
    expect_equal(as.integer(default_sim@n_valid), 2L)
    expect_identical(default_sim@metadata$demand_anchoring, "latent_quality")
    expect_identical(wrapped@metadata$demand_anchoring, "latent_quality")
    expect_equal(default_sim@price_pre, latent_sim@price_pre)
    expect_equal(default_sim@price_post, wrapped@price_post)
    expect_equal(default_sim@price_post, latent_sim@price_post)
    expect_equal(default_sim@price_post, default_sim@price_pre, tolerance = 1e-6)
    expect_gt(max(abs(default_sim@price_pre[[1L]] - observed_sim@price_pre[[1L]])), 1e-4)
    for (pos in seq_along(plan$draw_indices)) {
      model <- implicit$models[[1L]][[pos]]@model
      conditional <- observed$models[[1L]][[pos]]@model
      expect_equal(as.numeric(conditional@pricePre), fit@data$prices[1:3])
      expect_equal(model@mcPre, conditional@mcPre, tolerance = 1e-12)
      expect_equal(as.numeric(model@pricePre), default_sim@price_pre[[1L]][pos, ],
                   ignore_attr = TRUE, tolerance = 1e-8)
      state <- attr(implicit$models[[1L]][[pos]], "bayes_policy_draw")
      expect_identical(state$draw_index, plan$draw_indices[[pos]])
      expect_equal(state$quality, as.numeric(plan$quality_draws[pos, ]))
      if (demand == "logit") {
        alpha <- c(4, 5)[plan$draw_indices[[pos]]] / fit@price_sd
        shares <- as.numeric(antitrust::calcShares(model, TRUE))
        jacobian <- -alpha * (diag(shares) - tcrossprod(shares))
        ownership <- outer(fit@data$group[1:3], fit@data$group[1:3], "==")
        foc <- shares + (jacobian * ownership) %*%
          (as.numeric(model@pricePre) - as.numeric(model@mcPre))
        expect_lt(max(abs(foc)), 1e-7)
      }
    }
    expect_identical(serialize(fit, NULL), saved)
  }
})

test_that("native singleton conditional-k defaults forward one joint quality and role state", {
  for (structure in c("fringe", "stackelberg")) {
    fit <- .policy_default_fixture(structure = structure, k = 2L)
    plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 831L)
    implicit <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan, year = "1"))
    explicit <- suppressWarnings(posterior_policy_models(
      fit, draw_plan = plan, year = "1", demand_anchoring = "latent_quality"))
    expect_equal(nrow(implicit$failures), 0L)
    expect_equal(implicit$models, explicit$models)
    run <- function(...) suppressWarnings(bayes_sim(
      fit, mcDelta = rep(.02, 3L), year = "1", mc.cores = 1L, ...))
    automatic <- run(role_mapping = "map", fixed_k = 2L)
    planned <- run(draw_plan = plan)
    latent <- run(draw_plan = plan, demand_anchoring = "latent_quality")
    expect_equal(as.integer(automatic@n_valid), 2L)
    expect_identical(automatic@metadata$demand_anchoring, "latent_quality")
    expect_true(all(vapply(automatic@metadata$k_cache,
      function(state) identical(as.integer(state[[1L]]), 2L), logical(1))))
    expect_equal(as.integer(planned@n_valid), 2L)
    expect_identical(planned@metadata$demand_anchoring, "latent_quality")
    expect_equal(planned@price_pre, latent@price_pre)
    expect_equal(planned@price_post, latent@price_post)
    for (pos in 1:2) {
      state <- attr(implicit$models[[1L]][[pos]], "bayes_policy_draw")
      expect_identical(state$draw_index, plan$draw_indices[[pos]])
      expect_equal(state$realized_k, 2L)
      expect_equal(state$quality, as.numeric(plan$quality_draws[pos, ]))
      expect_equal(planned@metadata$latent_diagnostics[[1L]][[pos]]$quality,
                   as.numeric(plan$quality_draws[pos, ]))
      expect_equal(planned@metadata$k_cache[[pos]][[1L]], 2L)
      expect_equal(as.numeric(implicit$models[[1L]][[pos]]@pricePre),
                   planned@price_pre[[1L]][pos, ], ignore_attr = TRUE)
    }
  }
})

test_that("native model realization honors fixed-plan role semantics and checks integrity", {
  for (structure in c("fringe", "stackelberg")) {
    for (controls in list(list(role_mapping = "map", fixed_k = 2L),
                          list(role_mapping = "map", role_offset = 1L))) {
      fit <- .policy_default_fixture(structure = structure,
        k = if (is.null(controls$role_offset)) 2L else 1:2)
      plan <- do.call(posterior_policy_draws,
        c(list(fit, n_draws = 2L, seed = 718L), controls))
      realized <- suppressWarnings(posterior_policy_models(
        fit, draw_plan = plan, year = "1"))
      expect_equal(nrow(realized$failures), 0L)
      expect_identical(realized$draw_plan, plan)
      for (pos in 1:2) {
        state <- attr(realized$models[[1L]][[pos]], "bayes_policy_draw")
        expect_identical(state$realized_k, plan$k_cache[[pos]][[1L]])
        expect_equal(state$quality, as.numeric(plan$quality_draws[pos, ]))
      }
      changed <- plan
      changed$metadata$role_mapping <- "threshold"
      expect_error(posterior_policy_models(fit, draw_plan = changed, year = "1"),
                   "role metadata")
      changed <- plan
      changed$k_cache[[1L]][[1L]] <- 1L
      expect_error(posterior_policy_models(fit, draw_plan = changed, year = "1"),
                   "latent states")
    }
  }
})

test_that("estimated smooth-cutoff auctions fail explicitly instead of using posted-price quality", {
  fixture <- .relative_posterior_fixture("sampled")
  fit <- .policy_default_fixture("sampled")
  fit@supply_model <- 4L
  fit@metadata$structure <- "fringe"
  fit@metadata$partition_prior <- "cutoff"
  fit@stan_data$use_cutoff <- 1L
  fixture$draws$cutoff_share <- matrix(.2, 2L, 1L)
  fit@stanfit <- .relative_posterior_native(fixture$draws)
  run <- function(...) simulate(fit, c("a", "b", "c"), year = "1",
                                n_draws = 2L, mc.cores = 1L, ...)
  expect_error(run(), "smooth-cutoff Logit auctions have no implemented cost-based value-index")
  expect_error(run(demand_anchoring = "latent_quality"),
               "smooth-cutoff Logit auctions have no implemented cost-based value-index")
  observed <- suppressWarnings(run(demand_anchoring = "observed"))
  expect_identical(observed@metadata$demand_anchoring, "observed")
  expect_equal(as.integer(observed@n_valid), 2L)
  expect_equal(observed@price_pre[[1L]][1L, ], fit@data$prices[1:3],
               ignore_attr = TRUE)
})

test_that("posterior defaults hold quality fixed through sequential policies", {
  fit <- .policy_default_fixture()
  plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 402L)
  scenario <- antitrust::counterfactual(costs = rep(.02, 3L))
  scenario <- antitrust::add_step(scenario, ownership = c("a", "a", "c"))
  result <- suppressWarnings(simulate(
    fit, scenario, draw_plan = plan, year = "1", mc.cores = 1L))
  expect_equal(as.integer(result@n_valid), 2L)
  for (pos in 1:2) {
    path <- result@sims[[1L]][[pos]]
    expect_s4_class(path, "CounterfactualPath")
    first <- path@results[[1L]]
    second <- path@results[[2L]]
    expect_equal(first@slopes, path@initial@slopes)
    expect_equal(second@slopes, path@initial@slopes)
    expect_equal(second@pricePre, first@pricePost)
    expect_equal(as.numeric(second@mcPre), as.numeric(first@mcPost))
    expect_equal(result@metadata$latent_diagnostics[[1L]][[pos]]$quality,
                 as.numeric(plan$quality_draws[pos, ]))
    alpha <- c(4, 5)[plan$draw_indices[[pos]]] / fit@price_sd
    for (model in path@results) {
      shares <- as.numeric(antitrust::calcShares(model, FALSE))
      jacobian <- -alpha * (diag(shares) - tcrossprod(shares))
      foc <- shares + (jacobian * model@ownerPost) %*%
        (as.numeric(model@pricePost) - as.numeric(model@mcPost))
      expect_lt(max(abs(foc)), 1e-7)
    }
  }
})

test_that("relative quality recovery equals dense masked joint Gaussian conditioning", {
  fixture <- .relative_posterior_fixture()
  fit <- fixture$fit
  draws <- fixture$draws
  actual <- posterior_quality(fit, draws = draws, seed = 41L)
  q <- qr.Q(qr(stats::contr.helmert(3L)))
  contrast <- matrix(0, 6L, 4L)
  contrast[1:3, 1:2] <- q
  contrast[4:6, 3:4] <- q
  incidence <- diag(3L)[fit@stan_data$group, ]
  design <- crossprod(contrast, incidence %*% q)
  for (draw in 1:2) {
    sigma <- draws$sigma_share[[draw]]
    rho <- draws$rho[[draw]]
    idx <- fit@stan_data$margin_idx
    residual <- fit@stan_data$log_shareIn - draws$relative_share_predictor[draw, ]
    residual[idx] <- residual[idx] - rho * sigma / draws$sigma_margin[[draw]] *
      (fit@stan_data$log_marginNormInv[idx] - log(2 + 1e-9))
    variance <- rep(sigma^2, 6L)
    variance[idx] <- variance[idx] * (1 - rho^2)
    covariance <- crossprod(contrast, variance * contrast)
    precision <- solve(covariance)
    posterior <- solve(diag(1 / draws$sigma_b_strat[[draw]]^2, 2L) +
      crossprod(design, precision %*% design))
    mean <- q %*% posterior %*% crossprod(design,
      precision %*% crossprod(contrast, residual))
    expect_equal(as.numeric(actual$conditional_means[draw, ]), as.numeric(mean), tolerance = 1e-10)
    expect_equal(as.numeric(actual$conditional_sds[draw, ]),
      sqrt(diag(q %*% posterior %*% t(q))), tolerance = 1e-10)
  }
  shifted <- fit
  shifted@stan_data$log_shareIn <- shifted@stan_data$log_shareIn + rep(c(19, -7), each = 3L)
  shifted@stan_data$log_marginNormInv[-fit@stan_data$margin_idx] <- 900
  invariant <- posterior_quality(shifted, draws = draws, seed = 41L)
  expect_equal(actual$draws, invariant$draws, tolerance = 1e-12)
  expect_equal(rowSums(actual$draws), c(`1` = 0, `2` = 0), tolerance = 1e-12)
})

test_that("relative plans calibrate each market-year after the same quality draw", {
  for (quality in c("sampled", "collapsed")) for (demand in c("logit", "ces")) {
    fixture <- .relative_posterior_fixture(quality, demand = demand)
    fit <- fixture$fit
    plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 71L)
    expect_identical(.lq_validate_plan(fit, plan), plan)
    value <- fixture$draws$relative_share_predictor[plan$draw_indices, , drop = FALSE]
    if (quality == "collapsed") value <- value + plan$quality_draws[, fit@stan_data$group]
    for (position in 1:2) for (cell in 1:2) {
      rows <- which(fit@stan_data$market_year_idx == cell)
      utilities <- value[position, rows] + plan$derived_market_year_intercept_draws[position, cell]
      outside <- 1 / (1 + sum(exp(utilities)))
      expect_equal(outside, fixture$draws$s0_prob[plan$draw_indices[[position]], cell], tolerance = 1e-12)
    }
    altered <- plan
    altered$derived_market_year_intercept_draws[1L, 2L] <-
      altered$derived_market_year_intercept_draws[1L, 2L] + .1
    expect_error(.lq_validate_plan(fit, altered), "fingerprint|derived|state")
    altered$fit_signature$latent_state_fingerprint <- .lq_plan_state_fingerprint(altered)
    expect_error(.lq_validate_plan(fit, altered), "inconsistent derived")
    expect_false(any(c("share_mean", "share_mean_integrated", "mu_b_market_vec") %in%
      .posterior_contract(fit, "policy_plan")$canonical))
  }
})

test_that("relative predictive fixed shares give a dense conditional margin law", {
  count <- 3500L
  rho <- .7
  ss <- c(.3, .4, .5)
  sm <- c(.6, .7, .8)
  observed <- log(c(.5, .3, .2))
  predictor <- c(.1, -.2, .4)
  expected_precision <- diag(1 / ss^2) -
    (1 / ss^2) %o% (1 / ss^2) / sum(1 / ss^2)
  cross <- rho * ss * sm
  expected_mean <- cross * as.numeric(expected_precision %*% (observed - predictor))
  expected_covariance <- diag(sm^2) - (cross %o% cross) * expected_precision
  actual <- .posterior_predict_relative_residual_draws(
    matrix(rep(predictor, each = count), count), matrix(0, count, 3L),
    matrix(rep(ss, each = count), count), matrix(rep(sm, each = count), count),
    rho = rep(rho, count), rho_active = TRUE, likelihood_mode = 1L,
    cell = rep.int(1L, 3L), seed = 193L, share_conditioning = "observed",
    observed_share = observed)
  expect_lt(max(abs(colMeans(actual$log_marginNormInv) - expected_mean) /
    sqrt(diag(expected_covariance) / count)), 5)
  covariance_mcse <- sqrt((expected_covariance^2 +
    diag(expected_covariance) %o% diag(expected_covariance)) / (count - 1L))
  expect_lt(max(abs(stats::cov(actual$log_marginNormInv) - expected_covariance) /
    covariance_mcse), 5)
  expect_gt(expected_covariance[1L, 2L], .02)
  expect_equal(exp(actual$log_shareIn[1L, ]), exp(observed), tolerance = 1e-12)
  replicated <- .posterior_predict_relative_residual_draws(
    matrix(rep(predictor, each = 4L), 4L), matrix(0, 4L, 3L),
    matrix(rep(ss, each = 4L), 4L), matrix(rep(sm, each = 4L), 4L),
    rho = rep(rho, 4L), rho_active = TRUE, likelihood_mode = 1L,
    cell = rep.int(1L, 3L), seed = 193L)
  expect_equal(rowSums(exp(replicated$log_shareIn)), rep(1, 4L), tolerance = 1e-12)
})

test_that("relative latent policies retain the structural cost primitive", {
  for (quality in c("sampled", "collapsed")) for (demand in c("logit", "ces")) {
    fixture <- .relative_posterior_fixture(quality, demand = demand)
    fit <- fixture$fit
    data <- fit@data
    prepare <- if (demand == "ces") prepare_ces_stan_data else prepare_stan_data
    fit@stan_data <- prepare(prices = data$prices, quantities = data$quantities,
      margins = data$margins, market = data$market, year = data$year, group = data$group,
      market_side = "output", supply_model = 1L, quality_scope = "global", quality_parameterization = quality,
      share_basis = "relative", share_residual_model = "combined",
      use_rho = TRUE, use_cutoff = FALSE,
      use_hmt = FALSE, priors = list())
    fit@price_sd <- if (demand == "ces") fit@stan_data$log_price_sd else fit@stan_data$price_sd
    fit@metadata$market_year_labels <-
      attr(fit@stan_data, "market_year_size")$market_year_labels
    plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 11L)
    observed <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan,
      market = "m1", year = "1", insideSize = 1000, demand_anchoring = "observed"))
    latent <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan,
      market = "m1", year = "1", insideSize = 1000, demand_anchoring = "latent_quality"))
    expect_equal(nrow(observed$failures), 0L, info = paste(quality, demand))
    expect_equal(nrow(latent$failures), 0L, info = paste(quality, demand))
    for (draw in 1:2) {
      before <- observed$models[[1L]][[draw]]@model
      after <- latent$models[[1L]][[draw]]@model
      expect_equal(after@mcPre, before@mcPre, tolerance = 1e-10)
      expect_equal(as.numeric(antitrust::calcMC(after, preMerger = TRUE)),
        as.numeric(after@mcPre), tolerance = 1e-7)
    }
    merger <- suppressWarnings(simulate(fit, c("a", "a", "c"),
      draw_plan = plan, market = "m1", year = "1", insideSize = 1000,
      demand_anchoring = "latent_quality", mc.cores = 1L))
    expect_equal(as.integer(merger@n_valid), 2L, info = paste(quality, demand))
  }
})

test_that("relative posterior prediction returns normalized aligned compositions", {
  for (quality in c("sampled", "collapsed")) {
    fixture <- .relative_posterior_fixture(quality)
    fit <- fixture$fit
    saved <- serialize(fit, NULL)
    plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 29L)
    for (conditioning in c("replicate", "observed")) {
      predicted <- posterior_predict(fit, draw_plan = plan, seed = 83L,
        share_conditioning = conditioning)
      expect_identical(predicted$draw_indices, plan$draw_indices)
      expect_identical(predicted$metadata$share_basis, "relative")
      expect_identical(predicted$metadata$observation_basis, "ilr_quantity_shares_v1")
      expect_true(all(is.finite(predicted$replicated$log_marginNormInv)))
      for (rows in list(1:3, 4:6)) {
        expect_equal(rowSums(exp(predicted$replicated$log_shareIn[, rows])),
          rep(1, 2L), ignore_attr = TRUE, tolerance = 1e-12)
        if (conditioning == "observed")
          expect_equal(unname(predicted$replicated$log_shareIn[1L, rows]),
            fit@stan_data$log_shareIn[rows], tolerance = 1e-12)
      }
    }
    expect_identical(serialize(fit, NULL), saved)
  }
})

test_that("relative exact-role policies reuse their rho-conditioned quality and k state", {
  skip_if_not_installed("trade")
  for (demand in c("logit", "ces")) for (structure in c("fringe", "stackelberg"))
    for (supply in c(1L, 3L)) {
      fixture <- .relative_posterior_fixture("collapsed", demand = demand)
      fit <- fixture$fit
      data <- fit@data
      prepare <- if (demand == "ces") prepare_ces_stan_data else prepare_stan_data
      fit@stan_data <- prepare(prices = data$prices, quantities = data$quantities,
        margins = data$margins, market = data$market, year = data$year, group = data$group,
        market_side = "output", supply_model = supply, quality_scope = "market",
        quality_parameterization = "collapsed", share_basis = "relative",
        share_residual_model = "combined", use_rho = TRUE,
        use_cutoff = FALSE, use_hmt = FALSE, priors = list(), structure = structure,
        partition_prior = "ordered_k", k_support = 1L, k_prior = 1)
      fit@supply_model <- supply
      fit@price_sd <- if (demand == "ces") fit@stan_data$log_price_sd else fit@stan_data$price_sd
      fit@metadata$quality_scope <- "market"
      fit@metadata$structure <- structure
      fit@metadata$partition_prior <- "ordered_k"
      fit@metadata$k_support <- 1L
      fit@metadata$k_prior <- 1
      fit@metadata$firm_market_ranking <- attr(fit@stan_data, "firm_market_ranking")
      fit@metadata$market_year_labels <-
        attr(fit@stan_data, "market_year_size")$market_year_labels
      fixture$draws$posterior_component_prob <- matrix(1, 2L, fit@stan_data$N_partition_component)
      fixture$draws$prior_component_prob <- fixture$draws$posterior_component_prob
      fit@stanfit <- .relative_posterior_native(fixture$draws)
      plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 119L)
      expect_true(all(vapply(plan$k_cache, function(state) all(state == 1L), logical(1))))
      prepared <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan,
        market = "m1", year = "1", insideSize = 1000, demand_anchoring = "latent_quality"))
      expect_equal(nrow(prepared$failures), 0L, info = paste(demand, structure, supply))
      merger <- suppressWarnings(simulate(fit, c("a", "a", "c"),
        draw_plan = plan, market = "m1", year = "1", insideSize = 1000,
        demand_anchoring = "latent_quality", mc.cores = 1L))
      expect_equal(as.integer(merger@n_valid), 2L, info = paste(demand, structure, supply))
      for (position in 1:2) {
        source <- prepared$models[[1L]][[position]]
        expect_identical(attr(source, "bayes_policy_draw")$draw_index, plan$draw_indices[[position]])
        promoted <- trade::as_trade_fit(source, policy = "tariff", tariffPre = 0,
          cost_basis = "effective", margin_basis = "net_revenue")
        solved <- suppressWarnings(trade::simulate(promoted, tariffPost = c(0, .2, .2)))
        before <- trade::tariff_accounts(solved, preMerger = TRUE)
        after <- trade::tariff_accounts(solved, preMerger = FALSE)
        expect_true(all(is.finite(after$price)))
        expect_equal(after$physical_cost, before$physical_cost, tolerance = 1e-12)
        if (structure == "stackelberg")
          expect_identical(solved@tariff_state$leadersPost, solved@tariff_state$leadersPre)
        else expect_identical(solved@tariff_state$corePost, solved@tariff_state$corePre)
      }
    }
})

## Material failure: a constructor's placeholder shares rescale every quantity
## even when estimated slopes, costs and equilibrium prices remain identical.
test_that("estimated posterior baselines retain quantity mass and CES budget", {
  for (demand in c("logit", "ces")) {
    fit <- .policy_default_fixture("sampled", demand)
    plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 825L)
    observed <- posterior_policy_models(fit, draw_plan = plan, year = "1",
      insideSize = 90, demand_anchoring = "observed")
    latent <- posterior_policy_models(fit, draw_plan = plan, year = "1",
      insideSize = 90)
    expect_equal(nrow(observed$failures), 0L)
    expect_equal(nrow(latent$failures), 0L)
    expected_quantity <- 90 * c(.5, .3, .2)
    for (position in seq_along(plan$draw_indices)) {
      draw <- plan$draw_indices[[position]]
      s0 <- .posterior_extract(fit)$s0_prob[draw, 1L]
      expected_mass <- if (demand == "logit") 90 / (1 - s0) else
        sum(expected_quantity * fit@data$prices[1:3]) + 90 * s0 / (1 - s0)
      model <- observed$models[[1L]][[position]]@model
      latent_model <- latent$models[[1L]][[position]]@model
      expect_equal(model@mktSize, expected_mass, tolerance = 1e-10)
      expect_equal(unname(antitrust::calcQuantities(model, TRUE)),
        expected_quantity, tolerance = 1e-10)
      expect_equal(latent_model@mktSize, expected_mass, tolerance = 1e-10)
    }
  }
})

.demand_shock_fixture <- function(quality = "collapsed", rho = .6,
                                  likelihood_mode = 1L) {
  fixture <- .relative_posterior_fixture(quality, rho)
  fit <- fixture$fit
  attr(fit@stan_data, "model_spec_version") <- .antitrust_bayes_model_spec_version
  fit@stan_data$share_residual_model <- 2L
  fit@stan_data$likelihood_mode <- likelihood_mode
  fit@metadata$share_residual_model <- "demand_shock"
  fixture$draws$sigma_demand_shock <- c(.3, .35)
  fixture$draws$sigma_measurement <- c(.2, .25)
  fixture$draws$rho_demand_shock <- rep(rho, 2L)
  fixture$draws$sigma_share <- sqrt(fixture$draws$sigma_demand_shock^2 +
                                    fixture$draws$sigma_measurement^2)
  fixture$draws$rho <- rho * fixture$draws$sigma_demand_shock /
    fixture$draws$sigma_share
  fit@stanfit <- .relative_posterior_native(fixture$draws)
  fixture$fit <- fit
  fixture
}

## Material failure: period shocks can be lost or observation error carried into
## equilibrium demand; plan reuse can also silently change the fitted state.
test_that("new residual plans carry the same recovered demand state into policy", {
  for (quality in c("sampled", "collapsed")) {
    f <- .demand_shock_fixture(quality)
    plan <- posterior_policy_draws(f$fit, n_draws = 2L, seed = 915L)
    expect_identical(plan$fit_signature$fit_signature_version, 4L)
    expect_true(any(abs(plan$demand_shock_draws) > .01))
    expect_equal(unname(rowSums(plan$demand_shock_draws[, 1:3])), c(0, 0), tolerance = 1e-12)
    expect_equal(unname(rowSums(plan$demand_shock_draws[, 4:6])), c(0, 0), tolerance = 1e-12)
    expect_equal(plan$measurement_sd_draws,
      f$draws$sigma_measurement[plan$draw_indices])
    expect_identical(antitrustBayes:::.lq_validate_plan(f$fit, plan), plan)
    changed <- plan
    changed$demand_shock_draws[1L, 1L] <- changed$demand_shock_draws[1L, 1L] + .01
    expect_error(antitrustBayes:::.lq_validate_plan(f$fit, changed), "latent states")
    models <- suppressWarnings(posterior_policy_models(f$fit, draw_plan = plan, year = "1"))
    expect_equal(nrow(models$failures), 0L)
    scenario <- antitrust::counterfactual(ownership = c("a", "b", "c"))
    policy <- suppressWarnings(simulate(f$fit, scenario, draw_plan = plan,
      year = "1", mc.cores = 1L))
    expect_equal(as.integer(policy@n_valid), 2L)
    expect_equal(policy@price_post, policy@price_pre, tolerance = 1e-6)
    for (r in seq_along(plan$draw_indices)) {
      model <- models$models[[1L]][[r]]@model
      shares <- as.numeric(antitrust::calcShares(model, TRUE))
      alpha <- 4
      jacobian <- -alpha * (diag(shares) - tcrossprod(shares))
      foc <- shares + diag(diag(jacobian)) %*%
        (as.numeric(model@pricePre) - as.numeric(model@mcPre))
      expect_lt(max(abs(foc)), 1e-7)
      expected <- antitrustBayes:::.lq_relative_full_predictor(f$fit, f$draws,
        plan$draw_indices[[r]], plan$quality_draws[r, , drop = FALSE])[1L, 1:3] +
          plan$demand_shock_draws[r, 1:3]
      reference <- model
      reference@pricePre <- f$fit@data$prices[1:3]
      ref_shares <- as.numeric(antitrust::calcShares(reference, TRUE))
      expect_equal(as.numeric(ref_shares / sum(ref_shares)),
        as.numeric(exp(expected - max(expected)) / sum(exp(expected - max(expected)))),
        tolerance = 1e-7)
    }
  }
})

## Material failure: conditioning q on omega alone overstates information, and
## a wrong Matheron covariance leaks measurement into economic demand shocks.
test_that("conditional shocks match independent Gaussian contrast moments", {
  basis <- antitrustBayes:::.relative_share_basis(rep.int(1L, 3L))
  variance <- c(.3, .5, .2)^2
  residual <- c(.5, -.3, .1)
  shift <- c(.1, -.1, .05)
  m <- .2
  K <- crossprod(basis, (variance + m^2) * basis)
  cross <- diag(variance) %*% basis
  expected_mean <- as.numeric(basis %*% crossprod(basis,
    shift + cross %*% solve(K, crossprod(basis, residual - shift))))
  expected_cov <- basis %*% crossprod(basis,
    (diag(variance) - cross %*% solve(K, t(cross))) %*% basis) %*% t(basis)
  set.seed(91)
  states <- replicate(5000L, antitrustBayes:::.lq_conditional_demand_shock(
    residual, shift, variance, m, rep.int(1L, 3L)))
  expect_equal(rowMeans(states), expected_mean, tolerance = .008)
  covariance_mcse <- sqrt((expected_cov^2 + diag(expected_cov) %o%
    diag(expected_cov)) / (ncol(states) - 1L))
  expect_lt(max(abs(stats::cov(t(states)) - expected_cov) / covariance_mcse), 5)
  f <- .demand_shock_fixture()
  q <- posterior_quality(f$fit, draws = f$draws, seed = 94L)
  legacy <- f$fit
  legacy@metadata$share_residual_model <- "combined"
  legacy@stan_data$share_residual_model <- 1L
  reference <- posterior_quality(legacy, draws = f$draws, seed = 94L)
  expect_equal(q$draws, reference$draws)
})

## Material failure: omitted share equations must not constrain the economic
## shock or persistent q using data that the likelihood did not observe.
test_that("margin-only demand states ignore observed shares", {
  f <- .demand_shock_fixture(likelihood_mode = 3L)
  p <- posterior_policy_draws(f$fit, n_draws = 2L, seed = 31L)
  f$fit@stan_data$log_shareIn <- rev(f$fit@stan_data$log_shareIn) + 100
  q <- posterior_policy_draws(f$fit, n_draws = 2L, seed = 31L)
  expect_equal(p$quality_draws, q$quality_draws)
  expect_equal(p$demand_shock_draws, q$demand_shock_draws)
})

## Material failure: using a common total correlation with market shock scales
## gives the wrong share-margin covariance and hence the wrong posterior q.
test_that("market residual scales recover q using the economic shock covariance", {
  f <- .demand_shock_fixture()
  f$fit@stan_data$use_market_resid_scale <- 1L
  f$draws$sigma_demand_shock_market <- matrix(c(.6, .7), 2L, 1L)
  f$draws$sigma_share_market <- sqrt(f$draws$sigma_demand_shock_market^2 +
    f$draws$sigma_measurement^2)
  f$draws$sigma_margin_market <- matrix(c(.7, .8), 2L, 1L)
  f$draws$rho_market <- f$draws$rho_demand_shock *
    f$draws$sigma_demand_shock_market / f$draws$sigma_share_market
  f$fit@stanfit <- .relative_posterior_native(f$draws)
  q <- posterior_quality(f$fit, draws = f$draws, seed = 5L)
  for (r in 1:2) {
    residual <- f$fit@stan_data$log_shareIn - f$draws$relative_share_predictor[r, ]
    idx <- f$fit@stan_data$margin_idx
    omega <- f$draws$sigma_demand_shock_market[r, 1L]
    rho <- f$draws$rho_demand_shock[[r]]
    variance <- rep(omega^2 + f$draws$sigma_measurement[[r]]^2, 6L)
    residual[idx] <- residual[idx] - rho * omega / f$draws$sigma_margin_market[r, 1L] *
      (f$fit@stan_data$log_marginNormInv[idx] - log(2 + 1e-9))
    variance[idx] <- omega^2 * (1 - rho^2) + f$draws$sigma_measurement[[r]]^2
    expected <- .relative_share_gaussian(residual, variance,
      f$fit@stan_data$group, rep.int(1L, 3L), f$fit@stan_data$market_year_idx,
      tau = f$draws$sigma_b_strat[[r]], return_moments = TRUE)
    expect_equal(as.numeric(q$conditional_means[r, ]), expected$conditional_mean)
    expect_equal(as.numeric(q$conditional_sds[r, ]), sqrt(expected$conditional_variance))
  }
  p <- posterior_policy_draws(f$fit, n_draws = 2L, seed = 51L)
  expect_identical(.lq_validate_plan(f$fit, p), p)
  invalid <- f$fit
  invalid@metadata$share_residual_model <- "combined"
  expect_error(posterior_policy_draws(invalid), "Inconsistent")
  f$draws$sigma_demand_shock_market[1L, 1L] <- -.1
  f$fit@stanfit <- .relative_posterior_native(f$draws)
  expect_error(posterior_policy_draws(f$fit, n_draws = 2L), "market-specific")
})

## Material failure: exact-k marginalization can condition q or d on another
## component's FOC and then solve policy under a newly drawn role state.
test_that("demand-shock native policies keep one role and demand state", {
  for (demand in c("logit", "ces")) for (structure in c("fringe", "stackelberg")) {
    fit <- .policy_default_fixture(structure = structure, demand = demand, k = c(1L, 2L))
    fit@metadata$share_residual_model <- "demand_shock"
    fit@stan_data$share_residual_model <- 2L
    attr(fit@stan_data, "share_residual_model") <- "demand_shock"
    draws <- .bayes_backend_extract(fit@stanfit)
    draws$sigma_demand_shock <- c(.3, .35)
    draws$sigma_measurement <- c(.2, .25)
    draws$rho_demand_shock <- c(.6, .6)
    draws$sigma_share <- sqrt(draws$sigma_demand_shock^2 + draws$sigma_measurement^2)
    draws$rho <- draws$rho_demand_shock * draws$sigma_demand_shock / draws$sigma_share
    fit@stanfit <- .relative_posterior_native(draws)
    plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 209L)
    models <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan, year = "1"))
    expect_equal(nrow(models$failures), 0L)
    scenario <- antitrust::counterfactual(ownership = c("a", "b", "c"))
    policy <- suppressWarnings(simulate(fit, scenario, draw_plan = plan,
      year = "1", mc.cores = 1L))
    expect_equal(as.integer(policy@n_valid), 2L)
    expect_equal(policy@price_post, policy@price_pre, tolerance = 1e-6)
    for (r in seq_along(plan$draw_indices)) {
      state <- attr(models$models[[1L]][[r]], "bayes_policy_draw")
      expect_identical(state$k, as.integer(plan$k_cache[[as.character(plan$draw_indices[[r]])]][[1L]]))
      expect_equal(state$demand_shock, as.numeric(plan$demand_shock_draws[r, 1:3]))
    }
  }
})

## Material failure: a recovery-only economic-scale alias can invalidate a
## frozen v2/v3 plan despite identical fitted economic draws and latent states.
test_that("historical total-fraction plans retain their pre-alias fingerprint", {
  f <- .relative_posterior_fixture()
  f$fit@metadata$quality_scale_parameterization_effective <- "total_fraction"
  f$fit@stan_data$quality_scale_parameterization <- 2L
  f$draws$sigma_b_strat_economic <- f$draws$sigma_b_strat
  f$fit@stanfit <- .relative_posterior_native(f$draws)
  plan <- posterior_policy_draws(f$fit, n_draws = 2L, seed = 43L)
  for (version in c(2L, 3L)) {
    old <- plan
    old$fit_signature$fit_signature_version <- version
    old$fit_signature$fit_fingerprint <- .lq_fit_fingerprint(f$fit,
      version = version, omit_redundant_quality_alias = TRUE)
    old$fit_signature$latent_state_fingerprint <- .lq_plan_state_fingerprint(old)
    expect_identical(.lq_validate_plan(f$fit, old), old)
    corrupted <- old
    corrupted$quality_draws[1L, 1L] <- corrupted$quality_draws[1L, 1L] + .01
    expect_error(.lq_validate_plan(f$fit, corrupted), "latent states")
  }
  expect_error(.lq_fit_fingerprint(f$fit, version = 4L,
    omit_redundant_quality_alias = TRUE), "version 2 or 3")
})

## Material failure: sampled quality must not excuse missing new residual states,
## and share-only policies must not condition on omitted margin observations.
test_that("demand residual contracts respect omitted equations and required states", {
  f <- .demand_shock_fixture("sampled", likelihood_mode = 2L)
  p <- posterior_policy_draws(f$fit, n_draws = 2L, seed = 72L)
  changed <- f$fit
  changed@stan_data$log_marginNormInv <- rep(100, 6L)
  q <- posterior_policy_draws(changed, n_draws = 2L, seed = 72L)
  expect_equal(p$quality_draws, q$quality_draws)
  expect_equal(p$demand_shock_draws, q$demand_shock_draws)
  f$draws$sigma_measurement <- NULL
  f$fit@stanfit <- .relative_posterior_native(f$draws)
  expect_error(posterior_policy_draws(f$fit), "sigma_measurement")
})
