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
      share_basis = "relative", use_rho = TRUE, use_cutoff = FALSE,
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
        quality_parameterization = "collapsed", share_basis = "relative", use_rho = TRUE,
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
