## Deterministic posterior states with a nonzero sampled market intercept.
## These checks establish reconstruction identities, not sampler recovery.
.single_intercept_policy_fixture <- function(candidate = "B", shift = 0) {
  price <- c(1.2, 1.35, 1.5, 1.8)
  conditional_share <- c(.38, .29, .20, .13)
  outside_share <- .3
  alpha <- 6
  inside_share <- (1 - outside_share) * conditional_share
  k <- 2L
  h <- switch(candidate,
    B = 1 / (1 - inside_share),
    C = 1 + inside_share / outside_share,
    CF = c(1 + inside_share[1:2] / outside_share, 1, 1),
    CS = c(1 + inside_share[1:2] / (outside_share + sum(inside_share[3:4])),
           1 + inside_share[3:4] / outside_share))
  markup <- h / alpha
  structure <- if (candidate == "CF") "fringe" else
    if (candidate == "CS") "stackelberg" else "simultaneous"
  supply <- if (candidate == "B") 1L else 3L
  data <- data.frame(prices = price, quantities = 1000 * conditional_share,
    margins = markup / price, market = "single", year = "period",
    group = letters[1:4], outside_price = 1, margin_observed = TRUE)
  stan_data <- prepare_stan_data(share_basis = "absolute", data$prices, data$quantities, data$margins,
    data$market, data$year, data$group, supply_model = supply,
    market_side = "output", use_rho = FALSE, structure = structure,
    partition_prior = if (structure == "simultaneous") "cutoff" else "ordered_k",
    k_support = if (structure == "simultaneous") NULL else k,
    k_prior = if (structure == "simultaneous") NULL else 1,
    quality_parameterization = "sampled", quality_scope = "market",
    share_likelihood = "gaussian", margin_likelihood = "gaussian")
  intercept <- mean(log(conditional_share)) - stats::qlogis(outside_share)
  quality <- log(conditional_share) - stats::qlogis(outside_share) - intercept +
    alpha * stats::sd(price) * stan_data$price_std
  values <- c(`a_market[1]` = alpha * stats::sd(price), `s0_prob[1]` = outside_share,
    `mu_b_market_vec[1]` = intercept + shift, `b_market[1]` = intercept + shift,
    `market_intercept_location[1]` = intercept + shift,
    `gamma_market[1]` = intercept + shift,
    sigma_b_strat = .5, sigma_share = .2, sigma_margin = .3,
    `posterior_component_prob[1]` = 1, `prior_component_prob[1]` = 1)
  values <- c(values,
    stats::setNames(quality, paste0("b_quality_all[", 1:4, "]")),
    stats::setNames(log(conditional_share) + shift, paste0("share_mean[", 1:4, "]")),
    stats::setNames(1 / data$margins, paste0("pred_marginNormInv[", 1:4, "]")))
  posterior <- array(rep(values, each = 2L), c(2L, 1L, length(values)),
    dimnames = list(iteration = 1:2, chain = "1", variable = names(values)))
  native <- structure(list(draws_array = posterior), class = "NutpieFit")
  fit <- methods::new("BayesLogit", stanfit = native, stan_data = stan_data,
    data = data, supply_model = supply, demand_model = "logit",
    market_side = "output", price_sd = stats::sd(price),
    metadata = list(backend = "nutpie", inference = "hmc", structure = structure,
      partition_prior = attr(stan_data, "partition_prior"),
      firm_market_ranking = attr(stan_data, "firm_market_ranking"),
      k_support = if (structure == "simultaneous") NULL else k,
      k_prior = if (structure == "simultaneous") NULL else 1,
      quality_parameterization = "sampled", quality_scope = "market",
      share_likelihood = "gaussian", margin_likelihood = "gaussian",
      market_intercept_parameterization = "sampled", posterior_draws = 2L,
      model_spec_version = .antitrust_bayes_model_spec_version))
  list(fit = fit, intercept = intercept + shift, quality = quality,
    alpha = alpha, outside_share = outside_share, inside_share = inside_share,
    markup = markup, cost = price - markup, candidate = candidate)
}

test_that("single-market prediction and latent scores inherit a nonzero intercept once", {
  fixture <- .single_intercept_policy_fixture()
  fit <- fixture$fit
  expect_gt(abs(fixture$intercept), .1)
  expect_equal(sum(fixture$quality), 0, tolerance = 1e-12)
  draws <- .posterior_extract(fit, "posterior_predictive")
  plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 822L)
  unit <- .lq_all_units(fit)$info[[1L]]
  market <- bayes_market(fit@data$prices, fit@data$quantities, fit@data$margins,
    group = fit@data$group, outside_price = 1, market_side = "output")
  params <- .estimated_draw_params(fit, draws, plan$draw_indices[1L], unit)
  score <- .lq_reference_score(fit, draws, plan$draw_indices[1L], unit, plan, market, params)
  expected <- log(fit@data$quantities / sum(fit@data$quantities))
  expect_equal(score$score, expected, tolerance = 1e-12)
  means <- .posterior_predict_mean_state(fit, draws, plan, .posterior_predict_mask_and_mode(fit))
  expect_equal(unname(means$share_mean), matrix(rep(expected, 2L), 2L, byrow = TRUE),
    tolerance = 1e-12)
  expect_null(plan$market_intercept_draws)

  shifted <- .single_intercept_policy_fixture(shift = .45)$fit
  other_plan <- posterior_policy_draws(shifted, n_draws = 2L, seed = 822L)
  first <- posterior_predict(fit, draw_plan = plan, seed = 610L)
  second <- posterior_predict(shifted, draw_plan = other_plan, seed = 610L)
  expect_true(is.matrix(first$replicated$log_shareIn))
  expect_equal(unname(second$replicated$log_shareIn - first$replicated$log_shareIn), matrix(.45, 2L, 4L),
    tolerance = 1e-12)
  expect_equal(second$replicated$share_rep_normalized, first$replicated$share_rep_normalized, tolerance = 1e-12)
  expect_equal(second$replicated$log_marginNormInv, first$replicated$log_marginNormInv, tolerance = 1e-12)
})

test_that("nonzero single-market intercepts preserve observed equilibria costs and roles", {
  skip_if_not_installed("coordination", "0.0.0.9006")
  for (candidate in c("B", "C", "CF", "CS")) {
    fixture <- .single_intercept_policy_fixture(candidate)
    fit <- fixture$fit
    expect_true(all(fixture$cost > 0), info = candidate)
    plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 915L)
    routes <- lapply(c("observed", "latent_quality"), function(anchor)
      suppressWarnings(posterior_policy_models(fit, draw_plan = plan,
        demand_anchoring = anchor, insideSize = 1000)))
    for (route in routes) {
      expect_equal(nrow(route$failures), 0L, info = candidate)
      for (position in seq_along(plan$draw_indices)) {
        result <- route$models[[1L]][[position]]
        model <- if (methods::is(result, "AntitrustFit")) result@model else result
        expect_equal(as.numeric(model@pricePre), fit@data$prices, tolerance = 1e-7, info = candidate)
        expect_equal(as.numeric(antitrust::calcShares(model, TRUE, revenue = FALSE)),
          fixture$inside_share, tolerance = 1e-7, info = candidate)
        expect_equal(as.numeric(model@mcPre), fixture$cost, tolerance = 1e-7, info = candidate)
        expect_equal(as.numeric(antitrust::calcMargins(model, TRUE, level = TRUE)),
          fixture$markup, tolerance = 1e-7, info = candidate)
        expect_identical(attr(result, "bayes_policy_draw")$draw_index, plan$draw_indices[position])
        if (candidate %in% c("CF", "CS")) {
          expect_equal(attr(result, "bayes_policy_draw")$realized_k, 2L, info = candidate)
          expect_setequal(if (candidate == "CF") model@corePre else model@leadersPre, letters[1:2])
        }
      }
    }
    merger <- lapply(c("observed", "latent_quality"), function(anchor)
      suppressWarnings(simulate(fit, c("a", "a", "c", "d"), draw_plan = plan,
        demand_anchoring = anchor, insideSize = 1000, mc.cores = 1L)))
    expect_equal(as.integer(merger[[1L]]@n_valid), 2L, info = candidate)
    expect_equal(as.integer(merger[[2L]]@n_valid), 2L, info = candidate)
    expect_equal(merger[[1L]]@price_post, merger[[2L]]@price_post, tolerance = 1e-7, info = candidate)
  }
})

test_that("collapsed quality conditions on the single sampled intercept without refitting it", {
  fixture <- .single_intercept_policy_fixture()
  fit <- fixture$fit
  fit@data$group <- c("a", "a", "b", "c")
  fit@stan_data <- prepare_stan_data(share_basis = "absolute", fit@data$prices, fit@data$quantities, fit@data$margins,
    fit@data$market, fit@data$year, fit@data$group, supply_model = 1L,
    market_side = "output", use_rho = FALSE, quality_parameterization = "collapsed",
    quality_scope = "market", share_likelihood = "gaussian", margin_likelihood = "gaussian")
  fit@metadata$quality_parameterization <- "collapsed"
  intercept <- c(.45, -.35)
  offset <- c(-.9, -.3, .15, .4)
  means <- outer(intercept, rep(1, 4L)) + matrix(rep(offset, 2L), 2L, byrow = TRUE)
  draws <- list(share_mean = means, mu_b_market_vec = intercept,
    sigma_b_strat = rep(.5, 2L), sigma_share = rep(.3, 2L))
  result <- posterior_quality(fit, draws = draws, seed = 411L, conditioning = "share")
  incidence <- model.matrix(~ factor(fit@stan_data$obs_to_quality) - 1L)
  basis <- contr.helmert(ncol(incidence))
  basis <- sweep(basis, 2L, sqrt(colSums(basis^2)), "/")
  design <- incidence %*% basis
  covariance <- solve(diag(1 / .5^2, ncol(basis)) + crossprod(design) / .3^2)
  for (row in 1:2) {
    expected_mean <- basis %*% covariance %*% crossprod(design,
      fit@stan_data$log_shareIn - means[row, ]) / .3^2
    expected_sd <- sqrt(diag(basis %*% covariance %*% t(basis)))
    expect_equal(as.numeric(result$conditional_means[row, ]), as.numeric(expected_mean), tolerance = 1e-12)
    expect_equal(as.numeric(result$conditional_sds[row, ]), as.numeric(expected_sd), tolerance = 1e-12)
  }
  expect_gt(max(abs(result$conditional_means[1L, ] - result$conditional_means[2L, ])), .01)
  expect_null(result$market_intercept_draws)
})

test_that("nonzero single-market intercept states retain tariff costs and inherited roles", {
  skip_if_not_installed("trade", "0.8.9.9002")
  for (candidate in c("B", "C", "CF", "CS")) {
    fixture <- .single_intercept_policy_fixture(candidate)
    fit <- fixture$fit
    plan <- posterior_policy_draws(fit, n_draws = 1L, seed = 618L)
    prices <- lapply(c("observed", "latent_quality"), function(anchor) {
      realized <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan,
        demand_anchoring = anchor, insideSize = 1000))
      expect_equal(nrow(realized$failures), 0L, info = candidate)
      source <- realized$models[[1L]][[1L]]
      promoted <- trade::as_trade_fit(source, policy = "tariff", tariffPre = 0,
        cost_basis = "effective", margin_basis = "net_revenue")
      roles <- if (candidate == "CF") list(corePost = promoted@tariff_state$corePre) else
        if (candidate == "CS") list(leadersPost = promoted@tariff_state$leadersPre) else list()
      solved <- suppressWarnings(do.call(trade::simulate, c(list(promoted,
        ownerPost = fit@data$group, tariffPost = c(0, 0, .2, .2)), roles)))
      before <- trade::tariff_accounts(solved, preMerger = TRUE)
      after <- trade::tariff_accounts(solved, preMerger = FALSE)
      expect_equal(before$physical_cost, fixture$cost, tolerance = 1e-7, info = candidate)
      expect_equal(after$physical_cost, before$physical_cost, tolerance = 1e-10, info = candidate)
      expect_true(all(is.finite(after$price)), info = candidate)
      if (candidate == "CF") expect_setequal(solved@tariff_state$corePost, letters[1:2])
      if (candidate == "CS") expect_setequal(solved@tariff_state$leadersPost, letters[1:2])
      after$price
    })
    expect_equal(prices[[1L]], prices[[2L]], tolerance = 1e-7, info = candidate)
  }
})
