.mic_ns <- function(name) getFromNamespace(name, "antitrustBayes")

.mic_fixture <- function(n_market = 2L, overrides = list()) {
  if (n_market == 2L) {
    market <- rep(c("east", "west"), each = 4L)
    year <- rep(c("y1", "y1", "y2", "y2"), 2L)
    group <- rep(c("A", "B"), 4L)
  } else {
    market <- rep("solo", 4L)
    year <- rep(c("y1", "y1", "y2", "y2"), 1L)
    group <- rep(c("A", "B"), 2L)
  }
  n <- length(market)
  args <- list(
    prices = seq(1.1, 1.1 + 0.1 * (n - 1L), by = 0.1),
    quantities = rep(c(40, 60, 35, 65), length.out = n),
    margins = rep(c(.25, .35, .20, .30), length.out = n),
    market = market, year = year, group = group,
    supply_model = 1L, use_cutoff = FALSE, use_rho = TRUE,
    priors = list(), market_intercept_parameterization = "sampled"
  )
  args[names(overrides)] <- overrides
  data <- do.call(.mic_ns("prepare_stan_data"), args)
  list(data = data, prices = args$prices, quantities = args$quantities,
       margins = args$margins, market = market, year = year, group = group)
}

.mic_fit <- function(prepared, posterior_draws = 2L) {
  intercept_mode <- attr(
    prepared$data, "market_intercept_parameterization_effective",
    exact = TRUE
  )
  if (is.null(intercept_mode)) intercept_mode <- "sampled"
  df <- data.frame(
    prices = prepared$prices, quantities = prepared$quantities,
    margins = prepared$margins, market = factor(prepared$market),
    year = factor(prepared$year), group = factor(prepared$group),
    outside_price = rep(1, length(prepared$prices))
  )
  md <- list(
    quality_parameterization = "collapsed", quality_scope = "market",
    share_likelihood = "gaussian", margin_likelihood = "gaussian",
    likelihood_mode = 1L,
    market_intercept_parameterization = intercept_mode,
    market_intercept_parameterization_effective = intercept_mode,
    quality_scale_parameterization = "legacy",
    quality_scale_parameterization_effective = "legacy",
    structure = "simultaneous", partition_prior = "cutoff",
    posterior_draws = as.integer(posterior_draws),
    group_levels = attr(prepared$data, "group_levels"),
    market_levels = attr(prepared$data, "market_levels")
  )
  new("BayesLogit", spec = NULL, model = NULL, parameters = list(),
      observed = list(), diagnostics = list(), stanfit = NULL,
      data = df, stan_data = prepared$data, supply_model = 1L,
      demand_model = "logit", market_side = "input", price_sd = 1,
      call = quote(identity()), metadata = md)
}

.mic_dense_moments <- function(z, variance, group, tau, mu_b, sigma_b) {
  n_group <- max(group)
  if (n_group == 1L) {
    post_var <- 1 / (1 / sigma_b^2 + sum(1 / variance))
    post_mean <- post_var * sum(z / variance)
    return(list(intercept_mean = mu_b + post_mean,
                intercept_sd = sqrt(post_var), quality_mean = 0,
                quality_sd = 0))
  }
  basis <- qr.Q(qr(contr.helmert(n_group)))
  ## The columns of the orthonormal contrast basis span the centered subspace.
  basis <- basis[, seq_len(n_group - 1L), drop = FALSE]
  design <- cbind(1, basis[group, , drop = FALSE])
  prior_cov <- diag(c(sigma_b^2, rep(tau^2, n_group - 1L)))
  precision <- solve(prior_cov) + crossprod(design, design / variance)
  cov <- solve(precision)
  mean <- cov %*% crossprod(design, z / variance)
  output_map <- rbind(c(1, rep(0, n_group - 1L)), cbind(0, basis))
  out_mean <- as.numeric(output_map %*% mean)
  out_cov <- output_map %*% cov %*% t(output_map)
  list(intercept_mean = mu_b + out_mean[[1L]],
       intercept_sd = sqrt(out_cov[1L, 1L]),
       quality_mean = out_mean[-1L],
       quality_sd = sqrt(diag(out_cov)[-1L]))
}

.mic_exact_k_quality_fixture <- function(n_draw = 80L) {
  prepared <- .mic_fixture(overrides = list(
    market_intercept_parameterization = "collapsed"
  ))
  fit <- .mic_fit(prepared, posterior_draws = n_draw)
  fit@metadata$structure <- "stackelberg"
  fit@metadata$partition_prior <- "ordered_k"
  fit@metadata$k_scope <- "market"
  d <- fit@stan_data
  d$N_partition_unit <- 2L
  d$partition_unit_market <- 1:2
  d$partition_unit_start <- c(1L, 3L)
  d$partition_unit_count <- c(2L, 2L)
  d$partition_component_unit <- c(1L, 1L, 2L, 2L)
  d$partition_component_k <- c(1L, 2L, 1L, 2L)
  d$partition_prior_probability <- c(.55, .45, .60, .40)
  d$partition_prior_logprob <- log(d$partition_prior_probability)
  d$partition_prior_code <- 1L
  d$market_year_idx <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)
  d$N_market_year <- 4L
  d$market <- rep(1:2, each = 4L)
  d$N_market <- 2L
  d$firm_mky_start <- c(1L, 3L, 5L, 7L)
  d$firm_mky_n <- rep(2L, 4L)
  d$firm_mky_rank <- rep(c(1L, 2L), 4L)
  d$obs_to_firm_mky <- seq_len(8L)
  d$mky_to_market <- c(1L, 1L, 2L, 2L)
  d$firm_mky_shareIn <- rep(c(.60, .40), 4L)
  d$price <- prepared$prices
  d$price_sd <- 1
  d$log_marginNormInv <- log(c(.15, .42, .23, .51, .18, .47, .27, .58))
  d$margin_idx <- seq_along(prepared$prices)
  d$N_margin <- length(d$margin_idx)
  d$N_firm_mky <- 8L
  fit@stan_data <- d
  attr(fit@stan_data, "structure") <- "stackelberg"
  attr(fit@stan_data, "partition_prior") <- "ordered_k"
  draws <- list(
    share_mean_integrated = matrix(
      rep(seq(-1.1, -.65, length.out = length(prepared$prices)), n_draw),
      nrow = n_draw, byrow = TRUE
    ),
    sigma_b_strat = rep(.55, n_draw),
    sigma_share = rep(.35, n_draw),
    mu_b_market_vec = rep(c(.10, -.08), length.out = n_draw),
    sigma_b_market_vec = rep(c(.30, .36), length.out = n_draw),
    rho = rep(.35, n_draw),
    sigma_margin = rep(.60, n_draw),
    a_market = matrix(rep(c(.8, .9), each = n_draw), nrow = n_draw),
    s0_prob = matrix(rep(c(.20, .23, .25, .28), each = n_draw), nrow = n_draw),
    posterior_component_prob = matrix(
      rep(c(.55, .45, .60, .40), each = n_draw), nrow = n_draw
    ),
    prior_component_prob = matrix(
      rep(c(.50, .50, .50, .50), each = n_draw), nrow = n_draw
    )
  )
  list(fit = fit, draws = draws, prepared = prepared)
}

test_that("market-intercept API preserves sampled defaults and emits stable codes", {
  sampled <- .mic_fixture()
  expect_equal(sampled$data$market_intercept_parameterization, 1L)
  expect_identical(attr(sampled$data, "market_intercept_parameterization_effective"),
                   "sampled")

  collapsed <- .mic_fixture(overrides = list(
    market_intercept_parameterization = "collapsed"
  ))
  expect_equal(collapsed$data$market_intercept_parameterization, 2L)
  expect_identical(attr(collapsed$data, "market_intercept_parameterization_requested"),
                   "collapsed")
  expect_identical(attr(collapsed$data, "market_intercept_parameterization_effective"),
                   "collapsed")
  expect_equal(collapsed$data$N_market, 2L)
  expect_equal(collapsed$data$N_market_year, 4L)
  incidence <- attr(collapsed$data, "quality_incidence")
  expect_equal(incidence$N_quality, 4L)
  expect_equal(incidence$N_quality_block, 2L)
  expect_equal(incidence$quality_block_n, c(2L, 2L))

  single <- .mic_fixture(1L, list(
    market_intercept_parameterization = "collapsed",
    share_likelihood = "student_t", margin_likelihood = "student_t",
    quality_parameterization = "sampled", quality_scope = "global",
    likelihood_mode = 2L
  ))
  expect_equal(single$data$market_intercept_parameterization, 1L)
  expect_identical(attr(single$data, "market_intercept_parameterization_requested"),
                   "collapsed")
  expect_identical(attr(single$data, "market_intercept_parameterization_effective"),
                   "sampled")
  expect_true(attr(single$data, "market_intercept_parameterization_no_op"))
})

test_that("multi-market collapse validates the reviewed Gaussian coordinate contract", {
  common <- list(market_intercept_parameterization = "collapsed")
  expect_error(.mic_fixture(overrides = c(common, list(share_likelihood = "student_t"))),
               "Gaussian shares and Gaussian margins")
  expect_error(.mic_fixture(overrides = c(common, list(margin_likelihood = "student_t"))),
               "Gaussian shares and Gaussian margins")
  expect_error(.mic_fixture(overrides = c(common, list(quality_parameterization = "sampled"))),
               "quality_parameterization = 'collapsed'")
  expect_error(.mic_fixture(overrides = c(common, list(quality_scope = "global"))),
               "quality_scope = 'market'")
  expect_error(.mic_fixture(overrides = c(common, list(likelihood_mode = 2L))),
               "likelihood = 'joint'")
  expect_error(.mic_fixture(overrides = c(common, list(market_s0_parameterization = "rotated"))),
               "market_s0_parameterization = 'legacy'")
  expect_error(.mic_fixture(overrides = c(common, list(
    use_market_resid_scale = TRUE,
    quality_scale_parameterization = "total_fraction"
  ))), "quality_scale_parameterization = 'legacy'")

  market_scale <- .mic_fixture(overrides = c(common, list(
    use_market_resid_scale = TRUE,
    quality_scale_parameterization = "legacy"
  )))
  expect_equal(market_scale$data$market_intercept_parameterization, 2L)
  expect_equal(market_scale$data$use_market_resid_scale, 1L)
  total_fraction <- .mic_fixture(overrides = c(common, list(
    quality_scale_parameterization = "total_fraction"
  )))
  expect_equal(total_fraction$data$market_intercept_parameterization, 2L)
  expect_equal(total_fraction$data$quality_scale_parameterization, 2L)
})

test_that("collapsed intercept data accepts partially or wholly missing margins", {
  partial <- .mic_fixture(overrides = list(
    market_intercept_parameterization = "collapsed",
    margin_observed = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE)
  ))
  expect_equal(partial$data$market_intercept_parameterization, 2L)
  expect_equal(partial$data$N_margin, 3L)

  absent <- .mic_fixture(overrides = list(
    market_intercept_parameterization = "collapsed",
    margin_observed = rep(FALSE, 8L)
  ))
  expect_equal(absent$data$market_intercept_parameterization, 2L)
  expect_equal(absent$data$N_margin, 0L)
  expect_length(absent$data$margin_idx, 0L)
})

test_that("all four production kernels expose the collapsed-intercept contract", {
  kernels <- c("logit_supply.stan", "logit_stackelberg.stan",
               "ces_supply.stan", "ces_stackelberg.stan")
  for (kernel in kernels) {
    path <- testthat::test_path("..", "..", "inst", "stan", kernel)
    source <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_match(source, "market_intercept_parameterization", fixed = TRUE,
                 info = kernel)
    expect_match(source, "centered_intercept_gaussian_from_stats", fixed = TRUE,
                 info = kernel)
    expect_match(source, "share_mean_integrated", fixed = TRUE, info = kernel)
    expect_match(source, "log_lik_intercept_quality_integrated_market",
                 fixed = TRUE, info = kernel)
  }
})

test_that("joint intercept and centered quality moments match a dense Gaussian oracle", {
  z <- c(.4, -.15, .7, .05, -.3, .55, .12)
  variance <- c(.3, .6, .4, .25, .7, .5, .35)
  group <- c(1L, 2L, 3L, 1L, 2L, 3L, 1L)
  tau <- .42
  mu_b <- -.08
  sigma_b <- .31
  got <- .mic_ns(".posterior_quality_market_intercept_conditional")(
    z, variance, group, tau, mu_b, sigma_b
  )
  expected <- .mic_dense_moments(z, variance, group, tau, mu_b, sigma_b)
  expect_equal(got$intercept_mean, expected$intercept_mean, tolerance = 1e-10)
  expect_equal(got$intercept_sd, expected$intercept_sd, tolerance = 1e-10)
  expect_equal(got$quality_mean, expected$quality_mean, tolerance = 1e-10)
  expect_equal(got$quality_sd, expected$quality_sd, tolerance = 1e-10)
  expect_equal(sum(got$quality_mean), 0, tolerance = 1e-12)

  one <- .mic_ns(".posterior_quality_market_intercept_conditional")(
    z = c(.2, .5, -.1), variance = c(.2, .4, .3),
    group = rep(1L, 3L), tau = 4, mu_b = .15, sigma_b = .5
  )
  one_expected <- .mic_dense_moments(
    c(.2, .5, -.1), c(.2, .4, .3), rep(1L, 3L), 0, .15, .5
  )
  expect_equal(one$intercept_mean, one_expected$intercept_mean, tolerance = 1e-12)
  expect_equal(one$intercept_sd, one_expected$intercept_sd, tolerance = 1e-12)
  expect_equal(one$quality_mean, 0)
  expect_equal(one$quality_sd, 0)
  expect_equal(one$quality_draw, 0)
})

test_that("policy plans bind reconstructed intercepts and reuse their state in scores", {
  prepared <- .mic_fixture(overrides = list(
    market_intercept_parameterization = "collapsed",
    margin_observed = rep(FALSE, 8L), use_rho = FALSE
  ))
  fit <- .mic_fit(prepared)
  n <- length(prepared$prices)
  draws <- list(
    a_market = matrix(c(.8, .9, .85, .95), nrow = 2L, ncol = 2L),
    s0_prob = matrix(c(.20, .22, .24, .26), nrow = 2L, ncol = 2L),
    share_mean_integrated = matrix(rep(seq(-1.2, -.7, length.out = n), 2L),
                                   nrow = 2L, byrow = TRUE),
    sigma_b_strat = c(.45, .50), sigma_share = c(.22, .24),
    mu_b_market_vec = c(.10, -.08), sigma_b_market_vec = c(.35, .38)
  )
  plan <- .mic_ns(".lq_make_policy_plan")(fit, draws, n_draws = 2L, seed = 19L)
  expect_equal(plan$fit_signature$fit_signature_version, 3L)
  expect_identical(colnames(plan$market_intercept_draws), c("east", "west"))
  expect_identical(rownames(plan$market_intercept_draws), as.character(plan$draw_indices))
  expect_silent(.mic_ns(".lq_validate_plan")(fit, plan, draws = draws))

  tampered <- plan
  tampered$market_intercept_draws[1L, 1L] <-
    tampered$market_intercept_draws[1L, 1L] + .1
  expect_error(.mic_ns(".lq_validate_plan")(fit, tampered, draws = draws),
               "latent states or selected posterior rows were changed")

  row <- match(plan$draw_indices[[1L]], seq_len(nrow(draws$share_mean_integrated)))
  unit <- list(idx = 1:4, market_index = 1L)
  market_object <- list(prices = prepared$prices[1:4],
                        quantities = prepared$quantities[1:4])
  score <- .mic_ns(".lq_reference_score")(
    fit, draws, row, unit, plan, market_object, list(s0 = .2)
  )$score
  quality_index <- attr(prepared$data, "quality_incidence")$obs_to_quality[1:4]
  plan_row <- match(row, plan$draw_indices)
  expected <- draws$share_mean_integrated[row, 1:4] +
    plan$market_intercept_draws[plan_row, 1L] - draws$mu_b_market_vec[[row]] +
    plan$quality_draws[plan_row, quality_index]
  expect_equal(score, as.numeric(expected), tolerance = 1e-12)
})

test_that("exact-k joint quality and intercept recovery share sampled market states", {
  fixture <- .mic_exact_k_quality_fixture()
  fit <- fixture$fit
  draws <- fixture$draws
  n_draw <- nrow(draws$posterior_component_prob)
  seed <- 391L

  ## The component probabilities are normalized market by market, with a
  ## unique posterior mode at k = 1 in each market.
  expect_equal(rowSums(draws$posterior_component_prob[, 1:2, drop = FALSE]),
               rep(1, n_draw))
  expect_equal(rowSums(draws$posterior_component_prob[, 3:4, drop = FALSE]),
               rep(1, n_draw))
  k_summary <- posterior_k(fit, draws = draws)
  expect_equal(as.numeric(tapply(k_summary$posterior, k_summary$market, sum)),
               c(1, 1), tolerance = 1e-12)
  parsed <- .mic_ns(".posterior_quality_draws")(fit, draws, "auto")
  joint <- .mic_ns(".posterior_quality_joint_draws")(
    fit, parsed, draws, "auto", draw_indices = seq_len(n_draw), seed = seed + 1L
  )
  expect_equal(dim(joint$realized_k), c(n_draw, 2L))
  expect_true(all(joint$realized_k %in% c(1L, 2L)))
  ## A draw of k = 2 despite P(k = 1) > P(k = 2) establishes that role
  ## reconstruction samples the posterior mixture instead of substituting
  ## its MAP state.
  expect_true(any(joint$realized_k == 2L))
  expect_true(all(vapply(seq_len(2L), function(m) {
    length(unique(joint$realized_k[, m])) == 2L
  }, logical(1))))

  realized_prediction <- .mic_ns(".stack_pure_foc_margin_prediction")(
    fit, parsed, draws, seq_len(n_draw), joint$realized_k
  )
  expect_equal(joint$prediction, realized_prediction, tolerance = 1e-12)
  map_prediction <- .mic_ns(".stack_pure_foc_margin_prediction")(
    fit, parsed, draws, seq_len(n_draw), matrix(1L, nrow = n_draw, ncol = 2L)
  )
  expect_gt(max(abs(joint$prediction - map_prediction)), 1e-6)

  recovered <- posterior_quality(fit, draws = draws, seed = seed)
  expect_identical(recovered$realized_k, joint$realized_k)
  expect_identical(dim(recovered$market_intercept_draws), c(n_draw, 2L))
  expect_true(all(is.finite(recovered$market_intercept_draws)))
  expect_equal(unname(rowSums(recovered$draws[, 1:2, drop = FALSE])),
               rep(0, n_draw), tolerance = 1e-10)
  expect_equal(unname(rowSums(recovered$draws[, 3:4, drop = FALSE])),
               rep(0, n_draw), tolerance = 1e-10)

  ## Check both recovered blocks against the conditional Gaussian moments
  ## built from the same realized-k margin prediction. Forcing that sampled
  ## k through the cache reproduces exactly the same intercept and quality
  ## draws, so the two latent states cannot silently use different k values.
  row <- which(joint$realized_k[, 1L] == 2L)[[1L]]
  market <- 1L
  obs <- which(fit@stan_data$market == market)
  local_units <- which(parsed$quality_lookup$market_index == market)
  margin_obs <- joint$idx[joint$idx %in% obs]
  local_margin <- match(margin_obs, obs)
  sigma_share <- parsed$sigma[[row]]
  sigma_margin <- joint$sigma_by_market[row, market]
  rho <- joint$rho[[row]]
  z <- parsed$y[obs] - parsed$share_mean[row, obs]
  variance <- rep(sigma_share^2, length(obs))
  slope <- rho * sigma_share / sigma_margin
  residual_margin <- joint$y[margin_obs] -
    log(joint$prediction[row, margin_obs] + 1e-9)
  z[local_margin] <- z[local_margin] - slope * residual_margin
  variance[local_margin] <- sigma_share^2 * (1 - rho^2)
  expected <- .mic_ns(".posterior_quality_market_intercept_conditional")(
    z, variance, match(parsed$quality[obs], local_units),
    parsed$tau[[row]], parsed$mu_b_market[[row]], parsed$sigma_b_market[[row]]
  )
  expect_equal(recovered$market_intercept_conditional_means[row, market],
               expected$intercept_mean, tolerance = 1e-11)
  expect_equal(unname(recovered$conditional_means[row, local_units]),
               unname(expected$quality_mean), tolerance = 1e-11)

  k_cache <- stats::setNames(lapply(seq_len(n_draw), function(i) {
    as.integer(recovered$realized_k[i, ])
  }), as.character(seq_len(n_draw)))
  replay <- posterior_quality(fit, draws = draws, seed = seed, k_cache = k_cache)
  expect_identical(replay$realized_k, recovered$realized_k)
  expect_identical(replay$draws, recovered$draws)
  expect_identical(replay$market_intercept_draws,
                   recovered$market_intercept_draws)
})

test_that("intercept collapse preserves hypothesis identity but enters policy provenance", {
  prepared <- .mic_fixture(overrides = list(
    use_rho = FALSE, market_intercept_parameterization = "sampled"
  ))
  sampled <- .mic_fit(prepared, posterior_draws = 1L)
  collapsed <- sampled
  collapsed@stan_data$market_intercept_parameterization <- 2L
  attr(collapsed@stan_data, "market_intercept_parameterization") <- "collapsed"
  attr(collapsed@stan_data, "market_intercept_parameterization_requested") <- "collapsed"
  attr(collapsed@stan_data, "market_intercept_parameterization_effective") <- "collapsed"
  attr(collapsed@stan_data, "market_intercept_parameterization_no_op") <- FALSE
  collapsed@metadata$market_intercept_parameterization <- "collapsed"
  collapsed@metadata$market_intercept_parameterization_requested <- "collapsed"
  collapsed@metadata$market_intercept_parameterization_effective <- "collapsed"
  expect_identical(.mic_ns(".market_intercept_parameterization_from_fit")(sampled),
                   "sampled")
  expect_identical(.mic_ns(".market_intercept_parameterization_from_fit")(collapsed),
                   "collapsed")

  ## The collapse integrates a Gaussian nuisance effect without defining a
  ## second economic hypothesis for model-universe comparisons.
  expect_identical(.mic_ns(".mu_fingerprint")(sampled),
                   .mic_ns(".mu_fingerprint")(collapsed))

  common <- list(
    a_market = matrix(c(.8, .9), nrow = 1L),
    s0_prob = matrix(c(.2, .22, .24, .26), nrow = 1L),
    sigma_b_strat = .45, sigma_share = .25
  )
  sampled_draws <- c(common, list(
    share_mean = matrix(seq(-1, -.3, length.out = 8L), nrow = 1L)
  ))
  collapsed_draws <- c(common, list(
    share_mean_integrated = matrix(seq(-1, -.3, length.out = 8L), nrow = 1L),
    mu_b_market_vec = .1, sigma_b_market_vec = .3
  ))
  expect_false(identical(
    .mic_ns(".lq_fit_fingerprint")(sampled, sampled_draws, version = 3L),
    .mic_ns(".lq_fit_fingerprint")(collapsed, collapsed_draws, version = 3L)
  ))
  plan <- .mic_ns(".lq_make_policy_plan")(
    collapsed, collapsed_draws, n_draws = 1L, seed = 73L
  )
  expect_identical(plan$metadata$market_intercept_parameterization, "collapsed")
  expect_error(
    .mic_ns(".lq_validate_plan")(sampled, plan, draws = sampled_draws),
    "different fit or group map"
  )
})
