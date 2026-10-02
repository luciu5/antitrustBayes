.ppc_ns <- function(name) getFromNamespace(name, "antitrustBayes")

test_that("Gaussian rho replication matches the joint Stan residual law", {
  share_mean <- matrix(c(-1, -.7, -.2), nrow = 1L)
  margin_mean <- matrix(c(-2, -1.5, -1), nrow = 1L)
  share_scale <- matrix(c(.2, .3, .4), nrow = 1L)
  margin_scale <- matrix(c(.5, .6, .7), nrow = 1L)
  seed <- 501L
  actual <- .ppc_ns(".posterior_predict_residual_draws")(
    share_mean, margin_mean, share_scale, margin_scale,
    "gaussian", "gaussian", nu = 5, rho = .35, rho_active = TRUE,
    likelihood_mode = 1L, margin_mask = c(TRUE, FALSE, TRUE), seed = seed
  )

  set.seed(seed)
  err_share <- stats::rnorm(3L, sd = as.numeric(share_scale))
  expected_share <- as.numeric(share_mean) + err_share
  expected_margin_mean <- as.numeric(margin_mean) + .35 *
    as.numeric(margin_scale / share_scale) * err_share
  expected_margin <- stats::rnorm(
    3L, mean = expected_margin_mean,
    sd = as.numeric(margin_scale) * sqrt(1 - .35^2)
  )
  expect_equal(as.numeric(actual$log_shareIn), expected_share)
  expect_equal(as.numeric(actual$log_marginNormInv), expected_margin)
})

test_that("Student-t rho replication matches the conditional t law", {
  share_mean <- matrix(c(-1, -.7), nrow = 1L)
  margin_mean <- matrix(c(-2, -1.5), nrow = 1L)
  share_scale <- matrix(c(.2, .3), nrow = 1L)
  margin_scale <- matrix(c(.5, .6), nrow = 1L)
  nu <- 7
  rho <- -.25
  seed <- 812L
  actual <- .ppc_ns(".posterior_predict_residual_draws")(
    share_mean, margin_mean, share_scale, margin_scale,
    "student_t", "student_t", nu = nu, rho = rho, rho_active = TRUE,
    likelihood_mode = 1L, margin_mask = c(TRUE, FALSE), seed = seed
  )

  set.seed(seed)
  err_share <- stats::rt(2L, df = nu) * as.numeric(share_scale) *
    sqrt((nu - 2) / nu)
  expected_share <- as.numeric(share_mean) + err_share
  z_share <- err_share / as.numeric(share_scale)
  expected_mean <- as.numeric(margin_mean) + rho *
    as.numeric(margin_scale / share_scale) * err_share
  conditional_scale <- as.numeric(margin_scale) * sqrt(1 - rho^2) *
    sqrt((nu - 2 + z_share^2) / (nu + 1))
  expected_margin <- expected_mean + stats::rt(2L, df = nu + 1) *
    conditional_scale
  expect_equal(as.numeric(actual$log_shareIn), expected_share)
  expect_equal(as.numeric(actual$log_marginNormInv), expected_margin)
})

test_that("observed shares condition Gaussian margins at rho zero and nonzero", {
  residual <- .ppc_ns(".posterior_predict_residual_draws")
  share_mean <- matrix(c(-1, -.8, -.6), nrow = 1L)
  margin_mean <- matrix(c(-2, -1.8, -1.6), nrow = 1L)
  share_scale <- matrix(c(.2, .3, .4), nrow = 1L)
  margin_scale <- matrix(c(.5, .6, .7), nrow = 1L)
  observed <- c(-.9, -.95, -.45)
  mask <- c(TRUE, FALSE, TRUE)
  for (rho in c(0, .45)) {
    set.seed(951L)
    before <- .Random.seed
    actual <- residual(
      share_mean, margin_mean, share_scale, margin_scale,
      "gaussian", "gaussian", nu = 5, rho = rho, rho_active = TRUE,
      likelihood_mode = 1L, margin_mask = mask, seed = 72L,
      share_conditioning = "observed", observed_share = observed
    )
    expect_identical(.Random.seed, before)
    set.seed(72L)
    expected <- stats::rnorm(
      length(observed),
      mean = as.numeric(margin_mean) + rho *
        as.numeric(margin_scale / share_scale) *
        (observed - as.numeric(share_mean)),
      sd = as.numeric(margin_scale) * sqrt(1 - rho^2)
    )
    expect_equal(as.numeric(actual$log_shareIn), observed)
    expect_equal(as.numeric(actual$log_marginNormInv), expected)
    expect_true(is.finite(actual$log_marginNormInv[1L, 2L]))
  }
})

test_that("observed shares condition Student-t margins including rho zero", {
  residual <- .ppc_ns(".posterior_predict_residual_draws")
  share_mean <- matrix(c(-1, -.8), nrow = 1L)
  margin_mean <- matrix(c(-2, -1.8), nrow = 1L)
  share_scale <- matrix(c(.2, .3), nrow = 1L)
  margin_scale <- matrix(c(.5, .6), nrow = 1L)
  observed <- c(-.85, -1.05)
  nu <- 7
  for (rho in c(0, -.35)) {
    actual <- residual(
      share_mean, margin_mean, share_scale, margin_scale,
      "student_t", "student_t", nu = nu, rho = rho,
      rho_active = TRUE, likelihood_mode = 1L,
      margin_mask = c(TRUE, FALSE), seed = 712L,
      share_conditioning = "observed", observed_share = observed
    )
    error <- observed - as.numeric(share_mean)
    expected_mean <- as.numeric(margin_mean) + rho *
      as.numeric(margin_scale / share_scale) * error
    expected_scale <- as.numeric(margin_scale) * sqrt(1 - rho^2) *
      sqrt((nu - 2 + (error / as.numeric(share_scale))^2) / (nu + 1))
    set.seed(712L)
    expected <- expected_mean + stats::rt(2L, df = nu + 1) * expected_scale
    expect_equal(as.numeric(actual$log_shareIn), observed)
    expect_equal(as.numeric(actual$log_marginNormInv), expected)
  }
  inactive <- residual(
    share_mean, margin_mean, share_scale, margin_scale,
    "student_t", "student_t", nu = nu, rho = 0,
    rho_active = FALSE, likelihood_mode = 1L,
    margin_mask = c(TRUE, FALSE), seed = 712L,
    share_conditioning = "observed", observed_share = observed
  )
  set.seed(712L)
  independent <- as.numeric(margin_mean) + stats::rt(2L, df = nu) *
    as.numeric(margin_scale) * sqrt((nu - 2) / nu)
  expect_equal(as.numeric(inactive$log_marginNormInv), independent)
})

test_that("default residual replication is unchanged and holdout mask is not conditioning data", {
  residual <- .ppc_ns(".posterior_predict_residual_draws")
  args <- list(
    share_mean = matrix(c(-1, -.8), nrow = 1L),
    margin_mean = matrix(c(-2, -1.8), nrow = 1L),
    share_scale = matrix(.2, nrow = 1L, ncol = 2L),
    margin_scale = matrix(.5, nrow = 1L, ncol = 2L),
    share_family = "gaussian", margin_family = "gaussian",
    nu = 5, rho = .3, rho_active = TRUE, likelihood_mode = 1L,
    margin_mask = c(TRUE, FALSE), seed = 141L
  )
  expect_identical(do.call(residual, args), do.call(residual, c(args,
    list(share_conditioning = "replicate"))))
  conditional <- c(args, list(share_conditioning = "observed",
                              observed_share = c(-.9, -.75)))
  hidden <- do.call(residual, conditional)
  conditional$margin_mask <- c(TRUE, TRUE)
  expect_identical(do.call(residual, conditional), hidden)
  conditional$likelihood_mode <- 2L
  expect_error(do.call(residual, conditional), "requires the joint likelihood")
})

test_that("conditional prediction marks fixed shares and never reads held-out margins", {
  prices <- c(1.1, 1.2, 1.3, 1.4)
  quantities <- c(40, 30, 35, 25)
  margins <- c(.3, NA_real_, .25, .4)
  mask <- c(TRUE, FALSE, TRUE, TRUE)
  stan_data <- prepare_stan_data(
    prices, quantities, margins, rep("m1", 4L),
    rep(c("y1", "y2"), each = 2L), rep(c("a", "b"), 2L),
    supply_model = 1L, use_rho = TRUE, margin_observed = mask,
    likelihood_mode = 1L, quality_parameterization = "sampled"
  )
  fit <- methods::new(
    "BayesLogit", spec = NULL, model = NULL, parameters = list(),
    observed = list(), diagnostics = list(), stanfit = NULL,
    data = data.frame(prices = prices, quantities = quantities, margins = margins,
      market = rep("m1", 4L), year = rep(c("y1", "y2"), each = 2L),
      group = rep(c("a", "b"), 2L)),
    stan_data = stan_data, supply_model = 1L, demand_model = "logit",
    market_side = "input", price_sd = 1, call = quote(identity()),
    metadata = list(
      quality_parameterization = "sampled", quality_scope = "market",
      share_likelihood = "gaussian", margin_likelihood = "gaussian",
      likelihood_mode = 1L, use_rho = TRUE,
      market_intercept_parameterization = "sampled",
      market_intercept_parameterization_effective = "sampled",
      structure = "simultaneous", partition_prior = "cutoff",
      posterior_draws = 2L
    )
  )
  draws <- list(
    share_mean = matrix(rep(c(-1, -.8, -.9, -.7), 2L), nrow = 2L,
                        byrow = TRUE),
    pred_marginNormInv = matrix(c(2, 2.1, 2.2, 2.3,
                                  2.4, 2.5, 2.6, 2.7), nrow = 2L,
                                byrow = TRUE),
    sigma_share = c(.2, .25), sigma_margin = c(.3, .35), rho = c(.4, -.2)
  )
  plan <- list(draw_indices = 1:2,
               fit_signature = list(latent_state_fingerprint = "fixture"))
  predict_from_draws <- .ppc_ns(".posterior_predict_from_draws")
  actual <- predict_from_draws(fit, draws, plan, seed = 92L,
                               share_conditioning = "observed")
  expect_identical(actual$metadata$share_conditioning, "observed")
  expect_identical(actual$metadata$share_replicated, FALSE)
  expect_equal(actual$replicated$log_shareIn,
    matrix(rep(stan_data$log_shareIn, each = 2L), nrow = 2L,
           dimnames = dimnames(actual$replicated$log_shareIn)))
  expect_identical(unname(actual$replicated$log_shareIn[1L, ]),
                   as.numeric(stan_data$log_shareIn))
  expect_identical(actual$observed$margin_observed, mask)
  expect_true(all(is.finite(actual$replicated$log_marginNormInv)))
  altered <- fit
  altered@stan_data$log_marginNormInv[2L] <- 1e6
  hidden <- predict_from_draws(altered, draws, plan, seed = 92L,
                               share_conditioning = "observed")
  expect_identical(hidden$replicated, actual$replicated)
  default <- predict_from_draws(fit, draws, plan, seed = 92L)
  explicit <- predict_from_draws(fit, draws, plan, seed = 92L,
                                 share_conditioning = "replicate")
  expect_identical(default, explicit)
  fit@stan_data$likelihood_mode <- 2L
  expect_error(predict_from_draws(fit, draws, plan, seed = 92L,
    share_conditioning = "observed"), "requires the joint likelihood")
})

test_that("predictive masks keep off-likelihood outcomes NA and impute hidden margins", {
  share_mean <- matrix(c(-1, -.7, -.2), nrow = 1L)
  margin_mean <- matrix(c(-2, -1.5, -1), nrow = 1L)
  share_scale <- matrix(.2, nrow = 1L, ncol = 3L)
  margin_scale <- matrix(.5, nrow = 1L, ncol = 3L)
  mask <- c(TRUE, FALSE, TRUE)
  replicate <- .ppc_ns(".posterior_predict_residual_draws")

  joint <- replicate(share_mean, margin_mean, share_scale, margin_scale,
    "gaussian", "gaussian", nu = 5, rho = 0, rho_active = FALSE,
    likelihood_mode = 1L, margin_mask = mask, seed = 1L)
  margin_only <- replicate(share_mean, margin_mean, share_scale, margin_scale,
    "gaussian", "gaussian", nu = 5, rho = 0, rho_active = FALSE,
    likelihood_mode = 3L, margin_mask = mask, seed = 1L)
  share_only <- replicate(share_mean, margin_mean, share_scale, margin_scale,
    "gaussian", "gaussian", nu = 5, rho = 0, rho_active = FALSE,
    likelihood_mode = 2L, margin_mask = mask, seed = 1L)

  expect_true(all(is.finite(joint$log_marginNormInv)))
  expect_true(all(is.finite(margin_only$log_marginNormInv)))
  expect_true(all(is.na(margin_only$log_shareIn)))
  expect_true(all(is.finite(share_only$log_shareIn)))
  expect_true(all(is.na(share_only$log_marginNormInv)))
})

test_that("replicated Gaussian residuals recover the fitted scale and rho", {
  n_draw <- 40000L
  rho <- .4
  share_scale <- matrix(.3, nrow = n_draw, ncol = 2L)
  margin_scale <- matrix(.7, nrow = n_draw, ncol = 2L)
  z <- .ppc_ns(".posterior_predict_residual_draws")(
    matrix(0, n_draw, 2L), matrix(0, n_draw, 2L),
    share_scale, margin_scale, "gaussian", "gaussian", nu = 5,
    rho = rep(rho, n_draw), rho_active = TRUE, likelihood_mode = 1L,
    margin_mask = c(TRUE, TRUE), seed = 991L
  )
  share <- z$log_shareIn[, 1L]
  margin <- z$log_marginNormInv[, 1L]
  share_var_se <- .3^2 * sqrt(2 / (n_draw - 1))
  margin_var_se <- .7^2 * sqrt(2 / (n_draw - 1))
  correlation_se <- (1 - rho^2) / sqrt(n_draw - 1)
  expect_lte(abs(stats::var(share) - .3^2), 4 * share_var_se)
  expect_lte(abs(stats::var(margin) - .7^2), 4 * margin_var_se)
  expect_lte(abs(stats::cor(share, margin) - rho), 4 * correlation_se)
})

test_that("collapsed quality prediction reuses plan state without changing the fit", {
  prices <- c(1.1, 1.2, 1.3, 1.4)
  quantities <- c(40, 30, 35, 25)
  margins <- c(.3, .35, .25, .4)
  market <- rep("m1", 4L)
  year <- rep(c("y1", "y1", "y2", "y2"), 1L)
  group <- rep(c("a", "b"), 2L)
  stan_data <- prepare_stan_data(
    prices, quantities, margins, market, year, group, supply_model = 1L,
    use_cutoff = FALSE, use_rho = FALSE, likelihood_mode = 2L,
    quality_parameterization = "collapsed", quality_scope = "market"
  )
  fit <- methods::new(
    "BayesLogit", spec = NULL, model = NULL, parameters = list(),
    observed = list(), diagnostics = list(), stanfit = NULL,
    data = data.frame(prices = prices, quantities = quantities, margins = margins,
      market = market, year = year, group = group),
    stan_data = stan_data, supply_model = 1L, demand_model = "logit",
    market_side = "input", price_sd = 1, call = quote(identity()),
    metadata = list(
      quality_parameterization = "collapsed", quality_scope = "market",
      share_likelihood = "gaussian", margin_likelihood = "gaussian",
      likelihood_mode = 2L, market_intercept_parameterization = "sampled",
      market_intercept_parameterization_effective = "sampled",
      structure = "simultaneous", partition_prior = "cutoff",
      posterior_draws = 3L,
      market_levels = c("m1"), group_levels = c("a", "b")
    )
  )
  draws <- list(
    a_market = matrix(c(1, 1.1, .9), nrow = 3L, ncol = 1L),
    s0_prob = matrix(c(.2, .21, .22, .23, .24, .25), nrow = 3L),
    share_mean = matrix(rep(c(-1, -.8, -.9, -.7), 3L), nrow = 3L,
                        byrow = TRUE),
    sigma_b_strat = c(.4, .45, .5), sigma_share = c(.2, .22, .24),
    sigma_margin = c(.3, .32, .34)
  )
  plan <- .ppc_ns(".lq_make_policy_plan")(
    fit, draws, n_draws = 2L, seed = 72L, role_mapping = "mixture"
  )
  before <- fit
  ppc <- .ppc_ns(".posterior_predict_from_draws")(
    fit, draws, plan, seed = 81L
  )
  incidence <- .ppc_ns(".lq_quality_incidence")(fit)
  data_info <- .ppc_ns(".posterior_predict_mask_and_mode")(fit)
  means <- .ppc_ns(".posterior_predict_mean_state")(
    fit, draws, plan, data_info
  )
  expected_mean <- draws$share_mean[plan$draw_indices, , drop = FALSE] +
    plan$quality_draws[, incidence$obs_to_quality, drop = FALSE]
  expect_equal(means$share_mean, expected_mean)
  expect_equal(ppc$draw_plan$fit_signature$latent_state_fingerprint,
               plan$fit_signature$latent_state_fingerprint)
  expect_identical(ppc$draws$state_plan_fingerprint,
    rep(plan$fit_signature$latent_state_fingerprint, nrow(ppc$draws)))
  expect_true(all(is.na(ppc$replicated$log_marginNormInv)))
  for (idx in split(seq_len(4L), stan_data$market_year_idx)) {
    expect_equal(unname(rowSums(ppc$replicated$share_rep_normalized[, idx, drop = FALSE])),
                 rep(1, 2L), tolerance = 1e-12)
  }
  expect_identical(fit, before)
  expect_error(posterior_predict(fit, share_conditioning = "observed"),
               "requires the joint likelihood")
})

test_that("exact-k FOC means accept explicit market-year component states", {
  n <- 4L
  data <- data.frame(
    prices = rep(1, n), quantities = rep(1, n), margins = rep(.3, n),
    product_id = paste0("p", seq_len(n)), market = rep("m1", n),
    year = rep(c("y1", "y2"), each = 2L), group = letters[seq_len(n)]
  )
  d <- list(
    N = n, log_shareIn = rep(log(.25), n), market = rep(1L, n),
    N_market = 1L, market_year_idx = rep(1:2, each = 2L),
    N_market_year = 2L, mky_to_market = c(1L, 1L),
    firm_mky_start = c(1L, 3L), firm_mky_n = c(2L, 2L),
    firm_mky_rank = rep(c(1L, 2L), 2L), obs_to_firm_mky = 1:4,
    firm_mky_shareIn = rep(.5, n), price = rep(1, n), price_sd = 1
  )
  fit <- methods::new(
    "BayesLogit", spec = NULL, model = NULL, parameters = list(),
    observed = list(), diagnostics = list(), stanfit = NULL, data = data,
    stan_data = d, supply_model = 1L, demand_model = "logit",
    market_side = "input", price_sd = 1, call = quote(identity()),
    metadata = list(structure = "stackelberg")
  )
  draws <- list(
    a_market = matrix(5, nrow = 1L, ncol = 1L),
    s0_prob = matrix(c(.2, .2), nrow = 1L, ncol = 2L)
  )
  parsed <- list(n_draw = 1L, y = d$log_shareIn)
  predict <- .ppc_ns(".stack_pure_foc_margin_prediction")
  market_state <- predict(fit, parsed, draws, 1L, matrix(0L, 1L, 1L))
  expanded_state <- predict(fit, parsed, draws, 1L, matrix(c(0L, 0L), 1L),
                            realized_k_scope = "market_year")
  expect_equal(expanded_state, market_state)

  year_state <- predict(fit, parsed, draws, 1L, matrix(c(0L, 1L), 1L),
                        realized_k_scope = "market_year")
  expect_equal(year_state[, 1:2], market_state[, 1:2])
  expect_gt(max(abs(year_state[, 3:4] - market_state[, 3:4])), 1e-8)
  ids <- .ppc_ns(".posterior_predict_observation_ids")(
    fit, list(n = n, market = d$market), d$market_year_idx
  )
  expect_identical(ids$product_id, data$product_id)
  expect_identical(ids$year_id, data$year)
})

test_that("exact-k mixture realizations recover known conditional role probabilities", {
  component <- list(k = 0:3, probability = c(.1, .2, .3, .4), cutoff = NULL)
  n <- 20000L
  set.seed(405L)
  realized <- replicate(n, .ppc_ns(".stack_bridge_select_k")(
    component, params = list(), sizes = c(4, 3, 2),
    role_mapping = "mixture", fixed_k = NULL, role_offset = 0L
  ))
  frequency <- tabulate(realized + 1L, nbins = 4L) / n
  standard_error <- sqrt(component$probability *
                           (1 - component$probability) / n)
  expect_lte(max(abs(frequency - component$probability) /
                   standard_error), 4)
})
