rho_test_dt_scaled_log <- function(x, mu, s, nu) {
  stats::dt((x - mu) / s, df = nu, log = TRUE) - log(s)
}

# Sampling tests in this file exercise the compiled Stan models and are kept
# opt-in so the independent density checks remain part of the fast suite.
skip_unless_rho_stan_tests <- function() {
  enabled <- tolower(Sys.getenv("ANTITRUSTBAYES_RUN_STAN_TESTS", "false"))
  if (!identical(enabled, "true")) {
    testthat::skip(paste0(
      "Stan sampling tests are opt-in; set ",
      "ANTITRUSTBAYES_RUN_STAN_TESTS=true to run them."
    ))
  }
}

# Standard bivariate Student-t log density with correlation rho and Stan-style
# "scale" parameters s1, s2 (not SDs).
rho_test_dbvt_joint_log <- function(x1, x2, mu1, mu2, s1, s2, rho, nu) {
  z1 <- (x1 - mu1) / s1
  z2 <- (x2 - mu2) / s2
  q <- (z1^2 - 2 * rho * z1 * z2 + z2^2) / (1 - rho^2)
  log_det <- 2 * log(s1) + 2 * log(s2) + log(1 - rho^2)
  lgamma((nu + 2) / 2) - lgamma(nu / 2) - log(nu * pi) -
    0.5 * log_det - (nu + 2) / 2 * log1p(q / nu)
}

# Factorized marginal-share Student-t (df nu) + conditional-margin Student-t
# (df nu+1), using the exact SD-parameterized formulas implemented in
# logit_supply.stan / ces_supply.stan.
rho_test_factorized_log <- function(x1, x2, mu1, mu2, sigma_share, sigma_margin,
                                     rho, nu) {
  s1 <- sigma_share * sqrt((nu - 2) / nu)
  z_s <- (x1 - mu1) / sigma_share
  slope_cond <- rho * (sigma_margin / sigma_share)
  mu2_cond <- mu2 + slope_cond * (x1 - mu1)
  scale2_cond <- sigma_margin * sqrt(1 - rho^2) * sqrt((nu - 2 + z_s^2) / (nu + 1))
  rho_test_dt_scaled_log(x1, mu1, s1, nu) +
    rho_test_dt_scaled_log(x2, mu2_cond, scale2_cond, nu + 1)
}

make_rho_smoke_data <- function(seed = 123, N = 24) {
  set.seed(seed)
  list(
    prices = runif(N, 1, 5),
    quantities = runif(N, 10, 100),
    margins = runif(N, 0.05, 0.35),
    market = rep("m1", N),
    group = rep("g1", N)
  )
}

test_that("factorized bivariate Student-t matches the joint density across a grid", {
  grid <- expand.grid(
    nu = c(3, 5, 10, 50),
    rho = c(-0.8, -0.3, 0, 0.3, 0.8),
    sigma_share = c(0.3, 0.8, 1.5),
    sigma_margin = c(0.3, 0.8, 1.5),
    err_share = c(-2, -0.5, 0, 0.5, 2),
    err_margin = c(-1, 0, 1)
  )
  set.seed(42)
  grid <- grid[sample(nrow(grid), min(400, nrow(grid))), ]
  mu1 <- 0.7
  mu2 <- -0.3
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    x1 <- mu1 + g$err_share
    x2 <- mu2 + g$err_margin
    s1 <- g$sigma_share * sqrt((g$nu - 2) / g$nu)
    s2 <- g$sigma_margin * sqrt((g$nu - 2) / g$nu)
    joint <- rho_test_dbvt_joint_log(x1, x2, mu1, mu2, s1, s2, g$rho, g$nu)
    fact <- rho_test_factorized_log(
      x1, x2, mu1, mu2, g$sigma_share, g$sigma_margin, g$rho, g$nu
    )
    expect_equal(fact, joint, tolerance = 1e-6,
                 label = paste("row", i, "nu", g$nu, "rho", g$rho))
  }
})

test_that("nu -> Inf recovers the historical Normal conditional", {
  sigma_share <- 0.6
  sigma_margin <- 0.9
  rho <- 0.4
  nu <- 1e5
  z_s <- 1.3
  scale2_cond <- sigma_margin * sqrt(1 - rho^2) * sqrt((nu - 2 + z_s^2) / (nu + 1))
  expect_equal(scale2_cond, sigma_margin * sqrt(1 - rho^2), tolerance = 1e-4)
})

test_that("fixed and estimated df fits compile and produce finite draws", {
  skip_unless_rho_stan_tests()
  d <- make_rho_smoke_data()
  fast_opts <- list(chains = 1, iter = 60, warmup = 30, refresh = 0, seed = 1)

  fit_fixed <- suppressWarnings(bayes_est(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, group = d$group,
    demand_model = "logit", supply_model = "bertrand",
    df = 5,
    stan_options = fast_opts
  ))
  rho_draws <- rstan::extract(fit_fixed@stanfit, "rho")[[1]]
  s0_draws <- rstan::extract(fit_fixed@stanfit, "s0_prob")[[1]]
  alpha_draws <- rstan::extract(fit_fixed@stanfit, "a_market")[[1]]
  expect_true(all(is.finite(rho_draws)))
  expect_true(all(is.finite(s0_draws)))
  expect_true(all(is.finite(alpha_draws)))

  # Reconstruct one generated joint log-likelihood term independently from
  # the same Logit fit used by this smoke test.
  draws_fixed <- rstan::extract(fit_fixed@stanfit)
  i <- 1L
  draw_idx <- 1L
  expected_joint <- rho_test_factorized_log(
    fit_fixed@stan_data$log_shareIn[i],
    fit_fixed@stan_data$log_marginNormInv[i],
    log(as.numeric(draws_fixed$pred_shareIn[draw_idx, i])),
    log(as.numeric(draws_fixed$pred_marginNormInv[draw_idx, i]) + 1e-9),
    as.numeric(draws_fixed$sigma_share[draw_idx]),
    as.numeric(draws_fixed$sigma_margin[draw_idx]),
    as.numeric(draws_fixed$rho[draw_idx]),
    as.numeric(draws_fixed$nu[draw_idx])
  )
  expect_equal(
    as.numeric(draws_fixed$log_lik[draw_idx, i]),
    expected_joint,
    tolerance = 1e-5
  )

  fit_est <- suppressWarnings(bayes_est(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, group = d$group,
    demand_model = "logit", supply_model = "bertrand",
    df = "estimate",
    stan_options = fast_opts
  ))
  nu_draws <- rstan::extract(fit_est@stanfit, "nu")[[1]]
  expect_true(all(is.finite(nu_draws)))
  expect_true(all(nu_draws > 2))
})

test_that("rho is never instantiated when no margins are observed", {
  skip_unless_rho_stan_tests()
  d <- make_rho_smoke_data()
  fast_opts <- list(chains = 1, iter = 60, warmup = 30, refresh = 0, seed = 1)
  fit <- suppressWarnings(bayes_est(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, group = d$group,
    demand_model = "logit", supply_model = "bertrand",
    df = 5,
    margin_observed = rep(FALSE, length(d$prices)),
    stan_options = fast_opts
  ))
  expect_equal(fit@stanfit@par_dims[["rho_raw"]], 0)
  rho_draws <- rstan::extract(fit@stanfit, "rho")[[1]]
  expect_true(all(rho_draws == 0))
})

test_that("log_lik matches an independent Student-t computation for share/margin-only modes", {
  skip_unless_rho_stan_tests()
  d <- make_rho_smoke_data()
  fast_opts <- list(chains = 1, iter = 60, warmup = 30, refresh = 0, seed = 1)
  draw_idx <- 1L

  fit_share <- suppressWarnings(bayes_est(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, group = d$group,
    demand_model = "logit", supply_model = "bertrand",
    df = 5, likelihood = "share",
    stan_options = fast_opts
  ))
  draws <- rstan::extract(fit_share@stanfit)
  log_shareIn <- as.numeric(fit_share@stan_data$log_shareIn)
  sigma_share <- draws$sigma_share[draw_idx]
  nu <- draws$nu[draw_idx]
  pred_shareIn <- draws$pred_shareIn[draw_idx, ]
  scale_marginal <- sigma_share * sqrt((nu - 2) / nu)
  expected_ll <- rho_test_dt_scaled_log(log_shareIn, log(pred_shareIn), scale_marginal, nu)
  actual_ll <- draws$log_lik[draw_idx, ]
  expect_equal(actual_ll, expected_ll, tolerance = 1e-6)
  expect_equal(actual_ll, draws$log_lik_share_marginal[draw_idx, ], tolerance = 1e-8)

  fit_margin <- suppressWarnings(bayes_est(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, group = d$group,
    demand_model = "logit", supply_model = "bertrand",
    df = 5, likelihood = "margin",
    stan_options = fast_opts
  ))
  draws_m <- rstan::extract(fit_margin@stanfit)
  log_marginNormInv <- as.numeric(fit_margin@stan_data$log_marginNormInv)
  sigma_margin <- draws_m$sigma_margin[draw_idx]
  nu_m <- draws_m$nu[draw_idx]
  pred_marginNormInv <- draws_m$pred_marginNormInv[draw_idx, ]
  scale_marginal_m <- sigma_margin * sqrt((nu_m - 2) / nu_m)
  expected_ll_m <- rho_test_dt_scaled_log(
    log_marginNormInv, log(pred_marginNormInv), scale_marginal_m, nu_m
  )
  actual_ll_m <- draws_m$log_lik[draw_idx, ]
  expect_equal(actual_ll_m, expected_ll_m, tolerance = 1e-6)
  expect_equal(actual_ll_m, draws_m$log_lik_margin[draw_idx, ], tolerance = 1e-8)
})

test_that("CES rho and Student-t log_lik agree with independent calculations", {
  skip_unless_rho_stan_tests()
  d <- make_rho_smoke_data(seed = 321, N = 12)
  stan_options <- list(chains = 1, iter = 80, warmup = 40, refresh = 0, seed = 19)

  fit <- suppressWarnings(bayes_est(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, group = d$group,
    demand_model = "ces", supply_model = "bertrand",
    df = 5, use_rho = TRUE, likelihood = "joint",
    stan_options = stan_options
  ))
  draws <- rstan::extract(fit@stanfit)
  draw_idx <- 1L

  # The active CES fit must expose finite structural quantities, including
  # the correlation used by the factored bivariate Student-t likelihood.
  expect_equal(fit@stanfit@par_dims[["rho_raw"]], 1L)
  rho <- as.numeric(draws$rho[draw_idx])
  expect_true(is.finite(rho))
  expect_true(is.finite(as.numeric(draws$s0_prob[draw_idx, 1L])))
  expect_true(is.finite(as.numeric(draws$gamma_structural_market[draw_idx, 1L])))

  # Reconstruct the generated joint log likelihood for one paired observation
  # from the public Stan draws and the independent R Student-t density.
  i <- 1L
  log_share_obs <- fit@stan_data$log_shareIn[i]
  log_margin_obs <- fit@stan_data$log_marginNormInv[i]
  pred_log_q <- log(as.numeric(draws$pred_shareIn[draw_idx, i]))
  pred_log_margin <- log(as.numeric(draws$pred_marginNormInv[draw_idx, i]) + 1e-9)
  sigma_share <- as.numeric(draws$sigma_share[draw_idx])
  sigma_margin <- as.numeric(draws$sigma_margin[draw_idx])
  nu <- as.numeric(draws$nu[draw_idx])

  scale_share <- sigma_share * sqrt((nu - 2) / nu)
  ll_share <- rho_test_dt_scaled_log(
    log_share_obs, pred_log_q, scale_share, nu
  )
  err_share <- log_share_obs - pred_log_q
  z_s <- err_share / sigma_share
  mu_margin_cond <- pred_log_margin + rho *
    (sigma_margin / sigma_share) * err_share
  scale_margin_cond <- sigma_margin * sqrt(1 - rho^2) *
    sqrt((nu - 2 + z_s^2) / (nu + 1))
  ll_margin <- rho_test_dt_scaled_log(
    log_margin_obs, mu_margin_cond, scale_margin_cond, nu + 1
  )
  expected_ll <- ll_share + ll_margin
  actual_ll <- as.numeric(draws$log_lik[draw_idx, i])
  expect_equal(actual_ll, expected_ll, tolerance = 1e-5)

  # With rho disabled and no observed margin equation, Stan uses a zero-length
  # rho_raw parameter and the transformed correlation is exactly zero.
  fit_inactive <- suppressWarnings(bayes_est(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, group = d$group,
    demand_model = "ces", supply_model = "bertrand",
    df = 5, use_rho = FALSE, margin_observed = rep(FALSE, length(d$prices)),
    likelihood = "joint", stan_options = stan_options
  ))
  expect_equal(fit_inactive@stanfit@par_dims[["rho_raw"]], 0)
  inactive_rho <- rstan::extract(fit_inactive@stanfit, "rho")[[1]]
  expect_true(all(inactive_rho == 0))
})

test_that("predicted shares distinguish conditional and unconditional identities", {
  skip_unless_rho_stan_tests()
  d <- make_rho_smoke_data(N = 8)
  fast_opts <- list(chains = 1, iter = 60, warmup = 30, refresh = 0, seed = 7)
  fit <- suppressWarnings(bayes_est(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, group = d$group,
    demand_model = "logit", supply_model = "bertrand",
    df = 1e4, use_rho = FALSE,
    stan_options = fast_opts
  ))
  draws <- rstan::extract(fit@stanfit)
  draw_idx <- 1L
  s0_prob <- draws$s0_prob[draw_idx, ]
  mky <- fit@stan_data$market_year_idx
  a_market <- draws$a_market[draw_idx, ]
  b_market <- draws$b_market[draw_idx, ]
  price_std <- as.numeric(fit@stan_data$price_std)
  market_sign <- fit@stan_data$market_sign
  m <- fit@stan_data$market

  delta <- b_market[m] + market_sign * a_market[m] * price_std
  expected_pred_shareIn <- (s0_prob[mky] / (1 - s0_prob[mky])) * exp(delta)
  wrong_pred_shareIn <- s0_prob[mky] * exp(delta)

  actual_pred_shareIn <- draws$pred_shareIn[draw_idx, ]
  expect_equal(actual_pred_shareIn, expected_pred_shareIn, tolerance = 1e-6)
  expect_true(max(abs(actual_pred_shareIn - wrong_pred_shareIn)) > 1e-3)

  # s0_prob is the outside quantity share, so 1 - s0_prob is the total
  # inside quantity share. The generated unconditional product share must
  # therefore satisfy s_j = (1 - s0) * s_{j|I}.
  pred_shareUncond <- draws$pred_shareUncond[draw_idx, ]
  expect_equal(pred_shareUncond, (1 - s0_prob[mky]) * actual_pred_shareIn,
               tolerance = 1e-8)
})

test_that("share-only likelihood does not use the equilibrium-margin information that identifies s0", {
  skip_unless_rho_stan_tests()
  # likelihood_mode == 2 (share = "share") never touches margin_obs /
  # log_marginNormInv in the Stan model block, so s0's posterior should be
  # (up to MCMC noise) invariant to the observed margins, unlike the joint
  # fit where the margin FOC ties s0 to the markup. This tests the actual
  # mechanism behind the identification caveat directly, rather than
  # inferring it indirectly from a posterior-SD comparison (which is
  # confounded by market-intercept rotation and prior effects and does not
  # reliably move in one direction across data regimes).
  set.seed(99)
  N <- 40
  prices <- runif(N, 1, 5)
  quantities <- runif(N, 10, 100)
  market <- rep("m1", N)
  group <- rep("g1", N)
  margins_a <- runif(N, 0.05, 0.35)
  margins_b <- runif(N, 0.5, 0.9)
  opts <- list(chains = 2, iter = 300, warmup = 150, refresh = 0, seed = 11)

  fit_share_a <- suppressWarnings(bayes_est(
    prices = prices, quantities = quantities, margins = margins_a,
    market = market, group = group,
    demand_model = "logit", supply_model = "bertrand",
    df = 5, likelihood = "share",
    stan_options = opts
  ))
  fit_share_b <- suppressWarnings(bayes_est(
    prices = prices, quantities = quantities, margins = margins_b,
    market = market, group = group,
    demand_model = "logit", supply_model = "bertrand",
    df = 5, likelihood = "share",
    stan_options = opts
  ))
  fit_joint_a <- suppressWarnings(bayes_est(
    prices = prices, quantities = quantities, margins = margins_a,
    market = market, group = group,
    demand_model = "logit", supply_model = "bertrand",
    df = 5, likelihood = "joint",
    stan_options = opts
  ))
  fit_joint_b <- suppressWarnings(bayes_est(
    prices = prices, quantities = quantities, margins = margins_b,
    market = market, group = group,
    demand_model = "logit", supply_model = "bertrand",
    df = 5, likelihood = "joint",
    stan_options = opts
  ))

  s0_share_a <- rstan::extract(fit_share_a@stanfit, "s0_logit")[[1]]
  s0_share_b <- rstan::extract(fit_share_b@stanfit, "s0_logit")[[1]]
  s0_joint_a <- rstan::extract(fit_joint_a@stanfit, "s0_logit")[[1]]
  s0_joint_b <- rstan::extract(fit_joint_b@stanfit, "s0_logit")[[1]]

  share_shift <- abs(mean(s0_share_a) - mean(s0_share_b))
  joint_shift <- abs(mean(s0_joint_a) - mean(s0_joint_b))

  # Share-only s0 barely moves when only the (unused) margins change; the
  # joint fit's s0 moves substantially because the margin FOC pins it down.
  expect_lt(share_shift, joint_shift)
})
