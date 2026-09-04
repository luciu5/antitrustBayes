# Tests for hierarchical outside-share (s0) pooling, both the market-year pool
# (use_s0_market_year_pool: one s0 per market-year) and the market-level pool
# (use_s0_market_pool: one s0 per market, shared across the market's periods).
# These paths previously had no coverage; the market-level pool was added
# because market-year pooling is weakly identified (low BFMI) when markets span
# several periods with sparse per-cell information.

# Small multi-market, multi-period panel with single-product firms.
make_panel <- function(seed = 7, n_markets = 3, n_years = 3, n_firms = 4) {
  set.seed(seed)
  grid <- expand.grid(
    m = seq_len(n_markets), y = seq_len(n_years), f = seq_len(n_firms)
  )
  n <- nrow(grid)
  prices <- pmax(0.01, rnorm(n, 0.025, 0.006))
  quantities <- exp(rnorm(n, 0, 0.5))
  # A plausible positive markdown ratio for an input market.
  margins <- 1 / (80 * (1 - 0.3)) / prices * exp(rnorm(n, 0, 0.05))
  list(
    prices = prices, quantities = quantities, margins = margins,
    market = paste0("M", grid$m), year = paste0("Y", grid$y),
    group = paste0("M", grid$m, "_F", grid$f)
  )
}

test_that("use_s0_market_pool and use_s0_market_year_pool are mutually exclusive", {
  d <- make_panel()
  expect_error(
    bayes_est(
      prices = d$prices, quantities = d$quantities, margins = d$margins,
      market = d$market, year = d$year, group = d$group,
      demand_model = "logit", supply_model = "bertrand",
      use_s0_market_pool = TRUE, use_s0_market_year_pool = TRUE,
      stan_options = list(chains = 1, iter = 50, warmup = 25, refresh = 0)
    ),
    "at most one"
  )
})

test_that("market-level s0 pooling yields one outside share per market (constant across years)", {
  skip_on_cran()
  d <- make_panel(seed = 21)
  fit <- suppressWarnings(bayes_est(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, year = d$year, group = d$group,
    demand_model = "logit", supply_model = "bertrand",
    use_alpha_market_pool = TRUE,
    use_s0_market_pool = TRUE, use_centered_s0_pool = TRUE,
    priors = list(prior_sigma_s0_pool = 0.3),
    stan_options = list(chains = 1, iter = 200, warmup = 100, refresh = 0, seed = 3)
  ))
  expect_s4_class(fit, "BayesLogit")

  # s0_prob is a per-market-year transformed parameter; under market-level
  # pooling it must be identical across the years of a given market (the only
  # per-market-year variation, X_s0, is absent here).
  s0_prob <- rstan::extract(fit@stanfit, "s0_prob")[[1]]  # draws x N_market_year
  mky_levels <- attr(fit@stan_data, "mky_levels")
  market_of_mky <- sub("_.*$", "", mky_levels)  # mky labels are "market_year"
  post_mean <- colMeans(s0_prob)
  for (mk in unique(market_of_mky)) {
    cols <- which(market_of_mky == mk)
    if (length(cols) > 1) {
      # All market-years of this market share one s0 draw-for-draw, so their
      # posterior means coincide to numerical precision.
      expect_lt(max(post_mean[cols]) - min(post_mean[cols]), 1e-8)
    }
  }
})

test_that("market-year s0 pooling gives a distinct outside share per market-year", {
  skip_on_cran()
  d <- make_panel(seed = 22)
  fit <- suppressWarnings(bayes_est(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, year = d$year, group = d$group,
    demand_model = "logit", supply_model = "bertrand",
    use_alpha_market_pool = TRUE,
    use_s0_market_year_pool = TRUE, use_centered_s0_pool = TRUE,
    priors = list(prior_sigma_s0_pool = 0.3),
    stan_options = list(chains = 1, iter = 200, warmup = 100, refresh = 0, seed = 3)
  ))
  expect_s4_class(fit, "BayesLogit")
  # There is a free s0 per market-year, so within a market the posterior means
  # generally differ across years (not constrained equal as in market pooling).
  s0_prob <- rstan::extract(fit@stanfit, "s0_prob")[[1]]
  mky_levels <- attr(fit@stan_data, "mky_levels")
  market_of_mky <- sub("_.*$", "", mky_levels)
  multi <- names(which(table(market_of_mky) > 1))
  skip_if(length(multi) == 0, "no multi-year market in fixture")
  post_mean <- colMeans(s0_prob)
  cols <- which(market_of_mky == multi[1])
  expect_gt(max(post_mean[cols]) - min(post_mean[cols]), 1e-8)
})
