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

# CES fixture with strictly sub-unit Lerner indices and deliberately different
# inside price indices across years.  The latter makes the quantity-share s0
# versus revenue-share R0 distinction observable in the pooling smoke test.
make_ces_panel <- function(seed = 31, n_markets = 3, n_years = 3,
                           n_firms = 4) {
  set.seed(seed)
  grid <- expand.grid(
    m = seq_len(n_markets), y = seq_len(n_years), f = seq_len(n_firms)
  )
  n <- nrow(grid)
  prices <- 1 + 0.2 * grid$m + 0.4 * grid$y + 0.05 * grid$f +
    rnorm(n, 0, 0.005)
  quantities <- exp(rnorm(n, 0, 0.5))
  margins <- 0.30 * exp(rnorm(n, 0, 0.03))
  list(
    prices = prices, quantities = quantities, margins = margins,
    market = paste0("M", grid$m), year = paste0("Y", grid$y),
    group = paste0("M", grid$m, "_F", grid$f)
  )
}

skip_stan_pooling <- function() {
  skip_on_cran()
  skip_if(
    !identical(tolower(Sys.getenv("ANTITRUSTBAYES_RUN_STAN_TESTS")), "true"),
    "set ANTITRUSTBAYES_RUN_STAN_TESTS=true to run Stan sampling tests"
  )
}

test_that("multi-market intercept hierarchies are noncentered", {
  skip_stan_pooling()

  cases <- list(
    logit = list(
      data = make_panel(seed = 19),
      args = list(demand_model = "logit", share_basis = "absolute"),
      rotated_name = "gamma_market"
    ),
    ces = list(
      data = make_ces_panel(seed = 19),
      args = list(demand_model = "ces", market_side = "output", share_basis = "absolute"),
      rotated_name = "b_market_rot"
    )
  )

  for (case in cases) {
    d <- case$data
    fit <- suppressWarnings(do.call(
      bayes_est,
      c(
        list(
          prices = d$prices, quantities = d$quantities, margins = d$margins,
          market = d$market, year = d$year, group = d$group,
          supply_model = "bertrand", use_alpha_market_pool = TRUE,
          use_s0_market_pool = TRUE,
          priors = list(prior_sigma_s0_pool = 0.3),
          stan_options = list(
            chains = 1, iter = 160, warmup = 80, refresh = 0, seed = 3,
            sig_figs = 18L
          )
        ),
        case$args
      )
    ))

    expect_equal(
      as.integer(fit@stanfit$metadata()$stan_variable_sizes[["b_market_raw"]]),
      3L
    )
    draws <- antitrustBayes:::.cmdstan_extract(
      fit@stanfit,
      variables = c(
        "mu_b_market_vec", "sigma_b_market_vec", "b_market_raw",
        "b_market", case$rotated_name, "mean_logit_s0_market"
      )
    )
    mu <- as.numeric(draws$mu_b_market_vec)
    sigma <- as.numeric(draws$sigma_b_market_vec)
    for (m in seq_len(3L)) {
      expect_equal(
        draws$b_market[, m],
        mu + sigma * draws$b_market_raw[, m],
        tolerance = 1e-10
      )
      expect_equal(
        draws[[case$rotated_name]][, m],
        draws$b_market[, m] + draws$mean_logit_s0_market[, m],
        tolerance = 1e-10
      )
    }
  }
})

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
  d_ces <- make_ces_panel()
  expect_error(
    bayes_est(
      prices = d_ces$prices, quantities = d_ces$quantities,
      margins = d_ces$margins, market = d_ces$market,
      year = d_ces$year, group = d_ces$group,
      demand_model = "ces", market_side = "output", supply_model = "bertrand",
      use_s0_market_pool = TRUE, use_s0_market_year_pool = TRUE,
      stan_options = list(chains = 1, iter = 50, warmup = 25, refresh = 0)
    ),
    "at most one"
  )
})

test_that("CES market-level s0 pooling reaches the Stan data contract", {
  d <- make_ces_panel()
  for (centered in c(FALSE, TRUE)) {
    sd <- prepare_ces_stan_data(
      prices = d$prices, quantities = d$quantities, margins = d$margins,
      market = d$market, year = d$year, group = d$group,
      supply_model = 1L, use_cutoff = FALSE, use_rho = TRUE,
      use_hmt = FALSE, priors = list(),
      use_s0_market_year_pool = FALSE,
      use_s0_market_pool = TRUE,
      use_centered_s0_pool = centered
    )
    expect_equal(sd$use_s0_market_pool, 1L)
    expect_equal(sd$use_s0_market_year_pool, 0L)
    expect_equal(sd$N_market, 3L)
    expect_equal(sd$N_market_year, 9L)
    expect_equal(sd$K_s0, 0L)
    expect_equal(length(sd$inside_price_index), sd$N_market_year)
  }
})

test_that("market-level s0 pooling yields one outside share per market (constant across years)", {
  skip_stan_pooling()
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
  s0_prob <- antitrustBayes:::.cmdstan_extract(
    fit@stanfit, "s0_prob"
  )$s0_prob  # draws x N_market_year
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

test_that("CES market-level s0 pooling keeps quantity s0 fixed but transforms R0 by prices", {
  skip_stan_pooling()
  d <- make_ces_panel()
  for (centered in c(FALSE, TRUE)) {
    fit <- suppressWarnings(bayes_est(
      prices = d$prices, quantities = d$quantities, margins = d$margins,
      market = d$market, year = d$year, group = d$group,
      demand_model = "ces", market_side = "output", supply_model = "bertrand",
      use_alpha_market_pool = TRUE,
      use_s0_market_pool = TRUE, use_centered_s0_pool = centered,
      priors = list(prior_sigma_s0_pool = 0.3),
      stan_options = list(chains = 1, iter = 160, warmup = 80, refresh = 0, seed = 3,
                          sig_figs = 18L)
    ))
    expect_s4_class(fit, "BayesLogit")

    draws <- antitrustBayes:::.cmdstan_extract(
      fit@stanfit, c("s0_prob", "R0_prob")
    )
    s0_prob <- draws$s0_prob
    R0_prob <- draws$R0_prob
    mky_levels <- attr(fit@stan_data, "mky_levels")
    market_of_mky <- sub("_.*$", "", mky_levels)

    # Independently recover outside revenue share by expenditure accounting:
    # R0 = s0*p0 / (s0*p0 + (1-s0)*pbar_inside).
    for (j in seq_len(fit@stan_data$N_market_year)) {
      outside_revenue <- s0_prob[, j] * fit@stan_data$outside_price[j]
      inside_revenue <- (1 - s0_prob[, j]) *
        fit@stan_data$inside_price_index[j]
      expected_R0 <- outside_revenue / (outside_revenue + inside_revenue)
      expect_equal(R0_prob[, j], expected_R0, tolerance = 1e-10)
    }

    for (mk in unique(market_of_mky)) {
      cols <- which(market_of_mky == mk)
      if (length(cols) < 2L) next

      # A market-pool effect is shared draw-for-draw, not merely in posterior
      # expectation. X_s0 is NULL in this fit, so no year-specific term exists.
      for (j in cols[-1L]) {
        expect_equal(s0_prob[, j], s0_prob[, cols[1L]], tolerance = 1e-12)
      }

      # R0 can still vary because the CES conversion uses period-specific
      # inside and outside prices.
      r0_diffs <- vapply(cols[-1L], function(j) {
        max(abs(R0_prob[, j] - R0_prob[, cols[1L]]))
      }, numeric(1))
      expect_gt(max(r0_diffs), 1e-8)
    }
  }
})

test_that("market-year s0 pooling gives a distinct outside share per market-year", {
  skip_stan_pooling()
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
  s0_prob <- antitrustBayes:::.cmdstan_extract(
    fit@stanfit, "s0_prob"
  )$s0_prob
  mky_levels <- attr(fit@stan_data, "mky_levels")
  market_of_mky <- sub("_.*$", "", mky_levels)
  multi <- names(which(table(market_of_mky) > 1))
  skip_if(length(multi) == 0, "no multi-year market in fixture")
  post_mean <- colMeans(s0_prob)
  cols <- which(market_of_mky == multi[1])
  expect_gt(max(post_mean[cols]) - min(post_mean[cols]), 1e-8)
})
