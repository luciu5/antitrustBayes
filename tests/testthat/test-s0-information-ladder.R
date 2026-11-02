test_that("outside-share pooling resolver follows the information ladder", {
  resolve <- antitrustBayes:::.resolve_s0_pooling

  one_cell <- resolve(N_market = 1L, N_market_year = 1L)
  one_market_panel <- resolve(N_market = 1L, N_market_year = 3L)
  multi_market <- resolve(N_market = 3L, N_market_year = 3L)
  expect_equal(one_cell$resolved, "none")
  expect_equal(one_market_panel$resolved, "market")
  expect_equal(multi_market$resolved, "market")

  expect_equal(
    resolve(s0_pooling = "none", s0_pooling_supplied = TRUE,
            N_market = 3L, N_market_year = 3L)$resolved,
    "none"
  )
  expect_equal(
    resolve(s0_pooling = "market_year", s0_pooling_supplied = TRUE,
            N_market = 1L, N_market_year = 3L)$resolved,
    "market_year"
  )
  expect_equal(
    resolve(use_s0_market_pool = TRUE, legacy_supplied = TRUE,
            N_market = 1L, N_market_year = 1L)$resolved,
    "market"
  )
  expect_equal(
    resolve(use_s0_market_year_pool = TRUE, legacy_supplied = TRUE,
            N_market = 3L, N_market_year = 3L)$resolved,
    "market_year"
  )
  expect_equal(
    resolve(legacy_supplied = TRUE, N_market = 3L, N_market_year = 3L)$resolved,
    "none"
  )
  expect_error(
    resolve(use_s0_market_pool = TRUE, use_s0_market_year_pool = TRUE,
            legacy_supplied = TRUE, N_market = 2L, N_market_year = 2L),
    "at most one"
  )
  expect_error(
    resolve(s0_pooling = "auto", s0_pooling_supplied = TRUE,
            use_s0_market_pool = FALSE, legacy_supplied = TRUE,
            N_market = 1L, N_market_year = 1L),
    "cannot be supplied together"
  )
})

test_that("Logit and CES preparation retain pooling and X_s0 parity", {
  n <- 12L
  d <- list(
    prices = seq(1, 2, length.out = n),
    quantities = rep(seq_len(4), 3),
    margins = rep(0.3, n),
    market = rep(c("M1", "M2"), each = 6),
    year = rep(c("Y1", "Y2"), each = 3, times = 2),
    group = paste0("G", seq_len(n))
  )
  X_s0 <- matrix(seq_len(4), nrow = 4, ncol = 1)

  logit <- prepare_stan_data(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, year = d$year, group = d$group,
    supply_model = 1L, use_s0_market_pool = TRUE, X_s0 = X_s0,
    use_rho = FALSE, use_hmt = FALSE, priors = list()
  )
  ces <- prepare_ces_stan_data(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, year = d$year, group = d$group,
    supply_model = 1L, use_cutoff = FALSE, use_rho = FALSE,
    use_hmt = FALSE, priors = list(), use_s0_market_pool = TRUE,
    X_s0 = X_s0
  )
  for (sd in list(logit, ces)) {
    expect_equal(attr(sd, "s0_pooling"), "market")
    expect_equal(attr(sd, "s0_pooling_resolved"), "market")
    expect_equal(sd$use_s0_market_pool, 1L)
    expect_equal(sd$use_s0_market_year_pool, 0L)
    expect_equal(sd$K_s0, 1L)
    expect_equal(sd$X_s0, X_s0)
  }
})

test_that("single-market market pooling uses a direct prior source", {
  stan_dir <- .ab_test_path("inst", "stan")
  files <- file.path(stan_dir, c(
    "logit_supply.stan", "logit_stackelberg.stan",
    "ces_supply.stan", "ces_stackelberg.stan"
  ))
  for (path in files) {
    src <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_match(src, "use_s0_market_pool == 1 && is_single_market == 0",
                 fixed = TRUE)
    expect_match(src, "if (is_single_market == 1) base_logit = mu_s0_pool_vec[1];",
                 fixed = TRUE)
    expect_match(src, paste0(
      "if (is_single_market == 0) {\n",
      "      if (prior_sigma_s0_pool_family == 1) {\n",
      "        target += normal_lpdf(sigma_s0_pool_vec"
    ),
                 fixed = TRUE)
    expect_match(src, "target += normal_lpdf(mu_s0_pool_vec | s0_prior_mean, s0_prior_sd);",
                 fixed = TRUE)
  }
})

test_that("automatic one-market panels share one direct outside share", {
  skip_on_cran()
  skip_if(!identical(tolower(Sys.getenv("ANTITRUSTBAYES_RUN_STAN_TESTS")), "true"),
          "set ANTITRUSTBAYES_RUN_STAN_TESTS=true for the Stan dimension check")
  for (demand in c("logit", "ces")) {
    fit <- bayes_est(
      prices = c(1, 1.2, 1.4, 1.1, 1.3, 1.5),
      quantities = c(1, 2, 3, 2, 3, 4), margins = rep(.3, 6),
      market = rep("M", 6), year = rep(c("Y1", "Y2"), each = 3),
      group = rep(c("A", "B", "C"), 2), demand_model = demand,
      market_side = "output", use_rho = FALSE,
      share_likelihood = "gaussian", quality_parameterization = "collapsed",
      stan_options = list(chains = 1, iter = 40, warmup = 20,
        init = 0, seed = 16, refresh = 0)
    )
    expect_equal(fit@metadata$s0_pooling_requested, "auto")
    expect_equal(fit@metadata$s0_pooling_resolved, "market")
    expect_equal(attr(fit@stan_data, "s0_pooling_resolved"), "market")
    sizes <- fit@stanfit$metadata()$stan_variable_sizes
    expect_equal(as.integer(sizes[["mu_s0_pool_vec"]]), 1L)
    expect_false("sigma_s0_pool_vec" %in% names(sizes))
    expect_false("s0_market_raw" %in% names(sizes))
    expect_false("eta_s0_market" %in% names(sizes))
    outside <- antitrustBayes:::.cmdstan_extract(
      fit@stanfit, "s0_prob"
    )$s0_prob
    expect_equal(as.numeric(outside[, 1]), as.numeric(outside[, 2]),
                 tolerance = 1e-12)
  }
})
