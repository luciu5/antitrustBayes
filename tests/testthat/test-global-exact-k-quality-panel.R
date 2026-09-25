test_that("global collapsed quality reduces exact-k states to market scope", {
  n <- 4L
  data <- data.frame(
    prices = rep(1, n), quantities = rep(1, n), margins = rep(.3, n),
    market = rep("m1", n), year = rep(c("y1", "y2"), each = 2L),
    group = letters[seq_len(n)], outside_price = rep(1, n),
    stringsAsFactors = FALSE
  )
  stan_data <- list(
    N = n, N_market = 1L, N_market_year = 2L,
    market = rep(1L, n), market_year_idx = rep(1:2, each = 2L),
    group = seq_len(n), log_shareIn = rep(log(.25), n),
    likelihood_mode = 1L, use_rho = 1L, margin_idx = 1L,
    log_marginNormInv = c(log(2), rep(0, n - 1L)),
    N_partition_unit = 1L, partition_unit_market = 1L,
    partition_component_k = 1L, partition_component_unit = 1L,
    partition_unit_count = 1L
  )
  attr(stan_data, "model_spec_version") <-
    antitrustBayes:::.antitrust_bayes_model_spec_version
  attr(stan_data, "market_levels") <- "m1"
  attr(stan_data, "group_levels") <- letters[seq_len(n)]
  attr(stan_data, "quality_scope") <- "global"
  attr(stan_data, "quality_parameterization") <- "collapsed"
  attr(stan_data, "share_likelihood") <- "gaussian"

  fit <- methods::new(
    "BayesLogit", stanfit = NULL, data = data, stan_data = stan_data,
    supply_model = 1L, demand_model = "logit", market_side = "output",
    price_sd = 1, call = quote(fake()),
    metadata = list(
      model_spec_version = antitrustBayes:::.antitrust_bayes_model_spec_version,
      quality_scope = "global", quality_parameterization = "collapsed",
      share_likelihood = "gaussian", structure = "stackelberg",
      partition_prior = "ordered_k", k_scope = "market",
      posterior_draws = 1L
    )
  )

  incidence <- antitrustBayes:::.lq_quality_incidence(fit)
  year_states <- structure(
    list(c(m1_y1 = 1L, m1_y2 = 1L)), names = "1"
  )
  market_states <- antitrustBayes:::.lq_market_quality_k_cache(
    fit, incidence, year_states
  )

  expect_length(market_states[[1L]], 1L)
  expect_equal(unname(market_states[[1L]]), 1)
  recovered <- antitrustBayes:::.posterior_quality_exact_k(
    fit,
    parsed = list(n_draw = 1L, n_market = 1L, quality_scope = "global"),
    draws = list(posterior_component_prob = matrix(1, nrow = 1L)),
    draw_indices = 1L,
    k_cache = market_states
  )
  expect_identical(unname(recovered), matrix(1L, nrow = 1L, ncol = 1L))
})
