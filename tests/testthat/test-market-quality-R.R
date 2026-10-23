test_that("market quality incidence has a separate contiguous CSR map", {
  market <- c(1L, 2L, 1L, 2L, 1L)
  group <- c(2L, 1L, 2L, 2L, 1L)
  original_group <- group
  inc <- antitrustBayes:::.build_quality_incidence(
    market, group, c("east", "west"), c("alpha", "beta"),
    quality_scope = "market", margin_idx = c(2L, 3L, 5L)
  )
  expect_identical(group, original_group)
  expect_identical(inc$obs_to_quality, c(2L, 3L, 2L, 4L, 1L))
  expect_identical(inc$quality_block_start, c(1L, 3L))
  expect_identical(inc$quality_block_n, c(2L, 2L))
  expect_identical(inc$quality_coord_start, c(1L, 2L))
  expect_identical(inc$market_margin_order, c(2L, 3L, 1L))
  expect_identical(inc$market_margin_start, c(1L, 3L))
  expect_identical(inc$market_margin_n, c(2L, 1L))
  expect_identical(as.character(inc$quality_lookup$quality_unit),
                   c("m1_g1", "m1_g2", "m2_g1", "m2_g2"))
})

test_that("fully Gaussian likelihood has no active df coordinate", {
  modes <- antitrustBayes:::.resolve_likelihood_modes(
    share_likelihood = "gaussian", margin_likelihood = "gaussian",
    quality_parameterization = "collapsed", estimate_df = FALSE
  )
  expect_false(modes$df_active)
  expect_identical(modes$quality_scope, "market")
  expect_error(
    antitrustBayes:::.resolve_likelihood_modes(
      share_likelihood = "gaussian", margin_likelihood = "gaussian",
      quality_parameterization = "collapsed", estimate_df = TRUE
    ),
    "only available"
  )
})

test_that("missing historical quality scope is explicitly global", {
  stan_data <- list(group = c(1L, 2L), N_group = 2L)
  fit <- methods::new(
    "BayesLogit", stanfit = NULL,
    data = data.frame(group = c("a", "b")), stan_data = stan_data,
    price_sd = 1, supply_model = 1L, demand_model = "logit",
    market_side = "input", call = quote(fake()), metadata = list()
  )
  expect_identical(antitrustBayes:::.quality_scope_from_fit(fit), "global")
  attr(fit@stan_data, "model_spec_version") <- 13L
  expect_error(antitrustBayes:::.quality_scope_from_fit(fit), "missing quality_scope")
})

test_that("scope metadata must agree and economic group stays independent", {
  stan_data <- list(group = c(1L, 2L), quality_scope = 2L)
  attr(stan_data, "quality_scope") <- "global"
  fit <- methods::new(
    "BayesLogit", stanfit = NULL,
    data = data.frame(group = c("a", "b")), stan_data = stan_data,
    price_sd = 1, supply_model = 1L, demand_model = "logit",
    market_side = "input", call = quote(fake()), metadata = list()
  )
  expect_error(antitrustBayes:::.quality_scope_from_fit(fit), "disagree")
})

test_that("market quality persists by group across years and centers per market", {
  inc <- antitrustBayes:::.build_quality_incidence(
    market_int = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
    group_int = c(1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L),
    market_levels = c("east", "west"), group_levels = c("A", "B"),
    quality_scope = "market"
  )
  expect_identical(inc$obs_to_quality, c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L))
  expect_identical(inc$quality_block_start, c(1L, 3L))
  expect_identical(inc$quality_block_n, c(2L, 2L))
  expect_identical(inc$quality_coord_start, c(1L, 2L))
  expect_identical(inc$quality_lookup$market_index, c(1L, 1L, 2L, 2L))
  expect_identical(inc$quality_lookup$group_index, c(1L, 2L, 1L, 2L))
})

test_that("posterior quality uses scoped blocks and accepts missing margins", {
  inc <- antitrustBayes:::.build_quality_incidence(
    c(1L, 1L, 2L, 2L), c(1L, 2L, 1L, 2L),
    c("east", "west"), c("A", "B"), quality_scope = "market"
  )
  stan_data <- list(
    N = 4L, N_market = 2L, market = c(1L, 1L, 2L, 2L),
    N_group = 2L, group = c(1L, 2L, 1L, 2L),
    log_shareIn = c(-.1, .2, .3, .4),
    N_quality = 4L, quality_scope = 2L,
    quality_parameterization = 2L, share_likelihood = 2L,
    likelihood_mode = 2L, use_rho = 0L, margin_idx = integer()
  )
  attr(stan_data, "quality_incidence") <- inc
  attr(stan_data, "quality_scope") <- "market"
  attr(stan_data, "group_levels") <- c("A", "B")
  fit <- methods::new(
    "BayesLogit", stanfit = NULL,
    data = data.frame(market = c("east", "east", "west", "west"),
                      group = c("A", "B", "A", "B")),
    stan_data = stan_data, price_sd = 1, supply_model = 1L,
    demand_model = "logit", market_side = "input", call = quote(fake()),
    metadata = list(model_spec_version = 13L, quality_scope = "market",
                    quality_parameterization = "collapsed",
                    share_likelihood = "gaussian", posterior_draws = 2L,
                    group_levels = c("A", "B"))
  )
  got <- posterior_quality(
    fit,
    draws = list(
      share_mean = rbind(c(0, 0, 0, 0), c(.1, .2, .3, .4)),
      sigma_b_strat = c(1, 1), sigma_share = c(1, 1)
    ),
    seed = 11
  )
  expect_identical(got$quality_scope, "market")
  expect_identical(dim(got$draws), c(2L, 4L))
  expect_equal(unname(rowSums(got$draws[, 1:2, drop = FALSE])), c(0, 0), tolerance = 1e-8)
  expect_equal(unname(rowSums(got$draws[, 3:4, drop = FALSE])), c(0, 0), tolerance = 1e-8)
  expect_identical(got$metadata$conditional_moments,
                   "conditional_on_observed_quality_state")
  expect_true(all(c("market", "group", "quality_unit", "quality") %in%
                  names(got$draws_long)))
})

test_that("exact-k quality recovery maps one supported state per market", {
  inc <- antitrustBayes:::.build_quality_incidence(
    c(1L, 1L, 2L, 2L), c(1L, 2L, 1L, 2L),
    c("east", "west"), c("A", "B"), quality_scope = "market"
  )
  stan_data <- list(
    N = 4L, N_market = 2L, market = c(1L, 1L, 2L, 2L),
    N_group = 2L, group = c(1L, 2L, 1L, 2L),
    N_quality = 4L, quality_scope = 2L, quality_parameterization = 2L,
    share_likelihood = 2L, likelihood_mode = 1L, use_rho = 1L,
    margin_idx = c(1L, 3L), N_partition_unit = 2L,
    partition_unit_market = c(2L, 1L), partition_unit_count = c(2L, 2L),
    partition_component_unit = c(1L, 1L, 2L, 2L),
    partition_component_k = c(1L, 2L, 1L, 2L)
  )
  attr(stan_data, "quality_incidence") <- inc
  attr(stan_data, "quality_scope") <- "market"
  fit <- methods::new(
    "BayesLogit", stanfit = NULL,
    data = data.frame(market = c("east", "east", "west", "west"),
                      group = c("A", "B", "A", "B")),
    stan_data = stan_data, price_sd = 1, supply_model = 1L,
    demand_model = "logit", market_side = "input", call = quote(fake()),
    metadata = list(model_spec_version = 13L, quality_scope = "market",
                    structure = "stackelberg", partition_prior = "ordered_k")
  )
  draws <- list(posterior_component_prob = rbind(
    c(.25, .75, .6, .4), c(.5, .5, .2, .8)
  ))
  got <- antitrustBayes:::.posterior_quality_exact_k(
    fit, list(n_draw = 2L), draws, draw_indices = 1L,
    k_cache = list(`1` = c(2L, 1L))
  )
  expect_identical(as.integer(got[1, ]), c(1L, 2L))
  singleton <- stan_data
  singleton$N_partition_unit <- 2L
  singleton$partition_unit_market <- c(1L, 2L)
  singleton$partition_unit_count <- c(1L, 1L)
  singleton$partition_component_unit <- c(1L, 2L)
  singleton$partition_component_k <- c(2L, 1L)
  fit@stan_data <- singleton
  singleton_draw <- antitrustBayes:::.posterior_quality_exact_k(
    fit, list(n_draw = 1L),
    list(posterior_component_prob = matrix(c(1, 1), nrow = 1L)),
    draw_indices = 1L
  )
  expect_identical(as.integer(singleton_draw[1, ]), c(2L, 1L))
  bad <- stan_data
  bad$N_partition_unit <- 3L
  bad$partition_unit_market <- c(1L, 1L, 2L)
  bad$partition_unit_count <- c(2L, 2L, 2L)
  bad$partition_component_unit <- rep(1:3, each = 2L)
  bad$partition_component_k <- rep(c(1L, 2L), 3L)
  fit@stan_data <- bad
  expect_error(
    antitrustBayes:::.posterior_quality_exact_k(
      fit, list(n_draw = 1L), draws, draw_indices = 1L
    ),
    "one unit per market"
  )
})

test_that("current specification missing scope errors even through non-strict callers", {
  stan_data <- list(group = c(1L, 2L), N_group = 2L)
  attr(stan_data, "model_spec_version") <- 13L
  fit <- methods::new(
    "BayesLogit", stanfit = NULL,
    data = data.frame(group = c("a", "b")), stan_data = stan_data,
    price_sd = 1, supply_model = 1L, demand_model = "logit",
    market_side = "input", call = quote(fake()), metadata = list()
  )
  expect_error(
    antitrustBayes:::.quality_scope_from_fit(fit, strict = FALSE),
    "missing quality_scope"
  )
})

market_quality_fingerprint_fit <- function(quality_scope,
                                           market_levels = "m1") {
  market_int <- if (length(market_levels) == 1L) rep(1L, 4L) else
    rep(seq_along(market_levels), each = 2L)
  group_int <- rep(c(1L, 2L), length.out = length(market_int))
  incidence <- antitrustBayes:::.build_quality_incidence(
    market_int, group_int, market_levels, c("A", "B"), quality_scope,
    margin_idx = integer()
  )
  stan_data <- c(list(
    N = length(market_int), N_market = length(market_levels),
    N_market_year = length(market_levels), market = market_int,
    market_year_idx = market_int, N_group = 2L, group = group_int,
    shareIn = rep(.4, length(market_int)),
    marginNormInv = rep(.5, length(market_int)), margin_idx = integer(),
    use_cutoff = 0L, use_rho = 0L, use_alpha_market_pool = 0L,
    use_s0_market_pool = 0L, use_s0_market_year_pool = 0L,
    use_market_resid_scale = 0L, use_hmt = 0L, estimate_df = 0L,
    likelihood_mode = 2L, share_likelihood = 2L, margin_likelihood = 2L,
    quality_parameterization = 2L, quality_scope =
      if (quality_scope == "global") 1L else 2L,
    prior_alpha_sd = 1
  ), incidence[setdiff(names(incidence),
                       c("quality_scope", "quality_scope_code", "quality_lookup"))])
  attr(stan_data, "quality_incidence") <- incidence
  attr(stan_data, "quality_scope") <- quality_scope
  attr(stan_data, "market_levels") <- market_levels
  attr(stan_data, "model_spec_version") <- 13L
  data <- data.frame(
    prices = rep(1, length(market_int)), quantities = rep(1, length(market_int)),
    margins = rep(.5, length(market_int)),
    market = market_levels[market_int], year = "1",
    group = c("A", "B")[group_int], margin_observed = TRUE,
    outside_price = 1, stringsAsFactors = FALSE
  )
  methods::new(
    "BayesLogit", stanfit = NULL, data = data, stan_data = stan_data,
    price_sd = 1, supply_model = 1L, demand_model = "logit",
    market_side = "input", call = quote(fake()), metadata = list(
      model_spec_version = 13L, quality_scope = quality_scope,
      quality_parameterization = "collapsed", share_likelihood = "gaussian",
      margin_likelihood = "gaussian", structure = "simultaneous",
      partition_prior = "cutoff", posterior_draws = 1L,
      prior_normalization = list(model_spec_version = 13L),
      resolved_priors = list(prior_alpha_sd = 1)
    )
  )
}

test_that("one-market global and market quality share a model-universe identity", {
  global <- market_quality_fingerprint_fit("global")
  market <- market_quality_fingerprint_fit("market")
  expect_identical(antitrustBayes:::.mu_fingerprint(global),
                   antitrustBayes:::.mu_fingerprint(market))
  multi <- market_quality_fingerprint_fit("market", c("m1", "m2"))
  expect_false(identical(antitrustBayes:::.mu_fingerprint(global),
                         antitrustBayes:::.mu_fingerprint(multi)))
})

test_that("market exact-k policy cache is shared across panel years", {
  fit <- market_quality_fingerprint_fit("market", c("m1", "m2"))
  market_int <- rep(c(1L, 1L, 2L, 2L), 2L)
  group_int <- rep(c(1L, 2L), 4L)
  incidence <- antitrustBayes:::.build_quality_incidence(
    market_int, group_int, c("m1", "m2"), c("A", "B"), "market"
  )
  fit@data <- data.frame(
    prices = rep(1, 8), quantities = rep(1, 8), margins = rep(.5, 8),
    market = c(rep("m1", 2), rep("m2", 2), rep("m1", 2), rep("m2", 2)),
    year = rep(c("y1", "y2"), each = 4),
    group = c("A", "B", "A", "B", "A", "B", "A", "B"),
    outside_price = 1, stringsAsFactors = FALSE
  )
  fit@stan_data$market <- market_int
  fit@stan_data$group <- group_int
  fit@stan_data$N <- 8L
  fit@stan_data$N_market_year <- 4L
  fit@stan_data$market_year_idx <- rep(1:4, each = 2L)
  attr(fit@stan_data, "quality_incidence") <- incidence
  raw <- list(draw1 = c(m1_y1 = 1L, m1_y2 = 1L, m2_y1 = 2L, m2_y2 = 2L))
  got <- antitrustBayes:::.lq_market_quality_k_cache(fit, incidence, raw)
  expect_identical(as.integer(got[[1L]]), c(1L, 2L))
  raw_bad <- raw
  raw_bad[[1L]][["m1_y2"]] <- 2L
  expect_error(
    antitrustBayes:::.lq_market_quality_k_cache(fit, incidence, raw_bad),
    "one exact-k state"
  )
})
