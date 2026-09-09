test_that("validate_inputs catches mismatched lengths", {
  expect_error(
    validate_inputs(1:3, 1:2, c(0.2, 0.3, 0.4), 1:3, 1:3, 1:3),
    "'quantities' must have the same length"
  )
  expect_error(
    validate_inputs(1:3, 1:3, c(0.2, 0.3), 1:3, 1:3, 1:3),
    "'margins' must have the same length"
  )
})

test_that("validate_inputs rejects empty data and missing identifiers", {
  expect_error(
    validate_inputs(numeric(), numeric(), numeric(), character(),
                    character(), character()),
    "At least one observation"
  )
  expect_error(
    validate_inputs(1, 1, 0.2, NA_character_, "1", "g"),
    "identifiers may not be missing"
  )
})

test_that("validate_inputs catches bad margins", {
  expect_error(
    validate_inputs(
      c(1, 2, 3), c(10, 20, 30), c(-0.1, 0.2, 0.3),
      c(1, 1, 1), c(1, 1, 1), c(1, 2, 3)
    ),
    "'margins' must be positive and finite"
  )
  expect_error(
    validate_inputs(
      c(1, 2, 3), c(10, 20, 30), c(0.2, Inf, 0.3),
      c(1, 1, 1), c(1, 1, 1), c(1, 2, 3)
    ),
    "'margins' must be positive and finite"
  )
})

test_that("logit preparation allows input markdown ratios above one", {
  prices <- c(0.6, 0.8, 1.0, 1.2)
  margins <- c(0.5, 1.2, 2.0, 3.5)
  sd <- prepare_stan_data(
    prices       = prices,
    quantities   = c(10, 20, 30, 40),
    margins      = margins,
    market       = rep(1, 4),
    year         = rep(1, 4),
    group        = 1:4,
    supply_model = 1L,
    use_cutoff   = FALSE,
    use_rho      = TRUE,
    use_hmt      = FALSE,
    priors       = list()
  )

  expect_equal(sd$marginNormInv, 1 / margins)
  expect_equal(sd$log_marginNormInv, log(sd$marginNormInv))
})

test_that("validate_inputs catches bad quantities", {
  expect_error(
    validate_inputs(
      c(1, 2, 3), c(0, 20, 30), c(0.2, 0.3, 0.4),
      c(1, 1, 1), c(1, 1, 1), c(1, 2, 3)
    ),
    "'quantities' must be positive"
  )
})

test_that("prepare_stan_data returns correct structure", {
  N <- 10
  sd <- prepare_stan_data(
    prices       = runif(N, 1, 5),
    quantities   = runif(N, 10, 100),
    margins      = runif(N, 0.1, 0.5),
    market       = rep(1:2, each = 5),
    year         = rep(1, N),
    group        = sample(1:5, N, replace = TRUE),
    supply_model = 1L,
    use_cutoff   = FALSE,
    use_rho      = TRUE,
    use_hmt      = FALSE,
    priors       = list()
  )

  expect_true(is.list(sd))
  expect_equal(sd$N, N)
  expect_equal(sd$N_market, 2L)
  expect_equal(sd$supply_model, 1L)
  expect_equal(sd$use_rho, 1L)
  expect_equal(sd$use_cutoff, 0L)
  expect_equal(length(sd$shareIn), N)
  expect_equal(length(sd$marginNormInv), N)
  expect_equal(length(sd$log_shareIn), N)
  expect_equal(length(sd$clr_shareIn), N)
  expect_equal(length(sd$log_marginNormInv), N)
  expect_equal(sd$log_shareIn, log(sd$shareIn))
  expect_equal(sd$log_marginNormInv, log(sd$marginNormInv))
  expect_equal(length(sd$price_std), N)
  expect_true(sd$price_sd > 0)
  expect_false("implied_s0" %in% names(sd))
  expect_equal(attr(sd, "model_spec_version"), 6L)
  expect_equal(attr(sd, "s0_prior_reference"), 0.2, tolerance = 1e-8)
})

test_that("default outside-share prior has mean 20% and 5th percentile 10%", {
  sd <- prepare_stan_data(
    prices       = c(1, 2, 3, 4),
    quantities   = c(40, 30, 20, 10),
    margins      = rep(0.3, 4),
    market       = rep(1, 4),
    year         = rep(1, 4),
    group        = 1:4,
    supply_model = 1L,
    use_cutoff   = FALSE,
    use_rho      = FALSE,
    use_hmt      = FALSE,
    priors       = list()
  )

  prior_mean <- logit_normal_mean(sd$s0_prior_mean, sd$s0_prior_sd)
  prior_q05 <- plogis(sd$s0_prior_mean + qnorm(0.05) * sd$s0_prior_sd)
  expect_equal(prior_mean, 0.20, tolerance = 1e-8)
  expect_equal(prior_q05, 0.10, tolerance = 1e-8)
})

test_that("single-market CES vector data retain Stan dimensions", {
  sd <- prepare_ces_stan_data(
    prices = c(1, 1.2, 1.4), quantities = c(10, 20, 30),
    margins = rep(0.3, 3), market = rep("m", 3), year = rep(1, 3),
    group = 1:3, supply_model = 1L, use_cutoff = FALSE, use_rho = FALSE,
    use_hmt = FALSE, priors = list()
  )
  expect_equal(dim(sd$outside_price), 1L)
  expect_equal(dim(sd$inside_price_index), 1L)
  expect_equal(dim(sd$mky_n), 1L)
})

test_that("custom outside-share prior does not calibrate alpha from outcomes", {
  prices <- c(1, 2, 3, 4)
  quantities <- c(40, 30, 20, 10)
  margins <- rep(0.3, 4)
  custom_location <- qlogis(0.30)
  custom_scale <- 0.20
  sd <- prepare_stan_data(
    prices       = prices,
    quantities   = quantities,
    margins      = margins,
    market       = rep(1, 4),
    year         = rep(1, 4),
    group        = 1:4,
    supply_model = 1L,
    use_cutoff   = FALSE,
    use_rho      = FALSE,
    use_hmt      = FALSE,
    priors       = list(
      s0_prior_mean = custom_location,
      s0_prior_sd = custom_scale
    )
  )

  expected_reference <- logit_normal_mean(custom_location, custom_scale)
  expect_equal(attr(sd, "s0_prior_reference"), expected_reference)
  expect_equal(sd$prior_alpha_mean, 0)
  expect_equal(sd$prior_alpha_sd, 1)
})

test_that("prepare_stan_data computes quantity shares correctly", {
  N <- 6
  prices <- rep(2, N)
  quantities <- c(10, 20, 30, 40, 50, 60)
  sd <- prepare_stan_data(
    prices       = prices,
    quantities   = quantities,
    margins      = rep(0.3, N),
    market       = rep(1, N),
    year         = rep(1, N),
    group        = 1:6,
    supply_model = 1L,
    use_cutoff   = FALSE,
    use_rho      = TRUE,
    use_hmt      = FALSE,
    priors       = list()
  )

  expected_shares <- quantities / sum(quantities)
  expect_equal(sd$shareIn, expected_shares)
  expect_equal(sum(sd$clr_shareIn), 0, tolerance = 1e-12)
})

test_that("compositional share outcomes are scale invariant", {
  args <- list(
    prices = c(1, 2, 3, 4), quantities = c(10, 20, 30, 40),
    margins = rep(0.3, 4), market = rep("m", 4), year = rep(1, 4),
    group = letters[1:4], supply_model = 1L, use_cutoff = FALSE,
    use_rho = FALSE, use_hmt = FALSE, priors = list()
  )
  base <- do.call(prepare_stan_data, args)
  args$quantities <- args$quantities * 100
  scaled <- do.call(prepare_stan_data, args)
  expect_equal(base$clr_shareIn, scaled$clr_shareIn)
})

test_that("predictive clusters must respect equilibrium cells", {
  expect_error(
    prepare_stan_data(
      prices = c(1, 2), quantities = c(1, 1), margins = c(0.2, 0.2),
      market = c("m", "m"), year = c(1, 1), group = c("a", "b"),
      supply_model = 1L, use_cutoff = FALSE, use_rho = FALSE,
      use_hmt = FALSE, priors = list(), loo_cluster = c("x", "y")
    ),
    "constant within each market-year"
  )
})

test_that("strategic FOCs receive firm shares and market direction", {
  sd <- prepare_stan_data(
    prices       = c(1, 2, 3),
    quantities   = c(2, 3, 5),
    margins      = rep(0.2, 3),
    market       = rep("m", 3),
    year         = rep(1, 3),
    group        = c("a", "a", "b"),
    supply_model = 3L,
    use_cutoff   = FALSE,
    use_rho      = FALSE,
    use_hmt      = FALSE,
    market_side  = "output",
    priors       = list()
  )

  expect_equal(sd$market_sign, -1L)
  expect_equal(sd$firm_shareIn, c(0.5, 0.5, 0.5))
})

test_that("Logit alpha prior follows the standardized price scale", {
  base <- compute_alpha_prior(
    shares = c(0.2, 0.3, 0.5),
    margins_abs = c(1, 1, 1),
    s0_reference = 0.2,
    mky_int = rep(1L, 3),
    supply_model = 2L,
    price_sd = 1
  )
  scaled <- compute_alpha_prior(
    shares = c(0.2, 0.3, 0.5),
    margins_abs = c(1, 1, 1),
    s0_reference = 0.2,
    mky_int = rep(1L, 3),
    supply_model = 2L,
    price_sd = 2
  )
  expect_equal(scaled$mean, base$mean + log(2), tolerance = 1e-10)
})

test_that("prepare_stan_data handles constant prices and custom grainsize", {
  N <- 6
  sd <- prepare_stan_data(
    prices       = rep(2, N),
    quantities   = runif(N, 10, 100),
    margins      = rep(0.3, N),
    market       = rep(1, N),
    year         = rep(1, N),
    group        = 1:6,
    supply_model = 1L,
    use_cutoff   = FALSE,
    use_rho      = TRUE,
    use_hmt      = FALSE,
    priors       = list(),
    grainsize    = 3L
  )

  expect_equal(sd$price_sd, 1)
  expect_equal(sd$grainsize, 3L)
  expect_true(all(is.finite(sd$price_std)))
})

test_that("prepare_stan_data handles single market", {
  N <- 5
  sd <- prepare_stan_data(
    prices       = runif(N, 1, 5),
    quantities   = runif(N, 10, 100),
    margins      = runif(N, 0.1, 0.5),
    market       = rep(1, N),
    year         = rep(1, N),
    group        = 1:5,
    supply_model = 1L,
    use_cutoff   = FALSE,
    use_rho      = TRUE,
    use_hmt      = FALSE,
    priors       = list()
  )

  expect_equal(sd$is_single_market, 1L)
  expect_equal(sd$N_market, 1L)
})

test_that("prepare_stan_data respects custom priors", {
  N <- 5
  sd <- prepare_stan_data(
    prices       = runif(N, 1, 5),
    quantities   = runif(N, 10, 100),
    margins      = runif(N, 0.1, 0.5),
    market       = rep(1, N),
    year         = rep(1, N),
    group        = 1:5,
    supply_model = 1L,
    use_cutoff   = FALSE,
    use_rho      = TRUE,
    use_hmt      = FALSE,
    priors       = list(prior_alpha_mean = 2.0, prior_alpha_sd = 0.5)
  )

  expect_equal(sd$prior_alpha_mean, 2.0)
  expect_equal(sd$prior_alpha_sd, 0.5)
})



test_that("prepare_stan_data accepts s0 covariates", {
  N <- 6
  X_s0 <- cbind(intercept = 1, trend = c(0, 1))
  sd <- prepare_stan_data(
    prices       = runif(N, 1, 5),
    quantities   = runif(N, 10, 100),
    margins      = runif(N, 0.1, 0.5),
    market       = rep(1, N),
    year         = rep(1:2, each = 3),
    group        = 1:6,
    supply_model = 1L,
    use_cutoff   = FALSE,
    use_rho      = TRUE,
    use_hmt      = FALSE,
    priors       = list(),
    X_s0         = X_s0
  )

  expect_equal(sd$K_s0, 2L)
  expect_equal(sd$X_s0, X_s0)
})

test_that("prepare_stan_data validates s0 covariate rows", {
  N <- 6
  expect_error(
    prepare_stan_data(
      prices       = runif(N, 1, 5),
      quantities   = runif(N, 10, 100),
      margins      = runif(N, 0.1, 0.5),
      market       = rep(1, N),
      year         = rep(1:2, each = 3),
      group        = 1:6,
      supply_model = 1L,
      use_cutoff   = FALSE,
      use_rho      = TRUE,
      use_hmt      = FALSE,
      priors       = list(),
      X_s0         = matrix(1, nrow = 1, ncol = 2)
    ),
    "one row per market-year"
  )
})

test_that("preparation validates outside prices, priors, and covariates", {
  args <- list(
    prices = c(1, 2), quantities = c(1, 1), margins = c(0.2, 0.2),
    market = c("m", "m"), year = c(1, 1), group = c("a", "b"),
    supply_model = 1L, use_cutoff = FALSE, use_rho = FALSE,
    use_hmt = FALSE
  )
  expect_error(
    do.call(prepare_stan_data, c(args, list(
      priors = list(), outside_price = 0
    ))),
    "outside_price"
  )
  expect_error(
    do.call(prepare_stan_data, c(args, list(
      priors = list(prior_sigma_share = 0)
    ))),
    "strictly positive"
  )
  expect_error(
    do.call(prepare_stan_data, c(args, list(
      priors = list(), X_s0 = matrix(c(1, Inf), nrow = 1)
    ))),
    "finite numeric matrix"
  )
})

test_that("all-masked margins remain finite when HMT is disabled", {
  sd <- prepare_stan_data(
    prices = c(1, 2), quantities = c(1, 2), margins = c(NA, NA),
    market = c("m", "m"), year = c(1, 1), group = c("a", "b"),
    supply_model = 1L, use_cutoff = FALSE, use_rho = FALSE,
    use_hmt = FALSE, priors = list(), margin_observed = c(FALSE, FALSE)
  )
  expect_true(all(is.finite(sd$marginNormInv)))
  expect_equal(sd$avg_margin_hmt, 0)
  expect_error(
    bayes_est(
      prices = c(1, 2), quantities = c(1, 2), margins = c(NA_real_, NA_real_),
      market = c("m", "m"), group = c("a", "b"),
      use_hmt = TRUE, margin_observed = c(FALSE, FALSE)
    ),
    "requires at least one observed margin"
  )
})

test_that("prepare_stan_data warns on unknown prior", {
  N <- 5
  expect_warning(
    prepare_stan_data(
      prices       = runif(N, 1, 5),
      quantities   = runif(N, 10, 100),
      margins      = runif(N, 0.1, 0.5),
      market       = rep(1, N),
      year         = rep(1, N),
      group        = 1:5,
      supply_model = 1L,
      use_cutoff   = FALSE,
      use_rho      = TRUE,
      use_hmt      = FALSE,
      priors       = list(nonsense_param = 99)
    ),
    "Unknown prior parameter"
  )
})

test_that("prepare_ces_stan_data returns correct structure", {
  N <- 10
  sd <- prepare_ces_stan_data(
    prices       = runif(N, 1, 5),
    quantities   = runif(N, 10, 100),
    margins      = runif(N, 0.1, 0.5),
    market       = rep(1:2, each = 5),
    year         = rep(1, N),
    group        = sample(1:5, N, replace = TRUE),
    supply_model = 1L,
    use_cutoff   = FALSE,
    use_rho      = TRUE,
    use_hmt      = FALSE,
    priors       = list()
  )

  expect_true(is.list(sd))
  expect_equal(sd$N, N)
  expect_equal(sd$N_market, 2L)
  expect_equal(length(sd$revShare), N)
  expect_equal(length(sd$shareIn), N)
  expect_equal(length(sd$marginNormInv), N)
  expect_equal(length(sd$log_shareIn), N)
  expect_equal(length(sd$clr_shareIn), N)
  expect_equal(length(sd$log_marginNormInv), N)
  expect_equal(sd$log_shareIn, log(sd$shareIn))
  expect_equal(sd$log_marginNormInv, log(sd$marginNormInv))
  expect_equal(length(sd$log_price_std), N)
  expect_true(sd$log_price_sd > 0)
  expect_false("implied_s0" %in% names(sd))
  expect_equal(attr(sd, "model_spec_version"), 6L)
  expect_equal(attr(sd, "s0_prior_reference"), 0.2, tolerance = 1e-8)
  # Revenue shares should sum to ~1 within each market
  for (m in 1:sd$N_market_year) {
    idx <- which(attr(sd, "market_year_idx") == m)
    expect_equal(sum(sd$revShare[idx]), 1.0, tolerance = 1e-10)
  }
})

test_that("prepare_ces_stan_data computes revenue shares correctly", {
  N <- 4
  prices <- c(2, 3, 4, 5)
  quantities <- c(10, 20, 30, 40)
  sd <- prepare_ces_stan_data(
    prices       = prices,
    quantities   = quantities,
    margins      = rep(0.3, N),
    market       = rep(1, N),
    year         = rep(1, N),
    group        = 1:4,
    supply_model = 1L,
    use_cutoff   = FALSE,
    use_rho      = TRUE,
    use_hmt      = FALSE,
    priors       = list()
  )

  revenue <- prices * quantities
  expected_rev_shares <- revenue / sum(revenue)
  expect_equal(sd$revShare, expected_rev_shares)
  expect_equal(sd$shareIn, quantities / sum(quantities))
  expect_equal(sd$marginNormInv, rep(1 / 0.3, N))
})

test_that("CES uses the data-independent standardized prior", {
  N <- 6
  prices <- seq(1, 3, length.out = N)
  sd <- prepare_ces_stan_data(
    prices       = prices,
    quantities   = runif(N, 10, 100),
    margins      = rep(0.3, N),
    market       = rep(1, N),
    year         = rep(1, N),
    group        = 1:6,
    supply_model = 2L,
    use_cutoff   = FALSE,
    use_rho      = TRUE,
    use_hmt      = FALSE,
    priors       = list()
  )

  expect_equal(sd$prior_alpha_mean, 0)
  expect_equal(sd$prior_alpha_sd, 1)
})

test_that("CES FOC priors use unconditional revenue shares", {
  rev_share <- c(0.10, 0.25, 0.40)
  s0 <- 0.20
  r_unc <- rev_share * (1 - s0)
  gamma <- 4

  margins_by_model <- list(
    1 / (gamma - (gamma - 1) * r_unc),
    rep(1 / gamma, length(r_unc)),
    1 / gamma + (gamma - 1) * r_unc /
      (gamma * (1 + (gamma - 1) * s0)),
    1 - (1 - r_unc)^(1 / (gamma - 1))
  )

  for (supply_model in seq_along(margins_by_model)) {
    prior <- compute_gamma_prior(
      rev_share,
      margins_by_model[[supply_model]],
      s0,
      rep(1L, length(rev_share)),
      supply_model,
      log_price_sd = 1
    )
    expect_equal(prior$mean, log(gamma - 1), tolerance = 1e-8)
  }
})
