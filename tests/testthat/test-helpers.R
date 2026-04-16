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

test_that("validate_inputs catches bad margins", {
  expect_error(
    validate_inputs(
      c(1, 2, 3), c(10, 20, 30), c(-0.1, 0.2, 0.3),
      c(1, 1, 1), c(1, 1, 1), c(1, 2, 3)
    ),
    "'margins' \\(Lerner index \\(p-c\\)/p\\) must be positive and finite"
  )
  expect_error(
    validate_inputs(
      c(1, 2, 3), c(10, 20, 30), c(0.2, Inf, 0.3),
      c(1, 1, 1), c(1, 1, 1), c(1, 2, 3)
    ),
    "'margins' \\(Lerner index \\(p-c\\)/p\\) must be positive and finite"
  )
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
  expect_equal(length(sd$marginInv), N)
  expect_equal(length(sd$price_std), N)
  expect_true(sd$price_sd > 0)
  expect_true(all(sd$implied_s0 > 0 & sd$implied_s0 < 1))
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
  expect_equal(length(sd$lernerInv), N)
  expect_equal(length(sd$log_price_std), N)
  expect_true(sd$log_price_sd > 0)
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
  # lernerInv = 1/L
  expect_equal(sd$lernerInv, rep(1 / 0.3, N))
})
