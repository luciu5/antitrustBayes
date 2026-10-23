.relative_api_data <- function(markets = 1L, periods = 1L) {
  d <- expand.grid(firm = seq_len(4), year = seq_len(periods),
                   market = seq_len(markets))
  list(prices = 1 + d$firm / 10 + d$year / 30,
       quantities = c(4, 3, 2, 1)[d$firm],
       margins = 0.2 + d$firm / 100,
       market = d$market, year = d$year, group = d$firm,
       supply_model = 1L, use_cutoff = FALSE, use_rho = TRUE,
       use_hmt = FALSE, priors = list())
}

test_that("relative quantity-share basis is explicit and defaults in both demands", {
  expect_identical(eval(formals(bayes_est)$share_basis), c("relative", "absolute"))
  expect_identical(eval(formals(bayes_est)$backend), c("cmdstanr", "nutpie"))
  for (markets in c(1L, 4L)) for (periods in c(1L, 4L)) {
    f <- .relative_api_data(markets, periods)
    for (prepare in list(prepare_stan_data, prepare_ces_stan_data)) {
      d <- do.call(prepare, f)
      expect_identical(d$share_basis, 1L)
      expect_identical(attr(d, "share_basis"), "relative")
      expect_identical(attr(d, "share_observation_measure"), "ilr_quantity_shares_v1")
      expect_equal(attr(d, "share_dimension"), 3 * markets * periods)
      expect_identical(attr(d, "market_intercept_parameterization_effective"),
                       "derived_market_year")
      expect_equal(d$N_quality, 4L * markets)
      expect_equal(d$N_market_year, markets * periods)
      absolute <- do.call(prepare, c(f, list(share_basis = "absolute")))
      expect_identical(absolute$share_basis, 2L)
      expect_identical(attr(absolute, "market_intercept_parameterization_effective"),
                       "sampled")
    }
  }
})

test_that("unsupported relative-share coordinates and errors fail explicitly", {
  f <- .relative_api_data()
  expect_error(do.call(prepare_stan_data, c(f, list(share_likelihood = "student_t"))),
               "requires Gaussian")
  expect_error(do.call(prepare_stan_data, c(f, list(margin_likelihood = "student_t"))),
               "requires Gaussian")
  expect_error(do.call(prepare_stan_data, c(f, list(market_intercept_parameterization = "sampled"))),
               "applies only")
  expect_error(do.call(prepare_stan_data, c(f, list(market_s0_parameterization = "joint"))),
               "require.*legacy")
  expect_error(.resolve_share_basis(priors = list(prior_sigma_b_market_scale = 1)),
               "priors apply only")
  expect_no_error(do.call(prepare_stan_data, c(f, list(use_market_resid_scale = TRUE))))
})

test_that("conditional-k call generation forwards scientific basis unchanged", {
  f <- .relative_api_data()
  f$supply_model <- "cournot"
  f$structure <- "fringe"
  f$share_basis <- "relative"
  calls <- do.call(bayes_conditional_k_arguments, c(f, list(k_support = 0:3)))
  expect_true(all(vapply(calls, function(x) identical(x$share_basis, "relative"), logical(1))))
  expect_identical(attr(calls, "k_support"), 0:3)
})
