test_that("BayesLogit class can be instantiated", {
  obj <- new("BayesLogit",
    stanfit      = NULL,
    data         = data.frame(x = 1),
    stan_data    = list(N = 1L),
    supply_model = 1L,
    demand_model = "logit",
    market_side  = "input",
    price_sd     = 1.0,
    call         = call("bayes_est")
  )
  expect_s4_class(obj, "BayesLogit")
})

test_that("BayesSim class can be instantiated", {
  obj <- new("BayesSim",
    sims           = list(),
    price_pre      = list(),
    price_post     = list(),
    price_change   = list(),
    shares_pre     = list(),
    shares_post    = list(),
    market_ids     = character(0),
    product_labels = list(),
    n_draws        = 0L,
    n_valid        = integer(0),
    failures       = list(),
    call           = call("bayes_sim")
  )
  expect_s4_class(obj, "BayesSim")
})

test_that("fit specification validation rejects stale fits", {
  stale <- new("BayesLogit",
    stanfit      = NULL,
    data         = data.frame(x = 1),
    stan_data    = list(N = 1L),
    supply_model = 1L,
    demand_model = "logit",
    market_side  = "input",
    price_sd     = 1.0,
    call         = call("bayes_est")
  )
  expect_error(validate_fit_spec(stale), "Refit it with bayes_est")

  current_data <- list(N = 1L)
  attr(current_data, "model_spec_version") <- 5L
  current <- stale
  current@stan_data <- current_data
  expect_invisible(validate_fit_spec(current))
})
