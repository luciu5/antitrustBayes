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
  attr(current_data, "model_spec_version") <- 13L
  attr(current_data, "quality_scope") <- "global"
  current <- stale
  current@stan_data <- current_data
  expect_invisible(validate_fit_spec(current))

  old_normalization <- current
  attr(old_normalization@stan_data, "model_spec_version") <- 12L
  expect_error(validate_fit_spec(old_normalization),
               "current specification is 13.*Refit it with bayes_est")
})
