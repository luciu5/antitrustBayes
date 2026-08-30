make_test_sim <- function(prefix, markets = c("m1", "m2"), n = 3L) {
  sims <- price_pre <- price_post <- price_change <- shares_pre <-
    shares_post <- product_labels <- failures <- list()
  for (mk in markets) {
    sims[[mk]] <- lapply(seq_len(n), function(i) paste(prefix, mk, i, sep = "-"))
    price_pre[[mk]] <- matrix(seq_len(n), ncol = 1,
                              dimnames = list(NULL, "p"))
    price_post[[mk]] <- price_pre[[mk]] + 1
    price_change[[mk]] <- price_post[[mk]] / price_pre[[mk]] - 1
    shares_pre[[mk]] <- matrix(0.2, nrow = n, ncol = 1,
                                dimnames = list(NULL, "p"))
    shares_post[[mk]] <- matrix(0.2, nrow = n, ncol = 1,
                                 dimnames = list(NULL, "p"))
    product_labels[[mk]] <- "p"
    failures[[mk]] <- rep(NA_character_, n)
  }
  new("BayesSim",
    sims = sims,
    price_pre = price_pre,
    price_post = price_post,
    price_change = price_change,
    shares_pre = shares_pre,
    shares_post = shares_post,
    market_ids = markets,
    product_labels = product_labels,
    n_draws = as.integer(n),
    n_valid = stats::setNames(rep(as.integer(n), length(markets)), markets),
    failures = failures,
    call = call("bayes_sim")
  )
}

test_that("bayes_avg returns exactly the requested allocation", {
  set.seed(1)
  result <- bayes_avg(
    list(make_test_sim("a"), make_test_sim("b")),
    weights = c(0.5, 0.5),
    n_draws = 7L
  )
  expect_equal(result@n_draws, 7L)
  expect_equal(unname(result@n_valid), c(7L, 7L))
  expect_equal(length(result@sims$m1), 7L)
  prefixes <- substring(unlist(result@sims$m1), 1, 1)
  expect_equal(as.integer(table(prefixes)), c(4L, 3L))
})

test_that("bayes_avg preserves zero model weights", {
  result <- bayes_avg(
    list(make_test_sim("a"), make_test_sim("b")),
    weights = c(1, 0),
    n_draws = 5L
  )
  expect_true(all(startsWith(unlist(result@sims$m1), "a-")))
  expect_equal(result@n_draws, 5L)
})

test_that("bayes_weights falls back cleanly after failed bridges", {
  comp <- data.frame(
    model = c("a", "b"),
    weight_stacking = c(0.25, 0.75),
    post_prob = c(NA_real_, NA_real_)
  )
  expect_message(
    weights <- bayes_weights(comp, "bridge"),
    "using stacking"
  )
  expect_equal(weights, c(a = 0.25, b = 0.75))
})

make_contract_fit <- function(demand_model = "logit") {
  stan_data <- list(
    N = 2L,
    shareIn = c(0.4, 0.6),
    log_shareIn = log(c(0.4, 0.6)),
    marginNormInv = c(4, 5),
    log_marginNormInv = log(c(4, 5)),
    margin_obs = c(1L, 1L)
  )
  attr(stan_data, "model_spec_version") <- 4L
  new("BayesLogit",
    stanfit = NULL,
    data = data.frame(
      prices = c(1, 2), quantities = c(4, 6), margins = c(0.25, 0.2),
      market = c("m", "m"), year = c("1", "1"), group = c("a", "b"),
      margin_observed = c(TRUE, TRUE), outside_price = c(1, 1)
    ),
    stan_data = stan_data,
    supply_model = 1L,
    demand_model = demand_model,
    market_side = "output",
    price_sd = 1,
    call = call("bayes_est")
  )
}

test_that("comparison contract permits common-outcome Logit and CES fits", {
  expect_invisible(validate_comparison_contract(list(
    make_contract_fit("logit"), make_contract_fit("ces")
  )))
})

test_that("comparison contract rejects different observations", {
  fit1 <- make_contract_fit()
  fit2 <- make_contract_fit()
  fit2@data$prices[1] <- 9
  expect_error(
    validate_comparison_contract(list(fit1, fit2)),
    "same ordered observations"
  )
})
