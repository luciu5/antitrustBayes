test_that("posterior_k aggregates exact Stan components by mapped market", {
  stan_data <- list(
    N_partition_unit = 2L,
    partition_unit_market = c(2L, 1L),
    partition_unit_start = c(1L, 3L),
    partition_unit_count = c(2L, 2L),
    partition_component_unit = c(1L, 1L, 2L, 2L),
    # The support is intentionally not sorted in component order.
    partition_component_k = c(3L, 1L, 2L, 1L)
  )
  attr(stan_data, "market_levels") <- c("m1", "m2")
  fit <- methods::new(
    "BayesLogit", stanfit = NULL,
    data = data.frame(market = c("m1", "m2")),
    stan_data = stan_data, price_sd = 1,
    supply_model = 1L, demand_model = "logit", market_side = "input",
    call = quote(fake()),
    metadata = list(
      structure = "stackelberg", partition_prior = "ordered_k",
      posterior_draws = 2L
    )
  )
  draws <- list(
    posterior_component_prob = rbind(
      c(.2, .8, .6, .4),
      c(.4, .6, .8, .2)
    ),
    prior_component_prob = rbind(
      c(.5, .5, .25, .75),
      c(.5, .5, .25, .75)
    )
  )

  got <- posterior_k(fit, draws = draws)
  expect_s3_class(got, "posterior_k_summary")
  expect_identical(names(got)[1:5],
                   c("market", "k", "prior", "posterior", "role"))
  expect_identical(got$market, c("m2", "m2", "m1", "m1"))
  expect_identical(got$k, c(1L, 3L, 1L, 2L))
  expect_equal(got$prior, c(.5, .5, .75, .25))
  expect_equal(got$posterior, c(.7, .3, .3, .7))
  expect_true(all(got$role == "leaders"))
  expect_equal(unique(got$posterior_mean[got$market == "m2"]), 1.6)
  expect_equal(unique(got$posterior_median[got$market == "m2"]), 1)
  expect_equal(unique(got$posterior_lower[got$market == "m2"]), 1)
  expect_equal(unique(got$posterior_upper[got$market == "m2"]), 3)
  expect_equal(attr(got, "summary")$market, c("m2", "m1"))
  expect_equal(attr(got, "metadata")$mapping$component_k,
               stan_data$partition_component_k)
})

test_that("posterior_k distinguishes fringe core states and specified priors", {
  market <- bayes_market(
    prices = c(1, 1.1, 1.2), quantities = c(5, 3, 2),
    margins = rep(.3, 3), group = c("A", "B", "C"),
    outside_price = 1, market_side = "input", insideSize = 10
  )
  fit <- bayes_specify(
    bayes_params(alpha = .7, s0 = .2, supply_model = 1L,
                 structure = "fringe", k_support = c(0L, 2L, 3L),
                 k_prior = c(.2, .3, .5)), market
  )
  got <- posterior_k(fit)
  expect_true(all(got$role == "strategic_core"))
  expect_equal(got$k, c(0L, 2L, 3L))
  expect_equal(got$prior, c(.2, .3, .5))
  expect_equal(got$posterior, c(.2, .3, .5))
  expect_equal(attr(got, "summary")$posterior_mean, 2.1)
})

test_that("posterior_k fails clearly on missing or malformed generated draws", {
  stan_data <- list(
    N_partition_unit = 1L, partition_unit_market = 1L,
    partition_unit_start = 1L, partition_unit_count = 2L,
    partition_component_unit = c(1L, 1L),
    partition_component_k = c(1L, 2L)
  )
  attr(stan_data, "market_levels") <- "m"
  fit <- methods::new(
    "BayesLogit", stanfit = NULL,
    data = data.frame(market = "m"), stan_data = stan_data, price_sd = 1,
    supply_model = 1L, demand_model = "logit", market_side = "input",
    call = quote(fake()),
    metadata = list(structure = "stackelberg", partition_prior = "ordered_k",
                    posterior_draws = 2L)
  )
  expect_error(posterior_k(fit, draws = list()), "missing generated quantity")
  expect_error(
    posterior_k(fit, draws = list(
      posterior_component_prob = matrix(c(.5, .5), nrow = 1L),
      prior_component_prob = matrix(c(.5, .5), nrow = 1L)
    )),
    "has 1 draws; the fitted object records 2"
  )
  expect_error(
    posterior_k(fit, draws = list(
      posterior_component_prob = matrix(c(.5, NA, .5, .5), nrow = 2L),
      prior_component_prob = matrix(rep(c(.5, .5), 2), nrow = 2L)
    )),
    "contains invalid probabilities"
  )
})
