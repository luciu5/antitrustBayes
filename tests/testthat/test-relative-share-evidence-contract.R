.relative_evidence_fixture <- function(basis = "relative", k = 1L,
                                       structure = "fringe", backend = "nutpie") {
  data <- data.frame(prices = c(1.1, 1.3, 1.6), quantities = c(50, 30, 20),
    margins = c(.2, .25, .3), market = "m", year = "t", group = letters[1:3],
    loo_cluster = "m", margin_observed = TRUE, outside_price = 1)
  d <- prepare_stan_data(data$prices, data$quantities, data$margins,
    data$market, data$year, data$group, supply_model = 3L,
    market_side = "output", structure = structure,
    partition_prior = "ordered_k", k_support = k, k_prior = 1, share_basis = basis)
  methods::new("BayesLogit", data = data, stan_data = d, supply_model = 3L,
    demand_model = "logit", market_side = "output", price_sd = d$price_sd,
    metadata = list(backend = backend, share_basis = basis, structure = structure,
      partition_prior = "ordered_k", k_support = k, k_prior = 1,
      resolved_priors = bayes_prior_manifest(), prior_manifest = "baseline"))
}

test_that("evidence compares samplers within one share observation measure", {
  a <- .relative_evidence_fixture(backend = "cmdstanr")
  b <- .relative_evidence_fixture(backend = "nutpie")
  expect_invisible(validate_comparison_contract(list(a, b)))
  expect_error(validate_comparison_contract(list(a, .relative_evidence_fixture("absolute"))),
               "same share_basis")
  expect_error(.bayes_comparison_units(list(a), "observation"), "Observation LOO")
  expect_invisible(.bayes_comparison_units(list(a), "market"))
  expect_error(.bayes_comparison_units(list(a), "period"), "one ordered k")
})

test_that("conditional evidence keeps one family and basis while permitting mixed samplers", {
  fits <- list("1" = .relative_evidence_fixture(k = 1L, backend = "cmdstanr"),
               "2" = .relative_evidence_fixture(k = 2L))
  logml <- c("1" = -11, "2" = -12)
  result <- bayes_combine_conditional_k(logml, 1:2, k_prior = c(.25, .75),
                                       component_fits = fits)
  expect_identical(result$share_basis, "relative")
  expect_identical(unname(result$component_backends), c("cmdstanr", "nutpie"))
  expect_equal(exp(result$logml), sum(exp(logml) * c(.25, .75)))
  wrong <- fits; wrong[[2L]] <- .relative_evidence_fixture("absolute", k = 2L)
  expect_error(bayes_combine_conditional_k(logml, 1:2, component_fits = wrong), "same share_basis")
  wrong <- fits; wrong[[2L]] <- .relative_evidence_fixture(k = 2L, structure = "stackelberg")
  expect_error(bayes_combine_conditional_k(logml, 1:2, component_fits = wrong), "one conduct family")
  wrong <- fits; wrong[[2L]]@demand_model <- "ces"
  expect_error(bayes_combine_conditional_k(logml, 1:2, component_fits = wrong), "one conduct family")
  expect_error(bayes_combine_conditional_k(logml[-2L], 1:2, component_fits = fits[-2L]),
               "every original")
})
