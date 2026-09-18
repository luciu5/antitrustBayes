.stack_comparison_fixture <- function() {
  methods::new("BayesLogit",stanfit=NULL,supply_model=1L,demand_model="logit",
    market_side="output",price_sd=1,call=call("bayes_est"),
    metadata=list(structure="stackelberg", partition_prior="ordered_k", quality_scope="global",
      prior_normalization=list(
      mode="native_stan_normalized", model_spec_version=10L)),
    stan_data=list(structure=2L,partition_prior_code=1L,
      likelihood_mode=1L,use_hmt=0L,
      use_s0_market_pool=0L,use_s0_market_year_pool=0L,
      market_year_idx=c(1L,1L,2L,2L),mky_to_partition_unit=c(1L,1L),
      k_scope_code=1L,cluster=c(1L,1L,1L,1L)))
}
test_that("predictive units retain the marginalized equilibrium count", {
  fit <- .stack_comparison_fixture()
  expect_error(.bayes_comparison_units(list(fit),"observation"),"Observation LOO")
  expect_error(.bayes_comparison_units(list(fit),"period"),"Period LOO")
  expect_invisible(.bayes_comparison_units(list(fit),"market"))
  expect_invisible(.bayes_comparison_units(list(fit),"cluster"))
  fit@stan_data$cluster <- c(1L,1L,2L,2L)
  expect_error(.bayes_comparison_units(list(fit),"cluster"),"whole ordered-k")
  fit@stan_data$k_scope_code <- 0L
  fit@stan_data$mky_to_partition_unit <- c(1L,2L)
  expect_invisible(.bayes_comparison_units(list(fit),"period"))
  expect_invisible(.bayes_comparison_units(list(fit),"cluster"))
  fit@stan_data$cluster <- c(1L,2L,3L,3L)
  expect_error(.bayes_comparison_units(list(fit),"cluster"),"whole ordered-k")
})

test_that("exact fringe retains one k at the whole-market predictive unit", {
  fit <- .stack_comparison_fixture()
  fit@metadata$structure <- "fringe"
  fit@stan_data$structure <- 1L
  expect_false(.bayes_is_stackelberg(fit))
  expect_true(.bayes_is_exact_ordered_k(fit))
  expect_error(.bayes_comparison_units(list(fit), "observation"),
               "whole-market ordered-k")
  expect_error(.bayes_comparison_units(list(fit), "period"),
               "one ordered k")
  expect_invisible(.bayes_comparison_units(list(fit), "market"))
})
test_that("evidence rejects non-normalized factors and improper pools", {
  fit <- .stack_comparison_fixture()
  expect_invisible(.bayes_validate_bridge(list(fit)))
  fit@stan_data$use_hmt <- 1L
  expect_error(.bayes_validate_bridge(list(fit)),"HMT")
  fit@stan_data$likelihood_mode <- 3L
  expect_error(.bayes_validate_bridge(list(fit)),"HMT")
  fit@stan_data$likelihood_mode <- 2L
  expect_invisible(.bayes_validate_bridge(list(fit)))
  fit@stan_data$use_s0_market_pool <- fit@stan_data$use_s0_market_year_pool <- 1L
  expect_error(.bayes_validate_bridge(list(fit)),"without a prior")
})

.stack_prior_contract_fixture <- function() {
  methods::new("BayesLogit", stanfit = NULL,
    supply_model = 1L, demand_model = "logit", market_side = "output",
    price_sd = 1, call = call("bayes_est"),
    data = data.frame(prices = c(1, 2), quantities = c(4, 6),
      margins = c(.25, .2), market = "m", year = "1", group = c("a", "b"),
      loo_cluster = "market", margin_observed = TRUE, outside_price = 1),
    stan_data = list(shareIn = c(.4, .6), clr_shareIn = c(-.2, .2),
      marginNormInv = c(4, 5), log_marginNormInv = log(c(4, 5)),
      margin_obs = c(1L, 1L), market_year_idx = c(1L, 1L),
      cluster = c(1L, 1L), likelihood_mode = 1L),
    metadata = list(prior_manifest = "baseline",
      resolved_priors = bayes_prior_manifest("baseline")))
}

test_that("comparison ignores hyperparameters of absent role priors", {
  reference <- .stack_prior_contract_fixture()
  reference@stan_data$use_cutoff <- 1L
  reference@stan_data$structure <- 1L
  direct <- .stack_prior_contract_fixture()
  direct@metadata$structure <- "stackelberg"
  direct@stan_data$structure <- 2L
  direct@stan_data$partition_prior_code <- 1L
  direct@metadata$resolved_priors$prior_cutoff_alpha <- 2
  direct@metadata$resolved_priors$prior_cutoff_beta <- 8
  direct@metadata$resolved_priors$prior_sigma_meanval_fringe <- .8
  expect_invisible(validate_comparison_contract(list(reference, direct)))
  direct@metadata$resolved_priors$prior_alpha_sd <- .7
  expect_error(validate_comparison_contract(list(reference, direct)), "same prior manifest")
  direct@metadata$resolved_priors$prior_alpha_sd <- reference@metadata$resolved_priors$prior_alpha_sd
  direct@metadata$structure <- "fringe"
  direct@metadata$partition_prior <- "cutoff"
  direct@stan_data$structure <- 1L
  direct@stan_data$partition_prior_code <- 0L
  direct@stan_data$use_cutoff <- 1L
  expect_error(validate_comparison_contract(list(reference, direct)), "same prior manifest")
})


test_that("evidence requires verified native normalized fits", {
  fit <- .stack_comparison_fixture()
  expect_invisible(.bayes_validate_bridge(list(fit)))
  fit@metadata$prior_normalization$model_spec_version <- 6L
  expect_error(.bayes_validate_bridge(list(fit)), "Refit")
  fit@metadata$prior_normalization <- NULL
  expect_error(.bayes_validate_bridge(list(fit)), "Refit")
  fit@metadata$prior_normalization <- list(mode = "supplied_draws", model_spec_version = 8L)
  expect_error(.bayes_validate_bridge(list(fit)), "native normalized")
})
