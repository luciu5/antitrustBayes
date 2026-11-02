## These are identical-state tests, not evidence of independent sampler recovery.
## Scrambled Stan columns exercise both adapters before real economic solvers run.
.backend_policy_pair <- function(conduct, structure = "simultaneous") {
  data <- data.frame(prices = c(1.1, 1.2, 1.3), quantities = c(40, 30, 20),
    margins = c(.5, .4, .3), market = "m1", year = "t1",
    group = c("a", "b", "c"), outside_price = 1, margin_observed = TRUE)
  stan_data <- prepare_stan_data(share_basis = "absolute", data$prices, data$quantities, data$margins,
    data$market, data$year, data$group, supply_model = conduct,
    market_side = "output", use_rho = FALSE, structure = structure,
    partition_prior = if (structure == "simultaneous") "cutoff" else "ordered_k",
    k_support = if (structure == "simultaneous") NULL else 1L,
    k_prior = if (structure == "simultaneous") NULL else 1,
    quality_parameterization = "sampled")
  columns <- c("share_mean[3]", "b_quality_all[2]", "a_market[1]",
    "s0_prob[1]", "b_quality_all[1]", "share_mean[1]",
    "b_quality_all[3]", "share_mean[2]", "posterior_component_prob[1]",
    "prior_component_prob[1]")
  values <- array(0, c(2L, 1L, length(columns)),
    dimnames = list(iteration = 1:2, chain = "1", variable = columns))
  values[, , "a_market[1]"] <- c(.6, .7)
  values[, , "s0_prob[1]"] <- c(.2, .25)
  quality <- rbind(c(.1, -.1, 0), c(-.2, .1, .1))
  for (j in 1:3) {
    values[, , paste0("b_quality_all[", j, "]")] <- quality[, j]
    values[, , paste0("share_mean[", j, "]")] <-
      log(data$quantities[[j]] / sum(data$quantities)) + quality[, j]
  }
  values[, , "posterior_component_prob[1]"] <- 1
  values[, , "prior_component_prob[1]"] <- 1
  native <- structure(list(draws_array = values), class = "NutpieFit")
  cmdstan <- structure(list(
    metadata = function() list(stan_variables = unique(sub("\\[.*$", "", columns))),
    draws = function(variables = NULL, format = "draws_array")
      .bayes_backend_draws(native, variables, format)
  ), class = c("CmdStanMCMC", "CmdStanFit"))
  fit <- methods::new("BayesLogit", stanfit = cmdstan, stan_data = stan_data,
    data = data, supply_model = as.integer(conduct), demand_model = "logit",
    market_side = "output", price_sd = stats::sd(data$prices),
    metadata = list(backend = "cmdstanr", inference = "hmc",
      structure = structure, partition_prior = attr(stan_data, "partition_prior"),
      firm_market_ranking = attr(stan_data, "firm_market_ranking"),
      k_support = if (structure == "simultaneous") NULL else 1L,
      k_prior = if (structure == "simultaneous") NULL else 1,
      quality_parameterization = "sampled", quality_scope = "market",
      model_spec_version = .antitrust_bayes_model_spec_version))
  other <- fit
  other@stanfit <- native
  other@metadata$backend <- "nutpie"
  list(cmdstanr = fit, nutpie = other)
}

test_that("identical backend states give the same quality, roles and merger effects", {
  for (candidate in c("MC", "B", "C", "BF", "CF", "BS", "CS")) {
    structure <- if (candidate %in% c("BF", "CF")) "fringe" else
      if (candidate %in% c("BS", "CS")) "stackelberg" else "simultaneous"
    conduct <- if (candidate == "MC") 2L else
      if (candidate %in% c("C", "CF", "CS")) 3L else 1L
    pair <- .backend_policy_pair(conduct, structure)
    plans <- lapply(pair, posterior_policy_draws, n_draws = 2L, seed = 29L)
    expect_identical(plans[[1]]$draw_indices, plans[[2]]$draw_indices, info = candidate)
    expect_identical(plans[[1]]$quality_draws, plans[[2]]$quality_draws, info = candidate)
    expect_identical(plans[[1]]$k_cache, plans[[2]]$k_cache, info = candidate)
    for (anchoring in c("observed", "latent_quality")) {
      prepared <- lapply(seq_along(pair), function(i) suppressWarnings(
        posterior_policy_models(pair[[i]], draw_plan = plans[[i]],
          demand_anchoring = anchoring)))
      for (i in 1:2) expect_equal(nrow(prepared[[i]]$failures), 0L, info = candidate)
      for (draw in 1:2) {
        models <- lapply(prepared, function(x) {
          model <- x$models[[1]][[draw]]
          if (methods::is(model, "AntitrustFit")) model@model else model
        })
        expect_equal(models[[1]]@pricePre, models[[2]]@pricePre, info = candidate)
        expect_equal(antitrust::calcMC(models[[1]], TRUE),
          antitrust::calcMC(models[[2]], TRUE), info = candidate)
        expect_equal(antitrust::calcShares(models[[1]], TRUE),
          antitrust::calcShares(models[[2]], TRUE), info = candidate)
      }
      sims <- lapply(seq_along(pair), function(i) suppressWarnings(
        simulate(pair[[i]], c("a", "a", "c"), draw_plan = plans[[i]],
          demand_anchoring = anchoring, mc.cores = 1L)))
      expect_equal(as.integer(sims[[1]]@n_valid), 2L, info = candidate)
      expect_identical(sims[[1]]@n_valid, sims[[2]]@n_valid, info = candidate)
      expect_equal(sims[[1]]@price_post, sims[[2]]@price_post, info = candidate)
    }
  }
})

test_that("Nutpie conditional-k scheduling forwards backend and preserves failed support", {
  local_mocked_bindings(bayes_est = function(...) {
    args <- list(...)
    if (args$k_support == 1L) stop("retained Nutpie component failure")
    args
  }, .package = "antitrustBayes")
  prior <- c(.2, .3, .5)
  result <- bayes_est_conditional_k(market = rep("m", 3), group = letters[1:3],
    structure = "fringe", supply_model = "cournot", backend = "nutpie",
    nutpie_options = list(adaptation = "low_rank"), k_support = 0:2,
    k_prior = prior)
  expect_identical(result$k_support, 0:2)
  expect_equal(unname(result$k_prior), prior)
  expect_identical(names(result$fits), as.character(0:2))
  expect_s3_class(result$fits[[2]], "error")
  for (i in c(1L, 3L)) {
    expect_identical(result$fits[[i]]$backend, "nutpie")
    expect_identical(result$fits[[i]]$nutpie_options, list(adaptation = "low_rank"))
    expect_identical(result$fits[[i]]$k_prior, 1)
  }
  expect_error(bayes_combine_conditional_k(c(`0` = -1, `2` = -2),
    k_support = result$k_support, k_prior = result$k_prior), "Every original")
})
