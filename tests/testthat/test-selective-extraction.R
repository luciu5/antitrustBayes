## Unit coverage for the posterior draw contract.  These tests use a small
## estimated-fit shell and mock only rstan::extract; no Stan sampling occurs.

.selective_fit <- function(quality = "sampled", structure = "simultaneous") {
  data <- data.frame(
    prices = c(1.1, 1.2, 1.3), quantities = c(40, 30, 20),
    margins = c(.5, .4, .3), market = rep("m1", 3L), year = rep("1", 3L),
    group = c("a", "b", "c"), outside_price = 1,
    stringsAsFactors = FALSE
  )
  stan_data <- list(quality_parameterization = if (quality == "sampled") 1L else 2L,
                    N_group = 3L, group = 1:3)
  attr(stan_data, "group_levels") <- c("a", "b", "c")
  attr(stan_data, "quality_parameterization") <- quality
  methods::new(
    "BayesLogit", stanfit = NULL, stan_data = stan_data, data = data,
    supply_model = 1L, demand_model = "logit", market_side = "output",
    price_sd = 1, call = quote(fake_selective_fit()),
    metadata = list(quality_parameterization = quality,
                    posterior_draws = 2L, structure = structure)
  )
}

.selective_draws <- function(s0_name = "s0_prob") {
  out <- list(
    a_market = matrix(c(.8, .9), nrow = 2L),
    share_mean = matrix(rep(c(-.8, -.9, -1.0), 2L), nrow = 2L, byrow = TRUE),
    b_quality_all = rbind(c(.1, -.1, 0), c(-.2, .1, .1))
  )
  out[[s0_name]] <- matrix(c(.2, .25), nrow = 2L)
  out
}

if (!methods::isClass("selective_contract_stanfit")) {
  methods::setClass("selective_contract_stanfit",
                    slots = c(sim = "list", model_pars = "character"))
}

.selective_fit_with_schema <- function(pars_oi,
                                       model_pars = c("a_market", "s0_prob")) {
  fit <- .selective_fit()
  fit@stanfit <- methods::new(
    "selective_contract_stanfit",
    sim = list(pars_oi = as.character(pars_oi)),
    model_pars = as.character(model_pars)
  )
  fit
}

test_that("policy plans use one selected union and normalize s0 aliases", {
  fit <- .selective_fit()
  draws <- .selective_draws("s0")
  n_extract <- 0L
  requested <- NULL
  testthat::with_mocked_bindings(
    extract = function(object, pars = NULL, ...) {
      n_extract <<- n_extract + 1L
      requested <<- pars
      draws
    }, .package = "rstan", {
      plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 7L)
    }
  )
  expect_equal(n_extract, 1L)
  ## Without a model_pars registry the shell receives the canonical alias;
  ## the mocked response uses `s0` and is normalized back to `s0_prob`.
  expect_equal(requested, c("a_market", "s0_prob", "share_mean", "b_quality_all"))
  expect_s3_class(plan, "bayes_policy_draw_plan")
  expect_equal(plan$fit_signature$fit_signature_version, 2L)
  expect_equal(plan$fit_signature$posterior_draws, 2L)
})

test_that("missing active quality fields fail before a plan is constructed", {
  fit <- .selective_fit(quality = "sampled")
  draws <- .selective_draws("s0")
  draws$b_quality_all <- NULL
  testthat::with_mocked_bindings(
    extract = function(object, pars = NULL, ...) draws,
    .package = "rstan", {
      expect_error(posterior_policy_draws(fit, n_draws = 2L),
                   "missing required field")
    }
  )
})

test_that("one-group sampled quality retains the zero-centered fallback", {
  fit <- .selective_fit(quality = "sampled")
  fit@metadata$quality_parameterization <- NULL
  attr(fit@stan_data, "quality_parameterization") <- NULL
  fit@stan_data$quality_parameterization <- as.numeric(1L)
  fit@stan_data$N_group <- 1L
  fit@stan_data$group <- rep(1L, 3L)
  attr(fit@stan_data, "group_levels") <- "a"
  draws <- .selective_draws("s0")
  draws$b_quality_all <- NULL
  requested <- NULL
  testthat::with_mocked_bindings(
    extract = function(object, pars = NULL, ...) {
      requested <<- pars
      draws
    }, .package = "rstan", {
      plan <- posterior_policy_draws(fit, n_draws = 2L)
    }
  )
  expect_false("b_quality_all" %in% requested)
  expect_equal(unname(plan$quality_draws), matrix(0, nrow = 2L, ncol = 1L))
  expect_equal(rownames(plan$quality_draws), as.character(plan$draw_indices))
})

test_that("active rho quality conditioning extracts its union once", {
  fit <- .selective_fit(quality = "collapsed")
  fit@metadata$share_likelihood <- "gaussian"
  fit@stan_data$share_likelihood <- 2L
  fit@stan_data$likelihood_mode <- 1L
  fit@stan_data$use_rho <- 1L
  fit@stan_data$log_shareIn <- c(-.8, -.9, -1)
  fit@stan_data$margin_idx <- 1:3
  fit@stan_data$log_marginNormInv <- log(c(.4, .5, .6))
  draws <- list(
    share_mean = matrix(rep(c(-.8, -.9, -1), 2L), nrow = 2L, byrow = TRUE),
    sigma_b_strat = c(.5, .5), sigma_share = c(.2, .2),
    rho = c(.2, .2), sigma_margin = c(.3, .3),
    pred_marginNormInv = matrix(rep(c(.4, .5, .6), 2L), nrow = 2L, byrow = TRUE)
  )
  n_extract <- 0L
  requested <- NULL
  testthat::with_mocked_bindings(
    extract = function(object, pars = NULL, ...) {
      n_extract <<- n_extract + 1L
      requested <<- pars
      draws
    }, .package = "rstan", {
      got <- posterior_quality(fit, seed = 9L)
    }
  )
  expect_equal(n_extract, 1L)
  expect_equal(requested, c("share_mean", "sigma_b_strat", "sigma_share",
                            "rho", "sigma_margin", "pred_marginNormInv"))
  expect_equal(nrow(got$draws), 2L)
  bad <- draws
  bad$pred_marginNormInv <- NULL
  testthat::with_mocked_bindings(
    extract = function(object, pars = NULL, ...) bad,
    .package = "rstan", {
      expect_error(posterior_quality(fit, seed = 9L), "missing required field")
    }
  )
})

test_that("numeric Stan quality codes work without metadata while contradictions remain rejected", {
  sampled <- .selective_fit(quality = "sampled")
  sampled@metadata$quality_parameterization <- NULL
  attr(sampled@stan_data, "quality_parameterization") <- NULL
  sampled@stan_data$quality_parameterization <- as.numeric(1L)
  expect_identical(antitrustBayes:::.posterior_contract_quality(sampled), "sampled")
  expect_identical(antitrustBayes:::.lq_quality_mode(sampled), "sampled")

  collapsed <- .selective_fit(quality = "collapsed")
  collapsed@metadata$quality_parameterization <- NULL
  attr(collapsed@stan_data, "quality_parameterization") <- NULL
  collapsed@stan_data$quality_parameterization <- as.numeric(2L)
  expect_identical(antitrustBayes:::.posterior_contract_quality(collapsed), "collapsed")
  expect_identical(antitrustBayes:::.lq_quality_mode(collapsed), "collapsed")

  contradictory <- collapsed
  contradictory@metadata$quality_parameterization <- "sampled"
  contradictory@stan_data$group <- 1:3
  contradictory@stan_data$log_shareIn <- c(-.8, -.9, -1)
  contradictory@stan_data$likelihood_mode <- 2L
  contradictory@metadata$posterior_draws <- NULL
  draws <- list(share_mean = matrix(c(-.8, -.9, -1), nrow = 1L),
                sigma_b_strat = .5, sigma_share = .4)
  expect_error(posterior_quality(contradictory, draws = draws),
               "inconsistent quality parameterization")
})

test_that("s0 logit aliases become probabilities without changing native probabilities", {
  fit <- .selective_fit()
  raw_logit <- matrix(c(-20, -2, 0, 2, 20, 1), nrow = 2L, ncol = 3L,
                      byrow = TRUE)
  raw <- list(a_market = matrix(c(.8, .9), nrow = 2L), s0 = raw_logit)
  testthat::with_mocked_bindings(
    extract = function(object, pars = NULL, ...) raw,
    .package = "rstan", {
      got <- antitrustBayes:::.posterior_extract(fit, "simulation")
    }
  )
  expect_equal(got$s0_prob, stats::plogis(raw_logit), tolerance = 0)
  expect_identical(dim(got$s0_prob), dim(raw_logit))

  raw_matrix_prob <- matrix(c(.01, .2, .8, .99, .4, .6), nrow = 2L,
                            ncol = 3L, byrow = TRUE)
  raw_native <- list(a_market = matrix(c(.8, .9), nrow = 2L),
                     s0_prob = raw_matrix_prob,
                     s0 = raw_logit)
  testthat::with_mocked_bindings(
    extract = function(object, pars = NULL, ...) raw_native,
    .package = "rstan", {
      got_native <- antitrustBayes:::.posterior_extract(fit, "simulation")
    }
  )
  expect_identical(got_native$s0_prob, raw_matrix_prob)

  fit_logit <- .selective_fit_with_schema(c("a_market", "s0_logit"),
                                           model_pars = c("a_market", "s0_prob"))
  requested <- NULL
  testthat::with_mocked_bindings(
    extract = function(object, pars = NULL, ...) {
      requested <<- pars
      list(a_market = matrix(c(.8, .9), nrow = 2L), s0_logit = raw_logit)
    },
    .package = "rstan", {
      got_logit <- antitrustBayes:::.posterior_extract(fit_logit, "simulation")
    }
  )
  expect_equal(requested, c("a_market", "s0_logit"))
  expect_equal(got_logit$s0_prob, stats::plogis(raw_logit), tolerance = 0)
})

test_that("required NULL fields fail while absent known optional fields are skipped", {
  fit <- .selective_fit()
  bad <- .selective_draws("s0_prob")
  bad$a_market <- NULL
  testthat::with_mocked_bindings(
    extract = function(object, pars = NULL, ...) bad,
    .package = "rstan", {
      expect_error(antitrustBayes:::.posterior_extract(fit, "simulation"),
                   "missing required field")
    }
  )

  ordered <- .selective_fit_with_schema(
    c("posterior_component_prob", "prior_component_prob",
      "firm_leader_prob", "s0_prob"),
    model_pars = c("posterior_component_prob", "prior_component_prob",
                   "firm_leader_prob", "s0_prob", "prior_firm_leader_prob")
  )
  ordered@metadata$structure <- "stackelberg"
  ordered@metadata$partition_prior <- "ordered_k"
  contract <- antitrustBayes:::.posterior_contract(ordered, "cutoff_summary")
  expect_false("prior_firm_leader_prob" %in% contract$canonical)
  expect_equal(contract$actual, c("posterior_component_prob",
                                  "prior_component_prob", "firm_leader_prob",
                                  "s0_prob"))

  stack_cutoff <- .selective_fit_with_schema(
    c("posterior_component_prob", "prior_component_prob",
      "firm_leader_prob", "cutoff_share", "s0_prob"),
    model_pars = c("posterior_component_prob", "prior_component_prob",
                   "firm_leader_prob", "cutoff_share", "s0_prob",
                   "prior_firm_leader_prob")
  )
  stack_cutoff@metadata$structure <- "stackelberg"
  stack_cutoff@metadata$partition_prior <- "cutoff"
  cutoff_contract <- antitrustBayes:::.posterior_contract(stack_cutoff, "cutoff_summary")
  expect_equal(cutoff_contract$actual, c("posterior_component_prob",
                                         "prior_component_prob",
                                         "firm_leader_prob", "cutoff_share",
                                         "s0_prob"))
})

test_that("retained pars schema takes precedence over compiled model pars", {
  fit <- .selective_fit_with_schema(c("a_market", "s0_logit"),
                                    model_pars = c("a_market", "s0_prob"))
  contract <- antitrustBayes:::.posterior_contract(fit, "simulation")
  expect_equal(contract$actual, c("a_market", "s0_logit"))
  expect_equal(contract$canonical, c("a_market", "s0_prob"))
})

test_that("draw count follows retained dimensions and prioritizes a_market", {
  fit <- .selective_fit()
  fit@metadata$posterior_draws <- 4L
  draws <- list(
    share_mean = matrix(seq_len(12L), nrow = 4L),
    a_market = matrix(seq_len(6L), nrow = 2L)
  )
  expect_equal(antitrustBayes:::.posterior_draw_count_from_draws(draws, fit), 2L)
})

test_that("v2 canonical fingerprints ignore selected-field order and aliases", {
  fit <- .selective_fit()
  draws <- .selective_draws("s0")
  canonical <- list(
    a_market = draws$a_market, s0_prob = draws$s0,
    share_mean = draws$share_mean, b_quality_all = draws$b_quality_all
  )
  canonical$s0_prob <- stats::plogis(canonical$s0_prob)
  expect_identical(
    antitrustBayes:::.lq_fit_fingerprint(fit, draws),
    antitrustBayes:::.lq_fit_fingerprint(fit, canonical)
  )
  changed <- canonical
  changed$share_mean[1L, 1L] <- changed$share_mean[1L, 1L] + 1e-3
  expect_false(identical(
    antitrustBayes:::.lq_fit_fingerprint(fit, canonical),
    antitrustBayes:::.lq_fit_fingerprint(fit, changed)
  ))
})
