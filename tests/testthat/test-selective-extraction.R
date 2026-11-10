## Unit coverage for the posterior draw contract.  These tests use a small
## estimated-fit shell and mock only the CmdStanR draw adapter; no Stan sampling occurs.

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

test_that("selective total/fraction draws preserve economic quality across nested policy recovery", {
  ## Material failure: the adapter's economic SD used to lose its explicit
  ## alias, so nested recovery mistook it for the raw total scale and failed
  ## (or could condition quality on the wrong variance). Unlike a full-list
  ## mock, this backend double returns only requested Stan columns.
  fit <- .selective_fit(quality = "collapsed")
  fit@stan_data$N <- 3L
  fit@stan_data$log_shareIn <- c(-.8, -.9, -1)
  fit@stan_data$share_likelihood <- 2L
  fit@stan_data$likelihood_mode <- 2L
  fit@stan_data$quality_scale_parameterization <- 2L
  attr(fit@stan_data, "model_spec_version") <- .antitrust_bayes_model_spec_version
  attr(fit@stan_data, "quality_scope") <- "global"
  fit@metadata$quality_scale_parameterization_effective <- "total_fraction"
  fit@metadata$quality_scope <- "global"
  fit@metadata$share_likelihood <- "gaussian"
  fit@metadata$model_spec_version <- .antitrust_bayes_model_spec_version
  mu <- rbind(c(-.9, -.75, -.85), c(-.7, -.95, -1.2))
  tau <- c(.6, .3)
  sigma <- c(.4, .5)
  total <- sqrt(tau^2 + sigma^2)
  fraction_logit <- 2 * (log(sigma) - log(tau))
  columns <- c("sigma_b_strat", "share_mean[3]", "a_market[1]",
    "log_sigma_share", "s0_prob[1]", "sigma_b_strat_economic",
    "share_mean[1]", "sigma_share", "share_mean[2]")
  values <- array(0, c(2L, 1L, length(columns)),
    dimnames = list(iteration = 1:2, chain = "1", variable = columns))
  values[, , "sigma_b_strat"] <- total
  values[, , "log_sigma_share"] <- fraction_logit
  values[, , "sigma_b_strat_economic"] <- tau
  values[, , "sigma_share"] <- sigma
  values[, , "a_market[1]"] <- 6
  values[, , "s0_prob[1]"] <- .2
  for (j in 1:3) values[, , paste0("share_mean[", j, "]")] <- mu[, j]
  native <- structure(list(draws_array = values), class = "NutpieFit")
  requested <- list()
  fit@stanfit <- structure(list(
    metadata = function() list(stan_variables = unique(sub("\\[.*$", "", columns))),
    draws = function(variables = NULL, format = "draws_array") {
      requested[[length(requested) + 1L]] <<- variables
      .bayes_backend_draws(native, variables, format)
    }), class = c("CmdStanMCMC", "CmdStanFit"))

  extracted <- .posterior_extract(fit, "policy_plan")
  expect_equal(extracted$sigma_b_strat, tau)
  expect_equal(extracted$sigma_b_strat_economic, tau)
  expect_identical(.posterior_canonical_draws(extracted, fit), extracted)
  got <- posterior_quality(fit, seed = 47L)
  raw <- posterior_quality(fit, seed = 47L, draws = list(share_mean = mu,
    sigma_b_strat = total, log_sigma_share = fraction_logit, sigma_share = sigma))
  expect_equal(got$draws, raw$draws, tolerance = 1e-14)
  centered <- diag(3L) - matrix(1 / 3, 3L, 3L)
  for (i in 1:2) {
    expected_mean <- tau[i]^2 / (tau[i]^2 + sigma[i]^2) *
      as.numeric(centered %*% (fit@stan_data$log_shareIn - mu[i, ]))
    expected_sd <- sqrt(2 / 3 * tau[i]^2 * sigma[i]^2 / (tau[i]^2 + sigma[i]^2))
    expect_equal(as.numeric(got$conditional_means[i, ]), expected_mean, tolerance = 1e-12)
    expect_equal(as.numeric(got$conditional_sds[i, ]), rep(expected_sd, 3L), tolerance = 1e-12)
  }
  plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 47L)
  prepared <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan))
  result <- suppressWarnings(simulate(fit, fit@data$group, draw_plan = plan, mc.cores = 1L))
  expect_equal(nrow(prepared$failures), 0L)
  expect_equal(as.integer(result@n_valid), 2L)
  expect_gt(max(abs(plan$quality_draws)), 0)
  expect_equal(result@price_change$m1, matrix(0, 2L, 3L), ignore_attr = TRUE, tolerance = 1e-8)
  expect_true(all(vapply(requested, function(x) "sigma_b_strat_economic" %in% x,
                         logical(1))))
  expect_false(any(vapply(requested, function(x)
    any(c("sigma_b_strat", "log_sigma_share") %in% x), logical(1))))
})

test_that("policy plans use one selected union and normalize s0 aliases", {
  fit <- .selective_fit()
  draws <- .selective_draws("s0")
  n_extract <- 0L
  requested <- NULL
  testthat::with_mocked_bindings(
    .cmdstan_extract = function(fit, variables = NULL, ...) {
      n_extract <<- n_extract + 1L
      requested <<- variables
      draws
    }, .package = "antitrustBayes", {
      plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 7L)
    }
  )
  expect_equal(n_extract, 1L)
  ## Without a model_pars registry the shell receives the canonical alias;
  ## the mocked response uses `s0` and is normalized back to `s0_prob`.
  expect_equal(requested, c("a_market", "s0_prob", "share_mean", "b_quality_all"))
  expect_s3_class(plan, "bayes_policy_draw_plan")
  expect_equal(plan$fit_signature$fit_signature_version, 3L)
  expect_equal(plan$fit_signature$posterior_draws, 2L)
})

test_that("missing active quality fields fail before a plan is constructed", {
  fit <- .selective_fit(quality = "sampled")
  draws <- .selective_draws("s0")
  draws$b_quality_all <- NULL
  testthat::with_mocked_bindings(
    .cmdstan_extract = function(fit, variables = NULL, ...) draws,
    .package = "antitrustBayes", {
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
    .cmdstan_extract = function(fit, variables = NULL, ...) {
      requested <<- variables
      draws
    }, .package = "antitrustBayes", {
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
    sigma_b_strat_economic = c(.5, .5), sigma_share = c(.2, .2),
    rho = c(.2, .2), sigma_margin = c(.3, .3),
    pred_marginNormInv = matrix(rep(c(.4, .5, .6), 2L), nrow = 2L, byrow = TRUE)
  )
  n_extract <- 0L
  requested <- NULL
  testthat::with_mocked_bindings(
    .cmdstan_extract = function(fit, variables = NULL, ...) {
      n_extract <<- n_extract + 1L
      requested <<- variables
      draws
    }, .package = "antitrustBayes", {
      got <- posterior_quality(fit, seed = 9L)
    }
  )
  expect_equal(n_extract, 1L)
  expect_equal(requested, c("share_mean", "sigma_b_strat_economic", "sigma_share",
                            "rho", "sigma_margin", "pred_marginNormInv"))
  expect_equal(nrow(got$draws), 2L)
  bad <- draws
  bad$pred_marginNormInv <- NULL
  testthat::with_mocked_bindings(
    .cmdstan_extract = function(fit, variables = NULL, ...) bad,
    .package = "antitrustBayes", {
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
    .cmdstan_extract = function(fit, variables = NULL, ...) raw,
    .package = "antitrustBayes", {
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
    .cmdstan_extract = function(fit, variables = NULL, ...) raw_native,
    .package = "antitrustBayes", {
      got_native <- antitrustBayes:::.posterior_extract(fit, "simulation")
    }
  )
  expect_identical(got_native$s0_prob, raw_matrix_prob)

  fit_logit <- .selective_fit_with_schema(c("a_market", "s0_logit"),
                                           model_pars = c("a_market", "s0_prob"))
  requested <- NULL
  testthat::with_mocked_bindings(
    .cmdstan_extract = function(fit, variables = NULL, ...) {
      requested <<- variables
      list(a_market = matrix(c(.8, .9), nrow = 2L), s0_logit = raw_logit)
    },
    .package = "antitrustBayes", {
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
    .cmdstan_extract = function(fit, variables = NULL, ...) bad,
    .package = "antitrustBayes", {
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
