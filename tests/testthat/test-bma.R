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

test_that("BMA preserves weights and exposes joint source attrition", {
  a <- make_test_sim("a")
  a@sims$m1[2] <- list(NULL)
  a@sims$m2[1] <- list(NULL)
  a@metadata$draw_indices <- c(10L, 20L, 30L)
  expect_warning(result <- bayes_avg(list(a = a, b = make_test_sim("b")),
    weights = c(.5, .5), n_draws = 10L), "source attrition")
  expect_equal(unname(result@metadata$weights), c(.5, .5))
  expect_equal(result@metadata$source_attrition$a$common_eligible, 1L)
  expect_equal(unname(result@metadata$source_attrition$a$failed), c(1L, 1L))
  selected <- result@metadata$source_selection
  expect_true(all(selected$draw_id[selected$model_index == 1] == 30))
  expect_true(all(unlist(result@sims$m2)[selected$model_index == 1] == "a-m2-3"))
  a@sims$m1 <- rep(list(NULL), 3L)
  expect_error(bayes_avg(list(a, make_test_sim("b")), c(.0001, .9999), 10L),
    "no draw that is valid")
})

test_that("simulation summaries report failures without shifting draw IDs", {
  sim <- make_test_sim("a", markets = "m1", n = 4L)
  sim@sims$m1[1] <- list(NULL)
  sim@metadata$draw_indices <- c(10L, 20L, 30L, 40L)
  testthat::with_mocked_bindings(
    .bayes_sim_draw_summary = function(model, ...) {
      if (model == "a-m1-2") stop("deliberate summary failure")
      data.frame(price = if (model == "a-m1-3") NA_real_ else 2)
    }, {
      expect_warning(capture.output(result <- summary(sim)), "omitted draws or values")
      expect_equal(result$mean$price, 2)
      d <- attr(result, "diagnostics")$m1
      expect_equal(d$summary_draw_ids, c(30L, 40L))
      expect_equal(d$errors[2], "deliberate summary failure")
      expect_equal(d$estimands$valid, 1)
      expect_equal(d$estimands$missing, 3)
      expect_equal(d$nonfinite_values$draw_id, 30L)
      sim@sims$m1 <- rep(list("a-m1-2"), 4L)
      warning_data <- NULL
      all_failed <- withCallingHandlers(summary(sim), warning = function(w) {
        warning_data <<- w
        invokeRestart("muffleWarning")
      })
      expect_null(all_failed)
      expect_s3_class(warning_data, "bayes_summary_attrition")
      expect_match(conditionMessage(warning_data), "4 requested")
      expect_true(all(warning_data$diagnostics$m1$errors == "deliberate summary failure"))
    })
})

test_that("bridge weights never replace missing evidence with predictive stacking", {
  comp <- data.frame(
    model = c("a", "b"),
    weight_stacking = c(0.25, 0.75),
    post_prob = c(NA_real_, NA_real_)
  )
  expect_error(bayes_weights(comp, "bridge"), "Every positive-prior")
  expect_equal(bayes_weights(comp, "stacking"), c(a = 0.25, b = 0.75))
  comp$post_prob <- c(.3, NA_real_)
  expect_error(bayes_weights(comp, "bridge"), "no candidates can be dropped")
  comp$post_prob <- c(.3, .7)
  expect_equal(bayes_weights(comp, "bridge"), c(a = .3, b = .7))
  expect_error(bayes_weights(comp[1L, ], "bridge"), "complete model set")
  comp$model_prior <- c(1, 0); comp$post_prob <- c(1, NA_real_)
  expect_equal(bayes_weights(comp, "bridge"), c(a = 1, b = 0))
  comp$post_prob <- c(.8, .2)
  expect_error(bayes_weights(comp, "bridge"), "zero-prior")
  comp$post_prob <- NULL
  expect_error(bayes_weights(comp, "bridge"), "request method = 'stacking' explicitly")
})

make_contract_fit <- function(demand_model = "logit") {
  stan_data <- list(
    N = 2L,
    shareIn = c(0.4, 0.6),
    clr_shareIn = log(c(0.4, 0.6)) - mean(log(c(0.4, 0.6))),
    marginNormInv = c(4, 5),
    log_marginNormInv = log(c(4, 5)),
    margin_obs = c(1L, 1L),
    market_year_idx = c(1L, 1L),
    cluster = c(1L, 1L),
    likelihood_mode = 1L
  )
  attr(stan_data, "model_spec_version") <- 13L
  attr(stan_data, "quality_scope") <- "global"
  new("BayesLogit",
    stanfit = NULL,
    data = data.frame(
      prices = c(1, 2), quantities = c(4, 6), margins = c(0.25, 0.2),
      market = c("m", "m"), year = c("1", "1"), group = c("a", "b"),
      loo_cluster = c("case", "case"), margin_observed = c(TRUE, TRUE),
      outside_price = c(1, 1)
    ),
    stan_data = stan_data,
    supply_model = 1L,
    demand_model = demand_model,
    market_side = "output",
    price_sd = 1,
    call = call("bayes_est"),
    metadata = list(
      prior_manifest = "baseline",
      resolved_priors = bayes_prior_manifest("baseline"),
      model_prior = NA_real_
    )
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

test_that("bayes_compare dispatches exact-k restrictions through one fit", {
  stan_data <- list(
    N_partition_unit = 1L,
    partition_unit_market = 1L,
    partition_unit_start = 1L,
    partition_unit_count = 2L,
    partition_component_unit = c(1L, 1L),
    partition_component_k = c(1L, 2L)
  )
  attr(stan_data, "market_levels") <- "m1"
  fit <- new(
    "BayesLogit", stanfit = NULL, data = data.frame(market = "m1"),
    stan_data = stan_data, price_sd = 1, supply_model = 1L,
    demand_model = "logit", market_side = "output", call = call("fake"),
    metadata = list(structure = "stackelberg", partition_prior = "ordered_k",
                    posterior_draws = 2L)
  )
  draws <- list(
    posterior_component_prob = rbind(c(.8, .2), c(.6, .4)),
    prior_component_prob = matrix(.5, nrow = 2L, ncol = 2L)
  )
  got <- bayes_compare(
    fit, restrictions = list(k1 = 1L, k2 = 2L), reference = "k2",
    restriction_draws = draws
  )
  expect_s3_class(got, "bayes_restriction_comparison")
  expect_identical(got$method, rep("encompassing_prior", 2L))
  expect_equal(got$posterior_mass, c(.7, .3))
  expect_equal(got$bf_vs_reference, c(7 / 3, 1))
  expect_error(
    bayes_compare(fit, fit, restrictions = list(k1 = 1L),
                  restriction_draws = draws),
    "exactly one BayesFit"
  )
  expect_error(bayes_compare(fit, reference = "k1"), "require 'restrictions'")
  expect_error(
    bayes_compare(fit, restrictions = list(k1 = 1L),
                  restriction_draws = draws, prior_weights = 1),
    "fit's recorded k prior"
  )
})

test_that("evidence decomposition accepts common or model-specific share evidence", {
  comparison <- data.frame(model = c("a", "b"), logml = c(4, 6))
  common <- bayes_evidence_decomposition(comparison, 2)
  expect_equal(common$logml_margin_given_share, c(2, 4))
  specific <- bayes_evidence_decomposition(comparison, c(1, 4))
  expect_equal(specific$logml_margin_given_share, c(3, 2))
  expect_error(
    bayes_evidence_decomposition(comparison, c(1, 2, 3)),
    "one finite value or one value per model"
  )
})
