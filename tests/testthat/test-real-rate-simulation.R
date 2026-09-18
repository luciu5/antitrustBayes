.rate_domain_fixture <- function(structure = "fringe", demand = "logit",
                                 side = "input", supply = 3L) {
  args <- list(s0 = .1, supply_model = supply, structure = structure,
               partition_prior = "ordered_k", k_support = 1:2,
               k_prior = c(1, 0))
  args[[if (demand == "logit") "alpha" else "gamma"]] <-
    if (demand == "logit") .8 else 10
  params <- do.call(bayes_params, args)
  market <- bayes_market(
    prices = c(.2, .1, .3), quantities = c(45, 35, 20),
    margins = rep(.4, 3), group = c("a", "b", "c"),
    outside_price = 1, market_side = side, insideSize = 100,
    labels = c("a", "b", "c"))
  bayes_specify(params, market)
}

test_that("real rate domains are explicit and reject unsupported economics", {
  fit <- .rate_domain_fixture()
  expect_identical(.validate_simulation_price_domain(fit, c("positive", "real")),
                   "positive")
  expect_identical(.validate_simulation_price_domain(fit, "real"), "real")
  for (other in list(
    .rate_domain_fixture(side = "output"),
    .rate_domain_fixture(demand = "ces", side = "output"),
    .rate_domain_fixture(supply = 1L))) {
    expect_error(.validate_simulation_price_domain(other, "real"),
                 "only input Logit Cournot")
  }
  fit@metadata$structure <- "simultaneous"
  expect_error(.validate_simulation_price_domain(fit, "real"),
               "only input Logit Cournot")
})

test_that("failure metadata retains domain obstruction and solver stage", {
  error <- structure(list(message = "positive domain excludes root",
                          category = "positive_domain_obstruction",
                          domain = "positive", minimum_rate = -.2,
                          minimum_product = "b", preMerger = FALSE),
                     class = c("rate_domain_error", "error", "condition"))
  info <- .bayes_sim_failure_detail(error, "post_equilibrium", 17L)
  expect_identical(info$category, "positive_domain_obstruction")
  expect_identical(info$stage, "post_equilibrium")
  expect_identical(info$draw_index, 17L)
  expect_identical(info$price_domain, "positive")
  expect_identical(info$minimum_product, "b")
  expect_equal(info$minimum_rate, -.2)
})

test_that("assembly preserves failed draw positions and undefined ratios", {
  raw <- list(list(result = 1L, error = NA_character_),
              list(result = 2L, error = NA_character_),
              list(result = 3L, error = NA_character_))
  testthat::with_mocked_bindings(
    .result_values = function(result, fallback_labels) {
      if (result == 2L) stop("deliberate assembly failure")
      list(price_pre = c(0, -1), price_post = c(.5, -.5),
           shares_pre = c(.2, .3), shares_post = c(.3, .4),
           labels = c("a", "b"), latent_quality = NULL)
    },
    {
      got <- .assemble_unit(raw, list(labels = c("a", "b")), c(10L, 20L, 30L))
      expect_length(got$results, 3L)
      expect_null(got$results[[2L]])
      expect_identical(got$results[[3L]], 3L)
      expect_equal(got$n_valid, 2L)
      expect_true(all(is.na(got$price_post[2L, ])))
      expect_equal(got$price_post[3L, ], c(a = .5, b = -.5))
      expect_true(is.na(got$price_change[1L, 1L]))
      expect_true(got$undefined_relative_change[1L, 1L])
      expect_true(got$signed_relative_change[1L, 2L])
      expect_equal(unname(got$price_change[1L, 2L]), -.5)
      expect_identical(got$failure_details[[2L]]$stage, "result_assembly")
      expect_identical(got$failure_details[[2L]]$draw_index, 20L)
    })
})

test_that("native real rates reuse policy plans and preserve value primitives", {
  skip_if_not_installed("coordination", minimum_version = "0.0.0.9004")
  for (structure in c("fringe", "stackelberg")) {
    fit <- .rate_domain_fixture(structure)
    plan <- posterior_policy_draws(fit, seed = 34L)
    ownership <- c("a", "a", "c")
    merger <- suppressWarnings(simulate(
      fit, ownership, draw_plan = plan, mc.cores = 1L,
      price_domain = "real"))
    both <- suppressWarnings(simulate(
      fit, antitrust::counterfactual(ownership = ownership, costs = rep(.1, 3)),
      draw_plan = plan, mc.cores = 1L, price_domain = "real",
      demand_anchoring = "latent_quality"))
    expect_equal(as.integer(merger@n_valid), 1L)
    expect_equal(as.integer(both@n_valid), 1L)
    expect_identical(merger@metadata$price_domain, "real")
    expect_identical(merger@metadata$policy_distribution_conditioning,
                     "all_requested_draws")
    expect_identical(merger@metadata$k_cache, both@metadata$k_cache)
    expect_identical(merger@metadata$draw_indices, both@metadata$draw_indices)
    model <- both@sims$unit_1[[1L]]
    base <- merger@sims$unit_1[[1L]]
    expect_equal(model@mcPre, base@mcPre, tolerance = 1e-10)
    expect_equal(model@mcPost, 1.1 * model@mcPre, tolerance = 1e-10)
    expect_lt(max(abs(model@mcPost - model@pricePost -
                       antitrust::calcMargins(model, preMerger = FALSE,
                                              level = TRUE))), 1e-6)
    expect_equal(both@metadata$price_level_change$unit_1,
                 both@price_post$unit_1 - both@price_pre$unit_1)
    diag <- both@metadata$latent_diagnostics$unit_1[[1L]]
    expect_equal(diag$quality, unname(plan$quality_draws[1L, ]))
    expect_identical(diag$price_domain, "real")
    if (structure == "fringe") {
      expect_true(any(merger@price_post$unit_1 < 0))
      positive <- suppressWarnings(simulate(
        fit, ownership, draw_plan = plan, mc.cores = 1L,
        price_domain = "positive"))
      expect_equal(as.integer(positive@n_valid), 0L)
      detail <- positive@metadata$failure_details$unit_1[[1L]]
      expect_identical(detail$stage, "post_equilibrium")
      expect_identical(detail$category, "positive_domain_obstruction")
      expect_lt(detail$minimum_rate, 0)
      expect_identical(positive@metadata$policy_distribution_conditioning,
                       "successful_draws_only")
    }
  }
})

test_that("the positive native route remains the default", {
  skip_if_not_installed("coordination")
  fit <- .rate_domain_fixture()
  plan <- posterior_policy_draws(fit, seed = 19L)
  implicit <- simulate(fit, c("a", "b", "c"), draw_plan = plan, mc.cores = 1L)
  explicit <- simulate(fit, c("a", "b", "c"), draw_plan = plan, mc.cores = 1L,
                       price_domain = "positive")
  expect_identical(implicit@price_pre, explicit@price_pre)
  expect_identical(implicit@price_post, explicit@price_post)
  expect_identical(implicit@metadata$price_domain, "positive")
})
