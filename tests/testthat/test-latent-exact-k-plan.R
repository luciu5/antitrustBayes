make_latent_exact_k_fixture <- function() {
  groups <- c("a", "b", "c")
  market <- bayes_market(
    prices = c(1.1, 1.2, 1.3), quantities = c(40, 30, 20),
    margins = c(.50, .40, .30), group = groups, outside_price = 1,
    market_side = "input", insideSize = 90, labels = paste0("p", 1:3)
  )
  fit <- bayes_specify(
    bayes_params(
      alpha = 6, s0 = .2, supply_model = 3L, structure = "fringe",
      partition_prior = "ordered_k", k_support = 0:3,
      k_prior = c(.1, .2, .3, .4)
    ),
    market
  )
  list(fit = fit, owner = groups)
}

test_that("policy draw plans bind exact-k role semantics and validate reuse", {
  fx <- make_latent_exact_k_fixture()
  named_fixed_k <- .lq_validate_role_plan_args(
    fx$fit, "map", c(unit_1 = 3L, unit_2 = 2L), 0L
  )$fixed_k
  expect_identical(named_fixed_k, c(unit_1 = 3L, unit_2 = 2L))
  plan <- posterior_policy_draws(
    fx$fit, seed = 12L, role_mapping = "map", fixed_k = 3L,
    role_offset = 1L
  )
  draw <- as.character(plan$draw_indices[[1L]])

  expect_s3_class(plan, "bayes_policy_draw_plan")
  expect_identical(plan$role_mapping, "map")
  expect_identical(plan$fixed_k, 3L)
  expect_identical(plan$role_offset, 1L)
  expect_identical(plan$metadata$role_mapping, "map")
  expect_identical(plan$metadata$fixed_k, 3L)
  expect_identical(plan$metadata$role_offset, 1L)
  expect_identical(plan$k_cache[[draw]][["unit_1"]], 3L)
  expect_no_error(.lq_validate_plan(
    fx$fit, plan, role_mapping = "map", fixed_k = 3L, role_offset = 1L
  ))
  expect_error(.lq_validate_plan(
    fx$fit, plan, role_mapping = "map", fixed_k = 2L, role_offset = 1L
  ), "fixed_k.*conflicts")
  expect_error(.lq_validate_plan(
    fx$fit, plan, role_mapping = "mixture", fixed_k = 3L, role_offset = 1L
  ), "role_mapping.*conflicts")
  expect_error(.lq_validate_plan(
    fx$fit, plan, role_mapping = "map", fixed_k = 3L, role_offset = 0L
  ), "role_offset.*conflicts")

  changed <- plan
  changed$k_cache[[draw]][["unit_1"]] <- 2L
  expect_error(.lq_validate_plan(
    fx$fit, changed, role_mapping = "map", fixed_k = 3L, role_offset = 1L
  ), "latent states")
  changed <- plan
  changed$metadata$fixed_k <- 2L
  expect_error(.lq_validate_plan(
    fx$fit, changed, role_mapping = "map", fixed_k = 3L, role_offset = 1L
  ), "role metadata")
  changed <- plan
  changed$metadata$role_offset <- 1.5
  expect_error(.lq_validate_plan(
    fx$fit, changed, role_mapping = "map", fixed_k = 3L, role_offset = 1L
  ), "role metadata")

  legacy <- posterior_policy_draws(fx$fit, seed = 12L)
  expect_identical(legacy$role_mapping, "mixture")
  expect_null(legacy$fixed_k)
  expect_identical(legacy$role_offset, 0L)
})

test_that("an exact-k latent plan simulates only with matching role controls", {
  skip_if_not_installed("coordination", minimum_version = "0.0.0.9004")
  fx <- make_latent_exact_k_fixture()
  plan <- posterior_policy_draws(
    fx$fit, seed = 12L, role_mapping = "map", fixed_k = 3L,
    role_offset = 1L
  )

  result <- suppressWarnings(simulate(
    fx$fit, fx$owner, draw_plan = plan, role_mapping = "map",
    fixed_k = 3L, role_offset = 1L, demand_anchoring = "latent_quality",
    mc.cores = 1L
  ))
  expect_equal(as.integer(result@n_valid), 1L)
  expect_identical(result@metadata$draw_indices, plan$draw_indices)
  expect_identical(result@metadata$k_cache[[1L]][["unit_1"]], 3L)
  expect_identical(
    as.integer(result@metadata$latent_diagnostics$unit_1[[1L]]$k), 3L
  )
  result_real <- suppressWarnings(simulate(
    fx$fit, fx$owner, draw_plan = plan, role_mapping = "map",
    fixed_k = 3L, role_offset = 1L, demand_anchoring = "latent_quality",
    price_domain = "real", mc.cores = 1L
  ))
  expect_identical(result_real@metadata$draw_indices,
                   result@metadata$draw_indices)
  expect_identical(result_real@metadata$k_cache, result@metadata$k_cache)
  expect_identical(result_real@metadata$draw_plan, result@metadata$draw_plan)
  expect_error(simulate(
    fx$fit, fx$owner, draw_plan = plan, role_mapping = "map",
    fixed_k = 2L, role_offset = 1L, demand_anchoring = "latent_quality",
    mc.cores = 1L
  ), "fixed_k.*conflicts")
})
