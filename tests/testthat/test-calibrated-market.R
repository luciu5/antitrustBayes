test_that("calibrated Logit markets share one seeded geometry across conduct", {
  skip_if_not_installed("coordination")
  skip_if(utils::packageVersion("antitrust") < package_version("0.99.38"))

  had_seed <- exists(".Random.seed", envir = .GlobalEnv,
                     inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv)
  } else NULL
  set.seed(97)
  caller_seed <- get(".Random.seed", envir = .GlobalEnv)
  labels <- c("MC", "B", "C", "A", "BF", "CF", "BS", "CS")
  generated <- lapply(labels, function(conduct) {
    tryCatch(
      list(ok = TRUE, value = bayes_calibrated_market(
        conduct = conduct, J = 5L, dirichlet_concentration = 1,
        anchor_margin = 0.12, seed = 9101L
      )),
      error = function(e) list(ok = FALSE, error = conditionMessage(e))
    )
  })
  expect_identical(get(".Random.seed", envir = .GlobalEnv), caller_seed)
  if (had_seed) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  failures <- setNames(lapply(generated[!vapply(generated, `[[`, logical(1),
                                                "ok")], `[[`, "error"),
                     labels[!vapply(generated, `[[`, logical(1), "ok")])
  if (length(failures)) {
    fail(paste("conduct generation failures:",
               paste(names(failures), unlist(failures), sep = ": ",
                     collapse = "; ")))
  }
  generated <- lapply(generated, `[[`, "value")

  reference <- generated[[1L]]
  conduct_prices <- lapply(generated, function(x) unname(x$truth$data$prices))
  expect_true(any(vapply(conduct_prices[-1L], function(x)
    max(abs(x - conduct_prices[[1L]])) > 1e-5, logical(1))))
  for (market in generated) {
    expect_true(all(is.finite(market$truth$data$prices) &
                    market$truth$data$prices > 0))
    expect_equal(unname(market$truth$costs),
                 unname(reference$truth$costs), tolerance = 0)
    expect_equal(unname(market$truth$data$prices),
                 unname(market$truth$costs +
                          market$truth$unit_markup / market$truth$alpha),
                 tolerance = 1e-12)
    if (identical(market$design$conduct, "A")) {
      expect_equal(unname(market$truth$auction_value),
                   unname(market$truth$meanval +
                            market$truth$alpha * market$truth$costs),
                   tolerance = 1e-12)
    } else {
      expect_null(market$truth$auction_value)
    }
    expect_equal(unname(antitrust::calcMC(market$truth$model,
                                         preMerger = TRUE)),
                 unname(market$truth$costs), tolerance = 2e-6)
    expect_identical(unname(market$truth$data$shares),
                     unname(reference$truth$data$shares))
    expect_identical(market$design$anchor_margin,
                     reference$design$anchor_margin)
    expect_identical(market$design$anchor_margin, 0.12)
    expect_equal(market$design$realized_anchor_margin,
                 market$truth$margins[[market$design$anchor_index]])
    expect_equal(unname(market$truth$data$shares),
                 unname((1 - market$truth$s0) *
                          market$truth$conditional_shares),
                 tolerance = 0)
    expect_equal(unname(market$truth$data$quantities),
                 unname(10000 * market$truth$conditional_shares),
                 tolerance = 0)
    expected_rank <- market$truth$data$product[
      order(-market$truth$conditional_shares, seq_len(5L))
    ]
    expect_identical(unname(market$design$role_ranking), expected_rank)
    expect_equal(market$truth$s0,
                 stats::plogis(mean(log(market$truth$conditional_shares))),
                 tolerance = 1e-14)
    expect_equal(market$truth$s0_identity,
                 1 - sum(market$truth$shares), tolerance = 0)
    expect_equal(market$truth$s0, 1 - sum(market$truth$shares),
                 tolerance = 1e-14)
    expect_identical(market$design$anchor_index,
                     unname(which.max(market$truth$conditional_shares)))
    expect_equal(market$truth$margins[[market$design$anchor_index]],
                 0.12, tolerance = 1e-7)
    expect_identical(market$estimation$margin_observed,
                     seq_len(5L) == market$design$anchor_index)
    expect_equal(market$estimation$data$margins[[market$design$anchor_index]],
                 market$truth$margins[[market$design$anchor_index]])
    expect_true(all(is.na(market$estimation$data$margins[
      !market$estimation$margin_observed
    ])))
    expect_true(all(is.finite(market$truth$data$margins)))
    expect_true(all(market$truth$costs > 0))
    expect_true(market$truth$equilibrium$gates$perturbed_price_solve)
    expect_true(market$truth$equilibrium$gates$primitive_costs)
    expect_true(market$truth$equilibrium$gates$solved_share_reproduction)
    expect_true(market$truth$equilibrium$used_perturbed_start)
    expect_true(market$truth$equilibrium$passed)
    if (market$design$conduct %in% c("BF", "CF", "BS", "CS")) {
      expect_true(market$truth$equilibrium$gates$baseline_foc)
    } else {
      expect_true(is.na(market$truth$equilibrium$gates$baseline_foc))
    }
    expect_lt(market$fresh_oracle$max_relative_error, 2e-6)
    expect_true(market$fresh_oracle$used_perturbed_start)
    expect_true(any(abs(market$fresh_oracle$start -
                        market$truth$data$prices) > 0.03))
    expect_equal(market$fresh_oracle$shares,
                 market$truth$shares, tolerance = 2e-6)
  }
})

test_that("parameter-driven mode retains supplied alpha and one observed margin", {
  skip_if_not_installed("coordination")
  skip_if(utils::packageVersion("antitrust") < package_version("0.99.38"))

  market <- bayes_calibrated_market(
    conduct = "B", J = 5L, dirichlet_concentration = 1,
    mode = "parameter", alpha = 10,
    seed = 9102L
  )
  expect_equal(market$truth$alpha, 10)
  expect_identical(market$design$mode, "parameter")
  expect_identical(market$estimation$margin_observed,
                   seq_len(5L) == market$design$anchor_index)
  expect_equal(market$estimation$data$margins[[market$design$anchor_index]],
               market$truth$margins[[market$design$anchor_index]])
  expect_true(market$truth$equilibrium$gates$perturbed_price_solve)
  expect_true(is.na(market$truth$equilibrium$gates$baseline_foc))
  expect_true(market$truth$equilibrium$passed)
})

test_that("explicit costs and alpha are common primitives across conduct", {
  primitive_costs <- c(0.8, 1.1, 0.9, 1.2, 1.0)
  conducts <- c("MC", "B", "C", "A", "BF", "CF", "BS", "CS")
  markets <- lapply(conducts, function(conduct)
    suppressWarnings(bayes_calibrated_market(
      conduct = conduct, J = 5L, dirichlet_concentration = 1,
      mode = "parameter", alpha = 12, costs = primitive_costs,
      k = if (conduct %in% c("BF", "CF", "BS", "CS")) 2L else NULL,
      seed = 9103L
    )))
  for (market in markets) {
    expect_equal(unname(market$truth$costs), primitive_costs, tolerance = 0)
    expect_equal(market$truth$alpha, 12, tolerance = 0)
    expect_equal(unname(market$truth$shares),
                 unname(markets[[1L]]$truth$shares), tolerance = 0)
    expect_lt(market$truth$equilibrium$max_cost_error, 2e-6)
    expect_lt(market$truth$equilibrium$max_relative_error, 2e-6)
  }
  expect_true(any(abs(markets[[1L]]$truth$data$prices -
                      markets[[2L]]$truth$data$prices) > 1e-5))
})

test_that("calibrated market validates conduct-specific role counts", {
  expect_error(bayes_calibrated_market(
    conduct = "B", J = 5L, dirichlet_concentration = 1
  ), "margin_range")
  expect_error(bayes_calibrated_market(conduct = "B", J = 5L,
                                       dirichlet_concentration = 1,
                                       anchor_margin = 0.1, k = 2L),
               "only used for BF")
  expect_error(bayes_calibrated_market(conduct = "BS", J = 3L,
                                       dirichlet_concentration = 1,
                                       anchor_margin = 0.1, k = 3L),
               "1:\\(J - 1\\)")
  expect_error(bayes_calibrated_market(mode = "parameter", alpha = 1),
               "dirichlet_concentration.*required")
  expect_error(bayes_calibrated_market(
    dirichlet_concentration = -1
  ), "positive and have length one or J")
  expect_error(bayes_calibrated_market(
    J = 5L, dirichlet_concentration = 1, anchor_margin = .1,
    costs = c(1, 1, 0, 1, 1)
  ), "positive finite vector")
  expect_error(bayes_calibrated_market(
    J = 5L, dirichlet_concentration = 1, anchor_margin = .1,
    cost_range = c(1.3, 0.7)
  ), "two increasing positive")
})
