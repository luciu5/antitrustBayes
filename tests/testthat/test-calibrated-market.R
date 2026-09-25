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
  for (market in generated) {
    expect_identical(unname(market$truth$data$prices), rep(1, 5L))
    expect_identical(unname(market$truth$data$prices),
                     unname(reference$truth$data$prices))
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
    expect_true(market$truth$equilibrium$gates$solved_share_reproduction)
    expect_true(market$truth$equilibrium$used_perturbed_start)
    if (market$design$conduct %in% c("BF", "CF", "BS", "CS")) {
      expect_true(market$truth$equilibrium$gates$baseline_foc)
      expect_true(market$truth$equilibrium$passed)
    } else {
      expect_true(is.na(market$truth$equilibrium$gates$baseline_foc))
      expect_true(is.na(market$truth$equilibrium$passed))
    }
    expect_lt(market$fresh_oracle$max_relative_error, 2e-6)
    expect_true(market$fresh_oracle$used_perturbed_start)
    expect_false(all(market$fresh_oracle$start == 1))
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
  expect_true(is.na(market$truth$equilibrium$passed))
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
})
