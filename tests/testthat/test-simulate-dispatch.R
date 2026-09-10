simulation_fixture <- function(alpha = 0.65, n_draws = 1L) {
  data("pnb", package = "antitrustBayes", envir = environment())
  market <- bayes_market(
    prices = pnb$prices, quantities = pnb$quantities, margins = pnb$margins,
    group = pnb$group, outside_price = 1, market_side = "output",
    insideSize = attr(pnb, "insideSize"), labels = pnb$bank
  )
  params <- bayes_params(
    alpha = rep(alpha, n_draws), s0 = 0.015, supply_model = 1L
  )
  owner_post <- pnb$group
  owner_post[owner_post == attr(pnb, "target")] <- attr(pnb, "acquirer")
  list(
    fit = bayes_specify(params, market), params = params,
    market = market, ownerPost = owner_post, pnb = pnb
  )
}

test_that("BayesFit supplies the StructuralFit simulation contract", {
  fx <- simulation_fixture()

  expect_true(methods::isVirtualClass("BayesFit"))
  expect_s4_class(fx$fit, "BayesSpecifiedFit")
  expect_s4_class(fx$fit, "BayesFit")
  expect_s4_class(fx$fit, "StructuralFit")
  expect_identical(fx$fit@diagnostics$route, "bayes_specify")
  expect_s4_class(methods::selectMethod("simulate", "BayesFit"),
                  "MethodDefinition")
})

test_that("simulate handles point and multi-draw specified fits", {
  point <- simulation_fixture()
  point_result <- suppressWarnings(simulate(
    point$fit, antitrust::counterfactual(ownership = point$ownerPost),
    mc.cores = 1L
  ))
  expect_s4_class(point_result, "BayesSim")
  expect_equal(point_result@n_draws, 1L)
  expect_equal(as.integer(point_result@n_valid), 1L)
  expect_s4_class(point_result@sims[[1L]][[1L]], "Logit")

  draws <- simulation_fixture(n_draws = 3L)
  set.seed(10)
  draw_result <- suppressWarnings(simulate(
    draws$fit, antitrust::counterfactual(ownership = draws$ownerPost),
    n_draws = 3L, mc.cores = 1L
  ))
  expect_equal(draw_result@n_draws, 3L)
  expect_equal(dim(draw_result@price_post[[1L]]), c(3L, draws$market$n_prod))
  expect_setequal(draw_result@metadata$draw_indices, 1:3)
  expect_identical(draw_result@metadata$simulation_route,
                   "antitrust::simulate")
})

test_that("legacy ownership and Counterfactual interfaces agree", {
  fx <- simulation_fixture(n_draws = 3L)
  cf <- antitrust::counterfactual(ownership = fx$ownerPost)

  set.seed(42)
  grammar_result <- suppressWarnings(simulate(
    fx$fit, cf, n_draws = 3L, mc.cores = 1L
  ))
  set.seed(42)
  legacy_result <- suppressWarnings(simulate(
    fx$fit, fx$ownerPost, n_draws = 3L, mc.cores = 1L
  ))

  expect_equal(grammar_result@price_pre, legacy_result@price_pre)
  expect_equal(grammar_result@price_post, legacy_result@price_post,
               tolerance = 1e-8)
})

test_that("bayes_sim is a thin wrapper around simulate", {
  fx <- simulation_fixture(n_draws = 4L)

  set.seed(99)
  direct <- suppressWarnings(simulate(
    fx$fit, fx$ownerPost, n_draws = 3L, mc.cores = 1L
  ))
  set.seed(99)
  wrapped <- suppressWarnings(bayes_sim(
    fx$fit, fx$ownerPost, n_draws = 3L, mc.cores = 1L
  ))

  expect_equal(wrapped@price_pre, direct@price_pre)
  expect_equal(wrapped@price_post, direct@price_post)
  expect_equal(wrapped@metadata$draw_indices, direct@metadata$draw_indices)
  expect_identical(wrapped@call[[1L]], as.name("bayes_sim"))
})

test_that("parallel and serial simulation are deterministic after draw selection", {
  fx <- simulation_fixture(alpha = 1, n_draws = 6L)
  cores <- if (identical(.Platform$OS.type, "windows")) 1L else 2L

  set.seed(7)
  serial <- suppressWarnings(simulate(
    fx$fit, fx$ownerPost, n_draws = 6L, mc.cores = 1L
  ))
  set.seed(7)
  parallel_result <- suppressWarnings(simulate(
    fx$fit, fx$ownerPost, n_draws = 6L, mc.cores = cores
  ))

  expect_equal(serial@price_pre, parallel_result@price_pre)
  expect_equal(serial@price_post, parallel_result@price_post)
  expect_equal(serial@metadata$draw_indices,
               parallel_result@metadata$draw_indices)
})

test_that("simulate validates workers and isolates a bad draw", {
  fx <- simulation_fixture()
  expect_error(simulate(fx$fit, fx$ownerPost, mc.cores = 0), "mc.cores")
  expect_error(simulate(fx$fit, fx$ownerPost, mc.cores = 1.5), "mc.cores")

  params <- antitrustBayes:::.bayes_params_raw(
    alpha = rep(1, 6), s0 = c(0.015, 0.015, 1.5, 0.015, 0.015, 0.015),
    bargaining_power = rep(0, 6), cutoff_share = rep(0, 6),
    cutoff_scale = rep(0.2, 6), use_cutoff = FALSE, supply_model = 1L
  )
  bad_fit <- bayes_specify(params, fx$market)
  cores <- if (identical(.Platform$OS.type, "windows")) 1L else 2L

  expect_warning(
    result <- simulate(
      bad_fit, fx$ownerPost, n_draws = 6L, mc.cores = cores
    ),
    "1 of 6"
  )
  expect_equal(as.integer(result@n_valid), 5L)
  expect_equal(sum(!is.na(result@failures[[1L]])), 1L)
})

test_that("antitrust scenario machinery supports non-merger counterfactuals", {
  fx <- simulation_fixture()
  cost <- suppressWarnings(simulate(
    fx$fit,
    antitrust::counterfactual(costs = rep(0.02, fx$market$n_prod)),
    mc.cores = 1L
  ))
  exit <- suppressWarnings(simulate(
    fx$fit,
    antitrust::counterfactual(exit = fx$market$labels[[1L]]),
    mc.cores = 1L
  ))
  quality <- suppressWarnings(simulate(
    fx$fit,
    antitrust::counterfactual(
      quality = stats::setNames(0.1, fx$market$labels[[1L]])
    ),
    mc.cores = 1L
  ))

  expect_equal(as.integer(cost@n_valid), 1L)
  expect_equal(as.integer(exit@n_valid), 1L)
  expect_equal(as.integer(quality@n_valid), 1L)
  exit_model <- antitrustBayes:::.terminal_model(exit@sims[[1L]][[1L]])
  expect_false(exit_model@subset[match(fx$market$labels[[1L]],
                                       exit_model@labels)])
})

test_that("entry changes the assembled product dimension", {
  fx <- simulation_fixture()
  baseline <- build_conduct_object(
    fx$params, 1L, fx$market,
    ownerPre = fx$market$group, ownerPost = fx$market$group
  )
  baseline_fit <- suppressWarnings(as_antitrust_fit(baseline))
  entrant <- antitrust::entrant(
    "E1",
    meanval = mean(baseline_fit@model@slopes$meanval),
    cost = max(0.01, median(baseline_fit@model@mcPre)),
    priceStart = median(fx$market$prices)
  )

  result <- suppressWarnings(simulate(
    fx$fit, antitrust::counterfactual(entry = entrant), mc.cores = 1L
  ))
  expect_equal(ncol(result@price_post[[1L]]), fx$market$n_prod + 1L)
  expect_true("E1" %in% result@product_labels[[1L]])
})

test_that("multi-step paths remain stored while downstream methods terminalize", {
  fx <- simulation_fixture(n_draws = 2L)
  cf <- antitrust::counterfactual(ownership = fx$ownerPost)
  cf <- antitrust::add_step(cf, costs = rep(0.01, fx$market$n_prod))
  result <- suppressWarnings(simulate(
    fx$fit, cf, n_draws = 2L, mc.cores = 1L
  ))

  expect_true(all(vapply(result@sims[[1L]], methods::is,
                         logical(1), "CounterfactualPath")))
  elasticity <- elast(result, posterior = TRUE)
  expect_equal(dim(elasticity)[[3L]], 2L)
  expect_silent(invisible(capture.output(summary(result))))
})

test_that("antitrust's Counterfactual grammar is not copied into this namespace", {
  exports <- getNamespaceExports("antitrustBayes")
  expect_false(any(c(
    "bayes_counterfactual", "counterfactual", "add_step",
    "combine_counterfactuals", "entrant", "final_result", "result_at"
  ) %in% exports))
})

test_that("estimated fits dispatch through the same S4 method", {
  skip_if(tolower(Sys.getenv("ANTITRUSTBAYES_RUN_STAN_TESTS")) != "true",
          "Set ANTITRUSTBAYES_RUN_STAN_TESTS=true for the opt-in Stan test")
  skip_if_not_installed("rstan")

  set.seed(123)
  n <- 12
  prices <- runif(n, 1, 5)
  quantities <- runif(n, 10, 100)
  margins <- runif(n, 0.1, 0.4)
  group <- rep(letters[1:4], each = 3)
  fit <- suppressWarnings(bayes_est(
    prices = prices, quantities = quantities, margins = margins,
    market = rep("m1", n), group = group,
    demand_model = "logit", supply_model = "bertrand",
    stan_options = list(
      chains = 1, iter = 60, warmup = 30, refresh = 0, seed = 1
    )
  ))
  owner_post <- group
  owner_post[owner_post == "b"] <- "a"

  set.seed(4)
  direct <- simulate(fit, owner_post, n_draws = 2L, mc.cores = 1L)
  set.seed(4)
  wrapped <- bayes_sim(fit, owner_post, n_draws = 2L, mc.cores = 1L)

  expect_s4_class(fit, "BayesFit")
  expect_s4_class(fit, "StructuralFit")
  expect_equal(direct@price_post, wrapped@price_post, tolerance = 1e-8)
})
