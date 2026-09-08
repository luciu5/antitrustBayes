## The PNB dataset's prices/quantities/margins are jointly consistent (see
## bayes_est()'s own examples), so a hand-specified alpha/s0 in a plausible
## range solves without implying negative marginal costs -- unlike an
## arbitrary small synthetic market, which readily produces an inconsistent
## demand inversion for a hand-picked alpha.
pnb_fixture <- function() {
  data("pnb", package = "antitrustBayes", envir = environment())
  market <- bayes_market(
    prices = pnb$prices, quantities = pnb$quantities, margins = pnb$margins,
    group = pnb$group, outside_price = 1, market_side = "output",
    insideSize = attr(pnb, "insideSize"), labels = pnb$group
  )
  owner_post <- pnb$group
  owner_post[owner_post == attr(pnb, "target")] <- attr(pnb, "acquirer")
  list(market = market, ownerPost = owner_post, pnb = pnb)
}

test_that("bayes_counterfactual runs a merger from a point estimate with no Stan fit", {
  fx <- pnb_fixture()
  params <- bayes_params(alpha = 0.65, s0 = 0.015, supply_model = 1L)

  result <- suppressWarnings(
    bayes_counterfactual(params, fx$market, ownerPost = fx$ownerPost)
  )

  expect_s4_class(result, "BayesConductLogit")
  merging <- fx$pnb$group %in% c(attr(fx$pnb, "acquirer"), attr(fx$pnb, "target"))
  expect_true(all(result@pricePost[merging] > result@pricePre[merging] - 1e-8))
})

test_that("bayes_counterfactual accepts a draws matrix and returns a BayesSim", {
  fx <- pnb_fixture()
  params <- bayes_params(alpha = c(0.60, 0.65, 0.70), s0 = 0.015, supply_model = 1L)

  result <- suppressWarnings(
    bayes_counterfactual(params, fx$market, ownerPost = fx$ownerPost, n_draws = 3L)
  )

  expect_s4_class(result, "BayesSim")
  expect_equal(result@n_draws, 3L)
  expect_equal(nrow(result@price_pre[[1]]), 3L)
})

test_that(".default_mc_cores caps at 4 and returns 1 on Windows", {
  cores <- antitrustBayes:::.default_mc_cores()
  expect_true(cores >= 1L)
  expect_true(cores <= 4L)
  if (identical(.Platform$OS.type, "windows")) {
    expect_equal(cores, 1L)
  }
})

test_that("bayes_counterfactual validates mc.cores", {
  fx <- pnb_fixture()
  params <- bayes_params(alpha = 0.65, s0 = 0.015, supply_model = 1L)
  expect_error(
    bayes_counterfactual(params, fx$market, ownerPost = fx$ownerPost, mc.cores = 0),
    "mc.cores"
  )
  expect_error(
    bayes_counterfactual(params, fx$market, ownerPost = fx$ownerPost, mc.cores = 1.5),
    "mc.cores"
  )
})

test_that("bayes_counterfactual gives identical results with mc.cores = 1 and > 1", {
  fx <- pnb_fixture()
  params <- bayes_params(alpha = seq(0.9, 1.1, length.out = 6), s0 = 0.015,
                          supply_model = 1L)

  set.seed(42)
  serial <- suppressWarnings(bayes_counterfactual(
    params, fx$market, ownerPost = fx$ownerPost, n_draws = 6L, mc.cores = 1L
  ))
  set.seed(42)
  parallel_cores <- if (identical(.Platform$OS.type, "windows")) 1L else 2L
  parallel_result <- suppressWarnings(bayes_counterfactual(
    params, fx$market, ownerPost = fx$ownerPost, n_draws = 6L,
    mc.cores = parallel_cores
  ))

  expect_equal(serial@price_pre, parallel_result@price_pre)
  expect_equal(serial@price_post, parallel_result@price_post)
})

test_that("bayes_counterfactual isolates a single bad draw under mc.cores > 1", {
  fx <- pnb_fixture()
  ## Force one out-of-range s0 via the internal unvalidated constructor,
  ## bypassing bayes_params()'s upfront validation, to simulate a draw that
  ## fails deep inside the solver rather than at input validation.
  params <- antitrustBayes:::.bayes_params_raw(
    alpha = rep(1, 6), s0 = c(0.015, 0.015, 1.5, 0.015, 0.015, 0.015),
    bargaining_power = rep(0, 6), cutoff_share = rep(0, 6),
    cutoff_scale = rep(0.2, 6), use_cutoff = FALSE, supply_model = 1L
  )
  parallel_cores <- if (identical(.Platform$OS.type, "windows")) 1L else 2L

  warnings_seen <- character(0)
  result <- withCallingHandlers(
    bayes_counterfactual(params, fx$market, ownerPost = fx$ownerPost,
                        n_draws = 6L, mc.cores = parallel_cores),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("1 of 6", warnings_seen)))
  expect_equal(as.integer(result@n_valid), 5L)
  expect_equal(sum(!is.na(result@failures[[1]])), 1L)
})

test_that("bayes_counterfactual accepts an antitrust Counterfactual for cost shocks", {
  fx <- pnb_fixture()
  params <- bayes_params(alpha = 0.65, s0 = 0.015, supply_model = 1L)
  n_prod <- fx$market$n_prod

  cf <- counterfactual(ownership = fx$ownerPost, costs = rep(0, n_prod))
  result <- suppressWarnings(bayes_counterfactual(params, fx$market, ownerPost = cf))

  final <- if (is(result, "CounterfactualPath")) final_result(result) else result
  expect_true(all(is.finite(final@pricePost)))
})

test_that("bayes_counterfactual supports a multi-step path via add_step", {
  fx <- pnb_fixture()
  params <- bayes_params(alpha = 0.65, s0 = 0.015, supply_model = 1L)
  n_prod <- fx$market$n_prod

  cf <- counterfactual(ownership = fx$ownerPost)
  cf <- add_step(cf, costs = rep(0.02, n_prod))
  result <- suppressWarnings(bayes_counterfactual(params, fx$market, ownerPost = cf))

  expect_s4_class(result, "CounterfactualPath")
  final <- final_result(result)
  expect_true(all(is.finite(final@pricePost)))
})

test_that("the manual-parameter path is deterministic given the same inputs", {
  fx <- pnb_fixture()
  params <- bayes_params(alpha = 0.65, s0 = 0.015, supply_model = 1L)

  fit1 <- suppressWarnings(as_antitrust_fit(build_conduct_object(
    params, 1L, fx$market, ownerPre = fx$market$group, ownerPost = fx$ownerPost
  )))
  fit2 <- suppressWarnings(as_antitrust_fit(build_conduct_object(
    params, 1L, fx$market, ownerPre = fx$market$group, ownerPost = fx$ownerPost
  )))

  expect_equal(fit1@model@pricePre, fit2@model@pricePre)
  expect_equal(fit1@model@pricePost, fit2@model@pricePost, tolerance = 1e-8)
})

test_that("the manual-parameter path exactly reproduces one of bayes_sim()'s own draws", {
  ## Critical equivalence check: bayes_sim() was refactored to route through
  ## build_conduct_object()/as_antitrust_fit() (see R/simulate.R). Extract one
  ## posterior draw's raw parameters directly from the stanfit, rebuild the
  ## conduct object by hand through the public bayes_params()/bayes_market()
  ## API, and confirm the manually-produced post-merger prices appear
  ## verbatim among bayes_sim()'s own results for that fit -- proving the
  ## refactor did not change bayes_sim()'s behavior.
  skip_if_not_installed("rstan")

  set.seed(123)
  N <- 12
  prices <- runif(N, 1, 5)
  quantities <- runif(N, 10, 100)
  margins <- runif(N, 0.1, 0.4)
  group <- rep(letters[1:4], each = 3)

  fit <- suppressWarnings(bayes_est(
    prices = prices, quantities = quantities, margins = margins,
    market = rep("m1", N), group = group,
    demand_model = "logit", supply_model = "bertrand",
    stan_options = list(chains = 1, iter = 60, warmup = 30, refresh = 0, seed = 1)
  ))

  owner_post <- group
  owner_post[owner_post == "b"] <- "a"

  n_total <- length(rstan::extract(fit@stanfit)$a_market)
  sim_all <- bayes_sim(fit, ownerPost = owner_post, n_draws = n_total)

  draws <- rstan::extract(fit@stanfit)
  i <- 3L
  alpha_draw <- draws$a_market[i] / fit@price_sd
  s0_draw <- draws$s0_prob[i, 1]

  params <- bayes_params(alpha = alpha_draw, s0 = s0_draw, supply_model = 1L)
  market <- bayes_market(
    prices = prices, quantities = quantities, margins = margins,
    group = group, outside_price = unique(fit@data$outside_price),
    market_side = fit@market_side, labels = group
  )
  manual_obj <- build_conduct_object(params, 1L, market, ownerPre = group,
                                     ownerPost = owner_post)
  manual_fit <- as_antitrust_fit(manual_obj)
  manual_post <- as.numeric(manual_fit@model@pricePost)

  matches <- apply(sim_all@price_post[[1]], 1, function(row) {
    isTRUE(all.equal(as.numeric(row), manual_post, tolerance = 1e-8))
  })
  expect_true(any(matches))
})
