simple_market <- function() {
  bayes_market(
    prices = c(1, 1.2, 0.9),
    quantities = c(20, 15, 10),
    margins = rep(0.2, 3),
    group = c("a", "a", "b"),
    outside_price = 1,
    market_side = "output",
    insideSize = 1
  )
}

## The Philadelphia National Bank dataset shipped with the package (see
## bayes_est()'s own examples). Its prices/quantities/margins are jointly
## consistent, so a hand-specified alpha/s0 solves without implying negative
## marginal costs -- unlike the arbitrary 3-firm fixture above, which is only
## used for share-reproduction checks that never invoke the price solver.
pnb_market <- function() {
  data("pnb", package = "antitrustBayes", envir = environment())
  bayes_market(
    prices = pnb$prices, quantities = pnb$quantities, margins = pnb$margins,
    group = pnb$group, outside_price = 1, market_side = "output",
    insideSize = attr(pnb, "insideSize"), labels = pnb$group
  )
}

pnb_owner_post <- function() {
  data("pnb", package = "antitrustBayes", envir = environment())
  owner_post <- pnb$group
  owner_post[owner_post == attr(pnb, "target")] <- attr(pnb, "acquirer")
  owner_post
}

test_that("bayes_params requires exactly one of alpha/gamma", {
  expect_error(bayes_params(s0 = 0.5, supply_model = 1L),
               "exactly one")
  expect_error(bayes_params(alpha = -2, gamma = 3, s0 = 0.5, supply_model = 1L),
               "exactly one")
})

test_that("bayes_params validates gamma, s0, bargaining_power, supply_model", {
  expect_error(bayes_params(gamma = 0.5, s0 = 0.5, supply_model = 1L),
               "exceed one")
  expect_error(bayes_params(alpha = -2, s0 = 1.5, supply_model = 1L),
               "strictly between 0 and 1")
  expect_error(bayes_params(alpha = -2, s0 = 0.5, supply_model = 5L,
                             bargaining_power = 1.5),
               "\\[0, 1\\)")
  expect_error(bayes_params(alpha = -2, s0 = 0.5, supply_model = 9L),
               "supply_model")
})

test_that("bayes_params accepts a point estimate and a draws vector", {
  point <- bayes_params(alpha = -2, s0 = 0.5, supply_model = 1L)
  expect_equal(point$n_draws, 1L)
  expect_equal(point$demand_model, "logit")

  draws <- bayes_params(alpha = c(-2, -2.5, -1.8), s0 = 0.5, supply_model = 1L)
  expect_equal(draws$n_draws, 3L)
  expect_length(draws$s0, 3L)
})

test_that("bayes_market validates lengths and defaults", {
  expect_error(
    bayes_market(prices = c(1, 2), quantities = c(1, 2, 3), group = c("a", "b"),
                 outside_price = 1),
    "same length"
  )
  m <- simple_market()
  expect_equal(m$n_prod, 3L)
  expect_equal(m$labels, c("a", "a", "b"))
})

test_that("build_conduct_object reproduces the supplied quantity shares", {
  params <- bayes_params(alpha = -2, s0 = 0.5, supply_model = 1L)
  market <- simple_market()
  obj <- suppressWarnings(build_conduct_object(
    params, 1L, market, ownerPre = market$group, ownerPost = market$group
  ))
  expect_s4_class(obj, "BayesConductLogit")
  expect_equal(
    antitrust::calcShares(obj),
    market$quantities / sum(market$quantities) * (1 - 0.5),
    tolerance = 1e-10, ignore_attr = TRUE
  )
})

test_that("build_conduct_object builds CES objects with revenue shares", {
  params <- bayes_params(gamma = 3, s0 = 0.5, supply_model = 1L)
  market <- simple_market()
  obj <- suppressWarnings(build_conduct_object(
    params, 1L, market, ownerPre = market$group, ownerPost = market$group
  ))
  expect_s4_class(obj, "BayesConductCES")
})

test_that("as_antitrust_fit wraps a conduct object with a solved model", {
  params <- bayes_params(alpha = 0.65, s0 = 0.015, supply_model = 1L)
  market <- pnb_market()
  obj <- build_conduct_object(params, 1L, market, ownerPre = market$group,
                              ownerPost = pnb_owner_post())
  fit <- suppressWarnings(as_antitrust_fit(obj))
  expect_s4_class(fit, "AntitrustFit")
  expect_s4_class(fit@model, "BayesConductLogit")
  expect_true(all(is.finite(fit@model@pricePre)))
  expect_true(all(is.finite(fit@model@pricePost)))
})

test_that("as_antitrust_fit rejects non-conduct objects", {
  expect_error(as_antitrust_fit(list()), "BayesConductLogit or BayesConductCES")
})
