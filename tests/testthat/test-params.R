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
  expect_true(all(m$margin_observed))

  missing_margins <- bayes_market(
    prices = c(1, 2), quantities = c(2, 1), margins = NULL,
    group = c("a", "b"), outside_price = 1
  )
  expect_false(any(missing_margins$margin_observed))
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
  expect_s4_class(fit, "StructuralFit")
  expect_s4_class(fit@model, "Logit")
  expect_true(all(is.finite(fit@model@pricePre)))
  expect_true(all(is.finite(fit@model@pricePost)))
  expect_identical(fit@diagnostics$specification, "one posterior realization")
  expect_false(fit@diagnostics$demand_recalibrated)
  expect_identical(fit@diagnostics$cost_recovery, "model_foc")
})

test_that("as_antitrust_fit rejects non-conduct objects", {
  expect_error(as_antitrust_fit(list()), "BayesConductLogit or BayesConductCES")
})

test_that("no-cutoff posterior realizations use truthful registered conduct classes", {
  market <- pnb_market()
  owner_post <- pnb_owner_post()
  expected <- list(
    logit = c("Logit", "MonComLogit", "LogitCournot", "Auction2ndLogit",
              "BayesConductLogit"),
    ces = c("CES", "MonComCES", "CESCournot", "Auction2ndCES",
            "BargainingCES")
  )
  conducts <- c("bertrand", "moncom", "cournot", "auction2nd", "bargaining")

  for (supply_model in 1:5) {
    logit_params <- bayes_params(alpha = 0.65, s0 = 0.015,
                                 bargaining_power = 0.5,
                                 supply_model = supply_model)
    logit_obj <- build_conduct_object(logit_params, 1L, market,
                                      ownerPre = market$group,
                                      ownerPost = owner_post)
    logit_fit <- suppressWarnings(as_antitrust_fit(logit_obj))
    expect_s4_class(logit_fit, "AntitrustFit")
    expect_s4_class(logit_fit, "StructuralFit")
    expect_identical(logit_fit@spec$conduct, conducts[[supply_model]])
    expect_s4_class(logit_fit@model, expected$logit[[supply_model]])
    logit_margins <- unname(antitrust::calcMargins(logit_obj, TRUE,
                                                   level = TRUE))
    routed_logit_margins <- unname(antitrust::calcMargins(logit_fit@model,
                                                          TRUE, level = TRUE))
    if (supply_model < 5L) {
      expect_equal(routed_logit_margins, logit_margins, tolerance = 1e-8)
    }

    ces_params <- bayes_params(gamma = 3, s0 = 0.015,
                               bargaining_power = 0.5,
                               supply_model = supply_model)
    ces_obj <- build_conduct_object(ces_params, 1L, market,
                                    ownerPre = market$group,
                                    ownerPost = owner_post)
    ces_fit <- suppressWarnings(as_antitrust_fit(ces_obj))
    expect_s4_class(ces_fit, "AntitrustFit")
    expect_s4_class(ces_fit, "StructuralFit")
    expect_identical(ces_fit@spec$conduct, conducts[[supply_model]])
    expect_s4_class(ces_fit@model, expected$ces[[supply_model]])
    expect_equal(
      unname(antitrust::calcMargins(ces_fit@model, TRUE, level = TRUE)),
      unname(antitrust::calcMargins(ces_obj, TRUE, level = TRUE)),
      tolerance = 1e-8
    )
  }
})

test_that("multiproduct Logit bargaining keeps its explicit custom equation", {
  market <- simple_market()
  params <- bayes_params(alpha = 0.65, s0 = 0.015,
                         bargaining_power = 0.5, supply_model = 5L)
  object <- build_conduct_object(params, 1L, market,
                                 ownerPre = market$group,
                                 ownerPost = market$group)
  fit <- suppressWarnings(as_antitrust_fit(object))
  registered <- suppressWarnings(.bayes_core_model(
    object, .bayes_model_spec(params)
  ))

  expect_s4_class(fit@model, "BayesConductLogit")
  expect_identical(fit@diagnostics$custom_economics,
                   "multiproduct_logit_bargaining_equation")
  expect_gt(max(abs(
    unname(antitrust::calcMargins(fit@model, TRUE, level = TRUE)) -
      unname(antitrust::calcMargins(registered, TRUE, level = TRUE))
  )), 1e-6)
})

test_that("custom posterior fits are scenario-neutral baselines", {
  market <- pnb_market()
  params <- bayes_params(
    alpha = 0.65, s0 = 0.015, supply_model = 1L,
    use_cutoff = TRUE, cutoff_share = 0.2
  )
  object <- build_conduct_object(
    params, 1L, market,
    ownerPre = market$group,
    ownerPost = pnb_owner_post(),
    mcDelta = rep(0.02, market$n_prod),
    subset = seq_len(market$n_prod) != market$n_prod
  )
  fit <- suppressWarnings(as_antitrust_fit(object))
  model <- fit@model

  expect_s4_class(model, "BayesConductLogit")
  expect_equal(model@ownerPost, model@ownerPre)
  expect_equal(model@mcDelta, rep(0, market$n_prod))
  expect_true(all(model@subset))
  expect_equal(model@pricePre, object@prices)
  expect_equal(model@pricePost, object@prices)
  expect_equal(model@mcPost, model@mcPre, tolerance = 1e-10)
  expect_equal(
    antitrust::calcMC(model, preMerger = FALSE),
    model@mcPre,
    tolerance = 1e-10
  )

  scenario <- suppressWarnings(antitrust::simulate(
    fit,
    ownerPost = object@ownerPost,
    mcDelta = object@mcDelta,
    subset = object@subset
  ))
  expect_s4_class(scenario, "BayesConductLogit")
  expect_equal(scenario@ownerPost, object@ownerPost)
  expect_equal(scenario@mcDelta, object@mcDelta)
  expect_equal(scenario@subset, object@subset)
})

test_that("MonCom registration preserves ownership-invariant margins and prices", {
  market <- pnb_market()
  params <- bayes_params(alpha = 0.65, s0 = 0.015, supply_model = 2L)
  owner_post <- pnb_owner_post()
  obj <- build_conduct_object(params, 1L, market,
                              ownerPre = market$group,
                              ownerPost = owner_post)
  fit <- suppressWarnings(as_antitrust_fit(obj))

  expect_identical(fit@spec, antitrust::model_spec("logit", "moncom"))
  expect_equal(
    unname(antitrust::calcMargins(fit@model, TRUE, level = TRUE)),
    unname(antitrust::calcMargins(fit@model, FALSE, level = TRUE)),
    tolerance = 1e-8
  )
  expect_equal(unname(fit@model@pricePost), unname(fit@model@pricePre),
               tolerance = 1e-8)
})

test_that("FOCs recover model costs and positive cost shocks raise prices", {
  market <- pnb_market()
  owner <- market$group
  params <- bayes_params(alpha = 0.65, s0 = 0.015, supply_model = 1L)
  base <- build_conduct_object(params, 1L, market,
                               ownerPre = owner, ownerPost = owner)
  base_fit <- suppressWarnings(as_antitrust_fit(base))
  shocked <- suppressWarnings(antitrust::simulate(
    base_fit,
    ownerPost = owner,
    mcDelta = rep(0.02, market$n_prod)
  ))
  model <- base_fit@model
  foc <- model@pricePre - model@mcPre -
    antitrust::calcMargins(model, TRUE, level = TRUE)

  expect_lt(max(abs(foc)), 1e-8)
  expect_equal(shocked@mcPost, model@mcPre * 1.02, tolerance = 1e-10)
  expect_gt(sum(market$quantities *
                (shocked@pricePost - model@pricePre)), 0)
})

test_that("as_antitrust_fit does not mutate the conduct object", {
  market <- pnb_market()
  params <- bayes_params(alpha = 0.65, s0 = 0.015, supply_model = 1L)
  object <- build_conduct_object(params, 1L, market,
                                 ownerPre = market$group,
                                 ownerPost = pnb_owner_post())
  before <- object
  fit <- suppressWarnings(as_antitrust_fit(object))

  expect_equal(object@slopes, before@slopes)
  expect_equal(object@pricePre, before@pricePre)
  expect_equal(object@pricePost, before@pricePost)
  expect_equal(object@mcPre, before@mcPre)
  expect_s4_class(fit, "StructuralFit")
})

test_that("AntitrustFit slots survive an RDS round trip", {
  market <- pnb_market()
  params <- bayes_params(alpha = 0.65, s0 = 0.015, supply_model = 1L)
  object <- build_conduct_object(params, 1L, market,
                                 ownerPre = market$group,
                                 ownerPost = pnb_owner_post())
  fit <- suppressWarnings(as_antitrust_fit(object))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)

  expect_s4_class(restored, "AntitrustFit")
  expect_s4_class(restored, "StructuralFit")
  expect_equal(restored@spec, fit@spec)
  expect_equal(restored@model@pricePost, fit@model@pricePost,
               tolerance = 1e-10)
})
