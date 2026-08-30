make_conduct_object <- function(demand = c("logit", "ces"), supply_model = 1L,
                                use_cutoff = FALSE, cutoff_share = 0) {
  demand <- match.arg(demand)
  prices <- c(1, 1.2, 0.9)
  quantity_shares <- c(0.20, 0.15, 0.10)
  s0 <- 1 - sum(quantity_shares)
  owner_pre <- c("a", "a", "b")
  owner_post <- c("a", "a", "a")

  if (demand == "logit") {
    alpha <- -2
    slopes <- list(
      alpha = alpha,
      meanval = log(quantity_shares / s0) - alpha * (prices - 1)
    )
    class_name <- "BayesConductLogit"
  } else {
    gamma <- 3
    total_expenditure <- s0 + sum(prices * quantity_shares)
    revenue_shares <- prices * quantity_shares / total_expenditure
    revenue_outside <- s0 / total_expenditure
    slopes <- list(
      gamma = gamma,
      meanval = (quantity_shares / s0) * prices^gamma,
      alpha = revenue_outside / (1 - revenue_outside)
    )
    class_name <- "BayesConductCES"
  }

  object <- new(class_name,
    prices = prices,
    shares = if (demand == "ces") revenue_shares else quantity_shares,
    margins = rep(0.2, 3),
    weights = rep(1, 3),
    ownerPre = owner_pre,
    ownerPost = owner_post,
    mcDelta = rep(0, 3),
    subset = rep(TRUE, 3),
    insideSize = 1,
    priceOutside = 1,
    priceStart = prices,
    shareInside = if (demand == "ces") 1 - revenue_outside else
      sum(quantity_shares),
    normIndex = NA,
    output = TRUE,
    labels = letters[1:3],
    supply_model = as.integer(supply_model),
    use_cutoff = use_cutoff,
    cutoff_share = cutoff_share,
    cutoff_scale = 0.2
  )
  object@slopes <- slopes
  object@ownerPre <- antitrust::ownerToMatrix(object, TRUE)
  object@ownerPost <- antitrust::ownerToMatrix(object, FALSE)
  object@pricePre <- prices
  object@pricePost <- prices
  object
}

test_that("conduct classes reproduce fitted Logit and CES quantity shares", {
  for (demand in c("logit", "ces")) {
    object <- make_conduct_object(demand)
    expect_equal(
      antitrust::calcShares(object),
      c(0.20, 0.15, 0.10),
      tolerance = 1e-12,
      ignore_attr = TRUE
    )
  }
})

test_that("Logit conduct margins match all four Stan FOCs", {
  alpha <- 2
  prices <- c(1, 1.2, 0.9)
  firm_share <- c(0.35, 0.35, 0.10)
  s0 <- 0.55
  conduct <- list(
    1 - firm_share,
    rep(1, 3),
    s0 / (s0 + firm_share),
    firm_share / -log1p(-firm_share)
  )

  for (supply_model in 1:4) {
    object <- make_conduct_object("logit", supply_model)
    expected <- 1 / (alpha * conduct[[supply_model]] * prices)
    expect_equal(
      antitrust::calcMargins(object), expected,
      tolerance = 1e-12, ignore_attr = TRUE
    )
  }
})

test_that("CES conduct margins match all four Stan FOCs", {
  object <- make_conduct_object("ces")
  revenue_share <- antitrust::calcShares(object, revenue = TRUE)
  firm_share <- as.vector(object@ownerPre %*% revenue_share)
  gamma <- 3
  margin_inv <- list(
    gamma - (gamma - 1) * firm_share,
    rep(gamma, 3),
    gamma / (1 + (gamma - 1) * firm_share),
    1 / (1 - (1 - firm_share)^(1 / (gamma - 1)))
  )

  for (supply_model in 1:4) {
    object@supply_model <- as.integer(supply_model)
    expect_equal(
      antitrust::calcMargins(object), 1 / margin_inv[[supply_model]],
      tolerance = 1e-12, ignore_attr = TRUE
    )
  }
})

test_that("auction conduct supports antitrust ex-ante margin arguments", {
  logit <- make_conduct_object("logit", supply_model = 4L)
  ces <- make_conduct_object("ces", supply_model = 4L)
  expect_equal(
    antitrust::calcMargins(logit, exAnte = TRUE),
    antitrust::calcMargins(logit, level = FALSE) * antitrust::calcShares(logit),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_equal(
    antitrust::calcMargins(ces, exAnte = TRUE),
    antitrust::calcMargins(ces, level = FALSE) *
      antitrust::calcShares(ces, revenue = TRUE),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
})

test_that("cutoff conduct is recomputed from post-merger firm shares", {
  object <- make_conduct_object(
    "logit", supply_model = 1L, use_cutoff = TRUE, cutoff_share = 0.4
  )
  pre <- antitrust::calcMargins(object, preMerger = TRUE, level = TRUE)
  post <- antitrust::calcMargins(object, preMerger = FALSE, level = TRUE)
  expect_false(isTRUE(all.equal(pre, post)))
})

test_that("monopolistic conduct is ownership-independent", {
  for (demand in c("logit", "ces")) {
    object <- make_conduct_object(demand, supply_model = 2L)
    expect_equal(
      antitrust::calcMargins(object, TRUE),
      antitrust::calcMargins(object, FALSE),
      tolerance = 1e-12
    )
  }
})

test_that("CES conduct objects retain antitrust demand methods", {
  object <- make_conduct_object("ces", supply_model = 1L)
  object@mcPre <- antitrust::calcMC(object, TRUE)
  object@mcPost <- antitrust::calcMC(object, FALSE)
  object@pricePost <- antitrust::calcPrices(object, FALSE)
  expect_true(is.matrix(antitrust::elast(object)))
  expect_true(is.matrix(antitrust::diversion(object)))
  expect_true(is.numeric(antitrust::CV(object)))
  expect_type(antitrust::HypoMonTest(object, prodIndex = 1:2), "logical")
})
