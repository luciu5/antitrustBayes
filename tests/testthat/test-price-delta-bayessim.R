# calcPriceDelta,BayesSim-method: Paasche/Laspeyres deposit-rate-style index
# computed per posterior draw, delegating the arithmetic to
# antitrust::calcPriceDelta on each draw's terminal structural model.
#
# Fixtures below use the "Linear" demand class with slopes = 0, so
# calcQuantities()/calcShares() return the intercepts (normalized) exactly,
# regardless of price -- this lets a test pin down q_post exactly while
# pricePre/pricePost are set independently, without needing a solved
# equilibrium.

linear_fixture <- function(pricePre, pricePost, quantityPost, quantityPre = quantityPost,
                            ownerPre = NULL, ownerPost = NULL) {
  k <- length(pricePre)
  labels <- make.names(seq_len(k))
  if (is.null(ownerPre)) ownerPre <- diag(k)
  if (is.null(ownerPost)) ownerPost <- matrix(1, k, k)
  diversion <- matrix(1, k, k) - diag(2, k)
  methods::new("Linear",
    pricePre = pricePre, pricePost = pricePost,
    ownerPre = ownerPre, ownerPost = ownerPost,
    mcPre = pricePre / 2, mcPost = pricePost / 2,
    output = TRUE, labels = labels,
    shares = rep(1 / k, k), mcDelta = rep(0, k), subset = rep(TRUE, k),
    diversion = diversion,
    intercepts = quantityPost, slopes = matrix(0, k, k),
    prices = pricePre, quantities = quantityPost, margins = rep(0.5, k),
    priceStart = pricePre, symmetry = TRUE
  )
}

bayes_sim_fixture <- function(models, price_pre, price_post, shares_pre, shares_post,
                               market_id = "mkt1") {
  n_draws <- length(models)
  k <- length(price_pre[[1]])
  to_mat <- function(rows) {
    m <- do.call(rbind, rows)
    dimnames(m) <- NULL
    m
  }
  methods::new("BayesSim",
    sims = setNames(list(models), market_id),
    price_pre = setNames(list(to_mat(price_pre)), market_id),
    price_post = setNames(list(to_mat(price_post)), market_id),
    price_change = setNames(list(to_mat(Map(function(a, b) b / a - 1, price_pre, price_post))), market_id),
    shares_pre = setNames(list(to_mat(shares_pre)), market_id),
    shares_post = setNames(list(to_mat(shares_post)), market_id),
    market_ids = market_id,
    product_labels = setNames(list(make.names(seq_len(k))), market_id),
    n_draws = as.integer(n_draws),
    n_valid = setNames(as.integer(sum(!vapply(models, is.null, logical(1)))), market_id),
    failures = list(),
    call = quote(simulate(object)),
    metadata = list()
  )
}

test_that("Paasche matches the brief's hand calculation and antitrust's canonical method", {
  fixture <- linear_fixture(
    pricePre = c(1, 2), pricePost = c(0.9, 1.6), quantityPost = c(80, 20)
  )
  # Canonical antitrust reference: shares normalize (80,20) -> (0.8,0.2);
  # num = 0.8*0.9 + 0.2*1.6 = 1.04; den = 0.8*1 + 0.2*2 = 1.2;
  # 1.04/1.2 - 1 = -0.1333333333...
  canonical <- antitrust::calcPriceDelta(fixture, market = TRUE, index = "paasche")
  expect_equal(canonical, -0.1333333333, tolerance = 1e-9)

  sim <- bayes_sim_fixture(
    models = list(fixture),
    price_pre = list(c(1, 2)), price_post = list(c(0.9, 1.6)),
    shares_pre = list(c(0.8, 0.2)), shares_post = list(c(0.8, 0.2))
  )
  draws <- calcPriceDelta(sim, market = TRUE, index = "paasche", posterior = TRUE)
  expect_length(draws, 1L)
  expect_equal(unname(draws[1]), canonical, tolerance = 1e-12)
  expect_equal(attr(draws, "n_valid"), 1L)

  summary_out <- calcPriceDelta(sim, market = TRUE, index = "paasche")
  expect_equal(unname(summary_out["mean"]), canonical, tolerance = 1e-12)
  expect_equal(unname(summary_out["n_valid"]), 1)

  # This case must distinguish the Paasche index from the equal-bank mean
  # (-0.15) and the quantity-weighted mean of proportional changes (-0.12).
  equal_bank_mean <- mean(c(0.9 / 1 - 1, 1.6 / 2 - 1))
  weighted_mean_of_changes <- sum(c(0.8, 0.2) * c(0.9 / 1 - 1, 1.6 / 2 - 1))
  expect_equal(equal_bank_mean, -0.15, tolerance = 1e-12)
  expect_equal(weighted_mean_of_changes, -0.12, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(canonical, equal_bank_mean)))
  expect_false(isTRUE(all.equal(canonical, weighted_mean_of_changes)))
})

test_that("Paasche is invariant to a common rescaling of post-merger quantities", {
  small <- linear_fixture(c(1, 2), c(0.9, 1.6), quantityPost = c(80, 20))
  large <- linear_fixture(c(1, 2), c(0.9, 1.6), quantityPost = c(8000, 2000))
  sim_small <- bayes_sim_fixture(list(small), list(c(1, 2)), list(c(0.9, 1.6)),
    list(c(0.8, 0.2)), list(c(0.8, 0.2)))
  sim_large <- bayes_sim_fixture(list(large), list(c(1, 2)), list(c(0.9, 1.6)),
    list(c(0.8, 0.2)), list(c(0.8, 0.2)))
  r_small <- calcPriceDelta(sim_small, market = TRUE, index = "paasche", posterior = TRUE)
  r_large <- calcPriceDelta(sim_large, market = TRUE, index = "paasche", posterior = TRUE)
  expect_equal(unname(r_small[1]), unname(r_large[1]), tolerance = 1e-12)
})

test_that("Paasche is exactly zero with no rate change and matches a common proportional change", {
  no_change <- linear_fixture(c(1, 2), c(1, 2), quantityPost = c(60, 40))
  sim <- bayes_sim_fixture(list(no_change), list(c(1, 2)), list(c(1, 2)),
    list(c(0.6, 0.4)), list(c(0.6, 0.4)))
  r <- calcPriceDelta(sim, market = TRUE, index = "paasche", posterior = TRUE)
  expect_equal(unname(r[1]), 0, tolerance = 1e-14)

  proportional <- linear_fixture(c(1, 2), c(1.05, 2.10), quantityPost = c(60, 40))
  sim2 <- bayes_sim_fixture(list(proportional), list(c(1, 2)), list(c(1.05, 2.10)),
    list(c(0.6, 0.4)), list(c(0.6, 0.4)))
  r2 <- calcPriceDelta(sim2, market = TRUE, index = "paasche", posterior = TRUE)
  expect_equal(unname(r2[1]), 0.05, tolerance = 1e-12)
})

test_that("MC's exact zero-effect benchmark is not altered by this change", {
  identical_fixture <- linear_fixture(c(1, 1), c(1, 1), quantityPost = c(50, 50))
  sim <- bayes_sim_fixture(list(identical_fixture), list(c(1, 1)), list(c(1, 1)),
    list(c(0.5, 0.5)), list(c(0.5, 0.5)))
  r <- calcPriceDelta(sim, market = TRUE, index = "paasche", posterior = TRUE)
  expect_equal(unname(r[1]), 0, tolerance = 1e-14)
})

test_that("finite negative post-merger prices are retained, not filtered", {
  fixture <- linear_fixture(c(1, 2), c(-0.5, 1.6), quantityPost = c(80, 20))
  sim <- bayes_sim_fixture(list(fixture), list(c(1, 2)), list(c(-0.5, 1.6)),
    list(c(0.8, 0.2)), list(c(0.8, 0.2)))
  r <- calcPriceDelta(sim, market = TRUE, index = "paasche", posterior = TRUE)
  expect_true(is.finite(r[1]))
  # num = 0.8*(-0.5) + 0.2*1.6 = -0.08; den = 0.8*1+0.2*2 = 1.2
  expect_equal(unname(r[1]), -0.08 / 1.2 - 1, tolerance = 1e-9)
})

test_that("invalid draws are recorded as NA with a reason, never silently dropped", {
  valid <- linear_fixture(c(1, 2), c(0.9, 1.6), quantityPost = c(80, 20))
  fixtures <- list(valid, valid, valid)
  price_pre <- list(c(1, 2), c(1, 2), c(0, 0))       # draw 3: zero weighted denominator
  price_post <- list(c(0.9, 1.6), c(0.9, 1.6), c(0.9, 1.6))
  shares_pre <- list(c(0.8, 0.2), c(0.8, 0.2), c(0.8, 0.2))
  shares_post <- list(c(0.8, 0.2), c(-0.1, 1.1), c(0.8, 0.2))  # draw 2: negative weight

  sim <- bayes_sim_fixture(fixtures, price_pre, price_post, shares_pre, shares_post)
  draws <- calcPriceDelta(sim, market = TRUE, index = "paasche", posterior = TRUE)
  expect_equal(unname(is.finite(draws)), c(TRUE, FALSE, FALSE))
  expect_equal(attr(draws, "n_valid"), 1L)
  reasons <- attr(draws, "failure_reason")
  expect_true(is.na(reasons[1]))
  expect_false(is.na(reasons[2]))
  expect_false(is.na(reasons[3]))

  # A draw whose simulation already failed (NULL sim) is also NA, not dropped.
  fixtures_with_failure <- list(valid, NULL)
  price_pre2 <- list(c(1, 2), c(1, 2))
  price_post2 <- list(c(0.9, 1.6), c(0.9, 1.6))
  shares2 <- list(c(0.8, 0.2), c(0.8, 0.2))
  sim2 <- bayes_sim_fixture(fixtures_with_failure, price_pre2, price_post2, shares2, shares2)
  draws2 <- calcPriceDelta(sim2, market = TRUE, index = "paasche", posterior = TRUE)
  expect_equal(length(draws2), 2L)
  expect_true(is.finite(draws2[1]))
  expect_true(is.na(draws2[2]))
  expect_equal(attr(draws2, "n_valid"), 1L)
})

test_that("product/draw misalignment between stored matrices and sims is rejected", {
  valid <- linear_fixture(c(1, 2), c(0.9, 1.6), quantityPost = c(80, 20))
  sim <- bayes_sim_fixture(list(valid), list(c(1, 2)), list(c(0.9, 1.6)),
    list(c(0.8, 0.2)), list(c(0.8, 0.2)))
  # Corrupt the stored price_pre matrix to have an extra draw row.
  sim@price_pre[["mkt1"]] <- rbind(sim@price_pre[["mkt1"]], c(1, 2))
  expect_error(calcPriceDelta(sim, market = TRUE, index = "paasche"), "misaligned")
})

test_that("calcPriceDelta,BayesSim-method agrees with antitrust on laspeyres and levels", {
  fixture <- linear_fixture(c(1, 2), c(0.9, 1.6), quantityPost = c(80, 20), quantityPre = c(60, 40))
  sim <- bayes_sim_fixture(list(fixture), list(c(1, 2)), list(c(0.9, 1.6)),
    list(c(0.6, 0.4)), list(c(0.8, 0.2)))

  canonical_laspeyres <- antitrust::calcPriceDelta(fixture, market = TRUE, index = "laspeyres")
  r_laspeyres <- calcPriceDelta(sim, market = TRUE, index = "laspeyres", posterior = TRUE)
  expect_equal(unname(r_laspeyres[1]), canonical_laspeyres, tolerance = 1e-12)

  canonical_levels <- antitrust::calcPriceDelta(fixture, levels = TRUE, market = TRUE, index = "paasche")
  r_levels <- calcPriceDelta(sim, levels = TRUE, market = TRUE, index = "paasche", posterior = TRUE)
  expect_equal(unname(r_levels[1]), canonical_levels, tolerance = 1e-12)
})

test_that("product-level (market = FALSE) output matches antitrust::calcPriceDelta per product", {
  fixture <- linear_fixture(c(1, 2), c(0.9, 1.6), quantityPost = c(80, 20))
  sim <- bayes_sim_fixture(list(fixture), list(c(1, 2)), list(c(0.9, 1.6)),
    list(c(0.8, 0.2)), list(c(0.8, 0.2)))
  canonical <- antitrust::calcPriceDelta(fixture, market = FALSE)
  out <- calcPriceDelta(sim, market = FALSE)
  expect_equal(dim(out), c(1L, 2L))
  expect_equal(unname(out[1, ]), unname(canonical), tolerance = 1e-12)
})
