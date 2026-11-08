test_that("target-snapshot scaling preserves the economic slope prior across sets", {
  snapshot <- c(1, 2, 4)
  anchor <- list(center = mean(snapshot), scale = sd(snapshot))
  prepare <- function(prices, specification = NULL) prepare_stan_data(
    prices = prices, quantities = rep(1, length(prices)),
    margins = rep(.2, length(prices)), market = rep("m", length(prices)),
    year = rep(seq_len(length(prices)/3), each = 3), group = rep(letters[1:3], length(prices)/3),
    supply_model = 1L, market_side = "output", share_basis = "relative",
    price_standardization = specification)
  a <- prepare(snapshot, anchor)
  b <- prepare(c(snapshot, 2 * snapshot), anchor)
  default <- prepare(c(snapshot, 2 * snapshot))
  expect_equal(a$price_sd, b$price_sd)
  expect_equal(a$price_std, b$price_std[1:3])
  expect_equal(as.numeric(b$price_std), (c(snapshot, 2*snapshot) - anchor$center)/anchor$scale)
  expect_equal(default$price_sd, sd(c(snapshot, 2*snapshot)))
  expect_false(isTRUE(all.equal(default$price_sd, b$price_sd)))
  # In native units both the Logit index differences and inverse wedges use
  # alpha = a / scale; centering cancels from within-cell contrasts.
  alpha <- 1.7
  standardized_slope <- alpha * b$price_sd
  expect_equal(as.numeric(standardized_slope * diff(b$price_std[1:3])), alpha * diff(snapshot))
  expect_equal(standardized_slope/b$price_sd * (1-c(.1,.2,.3)),
    alpha * (1-c(.1,.2,.3)))
  expect_equal(attr(b, "price_standardization"),
    c(anchor, list(source = "supplied")))
  expect_equal(a$prior_alpha_mean, b$prior_alpha_mean)
  expect_equal(a$prior_alpha_sd, b$prior_alpha_sd)
  expect_false(identical(digest::digest(a), digest::digest(b)))
})

test_that("invalid or unsupported price anchors fail before sampling", {
  for (bad in list(list(center = 0, scale = 0), list(center = NA_real_, scale = 1),
      list(center = 0, scale = Inf), list(center = 0), c(center = 0, scale = 1)))
    expect_error(.resolve_price_standardization(1:3, bad), "price_standardization")
  expect_error(bayes_est(1:3, 1:3, rep(.2, 3), rep("m", 3), group = letters[1:3],
    demand_model = "ces", price_standardization = list(center = 2, scale = 1)),
    "relative-share Logit")
  expect_error(prepare_stan_data(1:3, 1:3, rep(.2, 3), rep("m", 3), rep(1,3),
    letters[1:3], supply_model = 1L, share_basis = "absolute",
    price_standardization = list(center = 2, scale = 1)), "relative-share Logit")
})
