smi_runner <- function() {
  path <- testthat::test_path("..", "..", "doc", "validation", "single_market_intercept_validation.R")
  if (!file.exists(path)) testthat::skip("Research runner is excluded from the source package.")
  environment <- new.env(parent = globalenv()); sys.source(path, envir = environment); environment
}

test_that("free Logit intercept supports outside shares beyond the old seven-firm ceiling", {
  r <- c(.25, .20, .16, .13, .11, .09, .06)
  prices <- c(1.2, 1.4, 1.6, 1.9, 2.1, 2.4, 2.8)
  for (s0 in c(.1, .2, .4)) {
    a <- mean(log(r)) - qlogis(s0)
    quality <- log(r) - mean(log(r)) + 2.4 * (prices - mean(prices))
    delta <- a + quality - 2.4 * (prices - mean(prices))
    shares <- exp(delta) / (1 + sum(exp(delta)))
    expect_equal(mean(quality), 0, tolerance = 1e-14)
    expect_equal(1 - sum(shares), s0, tolerance = 1e-14)
    expect_equal(shares / sum(shares), r, tolerance = 1e-14)
    expect_equal(shares, (1 - s0) * r, tolerance = 1e-14)
    expect_equal(mean(delta), a)
  }
  expect_lt(plogis(mean(log(r))), 1 / (length(r) + 1))
  expect_equal(plogis(mean(log(rep(1 / 7, 7)))), 1 / 8)
  expect_equal(mean(log(r)) - qlogis(plogis(mean(log(r)))), 0, tolerance = 1e-14)
})

test_that("MC margins leave outside-share identification to the prior and share level", {
  r <- c(.25, .20, .16, .13, .11, .09, .06)
  for (s0 in c(.1, .2, .4)) {
    inv_markup <- rep(2.4, length(r)) # MC FOC: inverse dollar markup is alpha.
    expect_equal(inv_markup, rep(2.4, 7))
  }
  conditional_mean <- function(prior_mean) {
    # Actual bounded link in the zero-share-noise conditional level equation.
    kernel <- function(eta) dnorm(mean(log(r)) - qlogis(1e-6 + (1 - 2e-6) * plogis(eta))) *
      dnorm(eta, prior_mean, .5)
    integrate(function(eta) eta * kernel(eta), -Inf, Inf)$value / integrate(kernel, -Inf, Inf)$value
  }
  expect_gt(conditional_mean(qlogis(.4)) - conditional_mean(qlogis(.1)), 1)
  expect_equal(integrate(function(a) exp(dnorm(a, log = TRUE)), -Inf, Inf)$value, 1)
  # A duplicated intercept prior changes the density ratio and integrated mass.
  expect_equal(integrate(function(a) exp(2 * (dnorm(a, log = TRUE) - dnorm(0, log = TRUE))) * dnorm(0),
    -Inf, Inf)$value, 1 / sqrt(2))
})

test_that("heterogeneous Cournot margins identify alpha and outside share independently of the intercept", {
  r <- c(.25, .20, .16, .13, .11, .09, .06)
  alpha <- 2.4; s0 <- .2
  inverse_markup <- alpha * s0 / (s0 + (1 - s0) * r)
  line <- stats::coef(stats::lm(I(1 / inverse_markup) ~ r))
  expect_equal(unname(1 / line[[1L]]), alpha, tolerance = 1e-12)
  expect_equal(unname(1 / (1 + line[[2L]] / line[[1L]])), s0, tolerance = 1e-12)
  profile <- function(candidate) {
    f <- candidate / (candidate + (1 - candidate) * r)
    fitted_alpha <- sum(f * inverse_markup) / sum(f^2)
    sum((inverse_markup - fitted_alpha * f)^2)
  }
  expect_equal(profile(s0), 0, tolerance = 1e-14)
  for (wrong in c(.1, .4)) expect_gt(profile(wrong), .01)
})

test_that("fresh-target protocol is bounded and preserves all four market-time layouts", {
  runner <- smi_runner(); design <- runner$smi_design()
  fits <- Filter(function(task) task$stage == "fit", design$tasks)
  expect_length(fits, 12L)
  expect_equal(design$controls$chains, 4L)
  expect_equal(design$controls$iter - design$controls$warmup, 500L)
  expect_setequal(design$fixtures$id[design$fixtures$conduct == "C"], c("C_M1T1", "C_M1T4", "C_M4T1", "C_M4T4"))
  expect_length(Filter(function(task) task$stage == "bridge", design$tasks), 6L)
  expect_equal(design$k, 3L)
  expect_equal(design$outside_share, .2)
})
