# No HMC. Named failures: outcome-dependent/case-repeated asset weighting,
# wrong heteroskedastic share-margin covariance, unnormalized contrast density,
# and conditional quality inconsistent with the estimation likelihood.
library(testthat)
pkgload::load_all(".", quiet = TRUE)
source("doc/validation/demand_shock_prototype.R")
source("doc/validation/asset_variance_prototype.R")

test_that("institution asset scaling is global and invariant to case duplication", {
  rows <- data.frame(tophold = c("a", "b", "a", "c"), year = c(1, 1, 1, 2),
    total_assets = c(10, 100, 10, 1000), market = c(1, 1, 2, 2),
    total_deposits = c(2, 20, 1, 200))
  x <- av_asset_index(rows)
  expect_equal(nrow(x$rows), 3L)
  expect_equal(x$rows$institution_log_assets_z, c(-1, 0, 1), tolerance = 1e-14)
  repeat_case <- rbind(rows, rows)
  repeat_case$total_deposits <- seq_len(nrow(repeat_case))
  expect_identical(av_asset_index(repeat_case), x)
  dollars <- rows; dollars$total_assets <- dollars$total_assets * 1000
  expect_equal(av_asset_index(dollars)$rows$institution_log_assets_z,
    x$rows$institution_log_assets_z, tolerance = 1e-14)
  conflicting <- rows; conflicting$total_assets[3L] <- 11
  expect_error(av_asset_index(conflicting), "disagree across cases")
})

test_that("heteroskedastic joint contrast likelihood equals a dense Gaussian", {
  cell <- rep(1:2, each = 3L); quality <- rep(1:3, 2L)
  basis <- antitrustBayes:::.relative_share_basis(cell)
  z <- diag(3)[quality, , drop = FALSE]
  tau <- .35; sigma <- .22 * exp(-.3 * c(-1, 0, 1, -.8, .1, 1.2))
  rho <- .6; sigma_m <- .17
  paired <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  margin_residual <- c(.03, 0, -.07, 0, .09, 0)
  residual <- c(-.2, .1, .3, -.1, .25, .05)
  shift <- rho * sigma / sigma_m * margin_residual
  variance <- sigma^2 * (1 - rho^2 * paired)
  q_covariance <- tau^2 * (diag(3) - matrix(1 / 3, 3, 3))
  covariance <- crossprod(basis, (z %*% q_covariance %*% t(z) + diag(variance)) %*% basis)
  y <- as.numeric(crossprod(basis, residual - shift))
  dense <- -.5 * (length(y) * log(2 * pi) + determinant(covariance, logarithm = TRUE)$modulus +
    crossprod(y, solve(covariance, y)))
  actual <- antitrustBayes:::.relative_share_gaussian(residual - shift, variance,
    quality, rep(1L, 3), cell, tau, return_moments = TRUE)
  expect_equal(actual$log_density, as.numeric(dense), tolerance = 1e-12)
  expected_q <- q_covariance %*% t(z) %*% basis %*% solve(covariance, y)
  expect_equal(actual$conditional_mean, as.numeric(expected_q), tolerance = 1e-12)
  data <- list(N = 6L, margin_obs = paired,
    log_marginNormInv = log(rep(2, 6) + 1e-9) + margin_residual,
    log_shareIn = residual, pred_margin = rep(2, 6), quality_block_n = 3L,
    obs_to_quality = quality, market_year_idx = cell)
  local_mocked_bindings(rnorm = function(n, ...) rep(0, n), .package = "stats")
  state <- av_quality_state(data, rep(0, 6), tau, sigma, sigma_m, rho)
  expect_equal(state$persistent_quality, as.numeric(expected_q), tolerance = 1e-12)
  expect_equal(state$structural_state, as.numeric(expected_q)[quality], tolerance = 1e-12)
  expect_identical(state$period_shock, rep(0, 6))
})

test_that("gamma zero retains original same-seed conditional quality draws", {
  data <- list(N = 6L, margin_obs = c(1L, 0L, 1L, 1L, 1L, 0L),
    log_marginNormInv = seq(-.2, .3, length.out = 6),
    log_shareIn = seq(-.5, .5, length.out = 6), pred_margin = rep(1, 6),
    quality_block_n = 3L, obs_to_quality = rep(1:3, 2), market_year_idx = rep(1:2, each = 3))
  set.seed(123); original <- ds_reference_state(data, rep(.02, 6), .35, .2, .17, .4)
  set.seed(123); candidate <- av_quality_state(data, rep(.02, 6), .35, rep(.2, 6), .17, .4)
  expect_identical(candidate, original)
})

test_that("asset kernel and policy evaluator retain the audited covariance routes", {
  file <- tempfile(fileext = ".stan")
  on.exit(unlink(file), add = TRUE)
  av_kernel("inst/stan/logit_supply.stan", file)
  text <- paste(readLines(file), collapse = "\n")
  count <- function(x) length(gregexpr(x, text, fixed = TRUE)[[1L]])
  expect_equal(count("relative_base_var[i] = square(sigma_share_observation[i])"), 2L)
  expect_equal(count("relative_cond_resid[i] -= rho * sigma_share_observation[i]"), 2L)
  expect_match(text, "normal_lpdf(asset_gamma | 0, asset_gamma_prior_sd)", fixed = TRUE)
  evaluator <- av_make_evaluator()
  code <- paste(deparse(evaluator), collapse = "\n")
  expect_match(code, "av_quality_state", fixed = TRUE)
  expect_false(grepl("ds_joint_state", code, fixed = TRUE))
  expect_match(code, "sigma_share_observation", fixed = TRUE)
})
