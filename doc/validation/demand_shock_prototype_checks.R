# No HMC. A material covariance failure would turn conditional shocks into
# independent firm noise or pair them with the wrong persistent-quality state.
library(testthat)
pkgload::load_all(".", quiet = TRUE)
e <- new.env(parent = globalenv())
source("doc/validation/demand_shock_prototype.R", local = e)

test_that("joint contrast recovery equals a dense posterior conditional mean", {
  cell <- rep(1:2, each = 3L); quality <- rep(1:3, 2L)
  basis <- antitrustBayes:::.relative_share_basis(cell)
  projection <- tcrossprod(basis)
  z <- diag(3)[quality, , drop = FALSE]
  tau <- .35; omega <- .22; rho <- .6; margin_sd <- .17
  paired <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  margin_residual <- c(.03, 0, -.07, 0, .09, 0)
  y <- c(-.2, .1, .3, -.1, .25, .05)
  shift <- rho * omega / margin_sd * margin_residual
  variance <- omega^2 * (1 - rho^2 * paired)
  q_covariance <- tau^2 * (diag(3) - matrix(1 / 3, 3, 3))
  covariance <- crossprod(basis, (z %*% q_covariance %*% t(z) +
    diag(variance + .01)) %*% basis)
  precision <- basis %*% solve(covariance, t(basis))
  expected_q <- q_covariance %*% t(z) %*% precision %*% (y - shift)
  expected_x <- projection %*% (shift + diag(variance) %*% precision %*% (y - shift))
  # Zero innovations select the exact conditional mean in the Gaussian
  # Matheron draw, allowing a deterministic independent identity test.
  local_mocked_bindings(rnorm = function(n, ...) rep(0, n), .package = "stats")
  e$rnorm <- function(n, ...) rep(0, n)
  data <- list(N = 6L, margin_obs = paired,
    log_marginNormInv = log(rep(2, 6) + 1e-9) + margin_residual,
    log_shareIn = array(y, 6L), quality_block_n = 3L,
    obs_to_quality = array(quality, 6L), market_year_idx = array(cell, 6L))
  state <- e$ds_joint_state(data, rep(0, 6), tau, omega, margin_sd, rho, rep(2, 6))
  expect_equal(state$persistent_quality, as.numeric(expected_q), tolerance = 1e-12)
  expect_equal(state$period_shock, as.numeric(expected_x), tolerance = 1e-12)
  expect_equal(state$structural_state,
    as.numeric(expected_q)[quality] + as.numeric(expected_x), tolerance = 1e-12)
  expect_equal(as.numeric(tapply(state$period_shock, cell, sum)), c(0, 0), tolerance = 1e-12)
  expect_gt(max(abs(state$period_shock)), .01)
  expect_gt(max(abs(state$structural_state - y)), .01)
})

test_that("fixed measurement and shock-margin covariance retain the intended target", {
  omega <- .22; sigma_margin <- .17; rho <- .6
  total <- sqrt(omega^2 + .1^2)
  effective <- rho * omega / total
  expect_equal(total * effective / sigma_margin, rho * omega / sigma_margin)
  expect_equal(total^2 * (1 - effective^2), omega^2 * (1 - rho^2) + .1^2)
  expect_lt(abs(effective), abs(rho))
})
