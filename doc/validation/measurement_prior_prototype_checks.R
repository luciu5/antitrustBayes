# Explicit extended-tier checks: variance allocation can change economic states
# despite an unchanged observable covariance. No sampling or production edits.
library(testthat)
root <- normalizePath(if (file.exists("DESCRIPTION")) "." else "../..", mustWork = TRUE)
pkgload::load_all(root, quiet = TRUE)
e <- new.env(parent = globalenv())
source(file.path(root, "doc/validation/demand_shock_prototype.R"), local = e)
source(file.path(root, "doc/validation/measurement_prior_prototype.R"), local = e)

test_that("estimated measurement is used in joint relative-contrast recovery", {
  cell <- rep(1:2, each = 3L); quality <- rep(1:3, 2L)
  basis <- antitrustBayes:::.relative_share_basis(cell)
  projection <- tcrossprod(basis); z <- diag(3)[quality, , drop = FALSE]
  tau <- .35; omega <- .22; rho <- .6; margin_sd <- .17
  paired <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  margin_residual <- c(.03, 0, -.07, 0, .09, 0)
  y <- c(-.2, .1, .3, -.1, .25, .05)
  shift <- rho * omega / margin_sd * margin_residual
  variance <- omega^2 * (1 - rho^2 * paired)
  q_covariance <- tau^2 * (diag(3) - matrix(1 / 3, 3, 3))
  local_mocked_bindings(rnorm = function(n, ...) rep(0, n), .package = "stats")
  e$rnorm <- function(n, ...) rep(0, n)
  data <- list(N = 6L, margin_obs = paired,
    log_marginNormInv = log(rep(2, 6) + 1e-9) + margin_residual,
    log_shareIn = array(y, 6L), quality_block_n = 3L,
    obs_to_quality = array(quality, 6L), market_year_idx = array(cell, 6L))
  states <- list()
  for (measurement in c(0, .04, .20)) {
    covariance <- crossprod(basis, (z %*% q_covariance %*% t(z) +
      diag(variance + measurement^2)) %*% basis)
    precision <- basis %*% solve(covariance, t(basis))
    expected_q <- q_covariance %*% t(z) %*% precision %*% (y - shift)
    expected_x <- projection %*% (shift + diag(variance) %*% precision %*% (y - shift))
    state <- e$mp_joint_state(data, rep(0, 6), tau, omega, margin_sd, rho,
      rep(2, 6), measurement_sd = measurement)
    expect_equal(state$persistent_quality, as.numeric(expected_q), tolerance = 1e-12)
    expect_equal(state$period_shock, as.numeric(expected_x), tolerance = 1e-12)
    expect_equal(state$structural_state, as.numeric(expected_q)[quality] +
      as.numeric(expected_x), tolerance = 1e-12)
    expect_equal(state$measurement_sd, measurement)
    expect_equal(as.numeric(tapply(state$period_shock, cell, sum)), c(0, 0), tolerance = 1e-12)
    if (measurement == 0) expect_equal(as.numeric(crossprod(basis,
      state$structural_state - y)), rep(0, ncol(basis)), tolerance = 1e-12)
    states[[as.character(measurement)]] <- state
  }
  expect_gt(max(abs(states[["0.04"]]$structural_state - states[["0.2"]]$structural_state)), .01)
  expect_error(e$mp_joint_state(data, rep(0, 6), tau, omega, margin_sd, rho,
    rep(2, 6), measurement_sd = -1))
})

test_that("free shock correlation leaves a variance-allocation likelihood ridge", {
  omega <- .14; measurement <- .10; rho <- .4
  total_variance <- omega^2 + measurement^2; covariance <- rho * omega
  other_measurement <- .12
  other_omega <- sqrt(total_variance - other_measurement^2)
  other_rho <- covariance / other_omega
  expect_lt(abs(other_rho), .95)
  expect_equal(omega^2 + measurement^2, other_omega^2 + other_measurement^2)
  expect_equal(rho * omega, other_rho * other_omega)
  # Independently integrate quality into a full dense Gaussian covariance.
  b <- antitrustBayes:::.relative_share_basis(rep(1:2, each = 3))
  z <- diag(3)[rep(1:3, 2), , drop = FALSE]
  q <- .3^2 * (diag(3) - matrix(1 / 3, 3, 3))
  y <- c(.3, -.2, .05, .1, -.4, .2); rm <- c(.2, 0, -.1, .15, 0, -.07)
  paired <- c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE)
  ll <- function(w, m, r) {
    shift <- r * w / .2 * rm
    v <- w^2 * (1 - r^2 * paired) + m^2
    c <- crossprod(b, (z %*% q %*% t(z) + diag(v)) %*% b)
    residual <- crossprod(b, y - shift)
    -.5 * (length(residual) * log(2 * pi) + determinant(c, logarithm = TRUE)$modulus +
      crossprod(residual, solve(c, residual)))
  }
  expect_equal(ll(omega, measurement, rho), ll(other_omega, other_measurement, other_rho), tolerance = 1e-12)
})

test_that("positive-scale prior is normalized on the primitive log parameter", {
  mu <- log(.1); sdlog <- .5
  expect_equal(integrate(function(m) exp(dnorm(log(m), mu, sdlog, log = TRUE) - log(m)),
    0, Inf)$value, 1, tolerance = 1e-8)
  m <- c(.04, .1, .2)
  expect_equal(dnorm(log(m), mu, sdlog, log = TRUE),
    dlnorm(m, mu, sdlog, log = TRUE) + log(m))
  expect_equal(qlnorm(c(.025, .975), mu, sdlog), c(.03753179, .2664408), tolerance = 1e-6)
})
