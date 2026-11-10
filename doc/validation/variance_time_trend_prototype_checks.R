# No HMC. Detect log-SD/log-variance confusion, shifted1960 reference priors,
# extra latent sampler states, and quality conditioned on the wrong covariance.
library(testthat)
pkgload::load_all(".", quiet = TRUE)
source("doc/validation/demand_shock_prototype.R")
source("doc/validation/asset_variance_prototype.R")
source("doc/validation/inverse_asset_variance_prototype.R")
source("doc/validation/variance_time_trend_prototype.R")

test_that("time coefficient is full-window log variance with1960 unchanged", {
  assets <- c(10, 40, 20, 80); time <- c(0, 0, 1, 1)
  sigma <- .3 * sqrt(30 / assets) * exp(.5 * log(4) * time)
  base <- .3 * sqrt(30 / assets)
  expect_equal(sigma[1:2], base[1:2])
  expect_equal((sigma[3:4] / base[3:4])^2, c(4, 4))
  expect_equal(sigma[3:4] / base[3:4], c(2, 2))
  expect_equal(sigma[3]^2 / sigma[4]^2, assets[4] / assets[3])
})

test_that("kernel adds exactly one trend parameter and a normalized prior", {
  original <- tempfile(fileext = ".stan"); trend <- tempfile(fileext = ".stan")
  on.exit(unlink(c(original, trend)), add = TRUE)
  iv_kernel("inst/stan/logit_supply.stan", original)
  tv_kernel("inst/stan/logit_supply.stan", trend)
  a <- paste(readLines(original), collapse = "\n")
  b <- paste(readLines(trend), collapse = "\n")
  block <- function(s) {
    first <- regexpr("\nparameters {", s, fixed = TRUE)[[1L]]
    last <- regexpr("\ntransformed parameters {", s, fixed = TRUE)[[1L]]
    gsub("[[:space:]]", "", substring(s, first, last - 1L))
  }
  removed <- sub("  vector[variance_time_trend_active] variance_time_gamma;", "", b, fixed = TRUE)
  expect_identical(block(a), block(removed))
  expect_match(b, "0.5 * variance_time_gamma[1] * variance_time01[i]", fixed = TRUE)
  expect_match(b, "normal_lpdf(variance_time_gamma | 0, variance_time_prior_sd)", fixed = TRUE)
  expect_match(b, "relative_cond_resid[i] -= rho * sigma_share_observation[i]", fixed = TRUE)
})

test_that("joint quality recovery uses the full time-and-size covariance", {
  cell <- rep(1:3, each = 3L); group <- rep(1:3, 3L)
  basis <- antitrustBayes:::.relative_share_basis(cell)
  z <- diag(3)[group, , drop = FALSE]
  assets <- c(10, 20, 80, 12, 25, 100, 15, 30, 120)
  time <- rep(c(0, .5, 1), each = 3)
  sigma <- .25 * sqrt(exp(mean(log(assets))) / assets) * exp(.5 * .7 * time)
  tau <- .4; rho <- -.35; sigma_m <- .3
  paired <- c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
  m <- c(.04, 0, -.05, .1, 0, -.02, 0, .08, -.07)
  residual <- c(-.3, .1, .4, -.2, .15, .5, -.1, .2, .3)
  shift <- rho * sigma / sigma_m * m
  v <- sigma^2 * (1 - rho^2 * paired)
  prior <- tau^2 * (diag(3) - matrix(1 / 3, 3, 3))
  covariance <- crossprod(basis, (z %*% prior %*% t(z) + diag(v)) %*% basis)
  expected <- prior %*% t(z) %*% basis %*% solve(covariance, crossprod(basis, residual - shift))
  data <- list(N = 9L, margin_obs = paired,
    log_marginNormInv = log(rep(2, 9) + 1e-9) + m,
    log_shareIn = residual, pred_margin = rep(2, 9), quality_block_n = 3L,
    obs_to_quality = group, market_year_idx = cell)
  local_mocked_bindings(rnorm = function(n, ...) rep(0, n), .package = "stats")
  q <- av_quality_state(data, rep(0, 9), tau, sigma, sigma_m, rho)
  expect_equal(q$persistent_quality, as.numeric(expected), tolerance = 1e-12)
  expect_equal(q$structural_state, as.numeric(expected)[group], tolerance = 1e-12)
  expect_identical(q$period_shock, rep(0, 9))
})
