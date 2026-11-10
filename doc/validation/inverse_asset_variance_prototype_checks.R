# No HMC. Failures: inverse SD mistaken for inverse variance, standardized
# asset exponent with wrong units, new sampler/prior coordinates, or a changed
# conditional quality covariance relative to the likelihood.
library(testthat)
pkgload::load_all(".", quiet = TRUE)
source("doc/validation/demand_shock_prototype.R")
source("doc/validation/asset_variance_prototype.R")
source("doc/validation/inverse_asset_variance_prototype.R")

test_that("fixed institution rule is inverse variance with unchanged units", {
  assets <- c(10, 20, 40, 160)
  center <- mean(log(assets)); scale <- sd(log(assets))
  z <- (log(assets) - center) / scale
  variance <- (.3 * exp(-.5 * scale * z))^2
  expect_equal(variance * assets, rep(.3^2 * exp(center), 4L), tolerance = 1e-14)
  expect_equal(variance[1] / variance[2], 2, tolerance = 1e-14)
  expect_equal(sqrt(variance[1] / variance[2]), sqrt(2), tolerance = 1e-14)
  dollars <- assets * 1000
  expect_equal(exp(mean(log(dollars))) / dollars, exp(center) / assets, tolerance = 1e-14)
})

test_that("inverse variance adds no HMC parameter and changes no prior", {
  file <- tempfile(fileext = ".stan"); on.exit(unlink(file), add = TRUE)
  iv_kernel("inst/stan/logit_supply.stan", file)
  original <- paste(readLines("inst/stan/logit_supply.stan"), collapse = "\n")
  text <- paste(readLines(file), collapse = "\n")
  block <- function(s) {
    start <- regexpr("\nparameters {", s, fixed = TRUE)[[1L]]
    end <- regexpr("\ntransformed parameters {", s, fixed = TRUE)[[1L]]
    gsub("[[:space:]]", "", substring(s, start, end - 1L))
  }
  expect_identical(block(text), block(original))
  expect_false(grepl("normal_lpdf(asset_gamma", text, fixed = TRUE))
  expect_false(grepl("asset_gamma_prior_sd", text, fixed = TRUE))
  expect_match(text, "-0.5 * institution_log_assets_log_sd", fixed = TRUE)
  expect_match(text, "asset_gamma * institution_log_assets_z[i]", fixed = TRUE)
  expect_match(text, "relative_base_var[i] = square(sigma_share_observation[i])", fixed = TRUE)
  expect_match(text, "relative_cond_resid[i] -= rho * sigma_share_observation[i]", fixed = TRUE)
})

test_that("inverse-asset conditional quality matches a dense Gaussian posterior", {
  cell <- rep(1:2, each = 3L); group <- rep(1:3, 2L)
  basis <- antitrustBayes:::.relative_share_basis(cell)
  z <- diag(3)[group, , drop = FALSE]
  tau <- .4; assets <- c(10, 20, 80, 12, 25, 100)
  sigma <- .25 * sqrt(exp(mean(log(assets))) / assets)
  rho <- -.45; sigma_m <- .3
  paired <- c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
  margin_resid <- c(.04, 0, -.05, .1, -.02, 0)
  residual <- c(-.3, .1, .4, -.2, .15, .5)
  shift <- rho * sigma / sigma_m * margin_resid
  variance <- sigma^2 * (1 - rho^2 * paired)
  prior <- tau^2 * (diag(3) - matrix(1 / 3, 3, 3))
  covariance <- crossprod(basis, (z %*% prior %*% t(z) + diag(variance)) %*% basis)
  expected <- prior %*% t(z) %*% basis %*% solve(covariance, crossprod(basis, residual - shift))
  data <- list(N = 6L, margin_obs = paired,
    log_marginNormInv = log(rep(2, 6) + 1e-9) + margin_resid,
    log_shareIn = residual, pred_margin = rep(2, 6), quality_block_n = 3L,
    obs_to_quality = group, market_year_idx = cell)
  local_mocked_bindings(rnorm = function(n, ...) rep(0, n), .package = "stats")
  q <- av_quality_state(data, rep(0, 6), tau, sigma, sigma_m, rho)
  expect_equal(q$persistent_quality, as.numeric(expected), tolerance = 1e-12)
  expect_equal(q$structural_state, as.numeric(expected)[group], tolerance = 1e-12)
  expect_identical(q$period_shock, rep(0, 6))
})
