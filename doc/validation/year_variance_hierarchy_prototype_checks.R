# No HMC. Detect incorrect contrast prior, redundant year-level coordinates,
# shifted 1960 base-scale prior and inconsistency with quality conditioning.
library(testthat)
pkgload::load_all(".", quiet = TRUE)
source("doc/validation/demand_shock_prototype.R")
source("doc/validation/asset_variance_prototype.R")
source("doc/validation/inverse_asset_variance_prototype.R")
source("doc/validation/variance_time_trend_prototype.R")
source("doc/validation/year_variance_hierarchy_prototype.R")

test_that("five Gaussian contrasts preserve the declared year-effect prior", {
  h <- contr.helmert(6); h <- sweep(h, 2, sqrt(colSums(h^2)), "/")
  expect_equal(crossprod(h), diag(5), tolerance = 1e-14)
  expect_equal(unname(tcrossprod(h)), diag(6) - matrix(1 / 6, 6, 6), tolerance = 1e-14)
  expect_equal(colSums(h), rep(0, 5), tolerance = 1e-14)
  difference <- h[6, ] - h[1, ]
  expect_equal(sum(difference^2) * (1 / sqrt(2))^2, 1, tolerance = 1e-14)
  u <- as.numeric(.7 * h %*% c(.4, -.2, .6, -.1, .3))
  offset <- u - u[1]
  expect_identical(offset[1], 0)
  expect_equal(exp(.5 * offset)^2, exp(offset), tolerance = 1e-14)
})

test_that("hierarchy adds only five noncentered effects and one proper SD prior", {
  base <- tempfile(fileext = ".stan"); path <- tempfile(fileext = ".stan")
  on.exit(unlink(c(base, path)), add = TRUE)
  iv_kernel("inst/stan/logit_supply.stan", base)
  yh_kernel("inst/stan/logit_supply.stan", path)
  a <- paste(readLines(base), collapse = "\n")
  b <- paste(readLines(path), collapse = "\n")
  block <- function(s) {
    first <- regexpr("\nparameters {", s, fixed = TRUE)[[1L]]
    last <- regexpr("\ntransformed parameters {", s, fixed = TRUE)[[1L]]
    gsub("[[:space:]]", "", substring(s, first, last - 1L))
  }
  removed <- sub("  vector[year_variance_active * (N_variance_year-1)] year_variance_raw;", "", b, fixed = TRUE)
  removed <- sub("  vector<lower=0>[year_variance_active] year_logvar_sd;", "", removed, fixed = TRUE)
  expect_identical(block(a), block(removed))
  expect_match(b, "std_normal_lpdf(year_variance_raw)", fixed = TRUE)
  expect_match(b, "normal_lccdf(0 | 0, year_logvar_prior_scale)", fixed = TRUE)
  expect_match(b, "0.5 * year_log_variance_offset[variance_year_idx[i]]", fixed = TRUE)
  expect_false(grepl("variance_time_gamma", b, fixed = TRUE))
})

test_that("arbitrary year variance effects retain dense Gaussian quality recovery", {
  cell <- rep(1:3, each = 3L); group <- rep(1:3, 3L)
  basis <- antitrustBayes:::.relative_share_basis(cell)
  z <- diag(3)[group, , drop = FALSE]
  asset <- rep(c(10, 20, 80), 3)
  offset <- rep(c(0, -.6, .9), each = 3)
  sigma <- .25 * sqrt(30 / asset) * exp(.5 * offset)
  tau <- .4; rho <- -.35; sigma_m <- .3
  paired <- c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
  m <- c(.04, 0, -.05, .1, 0, -.02, 0, .08, -.07)
  residual <- c(-.3, .1, .4, -.2, .15, .5, -.1, .2, .3)
  shift <- rho * sigma / sigma_m * m
  variance <- sigma^2 * (1 - rho^2 * paired)
  prior <- tau^2 * (diag(3) - matrix(1 / 3, 3, 3))
  covariance <- crossprod(basis, (z %*% prior %*% t(z) + diag(variance)) %*% basis)
  expected <- prior %*% t(z) %*% basis %*% solve(covariance, crossprod(basis, residual - shift))
  data <- list(N = 9L, margin_obs = paired, log_marginNormInv = log(rep(2, 9) + 1e-9) + m,
    log_shareIn = residual, pred_margin = rep(2, 9), quality_block_n = 3L,
    obs_to_quality = group, market_year_idx = cell)
  local_mocked_bindings(rnorm = function(n, ...) rep(0, n), .package = "stats")
  q <- av_quality_state(data, rep(0, 9), tau, sigma, sigma_m, rho)
  expect_equal(q$persistent_quality, as.numeric(expected), tolerance = 1e-12)
  expect_equal(q$structural_state, as.numeric(expected)[group], tolerance = 1e-12)
  expect_identical(q$period_shock, rep(0, 9))
})
