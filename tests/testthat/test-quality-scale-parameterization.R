test_that("quality-scale coordinate API is explicit and restricted", {
  resolve <- antitrustBayes:::.resolve_quality_scale_parameterization
  legacy <- resolve("legacy", "collapsed", "gaussian", FALSE)
  rotated <- resolve("total_fraction", "collapsed", "gaussian", FALSE)
  expect_equal(legacy$code, 1L)
  expect_equal(rotated$code, 2L)
  expect_match(rotated$coordinate_convention, "variance_fraction")
  expect_error(resolve("total_fraction", "sampled", "gaussian", FALSE),
               "requires Gaussian shares, collapsed quality")
  expect_error(resolve("total_fraction", "collapsed", "student_t", FALSE),
               "requires Gaussian shares, collapsed quality")
  expect_error(resolve("total_fraction", "collapsed", "gaussian", TRUE),
               "common residual scale")
})

test_that("Logit and CES data record total/fraction coordinates", {
  n <- 8L
  common <- list(
    prices = seq(1, 1.7, length.out = n),
    quantities = c(16, 13, 11, 9, 18, 14, 10, 7),
    margins = rep(.25, n),
    market = rep(1:2, each = 4), year = rep(1L, n),
    group = rep(1:4, 2), supply_model = 1L,
    use_cutoff = FALSE, use_rho = TRUE, use_hmt = FALSE,
    priors = list(), quality_parameterization = "collapsed",
    share_likelihood = "gaussian", margin_likelihood = "gaussian",
    quality_scope = "market",
    quality_scale_parameterization = "total_fraction"
  )
  logit <- do.call(antitrustBayes:::prepare_stan_data,
                   c(common, list(market_side = "input")))
  ces <- do.call(antitrustBayes:::prepare_ces_stan_data,
                 c(common, list(market_side = "output")))
  for (data in list(logit, ces)) {
    expect_equal(data$quality_scale_parameterization, 2L)
    expect_equal(attr(data, "quality_scale_parameterization"),
                 "total_fraction")
    expect_match(attr(data, "quality_scale_coordinate_convention"),
                 "total_scale")
  }
})

test_that("total/fraction map and Jacobian recover historical scales", {
  scenarios <- rbind(
    c(sigma = .4, tau = .8),
    c(sigma = .01, tau = 1.3),
    c(sigma = .7, tau = .01),
    c(sigma = .012, tau = .008)
  )
  for (i in seq_len(nrow(scenarios))) {
    sigma <- unname(scenarios[i, "sigma"])
    tau <- unname(scenarios[i, "tau"])
    x <- .5 * log(sigma^2 + tau^2)
    y <- 2 * (log(sigma) - log(tau))
    lambda <- plogis(y)
    sigma_new <- exp(x + .5 * plogis(y, log.p = TRUE))
    tau_new <- exp(x + .5 * plogis(-y, log.p = TRUE))
    expect_equal(sigma_new, sigma, tolerance = 1e-12)
    expect_equal(tau_new, tau, tolerance = 1e-12)

    h <- 1e-6
    map <- function(xy) c(
      ell = xy[1] + .5 * plogis(xy[2], log.p = TRUE),
      tau = exp(xy[1] + .5 * plogis(-xy[2], log.p = TRUE))
    )
    J <- cbind(
      (map(c(x + h, y)) - map(c(x - h, y))) / (2 * h),
      (map(c(x, y + h)) - map(c(x, y - h))) / (2 * h)
    )
    expect_equal(abs(det(J)), tau_new / 2, tolerance = 2e-7)
  }
})

test_that("all Stan kernels retain normalized legacy and total/fraction paths", {
  stan_dir <- .ab_test_path("inst", "stan")
  files <- file.path(stan_dir, c(
    "logit_supply.stan", "ces_supply.stan",
    "logit_stackelberg.stan", "ces_stackelberg.stan"
  ))
  for (file in files) {
    src <- paste(readLines(file, warn = FALSE), collapse = "\n")
    expect_match(src,
      "int<lower=1, upper=2> quality_scale_parameterization;", fixed = TRUE)
    expect_match(src, "real sigma_b_strat_economic", fixed = TRUE)
    expect_match(src, "real log_sigma_share_economic", fixed = TRUE)
    expect_match(src,
      "log(sigma_b_strat_economic) - log(2.0) - log(sigma_b_strat)",
      fixed = TRUE)
    expect_match(src,
      "normal_lpdf(log_sigma_share_economic |", fixed = TRUE)
  }
})
